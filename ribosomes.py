import mapping as mapping_tools
import os.path
import gtf
import pysam
import subprocess
import fastq
import mutations
import numpy as np
import matplotlib.pyplot as plt
import Serialize
from itertools import cycle
from collections import Counter, defaultdict

fraction_resolution = 10
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'BlueViolet', 'Gold']

def map_tophat(reads_file_name,
               bowtie2_index,
               gtf_file_name,
               transcriptome_index,
               tophat_dir,
              ):
    tophat_command = ['tophat2',
                      '--GTF', gtf_file_name,
                      '--no-novel-juncs',
                      '--num-threads', '1',
                      '--output-dir', tophat_dir,
                      '--transcriptome-index', transcriptome_index,
                      bowtie2_index,
                      reads_file_name,
                     ]
    # tophat maintains its own logs of everything that is writted to the
    # console, so discard output.
    subprocess.check_output(tophat_command, stderr=subprocess.STDOUT)

def get_aggregate_positions(clean_bam_fn,
                            simple_CDSs,
                            max_gene_length,
                            max_read_length,
                           ):
    edge_overlap = 50
    # To some extent, this is linked to the definition of simple_CDS, which
    # excludes genes that have another gene within 50 bp.

    # position_counts will hold, for each gene, for the lengths 28, 29, and 30,
    # for each position from -2 * edge_overlap to gene_length + edge_overlap,
    # the number of reads of that length starting at that position
    position_counts = []
    gene_names = []
    for gene in simple_CDSs:
        gene_name = gtf.parse_attribute(gene.attribute)['protein_id']
        gene_length = abs(gene.end - gene.start) + 1  
        gene_names.append((gene_name, gene_length))
        shape = (3, 3 * edge_overlap + gene_length)
        position_counts.append(np.zeros(shape, int))
    expression_counts = np.zeros((len(simple_CDSs), 2), int)

    # For from_starts, dimension 2 should be interpretted as going from 
    # -2 * edge_overlap to max_gene_length + edge_overlap.
    # For from_ends, dimension 2 should be interpretted as going from 
    # -(2 * edge_overlap + max_gene_length) to edge_overlap.
    shape = (max_read_length + 1, 3 * edge_overlap + max_gene_length)
    from_starts = np.zeros(shape, int)
    from_ends = np.zeros(shape, int)

    fraction_shape = (max_read_length + 1, fraction_resolution)
    fractions = np.zeros(fraction_shape, int)

    nonunique_mappings = 0

    bamfile = pysam.Samfile(clean_bam_fn, 'rb')
    for g, CDS in enumerate(simple_CDSs):
        max_from_start = CDS.end - CDS.start
        fraction_factor = float(fraction_resolution) / max_from_start

        reads = bamfile.fetch(CDS.seqname,
                              CDS.start - edge_overlap,
                              CDS.end + edge_overlap,
                             )
        for read in reads:
            if read.mapq != 50:
                nonunique_mappings += 1
            else:
                strand = '-' if read.is_reverse else '+'
                
                if strand != CDS.strand:
                    expression_counts[g, 1] += 1
                    continue
                else:
                    expression_counts[g, 0] += 1
                    if CDS.strand == '+':
                        from_start = read.pos - CDS.start
                        # CDS.end is the last base before the stop codon
                        from_end = read.pos - (CDS.end + 1)
                    elif CDS.strand == '-':
                        from_start = CDS.end - (read.aend - 1)
                        # CDS.start is the first base after the stop codon
                        from_end = (CDS.start - 1) - (read.aend - 1)
                
                start_index = 2 * edge_overlap + from_start 
                from_starts[read.qlen, start_index] += 1
                if 28 <= read.qlen <= 30:
                    position_counts[g][read.qlen - 28, start_index] += 1
                # For from_ends, index 2 * edge_overlap + max_gene_length
                # corresponds to a read that starts right at the stop codon.
                end_index = 2 * edge_overlap + max_gene_length + from_end
                if 0 <= end_index < shape[1]:
                    from_ends[read.qlen, end_index] += 1
                else:
                    print "bad from_end index"
                
                #if from_start >= 0:
                #    if from_start == max_from_start:
                #        fraction = fraction_resolution - 1  
                #    else:
                #        fraction = int(fraction_factor * from_start)
                #    if fraction < fraction_resolution:
                #        fractions[read.qlen, fraction] += 1

    return gene_names, position_counts, expression_counts, from_starts, from_ends

def get_extent_positions(clean_bam_fn, extent, genome_index):
    ''' position_counts: an array with rows for length 28, 29, and 30 reads
                         containing the number of uniquely mapped reads of that
                         length starting at that position
    '''
    edge_overlap = 50
    # To some extent, this is linked to the definition of simple_CDS, which
    # excludes genes that have another gene within 50 bp.

    seqname, gene_strand, start, end = extent
    gene_length = abs(end - start) + 1  
    # position_counts will hold, for the lengths 28, 29, and 30,
    # for each position from -2 * edge_overlap to gene_length + edge_overlap,
    # the number of uniquely mapping reads of that length starting at that
    # position
    shape = (3, 3 * edge_overlap + gene_length)
    position_counts = np.zeros(shape, int)
    expression_counts = np.zeros(2, int)
    
    bamfile = pysam.Samfile(clean_bam_fn, 'rb')
    reads = bamfile.fetch(seqname,
                          start - edge_overlap,
                          end + edge_overlap,
                         )
    for read in reads:
        if read.mapq != 50:
            pass
        else:
            strand = '-' if read.is_reverse else '+'
            
            if strand != gene_strand:
                expression_counts[1] += 1
                continue
            else:
                expression_counts[0] += 1
                if gene_strand == '+':
                    from_start = read.pos - start
                    # end is the last base before the stop codon
                    from_end = read.pos - (end + 1)
                elif gene_strand == '-':
                    from_start = end - (read.aend - 1)
                    # start is the first base after the stop codon
                    from_end = (start - 1) - (read.aend - 1)
            
            start_index = 2 * edge_overlap + from_start 
            if 28 <= read.qlen <= 30:
                position_counts[read.qlen - 28, start_index] += 1
                
    return position_counts, expression_counts

def determine_ambiguity_of_positions(extent, genome_index):
    edge_overlap = 50

    seqname, gene_strand, start, end = extent
    gene_length = abs(end - start) + 1  
    
    ref_seq = str(genome_index[seqname].seq)[start - 2 * edge_overlap:end + edge_overlap]
    
    shape = (3, 3 * edge_overlap + gene_length)
    position_ambiguity = np.empty(shape, int)
    for length in [28, 29, 30]:
        for start in range(3 * edge_overlap + gene_length - length):
            # Explanation of encoding:
            # 0 means a fragment of this length starting at this position would
            # end in an A and therefore can't be the result of poly-A trimming.
            # 1 means a fragment of this length starting at this position would
            # end in a non-A followed by an A, so the read produced can't be
            # distinguished from a read of greater length starting at the same
            # position.
            # 2 means a fragment of this length starting at this position would
            # end in a non-A followed by another non-A and is therefore
            # unambiguous.
            if ref_seq[start + length - 1].upper() == 'A':
                print length - 28, start 
                position_ambiguity[length - 28, start] = 0
            elif ref_seq[start + length].upper() == 'A':
                position_ambiguity[length - 28, start] = 1
            else:
                position_ambiguity[length - 28, start] = 2
                
    return position_ambiguity

def plot_starts_and_ends(from_starts_list, from_ends_list, names, fig_file_name):
    edge_overlap = 50

    fig, (start_ax, end_ax) = plt.subplots(2, 1, figsize=(12, 16))

    experiments = zip(from_starts_list, from_ends_list, names)
    for from_starts, from_ends, name in experiments:
        starts_xs = np.arange(-2 * edge_overlap, 1000)
        ends_xs = np.arange(-1000, edge_overlap)
        for length in [28, 29, 30]:
            starts = np.true_divide(from_starts[length], from_starts[length].sum())
            ends = np.true_divide(from_ends[length], from_starts[length].sum())
            start_ax.plot(starts_xs, starts[:len(starts_xs)],
                          '.-',
                          label=name + '_{0}'.format(length),
                         )
            end_ax.plot(ends_xs, ends[-len(ends_xs):],
                        '.-',
                        label=name + '_{0}'.format(length),
                       )

    start_ax.set_xlim(-20, 20)
    start_ax.set_xlabel('Position of read relative to start of CDS')
    start_ax.set_ylabel('Fraction of uniquely mapped reads of specific length')
    
    end_ax.set_xlim(-35, 5)
    end_ax.set_xlabel('Position of read relative to stop codon')
    end_ax.set_ylabel('Fraction of uniquely mapped readsof specific length')

    leg = start_ax.legend(loc='upper right', fancybox=True)
    leg.get_frame().set_alpha(0.5)
    leg = end_ax.legend(loc='upper right', fancybox=True)
    leg.get_frame().set_alpha(0.5)


    fig.savefig(fig_file_name)

def plot_aggregate(rpf_positions_list, names, fig_file_name):
    edge_overlap = 50

    min_codons = 200
    start_at_codon = -8
    end_at_codon = min_codons
    codons = np.arange(start_at_codon, end_at_codon)
    codon_starts = np.arange(2 * edge_overlap + 3 * start_at_codon,
                             2 * edge_overlap + 3 * end_at_codon,
                             3,
                            )
    
    fig, ax = plt.subplots(figsize=(12, 16))

    for name, rpf_positions in zip(names, rpf_positions_list):
        counts_in_long_genes = np.zeros((3, 2 * edge_overlap + min_codons * 3), float)
        for (gene_name, length), positions in zip(*rpf_positions):
            if length > min_codons * 3:
                gene_prefix = positions[:, :2 * edge_overlap + min_codons * 3]
                row_sums = gene_prefix.sum(axis=1)
                zeros_removed = np.maximum(row_sums, 1)
                normalized = np.true_divide(gene_prefix, zeros_removed[:, np.newaxis])
                counts_in_long_genes += normalized
    
        counts = [counts_in_long_genes[0][c] for c in codon_starts]
        ax.plot(codons, counts, '.-', label=name)

    ax.set_xlim(start_at_codon, end_at_codon)
    leg = ax.legend(loc='upper right', fancybox=True)
    leg.get_frame().set_alpha(0.5)
    fig.savefig(fig_file_name)

def plot_gene(rpf_positions_list, gene_name):
    edge_overlap = 50

    fig, ax = plt.subplots(figsize=(12, 16))

    for rpf_positions in rpf_positions_list:
        by_name = {name: counts for (name, length), counts in zip(*rpf_positions)}
        lengths = {name: length for (name, length), counts in zip(*rpf_positions)}
        length = lengths[gene_name]
        counts = by_name[gene_name]
        
        start_at_codon = -8
        end_at_codon = length / 3
        codons = np.arange(start_at_codon, end_at_codon)
        codon_starts = np.arange(2 * edge_overlap + 3 * start_at_codon,
                                 2 * edge_overlap + 3 * end_at_codon,
                                 3,
                                )
        counts = [counts[0][c] for c in codon_starts]
        ax.plot(codons, counts, '.-')
        ax.set_xlim(start_at_codon, end_at_codon)

def scatter_positions(rpf_positions_list):
    edge_overlap = 50

    fig, ax = plt.subplots(figsize=(12, 16))

    counts_list = []

    for rpf_positions in rpf_positions_list:
        experiment_counts = []
        for (name, length), counts in zip(*rpf_positions): 
            start_at_codon = -8
            end_at_codon = length / 3
            codons = np.arange(start_at_codon, end_at_codon)
            codon_starts = np.arange(2 * edge_overlap + 3 * start_at_codon,
                                     2 * edge_overlap + 3 * end_at_codon,
                                     3,
                                    )
            counts = [counts[0][c] for c in codon_starts]
            experiment_counts.extend(counts)
        counts_list.append(experiment_counts)
      
    ax.scatter(counts_list[0], counts_list[1], s=1)

def plot_density(from_starts_fns, from_ends_fns, names, read_lengths):
    from_starts_list = [np.loadtxt(fn) for fn in from_starts_fns]
    from_ends_list = [np.loadtxt(fn) for fn in from_ends_fns]

    fig, (start_ax, end_ax) = plt.subplots(2, 1, figsize=(12, 16))

    experiments = zip(from_starts_list, from_ends_list, names, read_lengths)
    for from_starts, from_ends, name, read_length in experiments:
        start_xs = np.arange(501)
        end_xs = np.arange(-500, 1)
        start_counts = from_starts[28][read_length + 1 - 12::3][:501]
        print len(start_xs), len(start_counts)
        end_counts = from_ends[28][-read_length - 15::-3][:501][::-1]
        start_ax.plot(start_xs, start_counts,
                      '.-',
                     )
        end_ax.plot(end_xs, end_counts,
                    '.-',
                   )

    #start_ax.set_xlim(-20, 20)
    #start_ax.set_xlabel('Position of read relative to start of CDS')
    #start_ax.set_ylabel('Fraction of uniquely mapped reads of specific length')
    #
    #end_ax.set_xlim(-35, 5)
    #end_ax.set_xlabel('Position of read relative to stop codon')
    #end_ax.set_ylabel('Fraction of uniquely mapped readsof specific length')

    #leg = start_ax.legend(loc='upper right', fancybox=True)
    #leg.get_frame().set_alpha(0.5)
    #leg = end_ax.legend(loc='upper right', fancybox=True)
    #leg.get_frame().set_alpha(0.5)

def plot_mRNA(from_starts_fns, from_ends_fns, fractions_fns, names, read_lengths, fragment_lengths):
    from_starts_list = [np.loadtxt(fn) for fn in from_starts_fns]
    from_ends_list = [np.loadtxt(fn) for fn in from_ends_fns]
    fractions_list = [np.loadtxt(fn) for fn in fractions_fns]

    #fig, (start_ax, end_ax, fraction_ax) = plt.subplots(3, 1, figsize=(12, 24))
    fig, fraction_ax = plt.subplots(1, 1, figsize=(12, 8))

    for from_starts, from_ends, fractions, name, read_length in zip(from_starts_list, from_ends_list, fractions_list, names, read_lengths):
        start_xs = np.arange(-read_length, 10000 + 1)
        end_xs = np.arange(-10000, read_length + 1)
        fraction_xs = np.linspace(0.5 / fraction_resolution,
                                  1 - 0.5 / fraction_resolution,
                                  fraction_resolution,
                                 )
        
        style = '.-'

        starts = np.true_divide(from_starts[fragment_lengths].sum(axis=0),
                                from_starts[fragment_lengths].sum(),
                               )
        #start_ax.plot(start_xs, starts, style, label=name)

        ends = np.true_divide(from_ends[fragment_lengths].sum(axis=0),
                              from_ends[fragment_lengths].sum(),
                             )
        #end_ax.plot(end_xs, ends, style, label=name)
        
        fractions = np.true_divide(fractions[fragment_lengths].sum(axis=0),
                                   fractions[fragment_lengths].sum(),
                                  )
        fraction_ax.plot(fraction_xs, fractions, 'o-', label=name)

    #start_ax.set_xlim(-max(read_lengths), 10000)
    #start_ax.legend()
    #
    #end_ax.set_xlim(-10000, max(read_lengths))
    #end_ax.legend()
    
    fraction_ax.set_xlim(0, 1)
    fraction_ax.set_ylim(ymin=0)
    fraction_ax.set_ylabel('Fraction of reads')
    fraction_ax.set_xlabel('Fraction of distance along CDS')
    fraction_ax.set_title('3\' bias in mRNA reads')
    leg = fraction_ax.legend(loc='upper left', fancybox=True)
    leg.get_frame().set_alpha(0.5)

def plot_frames(summary_fns, names):
    frames = {}
    for summary_fn, name in zip(summary_fns, names):
        counts = np.loadtxt(summary_fn, usecols=range(4, 13))
        frames[name, 28] = counts[:, 0:3].sum(axis=0)
        frames[name, 29] = counts[:, 3:6].sum(axis=0)
        frames[name, 30] = counts[:, 6:9].sum(axis=0)

    #fig, (ax_28, ax_29, ax_30) = plt.subplots(3, 1, figsize=(6, 12))
    #for ax, length in zip([ax_28, ax_29, ax_30], [28, 29, 30]):
    fig, (ax_28, ax_29) = plt.subplots(2, 1, figsize=(6, 12))
    for ax, length in zip([ax_28, ax_29], [28, 29]):
        rates_list = [frames[name, length] / frames[name, length].sum() for name in names]
        tick_labels = ['0', '1', '2']

        grouped_bar(rates_list, names, tick_labels, ax)
        ax.set_title('Length {0}'.format(length))
        ax.set_xlabel('Frame')
        ax.set_ylabel('Fraction of uniquely mapped reads')

    return counts

def grouped_bar(rates_list, names, tick_labels, ax):
    N = len(rates_list[0])

    group_starts = np.arange(N)
    gap = 1 # in multiples of bar width

    width = 1. / (len(rates_list) + gap)

    zipped = zip(rates_list, colors, names)
    for i, (rates, color, name) in enumerate(zipped):
        ax.bar(group_starts + i * width,
               rates,
               width,
               color=color,
               label=name,
              )

    ax.set_xlim(xmin=-width / 2, xmax=N - width / 2)

    ax.set_xticks(group_starts + (width * len(rates_list) / 2))
    ax.set_xticklabels(tick_labels)
    ax.xaxis.set_ticks_position('none')
    ax.set_ylim(0, 1)

    leg = ax.legend(loc='upper right', fancybox=True)
    leg.get_frame().set_alpha(0.5)

def compare(summary_fns, clean_length_fns, names):
    reads_list = [np.loadtxt(fn, dtype=int)[28:31].sum() for fn in clean_length_fns]

    sample_counts = defaultdict(lambda : defaultdict(float))

    for summary_fn, name, reads in zip(summary_fns, names, reads_list):
        for line in open(summary_fn):
            gene, length, count, _ = line.split('\t', 3)
            length = int(length)
            count = int(count)
            sample_counts[gene][name] = 1e3 * 1e6 * float(count) / (length * reads)
            #sample_counts[gene][name] = 1e6 * float(count) / (reads)

    fig = plt.figure(figsize=(12, 12))
    for r in range(len(names)):
        for c in range(r + 1, len(names)):
            ax = fig.add_subplot(len(names) - 1, len(names) - 1, r * (len(names) - 1) + (c - 1) + 1)
            first = names[c]
            second = names[r]
            genes = sample_counts.keys()
            xs = [sample_counts[gene][first] for gene in genes]
            ys = [sample_counts[gene][second] for gene in genes]

            print sum(xs)
            print sum(ys)
            print
            
            ax.scatter(xs, ys, s=1)
            ax.set_yscale('log')
            ax.set_xscale('log')
            ax.set_xlim(1e-1, 5e4)
            ax.set_ylim(1e-1, 5e4)
            ax.plot([1e-2, 1e5], [1e-2, 1e5], '-', color='red', alpha=0.2)
            ax.set_xlabel(first)
            ax.set_ylabel(second)

    return sample_counts

def plot_frameshifts(rpf_counts_list,
                     position_ambiguity_list,
                     edge_overlap,
                     gene_name,
                     gene_length,
                    ):
    ambiguity_to_color = {0: 'red',
                          1: 'blue',
                          2: 'green',
                         }

    length_data = zip([28, 29, 30], rpf_counts_list, position_ambiguity_list)
    for fragment_length, rpf_counts, position_ambiguity in length_data:
        start_at_codon = -8
        # codon_starts[i] is the index into a position array at which codon i
        # starts.
        codon_starts = np.arange(2 * edge_overlap + (start_at_codon * 3),
                                 2 * edge_overlap + gene_length,
                                 3,
                                )
        codon_numbers = start_at_codon + np.arange(len(codon_starts))
        
        # frame_counts_list[i, j] will be the number of RPF's starting at frame i of
        # codon j
        frame_counts_list = np.zeros((3, len(codon_starts)), int)
        # frame_colors_list will be used to visualize the ambiguity of each position
        frame_colors_list = [['']*len(codon_starts) for frame in range(3)]

        for c, codon_start in enumerate(codon_starts):
            for frame in range(3):
                frame_counts_list[frame, c] = rpf_counts[codon_start + frame]
                ambiguity = position_ambiguity[codon_start + frame]
                frame_colors_list[frame][c] = ambiguity_to_color[ambiguity]

        frames_so_far = frame_counts_list.cumsum(axis=1)
        fraction_frames_so_far = np.true_divide(frames_so_far, frames_so_far.sum(axis=0))

        frames_remaining = np.fliplr(np.fliplr(frame_counts_list).cumsum(axis=1))
        fraction_frames_remaining = np.true_divide(frames_remaining, frames_remaining.sum(axis=0))
        
        fig, axs = plt.subplots(4, 1, sharex=True)
        cumulative_ax = axs[0]
        frame_axs = axs[1:]

        for frame, (ax, frame_counts, frame_colors) in enumerate(zip(frame_axs, frame_counts_list, frame_colors_list)):
            ax.scatter(codon_numbers, frame_counts, s=10, c=frame_colors, linewidths=0)
            ax.set_ylim(-1, frame_counts_list.max() + 1)
            ax.set_xlim(codon_numbers[0], codon_numbers[-1])
            ax.set_title('Frame {0}'.format(frame))

        for so_far, remaining, color in zip(fraction_frames_so_far, fraction_frames_remaining, colors):
            cumulative_ax.plot(codon_numbers, so_far, color=color)
            cumulative_ax.plot(codon_numbers, remaining, color=color, linestyle='--')
            #difference = so_far - remaining
            #cumulative_ax.plot(codon_numbers, difference, color=color, linestyle=':')
            #difference = remaining - so_far
            #cumulative_ax.plot(codon_numbers, difference, color=color, linestyle=':')
        cumulative_ax.set_xlim(codon_numbers[0])
        
        fig.suptitle('{0} - length {1} fragments'.format(gene_name, fragment_length))

    return codon_numbers, frame_counts_list, fraction_frames_so_far, fraction_frames_remaining


def plot_all_aggregates():
    experiments = [
        ('geranshenko1', '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep1_foot/results/Initial_rep1_foot_rpf_positions.txt'),
        ('geranshenko2', '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep2_foot/results/Initial_rep2_foot_rpf_positions.txt'),
        ('ingolia1', '/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-1/results/Footprints-rich-1_rpf_positions.txt'),
        ('R98S', '/home/jah/projects/arlen/experiments/belgium_8_6_13/R98S_cDNA_sample/results/R98S_cDNA_sample_rpf_positions.txt'),
        ('suppressed', '/home/jah/projects/arlen/experiments/belgium_8_6_13/Suppressed_R98S_cDNA_sample/results/Suppressed_R98S_cDNA_sample_rpf_positions.txt'),
        ('WT',  '/home/jah/projects/arlen/experiments/belgium_8_6_13/WT_cDNA_sample/results/WT_cDNA_sample_rpf_positions.txt'),
    ]

    names = [name for name, _ in experiments]
    rpf_positions_list = [Serialize.read_file(fn, 'rpf_positions') for _, fn in experiments]
    fig_file_name = '/home/jah/projects/arlen/results/compare_aggregate_positions.pdf'
    plot_aggregate(rpf_positions_list, names, fig_file_name)
    
if __name__ == '__main__':
    genome_index = mapping_tools.get_genome_index('/home/jah/projects/arlen/data/organisms/saccharomyces_cerevisiae/genome', explicit_path=True)

    #clean_bam_fn = '/home/jah/projects/arlen/experiments/belgium_8_6_13/WT_cDNA_sample/results/WT_cDNA_sample_clean.bam'
    #clean_bam_fn = '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep1_foot/results/Initial_rep1_foot_clean.bam'
    clean_bam_fn = '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep2_foot/results/Initial_rep2_foot_clean.bam'
    
    gtf_fn = '/home/jah/projects/arlen/data/organisms/saccharomyces_cerevisiae/transcriptome/genes.gtf'
    
    gene_name = 'YOR239W'
    #gene_name = 'YPL052W'
    #gene_name = 'YAL053W'
    
    extent = gtf.get_extent_by_name(gtf_fn, gene_name)
    gene_length = extent[3] - extent[2] + 1
    position_counts, expression_counts = get_extent_positions(clean_bam_fn, extent, genome_index)
    position_ambiguity = determine_ambiguity_of_positions(extent, genome_index)
    codon_numbers, frames, so_far, remaining = plot_frameshifts(position_counts,
                                                                position_ambiguity,
                                                                50,
                                                                gene_name,
                                                                gene_length,
                                                               )
