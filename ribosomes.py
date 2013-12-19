import mapping as mapping_tools
import os.path
import gtf
import recycling
import pysam
import sam
import subprocess
import fastq
import mutations
from mutations import base_order, base_to_index, base_to_complement_index
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import Serialize
from itertools import cycle, izip
from collections import Counter, defaultdict

fraction_resolution = 10
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'BlueViolet', 'Gold']
    
# To some extent, this is linked to the definition of simple_CDS, which
# excludes genes that have another gene within 50 bp.
edge_overlap = 50

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

def get_read_positions(clean_bam_fn,
                       simple_CDSs,
                       relevant_lengths,
                       allow_opposite_strand):
    genes = {}
    for gene in simple_CDSs:
        gene_name = gtf.parse_attribute(gene.attribute)['protein_id']
        CDS_length = abs(gene.end - gene.start) + 1  
        shape = (3, 3 * edge_overlap + CDS_length)
        # position_counts will hold, for each valid length, for each position from 
        # -2 * edge_overlap to CDS_length + edge_overlap, the number of reads of
        # that length starting at that position.
        position_counts = {l: np.zeros(3 * edge_overlap + CDS_length, int)
                           for l in relevant_lengths}
        genes[gene_name] = {'CDS_length': CDS_length,
                            'position_counts': position_counts,
                            'expression': np.zeros(2, int),
                            'nonunique': 0,
                           }

    bamfile = pysam.Samfile(clean_bam_fn, 'rb')
    for CDS in simple_CDSs:
        gene_name = gtf.parse_attribute(CDS.attribute)['protein_id']
        position_counts = genes[gene_name]['position_counts']
        expression_counts = genes[gene_name]['expression']

        overlapping_reads = bamfile.fetch(CDS.seqname,
                                          CDS.start - edge_overlap,
                                          CDS.end + edge_overlap,
                                         )
        for read in overlapping_reads:
            if read.mapq != 50:
                genes[gene_name]['nonunique'] += 1
            else:
                strand = '-' if read.is_reverse else '+'

                if strand != CDS.strand:
                    expression_counts[1] += 1
                    if not allow_opposite_strand:
                        # For RPF reads, we don't want to record positions
                        # of reads on the opposite strand.
                        continue
                else:
                    expression_counts[0] += 1

                if read.qlen in relevant_lengths:
                    if CDS.strand == '+':
                        from_start = read.pos - CDS.start
                    elif CDS.strand == '-':
                        from_start = CDS.end - (read.aend - 1)

                    start_index = 2 * edge_overlap + from_start 
                    position_counts[read.qlen][start_index] += 1

    return genes

def old_aggregate(clean_bam_fn,
                  simple_CDSs,
                  max_gene_length,
                  max_read_length,
                  valid_lengths,
                 ):
    genes = {}
    for gene in simple_CDSs:
        gene_name = gtf.parse_attribute(gene.attribute)['protein_id']
        CDS_length = abs(gene.end - gene.start) + 1  
        shape = (3, 3 * edge_overlap + CDS_length)
        # counts will hold, for the lengths 28, 29, and 30,
        # for each position from -2 * edge_overlap to CDS_length + edge_overlap,
        # the number of reads of that length starting at that position
        genes[gene_name] = {'CDS_length': CDS_length,
                            'counts': dict(zip([28, 29, 30], np.zeros(shape, int))),
                            'expression': np.zeros(2, int),
                           }

    # For from_starts, dimension 2 should be interpretted as going from 
    # -2 * edge_overlap to max_gene_length + edge_overlap.
    # For from_ends, dimension 2 should be interpretted as going from 
    # -(2 * edge_overlap + max_gene_length) to edge_overlap.
    shape = (max_read_length + 1, 3 * edge_overlap + max_gene_length)
    from_starts = np.zeros(shape, int)
    from_ends = np.zeros(shape, int)

    nonunique_mappings = 0

    bamfile = pysam.Samfile(clean_bam_fn, 'rb')
    for CDS in simple_CDSs:
        gene_name = gtf.parse_attribute(CDS.attribute)['protein_id']
        position_counts = genes[gene_name]['counts']
        expression_counts = genes[gene_name]['expression']

        overlapping_reads = bamfile.fetch(CDS.seqname,
                                          CDS.start - edge_overlap,
                                          CDS.end + edge_overlap,
                                         )
        for read in overlapping_reads:
            if read.mapq != 50:
                nonunique_mappings += 1
            else:
                strand = '-' if read.is_reverse else '+'

                if strand != CDS.strand:
                    expression_counts[1] += 1
                    continue
                else:
                    expression_counts[0] += 1
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
                # For from_ends, index 2 * edge_overlap + max_gene_length
                # corresponds to a read that starts right at the stop codon.
                end_index = 2 * edge_overlap + max_gene_length + from_end
                if 0 <= end_index < shape[1]:
                    from_ends[read.qlen, end_index] += 1
                else:
                    print "bad from_end index"
                if 28 <= read.qlen <= 30:
                    position_counts[read.qlen][start_index] += 1

    return genes, from_starts, from_ends

def get_extent_positions(clean_bam_fn, extent, genome_index):
    ''' position_counts: an array with rows for length 28, 29, and 30 reads
                         containing the number of uniquely mapped reads of that
                         length starting at that position
    '''
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
    overlapping_reads = bamfile.fetch(seqname,
                                      start - edge_overlap,
                                      end + edge_overlap,
                                     )
    for read in overlapping_reads:
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

def get_codon_counts(gene, stringent=True):
    if stringent:
        length_28s = gene['position_counts'][28]
        A_site_offset = 15
        start_index = 2 * edge_overlap - A_site_offset
        end_index = -(edge_overlap + A_site_offset)
        in_frames = length_28s[start_index:end_index:3]
        return in_frames
    else:
        length_28s = gene['position_counts'][28]
        A_site_offset = 15
        start_index = 2 * edge_overlap - A_site_offset
        end_index = -(edge_overlap + A_site_offset)
        counts_28 = length_28s[start_index:end_index:3] + \
                    length_28s[start_index-1:end_index-1:3] + \
                    length_28s[start_index+1:end_index+1:3]
        
        length_29s = gene['position_counts'][29]
        A_site_offset = 15
        start_index = 2 * edge_overlap - A_site_offset
        end_index = -(edge_overlap + A_site_offset)
        counts_29 = length_29s[start_index:end_index:3] + \
                    length_29s[start_index-1:end_index-1:3] + \
                    length_29s[start_index-2:end_index-2:3]

        length_30s = gene['position_counts'][30]
        A_site_offset = 16
        start_index = 2 * edge_overlap - A_site_offset
        end_index = -(edge_overlap + A_site_offset)
        counts_30 = length_30s[start_index:end_index:3] + \
                    length_30s[start_index-1:end_index-1:3] + \
                    length_30s[start_index+1:end_index+1:3]

        return counts_28 + counts_29 + counts_30

def get_raw_counts(rpf_positions_dict):
    raw_counts = {}
    for gene_name in rpf_positions_dict:
        counts = get_codon_counts(rpf_positions_dict[gene_name])
        raw_counts[gene_name] = counts.sum()

    return raw_counts

def get_RPKMs(genes_dict):
    RPKMs = {}
    total_mapped_reads = 0
    for gene_name in genes_dict:
        counts = get_codon_counts(genes_dict[gene_name])
        length = genes_dict[gene_name]['CDS_length']
        reads = counts.sum()
        RPKMs[gene_name] = float(reads) / length
        total_mapped_reads += reads

    for gene_name in RPKMs:
        RPKMs[gene_name] = 1.e9 * RPKMs[gene_name] / total_mapped_reads

    return RPKMs

def get_TPMs(rpf_positions_dict):
    TPMs = {}
    T = 0
    for gene_name in rpf_positions_dict:
        counts = get_codon_counts(rpf_positions_dict[gene_name])
        length = rpf_positions_dict[gene_name]['CDS_length']
        reads = counts.sum()
        TPMs[gene_name] = float(reads) / length
        T += TPMs[gene_name]

    for gene_name in rpf_positions_dict:
        length = rpf_positions_dict[gene_name]['CDS_length']
        TPMs[gene_name] = 1.e6 * TPMs[gene_name] / T

    return TPMs

def determine_ambiguity_of_positions(extent, genome_index):
    edge_overlap = 50

    seqname, gene_strand, start, end = extent
    gene_length = abs(end - start) + 1  
    
    ref_seq = str(genome_index[seqname].seq)[start - 2 * edge_overlap:end + edge_overlap + 1]
    
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

def plot_starts_and_ends_new(from_starts_list, from_ends_list, names, fig_file_name):
    edge_overlap = 50

    start_fig, start_axs = plt.subplots(3, 1, figsize=(12, 16), sharex=True)
    end_fig, end_axs = plt.subplots(3, 1, figsize=(12, 16), sharex=True)

    experiments = zip(from_starts_list, names)
    for from_starts, name in experiments:
        starts_xs = np.arange(-2 * edge_overlap, 1000)
        for length, ax in zip([28, 29, 30], start_axs):
            starts = np.true_divide(from_starts[length], from_starts[length].sum())
            ax.plot(starts_xs, starts[:len(starts_xs)], '.-', label=name)
            ax.set_ylabel('Fraction of uniquely mapped reads of specific length')
            ax.set_title('Length {0} fragments'.format(length))
    
    experiments = zip(from_ends_list, names)
    for from_ends, name in experiments:
        ends_xs = np.arange(-1000, edge_overlap)
        for length, ax in zip([28, 29, 30], end_axs):
            ends = np.true_divide(from_ends[length], from_ends[length].sum())
            ax.plot(ends_xs, ends[-len(ends_xs):], '.-', label=name)
            ax.set_ylabel('Fraction of uniquely mapped reads of specific length')
            ax.set_title('Length {0} fragments'.format(length))

    leg = start_axs[0].legend(loc='upper right', fancybox=True)
    leg.get_frame().set_alpha(0.5)
    start_axs[-1].set_xlim(-20, 20)
    start_axs[-1].set_xlabel('Position of read relative to start of CDS')

    leg = end_axs[0].legend(loc='upper right', fancybox=True)
    leg.get_frame().set_alpha(0.5)
    end_axs[-1].set_xlim(-35, 5)
    end_axs[-1].set_xlabel('Position of read relative to start of CDS')

def read_codon_counts_file(fn, reports_P_site=False):
    genes = {}
    for line in open(fn):
        name, values = line.strip().split('\t', 1)
        values = np.fromstring(values, dtype=float, sep='\t')
        # If the P site is reported, shift over by one to convert to A site.
        if reports_P_site:
            values = np.concatenate(([0], values))
        genes[name] = values
            
    return genes

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

def compare():
    bartel_RPF, bartel_mRNA = read_bartel_file('/home/jah/projects/arlen/experiments/plotkin/RPKM_15aug30stop.txt') 
    names = [
        'Ingolia_1',
        'Ingolia_2',
        #'Brar',
        #'Gerashchenko',
        'UT_WT',
        'UT_R98S',
        'UT_Suppressed_R98S',
    ]

    #RPF_fns = ['/home/jah/projects/arlen/results/{0}_RPF_rpkm.txt'.format(name)
    #           for name in names]
    mRNA_fns = ['/home/jah/projects/arlen/results/{0}_mRNA_tpm.txt'.format(name)
                for name in names]
    #RPF_rpkms = [read_RPKMs_file(fn) for fn in RPF_fns]
    mRNA_rpkms = [read_RPKMs_file(fn) for fn in mRNA_fns]

    #names += ['Bartel']
    #RPF_rpkms += [bartel_RPF]
    #mRNA_rpkms += [bartel_mRNA]

    # Get the genes common to bartel's and mine.
    gene_names = list(set(bartel_RPF) & set(mRNA_rpkms[0]))
    
    fig = plt.figure(figsize=(12, 12))
    for r in range(len(names)):
        for c in range(r, len(names)):
            ax = fig.add_subplot(len(names),
                                 len(names),
                                 r * len(names) + c + 1,
                                )
            first = names[c]
            second = names[r]
            xs = [mRNA_rpkms[c][gene_name] for gene_name in gene_names]
            ys = [mRNA_rpkms[r][gene_name] for gene_name in gene_names]

            print first, second, scipy.stats.pearsonr(xs, ys)

            ax.scatter(xs, ys, s=1)
            ax.set_yscale('log')
            ax.set_xscale('log')
            ax.set_xlim(1e-1, 5e4)
            ax.set_ylim(1e-1, 5e4)
            ax.plot([1e-2, 1e5], [1e-2, 1e5], '-', color='red', alpha=0.2)
            ax.set_xlabel(first)
            ax.set_ylabel(second)

def get_ratios(first, second):
    assert set(first) == set(second)
    ratios = {key: np.divide(float(first[key]), second[key]) for key in first}
    return ratios

def plot_frameshifts(rpf_counts_list,
                     position_ambiguity_list,
                     edge_overlap,
                     gene_name,
                     gene_length,
                    ):
    ambiguity_to_color = {0: 'red',
                          1: 'green',
                          2: 'black',
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
            ax.scatter(codon_numbers, frame_counts, s=20, c=frame_colors, linewidths=0)
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

def plot_RPKMs():
    experiments = [
        ('geranshenko1', '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep1_foot/results/Initial_rep1_foot_rpf_positions.txt'),
        #('geranshenko2', '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep2_foot/results/Initial_rep2_foot_rpf_positions.txt'),
        ('ingolia1', '/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-1/results/Footprints-rich-1_rpf_positions.txt'),
        ('ingolia2', '/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-2/results/Footprints-rich-2_rpf_positions.txt'),
        ('ingolia1_mRNA', '/home/jah/projects/arlen/experiments/ingolia_science/mRNA-rich-1/results/mRNA-rich-1_rpf_positions.txt'),
        #('R98S', '/home/jah/projects/arlen/experiments/belgium_8_6_13/R98S_cDNA_sample/results/R98S_cDNA_sample_rpf_positions.txt'),
        #('suppressed', '/home/jah/projects/arlen/experiments/belgium_8_6_13/Suppressed_R98S_cDNA_sample/results/Suppressed_R98S_cDNA_sample_rpf_positions.txt'),
        #('WT',  '/home/jah/projects/arlen/experiments/belgium_8_6_13/WT_cDNA_sample/results/WT_cDNA_sample_rpf_positions.txt'),
    ]

    names = [name for name, _ in experiments]
    rpf_positions_lists = [Serialize.read_file(fn, 'rpf_positions') for _, fn in experiments]
    RPKMs_list = [get_TPMs(rpf_positions_list) for rpf_positions_list in rpf_positions_lists]
    vals = [[val for name, val in sorted(RPKMs.items())] for RPKMs in RPKMs_list]
    
    fig = plt.figure(figsize=(12, 12))
    for r in range(len(names)):
        for c in range(r + 1, len(names)):
            ax = fig.add_subplot(len(names) - 1, len(names) - 1, r * (len(names) - 1) + (c - 1) + 1)
            first = names[c]
            second = names[r]
            xs = vals[r]
            ys = vals[c]

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

def density_length_correlation():
    experiments = [
        ('ingolia1',
         '/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-1/results/Footprints-rich-1_rpf_positions.txt',
         '/home/jah/projects/arlen/experiments/ingolia_science/mRNA-rich-1/results/mRNA-rich-1_rpf_positions.txt',
        ),
        #('geranshenko1',
        # '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep1_foot/results/Initial_rep1_foot_rpf_positions.txt',
        # '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep1_mRNA/results/Initial_rep1_mRNA_rpf_positions.txt',
        #),
        ('R98S',
         '/home/jah/projects/arlen/experiments/belgium_8_6_13/R98S_cDNA_sample/results/R98S_cDNA_sample_rpf_positions.txt',
         '/home/jah/projects/arlen/experiments/belgium_8_6_13/R98S_cDNA_mRNA/results/R98S_cDNA_mRNA_rpf_positions.txt',
        ),
    ]

    for name, rpf_fn, mRNA_fn in experiments:
        rpfs = Serialize.read_file(rpf_fn, 'rpf_positions')
        rpf_TPMs = get_TPMs(rpfs)
        mRNAs = Serialize.read_file(mRNA_fn, 'rpf_positions')
        mRNA_TPMs = get_TPMs(mRNAs)
        ratios = get_ratios(rpf_TPMs, mRNA_TPMs)
        ratios_vals = [np.log2(val) for _, val in sorted(ratios.items())]
        lengths = [rpfs[gene_name]['CDS_length'] for gene_name in rpfs]
        fig, ax = plt.subplots()
        ax.scatter(lengths, ratios_vals, s=1)
        print name
        print scipy.stats.pearsonr(lengths, ratios_vals)
        print scipy.stats.spearmanr(lengths, ratios_vals)
        ax.set_xlabel('CDS length')
        ax.set_ylabel('Translational efficiency (log_2 ratio of TPMs)')
        ax.set_title(name)
    
def plot_all_starts_ends():
    experiments = [
        ('geranshenko1', 
         '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep1_foot/results/Initial_rep1_foot_from_starts.txt',
         '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep1_foot/results/Initial_rep1_foot_from_ends.txt',
        ),
        ('ingolia1',
         '/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-1/results/Footprints-rich-1_from_starts.txt',
         '/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-1/results/Footprints-rich-1_from_ends.txt',
        ),
        ('ingolia2',
         '/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-2/results/Footprints-rich-2_from_starts.txt',
         '/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-2/results/Footprints-rich-2_from_ends.txt',
        ),
        ('R98S',
         '/home/jah/projects/arlen/experiments/belgium_8_6_13/R98S_cDNA_sample/results/R98S_cDNA_sample_from_starts.txt',
         '/home/jah/projects/arlen/experiments/belgium_8_6_13/R98S_cDNA_sample/results/R98S_cDNA_sample_from_ends.txt',
        ),
        ('suppressed',
         '/home/jah/projects/arlen/experiments/belgium_8_6_13/Suppressed_R98S_cDNA_sample/results/Suppressed_R98S_cDNA_sample_from_starts.txt',
         '/home/jah/projects/arlen/experiments/belgium_8_6_13/Suppressed_R98S_cDNA_sample/results/Suppressed_R98S_cDNA_sample_from_ends.txt',
        ),
        ('WT',
         '/home/jah/projects/arlen/experiments/belgium_8_6_13/WT_cDNA_sample/results/WT_cDNA_sample_from_starts.txt',
         '/home/jah/projects/arlen/experiments/belgium_8_6_13/WT_cDNA_sample/results/WT_cDNA_sample_from_ends.txt',
        ),
    ]

    names = [name for name, _, _ in experiments]
    from_starts_list = [Serialize.read_file(fn, 'array') for _, fn, _ in experiments]
    from_ends_list = [Serialize.read_file(fn, 'array') for _, _, fn in experiments]
    fig_file_name = '/home/jah/projects/arlen/results/compare_starts_and_ends.pdf'
    plot_starts_and_ends_new(from_starts_list, from_ends_list, names, fig_file_name)

def error_profile(bam_file_name, simple_CDSs, fastq_type):
    edge_overlap = 50
    type_shape = (7,
                  50,
                  fastq.MAX_EXPECTED_QUAL + 1,
                  len(base_order),
                  len(base_order),
                 )
    type_counts = np.zeros(shape=type_shape, dtype=int)

    bamfile = pysam.Samfile(bam_file_name, 'rb')
    for CDS in simple_CDSs:
        reads = bamfile.fetch(CDS.seqname,
                              CDS.start - edge_overlap,
                              CDS.end + edge_overlap,
                             )
        for read in reads:
            if read.mapq != 50:
                # Non-unique mapping
                continue
            elif read.qlen < 25 or read.qlen > 31:
                continue
            elif sam.contains_indel_pysam(read):
                continue
            else:
                strand = '-' if read.is_reverse else '+'
                
                if strand != CDS.strand:
                    continue
                else:
                    alignment = sam.produce_alignment(read,
                                                      from_pysam=True,
                                                      fastq_type=fastq_type,
                                                     )

                    if strand == '+':
                        index_lookup = base_to_index
                    else:
                        index_lookup = base_to_complement_index
                   
                    for ref_char, read_char, qual, ref_pos, read_pos in alignment:
                        ref_index = index_lookup[ref_char]
                        read_index = index_lookup[read_char]
                        coords = (read.qlen - 25, read_pos, qual, ref_index, read_index)
                        type_counts[coords] += 1

    return type_counts

def recycling_ratios(rpf_positions_dict, simple_CDSs, genome):
    edge_overlap = 50
        
    ratio_lists = {amino_acid: defaultdict(list) for amino_acid in recycling.back_table}

    for gene in simple_CDSs:
        gene_name = gtf.parse_attribute(gene.attribute)['protein_id']
        codon_counts = get_codon_counts(rpf_positions_dict[gene_name], stringent=False)
        aa_locations = recycling.get_amino_acid_locations(gene, genome)
        if aa_locations == None:
            # MT or internal stop codon
            continue
        for amino_acid in aa_locations:
            location_pairs = izip(aa_locations[amino_acid], aa_locations[amino_acid][1:])
            for (first_position, first_codon), (second_position, second_codon) in location_pairs:
                first_count = codon_counts[first_position]
                second_count = codon_counts[second_position]
                ratio = (first_count, second_count)
                ratio_lists[amino_acid][first_codon, second_codon].append(ratio)
        
    return ratio_lists

def make_codon_counts_file(genes_dict, codon_counts_fn):
    with open(codon_counts_fn, 'w') as codon_counts_fh:
        for gene_name in sorted(genes_dict):
            counts = get_codon_counts(genes_dict[gene_name], stringent=False)
            counts_string = '\t'.join(str(count) for count in counts)
            line = '{0}\t{1}\n'.format(gene_name, counts_string)
            codon_counts_fh.write(line)

def make_expression_file(genes_dict, expression_fn, kind='RPKM'):
    if kind == 'RPKM':
        expression = get_RPKMs(genes_dict)
    elif kind == 'TPM':
        expression = get_TPMs(genes_dict)
    with open(expression_fn, 'w') as expression_fh:
        for gene_name in sorted(expression):
            line = '{0}\t{1:0.2f}\n'.format(gene_name, expression[gene_name])
            expression_fh.write(line)

def read_RPKMs_file(RPKMs_fn):
    def line_to_gene(line):
        name, value = line.strip().split()
        value = float(value)
        return name, value

    return dict(line_to_gene(line) for line in open(RPKMs_fn))

def read_bartel_file(bartel_file):
    RPF_dict = {}
    mRNA_dict = {}
    fh = open(bartel_file)
    # Skip column labels line
    fh.readline()
    for line in fh:
        name, mRNA_value, RPF_value = line.strip().split()
        mRNA_value = float(mRNA_value)
        RPF_value = float(RPF_value)
        mRNA_dict[name] = mRNA_value
        RPF_dict[name] = RPF_value

    return RPF_dict, mRNA_dict

def plot_metagene_averaged():
    # Generators that yields arrays of counts
    def counts_from_rpf_postiions_fn(rpf_positions_fn):
        rpf_positions_dict = Serialize.read_file(rpf_positions_fn, 'rpf_positions')
        for gene_name in rpf_positions_dict:
            counts = get_codon_counts(rpf_positions_dict[gene_name],
                                      stringent=False,
                                     )
            yield counts

    def counts_from_premal_fn(premal_fn):
        counts_dict = read_premal_file(premal_fn)
        for gene_name in counts_dict:
            counts = counts_dict[gene_name]
            yield counts
        
    rpf_experiments = [
        #('Geranshenko1', '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep1_foot/results/Initial_rep1_foot_rpf_positions.txt'),
        #('Geranshenko2', '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep2_foot/results/Initial_rep2_foot_rpf_positions.txt'),
        #('Geranshenko_mrna', '/home/jah/projects/arlen/experiments/gerashchenko_pnas/5min_rep1_mRNA/results/5min_rep1_mRNA_rpf_positions.txt'),

        ('Ingolia1', '/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-1/results/Footprints-rich-1_rpf_positions.txt'),
        #('Ingolia2', '/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-2/results/Footprints-rich-2_rpf_positions.txt'),
        #('Ingolia1_mRNA', '/home/jah/projects/arlen/experiments/ingolia_science/mRNA-rich-1/results/mRNA-rich-1_rpf_positions.txt'),

        ('Brar1' ,'/home/jah/projects/arlen/experiments/brar_science/s_tA-fp_100211_l3_sequence/results/s_tA-fp_100211_l3_sequence_rpf_positions.txt'),
        ('Brar2' ,'/home/jah/projects/arlen/experiments/brar_science/s_gb15exp_veg_-fp_100219_l4_sequence/results/s_gb15exp_veg_-fp_100219_l4_sequence_rpf_positions.txt'),
        ('Brar3' ,'/home/jah/projects/arlen/experiments/brar_science/s_14201exp_veg_-fp_100219_l6_sequence/results/s_14201exp_veg_-fp_100219_l6_sequence_rpf_positions.txt'),
        ('Brar4' ,'/home/jah/projects/arlen/experiments/brar_science/s_t1-fp_090807_l4/results/s_t1-fp_090807_l4_rpf_positions.txt'),

        #('suppressed', '/home/jah/projects/arlen/experiments/belgium_8_6_13/Suppressed_R98S_cDNA_sample/results/Suppressed_R98S_cDNA_sample_rpf_positions.txt'),
        #('R98S', '/home/jah/projects/arlen/experiments/belgium_8_6_13/R98S_cDNA_sample/results/R98S_cDNA_sample_rpf_positions.txt'),
        #('WT',  '/home/jah/projects/arlen/experiments/belgium_8_6_13/WT_cDNA_sample/results/WT_cDNA_sample_rpf_positions.txt'),
    ]

    premal_experiments = [
        ('Bartel', '/home/jah/projects/arlen/experiments/plotkin/genePosReads.txt'),
    ]

    all_experiments = [(name, counts_from_rpf_postiions_fn(fn)) for name, fn in rpf_experiments] + \
                      [(name, counts_from_premal_fn(fn)) for name, fn in premal_experiments]

    fig, ax = plt.subplots(figsize=(12, 12))

    plot_up_to = 500

    for name, counts_generator in all_experiments:
        sum_of_normalized = np.zeros(10000)
        long_enough_genes = np.zeros(10000)

        for counts in counts_generator:
            if counts.sum() < 64:
                continue

            num_codons = len(counts)
            density = counts.sum() / float(num_codons)
            normalized = counts / float(density)
            sum_of_normalized[:num_codons] += normalized
            long_enough_genes[:num_codons] += np.ones(num_codons)
        
        mean_densities = sum_of_normalized / long_enough_genes 

        ax.plot(mean_densities[:plot_up_to], '.-', label=name)
    
    leg = ax.legend(loc='upper right', fancybox=True)
    leg.get_frame().set_alpha(0.5)

    ax.set_xlabel('Position (codons)')
    ax.set_ylabel('Normalized mean reads')
    
    ax.plot(np.ones(plot_up_to), color='black', alpha=0.5)
    ax.set_ylim(0, 3)

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ax.set_aspect((xmax - xmin) / (ymax - ymin))

def plot_metagene_from_codon_counts():
    # Generators that yields arrays of counts
    def counts_from_codon_counts_fn(codon_counts_fn, report_P_site):
        counts_dict = read_codon_counts_file(codon_counts_fn)
        for gene_name in counts_dict:
            counts = counts_dict[gene_name]
            yield counts
        
    experiments = [
        ('Bartel', '/home/jah/projects/arlen/experiments/plotkin/genePosReads.txt', True),
        ('Ingolia', '/home/jah/projects/arlen/results/Ingolia_RPF_codon_counts.txt', False),
        ('Gerashchenko', '/home/jah/projects/arlen/results/Gerashchenko_RPF_codon_counts.txt', False),
        ('Brar', '/home/jah/projects/arlen/results/Brar_RPF_codon_counts.txt', False),
        ('Hussmann_WT', '/home/jah/projects/arlen/results/UT_WT_RPF_codon_counts.txt', False),
        ('Zinshteyn_1', '/home/jah/projects/arlen/results/Zinshteyn_1_RPF_codon_counts.txt', False),
        ('Zinshteyn_1_mRNA', '/home/jah/projects/arlen/results/Zinshteyn_1_mRNA_codon_counts.txt', False),
        #('Zinshteyn_2', '/home/jah/projects/arlen/results/Zinshteyn_2_RPF_codon_counts.txt', False),
    ]

    all_experiments = [(name, counts_from_codon_counts_fn(fn, reports_P_site))
                       for name, fn, reports_P_site in experiments]

    fig, ax = plt.subplots(figsize=(12, 12))

    plot_up_to = 500

    for name, counts_generator in all_experiments:
        sum_of_normalized = np.zeros(10000)
        long_enough_genes = np.zeros(10000)

        for counts in counts_generator:
            if counts.sum() < 64:
                continue

            num_codons = len(counts)
            density = counts.sum() / float(num_codons)
            normalized = counts / float(density)
            sum_of_normalized[:num_codons] += normalized
            long_enough_genes[:num_codons] += np.ones(num_codons)
        
        mean_densities = sum_of_normalized / long_enough_genes 

        ax.plot(mean_densities[:plot_up_to], '.-', label=name)
    
    leg = ax.legend(loc='upper right', fancybox=True)
    leg.get_frame().set_alpha(0.5)

    ax.set_xlabel('Position (codons)')
    ax.set_ylabel('Normalized mean reads')
    
    ax.plot(np.ones(plot_up_to), color='black', alpha=0.5)
    ax.set_ylim(0, 3)

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ax.set_aspect((xmax - xmin) / (ymax - ymin))

def plot_metagene_unaveraged(from_end=False):
    # Generators that yields arrays of counts
    def counts_from_rpf_postiions_fn(rpf_positions_fn):
        rpf_positions_dict = Serialize.read_file(rpf_positions_fn, 'rpf_positions')
        for gene_name in rpf_positions_dict:
            counts = get_codon_counts(rpf_positions_dict[gene_name],
                                      stringent=False,
                                     )
            yield counts

    def counts_from_premal_fn(premal_fn):
        counts_dict = read_premal_file(premal_fn)
        for gene_name in counts_dict:
            counts = counts_dict[gene_name]
            yield counts
        
    rpf_experiments = [
        #('Geranshenko1', '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep1_foot/results/Initial_rep1_foot_rpf_positions.txt'),
        #('Geranshenko2', '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep2_foot/results/Initial_rep2_foot_rpf_positions.txt'),
        #('Geranshenko_mrna', '/home/jah/projects/arlen/experiments/gerashchenko_pnas/5min_rep1_mRNA/results/5min_rep1_mRNA_rpf_positions.txt'),

        #('Ingolia1', '/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-1/results/Footprints-rich-1_rpf_positions.txt'),
        #('Ingolia2', '/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-2/results/Footprints-rich-2_rpf_positions.txt'),
        #('Ingolia1_mRNA', '/home/jah/projects/arlen/experiments/ingolia_science/mRNA-rich-1/results/mRNA-rich-1_rpf_positions.txt'),

        #('Brar1' ,'/home/jah/projects/arlen/experiments/brar_science/s_tA-fp_100211_l3_sequence/results/s_tA-fp_100211_l3_sequence_rpf_positions.txt'),
        #('Brar2' ,'/home/jah/projects/arlen/experiments/brar_science/s_gb15exp_veg_-fp_100219_l4_sequence/results/s_gb15exp_veg_-fp_100219_l4_sequence_rpf_positions.txt'),
        #('Brar3' ,'/home/jah/projects/arlen/experiments/brar_science/s_14201exp_veg_-fp_100219_l6_sequence/results/s_14201exp_veg_-fp_100219_l6_sequence_rpf_positions.txt'),
        #('Brar4' ,'/home/jah/projects/arlen/experiments/brar_science/s_t1-fp_090807_l4/results/s_t1-fp_090807_l4_rpf_positions.txt'),
        
        ('Nagalakshmi_RH_ori', '/home/jah/projects/arlen/experiments/nagalakshmi_science/RH_ori/results/RH_ori_rpf_positions.txt'),
        ('Nagalakshmi_RH_bio', '/home/jah/projects/arlen/experiments/nagalakshmi_science/RH_ori/results/RH_ori_rpf_positions.txt'),
        ('Nagalakshmi_dT_ori', '/home/jah/projects/arlen/experiments/nagalakshmi_science/dT_ori/results/dT_ori_rpf_positions.txt'),
        ('Nagalakshmi_dT_bio', '/home/jah/projects/arlen/experiments/nagalakshmi_science/dT_bio/results/dT_bio_rpf_positions.txt'),

        #('suppressed', '/home/jah/projects/arlen/experiments/belgium_8_6_13/Suppressed_R98S_cDNA_sample/results/Suppressed_R98S_cDNA_sample_rpf_positions.txt'),
        #('R98S', '/home/jah/projects/arlen/experiments/belgium_8_6_13/R98S_cDNA_sample/results/R98S_cDNA_sample_rpf_positions.txt'),
        #('WT',  '/home/jah/projects/arlen/experiments/belgium_8_6_13/WT_cDNA_sample/results/WT_cDNA_sample_rpf_positions.txt'),
    ]

    premal_experiments = [
        #('Bartel', '/home/jah/projects/arlen/experiments/plotkin/genePosReads.txt'),
    ]

    all_experiments = [(name, counts_from_rpf_postiions_fn(fn)) for name, fn in rpf_experiments] + \
                      [(name, counts_from_premal_fn(fn)) for name, fn in premal_experiments]

    plot_to = 500
    fig_cumulative, ax_cumulative = plt.subplots()

    if from_end:
        xs = np.arange(0, -plot_to, -1)
    else:
        xs = np.arange(plot_to)

    for name, counts_generator in all_experiments:
        print name

        expected_counts = np.zeros(5000)
        actual_counts = np.zeros(5000)

        for counts in counts_generator:
            num_codons = len(counts)
            r_g = counts.sum()
            uniform_counts = np.ones(num_codons) * r_g / num_codons
            
            actual_counts[:num_codons] += counts
            expected_counts[:num_codons] += uniform_counts

        normalized_cumulative = np.cumsum(actual_counts - expected_counts) / actual_counts.sum()

        #ax_cumulative.plot(xs, normalized_cumulative[:plot_to], label=name)
        ax_cumulative.plot(xs, ((actual_counts - expected_counts) / expected_counts)[:plot_to], label=name)

    ax_cumulative.plot(xs, np.zeros(plot_to), 'k--')
    ax_cumulative.legend()
    if from_end:
        xlabel = 'Codon position relative to end'
    else:
        xlabel = 'Codon position relative to start'
    ax_cumulative.set_xlabel(xlabel)
    ax_cumulative.set_ylabel('Cumulative sum of (actual - expected)')


#if __name__ == '__main__':
#    genome_index = mapping_tools.get_genome_index('/home/jah/projects/arlen/data/organisms/saccharomyces_cerevisiae/genome', explicit_path=True)
#
#    #clean_bam_fn = '/home/jah/projects/arlen/experiments/belgium_8_6_13/WT_cDNA_sample/results/WT_cDNA_sample_clean.bam'
#    #clean_bam_fn = '/home/jah/projects/arlen/experiments/belgium_8_6_13/R98S_cDNA_sample/results/R98S_cDNA_sample_clean.bam'
#    #clean_bam_fn = '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep1_foot/results/Initial_rep1_foot_clean.bam'
#    clean_bam_fn = '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep2_foot/results/Initial_rep2_foot_clean.bam'
#    
#    gtf_fn = '/home/jah/projects/arlen/data/organisms/saccharomyces_cerevisiae/transcriptome/genes.gtf'
#    
#    gene_name = 'YOR239W'
#    #gene_name = 'YPL052W'
#    #gene_name = 'YAL053W'
#    
#    extent = gtf.get_extent_by_name(gtf_fn, gene_name)
#    gene_length = extent[3] - extent[2] + 1
#    position_counts, expression_counts = get_extent_positions(clean_bam_fn, extent, genome_index)
#    position_ambiguity = determine_ambiguity_of_positions(extent, genome_index)
#    codon_numbers, frames, so_far, remaining = plot_frameshifts(position_counts,
#                                                                position_ambiguity,
#                                                                50,
#                                                                gene_name,
#                                                                gene_length,
#                                                               )

#if __name__ == '__main__':
#    density_length_correlation()
