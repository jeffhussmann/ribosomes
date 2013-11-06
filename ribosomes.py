import mapping
import fasta
import random
import os.path
import gtf
import pysam
import sam
import subprocess
import sys
import fastq
import mutations
import numpy as np
import matplotlib.pyplot as plt
from itertools import izip, cycle, chain
from collections import Counter, defaultdict

fraction_resolution = 10
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'BlueViolet', 'Gold']

def pre_filter(contaminant_index,
               trimmed_reads_fn,
               filtered_reads_fn,
               sam_fn,
               bam_fn,
              ):
    ''' Maps reads in trimmed_reads_fn to contaminant_index. Records reads that
        don't map in filtered_reads_fn. 
    '''
    mapping.map_bowtie2(trimmed_reads_fn,
                        contaminant_index,
                        sam_fn,
                        unaligned_reads_file_name=filtered_reads_fn,
                        threads=1,
                        report_all=True,
                        explicit_path=True,
                       )
    sam.make_sorted_indexed_bam(sam_fn, bam_fn)

def post_filter(input_bam_fn,
                gtf_fn,
                clean_bam_fn,
                more_rRNA_bam_fn,
                tRNA_bam_fn,
               ):
    ''' Removes any remaining mappings to tRNA or rRNA genes.
        Removes secondary mappings so that each read is represented at most
        once, but this DOES NOT mean that only unique mappings are retained.
    '''
    qnames_already_used = set()

    rRNA_genes = gtf.get_rRNA_genes(gtf_fn)
    tRNA_genes = gtf.get_tRNA_genes(gtf_fn)

    input_bam_file = pysam.Samfile(input_bam_fn, 'rb')
   
    # Find reads with any mappings that overlap rRNA or tRNA genes and write one
    # mapping for each such read to a contaminant bam file.
    for genes, bam_fn in [(rRNA_genes, more_rRNA_bam_fn),
                          (tRNA_genes, tRNA_bam_fn)]:
        contaminant_reads = defaultdict(list)
        for gene in genes:
            overlapping_mappings = input_bam_file.fetch(gene.seqname,
                                                        gene.start,
                                                        gene.end,
                                                       )
            for mapping in overlapping_mappings:
                contaminant_reads[read.qname].append(mapping)

        with pysam.Samfile(bam_fn, 'wb', header=input_bam_file.header) as bam_file:
            for qname in contaminant_reads:
                if qname not in qnames_already_used:
                    mapping = random.choice(contaminant_reads[qname])
                    bam_file.write(mapping)
                    qnames_already_used.add(qname)

    input_bam_file.close()
         
    # Create a new clean bam file consisting of the primary mapping of each
    # read that wasn't flagged as a contaminant.
    input_bam_file = pysam.Samfile(input_bam_fn, 'rb')
    with pysam.Samfile(clean_bam_fn, 'wb', header=input_bam_file.header) as clean_bam_file:
        for mapping in input_bam_file:
            if mapping.qname not in qnames_already_used and not mapping.is_secondary:
                clean_bam_file.write(mapping)

    pysam.index(clean_bam_fn)

def map_tophat(reads_file_name,
               bowtie2_index,
               gtf_file_name,
               transcriptome_index,
               tophat_dir,
              ):
    tophat_command = ['tophat2',
                      '--GTF', gtf_file_name,
                      '--no-novel-juncs',
                      '--num-threads', '8',
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

    # position_counts will hold, for each gene, for each position from
    # -edge_overlap to gene_length + edge_overlap
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

def frameshifts(rpf_positions, edge_overlap):
    rpf_positions = dict(zip(*rpf_positions))
    gene_name = 'YLR233C'
    gene_length = 2097
    positions = rpf_positions[gene_name, gene_length][0]
    start_at_codon = -4
    codon_starts = np.arange(edge_overlap + (start_at_codon * 3),
                             edge_overlap + gene_length,
                             3,
                            )
    print codon_starts
    print len(positions)
    codon_numbers = [start_at_codon + i for i in range(len(codon_starts))]
    frames = np.zeros((3, len(codon_starts)), int)
    for c, codon_start in enumerate(codon_starts):
        for frame in range(3):
            frames[frame, c] = positions[codon_start + frame]

    frames_so_far = frames.cumsum(axis=1)
    fraction_frames_so_far = np.true_divide(frames_so_far, frames_so_far.sum(axis=0))

    frames_remaining = np.fliplr(np.fliplr(frames).cumsum(axis=1))
    fraction_frames_remaining = np.true_divide(frames_remaining, frames_remaining.sum(axis=0))
    
    print frames
    return codon_numbers, fraction_frames_so_far, fraction_frames_remaining

def plot_RPF(from_starts_fns, from_ends_fns, names):
    edge_overlap = 50
    from_starts_list = [np.loadtxt(fn) for fn in from_starts_fns]
    from_ends_list = [np.loadtxt(fn) for fn in from_ends_fns]

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

    colors = cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])

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

def produce_rRNA_coverage(bam_file_names):
    bam_files = [pysam.Samfile(fn, 'rb') for fn in bam_file_names]
    all_reads = chain.from_iterable(bam_files)

    seq_lengths = bam_files[0].lengths
    counts = [np.zeros(length, int) for length in seq_lengths]
    names = [bam_files[0].getrname(tid) for tid in range(len(seq_lengths))]
    
    read_names = set()

    for aligned_read in all_reads:
        read_names.add(aligned_read.qname)
        tid = aligned_read.tid
        for position in aligned_read.positions:
            counts[tid][position] += 1

    total_reads = len(read_names)
    return total_reads, names, counts

def plot_rRNA_coverage(names, data_list, oligos_sam_fn, fig_fn_template):
    oligos_sam_file = pysam.Samfile(oligos_sam_fn, 'r')
    rnames = oligos_sam_file.references
    lengths = oligos_sam_file.lengths
    oligo_mappings = load_oligo_mappings(oligos_sam_fn)
    
    total_reads = {}
    counts = {}
    for name, (these_reads, _, these_counts) in zip(names, data_list):
        total_reads[name] = these_reads
        counts[name] = these_counts

    figs = {}
    axs = {}
    for i, (rname, length) in enumerate(zip(rnames, lengths)):
        figs[rname], axs[rname] = plt.subplots(figsize=(18, 12))
        axs[rname].set_title('rRNA identity - ' + rname)
        axs[rname].set_xlim(0, length)

        for name in names:
            normalized_counts = np.true_divide(counts[name][i], total_reads[name])
            axs[rname].plot(normalized_counts, label=name)

        leg = axs[rname].legend(loc='upper right', fancybox=True)
        leg.get_frame().set_alpha(0.5)
        axs[rname].set_xlabel('Position in rRNA')
        axs[rname].set_ylabel('Fraction of all reads mapping to position')

    
    for oligo, color in izip(oligo_mappings, colors):
        for rname, start, end in oligo_mappings[oligo]:
            axs[rname].axvspan(start, end, color=color, alpha=0.12, linewidth=0)
            _, y_max = axs[rname].get_ylim()
            axs[rname].text(float(start + end) / 2,
                            y_max,
                            oligo, 
                            horizontalalignment='center',
                            verticalalignment='top',
                           )
    for rname in rnames:
        figs[rname].savefig(fig_fn_template.format(rname))
        plt.close(figs[rname])

def load_oligo_mappings(oligos_sam_fn):
    oligos_sam_file = pysam.Samfile(oligos_sam_fn, 'r')
    oligo_mappings = defaultdict(list)
    for aligned_read in oligos_sam_file:
        positions = aligned_read.positions
        rname = oligos_sam_file.getrname(aligned_read.tid)
        extent = (rname, min(positions), max(positions))
        oligo_mappings[aligned_read.qname].append(extent)
    return oligo_mappings

def get_oligo_hit_lengths(bam_fn,
                          oligos_fasta_fn,
                          oligos_sam_fn,
                          max_read_length):
    oligo_mappings = load_oligo_mappings(oligos_sam_fn)
    bam_file = pysam.Samfile(bam_fn, 'rb')

    oligo_names = [read.name for read in fasta.reads(oligos_fasta_fn)]
    lengths = np.zeros((len(oligo_names), max_read_length + 1), int)

    for oligo_number, oligo_name in enumerate(oligo_names):
        for rname, start, end in oligo_mappings[oligo_name]:
            reads = bam_file.fetch(rname, start, end)
            for aligned_read in reads:
                if not aligned_read.is_secondary:
                    lengths[oligo_number][aligned_read.qlen] += 1
    
    return lengths

def plot_oligo_hit_lengths(oligos_fasta_fn, lengths, fig_fn):
    oligo_names = [read.name for read in fasta.reads(oligos_fasta_fn)]
    
    fig, ax = plt.subplots(figsize=(18, 12))
    for oligo_name, oligo_lengths, color in zip(oligo_names, lengths, colors):
        normalized_lengths = np.true_divide(oligo_lengths, oligo_lengths.sum())
        ax.plot(normalized_lengths, 'o-', color=color, label=oligo_name)
    
    leg = ax.legend(loc='upper right', fancybox=True)
    leg.get_frame().set_alpha(0.5)
    
    ax.set_xlim(0, lengths.shape[1] - 1)

    ax.set_xlabel('Length of original RNA fragment')
    ax.set_ylabel('Number of fragments')
    ax.set_title('Distribution of fragment lengths overlapping each oligo')
    
    fig.savefig(fig_fn)
    plt.close(fig)
