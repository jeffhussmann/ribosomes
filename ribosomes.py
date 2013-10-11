import mapping
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
import itertools
from collections import Counter, defaultdict

fraction_resolution = 10

def pre_filter(contaminant_index,
               trimmed_reads_fn,
               filtered_reads_fn,
               sam_fn,
               bam_fn,
               filtered_lengths_fn,
              ):
    ''' Maps reads in trimmed_reads_fn to contaminant_index. Records reads that
        don't map in filtered_reads_fn. Records distribution of lengths of reads
        that don't map in filtered_lengths_fn.
    '''
    mapping.map_bowtie2(trimmed_reads_fn,
                        contaminant_index,
                        sam_fn,
                        unaligned_reads_file_name=filtered_reads_fn,
                        threads=8,
                        report_all=True,
                        explicit_path=True,
                       )
    sam.sam_to_sorted_indexed_bam(sam_fn, bam_fn)
        
    filtered_reads = fastq.reads(filtered_reads_fn)
    lengths = Counter(len(read.seq) for read in filtered_reads)
    lengths = mutations.counts_to_array(lengths)
    np.savetxt(filtered_lengths_fn, lengths, fmt='%d')

def filter(input_bam_fn,
           gtf_fn,
           clean_bam_fn,
           read_length,
          ):
    ''' Removes any mappings to tRNA or rRNA genes in input_bam_fn and records
        all other primary mappings to clean_bam_fn.
    '''
    lengths = {'tRNA': np.zeros(read_length + 1, int),
               'rRNA': np.zeros(read_length + 1, int),
              }

    bamfile = pysam.Samfile(input_bam_fn, 'rb')
    primary_reads = (read for read in bamfile if not read.is_secondary)
    
    contaminant_genes = iter(gtf.get_contaminant_genes(gtf_fn))
    gene = contaminant_genes.next()

    with pysam.Samfile(clean_bam_fn, 'wb', header=bamfile.header) as clean_bam_fh:
        for read in primary_reads:
            seqname = bamfile.getrname(read.tid)
            while (gene.seqname, gene.end) < (seqname, read.pos):
                try:
                    gene = contaminant_genes.next()
                except StopIteration:
                    break
            if read.overlap(gene.start, gene.end) > 0:
                lengths[gene.source][read.qlen] += 1
            else:
                clean_bam_fh.write(read)

    pysam.index(clean_bam_fn)

    return lengths['tRNA'], lengths['rRNA']

def filter_fetch(input_bam_fn,
                 gtf_fn,
                 clean_bam_fn,
                 dirty_bam_fn,
                ):
    ''' Removes any m'''
    pysam.index(input_bam_fn)
    dirty_mappings = []
    dirty_names = set()
    contaminant_genes = gtf.get_contaminant_genes(gtf_fn)

    gene_counts = Counter()

    bamfile = pysam.Samfile(input_bam_fn, 'rb')
    with pysam.Samfile(dirty_bam_fn, 'wb', header=bamfile.header) as dirty_bam_fh:
        for gene in contaminant_genes:
            gene_name = gtf.parse_attribute(gene.attribute)['gene_name']
            overlapping_reads = bamfile.fetch(gene.seqname, gene.start, gene.end)
            for read in overlapping_reads:
                dirty_names.add(read.qname)
                if not read.is_secondary:
                    dirty_bam_fh.write(read)
                    gene_counts[gene_name] += 1
    bamfile.close()
    
    bamfile = pysam.Samfile(input_bam_fn, 'rb')
    with pysam.Samfile(clean_bam_fn, 'wb', header=bamfile.header) as clean_bam_fh:
        for read in bamfile:
            if read.qname not in dirty_names and not read.is_secondary:
                clean_bam_fh.write(read)

    for thing in gene_counts.most_common():
        print thing

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
    subprocess.check_call(tophat_command, stdout=sys.stdin, stderr=sys.stderr)

def get_aggregate_positions(gtf_fn, filtered_bam_fn):
    max_length = 50 
    max_position = 10000
    edge_overlap = max_length

    # Dimension 2 goes from -max_length to max_position for from_starts
    # and -max_position to max_length for from_ends.
    shape = (max_length + 1, max_length + max_position + 1)
    from_starts = np.zeros(shape, int)
    from_ends = np.zeros(shape, int)

    bamfile = pysam.Samfile(filtered_bam_fn, 'rb')

    counts = defaultdict(lambda : {'-': 0, '+': 0})

    simple_CDSs = iter(gtf.get_simple_CDSs(gtf_fn))
    gene = simple_CDSs.next()
    
    for read in bamfile:
        seqname = bamfile.getrname(read.tid)
        strand = '-' if read.is_reverse else '+'
        
        while (gene.seqname, gene.end + edge_overlap) < (seqname, read.pos):
            try:
                gene = simple_CDSs.next()
            except StopIteration:
                break

        if gene.strand == '+':
            overlaps = read.overlap(gene.start, gene.end + edge_overlap) > 0
        elif gene.strand == '-':
            overlaps = read.overlap(gene.start - edge_overlap, gene.end) > 0

        if overlaps:
            counts[gene][strand] += 1
            if strand == '+' and gene.strand == '+':
                from_start = read.pos - gene.start
                if from_start < max_position:
                    from_starts[read.qlen, max_length + from_start] += 1

                from_end = read.pos - gene.end
                if from_end > -max_position:
                    from_ends[read.qlen, max_position + from_end] += 1

            elif strand == '-' and gene.strand == '-':
                from_start = gene.end - (read.aend - 1)
                if from_start < max_position:
                    from_starts[read.qlen, max_length + from_start] += 1

                from_end = gene.start - (read.aend - 1)
                if from_end > -max_position:
                    from_ends[read.qlen, max_position + from_end] += 1
    
    return from_starts, from_ends

def get_aggregate_positions_fetch(gtf_fn,
                                  filtered_bam_fn,
                                  summary_fn,
                                  read_length,
                                 ):
    max_length = read_length 
    max_position = 10000
    edge_overlap = max_length

    # Dimension 2 goes from -max_length to max_position for from_starts
    # and -max_position to max_length for from_ends.
    shape = (max_length + 1, edge_overlap + max_position + 1)
    from_starts = np.zeros(shape, int)
    from_ends = np.zeros(shape, int)

    fraction_shape = (max_length + 1, fraction_resolution)
    fractions = np.zeros(fraction_shape, int)

    bamfile = pysam.Samfile(filtered_bam_fn, 'rb')

    read_counts = defaultdict(lambda : defaultdict(int))

    simple_CDSs = gtf.get_simple_CDSs(gtf_fn)
   
    for CDS in simple_CDSs:
        max_from_start = CDS.end - CDS.start
        fraction_factor = float(fraction_resolution) / max_from_start

        reads = bamfile.fetch(CDS.seqname, CDS.start - edge_overlap, CDS.end + edge_overlap)
        primary_reads = (read for read in reads if not read.is_secondary)

        for read in primary_reads:
            strand = '-' if read.is_reverse else '+'
        
            if 28 <= len(read.seq) <= 30:
                read_counts[CDS][strand] += 1
            
            if strand == '+' and CDS.strand == '+':
                from_start = read.pos - CDS.start
                # CDS.end is the last base before the stop codon
                from_end = read.pos - (CDS.end + 1)
            elif strand == '-' and CDS.strand == '-':
                from_start = CDS.end - (read.aend - 1)
                # CDS.start is the first base after the stop codon
                from_end = (CDS.start - 1) - (read.aend - 1)
            else:
                continue

            if -edge_overlap <= from_start < max_position:
                from_starts[read.qlen, edge_overlap + from_start] += 1
                if 28 <= len(read.seq) <= 30:
                    read_counts[CDS][len(read.seq), from_start % 3] += 1
            
            if from_end > -max_position:
                from_ends[read.qlen, max_position + from_end] += 1
            
            if from_start >= 0:
                if from_start == max_from_start:
                    fraction = fraction_resolution - 1  
                else:
                    fraction = int(fraction_factor * from_start)
                if fraction < fraction_resolution:
                    fractions[read.qlen, fraction] += 1

    with open(summary_fn, 'w') as summary_fh:
        for CDS, counts in read_counts.items():
            name = gtf.parse_attribute(CDS.attribute)['gene_name']
            right_strand = CDS.strand
            wrong_strand = '-' if CDS.strand == '+' else '+'
            length = CDS.end - CDS.start + 1
            keys = [right_strand, wrong_strand,
                    (28, 0), (28, 1), (28, 2),
                    (29, 0), (29, 1), (29, 2),
                    (30, 0), (30, 1), (30, 2),
                   ]
            info = '\t'.join(str(counts[key]) for key in keys)
            summary_fh.write('{0}\t{1}\t{2}\n'.format(name, length, info))

    return from_starts, from_ends, fractions

def plot_RPF(from_starts_fns, from_ends_fns, names, read_lengths):
    from_starts_list = [np.loadtxt(fn) for fn in from_starts_fns]
    from_ends_list = [np.loadtxt(fn) for fn in from_ends_fns]

    fig, (start_ax, end_ax) = plt.subplots(2, 1, figsize=(12, 16))

    for from_starts, from_ends, name, read_length in zip(from_starts_list, from_ends_list, names, read_lengths):
        starts_xs = np.arange(-read_length, 10000 + 1)
        ends_xs = np.arange(-10000, read_length + 1)
        for length in [22, 23, 24, 25, 28]:
            starts = np.true_divide(from_starts[length], from_starts[length].sum())
            ends = np.true_divide(from_ends[length], from_starts[length].sum())
            start_ax.plot(starts_xs, starts, '.-', label=name + '_{0}'.format(length))
            end_ax.plot(ends_xs, ends, '.-', label=name + '_{0}'.format(length))

    start_ax.set_xlim(-20, 20)
    start_ax.set_xlabel('Position of read relative to start of CDS')
    start_ax.set_ylabel('Fraction of uniquely mapped reads')
    #start_ax.set_ylim(ymax=0.008)
    end_ax.set_xlim(-35, 5)
    end_ax.set_xlabel('Position of read relative to stop codon')
    end_ax.set_ylabel('Fraction of uniquely mapped reads')

    leg = start_ax.legend(loc='upper right', fancybox=True)
    leg.get_frame().set_alpha(0.5)
    leg = end_ax.legend(loc='upper right', fancybox=True)
    leg.get_frame().set_alpha(0.5)

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

    colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])

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

def plot_clean_lengths(clean_length_fns, names):
    fig, ax = plt.subplots(figsize=(12, 8))
    colors = itertools.cycle(['g', 'b', 'r', 'c', 'm', 'y', 'k'])

    for clean_fn, name, color in zip(clean_length_fns, names, colors):
        clean_lengths = np.loadtxt(clean_fn)
        normalized = np.true_divide(clean_lengths,
                                    clean_lengths.sum(),
                                    #1,
                                   )
        ax.plot(normalized, '.-', label=name, color=color)
        #ax.axvspan(27.5, 28.5, color='green', alpha=0.1)
        ax.set_xlim(0, 50)
        ax.legend()
        ax.set_title('\'Clean\' fragment length distribution by sample')
        ax.set_xlabel('Length of original RNA fragment')
        ax.set_ylabel('Fraction of reads')

    leg = ax.legend(loc='upper right', fancybox=True)
    leg.get_frame().set_alpha(0.5)

def plot_trimmed_lengths(trimmed_length_fns, names):
    fig, ax = plt.subplots(figsize=(12, 8))

    for trimmed_fn, name in zip(trimmed_length_fns, names):
        trimmed_lengths = np.loadtxt(trimmed_fn)
        normalized = np.true_divide(trimmed_lengths,
                                    #trimmed_lengths.sum(),
                                    1,
                                   )
        ax.plot(normalized, '.-', label=name, color='black')

    leg = ax.legend(loc='upper left', fancybox=True)
    leg.get_frame().set_alpha(0.5)

    ax.set_title('Fragment length distribution')
    ax.set_ylabel('Number of reads')
    ax.set_xlabel('Length of original RNA fragment') 

def get_all_lengths(bam_fn):
    bamfile = pysam.Samfile(bam_fn, 'rb')
    qlens = Counter(ar.qlen for ar in bamfile)
    return mutations.counts_to_array(qlens)

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

def produce_rRNA_coverage(bam_fn, oligos_sam_fn, coverage_fn):
    bam_file = pysam.Samfile(bam_fn, 'rb')

    counts = {i: np.zeros(length, int) for i, length in enumerate(bam_file.lengths)}
    
    read_names = set()

    progress = mutations.progress_bar(bam_file.mapped + bam_file.unmapped)
    for i, aligned_read in itertools.izip(progress, bam_file):
        read_names.add(aligned_read.qname)
        tid = aligned_read.tid
        for position in aligned_read.positions:
            counts[tid][position] += 1

    with open(coverage_fn, 'w') as output_fh:
        output_fh.write('>Total reads\n')
        output_fh.write('{0}\n'.format(len(read_names)))
        for tid in counts:
            output_fh.write('>{0}\n'.format(bam_file.getrname(tid)))
            output_fh.write('{0}\n'.format(' '.join(str(c) for c in counts[tid])))

def load_rRNA_coverage(input_fn):
    input_fh = open(input_fn)
    line_pairs = itertools.izip(*[input_fh]*2)
    _, total_reads  = line_pairs.next()
    total_reads = int(total_reads.strip())
    counts = {}
    for name_line, counts_line in line_pairs:
        name = name_line.strip().lstrip('>')
        counts[name] = np.fromstring(counts_line, sep=' ', dtype=int)
    return total_reads, counts

def plot_rRNA_coverage(names, coverage_fns, oligos_sam_fn):
    oligos_sam_file = pysam.Samfile(oligos_sam_fn, 'r')
    
    total_reads = {}
    counts = {}
    for name, coverage_fn in zip(names, coverage_fns):
        total_reads[name], counts[name] = load_rRNA_coverage(coverage_fn)

    rnames = oligos_sam_file.references
    lengths = oligos_sam_file.lengths

    figs = {}
    axs = {}
    for rname, length in zip(rnames, lengths):
        figs[rname], axs[rname] = plt.subplots(figsize=(12, 8))
        axs[rname].set_title(rname)
        axs[rname].set_xlim(0, length)

        for name in names:
            normalized_counts = np.true_divide(counts[name][rname], total_reads[name])
            axs[rname].plot(normalized_counts, label=name)

        leg = axs[rname].legend(loc='upper right', fancybox=True)
        leg.get_frame().set_alpha(0.5)
    
    oligo_mappings = load_oligo_mappings(oligos_sam_fn)

    colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])
    for (oligo, mappings), color in itertools.izip(oligo_mappings.items(), colors):
        for rname, start, end in mappings:
            axs[rname].axvspan(start, end, color=color, alpha=0.2, linewidth=0)
            _, y_max = axs[rname].get_ylim()
            axs[rname].text(float(start + end) / 2,
                            y_max,
                            oligo, 
                            horizontalalignment='center',
                            verticalalignment='top',
                           )

def load_oligo_mappings(oligos_sam_fn):
    oligos_sam_file = pysam.Samfile(oligos_sam_fn, 'r')
    oligo_mappings = defaultdict(list)
    for aligned_read in oligos_sam_file:
        positions = aligned_read.positions
        rname = oligos_sam_file.getrname(aligned_read.tid)
        extent = (rname, min(positions), max(positions))
        oligo_mappings[aligned_read.qname].append(extent)
    return oligo_mappings

def get_oligo_hit_lengths(bam_fn, oligos_sam_fn):
    oligo_mappings = load_oligo_mappings(oligos_sam_fn)
    bam_file = pysam.Samfile(bam_fn, 'rb')

    lengths = {qname: np.zeros(51, int) for qname in oligo_mappings}
    for qname, extents in oligo_mappings.items():
        print qname
        for rname, start, end in extents:
            reads = bam_file.fetch(rname, start, end)
            for aligned_read in reads:
                if not aligned_read.is_secondary:
                    lengths[qname][aligned_read.qlen] += 1
    
    return lengths

def plot_oligo_hit_lengths(bam_fn, oligos_sam_fn, fig_fn, name):
    fig, ax = plt.subplots(figsize=(18, 12))
    lengths = get_oligo_hit_lengths(bam_fn, oligos_sam_fn)
    colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'BlueViolet', 'Gold'])
    qnames = sorted(lengths.keys())
    for qname, color in zip(qnames, colors):
        ax.plot(lengths[qname], 'o-', color=color, label=qname)
    leg = ax.legend(loc='upper right', fancybox=True)
    leg.get_frame().set_alpha(0.5)
    ax.set_title(name)
    plt.savefig(fig_fn)
    plt.close()
