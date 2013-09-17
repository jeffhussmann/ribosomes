import mapping
import gtf
import pysam
import subprocess
import sys
import fastq
import mutations
import numpy as np
import matplotlib.pyplot as plt
import itertools
from collections import Counter, defaultdict

fraction_resolution = 10

def filter_contaminant(contaminant_name,
                       trimmed_reads_fn,
                       no_contaminant_reads_fn,
                       sam_fn,
                       log_fn,
                       no_contaminant_lengths_fn,
                      ):
    ''' Maps reads in trimmed_reads_fn to the contaminant_name genome and
        records reads that don't map in no_contaminant_reads_fn. Records
        distribution of lengths of reads that don't map in
        no_contaminant_lengths_fn.
    '''
    mapping.map_bowtie2(trimmed_reads_fn,
                        contaminant_name,
                        sam_fn,
                        error_file_name=log_fn,
                        unaligned_reads_file_name=no_contaminant_reads_fn,
                        threads=8,
                        report_timing=True,
                       )
        
    no_rRNA_reads = fastq.reads(no_contaminant_reads_fn)
    lengths = Counter(len(read.seq) for read in no_rRNA_reads)
    lengths = mutations.counts_to_array(lengths)
    np.savetxt(no_contaminant_lengths_fn, lengths, fmt='%d')

def further_filtering(input_bam_fn,
                      gtf_fn,
                      filtered_bam_fn,
                      tRNA_lengths_fn,
                      rRNA_lengths_fn,
                      read_length,
                     ):
    ''' Removes any mappings to tRNA or rRNA genes in input_bam_fn and records
        all other primary mappings to filtered_bam_fn.
    '''
    lengths = {'tRNA': np.zeros(read_length + 1, int),
               'rRNA': np.zeros(read_length + 1, int),
              }

    bamfile = pysam.Samfile(input_bam_fn, 'rb')
    primary_reads = (read for read in bamfile if not read.is_secondary)
    
    contaminant_genes = iter(gtf.get_contaminant_genes(gtf_fn))
    gene = contaminant_genes.next()

    with pysam.Samfile(filtered_bam_fn, 'wb', header=bamfile.header) as filtered_bam_fh:
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
                filtered_bam_fh.write(read)

    pysam.index(filtered_bam_fn)

    np.savetxt(tRNA_lengths_fn, lengths['tRNA'], fmt='%d')
    np.savetxt(rRNA_lengths_fn, lengths['rRNA'], fmt='%d')

def map_tophat(file_names):
    tophat_command = ['tophat2',
                      '--GTF', file_names['gtf'],
                      '--no-novel-juncs',
                      '--num-threads', '8',
                      '--output-dir', file_names['tophat'],
                      file_names['bowtie2_index'],
                      file_names['no_rRNA_reads'],
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
                from_end = read.pos - CDS.end
            elif strand == '-' and CDS.strand == '-':
                from_start = CDS.end - (read.aend - 1)
                from_end = CDS.start - (read.aend - 1)
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

            keys = [right_strand, wrong_strand,
                    (28, 0), (28, 1), (28, 2),
                    (29, 0), (29, 1), (29, 2),
                    (30, 0), (30, 1), (30, 2),
                   ]
            info = '\t'.join(str(counts[key]) for key in keys)
            summary_fh.write('{0}\t{1}\n'.format(name, info))

    return from_starts, from_ends, fractions

def plot_RPF(from_starts_fns, from_ends_fns, names, read_lengths):
    from_starts_list = [np.loadtxt(fn) for fn in from_starts_fns]
    from_ends_list = [np.loadtxt(fn) for fn in from_ends_fns]

    fig, (start_ax, end_ax) = plt.subplots(2, 1, figsize=(12, 16))

    for from_starts, from_ends, name, read_length in zip(from_starts_list, from_ends_list, names, read_lengths):
        starts_xs = np.arange(-read_length, 10000 + 1)
        ends_xs = np.arange(-10000, read_length + 1)
        for length in range(28, 30):
            starts = np.true_divide(from_starts[length], from_starts[length].sum())
            ends = np.true_divide(from_ends[length], from_starts[length].sum())
            start_ax.plot(starts_xs, starts, '.-', label=name + '_{0}'.format(length))
            end_ax.plot(ends_xs, ends, '.-', label=name + '_{0}'.format(length))

    start_ax.set_xlim(-20, 20)
    end_ax.set_xlim(-30, 10)

    start_ax.legend()
    end_ax.legend()

def plot_mRNA(from_starts_fns, from_ends_fns, fractions_fns, names, read_lengths, fragment_lengths):
    from_starts_list = [np.loadtxt(fn) for fn in from_starts_fns]
    from_ends_list = [np.loadtxt(fn) for fn in from_ends_fns]
    fractions_list = [np.loadtxt(fn) for fn in fractions_fns]

    fig, (start_ax, end_ax, fraction_ax) = plt.subplots(3, 1, figsize=(12, 24))

    for from_starts, from_ends, fractions, name, read_length in zip(from_starts_list, from_ends_list, fractions_list, names, read_lengths):
        start_xs = np.arange(-read_length, 10000 + 1)
        end_xs = np.arange(-10000, read_length + 1)
        fraction_xs = np.linspace(0.5 / fraction_resolution,
                                  1 - 0.5 / fraction_resolution,
                                  fraction_resolution,
                                 )
        
        style = '.'

        starts = np.true_divide(from_starts[fragment_lengths].sum(axis=0), 1)
        start_ax.plot(start_xs, starts, style, label=name)

        ends = np.true_divide(from_ends[fragment_lengths].sum(axis=0), 1)
        end_ax.plot(end_xs, ends, style, label=name)
        
        fractions = np.true_divide(fractions[fragment_lengths].sum(axis=0),
                                   fractions.sum(),
                                  )
        fraction_ax.plot(fraction_xs, fractions, 'o-', label=name)

    start_ax.set_xlim(-max(read_lengths), 10000)
    start_ax.legend()
    
    end_ax.set_xlim(-10000, max(read_lengths))
    end_ax.legend()
    
    fraction_ax.set_xlim(0, 1)
    fraction_ax.set_ylim(ymin=0)
    end_ax.legend()

def plot_frames(summary_fns, names):
    frames = {}
    for summary_fn, name in zip(summary_fns, names):
        counts = np.loadtxt(summary_fn, usecols=range(3, 12))
        frames[name, 28] = counts[:, 0:3].sum(axis=0)
        frames[name, 29] = counts[:, 3:6].sum(axis=0)
        frames[name, 30] = counts[:, 6:9].sum(axis=0)

    fig, (ax_28, ax_29, ax_30) = plt.subplots(3, 1, figsize=(12, 24))
    for ax, length in zip([ax_28, ax_29, ax_30], [28, 29, 30]):
        rates_list = [frames[name, length] / frames[name, length].sum() for name in names]
        tick_labels = ['0', '1', '2']

        grouped_bar(rates_list, names, tick_labels, ax)
        ax.set_title(str(length))

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
