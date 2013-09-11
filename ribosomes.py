import mapping
import gtf
import pysam
import subprocess
import sys
import fastq
import mutations
import numpy as np
from collections import Counter, defaultdict

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
                     ):
    ''' Removes any mappings to tRNA or rRNA genes in input_bam_fn and records
        all other primary mappings to filtered_bam_fn.
    '''
    lengths = {'tRNA': np.zeros(51, int),
               'rRNA': np.zeros(51, int),
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
        
        while (gene.seqname, gene.end) < (seqname, read.pos):
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
