from adapters_cython import *
import fastq
import sam
from collections import Counter, namedtuple, defaultdict
import glob
import os.path
import mapping
import mutations
import numpy as np
import time
import pysam
import matplotlib.pyplot as plt
import subprocess
import sys

truseq_R2_rc_up_to_index = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'

def get_read_length(fn):
    return len(fastq.reads(fn).next().seq)

def characterize(reads_fn, distances_fn):
    read_length = get_read_length(reads_fn)
    d_array = characterize_adapters(fastq.reads(reads_fn),
                                    truseq_R2_rc_up_to_index,
                                    read_length,
                                   )
    np.savetxt(distances_fn, d_array, fmt='%d')

def trim_adapters(in_fn, out_fn, lengths_fn):
    lengths = np.zeros(get_read_length(in_fn) + 1, int)
    adapter_prefix = truseq_R2_rc_up_to_index[:10]
    
    with open(out_fn, 'w') as out_fh:
        for read in fastq.reads(in_fn):
            p = find_adapter_position(read.seq, truseq_R2_rc_up_to_index, 10, 3)
            lengths[p] += 1
            
            trimmed_slice = slice(None, p)
            seq = read.seq[trimmed_slice]
            qual = read.qual[trimmed_slice]
            
            if len(seq) > 5:
                out_fh.write(fastq.make_record(read.name, seq, qual))

    np.savetxt(lengths_fn, lengths, fmt='%d')

def filter_contaminant(contaminant, trimmed_fn, no_contaminant_fn, sam_fn, log_fn, lengths_fn):
    mapping.map_bowtie2(trimmed_fn,
                        contaminant,
                        sam_fn,
                        error_file_name=log_fn,
                        unaligned_reads_file_name=no_contaminant_fn,
                        threads=8,
                        report_timing=True,
                       )
    
    lengths = Counter(len(read.seq) for read in fastq.reads(no_contaminant_fn))
    lengths = mutations.counts_to_array(lengths)
    np.savetxt(lengths_fn, lengths, fmt='%d')

def make_file_names(sample):
    data_dir = '/home/jah/projects/arlen/data'
    results_dir = '/home/jah/projects/arlen/results'
    
    if not os.path.isdir('{0}/{1}/'.format(results_dir, sample)):
        os.makedirs('{0}/{1}/'.format(results_dir, sample))

    data_base = '{0}/{1}/{1}'.format(data_dir, sample)
    results_base = '{0}/{1}/{1}'.format(results_dir, sample)
    tophat_base = '{0}/{1}/tophat/'.format(results_dir, sample)

    file_names = {
        'results': results_dir,
        'data': data_dir,
        'tophat': tophat_base,
        'R1_reads': glob.glob(data_base + '*.fastq')[0],
        'gtf': data_dir + '/igenomes/Saccharomyces_cerevisiae/Ensembl/EF2/Annotation/Genes/genes.gtf'
    }

    tails = {
        'trimmed_reads':     '_trimmed.fastq',
        'trimmed_lengths':   '_trimmed_lengths.txt',
        'adapter_distances': '_adapter_distances.txt',
        'no_rRNA_reads':     '_no_rRNA.fastq',
        'rRNA_mapping_log':  '_rRNA_mapping_log.txt',
        'rRNA_sam':          '_rRNA.sam',
        'no_rRNA_lengths':   '_no_rRNA_lengths.txt',
        'filtered_bam':      '_filtered.bam',
        'tRNA_lengths':      '_tRNA_lengths.txt',
        'more_rRNA_lengths': '_more_rRNA_lengths.txt',
        'lengths_fig':       '_lengths.png',
    }

    tophat_tails = {
        'accepted_hits': 'accepted_hits.bam',
        'unmapped': 'unmapped.bam',
    }

    results = {name: results_base + tail for name, tail in tails.items()}
    tophat_results = {name: tophat_base + tail for name, tail in tophat_tails.items()}
    file_names.update(results)
    file_names.update(tophat_results)
    return file_names

def map_tophat(file_names):
    igenomes_dir = '/home/jah/projects/arlen/data/igenomes/Saccharomyces_cerevisiae/Ensembl/EF2'

    tophat_command = ['tophat2',
                      '--GTF', '{0}/Annotation/Genes/genes.gtf'.format(igenomes_dir),
                      '--no-novel-juncs',
                      '--num-threads', '8',
                      '--output-dir', file_names['tophat'],
                      '{0}/Sequence/Bowtie2Index/genome'.format(igenomes_dir),
                      file_names['no_rRNA_reads'],
                     ]
    subprocess.check_call(tophat_command, stdout=sys.stdin, stderr=sys.stderr)

gtf_fields = ['seqname',
              'source',
              'feature',
              'start',
              'end',
              'score',
              'strand',
              'frame',
              'attribute',
             ]
Gene = namedtuple('Gene', gtf_fields)

def parse_gtf_line(line):
    gene = Gene._make(line.strip().split('\t'))
    start = int(gene.start)
    end = int(gene.end)
    if gene.frame != '.':
        frame = int(gene.frame)
    else:
        frame = gene.frame
    gene = gene._replace(start=start, end=end, frame=frame)
    return gene

def get_all_genes(gtf_fn):
    all_genes = [parse_gtf_line(line) for line in open(gtf_fn)]
    return all_genes

def get_contaminant_genes(gtf_fn):
    all_genes = get_all_genes(gtf_fn)
    contaminant_genes = [gene for gene in all_genes if gene.source == 'rRNA' or gene.source == 'tRNA']
    return contaminant_genes

def get_easy_CDSs(gtf_fn):
    all_genes = get_all_genes(gtf_fn)
    CDSs = []
    for i, gene in enumerate(all_genes):
        if gene.source == 'protein_coding' and gene.feature == 'CDS':
            CDSs.append(gene)

    return CDSs

def further_filtering(input_bam_fn,
                      filtered_bam_fn,
                      tRNA_lengths_fn,
                      rRNA_lengths_fn,
                      gtf_fn,
                     ):
    ''' Removes any mappings to tRNA or rRNA genes. '''
    lengths = {'tRNA': np.zeros(51, int),
               'rRNA': np.zeros(51, int),
              }

    bamfile = pysam.Samfile(input_bam_fn, 'rb')
    primary_reads = (aligned_read for aligned_read in bamfile if not aligned_read.is_secondary)
    
    contaminant_genes = iter(get_contaminant_genes(gtf_fn))
    
    gene = contaminant_genes.next()

    with pysam.Samfile(filtered_bam_fn, 'wb', header=bamfile.header) as filtered_bam_fh:
        for aligned_read in primary_reads:
            seqname = bamfile.getrname(aligned_read.tid)
            while (gene.seqname, gene.end) < (seqname, aligned_read.pos):
                try:
                    gene = contaminant_genes.next()
                except StopIteration:
                    break
            if aligned_read.overlap(gene.start, gene.end) > 0:
                lengths[gene.source][aligned_read.qlen] += 1
            else:
                filtered_bam_fh.write(aligned_read)

    np.savetxt(tRNA_lengths_fn, lengths['tRNA'], fmt='%d')
    np.savetxt(rRNA_lengths_fn, lengths['rRNA'], fmt='%d')

def get_all_lengths(bam_fn):
    bamfile = pysam.Samfile(bam_fn, 'rb')
    qlens = Counter(ar.qlen for ar in bamfile)
    return mutations.counts_to_array(qlens)

def get_gaps(genes):
    gaps = Counter()
    for i in range(len(genes) - 1):
        if genes[i].seqname == genes[i + 1].seqname:
            gap = genes[i + 1].start - genes[i].end
            if gap == -1304:
                print genes[i]
                print genes[i + 1]
                print
            gaps[gap] += 1

    xs = np.arange(min(gaps), max(gaps) + 1)
    ys = np.asarray([gaps[x] for x in xs])

    return xs, ys

def get_gaps(genes):
    gaps = Counter()
    for i in range(len(genes) - 1):
        if genes[i].seqname == genes[i + 1].seqname:
            gap = genes[i + 1].start - genes[i].end
            if gap == -1304:
                print genes[i]
                print genes[i + 1]
                print
            gaps[gap] += 1

    xs = np.arange(min(gaps), max(gaps) + 1)
    ys = np.asarray([gaps[x] for x in xs])

    return xs, ys

def get_nonoverlapping(genes):
    overlapping = set()
    nonoverlapping = []
    for i, gene in enumerate(genes):
        j = i + 1
        while j < len(genes) and (genes[j].seqname, genes[j].start) <= (gene.seqname, gene.end):
            overlapping.add(gene)
            overlapping.add(genes[j])
            j += 1
        if gene not in overlapping:
            nonoverlapping.append(gene)

    return nonoverlapping

def plot_lengths(file_names):
    trimmed_lengths = np.loadtxt(file_names['trimmed_lengths'])
    no_rRNA_lengths = np.loadtxt(file_names['no_rRNA_lengths'])
    unmapped_lengths = get_all_lengths(file_names['unmapped'])
    tRNA_lengths = np.loadtxt(file_names['tRNA_lengths'])
    more_rRNA_lengths = np.loadtxt(file_names['more_rRNA_lengths'])

    rRNA_lengths = trimmed_lengths - no_rRNA_lengths + more_rRNA_lengths
    rpf_lengths = no_rRNA_lengths - more_rRNA_lengths - tRNA_lengths - unmapped_lengths

    plt.plot(rRNA_lengths, '.-', label='rRNA')
    plt.plot(tRNA_lengths, '.-', label='tRNA')
    plt.plot(rpf_lengths, '.-', label='RPF')
    plt.legend()
    plt.savefig(file_names['lengths_fig'])
    plt.close()

def count_reads(nonoverlapping_genes, filtered_bam_fn):
    positions = Counter()
    bamfile = pysam.Samfile(filtered_bam_fn, 'rb')

    counts = defaultdict(lambda : {'-': 0, '+': 0})

    nonoverlapping_genes = iter(nonoverlapping_genes)
    gene = nonoverlapping_genes.next()
    
    for aligned_read in bamfile:
        seqname = bamfile.getrname(aligned_read.tid)
        strand = '-' if aligned_read.is_reverse else '+'
        while (gene.seqname, gene.end) < (seqname, aligned_read.pos):
            try:
                gene = nonoverlapping_genes.next()
            except StopIteration:
                break
        if aligned_read.overlap(gene.start, gene.end) > 0:
            counts[gene][strand] += 1
            if strand == '+':
                position = aligned_read.pos - gene.start
                if position < -50:
                    print (gene.seqname, gene.start)
                    print (seqname, aligned_read.pos)
                positions[position] += 1
    
    ps = np.arange(min(positions), max(positions) + 1)
    cs = np.asarray([positions[p] for p in ps])

    return ps, cs
    
if __name__ == '__main__':
    samples = [
        #'R98S_cDNA_mRNA',
        #'R98S_cDNA_sample',
        #'Suppressed_R98S_cDNA_mRNA',
        'Suppressed_R98S_cDNA_sample',
        #'WT_cDNA_mRNA',
        #'WT_cDNA_sample',
    ]

    for sample in samples:
        print sample
        file_names = make_file_names(sample)

        start_time = time.time()
        #trim_adapters(file_names['R1_reads'],
        #              file_names['trimmed_reads'],
        #              file_names['trimmed_lengths'],
        #             )
        #filter_contaminant('yeast_rRNA',
        #                   file_names['trimmed_reads'],
        #                   file_names['no_rRNA_reads'],
        #                   file_names['rRNA_sam'],
        #                   file_names['rRNA_mapping_log'],
        #                   file_names['no_rRNA_lengths'],
        #                  )
        #map_tophat(file_names)
        further_filtering(file_names['accepted_hits'],
                          file_names['filtered_bam'],
                          file_names['tRNA_lengths'],
                          file_names['more_rRNA_lengths'],
                          file_names['gtf'],
                         )
        plot_lengths(file_names)

        end_time = time.time()
        print 'Processed in {0:0.2f} s'.format(end_time - start_time)
