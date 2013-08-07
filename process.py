from adapters_cython import *
import fastq
import sam
from collections import Counter
import glob
import os.path
import mapping
import mutations
import numpy as np
import time

def trim_adapters(in_fn, out_fn, lengths_fn):
    up_to_index_rc = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
    lengths = Counter()
    with open(out_fn, 'w') as out_fh:
        for read in fastq.reads(in_fn):
            p = find_adapter_position(read.seq, up_to_index_rc, 10, 3)
            trimmed_slice = slice(None, p)
            seq = read.seq[trimmed_slice]
            qual = read.qual[trimmed_slice]
            lengths[len(seq)] += 1
            if len(seq) > 5:
                out_fh.write(fastq.make_record(read.name, seq, qual))
    
    lengths = mutations.counts_to_array(lengths)
    np.savetxt(lengths_fn, lengths, fmt='%d')

def filter_rRNA(trimmed_fn, no_rRNA_fn, sam_fn, log_fn, lengths_fn):
    mapping.map_bowtie2(trimmed_fn,
                        'yeast_rRNA',
                        sam_fn,
                        error_file_name=log_fn,
                        unaligned_reads_file_name=no_rRNA_fn,
                        threads=8,
                        report_timing=True,
                       )
    
    lengths = Counter(len(read.seq) for read in fastq.reads(no_rRNA_fn))
    lengths = mutations.counts_to_array(lengths)
    np.savetxt(lengths_fn, lengths, fmt='%d')

samples = ['R98S_cDNA_mRNA',
           'R98S_cDNA_sample',
           'Suppressed_R98S_cDNA_mRNA',
           'Suppressed_R98S_cDNA_sample',
           'WT_cDNA_mRNA',
           'WT_cDNA_sample',
          ]

data_dir = '/home/jah/projects/arlen/data'
results_dir = '/home/jah/projects/arlen/results'

for sample in samples:
    print sample
    if not os.path.isdir('{0}/{1}/'.format(results_dir, sample)):
        os.makedirs('{0}/{1}/'.format(results_dir, sample))

    R1_fn = glob.glob('{0}/{1}/{1}*.fastq'.format(data_dir, sample))[0]
    trimmed_fn = '{0}/{1}/{1}_trimmed.fastq'.format(results_dir, sample)
    trimmed_lengths_fn = '{0}/{1}/{1}_trimmed_lengths.txt'.format(results_dir, sample)
    no_rRNA_fn = '{0}/{1}/{1}_no_rRNA.fastq'.format(results_dir, sample)
    rRNA_log_fn = '{0}/{1}/{1}_rRNA.log'.format(results_dir, sample)
    rRNA_sam_fn = '{0}/{1}/{1}_rRNA.sam'.format(results_dir, sample)
    no_rRNA_lengths_fn = '{0}/{1}/{1}_no_rRNA_lengths.txt'.format(results_dir, sample)
    
    start_time = time.time()
    trim_adapters(R1_fn, trimmed_fn, trimmed_lengths_fn)
    filter_rRNA(trimmed_fn, no_rRNA_fn, rRNA_sam_fn, rRNA_log_fn, no_rRNA_lengths_fn)
    end_time = time.time()
    print 'Processed in {0:0.2f} s'.format(end_time - start_time)
