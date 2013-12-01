import numpy as np
import mapping
import mutations
import fastq
import functools
import trim_cython
import pysam
from itertools import chain

def trim(reads, min_length, max_read_length, trimmed_reads_fn, find_position):
    ''' Wrapper that handles the logistics of trimming reads given a function
        find_postion that takes a sequence and returns a position that trimming
        should occur at.
    '''
    trimmed_lengths = np.zeros(max_read_length + 1, int)
    too_short_lengths = np.zeros(max_read_length + 1, int)
    
    with open(trimmed_reads_fn, 'w') as trimmed_reads_fh:
        for read in reads:
            position = find_position(read.seq)
            if position < min_length:
                too_short_lengths[position] += 1
            else:
                trimmed_lengths[position] += 1
                trimmed_record = fastq.make_record(read.name, 
                                                   read.seq[:position],
                                                   read.qual[:position],
                                                  )
                trimmed_reads_fh.write(trimmed_record)

    return trimmed_lengths, too_short_lengths

truseq_R2_rc = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
linker = 'CTGTAGGCACCATCAAT'

adapter_prefix_length = 10
max_distance = 1

find_truseq = functools.partial(trim_cython.find_adapter,
                                truseq_R2_rc[:adapter_prefix_length],
                                max_distance,
                               )
trim_truseq = functools.partial(trim, find_position=find_truseq)

find_linker = functools.partial(trim_cython.find_adapter,
                                linker[:adapter_prefix_length],
                                max_distance,
                               )
trim_linker = functools.partial(trim, find_position=find_linker)

trim_poly_A = functools.partial(trim, find_position=trim_cython.find_poly_A)

trim_nothing = functools.partial(trim, find_position=len)

def unambiguously_trimmed(filtered_bam_fn, unambiguous_bam_fn, genome_dir):
    ''' Reads that have had poly-As trimmed may have had some real RPF A's
        trimmed as well. Retains only mapped reads for which the next base in
        the reference is a non-A.
    '''
    def first_non_A(seq):
        for i, c in enumerate(seq):
            if c != 'A':
                return i
        return len(seq)

    loaded_genome = mapping.load_genome('saccharomyces_cerevisiae',
                                        explicit_path=genome_dir,
                                       )
    
    bamfile = pysam.Samfile(filtered_bam_fn, 'rb')
    with pysam.Samfile(unambiguous_bam_fn, 'wb', header=bamfile.header) as unambiguous_bam_fh:
        for read in bamfile:
            seqname = bamfile.getrname(read.tid)
            if not read.is_reverse:
                start = read.pos + len(read.seq)
                end = start + 10
                after = loaded_genome[seqname][start:end]
            else:
                start = read.pos - 10
                end = read.pos
                after = loaded_genome[seqname][start:end]
                after = mutations.reverse_complement(after)

            if first_non_A(after) == 0:
                unambiguous_bam_fh.write(read)

    pysam.index(unambiguous_bam_fn)
