import numpy as np
import mapping
import mutations
import fastq
import functools
import trim_cython
import pysam
from itertools import chain

def get_read_length(fns):
    lengths = [len(fastq.reads(fn).next().seq) for fn in fns]
    return max(lengths)

def get_all_reads(fns):
    reads_list = [fastq.reads(fn) for fn in fns]
    return chain.from_iterable(reads_list)

def trim(R1_reads_fns, trimmed_reads_fn, trimmed_lengths_fn, find_position):
    read_length = get_read_length(R1_reads_fns)
    trimmed_lengths = np.zeros(read_length + 1, int)
    
    with open(trimmed_reads_fn, 'w') as trimmed_reads_fh:
        for read in get_all_reads(R1_reads_fns):
            p = find_position(read.seq)
            trimmed_lengths[p] += 1
            if p > 5:
                trimmed_record = fastq.make_record(read.name, 
                                                   read.seq[:p],
                                                   read.qual[:p],
                                                  )
                trimmed_reads_fh.write(trimmed_record)

    np.savetxt(trimmed_lengths_fn, trimmed_lengths, fmt='%d')

truseq_R2_rc = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
adapter_prefix_length = 10
max_distance = 3
find_adapter = functools.partial(trim_cython.find_adapter,
                                 truseq_R2_rc[:adapter_prefix_length],
                                 max_distance,
                                )
trim_adapters = functools.partial(trim, find_position=find_adapter)

trim_poly_A = functools.partial(trim, find_position=trim_cython.find_poly_A)

trim_nothing = functools.partial(trim, find_position=len)

def unambiguously_trimmed(filtered_bam_fn, unambiguous_bam_fn, genome_dir):
    ''' Reads that have had poly-As trimmed may have had some real RPF As
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
