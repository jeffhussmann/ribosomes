import numpy as np
import fastq
from trim_cython import find_adapter_position

truseq_R2_rc_up_to_index = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
adapter_prefix_length = 10
max_distance = 3

def get_read_length(fn):
    return len(fastq.reads(fn).next().seq)

def trim_adapters(R1_reads_fn, trimmed_reads_fn, trimmed_lengths_fn):
    read_length = get_read_length(R1_reads_fn)
    trimmed_lengths = np.zeros(read_length + 1, int)
    
    adapter_prefix = truseq_R2_rc_up_to_index[:adapter_prefix_length]
    
    with open(trimmed_reads_fn, 'w') as trimmed_reads_fh:
        for read in fastq.reads(R1_reads_fn):
            p = find_adapter_position(read.seq, adapter_prefix, max_distance)
            trimmed_lengths[p] += 1
            
            if p > 5:
                trimmed_record = fastq.make_record(read.name, 
                                                   read.seq[:p],
                                                   read.seq[:p],
                                                  )
                trimmed_reads_fh.write(trimmed_record)

    np.savetxt(trimmed_lengths_fn, trimmed_lengths, fmt='%d')
