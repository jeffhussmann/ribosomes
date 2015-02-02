import numpy as np
from Sequencing import fastq
cimport cython

def length_stratified_composition(seq_info_pairs, int max_seq_length):
    ''' Unnecessary code duplication. '''
    cdef int i, b, seq_length
    cdef char* seq_typed
    relevant_bases = fastq.base_order[:5]
    cdef char* bases = relevant_bases

    shape = (max_seq_length + 1, max_seq_length + 1, 256)
    all_array = np.zeros(shape, int)
    perfect_array = np.zeros(shape, int)

    cdef long[:, :, ::1] all_array_view = all_array
    cdef long[:, :, ::1] perfect_array_view = perfect_array
    
    for seq, perfect_and_unique in seq_info_pairs:
        seq_length = len(seq)
        seq_typed = seq # To avoid 'Obtaining char* from temporary Python value'
        for i in range(seq_length):
            # Automatic type conversion means ord() is unneccesary
            b = seq_typed[i]
            all_array_view[seq_length, i, b] += 1
            if perfect_and_unique:
                perfect_array_view[seq_length, i, b] += 1
        
    # This pulls out only the columns corresponding to possible base
    # identities. 
    all_array = np.dstack([all_array[:, :, b] for b in bases])
    perfect_array = np.dstack([perfect_array[:, :, b] for b in bases])
    
    return all_array, perfect_array
