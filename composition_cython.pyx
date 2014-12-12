import numpy as np
from Sequencing import fastq
cimport numpy as np
cimport cython

DTYPEINT = np.int
ctypedef np.int_t DTYPEINT_t

def length_stratified_composition(seq_info_pairs, int max_seq_length):
    ''' Unnecessary code duplication. '''
    cdef int i, q, b, seq_length
    cdef char* seq_typed
    cdef char* bases = fastq.base_order[:5]

    shape = (max_seq_length + 1, max_seq_length + 1, 256)
    cdef np.ndarray[DTYPEINT_t, ndim=3] all_array = np.zeros(shape, dtype=DTYPEINT)
    cdef np.ndarray[DTYPEINT_t, ndim=3] perfect_array = np.zeros(shape, dtype=DTYPEINT)
    
    for seq, perfect_and_unique in seq_info_pairs:
        seq_length = len(seq)
        seq_typed = seq # To avoid 'Obtaining char* from temporary Python value'
        for i in range(seq_length):
            # Automatic type conversion means ord() is unneccesary
            b = seq_typed[i]
            all_array[seq_length, i, b] += 1
            if perfect_and_unique:
                perfect_array[seq_length, i, b] += 1
        
    # This pulls out only the columns corresponding to possible base
    # identities. 
    all_array = np.dstack([all_array[:, :, b] for b in bases])
    perfect_array = np.dstack([perfect_array[:, :, b] for b in bases])
    
    return all_array, perfect_array
