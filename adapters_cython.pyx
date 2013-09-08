import numpy as np
cimport numpy as np
cimport cython

DTYPEINT = np.int
ctypedef np.int_t DTYPEINT_t

cdef int hamming_distance(char *read,
                          char *adapter,
                          int read_length,
                          int adapter_length,
                          int start,
                         ):
    cdef int compare_length = min(adapter_length, read_length - start)
    cdef int mismatches = 0
    cdef int i

    for i in range(compare_length):
        if read[start + i] != adapter[i]:
            mismatches += 1

    return mismatches

def find_adapter_position(read, char *adapter, int min_comparison_length, int max_distance):
    cdef int read_length = len(read)
    cdef int adapter_length = len(adapter)
    cdef int max_start = len(read) - min_comparison_length
    cdef int distance, start
        
    for start in range(max_start + 1):
        distance = hamming_distance(read, adapter, read_length, adapter_length, start)
        if distance <= max_distance:
            return start
    # Convention: position of read_length means no position was found
    return read_length
    
def characterize_adapters(reads, char *adapter, int read_length):
    cdef int adapter_length = len(adapter)
    cdef int distance, start
    shape = (read_length, adapter_length + 1)
    cdef np.ndarray[DTYPEINT_t, ndim=2] d_array = np.zeros(shape, dtype=DTYPEINT)

    for _, seq, _ in reads:
        for start in range(read_length):
            d = hamming_distance(seq,
                                 adapter,
                                 read_length,
                                 adapter_length,
                                 start,
                                )
            d_array[start, d] += 1

    return d_array
