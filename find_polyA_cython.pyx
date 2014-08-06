import numpy as np
cimport numpy as np
cimport cython

def find_polyA(char* seq, int min_length):
    ''' Find the first polyA stretch at least min_length long in seq. 
        If there is no such stretch, return the index and length of 
        the longest stretch.
    '''
    cdef int i = 0
    cdef int current_length = 0
    cdef int longest_length = 0
    cdef int longest_start = 0
    cdef int current_start = 0
    cdef int seq_length = len(seq)

    for i in range(seq_length):
        if seq[i] != 'A':
            if current_length >= min_length:
                return current_start, current_length
            else:
                if current_length > longest_length:
                    longest_length = current_length
                    longest_start = current_start
                current_start = i + 1
                current_length = 0
        else:
            current_length += 1
    if current_length >= min_length:
        return current_start, current_length
    else:
        if current_length > longest_length:
            longest_length = current_length
            longest_start = current_start
        return longest_start, longest_length
