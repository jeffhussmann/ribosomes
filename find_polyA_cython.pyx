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

def predominantly_A(char* seq):
    ''' Returns True if at least 19 of the first first 20 bases in seq are A if
    seq is at least 20 bases long, or if seq is all A's if it is less than 20
    bases long.
    '''
    
    cdef int max_i = min(20, len(seq))
    cdef int A_count = 0

    for i in range(max_i):
        if seq[i] == 'A':
            A_count += 1

    if max_i < 20:
        if A_count == max_i:
            return True
        else:
            return False
    else:
        if A_count >= 19:
            return True
        else:
            return False
