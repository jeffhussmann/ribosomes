def find_poly_A(char *seq):
    cdef int seq_length = len(seq)
    cdef int start

    for start in range(seq_length, 0, -1):
        if seq[start - 1] != 'A' and seq[start - 1] != 'N':
            return start
    return 0

def find_jeff_start(seq):
    start = seq.find('CAGTA')
    if start == 11 or start == 12 or start == 13:
        return start + 5
    else:
        return 0
