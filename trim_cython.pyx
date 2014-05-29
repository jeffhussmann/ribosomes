cdef int hamming_distance(char *seq,
                          char *adapter,
                          int seq_length,
                          int adapter_length,
                          int start,
                         ):
    ''' Returns the hamming distance between the overlap of seq[start:] and
        adapter.
    '''
    cdef int compare_length = min(adapter_length, seq_length - start)
    cdef int mismatches = 0
    cdef int i

    for i in range(compare_length):
        if seq[start + i] != adapter[i]:
            mismatches += 1

    return mismatches

def find_adapter(char *adapter, int max_distance, char *seq):
    ''' Returns the leftmost position in seq for which seq[position:] is within
        hamming distance max_distance of adapter. Only checks positions for
        which adapter is completely contained in seq[position:].
    '''
    cdef int seq_length = len(seq)
    cdef int adapter_length = len(adapter)
    cdef int max_start = len(seq) - adapter_length
    cdef int distance, start
        
    for start in range(max_start + 1):
        distance = hamming_distance(seq, adapter, seq_length, adapter_length, start)
        if distance <= max_distance:
            return start
    
    # Convention: position of seq_length means no position was found
    return seq_length

def find_medium_adapter(char *adapter, char *seq):
    ''' Searches for exact matches of a suffix of seq with a prefix of adapter.
        Allows 1 mismatch if the compared regions are long enough.
    '''
    cdef int seq_length = len(seq)
    cdef int adapter_length = len(adapter)
    cdef int start

    for start in range(seq_length):
        distance = hamming_distance(seq, adapter, seq_length, adapter_length, start)
        if distance == 0:
            return start
        elif seq_length - start >= 5 and distance <= 1:
            return start

    return seq_length

def find_short_adapter(char *adapter, char *seq):
    ''' Searches for exact matches of a suffix of seq with a prefix of adapter.
    '''
    cdef int seq_length = len(seq)
    cdef int adapter_length = len(adapter)
    cdef int start

    for start in range(seq_length):
        distance = hamming_distance(seq, adapter, seq_length, adapter_length, start)
        if distance == 0:
            return start

    return seq_length

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
