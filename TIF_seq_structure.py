from Circles.adapters import find_adapter_positions, simple_hamming_distance
from Sequencing import utilities, fastq
import find_polyA_cython

common_right_reverse = 'GGTATTGCTCAGAGTGATA'
after_right = {'A': 'GCGGCCGCCT',
               'B': 'TAGCGGCCGCCT',
              }
after_left = {'A': 'GGAG',
              'B' :'CCAC',
             }
after_right_length = 4
after_right_prefix = {key: seq[:after_right_length] for key, seq in after_right.items()}

def find_boundary_sequences(R1, R2, counters):
    R1_ps = find_adapter_positions(R1.seq, common_right_reverse, 16, 1)
    R2_ps = find_adapter_positions(R2.seq, common_right_reverse, 16, 1)

    if len(R1_ps) + len(R2_ps) != 1:
        return None, None
    elif len(R1_ps) == 1:
        reverse_read = R1
        forward_read = R2
        polyA_read = 'R2_forward'
        p = R1_ps.pop()
        counters['positions']['R1_reverse'][p] += 1
    elif len(R2_ps) == 1:
        reverse_read = R2
        forward_read = R1
        polyA_read = 'R1_forward'
        p = R2_ps.pop()
        counters['positions']['R2_reverse'][p] += 1

    right_id = 'N'
    left_id = 'N'

    five_slice = slice(None, p)
    five_seq = utilities.reverse_complement(reverse_read.seq[five_slice])
    five_qual = reverse_read.qual[five_slice][::-1]

    after_p = p + len(common_right_reverse)
    if after_p < len(reverse_read.seq) - after_right_length:
        right_id_seq = reverse_read.seq[after_p:after_p + after_right_length]
        for key, prefix in after_right_prefix.items():
            if right_id_seq == prefix:
                right_id = key

    if right_id != 'N':
        after_p += len(after_right[right_id])
        if after_p < len(reverse_read.seq) - 4:
            left_id_seq = reverse_read.seq[after_p: after_p + 4]
            for key, sequence in after_left.items():
                if left_id_seq == sequence:
                    left_id = key

    polyA_start, poly_A_length = find_polyA_cython.find_polyA(forward_read.seq, 15)
    three_slice = slice(None, polyA_start)
    three_seq = forward_read.seq[three_slice]
    three_qual = forward_read.qual[three_slice]

    # Remove trailing read number identifier to allow IGV 'View as pairs'
    common_name, _ = R1.name.rsplit('.', 1)
    control_ids_string = '{0}-{1}'.format(left_id, right_id)
    name = '{0}_{1}'.format(common_name, control_ids_string)
    five_record = fastq.make_record(name, five_seq, five_qual)
    three_record = fastq.make_record(name, three_seq, three_qual)

    counters['positions'][polyA_read][polyA_start] += 1
    counters['polyA_lengths'][poly_A_length] += 1
    counters['control_ids'][control_ids_string] += 1
    counters['well_formed'] += 1

    if p > 12 and polyA_start > 12: 
        return five_record, three_record
    else:
        return None, None
