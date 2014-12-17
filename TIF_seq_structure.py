'''
After ligations, RT, and PCR, the ds product will look like

A                                                                              A

TATAGCGGCCGCTATCACTCTGAGCAATACC-payload-NBAAAAAAAAAAAAAAAACTCCAGGCGGCCGCTATA
ATATCGCCGGCGATAGTGAGACTCGTTATGG-payload-NVTTTTTTTTTTTTTTTTGAGGTCCGCCGGCGATATGTAC

B                                                                              B

TATAGCGGCCGCTATATCACTCTGAGCAATACC-payload-NBAAAAAAAAAAAAAAAAGTGGAGGCGGCCGCTATA
ATATCGCCGGCGATATAGTGAGACTCGTTATGG-payload-NVTTTTTTTTTTTTTTTTCACCTCCGCCGGCGATATGTAC

After digestion with NotI, 

A                                                                              A

      GGCCGCTATCACTCTGAGCAATACC-payload-NBAAAAAAAAAAAAAAAACTCCAGGC          
          CGATAGTGAGACTCGTTATGG-payload-NVTTTTTTTTTTTTTTTTGAGGTCCGCCGG      

After circularization,
A                                                                              A

5'(3' payload)-NBAAAAAAAAAAAAAAAACTCCAGGCGGCCGCTATCACTCTGAGCAATACC-(5' payload)3'
3'(3' payload)-NVTTTTTTTTTTTTTTTTGAGGTCCGCCGGCGATAGTGAGACTCGTTATGG-(5' payload)5'
                                               <------------------
                                               common_right_reverse
                                     <---------
                                  after_right['A']
B                                                                              B

5'(3' payload)-NBAAAAAAAAAAAAAAAAGTGGAGGCGGCCGCTATATCACTCTGAGCAATACC-(5' payload)3'
3'(3' payload)-NVTTTTTTTTTTTTTTTTCACCTCCGCCGGCGATATAGTGAGACTCGTTATGG-(5' payload)5'
                                                 <------------------
                                                 common_right_reverse
                                     <-----------
                                   after_right['B']
'''

from Sequencing.adapters_cython import *
from Sequencing import utilities, fastq
from Sequencing.annotation import Annotation_factory
import find_polyA_cython
import trim

barcodes = {'mp1': 'AGCGCTT',
            'mp5': 'CACTGTT',
            'mp19': 'ATTCCGT',
            'mp22': 'GTATAGT',
            'mp34': 'GCTACCT',
            'mp37': 'CGAAACT',
           }

common_right_reverse = 'GGTATTGCTCAGAGTGATA'
after_right = {'A': 'GCGGCCGCCT',
               'B': 'TAGCGGCCGCCT',
              }
after_left = {'A': 'GGAG',
              'B' :'CCAC',
             }
after_right_length = 4
after_right_prefix = {key: seq[:after_right_length] for key, seq in after_right.items()}

def all_adapter_possibilites(read, adapter):
    completely_contained = find_adapter_positions(read, adapter, len(adapter), 1)
    prefix = find_adapter(adapter, 1, read[-len(adapter):])
    suffix = find_adapter(adapter[::-1], 1, read[:len(adapter)][::-1])
    return completely_contained, prefix, suffix

def find_boundary_sequences(R1, R2, counters):
    # Find which read in the read pair is from the reverse strand by looking for
    # common_right_reverse.
    # First try to find a unique position entirely contained within R1 or R2
    # that is close to common_right_reverse.
    # Failing this, find the longest of (the longest suffix of R1 or R2 that
    # matches a prefix of common_right_reverse) or (the longest prefix of R1 or
    # R2 that matches a suffix of common_right_reverse).

    R1_contained, R1_prefix, R1_suffix = all_adapter_possibilites(R1.seq, common_right_reverse)
    R2_contained, R2_prefix, R2_suffix = all_adapter_possibilites(R2.seq, common_right_reverse)

    if len(R1_contained) + len(R2_contained) > 1:
        # Only one of occurence of common_right_reverse should exist between R1
        # and R2.
        return None, None
    elif len(R1_contained) + len(R2_contained) == 0:
        possiblities = [(len(common_right_reverse) - R1_prefix, 'R1_prefix'),
                        (len(common_right_reverse) - R2_prefix, 'R2_prefix'),
                        (len(common_right_reverse) - R1_suffix, 'R1_suffix'),
                        (len(common_right_reverse) - R2_suffix, 'R2_suffix'),
                       ]
        length, kind = max(possiblities)
        if length > 5:
            if 'R1' in kind:
                reverse_read = R1
                forward_read = R2
                polyA_read = 'R2_forward'
                polyT_read = 'R1_reverse'
            elif 'R2' in kind:
                reverse_read = R2
                forward_read = R1
                polyA_read = 'R1_forward'
                polyT_read = 'R2_reverse'
            if 'prefix' in kind:
                common_right_reverse_start = len(reverse_read.seq) - length
            elif 'suffix' in kind:
                common_right_reverse_start = -length
        else:
            return None, None

    elif len(R1_contained) == 1:
        reverse_read = R1
        forward_read = R2
        polyA_read = 'R2_forward'
        polyT_read = 'R1_reverse'
        common_right_reverse_start = R1_contained.pop()
    elif len(R2_contained) == 1:
        reverse_read = R2
        forward_read = R1
        polyA_read = 'R1_forward'
        polyT_read = 'R2_reverse'
        common_right_reverse_start = R2_contained.pop()

    # '*' means that there was no opportunity to see this id.
    # 'X' means that there was an opportunity and it was neither A nor B.
    right_id = '*'
    left_id = '*'

    five_payload_slice = slice(None, max(0, common_right_reverse_start))
    five_payload_seq = utilities.reverse_complement(reverse_read.seq[five_payload_slice])
    five_payload_qual = reverse_read.qual[five_payload_slice][::-1]

    current_p = common_right_reverse_start + len(common_right_reverse)
    if current_p < len(reverse_read.seq) - after_right_length:
        right_id_seq = reverse_read.seq[current_p:current_p + after_right_length]
        for key, prefix in after_right_prefix.items():
            if right_id_seq == prefix:
                right_id = key
        if right_id == '*':
            right_id = 'X'

        counters['right_ids'][right_id_seq] += 1

        if right_id != 'X':
            current_p += len(after_right[right_id])
            if current_p < len(reverse_read.seq) - 4:
                left_id_seq = reverse_read.seq[current_p:current_p + 4]
                for key, sequence in after_left.items():
                    if left_id_seq == sequence:
                        left_id = key
                if left_id == '*':
                    left_id = 'X'
            
                counters['left_ids'][left_id_seq] += 1

    polyA_start, polyA_length = find_polyA_cython.find_polyA(forward_read.seq, 15)
    polyA_slice = slice(polyA_start, polyA_start + polyA_length)
    polyA_seq = forward_read.seq[polyA_slice]
    polyA_qual = fastq.sanitize_qual(forward_read.qual[polyA_slice])
    three_payload_slice = slice(None, polyA_start)
    three_payload_seq = forward_read.seq[three_payload_slice]
    three_payload_qual = forward_read.qual[three_payload_slice]

    common_name, _ = R1.name.rsplit(':', 1)
    control_ids_string = '{0}-{1}'.format(left_id, right_id)
    five_annotation = trim.PayloadAnnotation(original_name=common_name,
                                             barcode=control_ids_string,
                                             trimmed_seq='',
                                             trimmed_qual='',
                                            )
    three_annotation = trim.PayloadAnnotation(original_name=common_name,
                                              barcode=control_ids_string,
                                              trimmed_seq=polyA_seq,
                                              trimmed_qual=polyA_qual,
                                             )
    five_payload_read = fastq.Read(five_annotation.identifier, five_payload_seq, five_payload_qual)
    three_payload_read = fastq.Read(three_annotation.identifier, three_payload_seq, three_payload_qual)

    counters['positions'][polyT_read][max(0, common_right_reverse_start)] += 1
    counters['positions'][polyA_read][polyA_start] += 1
    counters['joint_lengths'][max(0, common_right_reverse_start), polyA_start] += 1
    counters['polyA_lengths'][polyA_length] += 1
    counters['control_ids'][control_ids_string] += 1

    if polyA_length < 13:
        return None, None

    return five_payload_read, three_payload_read
