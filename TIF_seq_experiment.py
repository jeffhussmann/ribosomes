import matplotlib
matplotlib.use('Agg', warn=False)
import matplotlib.pyplot as plt
from Circles import fastq, fasta, utilities, adapters, mapping_tools, sam
from Circles.Parallel import map_reduce, split_file
from collections import Counter
from Circles.adapters import find_adapter_positions, simple_hamming_distance
from Circles.utilities import reverse_complement, counts_to_array
import trim
import os
import pysam
import glob
from itertools import izip, chain
from collections import defaultdict

mp1  = 'AGCGCTT'
mp5  = 'CACTGTT'
mp19 = 'ATTCCGT'
mp22 = 'GTATAGT'
mp34 = 'GCTACCT'
mp37 = 'CGAAACT'

forwards = {'middle': ('AGGCGGCCGCTAT', 1),
            'left_A': ('CTCC', 0),
            'left_B': ('GTGG', 0),
            'right_A': ('CACTCTGAGCAATACC', 2),
            'right_B': ('ATCACTCTGAGCAATACC', 2),
           }

reverses = {}
for key, (string, mismatches) in forwards.items():
    if 'left' in key:
        new_key = key.replace('left', 'right')
    elif 'right' in key:
        new_key = key.replace('right', 'left')
    else:
        new_key = key
    reverses[new_key] = (utilities.reverse_complement(string), mismatches)
        
orientations = ['R1_forward', 'R1_reverse', 'R2_forward', 'R2_reverse']

def valid_match(possible_seq, expected_seq, allowed_mismatches):
    return (len(possible_seq) == len(expected_seq) and
            simple_hamming_distance(possible_seq, expected_seq) <= allowed_mismatches)

def find_structure(seq, adapters):
    ''' Strategy: look for common middle or its RC, then see which of A or B is
        on right and left.
    '''
    possibilities = []

    middle_seq, middle_mismatches = adapters['middle']
    ps = find_adapter_positions(seq, middle_seq, len(middle_seq), middle_mismatches)
    for p in ps:
        bad_left = False
        bad_right = False

        left_A_seq, left_A_mismatches = adapters['left_A']
        left_B_seq, left_B_mismatches = adapters['left_B']
        possible_left_A = seq[p - len(left_A_seq):p]
        possible_left_B = seq[p - len(left_B_seq):p]

        if valid_match(possible_left_A, left_A_seq, left_A_mismatches):
            left_control_id = 'A'
            left_p = p - len(left_A_seq)
        elif valid_match(possible_left_B, left_B_seq, left_B_mismatches):
            left_control_id = 'B'
            left_p = p - len(left_B_seq)
        else:
            bad_left = True
            left_control_id = possible_left_A

        right_A_seq, right_A_mismatches = adapters['right_A']
        right_B_seq, right_B_mismatches = adapters['right_B']
        after_middle = p + len(middle_seq)
        possible_right_A = seq[after_middle:after_middle + len(right_A_seq)]
        possible_right_B = seq[after_middle:after_middle + len(right_B_seq)]
        
        # right_p is the index of the last base in 
        if valid_match(possible_right_A, right_A_seq, right_A_mismatches):
            right_control_id = 'A'
            right_p = after_middle + len(right_A_seq) - 1
        elif valid_match(possible_right_B, right_B_seq, right_B_mismatches):
            right_control_id = 'B'
            right_p = after_middle + len(right_B_seq) - 1
        else:
            bad_right = True
            right_control_id = possible_right_A

        if not (bad_left or bad_right):
            possibility = [(left_p, right_p),
                           (left_control_id, right_control_id),
                          ]
            possibilities.append(possibility)

    if len(possibilities) > 1:
        #raise ValueError(seq, possibilities)
        return None
    elif len(possibilities) == 0:
        return None
    else:
        return possibilities.pop()

def find_boundary_sequences(R1, R2, counters):
    possibilities = {'R1_forward': find_structure(R1.seq, forwards),
                     'R1_reverse': find_structure(R1.seq, reverses),
                     'R2_forward': find_structure(R2.seq, forwards),
                     'R2_reverse': find_structure(R2.seq, reverses),
                    }
    #TODO: use instances of R1/R2 overlap to extend UTR sequences
    # For now, only valid if possible in only one.

    valids = [key for key, possibility in possibilities.items() if possibility]

    if len(valids) == 1:
        orientation = valids.pop()
        (left_p, right_p), (left_control_id, right_control_id) = possibilities[orientation]
        if 'R1' in orientation:
            seq = R1.seq
            qual = R1.qual
        elif 'R2' in orientation:
            seq = R2.seq
            qual = R2.qual

        counters['positions'][orientation][left_p] += 1

        left_slice = slice(None, left_p)
        left_seq = seq[left_slice]
        left_qual = qual[left_slice]
        right_slice = slice(right_p + 1, None)
        right_seq = seq[right_slice]
        right_qual = qual[right_slice]

        if 'forward' in orientation:
            three_seq, three_qual = left_seq, left_qual
            five_seq, five_qual = right_seq, right_qual
            control_id_1, control_id_2 = left_control_id, right_control_id
        elif 'reverse' in orientation:
            three_seq, three_qual = reverse_complement(right_seq), right_qual[::-1]
            five_seq, five_qual = reverse_complement(left_seq), left_qual[::-1]
            control_id_1, control_id_2 = right_control_id, left_control_id
            if control_id_1 != 'A' and control_id_1 != 'B':
                control_id_1 = reverse_complement(control_id_1)
            if control_id_2 != 'A' and control_id_2 != 'B':
                control_id_2 = reverse_complement(control_id_2)

        control_ids_string = '{0}-{1}'.format(control_id_1, control_id_2)
        counters['control_ids'][control_ids_string] += 1

        remove_poly_A_slice = slice(None, trim.find_poly_A(three_seq))
        three_seq_trimmed = three_seq[remove_poly_A_slice]
        three_qual_trimmed = three_qual[remove_poly_A_slice]
        polyA_length = len(three_seq) - len(three_seq_trimmed)
        counters['polyA_lengths'][orientation][polyA_length] += 1
        
        if len(five_seq) > 12 and len(three_seq_trimmed) > 12:
            # Remove trailing read number identifier to allow IGV 'View as pairs'
            common_name, _ = R1.name.rsplit('.', 1)
            name = '{0}_{1}'.format(common_name, control_ids_string)
            five_record = fastq.make_record(name, five_seq, five_qual)
            three_record = fastq.make_record(name, three_seq_trimmed, three_qual_trimmed)
        else:
            five_record = None
            three_record = None

        return orientation, possibilities[orientation], five_record, three_record
    else:
        return None

def trim_read_pairs(read_pairs):
    num_to_trim = len(mp1)
    def trim_read(read):
        return fastq.Read(read.name, read.seq[num_to_trim:], read.qual[num_to_trim:])
    for R1, R2 in read_pairs:
        yield trim_read(R1), trim_read(R2)

class TIFSeqExperiment(map_reduce.MapReduceExperiment):
    num_stages = 1

    def __init__(self, **kwargs):
        map_reduce.MapReduceExperiment.__init__(self, **kwargs)

        self.data_dir = kwargs['data_dir']

        self.bowtie2_index = kwargs['bowtie2_index']

        specific_results_files = [
            ('five_prime_boundaries', 'fastq', '{name}_five_prime_boundaries.fastq'),
            ('three_prime_boundaries', 'fastq', '{name}_three_prime_boundaries.fastq'),
            ('R1_forward_positions', 'array_1d', '{name}_R1_forward_positions.txt'),
            ('R1_reverse_positions', 'array_1d', '{name}_R1_reverse_positions.txt'),
            ('R2_forward_positions', 'array_1d', '{name}_R2_forward_positions.txt'),
            ('R2_reverse_positions', 'array_1d', '{name}_R2_reverse_positions.txt'),
            ('R1_forward_polyA_lengths', 'array_1d', '{name}_R1_forward_polyA_lengths.txt'),
            ('R1_reverse_polyA_lengths', 'array_1d', '{name}_R1_reverse_polyA_lengths.txt'),
            ('R2_forward_polyA_lengths', 'array_1d', '{name}_R2_forward_polyA_lengths.txt'),
            ('R2_reverse_polyA_lengths', 'array_1d', '{name}_R2_reverse_polyA_lengths.txt'),
            ('control_ids', 'counts', '{name}_control_ids.txt'),
            ('mapped_sam', 'sam_unsorted', '{name}_mapped.sam'),
            ('mapped_bam', 'bam', '{name}_mapped.bam'),
            ('mapped_bam_sorted', 'bam', '{name}_mapped_sorted.bam'),
        ]

        specific_figure_files = [
            ('positions', '{name}_positions.pdf'),
            ('polyA_lengths', '{name}_polyA_lengths.pdf'),
        ]

        specific_outputs = [
            ['R1_forward_positions',
             'R1_reverse_positions',
             'R2_forward_positions',
             'R2_reverse_positions',
             'R1_forward_polyA_lengths',
             'R1_reverse_polyA_lengths',
             'R2_forward_polyA_lengths',
             'R2_reverse_polyA_lengths',
             'control_ids',
             'mapped_bam_sorted',
            ],
        ]

        specific_work = [
            [#(self.extract_boundary_sequences, 'Extracting boundary sequences'),
             #(self.map, 'Mapping'),
             #(self.filter_mappings, 'Filtering mappings'),
            ],
        ]

        specific_cleanup = [
            [self.plot_positions,
             self.plot_polyA_lengths,
             self.index_bam,
            ],
        ]

        self.results_files.extend(specific_results_files)
        self.figure_files.extend(specific_figure_files)
        map_reduce.extend_stages(self.outputs, specific_outputs)
        map_reduce.extend_stages(self.work, specific_work)
        map_reduce.extend_stages(self.cleanup, specific_cleanup)

        self.make_file_names()
    
    def get_read_pairs(self):
        data_fns = glob.glob(self.data_dir + '/*.fastq') + glob.glob(self.data_dir + '/*.fq')
        fn_pairs = defaultdict(lambda: {1: None, 2: None})
        for data_fn in data_fns:
            head, tail = os.path.split(data_fn)
            root, ext = os.path.splitext(tail)
            # Expect root to end in either 1 or 2
            prefix, which_member = root[:-1], int(root[-1])
            fn_pairs[prefix][which_member] = data_fn 

        read_pairs_list = []

        for prefix in sorted(fn_pairs):
            R1_fn = fn_pairs[prefix][1]
            R2_fn = fn_pairs[prefix][2]

            if R1_fn == None or R2_fn == None:
                raise ValueError('unpaired file names in data_dir')

            R1_lines = split_file.piece(R1_fn, self.num_pieces, self.which_piece, 'fastq')
            R2_lines = split_file.piece(R2_fn, self.num_pieces, self.which_piece, 'fastq')
            read_pairs = fastq.read_pairs(R1_lines, R2_lines, ensure_sanger_encoding=True)
            read_pairs_list.append(read_pairs)
        
        all_read_pairs = chain.from_iterable(read_pairs_list)

        return all_read_pairs

    def extract_boundary_sequences(self):
        read_pairs = self.get_read_pairs()
        trimmed_read_pairs = trim_read_pairs(read_pairs)

        total_reads = 0
    
        counters = {'positions': {orientation: Counter() for orientation in orientations},
                    'control_ids': Counter(),
                    'polyA_lengths': {orientation: Counter() for orientation in orientations},
                   }

        with open(self.file_names['five_prime_boundaries'], 'w') as fives_fh, \
             open(self.file_names['three_prime_boundaries'], 'w') as threes_fh:

            for R1, R2 in trimmed_read_pairs:
                total_reads += 1
                boundary_sequences = find_boundary_sequences(R1, R2, counters)
                if boundary_sequences:
                    orientation, structure, five_record, three_record = boundary_sequences
                    if five_record and three_record:
                        fives_fh.write(five_record)
                        threes_fh.write(three_record)

        for suffix in ['positions', 'polyA_lengths']:
            for orientation in orientations:
                key = '{0}_{1}'.format(orientation, suffix)
                array = counts_to_array(counters[suffix][orientation])
                self.write_file(key, array)

        self.write_file('control_ids', counters['control_ids'])

        self.log.extend(
            [('Total read pairs', total_reads),
            ],
        )

    def map(self):
        mapping_tools.map_bowtie2_paired(self.file_names['five_prime_boundaries'],
                                         self.file_names['three_prime_boundaries'],
                                         self.bowtie2_index,
                                         self.file_names['mapped_sam'],
                                         threads=1,
                                         forward_forward=True,
                                         max_insert_size=10000,
                                        )

    def filter_mappings(self):
        sam_file = pysam.Samfile(self.file_names['mapped_sam'])
        with pysam.Samfile(self.file_names['mapped_bam'], 'wb', template=sam_file) as bam_file:
            read_pairs = izip(*[sam_file]*2)
            for R1, R2 in read_pairs:
                if R1.qname != R2.qname:
                    # Ensure that the iteration through pairs is in sync.
                    raise ValueError
                if R1.is_unmapped or \
                   R2.is_unmapped or \
                   R1.mapq < 40 or \
                   R2.mapq < 40 or \
                   R1.tid != R2.tid or \
                   abs(R1.tlen) > 10000:
                    continue

                bam_file.write(R1)
                bam_file.write(R2)

        sam.sort_bam(self.file_names['mapped_bam'], self.file_names['mapped_bam_sorted'])

    def index_bam(self):
        pysam.index(self.file_names['mapped_bam_sorted'])

    def plot_positions(self):
        fig, ax = plt.subplots()

        max_length = 0
        for orientation in sorted(orientations):
            key = '{0}_positions'.format(orientation)
            array = self.read_file(key)
            max_length = max(max_length, len(array))
            ax.plot(array, '.-', label=orientation)

        ax.set_xlim(right=max_length - 1)
        ax.legend(loc='upper right', framealpha=0.5)
        fig.savefig(self.figure_file_names['positions'])

    def plot_polyA_lengths(self):
        fig, ax = plt.subplots()

        for orientation in sorted(orientations):
            key = '{0}_polyA_lengths'.format(orientation)
            array = self.read_file(key)
            ax.plot(array, '.-', label=orientation)

        ax.legend(loc='upper right', framealpha=0.5)
        fig.savefig(self.figure_file_names['polyA_lengths'])

if __name__ == '__main__':
    script_path = os.path.realpath(__file__)
    map_reduce.controller(TIFSeqExperiment, script_path)
