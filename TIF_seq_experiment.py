import matplotlib
matplotlib.use('Agg', warn=False)
import matplotlib.pyplot as plt
import trim
import os
import pysam
import glob
from collections import Counter, defaultdict
from itertools import izip, chain
from Sequencing import fastq, fasta, utilities, mapping_tools, sam
from Sequencing.Parallel import map_reduce, split_file
from Sequencing.utilities import reverse_complement, counts_to_array
from Circles import adapters
from Circles.adapters import find_adapter_positions, simple_hamming_distance
import TIF_seq_structure

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

# For reference:
#CTCC AGGCGGCCGCTAT CACTCTGAGCAATACC
#GTGG AGGCGGCCGCTAT ATCACTCTGAGCAATACC

#  GGTATTGCTCAGAGTG ATAGCGGCCGCCT GGAG
#GGTATTGCTCAGAGTGAT ATAGCGGCCGCCT CCAC

#  GGTATTGCTCAGAGTG   ATAGCGGCCGCCT GGAG
#  GGTATTGCTCAGAGTG   ATATAGCGGCCGCCT CCAC

orientations = ['R1_forward', 'R1_reverse', 'R2_forward', 'R2_reverse']

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
            ('polyA_lengths', 'array_1d', '{name}_polyA_lengths.txt'),
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
             'polyA_lengths',
             'control_ids',
             'mapped_bam_sorted',
            ],
        ]

        specific_work = [
            [(self.extract_boundary_sequences, 'Extracting boundary sequences'),
             (self.map, 'Mapping'),
             (self.filter_mappings, 'Filtering mappings'),
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
        long_enough = 0
    
        counters = {'positions': {orientation: Counter() for orientation in orientations},
                    'control_ids': Counter(),
                    'polyA_lengths': Counter(),
                    'well_formed': 0,
                   }

        with open(self.file_names['five_prime_boundaries'], 'w') as fives_fh, \
             open(self.file_names['three_prime_boundaries'], 'w') as threes_fh:

            for R1, R2 in trimmed_read_pairs:
                total_reads += 1
                five_record, three_record = TIF_seq_structure.find_boundary_sequences(R1, R2, counters)
                if five_record and three_record:
                    long_enough += 1
                    fives_fh.write(five_record)
                    threes_fh.write(three_record)

        for orientation in orientations:
            key = '{0}_{1}'.format(orientation, 'positions')
            array = counts_to_array(counters['positions'][orientation])
            self.write_file(key, array)

        self.write_file('control_ids', counters['control_ids'])
        self.write_file('polyA_lengths', counts_to_array(counters['polyA_lengths']))

        self.log.extend(
            [('Total read pairs', total_reads),
             ('Well-formed read pairs', counters['well_formed']),
             ('Long enough', long_enough),
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
        for orientation in orientations:
            key = '{0}_positions'.format(orientation)
            array = self.read_file(key)
            max_length = max(max_length, len(array))
            ax.plot(array, '.-', label=orientation)

        ax.set_xlim(right=max_length - 1)
        ax.legend(loc='upper right', framealpha=0.5)
        fig.savefig(self.figure_file_names['positions'])

    def plot_polyA_lengths(self):
        fig, ax = plt.subplots()

        array = self.read_file('polyA_lengths')
        ax.plot(array, '.-', label='polyA_lengths')

        ax.legend(loc='upper right', framealpha=0.5)
        fig.savefig(self.figure_file_names['polyA_lengths'])

if __name__ == '__main__':
    script_path = os.path.realpath(__file__)
    map_reduce.controller(TIFSeqExperiment, script_path)
