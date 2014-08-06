import matplotlib
matplotlib.use('Agg', warn=False)
import matplotlib.pyplot as plt
import trim
import os
import pysam
import glob
from collections import Counter
from itertools import izip, chain
from Sequencing import fastq, utilities, mapping_tools, sam, genomes
from Circles.annotation import Annotation_factory
from Sequencing.Parallel import map_reduce, split_file
from Sequencing.Serialize import array_1d, array_2d, counts, sparse_joint_counts
from Serialize import read_positions
from Sequencing.utilities import counts_to_array
import TIF_seq_structure
import rna_experiment
import positions
import visualize
import three_p_experiment

orientations = ['R1_forward', 'R1_reverse', 'R2_forward', 'R2_reverse']

class TIFSeqExperiment(rna_experiment.RNAExperiment):
    num_stages = 2

    specific_results_files = [
        ('five_prime_boundaries', 'fastq', '{name}_five_prime_boundaries.fastq'),
        ('three_prime_boundaries', 'fastq', '{name}_three_prime_boundaries.fastq'),
        ('R1_forward_positions', array_1d, '{name}_R1_forward_positions.txt'),
        ('R1_reverse_positions', array_1d, '{name}_R1_reverse_positions.txt'),
        ('R2_forward_positions', array_1d, '{name}_R2_forward_positions.txt'),
        ('R2_reverse_positions', array_1d, '{name}_R2_reverse_positions.txt'),
        ('polyA_lengths', array_1d, '{name}_polyA_lengths.txt'),
        ('joint_lengths', array_2d, '{name}_joint_lengths.txt'),
        ('control_ids', counts, '{name}_control_ids.txt'),
        ('left_ids', counts, '{name}_left_ids.txt'),
        ('right_ids', counts, '{name}_right_ids.txt'),

        ('five_prime_tophat_dir', 'dir', 'tophat_five_prime'),
        ('five_prime_accepted_hits', 'bam', 'tophat_five_prime/accepted_hits.bam'),
        ('five_prime_unmapped', 'bam', 'tophat_five_prime/unmapped.bam'),
        ('five_prime_sorted_by_name', 'bam', '{name}_five_prime_by_name.bam'),
        ('five_prime_read_positions', read_positions, '{name}_five_prime_read_positions.hdf5'),

        ('three_prime_tophat_dir', 'dir', 'tophat_three_prime'),
        ('three_prime_accepted_hits', 'bam', 'tophat_three_prime/accepted_hits.bam'),
        ('three_prime_unmapped', 'bam', 'tophat_three_prime/unmapped.bam'),
        ('three_prime_sorted_by_name', 'bam', '{name}_three_prime_by_name.bam'),
        ('three_prime_read_positions', read_positions, '{name}_three_prime_read_positions.hdf5'),

        ('combined', 'bam', '{name}_combined.bam'),
        ('combined_extended', 'bam', '{name}_combined_extended.bam'),
        ('combined_extended_sorted', 'bam', '{name}_combined_extended_sorted.bam'),
        
        ('nongenomic_lengths', array_1d, '{name}_nongenomic_lengths.txt'),

        ('joint_positions', sparse_joint_counts, '{name}_joint_positions.txt'),
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
         'left_ids',
         'right_ids',
         'joint_lengths',
         'combined_extended_sorted',
        ],
        ['five_prime_read_positions',
         'three_prime_read_positions',
         'from_starts_and_ends',
         'joint_positions',
        ],
    ]

    specific_work = [
        ['extract_boundary_sequences',
         'map_tophat',
         'sort_by_name',
         'combine_mappings',
         'extend_mappings',
        ],
        ['get_read_positions',
         'get_metagene_positions',
        ],
    ]

    specific_cleanup = [
        ['plot_positions',
         'plot_polyA_lengths',
        ],
        ['plot_starts_and_ends',
        ],
    ]

    def __init__(self, **kwargs):
        super(TIFSeqExperiment, self).__init__(**kwargs)

        self.min_payload_length = 12
        
    def trim_barcodes(self, read_pairs):
        num_to_trim = len(TIF_seq_structure.barcodes['mp1'])
        def trim_read(read):
            return fastq.Read(read.name, read.seq[num_to_trim:], read.qual[num_to_trim:])
        for R1, R2 in read_pairs:
            yield trim_read(R1), trim_read(R2)

    def extract_boundary_sequences(self):
        read_pairs = self.get_read_pairs()
        trimmed_read_pairs = self.trim_barcodes(read_pairs)

        total_reads = 0
        well_formed = 0
        long_enough = 0
    
        counters = {'positions': {orientation: Counter() for orientation in orientations},
                    'control_ids': Counter(),
                    'polyA_lengths': Counter(),
                    'left_ids': Counter(),
                    'right_ids': Counter(),
                    'joint_lengths': Counter(),
                   }

        with open(self.file_names['five_prime_boundaries'], 'w') as fives_fh, \
             open(self.file_names['three_prime_boundaries'], 'w') as threes_fh:

            for R1, R2 in trimmed_read_pairs:
                total_reads += 1
                five_payload_read, three_payload_read = TIF_seq_structure.find_boundary_sequences(R1, R2, counters)
                if five_payload_read and three_payload_read:
                    well_formed += 1
                    if len(five_payload_read.seq) >= self.min_payload_length and \
                       len(three_payload_read.seq) >= self.min_payload_length:
                        long_enough += 1
                        fives_fh.write(fastq.make_record(*five_payload_read))
                        threes_fh.write(fastq.make_record(*three_payload_read))

        for orientation in orientations:
            key = '{0}_{1}'.format(orientation, 'positions')
            array = counts_to_array(counters['positions'][orientation])
            self.write_file(key, array)

        self.write_file('left_ids', counters['left_ids'])
        self.write_file('right_ids', counters['right_ids'])
        self.write_file('control_ids', counters['control_ids'])
        self.write_file('polyA_lengths', counts_to_array(counters['polyA_lengths']))
        self.write_file('joint_lengths', counts_to_array(counters['joint_lengths'], dim=2))

        self.log.extend(
            [('Total read pairs', total_reads),
             ('Well-formed', well_formed),
             ('Long enough', long_enough),
            ],
        )

    def map_tophat(self):
        mapping_tools.map_tophat([self.file_names['five_prime_boundaries']],
                                 self.file_names['bowtie2_index_prefix'],
                                 self.file_names['genes'],
                                 self.file_names['transcriptome_index'],
                                 self.file_names['five_prime_tophat_dir'],
                                )
        mapping_tools.map_tophat([self.file_names['three_prime_boundaries']],
                                 self.file_names['bowtie2_index_prefix'],
                                 self.file_names['genes'],
                                 self.file_names['transcriptome_index'],
                                 self.file_names['three_prime_tophat_dir'],
                                )

    def sort_by_name(self):
        for edge in ['five', 'three']:
            accepted_fn = self.file_names['{}_prime_accepted_hits'.format(edge)]
            accepted_prefix, _ = os.path.splitext(accepted_fn)
            accepted_sorted_fn = '{}_sorted.bam'.format(accepted_prefix)

            unmapped_fn = self.file_names['{}_prime_unmapped'.format(edge)]
            unmapped_prefix, _ = os.path.splitext(unmapped_fn)
            unmapped_sorted_fn = '{}_sorted.bam'.format(unmapped_prefix)

            merged_fn = self.file_names['{}_prime_sorted_by_name'.format(edge)]

            sam.sort_bam(accepted_fn, accepted_sorted_fn, by_name=True)
            sam.sort_bam(unmapped_fn, unmapped_sorted_fn, by_name=True)
            pysam.merge('-nf', merged_fn, accepted_sorted_fn, unmapped_sorted_fn)

    def combine_mappings(self):
        num_unmapped = 0
        num_five_unmapped = 0
        num_three_unmapped = 0
        num_nonunique = 0
        num_discordant = 0
        num_concordant = 0

        five_prime_mappings = pysam.Samfile(self.file_names['five_prime_sorted_by_name'])
        three_prime_mappings = pysam.Samfile(self.file_names['three_prime_sorted_by_name'])

        five_prime_grouped = utilities.group_by(five_prime_mappings, lambda m: m.qname)
        three_prime_grouped = utilities.group_by(three_prime_mappings, lambda m: m.qname)
        with pysam.Samfile(self.file_names['combined'], 'wb', template=five_prime_mappings) as combined_bam_file:
            group_pairs = izip(five_prime_grouped, three_prime_grouped)
            for (five_qname, five_group), (three_qname, three_group) in group_pairs:
                five_annotation = trim.PayloadAnnotation.from_identifier(five_qname)
                three_annotation = trim.PayloadAnnotation.from_identifier(three_qname)
                if five_annotation['original_name'] != three_annotation['original_name']:
                    # Ensure that the iteration through pairs is in sync.
                    print five_qname, three_qname
                    raise ValueError

                five_unmapped = any(m.is_unmapped for m in five_group)
                three_unmapped = any(m.is_unmapped for m in three_group)
                if five_unmapped:
                    num_five_unmapped += 1
                if three_unmapped:
                    num_three_unmapped += 1
                if five_unmapped or three_unmapped:
                    num_unmapped += 1
                    continue

                five_nonunique = len(five_group) > 1 or any(m.mapq < 40 for m in five_group)
                three_nonunique = len(three_group) > 1 or any(m.mapq < 40 for m in three_group)
                if five_nonunique or three_nonunique:
                    num_nonunique += 1
                    continue
                
                five_m = five_group.pop()
                three_m = three_group.pop()

                five_strand = '-' if five_m.is_reverse else '+'
                three_strand = '-' if three_m.is_reverse else '+'

                tlen = max(five_m.aend, three_m.aend) - min(five_m.pos, three_m.pos)
                discordant = (five_m.tid != three_m.tid) or (five_strand) != (three_strand) or (tlen > 10000) 
                if discordant:
                    num_discordant += 1
                    continue
                
                if five_strand == '+':
                    first_read = five_m
                    second_read = three_m
                elif five_strand == '-':
                    first_read = three_m
                    second_read = five_m
                
                gap = second_read.pos - first_read.aend
                if gap < 0:
                    num_discordant += 1
                    continue
                
                combined_read = pysam.AlignedRead()
                # qname needs to come from three_m to include trimmed As
                combined_read.qname = three_m.qname
                combined_read.tid = five_m.tid
                combined_read.seq = first_read.seq + second_read.seq
                combined_read.qual = first_read.qual + second_read.qual
                combined_read.cigar = first_read.cigar + [(3, gap)] + second_read.cigar
                combined_read.pos = first_read.pos
                combined_read.is_reverse = first_read.is_reverse
                combined_read.mapq = min(first_read.mapq, second_read.mapq)
                combined_read.rnext = -1
                combined_read.pnext = -1
                
                num_concordant += 1

                combined_bam_file.write(combined_read)

        self.log.extend(
            [('Unmapped', num_unmapped),
             ('Five prime unmapped', num_five_unmapped),
             ('Three prime unmapped', num_three_unmapped),
             ('Nonunique', num_nonunique),
             ('Discordant', num_discordant),
             ('Concordant', num_concordant),
            ],
        )

    def extend_mappings(self):
        trim.extend_polyA_ends(self.file_names['combined'],
                               self.file_names['combined_extended'],
                               self.file_names['genome'],
                              )

        sam.sort_bam(self.file_names['combined_extended'],
                     self.file_names['combined_extended_sorted'],
                    )
    
    def get_read_positions(self):
        piece_CDSs, max_gene_length = self.get_CDSs()
        five_prime_gene_infos = positions.get_Transcript_position_counts(self.merged_file_names['combined_extended_sorted'],
                                                                         piece_CDSs,
                                                                         relevant_lengths=[],
                                                                         left_buffer=500,
                                                                         right_buffer=500,
                                                                        )
        
        self.five_prime_read_positions = {name: info['position_counts']
                                          for name, info in five_prime_gene_infos.iteritems()}

        self.write_file('five_prime_read_positions', self.five_prime_read_positions)
        

        three_prime_gene_infos = positions.get_Transcript_polyA_position_counts(self.merged_file_names['combined_extended_sorted'],
                                                                                piece_CDSs,
                                                                                max_relevant_length=5,
                                                                                left_buffer=500,
                                                                                right_buffer=500,
                                                                               )
        self.three_prime_read_positions = {name: info['position_counts']
                                          for name, info in three_prime_gene_infos.iteritems()}

        self.write_file('three_prime_read_positions', self.three_prime_read_positions)

        joint_position_counts = {}
        for transcript in piece_CDSs:
            counts = positions.get_joint_position_counts_sparse(self.merged_file_names['combined_extended_sorted'],
                                                                transcript,
                                                                left_buffer=500,
                                                                right_buffer=500,
                                                               )
            joint_position_counts[transcript.name] = counts

        self.write_file('joint_positions', joint_position_counts)

    def get_metagene_positions(self):
        piece_CDSs, max_gene_length = self.get_CDSs()

        five_prime_read_positions = self.load_read_positions(modifier='five_prime')
        three_prime_read_positions = self.load_read_positions(modifier='three_prime')

        processed_read_positions = {}
        for name in five_prime_read_positions:
            gene = {'five_prime': five_prime_read_positions[name]['all'],
                    'three_prime_genomic': three_prime_read_positions[name][0],
                    'three_prime_nongenomic': three_prime_read_positions[name]['all'] - three_prime_read_positions[name][0],
                   }
            processed_read_positions[name] = gene

        from_starts_and_ends = positions.compute_metagene_positions(processed_read_positions, max_gene_length)
        self.write_file('from_starts_and_ends', from_starts_and_ends)
    
    def plot_starts_and_ends(self):
        from_starts_and_ends = self.read_file('from_starts_and_ends')

        visualize.plot_metagene_positions(from_starts_and_ends['from_starts'],
                                          from_starts_and_ends['from_ends'],
                                          self.figure_file_names['starts_and_ends_zoomed_out'],
                                          zoomed_out=True,
                                         )
        
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


    def get_joint_position_counts(self, gene_name):
        CDSs, _ = self.get_CDSs()
        CDS_dict = {t.name: t for t in CDSs}
        transcript = CDS_dict[gene_name]

        joint_position_counts = positions.get_joint_position_counts_sparse(self.file_names['combined_extended_sorted'],
                                                                           transcript,
                                                                          )
        return joint_position_counts, transcript

    def get_total_eligible_reads(self):
        log_pairs = self.read_file('log')
        log_dict = {name: values[0] for name, values in log_pairs}
        total_mapped_reads = log_dict['Nonunique'] + log_dict['Concordant']
        return total_mapped_reads

if __name__ == '__main__':
    script_path = os.path.realpath(__file__)
    map_reduce.controller(TIFSeqExperiment, script_path)
