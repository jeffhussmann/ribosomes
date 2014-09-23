import matplotlib
matplotlib.use('Agg', warn=False)
import os
import positions
import visualize
import trim
import pysam
from Sequencing import fastq, sam, mapping_tools, genomes, utilities
from Sequencing.Parallel import map_reduce
from Serialize import read_positions
from Sequencing.Serialize import array_1d
import rna_experiment
from Sequencing.annotation import Annotation_factory
from collections import Counter

class ThreePExperiment(rna_experiment.RNAExperiment):
    num_stages = 2
        
    specific_results_files = [
        ('trimmed_reads', 'fastq', '{name}_trimmed.fastq'),
        ('too_short_lengths', array_1d, '{name}_too_short_lengths.txt'),
        ('trimmed_lengths', array_1d, '{name}_trimmed_lengths.txt'),
        ('nongenomic_lengths', array_1d, '{name}_nongenomic_lengths.txt'),

        ('sorted_by_name', 'bam', '{name}_by_name.bam'),
        ('extended', 'bam', '{name}_extended.bam'),
        ('extended_sorted', 'bam', '{name}_extended_sorted.bam'),
        ('extended_filtered', 'bam', '{name}_extended_filtered.bam'),
        ('extended_filtered_sorted', 'bam', '{name}_extended_filtered_sorted.bam'),

        ('three_prime_read_positions', read_positions, '{name}_three_prime_read_positions.hdf5'),
    ]

    specific_figure_files = []

    specific_outputs = [
        ['trimmed_lengths',
         'too_short_lengths',
         #'extended_sorted',
         #'extended_filtered_sorted',
         'nongenomic_lengths',
        ],
        [#'from_starts_and_ends',
         #'three_prime_read_positions',
        ],
    ]

    specific_work = [
        [#'trim_reads',
         #'map_tophat',
         #'sort_by_name',
         #'extend_mappings',
         'filter_mappings',
        ],
        [#'get_polyA_positions',
         #'get_metagene_positions',
        ],
    ]

    specific_cleanup = [
        [],
        ['plot_starts_and_ends',
        ],
    ]

    def __init__(self, **kwargs):
        super(ThreePExperiment, self).__init__(**kwargs)

        self.max_read_length = self.get_max_read_length()
        self.min_length = 12

    def trim_reads(self):
        trim_function = trim.bound_trim['polyA']

        trimmed_lengths, too_short_lengths, barcode_counts = trim_function(self.get_reads(),
                                                                           self.file_names['trimmed_reads'],
                                                                           self.min_length,
                                                                           self.max_read_length,
                                                                          )
        self.write_file('trimmed_lengths', trimmed_lengths)
        self.write_file('too_short_lengths', too_short_lengths)
    
    def map_tophat(self):
        mapping_tools.map_tophat([self.file_names['trimmed_reads']],
                                 self.file_names['bowtie2_index_prefix'],
                                 self.file_names['genes'],
                                 self.file_names['transcriptome_index'],
                                 self.file_names['tophat_dir'],
                                )
    def sort_by_name(self):
        accepted_fn = self.file_names['accepted_hits']
        accepted_prefix, _ = os.path.splitext(accepted_fn)
        accepted_sorted_fn = '{}_sorted.bam'.format(accepted_prefix)

        unmapped_fn = self.file_names['unmapped_bam']
        unmapped_prefix, _ = os.path.splitext(unmapped_fn)
        unmapped_sorted_fn = '{}_sorted.bam'.format(unmapped_prefix)

        merged_fn = self.file_names['sorted_by_name']

        sam.sort_bam(accepted_fn, accepted_sorted_fn, by_name=True)
        sam.sort_bam(unmapped_fn, unmapped_sorted_fn, by_name=True)
        pysam.merge('-nf', merged_fn, accepted_sorted_fn, unmapped_sorted_fn)

    def extend_mappings(self):
        trim.extend_polyA_ends(self.file_names['sorted_by_name'],
                               self.file_names['extended'],
                               self.file_names['genome'],
                              )
        
        sam.sort_bam(self.file_names['extended'],
                     self.file_names['extended_sorted'],
                    )

    def filter_mappings(self):
        num_unmapped = 0
        num_entirely_genomic = 0
        num_nonunique = 0
        num_unique = 0

        nongenomic_lengths = Counter()

        extended_mappings = pysam.Samfile(self.file_names['extended'])
        groups = utilities.group_by(extended_mappings, lambda m: m.qname)

        with pysam.Samfile(self.file_names['extended_filtered'], 'wb', template=extended_mappings) as filtered_bam_file:
            for qname, group in groups:
                unmapped = any(m.is_unmapped for m in group)
                if unmapped:
                    num_unmapped += 1
                    continue

                min_nongenomic_length = min(trim.get_nongenomic_length(m) for m in group)
                nongenomic_lengths[min_nongenomic_length] += 1
                if min_nongenomic_length == 0:
                    num_entirely_genomic += 1
                    continue
                
                for m in group:
                    filtered_bam_file.write(m)

                nonunique = len(group) > 1 or any(m.mapq < 40 for m in group)
                if nonunique:
                    num_nonunique += 1
                    continue

                num_unique += 1
            
        self.log.extend(
            [('Unmapped', num_unmapped),
             ('Entirely genomic', num_entirely_genomic),
             ('Nonunique', num_nonunique),
             ('Unique', num_unique),
            ],
        )

        sam.sort_bam(self.file_names['extended_filtered'],
                     self.file_names['extended_filtered_sorted'],
                    )
        
        nongenomic_lengths = utilities.counts_to_array(nongenomic_lengths)
        self.write_file('nongenomic_lengths', nongenomic_lengths)

    def get_polyA_positions(self):
        piece_CDSs, max_gene_length = self.get_CDSs()
        gene_infos = positions.get_Transcript_polyA_position_counts(self.merged_file_names['extended_sorted'],
                                                                    piece_CDSs,
                                                                    max_relevant_length=5,
                                                                    left_buffer=500,
                                                                    right_buffer=500,
                                                                   )

        self.three_prime_read_positions = {name: info['position_counts']
                                           for name, info in gene_infos.iteritems()}

        self.write_file('three_prime_read_positions', self.three_prime_read_positions)
        
    def get_metagene_positions(self):
        piece_CDSs, max_gene_length = self.get_CDSs()
        read_positions = self.load_read_positions(modifier='three_prime')
        
        processed_read_positions = {}
        for name in read_positions:
            gene = {'three_prime_genomic': read_positions[name][0],
                    'three_prime_nongenomic': read_positions[name]['all'] - read_positions[name][0],
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

    def get_total_eligible_reads(self):
        log_pairs = self.read_file('log')
        log_dict = {name: values[0] for name, values in log_pairs}
        total_mapped_reads = log_dict['Nonunique'] + log_dict['Unique']
        return total_mapped_reads

if __name__ == '__main__':
    script_path = os.path.realpath(__file__)
    map_reduce.controller(ThreePExperiment, script_path)
