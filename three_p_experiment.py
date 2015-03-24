import matplotlib
matplotlib.use('Agg', warn=False)
import os
import positions
import visualize
import trim
import pysam
from Sequencing import fastq, sam, mapping_tools, genomes, utilities
from Sequencing.Parallel import map_reduce
import Serialize.read_positions
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

        ('extended', 'bam', '{name}_extended.bam'),
        ('extended_filtered', 'bam', '{name}_extended_filtered.bam'),
    ]

    specific_figure_files = []

    specific_outputs = [
        ['trimmed_lengths',
         'too_short_lengths',
         'extended',
         'extended_filtered',
         'nongenomic_lengths',
        ],
        ['metagene_positions',
         'read_positions',
        ],
    ]

    specific_work = [
        ['preprocess',
         'map_tophat',
         'extend_mappings',
         'filter_mappings',
        ],
        ['get_polyA_positions',
         'get_metagene_positions',
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
        
        self.trim_function = trim.bound_trim['polyA']

    def preprocess(self):
        reads = self.get_reads()
        trimmed_reads = self.trim_reads(reads)

        with open(self.file_names['trimmed_reads'], 'w') as trimmed_fh:
            for read in trimmed_reads:
                trimmed_fh.write(str(read))
        
    def map_tophat(self):
        mapping_tools.map_tophat([self.file_names['trimmed_reads']],
                                 self.file_names['bowtie2_index_prefix'],
                                 self.file_names['genes'],
                                 self.file_names['transcriptome_index'],
                                 self.file_names['tophat_dir'],
                                 no_sort=True,
                                )

    def extend_mappings(self):
        trim.extend_polyA_ends(self.file_names['accepted_hits'],
                               self.file_names['extended'],
                               self.file_names['genome'],
                              )
        
    def filter_mappings(self):
        num_unmapped = 0
        num_entirely_genomic = 0
        num_nonunique = 0
        num_unique = 0

        nongenomic_lengths = Counter()

        sam_file = pysam.Samfile(self.file_names['accepted_hits'])
    
        region_fetcher = genomes.build_region_fetcher(self.file_names['genome'],
                                                      load_references=True,
                                                      sam_file=sam_file,
                                                     )

        extended_sorter = sam.AlignmentSorter(sam_file.references,
                                              sam_file.lengths,
                                              self.file_names['extended'],
                                             )
        filtered_sorter = sam.AlignmentSorter(sam_file.references,
                                              sam_file.lengths,
                                              self.file_names['extended_filtered'],
                                             )

        extended_mappings = (trim.extend_polyA_end(mapping, region_fetcher) for mapping in sam_file)
        mapping_groups = utilities.group_by(extended_mappings, lambda m: m.qname)

        with extended_sorter, filtered_sorter:
            for qname, group in mapping_groups:
                for m in group:
                    extended_sorter.write(m)

                min_nongenomic_length = min(trim.get_nongenomic_length(m) for m in group)
                nongenomic_lengths[min_nongenomic_length] += 1
                if min_nongenomic_length == 0:
                    num_entirely_genomic += 1
                    continue
                
                nonunique = len(group) > 1 or any(m.mapq < 40 for m in group)
                if nonunique:
                    num_nonunique += 1
                    continue
                
                num_unique += 1
                
                for m in group:
                    filtered_sorter.write(m)

        self.summary.extend(
            [('Mapped with no non-genomic A\'s', num_entirely_genomic),
             ('Nonunique', num_nonunique),
             ('Unique', num_unique),
            ],
        )

        nongenomic_lengths = utilities.counts_to_array(nongenomic_lengths)
        self.write_file('nongenomic_lengths', nongenomic_lengths)

    def get_polyA_positions(self):
        piece_CDSs, max_gene_length = self.get_CDSs()
        gene_infos = positions.get_Transcript_position_counts(self.merged_file_names['extended'],
                                                              piece_CDSs,
                                                              [],
                                                              left_buffer=500,
                                                              right_buffer=500,
                                                             )

        self.read_positions = {name: info['three_prime_positions']
                               for name, info in gene_infos.iteritems()}

        self.write_file('read_positions', self.read_positions)
        
    def get_metagene_positions(self):
        piece_CDSs, max_gene_length = self.get_CDSs()
        read_positions = self.load_read_positions()
        
        processed_read_positions = {}
        for name in read_positions:
            gene = {'three_prime_genomic': read_positions[name][0],
                    'three_prime_nongenomic': read_positions[name]['all'] - read_positions[name][0],
                    'three_prime_nonunique': three_prime_counts['all_nonunique'],
                    'sequence': read_positions[name]['sequence'],
                   }
            processed_read_positions[name] = gene
    
        metagene_positions = positions.compute_metagene_positions(piece_CDSs,
                                                                  processed_read_positions,
                                                                  max_gene_length,
                                                                 )

        self.write_file('metagene_positions', metagene_positions)
    
    def plot_starts_and_ends(self):
        metagene_positions = self.read_file('metagene_positions')

        visualize.plot_metagene_positions(metagene_positions,
                                          self.figure_file_names['starts_and_ends'],
                                          ['three_prime_genomic', 'three_prime_nongenomic'],
                                         )

    def get_total_eligible_reads(self):
        summary_pairs = self.read_file('summary')
        summary_dict = {name: values[0] for name, values in summary_pairs}
        total_mapped_reads = summary_dict['Nonunique'] + summary_dict['Unique']
        return total_mapped_reads

if __name__ == '__main__':
    script_path = os.path.realpath(__file__)
    map_reduce.controller(ThreePExperiment, script_path)
