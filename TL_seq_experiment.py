import os
import pysam
from Sequencing.Parallel import map_reduce
from Sequencing import mapping_tools, fastq, sam, utilities
from Serialize import read_positions
import positions
import visualize
import rna_experiment

class TLSeqExperiment(rna_experiment.RNAExperiment):
    num_stages = 2

    specific_results_files = [
        ('sorted_by_name', 'bam', '{name}_by_name.bam'),
        ('bam', 'bam', '{name}.bam'),

        ('five_prime_read_positions', read_positions, '{name}_five_prime_read_positions.hdf5'),
    ]

    specific_figure_files = [
    ]

    specific_outputs = [
        ['bam',
        ],
        ['five_prime_read_positions',
         'from_starts_and_ends',
        ],
    ]

    specific_work = [
        ['preprocess',
         'map_tophat',
         'sort_by_name',
         'count_mappings',
        ],
        ['get_read_positions',
         'get_metagene_positions',
        ]
    ]

    specific_cleanup = [
        [],
        ['plot_starts_and_ends',
        ],
    ]

    def __init__(self, **kwargs):
        super(TLSeqExperiment, self).__init__(**kwargs)
        
    def map_tophat(self):
        mapping_tools.map_tophat([self.file_names['preprocessed_reads']],
                                 self.file_names['bowtie2_index_prefix'],
                                 self.file_names['genes'],
                                 self.file_names['transcriptome_index'],
                                 self.file_names['tophat_dir'],
                                )

    def get_read_positions(self):
        piece_CDSs, max_gene_length = self.get_CDSs()
        gene_infos = positions.get_Transcript_position_counts(self.merged_file_names['bam'],
                                                              piece_CDSs,
                                                              relevant_lengths=[],
                                                              left_buffer=500,
                                                              right_buffer=500,
                                                             )

        self.five_prime_read_positions = {name: info['position_counts']
                                          for name, info in gene_infos.iteritems()}

        self.write_file('five_prime_read_positions', self.five_prime_read_positions)
    
    def get_metagene_positions(self):
        piece_CDSs, max_gene_length = self.get_CDSs()
        read_positions = self.load_read_positions(modifier='five_prime')
        processed_read_positions = {}
        for name in read_positions:
            gene = {'five_prime': read_positions[name]['all']}
            processed_read_positions[name] = gene
        from_starts_and_ends = positions.compute_metagene_positions(read_positions, max_gene_length)

        self.write_file('from_starts_and_ends', from_starts_and_ends)
    
    def plot_starts_and_ends(self):
        from_starts_and_ends = self.read_file('from_starts_and_ends')

        visualize.plot_metagene_positions(from_starts_and_ends['from_starts'],
                                          from_starts_and_ends['from_ends'],
                                          self.figure_file_names['starts_and_ends_zoomed_out'],
                                          zoomed_out=True,
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

        os.rename(self.file_names['accepted_hits'], self.file_names['bam'])

    def count_mappings(self):
        num_unmapped = 0
        num_nonunique = 0
        num_unique = 0

        mappings = pysam.Samfile(self.file_names['sorted_by_name'])

        groups = utilities.group_by(mappings, lambda m: m.qname)
        for qname, group in groups:
            unmapped = any(m.is_unmapped for m in group)
            if unmapped:
                num_unmapped += 1
                continue

            nonunique = len(group) > 1 or any(m.mapq < 40 for m in group)
            if nonunique:
                num_nonunique += 1
                continue

            num_unique += 1
            
        self.log.extend(
            [('Unmapped', num_unmapped),
             ('Nonunique', num_nonunique),
             ('Unique', num_unique),
            ],
        )
    
    def get_total_eligible_reads(self):
        log_pairs = self.read_file('log')
        log_dict = {name: values[0] for name, values in log_pairs}
        total_mapped_reads = log_dict['Nonunique'] + log_dict['Unique']
        return total_mapped_reads

if __name__ == '__main__':
    script_path = os.path.realpath(__file__)
    map_reduce.controller(TLSeqExperiment, script_path)
