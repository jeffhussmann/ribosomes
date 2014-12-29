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
         'combine_mappings',
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
        
    def preprocess(self):                                                                                 
        ''' tophat can't handle named pipes, so need to make a file. '''                                                                                               
        reads = self.get_reads()                                                                          
        total_reads = 0                                                                                   
                                                                                                          
        with open(self.file_names['preprocessed_reads'], 'w') as preprocessed_file:                       
            for read in reads:                                                                            
                total_reads += 1                                                                          
                record = fastq.make_record(*read)                                                         
                preprocessed_file.write(record)                                                           
                                                                                                          
        self.summary.extend(                                                                                  
            [('Total reads', total_reads),                                                                
            ],                                                                                            
        )                                                                                                 
    
    def map_tophat(self):
        mapping_tools.map_tophat([self.file_names['preprocessed_reads']],
                                 self.file_names['bowtie2_index_prefix'],
                                 self.file_names['genes'],
                                 self.file_names['transcriptome_index'],
                                 self.file_names['tophat_dir'],
                                 no_sort=True,
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

    def combine_mappings(self):
        num_unmapped = 0
        num_nonunique = 0
        num_unique = 0

        mappings = pysam.Samfile(self.file_names['accepted_hits'])
        unmapped = pysam.Samfile(self.file_names['unmapped_bam'])
        merged = sam.merge_by_name(mappings, unmapped)
        grouped = utilities.group_by(merged, lambda m: m.qname)

        alignment_sorter = sam.AlignmentSorter(mappings.references,
                                               mappings.lengths,
                                               self.file_names['bam'],
                                              )
        with alignment_sorter:
            for qname, group in grouped:
                unmapped = any(m.is_unmapped for m in group)
                if unmapped:
                    num_unmapped += 1
                    continue

                nonunique = len(group) > 1 or any(m.mapq < 40 for m in group)
                if nonunique:
                    num_nonunique += 1
                else:
                    num_unique += 1

                for mapping in group:
                    alignment_sorter.write(mapping)
            
        self.summary.extend(
            [('Unmapped', num_unmapped),
             ('Nonunique', num_nonunique),
             ('Unique', num_unique),
            ],
        )
    
    def get_total_eligible_reads(self):
        summary_pairs = self.read_file('summary')
        summary_dict = {name: values[0] for name, values in summary_pairs}
        total_mapped_reads = summary_dict['Nonunique'] + summary_dict['Unique']
        return total_mapped_reads

if __name__ == '__main__':
    script_path = os.path.realpath(__file__)
    map_reduce.controller(TLSeqExperiment, script_path)
