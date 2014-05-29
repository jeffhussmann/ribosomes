import matplotlib
matplotlib.use('Agg', warn=False)
import os
import glob
import ribosomes
import gtf
import positions
import visualize
from itertools import chain
import Sequencing.fastq as fastq
import Sequencing.sam as sam
from Sequencing.Parallel import map_reduce, split_file, piece_of_list

class ThreePExperiment(map_reduce.MapReduceExperiment):
    num_stages = 2

    def __init__(self, **kwargs):
        map_reduce.MapReduceExperiment.__init__(self, **kwargs)

        self.data_dir = kwargs['data_dir'].rstrip('/')
        self.organism_dir = kwargs['organism_dir'].rstrip('/')
        self.transcripts_file_name = kwargs.get('transcripts_file_name')

        specific_results_files = [
            ('preprocessed', 'fastq', '{name}_preprocessed.fastq'),

            ('read_positions', 'read_positions', '{name}_read_positions.txt'),
            ('read_positions_remapped', 'read_positions', '{name}_read_positions_remapped.txt'),
            ('from_starts_and_ends', 'read_positions', '{name}_from_starts_and_ends.txt'),
            ('from_starts_and_ends_remapped', 'read_positions', '{name}_from_starts_and_ends_remapped.txt'),

            ('tophat_dir', 'dir', 'tophat'),
            ('accepted_hits', 'bam', 'tophat/accepted_hits.bam'),
            ('clean_bam', 'bam', '{name}_clean.bam'),
            ('unmapped_bam', 'bam', 'tophat/unmapped.bam'),
            ('unmapped_trimmed_fastq', 'fastq', '{name}_unmapped_trimmed.fastq'),
            ('tophat_remapped_polyA_dir', 'dir', 'tophat_remapped_polyA'),
            ('remapped_accepted_hits', 'bam', 'tophat_remapped_polyA/accepted_hits.bam'),
            ('remapped_clean_bam', 'bam', '{name}_remapped_clean.bam'),
        
            ('yield', '', '{name}_yield.txt'),
        ]


        specific_figure_files = [
            ('starts_and_ends', '{name}_starts_and_ends.pdf'),
            ('starts_and_ends_zoomed_out', '{name}_starts_and_ends_zoomed_out.pdf'),
            ('starts_and_ends_remapped', '{name}_starts_and_ends_remapped.pdf'),
            ('starts_and_ends_remapped_zoomed_out', '{name}_starts_and_ends_remapped_zoomed_out.pdf'),
        ]

        self.organism_files = [
            ('bowtie2_index_prefix', 'genome/genome'),
            ('genome', 'genome'),
            ('genes', 'transcriptome/genes.gtf'),
            ('transcriptome_index', 'transcriptome/bowtie2_index/genes'),
        ]

        specific_outputs = [
            ['clean_bam',
             'remapped_clean_bam',
            ],
            ['from_starts_and_ends',
             'from_starts_and_ends_remapped',
            ],
        ]

        specific_work = [
            [(self.preprocess, 'Preprocessing'),
             (self.map_tophat, 'Mapping'),
             (self.remap_trimmed, 'Remapping trimmed'),
            ],
            [(self.get_read_positions, 'Counting mapping positions'),
             (self.get_metagene_positions, 'Aggregating metagene positions'),
            ],
        ]

        specific_cleanup = [
            [self.index_bams,
            ],
            [self.plot_starts_and_ends,
            ],
        ]

        self.results_files.extend(specific_results_files)
        self.figure_files.extend(specific_figure_files)
        map_reduce.extend_stages(self.outputs, specific_outputs)
        map_reduce.extend_stages(self.work, specific_work)
        map_reduce.extend_stages(self.cleanup, specific_cleanup)

        self.make_file_names()
        
        self.data_fns = glob.glob(self.data_dir + '/*.fastq') + glob.glob(self.data_dir + '/*.fq')

        self.max_read_length = self.get_max_read_length()

        for key, tail in self.organism_files:
            self.file_names[key] = '{0}/{1}'.format(self.organism_dir, tail)
    
    def get_reads(self):
        ''' Returns a generator over the reads in a piece of each data file.
            Can handle a mixture of different fastq encodings across (but not
            within) files.
        '''
        file_pieces = [split_file.piece(file_name,
                                        self.num_pieces,
                                        self.which_piece,
                                        'fastq',
                                       )
                       for file_name in self.data_fns]
        read_pieces = [fastq.reads(piece, standardize_names=True, ensure_sanger_encoding=True)
                       for piece in file_pieces]
        reads = chain.from_iterable(read_pieces)
        return reads

    def preprocess(self):
        reads = self.get_reads()
        with open(self.file_names['preprocessed'], 'w') as preprocessed_file:
            for read in reads:
                record = fastq.make_record(*read)
                preprocessed_file.write(record)
    
    def get_max_read_length(self):
        def length_from_file_name(file_name):
            length = len(fastq.reads(file_name).next().seq)
            return length
        
        max_length = max(length_from_file_name(fn) for fn in self.data_fns)
        return max_length

    def map_tophat(self):
        ribosomes.map_tophat([self.file_names['preprocessed']],
                             self.file_names['bowtie2_index_prefix'],
                             self.file_names['genes'],
                             self.file_names['transcriptome_index'],
                             self.file_names['tophat_dir'],
                            )
        os.rename(self.file_names['accepted_hits'],
                  self.file_names['clean_bam'],
                 )
    
    def remap_trimmed(self):
        ribosomes.trim_polyA_from_unmapped(self.file_names['unmapped_bam'],
                                          self.file_names['unmapped_trimmed_fastq'],
                                          15,
                                          self.max_read_length,
                                         )
        ribosomes.map_tophat([self.file_names['unmapped_trimmed_fastq']],
                             self.file_names['bowtie2_index_prefix'],
                             self.file_names['genes'],
                             self.file_names['transcriptome_index'],
                             self.file_names['tophat_remapped_polyA_dir'],
                            )
        os.rename(self.file_names['remapped_accepted_hits'],
                  self.file_names['remapped_clean_bam'],
                 )
    
    def index_bams(self):
        sam.index_bam(self.merged_file_names['clean_bam'])
        sam.index_bam(self.merged_file_names['remapped_clean_bam'])
    
    def get_CDSs(self):
        all_CDSs = gtf.get_CDSs(self.file_names['genes'])
        transcripts = {line.strip() for line in open(self.transcripts_file_name)}
        CDSs = [t for t in all_CDSs if t.name in transcripts]
        
        max_gene_length = 0
        for CDS in CDSs:
            CDS.build_coordinate_maps()
            max_gene_length = max(max_gene_length, CDS.CDS_length)
            CDS.delete_coordinate_maps()
        
        piece_CDSs = piece_of_list(CDSs, self.num_pieces, self.which_piece)
        return piece_CDSs, max_gene_length

    def get_read_positions(self):
        piece_CDSs, max_gene_length = self.get_CDSs()
        gene_infos = positions.get_Transcript_position_counts(self.merged_file_names['clean_bam'],
                                                              piece_CDSs,
                                                              relevant_lengths=[],
                                                             )

        self.read_positions = {name: info['position_counts']
                               for name, info in gene_infos.iteritems()}

        self.write_file('read_positions', self.read_positions)
        
        remapped_gene_infos = positions.get_Transcript_position_counts(self.merged_file_names['remapped_clean_bam'],
                                                                       piece_CDSs,
                                                                       relevant_lengths=[],
                                                                      )

        self.remapped_read_positions = {name: info['position_counts']
                                        for name, info in remapped_gene_infos.iteritems()}

    def get_metagene_positions(self):
        piece_CDSs, max_gene_length = self.get_CDSs()
        read_positions = self.load_read_positions()
        from_starts_and_ends = positions.compute_metagene_positions(read_positions, max_gene_length)
        remapped_read_positions = self.load_read_positions(modifier='remapped')
        from_starts_and_ends_remapped = positions.compute_metagene_positions(remapped_read_positions, max_gene_length)

        self.write_file('from_starts_and_ends', from_starts_and_ends)
        self.write_file('from_starts_and_ends_remapped', from_starts_and_ends_remapped)
    
    def load_read_positions(self, modifier=None):
        ''' Read position arrays can be expensive to read from disk. This is a
            a wrapper to make sure this is only done once.
        '''
        if modifier != None:
            attribute_name = '{0}_read_positions'.format(modifier)
        else:
            attribute_name = 'read_positions'

        if hasattr(self, attribute_name):
            pass
        else:
            setattr(self, attribute_name, self.read_file(attribute_name))

        return getattr(self, attribute_name)

    def plot_starts_and_ends(self):
        from_starts_and_ends = self.read_file('from_starts_and_ends')

        visualize.plot_metagene_positions(from_starts_and_ends['from_starts'],
                                          from_starts_and_ends['from_ends'],
                                          self.figure_file_names['starts_and_ends'],
                                         )
        
        visualize.plot_metagene_positions(from_starts_and_ends['from_starts'],
                                          from_starts_and_ends['from_ends'],
                                          self.figure_file_names['starts_and_ends_zoomed_out'],
                                          zoomed_out=True,
                                         )
        
        from_starts_and_ends_remapped = self.read_file('from_starts_and_ends_remapped')

        visualize.plot_metagene_positions(from_starts_and_ends_remapped['from_starts'],
                                          from_starts_and_ends_remapped['from_ends'],
                                          self.figure_file_names['starts_and_ends_remapped'],
                                         )
        
        visualize.plot_metagene_positions(from_starts_and_ends_remapped['from_starts'],
                                          from_starts_and_ends_remapped['from_ends'],
                                          self.figure_file_names['starts_and_ends_remapped_zoomed_out'],
                                          zoomed_out=True,
                                         )

if __name__ == '__main__':
    script_path = os.path.realpath(__file__)
    map_reduce.controller(ThreePExperiment, script_path)
