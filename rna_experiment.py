import glob
from Sequencing.Parallel import map_reduce, split_file, piece_of_list
from Serialize import read_positions
from Sequencing import fastq
from itertools import chain
import gtf
import os
from collections import defaultdict

class RNAExperiment(map_reduce.MapReduceExperiment):
    specific_results_files = [ 
        ('preprocessed_reads', 'fastq', '{name}_preprocessed.fastq'),
        ('tophat_dir', 'dir', 'tophat'),
        ('accepted_hits', 'bam', 'tophat/accepted_hits.bam'),
        ('unmapped_bam', 'bam', 'tophat/unmapped.bam'),
        ('read_positions', read_positions, '{name}_read_positions.hdf5'),
        ('from_starts_and_ends', read_positions, '{name}_from_starts_and_ends.hdf5'),
    ]

    specific_figure_files = [
        ('starts_and_ends', '{name}_starts_and_ends.pdf'),
        ('starts_and_ends_zoomed_out', '{name}_starts_and_ends_zoomed_out.pdf'),
    ]
    
    specific_outputs = []
    specific_work = []
    specific_cleanup = []
    
    def __init__(self, **kwargs):
        super(RNAExperiment, self).__init__(**kwargs)
        
        self.group = kwargs['group']
        self.data_dir = kwargs['data_dir'].rstrip('/')
        self.organism_dir = kwargs['organism_dir'].rstrip('/')
        self.transcripts_file_name = kwargs.get('transcripts_file_name', None)
        
        self.organism_files = [
            ('bowtie2_index_prefix', 'genome/genome'),
            ('genome', 'genome'),
            ('genes', 'transcriptome/genes.gtf'),
            ('transcriptome_index', 'transcriptome/bowtie2_index/genes'),
            ('rRNA_index', 'contaminant/bowtie2_index/rRNA'),
            ('oligos', 'contaminant/subtraction_oligos.fasta'),
            ('oligos_sam', 'contaminant/subtraction_oligos.sam'),
        ]
        
        self.data_fns = glob.glob(self.data_dir + '/*.fastq') + glob.glob(self.data_dir + '/*.fq')
        
        for key, tail in self.organism_files:
            self.file_names[key] = '{0}/{1}'.format(self.organism_dir, tail)

    def get_max_read_length(self):
        def length_from_file_name(file_name):
            length = len(fastq.reads(file_name).next().seq)
            return length
        
        max_length = max(length_from_file_name(fn) for fn in self.data_fns)
        return max_length

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
            read_pairs = fastq.read_pairs(R1_lines, R2_lines, standardize_names=True, ensure_sanger_encoding=True)
            read_pairs_list.append(read_pairs)
        
        all_read_pairs = chain.from_iterable(read_pairs_list)

        return all_read_pairs

    def preprocess(self):
        ''' Needs to exist to make files to feed to mapping programs. Would be
            better off as named pipes.
        '''
        reads = self.get_reads()
        total_reads = 0

        with open(self.file_names['preprocessed_reads'], 'w') as preprocessed_file:
            for read in reads:
                total_reads += 1
                record = fastq.make_record(*read)
                preprocessed_file.write(record)

        self.log.extend(
            [('Total reads', total_reads),
            ],
        )
    
    def get_CDSs(self):
        if self.transcripts_file_name == None:
            CDSs = gtf.get_CDSs(self.file_names['genes'])
        else:
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
