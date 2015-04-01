import glob
import numpy as np
from Sequencing.Parallel import map_reduce, split_file, piece_of_list
from Serialize import read_positions, lengths
from Sequencing import fastq
from itertools import chain
import gtf
import gff
import os
from collections import defaultdict
import logging

class RNAExperiment(map_reduce.MapReduceExperiment):
    specific_results_files = [ 
        ('preprocessed_reads', 'fastq', '{name}_preprocessed.fastq'),
        ('tophat_dir', 'dir', 'tophat'),
        ('accepted_hits', 'bam', 'tophat/accepted_hits.bam'),
        ('unmapped_bam', 'bam', 'tophat/unmapped.bam'),
        ('read_positions', read_positions, '{name}_read_positions.hdf5'),
        ('metagene_positions', read_positions, '{name}_metagene_positions.hdf5'),
        ('lengths', lengths, '{name}_lengths.hdf5'),
    ]

    specific_figure_files = [
        ('starts_and_ends', '{name}_starts_and_ends.pdf'),
        ('three_prime_starts_and_ends', '{name}_three_prime_starts_and_ends.pdf'),
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
        
        # As part of the process of determining gene annotation boundaries,
        # experiments need to be run on placeholder gene models.
        # This is controlled by setting bootstrap to a file name suffix.
        self.bootstrap = kwargs.get('bootstrap', '')
        
        self.organism_files = [
            ('bowtie2_index_prefix', 'genome/genome'),
            ('genome', 'genome'),
            ('genes', 'transcriptome/genes{0}.gff'.format(self.bootstrap)),
            ('transcriptome_index', 'transcriptome/bowtie2_index/genes{0}'.format(self.bootstrap)),
            ('rRNA_index', 'contaminant/bowtie2_index/rRNA'),
            ('oligos', 'contaminant/subtraction_oligos.fasta'),
            ('oligos_sam', 'contaminant/subtraction_oligos.sam'),
        ]

        self.data_fns = glob.glob(self.data_dir + '/*.fastq') + glob.glob(self.data_dir + '/*.fq')
        
        for key, tail in self.organism_files:
            self.file_names[key] = '{0}/{1}'.format(self.organism_dir, tail)
        
        self.min_length = 12
        self.max_read_length = kwargs.get('max_read_length', None)
        if self.max_read_length == None:
            self.max_read_length = self.get_max_read_length()
        else:
            self.max_read_length = int(self.max_read_length)
        
    def get_max_read_length(self):
        def length_from_file_name(file_name):
            length = len(fastq.reads(file_name).next().seq)
            return length
        
        if self.data_fns:
            max_length = max(length_from_file_name(fn) for fn in self.data_fns)
        else:
            max_length = 0

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
        ''' A generator over the reads in a piece of each data file.
            Can handle a mixture of different fastq encodings across (but not
            within) files.
        '''
        total_reads = 0
        for file_name in self.data_fns:
            total_reads_from_file = 0
            file_piece = split_file.piece(file_name,
                                          self.num_pieces,
                                          self.which_piece,
                                          'fastq',
                                         )
            for read in fastq.reads(file_piece, standardize_names=True, ensure_sanger_encoding=True):
                yield read
                
                total_reads += 1
                total_reads_from_file += 1
                if total_reads % 10000 == 0:
                    logging.info('{0:,} reads processed'.format(total_reads))

            head, tail = os.path.split(file_name)
            self.summary.append(('Reads in {0}'.format(tail), total_reads_from_file))

        logging.info('{0:,} total reads processed'.format(total_reads))
        
        self.summary.append(('Total reads', total_reads))
    
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

    def trim_reads(self, reads):
        trimmed_lengths = np.zeros(self.max_read_length + 1, int)
        too_short_lengths = np.zeros(self.max_read_length + 1, int)
    
        for trimmed_read in self.trim_function(reads):
            length = len(trimmed_read.seq)
            if length < self.min_length:
                too_short_lengths[length] += 1
            else:
                trimmed_lengths[length] += 1
                yield trimmed_read

        self.write_file('lengths', {'trimmed': trimmed_lengths,
                                    'too_short': too_short_lengths,
                                   },
                       )

    def get_CDSs(self, force_all=False):
        all_CDSs = gff.get_CDSs(self.file_names['genes'],
                                self.file_names['genome'],
                               )

        if self.transcripts_file_name == None:
            CDSs = all_CDSs
        else:
            transcripts = {line.strip() for line in open(self.transcripts_file_name)}
            CDSs = [t for t in all_CDSs if t.name in transcripts]
        
        max_gene_length = 0
        for CDS in CDSs:
            CDS.build_coordinate_maps()
            max_gene_length = max(max_gene_length, CDS.transcript_length)
            CDS.delete_coordinate_maps()
        
        if force_all:
            piece_CDSs = CDSs
        else:
            piece_CDSs = piece_of_list(CDSs, self.num_pieces, self.which_piece)

        return piece_CDSs, max_gene_length
