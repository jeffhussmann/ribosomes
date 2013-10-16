import matplotlib
matplotlib.use('Agg', warn=False)
import argparse
import glob
import time
import trim
import os
import Serialize
import ribosomes
import subprocess
import fastq
from itertools import chain
from collections import Counter
import Parallel.split_file
import mutations
import pysam
import sam
import numpy as np
import Serialize

class RibosomeProfilingExperiment(object):
    def __init__(self, **kwargs):
        self.name = kwargs['name']
        self.work_prefix = kwargs['work_prefix']
        self.scratch_prefix = kwargs['scratch_prefix']
        self.data_dir = kwargs['data_dir']
        self.organism_dir = kwargs['organism_dir']
        self.relative_results_dir = kwargs['relative_results_dir']
        self.adapter_type = kwargs['adapter_type']
        self.num_pieces = kwargs['num_pieces']
        self.which_piece = kwargs['which_piece']
        
        suffix = Parallel.split_file.generate_suffix(self.num_pieces, self.which_piece)

        self.scratch_results_dir = '{0}/{1}{2}'.format(self.scratch_prefix,
                                                       self.relative_results_dir,
                                                       suffix,
                                                      )
        self.work_results_dir = '{0}/{1}'.format(self.work_prefix,
                                                 self.relative_results_dir,
                                                )
        
        if not os.path.isdir(self.scratch_results_dir):
            os.makedirs(self.scratch_results_dir)
        if self.which_piece == -1:
            # Sentinel value that indicates this is the merged experiment.
            # Only create the directory in this instance to avoid race
            # conditions.
            if not os.path.isdir(self.work_results_dir):
                os.makedirs(self.work_results_dir)

        self.results_files = [
            ('trimmed_reads', 'fastq', '{name}_trimmed.fastq'),
            ('filtered_reads', 'fastq', '{name}_filtered.fastq'), 

            ('trimmed_lengths', 'array', '{name}_trimmed_lengths.txt'),
            ('filtered_lengths', 'array', '{name}_filtered_lengths.txt'),
            ('tRNA_lengths', 'array', '{name}_tRNA_lengths.txt'),
            ('rRNA_lengths', 'array', '{name}_rRNA_lengths.txt'),
            ('clean_lengths', 'array', '{name}_clean_lengths.txt'),
            ('unmapped_lengths', 'array', '{name}_unmapped_lengths.txt'),
            ('unambiguous_lengths', 'array', '{name}_unambiguous_lengths.txt'),

            ('rRNA_sam', 'sam', '{name}_rRNA.sam'),
            ('rRNA_bam', 'bam', '{name}_rRNA.bam'),
            ('clean_bam', 'bam', '{name}_clean.bam'),
            ('more_rRNA_bam', 'bam', '{name}_more_rRNA.bam'),
            ('tRNA_bam', 'bam', '{name}_tRNA.bam'),
            ('unambiguous_bam', 'bam', '{name}_unambiguous.bam'),

            ('from_starts', 'array', '{name}_from_starts.txt'),
            ('from_ends', 'array', '{name}_from_ends.txt'),
            ('from_starts_unambiguous', 'array', '{name}_from_starts_unambiguous.txt'),
            ('from_ends_unambiguous', 'array', '{name}_from_ends_unambiguous.txt'),
            ('binned_positions', 'array', '{name}_binned_positions.txt'),

            ('frames', 'frames', '{name}_frames.txt'),

            ('tophat_dir', 'dir', 'tophat'),
            ('accepted_hits', 'bam', 'tophat/accepted_hits.bam'),
            ('unmapped_bam', 'bam', 'tophat/unmapped.bam'),

            ('log', '', '{name}_log.txt'),
        ]

        self.organism_files = [
            ('index', 'genome/genome'),
            ('genome', 'genome/genome.fa'),
            ('genes', 'transcriptome/genes.gtf'),
            ('transcriptome_index', 'transcriptome/bowtie2_index/genes'),
            ('rRNA_index', 'contaminant/bowtie2_index/rRNA'),
            ('oligos', 'contaminant/subtraction_oligos.fastq'),
            ('oligos_sam', 'contaminant/subtraction_oligos.sam'),
        ]

        self.outputs = [('trimmed_lengths',
                         'filtered_lengths',
                         'tRNA_lengths',
                         'rRNA_lengths',
                         'clean_lengths',
                         'unmapped_lengths',
                        ),
                       ]

        self.work = [[(self.trim_reads, 'Trim reads'),
                      (self.pre_filter_rRNA, 'Filter contaminants'),
                      (self.tophat, 'Map with tophat'),
                      (self.post_filter_contaminants, 'Further contaminant filtering'),
                     ],
                    ]

        self.cleanup = [(self.make_log,
                        ),
                       ]

        self.data_fns = glob.glob(self.data_dir + '/*.fastq')
        self.max_read_length = self.get_max_read_length()
        
        if self.adapter_type == 'truseq':
            self.trim_function = trim.trim_adapters
        elif self.adapter_type == 'polyA':
            self.trim_function = trim.trim_poly_A
        elif self.adapter_type == 'nothing':
            self.trim_function = trim.trim_nothing

        num_stages = len(self.outputs)
        timing_names = [('timing_{0}'.format(i), 'log', '{{name}}_timing_{0}.txt'.format(i))
                        for i in range(num_stages)]
        self.results_files.extend(timing_names)
        
        self.file_names = {}
        self.merged_file_names = {}
        self.file_types = {}
        for key, serialize_type, tail_template in self.results_files:
            file_tail = tail_template.format(name=self.name)
            self.file_names[key] = '{0}/{1}'.format(self.scratch_results_dir, file_tail)
            self.merged_file_names[key] = '{0}/{1}'.format(self.work_results_dir, file_tail)
            self.file_types[key] = serialize_type

        if self.which_piece == -1:
            self.file_names = self.merged_file_names

        for key, tail in self.organism_files:
            self.file_names[key] = '{0}/{1}'.format(self.organism_dir, tail)

    def write_file(self, key, data):
        Serialize.write_file(data, self.file_names[key], self.file_types[key])

    def read_file(self, key):
        data = Serialize.read_file(self.file_names[key], self.file_types[key])
        return data

    def do_work(self, stage):
        times = []
        for function, description in self.work[stage]:
            start_time = time.time()
            function()
            end_time = time.time()
            times.append((description, end_time - start_time))

        self.write_file('timing_{0}'.format(stage), times)

    def do_cleanup(self, stage):
        for function in self.cleanup[stage]:
            function()

    def trim_reads(self):
        trimmed_lengths = self.trim_function(self.get_reads(),
                                             self.max_read_length,
                                             self.file_names['trimmed_reads'],
                                             self.file_names['trimmed_lengths'],
                                            )
        self.write_file('trimmed_lengths', trimmed_lengths)

    def pre_filter_rRNA(self):
        ribosomes.pre_filter(self.file_names['rRNA_index'],
                             self.file_names['trimmed_reads'],
                             self.file_names['filtered_reads'],
                             self.file_names['rRNA_sam'],
                             self.file_names['rRNA_bam'],
                            )
        
        filtered_reads = fastq.reads(self.file_names['filtered_reads'])
        filtered_lengths = Counter(len(read.seq) for read in filtered_reads)
        filtered_lengths = self.zero_padded_array(filtered_lengths)
        self.write_file('filtered_lengths', filtered_lengths)

    def tophat(self):
        ribosomes.map_tophat(self.file_names['filtered_reads'],
                             self.file_names['index'],
                             self.file_names['genes'],
                             self.file_names['transcriptome_index'],
                             self.file_names['tophat_dir'],
                            )
        pysam.index(self.file_names['accepted_hits'])
    
    def post_filter_contaminants(self):
        ribosomes.filter_fetch(self.file_names['accepted_hits'],
                               self.file_names['genes'],
                               self.file_names['clean_bam'],
                               self.file_names['more_rRNA_bam'],
                               self.file_names['tRNA_bam'],
                              )
        
        tRNA_length_counts = sam.get_length_counts(self.file_names['tRNA_bam'])
        tRNA_lengths = self.zero_padded_array(tRNA_length_counts)
        self.write_file('tRNA_lengths', tRNA_lengths)

        # Anything that was in trimmed_reads and didn't make it to
        # filtered_reads was an rRNA read.
        rRNA_length_counts = sam.get_length_counts(self.file_names['more_rRNA_bam'])
        rRNA_lengths = self.zero_padded_array(rRNA_length_counts)
        trimmed_lengths = self.read_file('trimmed_lengths')
        filtered_lengths = self.read_file('filtered_lengths')
        rRNA_lengths += trimmed_lengths - filtered_lengths
        self.write_file('rRNA_lengths', rRNA_lengths)
        
        unmapped_length_counts = sam.get_length_counts(self.file_names['unmapped_bam'])
        unmapped_lengths = self.zero_padded_array(unmapped_length_counts)
        self.write_file('unmapped_lengths', unmapped_lengths)

        clean_length_counts = sam.get_length_counts(self.file_names['clean_bam'])
        clean_lengths = self.zero_padded_array(clean_length_counts)
        self.write_file('clean_lengths', clean_lengths)

    def make_log(self):
        trimmed_lengths = self.read_file('trimmed_lengths')
        filtered_lengths = self.read_file('filtered_lengths')
        tRNA_lengths = self.read_file('tRNA_lengths')
        rRNA_lengths = self.read_file('rRNA_lengths')
        unmapped_lengths = self.read_file('unmapped_lengths')
        clean_lengths = self.read_file('clean_lengths')
        
        total_reads = int(trimmed_lengths.sum())
        rRNA_reads = int(rRNA_lengths.sum())
        tRNA_reads = int(tRNA_lengths.sum())
        unmapped_reads = int(unmapped_lengths.sum())
        clean_reads = int(clean_lengths.sum())

        with open(self.file_names['log'], 'w') as log_file:
            log_file.write('Total reads: {0:,}\n'.format(total_reads))
            for category, count in [('rRNA reads', rRNA_reads),
                                    ('tRNA reads', tRNA_reads),
                                    ('Unmapped reads', unmapped_reads),
                                    ('Clean reads', clean_reads),
                                   ]:
                fraction = float(count) / total_reads
                line = '{0}: {1:,} ({2:.2%})\n'.format(category,
                                                         count,
                                                         fraction,
                                                        )
                log_file.write(line)

    def get_max_read_length(self):
        def length_from_file_name(file_name):
            length = len(fastq.reads(file_name).next().seq)
            return length
        
        max_length = max(length_from_file_name(fn) for fn in self.data_fns)
        return max_length

    def zero_padded_array(self, counts):
        array = mutations.counts_to_array(counts)
        pad_length = self.max_read_length + 1 - len(array)
        if pad_length > 0:
            padding = np.zeros(pad_length, int)
            array = np.append(array, padding)
        return array

    def get_reads(self):
        file_pieces = [Parallel.split_file.piece(file_name,
                                                 self.num_pieces,
                                                 self.which_piece,
                                                 'fastq',
                                                )
                       for file_name in self.data_fns]
        lines = chain.from_iterable(file_pieces)
        reads = fastq.reads(lines)
        return reads
    
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--job_dir',
                        required=True,
                       )

    subparsers = parser.add_subparsers(dest='subparser_name')

    parser_process = subparsers.add_parser('process')
    parser_process.add_argument('--num_pieces', type=int,
                                help='how many pieces',
                                default=1,
                               )
    parser_process.add_argument('--which_piece', type=int,
                                help='which piece this is',
                                default=0,
                               )
    parser_process.add_argument('--stage', type=int,
                                help='stage',
                                default=0,
                               )
    
    parser_launch = subparsers.add_parser('launch')
    parser_launch.add_argument('--num_pieces', type=int,
                               default=1,
                              )

    parser_finish = subparsers.add_parser('finish')
    parser_finish.add_argument('--num_pieces', type=int,
                                help='how many pieces',
                                default=1,
                               )
    parser_finish.add_argument('--stage', type=int,
                                help='stage',
                                default=0,
                               )

    args = parser.parse_args()
    return args

def parse_experiment_description(description_fn):
    description = dict(line.strip().split() for line in open(description_fn))
    return description

def launch(args):
    if not os.path.isdir(args.job_dir):
        os.makedirs(args.job_dir)
    
    description_file_name = '{0}/description.txt'.format(args.job_dir)
    description = parse_experiment_description(description_file_name)

    script_path = os.path.realpath(__file__)

    def make_process_command(args, which_piece, stage):
        command = ['python', script_path,
                   '--job_dir', args.job_dir,
                   'process',
                   '--num_pieces', str(args.num_pieces),
                   '--which_piece', str(which_piece),
                   '--stage', str(stage),
                  ]
        command_string = ' '.join(command) + '\n'
        return command_string
    
    def make_finish_command(args, stage):
        command = ['python', script_path,
                   '--job_dir', args.job_dir,
                   'finish',
                   '--num_pieces', str(args.num_pieces),
                   '--stage', str(stage),
                  ]
        command_string = ' '.join(command) + '\n'
        return command_string

    for stage in range(1):
        process_file_name = '{0}/process_{1}_stage_{2}'.format(args.job_dir,
                                                               args.num_pieces,
                                                               stage,
                                                              )
        finish_file_name = '{0}/finish_{1}_stage_{2}'.format(args.job_dir,
                                                             args.num_pieces,
                                                             stage,
                                                            )
                              
        with open(process_file_name, 'w') as process_file:
            for which_piece in range(args.num_pieces):
                line = make_process_command(args, which_piece, stage)
                process_file.write(line)

        with open(finish_file_name, 'w') as finish_file:
            line = make_finish_command(args, stage)
            finish_file.write(line)
    
        print 'Launched {0} with parallel'.format(stage)
        subprocess.check_call('parallel < {0}'.format(process_file_name), shell=True)
        subprocess.check_call('bash {0}'.format(finish_file_name), shell=True)

def process(args):
    description_file_name = '{0}/description.txt'.format(args.job_dir)
    description = parse_experiment_description(description_file_name)

    experiment = RibosomeProfilingExperiment(num_pieces=args.num_pieces,
                                             which_piece=args.which_piece,
                                             **description)
    experiment.do_work(args.stage)

def finish(args):
    description_file_name = '{0}/description.txt'.format(args.job_dir)
    description = parse_experiment_description(description_file_name)
    
    merged = RibosomeProfilingExperiment(num_pieces=args.num_pieces,
                                         which_piece=-1,
                                         **description)
    pieces = [RibosomeProfilingExperiment(num_pieces=args.num_pieces,
                                          which_piece=which_piece,
                                          **description)
              for which_piece in range(args.num_pieces)]

    for key in merged.outputs[args.stage]:
        piece_file_names = [piece.file_names[key] for piece in pieces]
        merged_file_name = merged.merged_file_names[key]
        Serialize.merge_files(piece_file_names,
                             merged_file_name,
                             merged.file_types[key],
                            )

    merged.do_cleanup(args.stage)

if __name__ == '__main__':
    args = parse_arguments()
    if args.subparser_name == 'launch':
        launch(args)
    elif args.subparser_name == 'process':
        process(args)
    elif args.subparser_name == 'finish':
        finish(args)
