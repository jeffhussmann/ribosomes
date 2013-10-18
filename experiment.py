import matplotlib
matplotlib.use('Agg', warn=False)
import glob
import trim
import os
import ribosomes
import fastq
from itertools import chain
from collections import Counter
import mutations
import pysam
import sam
import numpy as np
import mapreduce
import Parallel.split_file

class RibosomeProfilingExperiment(mapreduce.MapReduceExperiment):
    def __init__(self, **kwargs):
        self.data_dir = kwargs['data_dir'].rstrip('/')
        self.adapter_type = kwargs['adapter_type']
        self.organism_dir = kwargs['organism_dir'].rstrip('/')
        
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

            ('rRNA_coverage', 'coverage', '{name}_rRNA_coverage.txt'),
            ('frames', 'frames', '{name}_frames.txt'),

            ('tophat_dir', 'dir', 'tophat'),
            ('accepted_hits', 'bam', 'tophat/accepted_hits.bam'),
            ('unmapped_bam', 'bam', 'tophat/unmapped.bam'),

            ('log', '', '{name}_log.txt'),
            ('rRNA_coverage_figure_template', '', '{name}_rRNA_coverage_{{0}}.pdf'),
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
                         'rRNA_coverage',
                        ),
                       ]

        self.work = [[(self.trim_reads, 'Trim reads'),
                      (self.pre_filter_rRNA, 'Filter contaminants'),
                      (self.tophat, 'Map with tophat'),
                      (self.post_filter_contaminants, 'Further contaminant filtering'),
                      (self.get_rRNA_coverage, 'Counting rRNA coverage'),
                     ],
                    ]

        self.cleanup = [(self.make_log,
                         self.plot_rRNA_coverage,
                        ),
                       ]
        
        mapreduce.MapReduceExperiment.__init__(self, **kwargs)

        self.data_fns = glob.glob(self.data_dir + '/*.fastq')
        self.max_read_length = self.get_max_read_length()
        
        if self.adapter_type == 'truseq':
            self.trim_function = trim.trim_adapters
        elif self.adapter_type == 'polyA':
            self.trim_function = trim.trim_poly_A
        elif self.adapter_type == 'nothing':
            self.trim_function = trim.trim_nothing

        for key, tail in self.organism_files:
            self.file_names[key] = '{0}/{1}'.format(self.organism_dir, tail)
        
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

    def get_rRNA_coverage(self):
        bam_file_names = [self.file_names['rRNA_bam'],
                          #self.file_names['more_rRNA_bam'],
                         ]
        data = ribosomes.produce_rRNA_coverage(bam_file_names)
        self.write_file('rRNA_coverage', data)
    
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

    def plot_rRNA_coverage(self):
        ribosomes.plot_rRNA_coverage([self.name],
                                     [self.read_file('rRNA_coverage')],
                                     self.file_names['oligos_sam'],
                                     self.file_names['rRNA_coverage_figure_template'],
                                    )

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

if __name__ == '__main__':
    script_path = os.path.realpath(__file__)
    mapreduce.controller(RibosomeProfilingExperiment, script_path)
