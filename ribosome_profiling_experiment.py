import matplotlib
matplotlib.use('Agg', warn=False)
import matplotlib.pyplot as plt
import glob
import trim
import os
import ribosomes
import contaminants
import fastq
from itertools import chain
from collections import Counter
import mutations
import pysam
import sam
import numpy as np
import mapreduce
import Parallel.split_file
import gtf
import Visualize.mismatches
import mapping

class RibosomeProfilingExperiment(mapreduce.MapReduceExperiment):
    num_stages = 2

    def __init__(self, **kwargs):
        mapreduce.MapReduceExperiment.__init__(self, **kwargs)

        self.data_dir = kwargs['data_dir'].rstrip('/')
        self.adapter_type = kwargs['adapter_type']
        self.organism_dir = kwargs['organism_dir'].rstrip('/')
        self.fastq_type = kwargs.get('fastq_type', 'sanger')
        self.max_read_length = kwargs.get('max_read_length', None)
        
        self.min_length = 10
        
        specific_results_files = [
            ('bowtie_error', '', '{name}_bowtie_error.txt'),
            ('trimmed_reads', 'fastq', '{name}_trimmed.fastq'),
            ('filtered_reads', 'fastq', '{name}_filtered.fastq'), 

            ('too_short_lengths', 'array', '{name}_too_short_lengths.txt'),
            ('trimmed_lengths', 'array', '{name}_trimmed_lengths.txt'),
            ('filtered_lengths', 'array', '{name}_filtered_lengths.txt'),
            ('tRNA_lengths', 'array', '{name}_tRNA_lengths.txt'),
            ('rRNA_lengths', 'array', '{name}_rRNA_lengths.txt'),
            ('clean_lengths', 'array', '{name}_clean_lengths.txt'),
            ('unmapped_lengths', 'array', '{name}_unmapped_lengths.txt'),
            ('unambiguous_lengths', 'array', '{name}_unambiguous_lengths.txt'),
            ('oligo_hit_lengths', 'array', '{name}_oligo_hit_lengths.txt'),

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
            ('rpf_positions', 'rpf_positions', '{name}_rpf_positions.txt'),
            ('expression', 'expression', '{name}_expression.txt'),

            ('rRNA_coverage', 'coverage', '{name}_rRNA_coverage.txt'),
            ('frames', 'frames', '{name}_frames.txt'),

            ('tophat_dir', 'dir', 'tophat'),
            ('accepted_hits', 'bam', 'tophat/accepted_hits.bam'),
            ('unmapped_bam', 'bam', 'tophat/unmapped.bam'),

            ('mismatches_25', 'mismatches', '{name}_mismatches_25.txt'),
            ('mismatches_26', 'mismatches', '{name}_mismatches_26.txt'),
            ('mismatches_27', 'mismatches', '{name}_mismatches_27.txt'),
            ('mismatches_28', 'mismatches', '{name}_mismatches_28.txt'),
            ('mismatches_29', 'mismatches', '{name}_mismatches_29.txt'),
            ('mismatches_30', 'mismatches', '{name}_mismatches_30.txt'),
            ('mismatches_31', 'mismatches', '{name}_mismatches_31.txt'),

            ('recycling_ratios', 'ratios', '{name}_recycling_ratios.txt'),

            ('yield', '', '{name}_yield.txt'),
        ]

        specific_figure_files = [
            ('all_lengths', '{name}_all_lengths.pdf'),
            ('clean_lengths', '{name}_clean_lengths.pdf'),
            ('rRNA_coverage_template', '{name}_rRNA_coverage_{{0}}.pdf'),
            ('oligo_hit_lengths', '{name}_oligo_hit_lengths.pdf'),
            ('starts_and_ends', '{name}_starts_and_ends.png'),
            ('mismatch_positions', '{name}_mismatch_positions.png'),
            ('first_mismatch_types', '{name}_first_mismatch_types.png'),
            ('last_mismatch_types', '{name}_last_mismatch_types.png'),
        ]

        self.organism_files = [
            ('index', 'genome/genome'),
            ('genome', 'genome'),
            ('genes', 'transcriptome/genes.gtf'),
            ('transcriptome_index', 'transcriptome/bowtie2_index/genes'),
            ('rRNA_index', 'contaminant/bowtie2_index/rRNA'),
            ('oligos', 'contaminant/subtraction_oligos.fasta'),
            ('oligos_sam', 'contaminant/subtraction_oligos.sam'),
        ]

        specific_outputs = [
            ['too_short_lengths',
             'trimmed_lengths',
             'filtered_lengths',
             'tRNA_lengths',
             'rRNA_lengths',
             'clean_lengths',
             'unmapped_lengths',
             'rRNA_coverage',
             'oligo_hit_lengths',
             'clean_bam',
            ],
            ['from_starts',
             'from_ends',
             'rpf_positions',
             'expression',
             'mismatches_25',
             'mismatches_26',
             'mismatches_27',
             'mismatches_28',
             'mismatches_29',
             'mismatches_30',
             'mismatches_31',
             'recycling_ratios',
            ],
        ]

        specific_work = [
            [(self.trim_reads, 'Trim reads'),
             (self.pre_filter_rRNA, 'Filter contaminants'),
             (self.tophat, 'Map with tophat'),
             (self.post_filter_contaminants, 'Further contaminant filtering'),
             (self.get_rRNA_coverage, 'Counting rRNA coverage'),
             (self.get_oligo_hit_lengths, 'Counting oligo hit length distributions'),
            ],
            [(self.get_aggregate_positions, 'Counting mapping positions'),
             (self.get_error_profile, 'Getting error profile'),
             (self.get_recycling_ratios, 'Getting recycling ratios'),
            ],
        ]

        specific_cleanup = [
            [self.compute_yield,
             self.plot_lengths,
             self.plot_rRNA_coverage,
             self.plot_oligo_hit_lengths,
             self.index_clean_bam,
            ],
            [self.plot_starts_and_ends,
             self.plot_mismatch_positions,
             self.plot_mismatch_types,
            ],
        ]

        self.results_files.extend(specific_results_files)
        self.figure_files.extend(specific_figure_files)
        mapreduce.extend_stages(self.outputs, specific_outputs)
        mapreduce.extend_stages(self.work, specific_work)
        mapreduce.extend_stages(self.cleanup, specific_cleanup)

        self.make_file_names()
        self.data_fns = glob.glob(self.data_dir + '/*.fastq') + glob.glob(self.data_dir + '/*.fq')
        if self.max_read_length == None:
            self.max_read_length = self.get_max_read_length()
        else:
            self.max_read_length = int(self.max_read_length)
        
        if self.adapter_type == 'truseq':
            self.trim_function = trim.trim_truseq
        elif self.adapter_type == 'linker':
            self.trim_function = trim.trim_linker
        elif self.adapter_type == 'polyA':
            self.trim_function = trim.trim_poly_A
        elif self.adapter_type == 'nothing':
            self.trim_function = trim.trim_nothing

        for key, tail in self.organism_files:
            self.file_names[key] = '{0}/{1}'.format(self.organism_dir, tail)

    def get_simple_CDSs(self):
        all_simple_CDSs = gtf.get_simple_CDSs(self.file_names['genes'])
        max_gene_length = max(abs(gene.end - gene.start) + 1 for gene in all_simple_CDSs)
        piece_simple_CDSs = Parallel.split_file.piece_of_list(all_simple_CDSs,
                                                              self.num_pieces,
                                                              self.which_piece,
                                                             )
        return piece_simple_CDSs, max_gene_length
        
    def trim_reads(self):
        trimmed_lengths, too_short_lengths = self.trim_function(self.get_reads(),
                                                                self.min_length,
                                                                self.max_read_length,
                                                                self.file_names['trimmed_reads'],
                                                               )
        self.write_file('trimmed_lengths', trimmed_lengths)
        self.write_file('too_short_lengths', too_short_lengths)

    def pre_filter_rRNA(self):
        contaminants.pre_filter(self.file_names['rRNA_index'],
                                self.file_names['trimmed_reads'],
                                self.file_names['filtered_reads'],
                                self.file_names['rRNA_sam'],
                                self.file_names['rRNA_bam'],
                                self.file_names['bowtie_error'],
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
        contaminants.post_filter(self.file_names['accepted_hits'],
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
        data = contaminants.produce_rRNA_coverage(bam_file_names)
        self.write_file('rRNA_coverage', data)

    def get_oligo_hit_lengths(self):
        lengths = contaminants.get_oligo_hit_lengths(self.file_names['rRNA_bam'],
                                                     self.file_names['oligos'],
                                                     self.file_names['oligos_sam'],
                                                     self.max_read_length,
                                                    )
        self.write_file('oligo_hit_lengths', lengths)

    def compute_yield(self):
        trimmed_lengths = self.read_file('trimmed_lengths')
        too_short_lengths = self.read_file('too_short_lengths')
        tRNA_lengths = self.read_file('tRNA_lengths')
        rRNA_lengths = self.read_file('rRNA_lengths')
        unmapped_lengths = self.read_file('unmapped_lengths')
        clean_lengths = self.read_file('clean_lengths')

        total_reads = trimmed_lengths.sum() + too_short_lengths.sum()
        long_enough_reads = trimmed_lengths.sum()
        rRNA_reads = rRNA_lengths.sum()
        tRNA_reads = tRNA_lengths.sum()
        unmapped_reads = unmapped_lengths.sum()
        clean_reads = clean_lengths.sum()

        with open(self.file_names['yield'], 'w') as yield_file:
            yield_file.write('Total reads: {0:,}\n'.format(total_reads))
            for category, count in [('Long enough reads', long_enough_reads),
                                    ('rRNA reads', rRNA_reads),
                                    ('tRNA reads', tRNA_reads),
                                    ('Unmapped reads', unmapped_reads),
                                    ('Clean reads', clean_reads),
                                   ]:
                fraction = float(count) / total_reads
                line = '{0}: {1:,} ({2:.2%})\n'.format(category,
                                                       count,
                                                       fraction,
                                                      )
                yield_file.write(line)

    def index_clean_bam(self):
        sam.index_bam(self.file_names['clean_bam'])

    def plot_lengths(self):
        too_short_lengths = self.read_file('too_short_lengths')
        tRNA_lengths = self.read_file('tRNA_lengths')
        rRNA_lengths = self.read_file('rRNA_lengths')
        unmapped_lengths = self.read_file('unmapped_lengths')
        clean_lengths = self.read_file('clean_lengths')

        fig_all, ax_all = plt.subplots(figsize=(12, 8))
    
        ax_all.plot(rRNA_lengths, '.-', label='rRNA', color='red')
        ax_all.plot(tRNA_lengths, '.-', label='tRNA', color='blue')
        ax_all.plot(unmapped_lengths, '.-', label='unmapped', color='cyan')
        ax_all.plot(too_short_lengths, '.-', label='too short', color='purple')
        ax_all.plot(clean_lengths, '.-', label='clean', color='green')
        ax_all.set_xlim(0, self.max_read_length)
        ax_all.set_title('Fragment length distribution by source')
        ax_all.set_xlabel('Length of original RNA fragment')
        ax_all.set_ylabel('Number of reads')
        leg = ax_all.legend(loc='upper right', fancybox=True)
        leg.get_frame().set_alpha(0.5)
        fig_all.savefig(self.figure_file_names['all_lengths'])
        
        fig_clean, ax_clean = plt.subplots(figsize=(12, 8))
        
        ax_clean.plot(clean_lengths, '.-', label='clean', color='green')
        ax_clean.axvspan(27.5, 28.5, color='green', alpha=0.2)
        ax_clean.set_xlim(0, 50)
        ax_clean.legend()
        ax_clean.set_title('Fragment length distribution by source')
        ax_clean.set_xlabel('Length of original RNA fragment')
        ax_clean.set_ylabel('Number of reads')
        leg = ax_clean.legend(loc='upper right', fancybox=True)
        leg.get_frame().set_alpha(0.5)
    
        fig_clean.savefig(self.figure_file_names['clean_lengths'])

    def plot_rRNA_coverage(self):
        coverage_data = {self.name: self.read_file('rRNA_coverage')}
        contaminants.plot_rRNA_coverage(coverage_data,
                                        self.file_names['oligos_sam'],
                                        self.figure_file_names['rRNA_coverage_template'],
                                       )

    def plot_oligo_hit_lengths(self):
        lengths = self.read_file('oligo_hit_lengths')
        contaminants.plot_oligo_hit_lengths(self.file_names['oligos'],
                                            lengths,
                                            self.figure_file_names['oligo_hit_lengths'],
                                           )


    def plot_starts_and_ends(self):
        from_starts = self.read_file('from_starts')
        from_ends = self.read_file('from_ends')

        ribosomes.plot_starts_and_ends([from_starts],
                                       [from_ends],
                                       [self.name],
                                       self.figure_file_names['starts_and_ends'],
                                      )

    def plot_mismatch_positions(self):
        fig, ax = plt.subplots(figsize=(16, 12))
        for length in range(25, 32):
            type_counts = self.read_file('mismatches_{0}'.format(length))
            Visualize.mismatches.plot_read_positions(type_counts[:length], str(length), ax=ax, min_q=30)
        ax.set_xlim(-1, 31)
        leg = ax.legend(loc='upper right', fancybox=True)  
        leg.get_frame().set_alpha(0.5)
        fig.savefig(self.figure_file_names['mismatch_positions'])

    def plot_mismatch_types(self):
        first_data_sets = []
        last_data_sets = []
        for length in range(27, 32):
            type_counts = self.read_file('mismatches_{0}'.format(length))
            first_counts = type_counts[0]
            last_counts = type_counts[length - 1]
            first_data_set = (str(length),
                              first_counts,
                              (30, fastq.MAX_EXPECTED_QUAL),
                             )
            first_data_sets.append(first_data_set)
            last_data_set = (str(length),
                             last_counts,
                             (30, fastq.MAX_EXPECTED_QUAL),
                            )
            last_data_sets.append(last_data_set)

        Visualize.mismatches.plot_rates(first_data_sets,
                                        save_as=self.figure_file_names['first_mismatch_types'],
                                       )

        Visualize.mismatches.plot_rates(last_data_sets,
                                        save_as=self.figure_file_names['last_mismatch_types'],
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

    def get_aggregate_positions(self):
        piece_simple_CDSs, max_gene_length = self.get_simple_CDSs()
        genes, from_starts, from_ends = ribosomes.get_aggregate_positions(self.merged_file_names['clean_bam'],
                                                                          piece_simple_CDSs,
                                                                          max_gene_length,
                                                                          self.max_read_length,
                                                                         )
        self.write_file('rpf_positions', genes)
        self.write_file('expression', genes)
        self.write_file('from_starts', from_starts)
        self.write_file('from_ends', from_ends)

    def get_error_profile(self):
        piece_simple_CDSs, _ = self.get_simple_CDSs()
        type_counts = ribosomes.error_profile(self.merged_file_names['clean_bam'],
                                              piece_simple_CDSs,
                                              self.fastq_type,
                                             )
        for length, length_counts in zip(range(25, 32), type_counts): 
            self.write_file('mismatches_{0}'.format(length), length_counts)

    def get_recycling_ratios(self):
        piece_simple_CDSs, _ = self.get_simple_CDSs()
        rpf_positions_dict = self.read_file('rpf_positions')
        genome = mapping.load_genome(self.file_names['genome'], explicit_path=True)
        ratio_lists = ribosomes.recycling_ratios(rpf_positions_dict,
                                                 piece_simple_CDSs,
                                                 genome,
                                                )
        self.write_file('recycling_ratios', ratio_lists)

if __name__ == '__main__':
    script_path = os.path.realpath(__file__)
    mapreduce.controller(RibosomeProfilingExperiment, script_path)
