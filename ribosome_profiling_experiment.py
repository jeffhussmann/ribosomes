import matplotlib
matplotlib.use('Agg', warn=False)
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import sys
import pysam
from itertools import chain
from collections import Counter
import trim
import ribosomes
import positions
import contaminants
import gtf
from Circles import fastq
from Circles import utilities
from Circles import sam
from Circles.Parallel import map_reduce
from Circles.Parallel import piece_of_list
from Circles.Parallel import split_file
from Circles import Visualize
from Circles import mapping_tools

class RibosomeProfilingExperiment(map_reduce.MapReduceExperiment):
    num_stages = 2

    def __init__(self, **kwargs):
        map_reduce.MapReduceExperiment.__init__(self, **kwargs)

        self.data_dir = kwargs['data_dir'].rstrip('/')
        self.adapter_type = kwargs['adapter_type']
        self.organism_dir = kwargs['organism_dir'].rstrip('/')
        self.transcripts_file_name = kwargs.get('transcripts_file_name', None)
        
        # A bowtie2 index of synthetic sequences (markers, adapters) that need
        # to be filtered out can be provided. 
        self.synthetic_index_prefix = kwargs.get('synthetic_index_prefix', None)
        
        # Which mapping of length to A-site offset to use. Varies by experiment
        # because of different digestion and library preparation strategies.
        self.offset_type = kwargs['offset_type']
        
        self.max_read_length = kwargs.get('max_read_length', None)
        length_range = kwargs.get('relevant_lengths', None)
        if length_range == None:
            self.relevant_lengths = range(27, 32)
        else:
            start, stop = map(int, length_range.split(','))
            self.relevant_lengths = range(start, stop + 1)
        
        self.min_length = 10
        
        specific_results_files = [
            ('quality', 'array_2d', '{name}_quality.txt'),
            ('complexity', 'array_2d', '{name}_complexity.txt'),

            ('bowtie_error', '', '{name}_bowtie_error.txt'),
            ('trimmed_reads', 'fastq', '{name}_trimmed.fastq'),
            ('rRNA_filtered_reads', 'fastq', '{name}_rRNA_filtered.fastq'), 
            ('synthetic_filtered_reads', 'fastq', '{name}_synthetic_filtered.fastq'), 

            ('too_short_lengths', 'array_1d', '{name}_too_short_lengths.txt'),
            ('trimmed_lengths', 'array_1d', '{name}_trimmed_lengths.txt'),
            ('tRNA_lengths', 'array_1d', '{name}_tRNA_lengths.txt'),
            ('rRNA_lengths', 'array_1d', '{name}_rRNA_lengths.txt'),
            ('synthetic_lengths', 'array_1d', '{name}_synthetic_lengths.txt'),
            ('other_ncRNA_lengths', 'array_1d', '{name}_ncRNA_lengths.txt'),
            ('clean_lengths', 'array_1d', '{name}_clean_lengths.txt'),
            ('unmapped_lengths', 'array_1d', '{name}_unmapped_lengths.txt'),
            ('unambiguous_lengths', 'array_1d', '{name}_unambiguous_lengths.txt'),
            ('oligo_hit_lengths', 'array_2d', '{name}_oligo_hit_lengths.txt'),

            ('rRNA_sam', 'sam', '{name}_rRNA.sam'),
            ('rRNA_bam', 'bam', '{name}_rRNA.bam'),
            ('synthetic_sam', 'sam', '{name}_synthetic.sam'),
            ('synthetic_bam', 'bam', '{name}_synthetic.bam'),
            ('clean_bam', 'bam', '{name}_clean.bam'),
            ('more_rRNA_bam', 'bam', '{name}_more_rRNA.bam'),
            ('tRNA_bam', 'bam', '{name}_tRNA.bam'),
            ('other_ncRNA_bam', 'bam', '{name}_other_ncRNA.bam'),
            ('unambiguous_bam', 'bam', '{name}_unambiguous.bam'),

            ('read_positions', 'read_positions', '{name}_read_positions.txt'),
            ('from_starts', 'read_positions', '{name}_from_starts.txt'),
            ('from_ends', 'read_positions', '{name}_from_ends.txt'),
            ('codon_counts', 'codon_counts', '{name}_codon_counts.txt'),
            #('codon_counts_strict', 'codon_counts', '{name}_codon_counts_strict.txt'),
            #('codon_counts_relaxed', 'codon_counts', '{name}_codon_counts_relaxed.txt'),
            ('RPKMs', 'RPKMs', '{name}_RPKMs.txt'),
            ('read_counts', 'expression', '{name}_read_counts.txt'),
            ('strand_counts', 'expression', '{name}_strand_counts.txt'),

            ('rRNA_coverage', 'coverage', '{name}_rRNA_coverage.txt'),
            ('rRNA_coverage_length_28', 'coverage', '{name}_rRNA_coverage_length_28.txt'),

            ('tophat_dir', 'dir', 'tophat'),
            ('accepted_hits', 'bam', 'tophat/accepted_hits.bam'),
            ('unmapped_bam', 'bam', 'tophat/unmapped.bam'),

            ('recycling_ratios', 'ratios', '{name}_recycling_ratios.txt'),

            ('yield', '', '{name}_yield.txt'),
        ]

        for length in self.relevant_lengths:
            entry = ('mismatches_{0}'.format(length),
                     'mismatches',
                     '{{name}}_mismatches_{0}.txt'.format(length),
                    )
            specific_results_files.append(entry)

        specific_figure_files = [
            ('quality', '{name}_quality.png'),
            ('complexity', '{name}_complexity.png'),
            ('all_lengths', '{name}_all_lengths.pdf'),
            ('clean_lengths', '{name}_clean_lengths.pdf'),
            ('rRNA_coverage_template', '{name}_rRNA_coverage_{{0}}.pdf'),
            ('rRNA_coverage_length_28_template', '{name}_rRNA_coverage_length_28_{{0}}.pdf'),
            ('oligo_hit_lengths', '{name}_oligo_hit_lengths.pdf'),
            ('starts_and_ends', '{name}_starts_and_ends.pdf'),
            ('mismatch_positions', '{name}_mismatch_positions.png'),
            ('first_mismatch_types', '{name}_first_mismatch_types.png'),
            ('last_mismatch_types', '{name}_last_mismatch_types.png'),
            ('frames', '{name}_frames.pdf'),
        ]

        self.organism_files = [
            ('bowtie2_index_prefix', 'genome/genome'),
            ('genome', 'genome'),
            ('genes', 'transcriptome/genes.gtf'),
            ('dubious_ORFs', 'transcriptome/dubious_ORF_gene_ids.txt'),
            ('transcriptome_index', 'transcriptome/bowtie2_index/genes'),
            ('rRNA_index', 'contaminant/bowtie2_index/rRNA'),
            ('oligos', 'contaminant/subtraction_oligos.fasta'),
            ('oligos_sam', 'contaminant/subtraction_oligos.sam'),
        ]

        specific_outputs = [
            ['quality',
             'complexity',
             'too_short_lengths',
             'trimmed_lengths',
             'tRNA_lengths',
             'rRNA_lengths',
             'synthetic_lengths',
             'other_ncRNA_lengths',
             'clean_lengths',
             'unmapped_lengths',
             'rRNA_coverage',
             'rRNA_coverage_length_28',
             'oligo_hit_lengths',
             'clean_bam',
            ],
            ['read_positions',
             'codon_counts',
             'from_starts',
             'from_ends',
             'strand_counts',
             'read_counts'
             #'recycling_ratios',
            ],
        ]

        specific_outputs[1].extend(['mismatches_{0}'.format(length) for length in self.relevant_lengths])

        specific_work = [
            [(self.collect_data_statistics, 'Collecting data statistics'),
             (self.trim_reads, 'Trim reads'),
             (self.pre_filter_contaminants, 'Filter contaminants'),
             (self.map_tophat, 'Map with tophat'),
             (self.post_filter_contaminants, 'Further contaminant filtering'),
             (self.get_rRNA_coverage, 'Counting rRNA coverage'),
             (self.get_oligo_hit_lengths, 'Counting oligo hit length distributions'),
            ],
            [#(self.get_read_positions, 'Counting mapping positions'),
             (self.get_read_positions_splicing, 'Counting mapping positions, splicing'),
             (self.get_metagene_positions, 'Aggregating metagene positions'),
             (self.compute_codon_occupancy_counts, 'Computing codon occupancies'),
             (self.get_error_profile, 'Getting error profile'),
             #(self.get_recycling_ratios, 'Getting recycling ratios'),
            ],
        ]

        specific_cleanup = [
            [self.plot_data_statistics,
             self.compute_yield,
             self.plot_lengths,
             self.plot_rRNA_coverage,
             self.plot_oligo_hit_lengths,
             self.index_clean_bam,
            ],
            [self.compute_RPKMs,
             self.plot_starts_and_ends,
             self.plot_frames,
             self.plot_mismatch_positions,
             self.plot_mismatch_types,
            ],
        ]

        self.results_files.extend(specific_results_files)
        self.figure_files.extend(specific_figure_files)
        map_reduce.extend_stages(self.outputs, specific_outputs)
        map_reduce.extend_stages(self.work, specific_work)
        map_reduce.extend_stages(self.cleanup, specific_cleanup)

        self.make_file_names()
        self.data_fns = glob.glob(self.data_dir + '/*.fastq') + glob.glob(self.data_dir + '/*.fq')
        if self.max_read_length == None:
            self.max_read_length = self.get_max_read_length()
        else:
            self.max_read_length = int(self.max_read_length)
        
        self.trim_function = trim.bound_trim[self.adapter_type]

        for key, tail in self.organism_files:
            self.file_names[key] = '{0}/{1}'.format(self.organism_dir, tail)

    def get_simple_CDSs(self):
        all_simple_CDSs = gtf.get_simple_CDSs(self.file_names['genes'],
                                              #exclude_from=(self.file_names['dubious_ORFs'],),
                                             )
        max_gene_length = max(abs(gene.end - gene.start) + 1 for gene in all_simple_CDSs)
        piece_simple_CDSs = piece_of_list(all_simple_CDSs,
                                          self.num_pieces,
                                          self.which_piece,
                                         )
        return piece_simple_CDSs, max_gene_length

    def get_CDSs(self):
        if self.transcripts_file_name == None:
            CDSs = gtf.get_simple_CDSs(self.file_names['genes'])
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

    def collect_data_statistics(self):
        q_array, c_array = fastq.quality_and_complexity(self.get_reads(), self.max_read_length)
        
        self.write_file('quality', q_array)
        self.write_file('complexity', c_array)

    def plot_data_statistics(self):
        quality_counts = self.read_file('quality', merged=True)
        Visualize.fastq.plot_quality_histograms(quality_counts,
                                                self.figure_file_names['quality'],
                                               )
        base_counts = self.read_file('complexity', merged=True)
        Visualize.fastq.plot_sequence_complexities(base_counts,
                                                   self.figure_file_names['complexity'],
                                                  )

    def trim_reads(self):
        trimmed_lengths, too_short_lengths, barcode_counts = self.trim_function(self.get_reads(),
                                                                                self.file_names['trimmed_reads'],
                                                                                self.min_length,
                                                                                self.max_read_length,
                                                                               )
        self.write_file('trimmed_lengths', trimmed_lengths)
        self.write_file('too_short_lengths', too_short_lengths)

    def pre_filter_contaminants(self):
        contaminants.pre_filter(self.file_names['rRNA_index'],
                                self.file_names['trimmed_reads'],
                                self.file_names['rRNA_filtered_reads'],
                                self.file_names['rRNA_sam'],
                                self.file_names['rRNA_bam'],
                                self.file_names['bowtie_error'],
                               )

        if self.synthetic_index_prefix:
            contaminants.pre_filter(self.synthetic_index_prefix,
                                    self.file_names['rRNA_filtered_reads'],
                                    self.file_names['synthetic_filtered_reads'],
                                    self.file_names['synthetic_sam'],
                                    self.file_names['synthetic_bam'],
                                    self.file_names['bowtie_error'],
                                   )
            synthetic_length_counts = sam.get_length_counts(self.file_names['synthetic_bam'])
            synthetic_lengths = self.zero_padded_array(synthetic_length_counts)
            self.file_names['filtered_reads'] = self.file_names['synthetic_filtered_reads']
        else:
            synthetic_lengths = np.zeros(self.max_read_length + 1, int)
            self.file_names['filtered_reads'] = self.file_names['rRNA_filtered_reads']

        self.write_file('synthetic_lengths', synthetic_lengths)

    def map_tophat(self):
        ribosomes.map_tophat(self.file_names['filtered_reads'],
                             self.file_names['bowtie2_index_prefix'],
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
                                 self.file_names['other_ncRNA_bam'],
                                )

        tRNA_length_counts = sam.get_length_counts(self.file_names['tRNA_bam'])
        tRNA_lengths = self.zero_padded_array(tRNA_length_counts)
        self.write_file('tRNA_lengths', tRNA_lengths)
        
        other_ncRNA_length_counts = sam.get_length_counts(self.file_names['other_ncRNA_bam'])
        other_ncRNA_lengths = self.zero_padded_array(other_ncRNA_length_counts)
        self.write_file('other_ncRNA_lengths', other_ncRNA_lengths)

        rRNA_length_counts = sam.get_length_counts(self.file_names['rRNA_bam'])
        rRNA_length_counts += sam.get_length_counts(self.file_names['more_rRNA_bam'])
        rRNA_lengths = self.zero_padded_array(rRNA_length_counts)
        self.write_file('rRNA_lengths', rRNA_lengths)
        
        unmapped_length_counts = sam.get_length_counts(self.file_names['unmapped_bam'], only_primary=False)
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

        data_length_28 = contaminants.produce_rRNA_coverage(bam_file_names, specific_length=28)
        self.write_file('rRNA_coverage_length_28', data_length_28)

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
        rRNA_lengths = self.read_file('rRNA_lengths')
        synthetic_lengths = self.read_file('synthetic_lengths')
        tRNA_lengths = self.read_file('tRNA_lengths')
        other_ncRNA_lengths = self.read_file('other_ncRNA_lengths')
        unmapped_lengths = self.read_file('unmapped_lengths')
        clean_lengths = self.read_file('clean_lengths')

        total_reads = trimmed_lengths.sum() + too_short_lengths.sum()
        long_enough_reads = trimmed_lengths.sum()
        rRNA_reads = rRNA_lengths.sum()
        synthetic_reads = synthetic_lengths.sum()
        tRNA_reads = tRNA_lengths.sum()
        other_ncRNA_reads = other_ncRNA_lengths.sum()
        unmapped_reads = unmapped_lengths.sum()
        clean_reads = clean_lengths.sum()

        with open(self.file_names['yield'], 'w') as yield_file:
            yield_file.write('Total reads: {0:,}\n'.format(total_reads))
            for category, count in [('Long enough reads', long_enough_reads),
                                    ('rRNA reads', rRNA_reads),
                                    ('synthetic reads', synthetic_reads),
                                    ('tRNA reads', tRNA_reads),
                                    ('Other ncRNA reads', other_ncRNA_reads),
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
        sam.index_bam(self.merged_file_names['clean_bam'])

    def plot_lengths(self):
        too_short_lengths = self.read_file('too_short_lengths')
        rRNA_lengths = self.read_file('rRNA_lengths')
        synthetic_lengths = self.read_file('synthetic_lengths')
        tRNA_lengths = self.read_file('tRNA_lengths')
        other_ncRNA_lengths = self.read_file('other_ncRNA_lengths')
        unmapped_lengths = self.read_file('unmapped_lengths')
        clean_lengths = self.read_file('clean_lengths')

        fig_all, ax_all = plt.subplots(figsize=(12, 8))
    
        ax_all.plot(rRNA_lengths, '.-', label='rRNA', color='red')
        ax_all.plot(synthetic_lengths, '.-', label='synthetic', color='orange')
        ax_all.plot(tRNA_lengths, '.-', label='tRNA', color='blue')
        ax_all.plot(other_ncRNA_lengths, '.-', label='other_ncRNA', color='black')
        ax_all.plot(unmapped_lengths, '.-', label='unmapped', color='cyan')
        ax_all.plot(too_short_lengths, '.-', label='too short', color='purple')
        ax_all.plot(clean_lengths, '.-', label='clean', color='green')
        ax_all.set_xlim(0, self.max_read_length)
        ax_all.set_title('Fragment length distribution by source')
        ax_all.set_xlabel('Length of original RNA fragment')
        ax_all.set_ylabel('Number of reads')
        ax_all.legend(loc='upper left', framealpha=0.5)
        fig_all.savefig(self.figure_file_names['all_lengths'])
        
        fig_clean, ax_clean = plt.subplots(figsize=(12, 8))
        
        ax_clean.plot(clean_lengths, '.-', label='clean', color='green')
        ax_clean.axvspan(27.5, 28.5, color='green', alpha=0.2)
        ax_clean.set_xlim(0, 50)
        ax_clean.legend()
        ax_clean.set_title('Fragment length distribution by source')
        ax_clean.set_xlabel('Length of original RNA fragment')
        ax_clean.set_ylabel('Number of reads')
        ax_clean.legend(loc='upper left', framealpha=0.5)
    
        fig_clean.savefig(self.figure_file_names['clean_lengths'])

    def plot_rRNA_coverage(self):
        trimmed_lengths = self.read_file('trimmed_lengths')
        too_short_lengths = self.read_file('too_short_lengths')

        total_reads = trimmed_lengths.sum() + too_short_lengths.sum()
        _, counts = self.read_file('rRNA_coverage')
        coverage_data = {self.name: (total_reads, counts, 'blue')}
        contaminants.plot_rRNA_coverage(coverage_data,
                                        self.file_names['oligos_sam'],
                                        self.figure_file_names['rRNA_coverage_template'],
                                       )
        
        _, counts = self.read_file('rRNA_coverage_length_28')
        coverage_data = {self.name: (total_reads, counts, 'blue')}
        contaminants.plot_rRNA_coverage(coverage_data,
                                        self.file_names['oligos_sam'],
                                        self.figure_file_names['rRNA_coverage_length_28_template'],
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

        positions.plot_metagene_positions(from_starts['from_starts']['position_counts'],
                                          from_ends['from_ends']['position_counts'],
                                          self.figure_file_names['starts_and_ends'],
                                         )

    def plot_frames(self):
        from_starts = self.read_file('from_starts')

        positions.plot_frames(from_starts['from_starts']['position_counts'],
                              self.figure_file_names['frames'],
                             )

    def plot_mismatch_positions(self):
        fig, ax = plt.subplots(figsize=(16, 12))
        for length in self.relevant_lengths:
            type_counts = self.read_file('mismatches_{0}'.format(length))
            Visualize.mismatches.plot_read_positions(type_counts[:length], str(length), ax=ax, min_q=30)
        ax.set_xlim(-1, max(self.relevant_lengths))
        ax.legend(loc='upper right', framealpha=0.5)
        fig.savefig(self.figure_file_names['mismatch_positions'])

    def plot_mismatch_types(self):
        first_data_sets = []
        last_data_sets = []
        for length in self.relevant_lengths:
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
        array = utilities.counts_to_array(counts)
        if len(array) < self.max_read_length + 1:
            padded_array = np.zeros(self.max_read_length + 1, int)
            padded_array[:len(array)] += array
        else:
            padded_array = array
        return padded_array

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

    def get_read_positions(self):
        piece_simple_CDSs, max_gene_length = self.get_simple_CDSs()
        genes = positions.get_CDS_position_counts(self.merged_file_names['clean_bam'],
                                                  piece_simple_CDSs,
                                                  relevant_lengths=self.relevant_lengths,
                                                 )
        self.write_file('read_positions', genes)
        self.write_file('strand_counts', genes)
    
    def get_read_positions_splicing(self):
        piece_CDSs, max_gene_length = self.get_CDSs()
        genes = positions.get_Transcript_position_counts(self.merged_file_names['clean_bam'],
                                                         piece_CDSs,
                                                         relevant_lengths=self.relevant_lengths,
                                                        )
        self.write_file('read_positions', genes)
        self.write_file('strand_counts', genes)

    def get_metagene_positions(self):
        #piece_simple_CDSs, max_gene_length = self.get_simple_CDSs()
        piece_CDSs, max_gene_length = self.get_CDSs()
        gene_infos = self.read_file('read_positions')
        from_starts, from_ends = positions.compute_metagene_positions(gene_infos, max_gene_length)

        from_starts = {'from_starts': {'CDS_length': max_gene_length,
                                       'position_counts': from_starts,
                                      }
                      }
        from_ends = {'from_ends': {'CDS_length': max_gene_length,
                                   'position_counts': from_ends,
                                  }
                    }

        self.write_file('from_starts', from_starts)
        self.write_file('from_ends', from_ends)
        
    def get_error_profile(self):
        #piece_simple_CDSs, _ = self.get_simple_CDSs()
        piece_simple_CDSs, _ = self.get_CDSs()
        type_counts = ribosomes.error_profile(self.merged_file_names['clean_bam'],
                                              piece_simple_CDSs,
                                              self.relevant_lengths,
                                             )
        for length in type_counts:
            self.write_file('mismatches_{0}'.format(length), type_counts[length])

    def get_recycling_ratios(self):
        #piece_simple_CDSs, _ = self.get_simple_CDSs()
        piece_simple_CDSs, _ = self.get_CDSs()
        rpf_positions_dict = self.read_file('rpf_positions')
        genome = mapping_tools.load_genome(self.file_names['genome'], explicit_path=True)
        ratio_lists = ribosomes.recycling_ratios(rpf_positions_dict,
                                                 piece_simple_CDSs,
                                                 genome,
                                                )
        self.write_file('recycling_ratios', ratio_lists)

    def compute_codon_occupancy_counts(self):
        transcripts = {line.strip() for line in open(self.transcripts_file_name)}
        gene_infos = self.read_file('read_positions')
        gene_infos = {name: info for name, info in gene_infos.iteritems() if name in transcripts}

        ribosomes.make_codon_counts_file(gene_infos,
                                         self.file_names['codon_counts'],
                                         self.offset_type,
                                        )
        gene_infos = ribosomes.compute_read_counts(gene_infos, stringency='everything')
        self.write_file('read_counts', gene_infos)

    def compute_RPKMs(self):
        gene_infos = self.read_file('read_counts', merged=True)
        ribosomes.make_RPKMs_file(gene_infos, self.file_names['RPKMs'])

if __name__ == '__main__':
    script_path = os.path.realpath(__file__)
    map_reduce.controller(RibosomeProfilingExperiment, script_path)
