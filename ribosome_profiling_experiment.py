import matplotlib
matplotlib.use('Agg', warn=False)
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import sys
import pysam
from itertools import chain, product
from collections import Counter
import trim
import ribosomes
import positions
import contaminants
import composition
import gtf
import codons
import visualize
from Sequencing import mapping_tools, fastq, sam, utilities
from Sequencing.Parallel import map_reduce, piece_of_list, split_file
from Circles import Visualize

class RibosomeProfilingExperiment(map_reduce.MapReduceExperiment):
    num_stages = 2

    def __init__(self, **kwargs):
        map_reduce.MapReduceExperiment.__init__(self, **kwargs)

        self.data_dir = kwargs['data_dir'].rstrip('/')
        self.adapter_type = kwargs['adapter_type']
        self.organism_dir = kwargs['organism_dir'].rstrip('/')
        self.transcripts_file_name = kwargs.get('transcripts_file_name', None)
        self.possibly_misannotated_file_name = kwargs.get('possibly_misannotated_file_name', None)
        
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
            ('clean_composition', 'array_3d', '{name}_clean_composition.npy'),
            ('clean_composition_perfect', 'array_3d', '{name}_clean_composition_perfect.npy'),

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
            ('more_rRNA_bam_sorted', 'bam', '{name}_more_rRNA_sorted.bam'),
            ('tRNA_bam', 'bam', '{name}_tRNA.bam'),
            ('tRNA_bam_sorted', 'bam', '{name}_tRNA_sorted.bam'),
            ('other_ncRNA_bam', 'bam', '{name}_other_ncRNA.bam'),
            ('other_ncRNA_bam_sorted', 'bam', '{name}_other_ncRNA_sorted.bam'),
            ('unambiguous_bam', 'bam', '{name}_unambiguous.bam'),

            ('mismatches', 'mismatches', '{name}_mismatches.npy'),

            ('read_positions', 'read_positions', '{name}_read_positions.txt'),
            ('read_positions_remapped', 'read_positions', '{name}_read_positions_remapped.txt'),
            ('unambiguous_read_positions', 'read_positions', '{name}_unambiguous_read_positions.txt'),
            ('from_starts_and_ends', 'read_positions', '{name}_from_starts_and_ends.txt'),
            ('unambiguous_from_starts', 'read_positions', '{name}_unambiguous_from_starts.txt'),
            ('from_starts_and_ends_remapped', 'read_positions', '{name}_from_starts_and_ends_remapped.txt'),
            ('codon_counts', 'codon_counts', '{name}_codon_counts.txt'),
            ('codon_counts_stringent', 'codon_counts', '{name}_codon_counts_stringent.txt'),
            ('buffered_codon_counts', 'read_positions', '{name}_buffered_codon_counts.txt'),
            ('metacodon_counts', 'read_positions', '{name}_metacodon_counts.txt'),
            ('metacodon_counts_stringent', 'read_positions', '{name}_metacodon_counts_stringent.txt'),
            ('metacodon_counts_nucleotide_resolution', 'read_positions', '{name}_metacodon_counts_nucleotide_resolution.txt'),
            ('mean_densities', 'read_positions', '{name}_mean_densities.txt'),
            ('mean_densities_no_misannotated', 'read_positions', '{name}_mean_densities_no_misannotated.txt'),
            ('RPKMs', 'RPKMs', '{name}_RPKMs.txt'),
            ('RPKMs_exclude_edges', 'RPKMs', '{name}_RPKMs_exclude_edges.txt'),
            ('read_counts', 'expression', '{name}_read_counts.txt'),
            ('read_counts_exclude_edges', 'expression', '{name}_read_counts_exclude_edges.txt'),

            ('rRNA_coverage', 'coverage', '{name}_rRNA_coverage.txt'),
            ('rRNA_coverage_length_28', 'coverage', '{name}_rRNA_coverage_length_28.txt'),
            ('dominant_stretches', 'PH', '{name}_dominant_stretches.txt'),

            ('tophat_dir', 'dir', 'tophat'),
            ('accepted_hits', 'bam', 'tophat/accepted_hits.bam'),
            ('unmapped_bam', 'bam', 'tophat/unmapped.bam'),
            ('unmapped_trimmed_fastq', 'fastq', '{name}_unmapped_trimmed.fastq'),
            ('tophat_remapped_polyA_dir', 'dir', 'tophat_remapped_polyA'),
            ('remapped_accepted_hits', 'bam', 'tophat_remapped_polyA/accepted_hits.bam'),
            ('remapped_clean_bam', 'bam', '{name}_remapped_clean.bam'),
        
            ('recycling_ratios', 'ratios', '{name}_recycling_ratios.txt'),

            ('reciprocal_rates', 'pickle', '{name}_reciprocal_rates.pkl'),

            ('yield', '', '{name}_yield.txt'),
        ]

        specific_figure_files = [
            ('quality', '{name}_quality.png'),
            ('clean_composition', '{name}_clean_composition.pdf'),
            ('clean_composition_perfect', '{name}_clean_composition_perfect.pdf'),
            ('all_lengths', '{name}_all_lengths.pdf'),
            ('clean_lengths', '{name}_clean_lengths.pdf'),
            ('rRNA_coverage_template', '{name}_rRNA_coverage_{{0}}.pdf'),
            ('rRNA_coverage_length_28_template', '{name}_rRNA_coverage_length_28_{{0}}.pdf'),
            ('oligo_hit_lengths', '{name}_oligo_hit_lengths.pdf'),
            ('starts_and_ends', '{name}_starts_and_ends.pdf'),
            ('starts_and_ends_zoomed_out', '{name}_starts_and_ends_zoomed_out.pdf'),
            ('starts_and_ends_remapped', '{name}_starts_and_ends_remapped.pdf'),
            ('starts_and_ends_remapped_zoomed_out', '{name}_starts_and_ends_remapped_zoomed_out.pdf'),
            ('mean_densities', '{name}_mean_densities.pdf'),
            ('mismatch_positions', '{name}_mismatch_positions.png'),
            ('first_mismatch_types', '{name}_first_mismatch_types.png'),
            ('last_mismatch_types', '{name}_last_mismatch_types.png'),
            ('frames', '{name}_frames.pdf'),
            ('unambiguous_frames', '{name}_unambiguous_frames.pdf'),
            ('metacodon_counts', '{name}_metacodon_counts.pdf'),
            ('metacodon_counts_stringent', '{name}_metacodon_counts_stringent.pdf'),
            ('metacodon_mean_enrichments', '{name}_metacodon_mean_enrichments.pdf'),
            ('metacodon_counts_nucleotide_resolution', '{name}_metacodon_counts_nucleotide_resolution.pdf'),
            ('metanucleotide_counts', '{name}_metanucleotide_counts.pdf'),
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
            ['clean_composition',
             'clean_composition_perfect',
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
             'remapped_clean_bam',
             'rRNA_bam',
             'more_rRNA_bam',
             'tRNA_bam_sorted',
             'other_ncRNA_bam_sorted',
            ],
            ['read_positions',
             'buffered_codon_counts',
             'codon_counts',
             'metacodon_counts',
             #'metacodon_counts_stringent',
             #'metacodon_counts_nucleotide_resolution',
             'from_starts_and_ends',
             'from_starts_and_ends_remapped',
             'read_counts',
             'read_counts_exclude_edges',
             'mismatches',
            ],
        ]

        if self.adapter_type == 'polyA':
            specific_outputs[0].extend(['unambiguous_lengths',
                                        'unambiguous_bam',
                                       ])
            specific_outputs[1].extend(['unambiguous_read_positions',
                                        'unambiguous_from_starts'
                                       ])

        specific_work = [
            [(self.trim_reads, 'Trim reads'),
             (self.pre_filter_contaminants, 'Filter contaminants'),
             (self.map_tophat, 'Map with tophat'),
             (self.remap_trimmed, 'Remapping trimmed unmapped read'),
             (self.post_filter_contaminants, 'Further contaminant filtering'),
             (self.compute_base_composition, 'Computing base composition'),
             (self.find_unambiguous_lengths, 'Finding unambiguous lengths'),
             (self.get_rRNA_coverage, 'Counting rRNA coverage'),
             (self.get_oligo_hit_lengths, 'Counting oligo hit length distributions'),
            ],
            [(self.get_read_positions, 'Counting mapping positions'),
             (self.get_metagene_positions, 'Aggregating metagene positions'),
             (self.compute_total_read_counts, 'Computing total read counts'),
             (self.compute_codon_occupancy_counts, 'Computing codon occupancies'),
             (self.compute_metacodon_counts, 'Computing metacodon counts'),
             (self.get_error_profile, 'Getting error profile'),
            ],
        ]

        specific_cleanup = [
            [self.compute_yield,
             self.plot_base_composition,
             self.plot_lengths,
             self.plot_rRNA_coverage,
             self.plot_oligo_hit_lengths,
            ],
            [self.compute_RPKMs,
             self.compute_mean_densities,
             self.plot_starts_and_ends,
             self.plot_frames,
             self.plot_mismatch_positions,
             self.plot_mismatch_types,
             self.plot_metacodon_counts,
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

    def compute_base_composition(self):
        seq_info_pairs = composition.get_seq_info_pairs(self.file_names['clean_bam'])
        all_array, perfect_array = composition.length_stratified_composition(seq_info_pairs, self.max_read_length)
        
        self.write_file('clean_composition', all_array)
        self.write_file('clean_composition_perfect', perfect_array)

    def plot_base_composition(self):
        all_array = self.read_file('clean_composition', merged=True)
        composition.length_stratified_plot(all_array,
                                           self.figure_file_names['clean_composition'],
                                          )
        
        perfect_array = self.read_file('clean_composition_perfect', merged=True)
        composition.length_stratified_plot(perfect_array,
                                           self.figure_file_names['clean_composition_perfect'],
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
        ribosomes.map_tophat([self.file_names['filtered_reads']],
                             self.file_names['bowtie2_index_prefix'],
                             self.file_names['genes'],
                             self.file_names['transcriptome_index'],
                             self.file_names['tophat_dir'],
                            )
        pysam.index(self.file_names['accepted_hits'])

    def remap_trimmed(self):
        ribosomes.trim_polyA_from_unmapped(self.file_names['unmapped_bam'],
                                           self.file_names['unmapped_trimmed_fastq'],
                                           self.min_length,
                                           self.max_read_length,
                                           second_time=True,
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

    def post_filter_contaminants(self):
        contaminants.post_filter(self.file_names['accepted_hits'],
                                 self.file_names['genes'],
                                 self.file_names['clean_bam'],
                                 self.file_names['more_rRNA_bam'],
                                 self.file_names['more_rRNA_bam_sorted'],
                                 self.file_names['tRNA_bam'],
                                 self.file_names['tRNA_bam_sorted'],
                                 self.file_names['other_ncRNA_bam'],
                                 self.file_names['other_ncRNA_bam_sorted'],
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

    def find_unambiguous_lengths(self):
        if self.adapter_type == 'polyA':
            trim.unambiguously_trimmed(self.file_names['clean_bam'],
                                       self.file_names['unambiguous_bam'],
                                       self.file_names['genome'],
                                      )
            unambiguous_length_counts = sam.get_length_counts(self.file_names['unambiguous_bam'])
            unambiguous_lengths = self.zero_padded_array(unambiguous_length_counts)
        else:
            # Need to write the file so that there is something to merge.
            unambiguous_lengths = np.zeros(self.max_read_length + 1, int)
        self.write_file('unambiguous_lengths', unambiguous_lengths)
    
    def get_rRNA_coverage(self):
        data = contaminants.produce_rRNA_coverage(self.file_names['rRNA_bam'])
        self.write_file('rRNA_coverage', data)

        data_length_28 = contaminants.produce_rRNA_coverage(self.file_names['rRNA_bam'], specific_length=28)
        self.write_file('rRNA_coverage_length_28', data_length_28)

    def get_oligo_hit_lengths(self):
        lengths = contaminants.get_oligo_hit_lengths(self.file_names['rRNA_bam'],
                                                     self.file_names['oligos'],
                                                     self.file_names['oligos_sam'],
                                                     self.max_read_length,
                                                    )
        self.write_file('oligo_hit_lengths', lengths)

    def compute_yield(self):
        trimmed_lengths = self.read_file('trimmed_lengths', merged=True)
        too_short_lengths = self.read_file('too_short_lengths', merged=True)
        rRNA_lengths = self.read_file('rRNA_lengths', merged=True)
        synthetic_lengths = self.read_file('synthetic_lengths', merged=True)
        tRNA_lengths = self.read_file('tRNA_lengths', merged=True)
        other_ncRNA_lengths = self.read_file('other_ncRNA_lengths', merged=True)
        unmapped_lengths = self.read_file('unmapped_lengths', merged=True)
        clean_lengths = self.read_file('clean_lengths', merged=True)

        total_reads = trimmed_lengths.sum() + too_short_lengths.sum()
        long_enough_reads = trimmed_lengths.sum()
        rRNA_reads = rRNA_lengths.sum()
        synthetic_reads = synthetic_lengths.sum()
        tRNA_reads = tRNA_lengths.sum()
        other_ncRNA_reads = other_ncRNA_lengths.sum()
        unmapped_reads = unmapped_lengths.sum()
        clean_reads = clean_lengths.sum()

        dominant_reads, boundaries = contaminants.identify_dominant_stretches(self.read_file('rRNA_coverage', merged=True),
                                                                              total_reads,
                                                                              self.merged_file_names['rRNA_bam'],
                                                                             )
        other_reads = rRNA_reads - dominant_reads

        with open(self.merged_file_names['dominant_stretches'], 'w') as dominant_stretches_file:
            for rname in sorted(boundaries):
                for start, stop in boundaries[rname]:
                    dominant_stretches_file.write('{0}:{1}-{2}\n'.format(rname, start, stop))

        with open(self.file_names['yield'], 'w') as yield_file:
            yield_file.write('Total reads: {0:,}\n'.format(total_reads))
            for category, count in [('Long enough reads', long_enough_reads),
                                    ('rRNA reads', rRNA_reads),
                                    ('(rRNA reads from non-dominant stetches)', other_reads),
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

    def plot_lengths(self):
        too_short_lengths = self.read_file('too_short_lengths')
        rRNA_lengths = self.read_file('rRNA_lengths')
        synthetic_lengths = self.read_file('synthetic_lengths')
        tRNA_lengths = self.read_file('tRNA_lengths')
        other_ncRNA_lengths = self.read_file('other_ncRNA_lengths')
        unmapped_lengths = self.read_file('unmapped_lengths')
        clean_lengths = self.read_file('clean_lengths')

        if self.adapter_type == 'polyA':
            unambiguous_lengths = self.read_file('unambiguous_lengths')

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
        
        normalized_clean_lengths = np.true_divide(clean_lengths, clean_lengths.sum())
        ax_clean.plot(normalized_clean_lengths, '.-', label='clean', color='green')
        if self.adapter_type == 'polyA':
            normalized_unambiguous_lengths = np.true_divide(unambiguous_lengths, unambiguous_lengths.sum())
            ax_clean.plot(normalized_unambiguous_lengths, '.-', label='unambiguous', color='blue')
        ax_clean.axvspan(27.5, 28.5, color='green', alpha=0.2)
        ax_clean.set_xlim(0, self.max_read_length)
        ax_clean.legend()
        ax_clean.set_title('Fragment length distribution by source')
        ax_clean.set_xlabel('Length of original RNA fragment')
        ax_clean.set_ylabel('Fraction of clean reads')
        ax_clean.legend(loc='upper left', framealpha=0.5)
    
        fig_clean.savefig(self.figure_file_names['clean_lengths'])

    def plot_rRNA_coverage(self):
        trimmed_lengths = self.read_file('trimmed_lengths')
        too_short_lengths = self.read_file('too_short_lengths')

        total_reads = trimmed_lengths.sum() + too_short_lengths.sum()
        counts = self.read_file('rRNA_coverage')
        coverage_data = {self.name: (total_reads, counts, 'blue')}
        contaminants.plot_rRNA_coverage(coverage_data,
                                        self.file_names['oligos_sam'],
                                        self.figure_file_names['rRNA_coverage_template'],
                                       )
        
        counts = self.read_file('rRNA_coverage_length_28')
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

        visualize.plot_averaged_codon_densities([(self.name, self.read_file('mean_densities'), 0)],
                                                self.figure_file_names['mean_densities'],
                                                past_edge=10,
                                                plot_up_to=100,
                                                smooth=False,
                                               )

    def plot_frames(self):
        from_starts_and_ends = self.read_file('from_starts_and_ends')

        visualize.plot_frames(from_starts_and_ends['from_starts'],
                              self.figure_file_names['frames'],
                             )

        if self.adapter_type == 'polyA':
            from_starts = self.read_file('unambiguous_from_starts')

            visualize.plot_frames(from_starts['from_starts'],
                                  self.figure_file_names['unambiguous_frames'],
                                 )

    def plot_metacodon_counts(self):
        metacodon_counts = self.read_file('metacodon_counts', merged=True)
        visualize.plot_metacodon_counts(metacodon_counts, self.figure_file_names['metacodon_counts'])
        #positions.plot_metacodon_counts(metacodon_counts, self.figure_file_names['metacodon_mean_enrichments'], enrichment=True)
        
        #metacodon_counts_stringent = self.read_file('metacodon_counts_stringent', merged=True)
        #positions.plot_metacodon_counts(metacodon_counts_stringent, self.figure_file_names['metacodon_counts_stringent'])
        #
        #metacodon_counts_nucleotide_resolution = self.read_file('metacodon_counts_nucleotide_resolution', merged=True)
        #positions.plot_metacodon_counts(metacodon_counts_nucleotide_resolution,
        #                                self.figure_file_names['metacodon_counts_nucleotide_resolution'],
        #                                keys_to_plot=[28, 29, 30, 31, 32],
        #                               )
        #all_codon_ids = metacodon_counts_nucleotide_resolution.keys()
        #positions.plot_metacodon_counts(metacodon_counts_nucleotide_resolution,
        #                                self.figure_file_names['metanucleotide_counts'],
        #                                #codon_ids=['T', 'C', 'A', 'G'] + list(''.join(pair) for pair in product('TCAG', repeat=2)),
        #                                codon_ids=[c for c in all_codon_ids if c.endswith('C') and len(c) < 3],
        #                                #codon_ids=[c for c in codons.non_stop_codons if c.endswith('CG') or c.endswith('CA')],
        #                                #codon_ids=[c for c in codons.non_stop_codons if c.startswith('AC')],
        #                                keys_to_plot=[28, 29, 30, 31, 32],
        #                               )

    def plot_mismatch_positions(self):
        fig, ax = plt.subplots(figsize=(16, 12))
        length_to_index = {length: i for i, length in enumerate(self.relevant_lengths)}
        type_counts = self.read_file('mismatches')
        for length in self.relevant_lengths:
            Visualize.mismatches.plot_read_positions(type_counts[length_to_index[length]][:length],
                                                     str(length),
                                                     ax=ax,
                                                     min_q=30,
                                                    )
        ax.set_xlim(-1, max(self.relevant_lengths))
        ax.legend(loc='upper right', framealpha=0.5)
        fig.savefig(self.figure_file_names['mismatch_positions'])

    def plot_mismatch_types(self):
        first_data_sets = []
        last_data_sets = []
        type_counts = self.read_file('mismatches')
        length_to_index = {length: i for i, length in enumerate(self.relevant_lengths)}
        for length in self.relevant_lengths:
            first_counts = type_counts[length_to_index[length]][0]
            last_counts = type_counts[length_to_index[length]][length - 1]
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
        piece_CDSs, max_gene_length = self.get_CDSs()
        gene_infos = positions.get_Transcript_position_counts(self.merged_file_names['clean_bam'],
                                                              piece_CDSs,
                                                              relevant_lengths=self.relevant_lengths,
                                                             )

        self.read_positions = {name: info['position_counts']
                               for name, info in gene_infos.iteritems()}

        self.write_file('read_positions', self.read_positions)
        
        remapped_gene_infos = positions.get_Transcript_position_counts(self.merged_file_names['remapped_clean_bam'],
                                                                       piece_CDSs,
                                                                       relevant_lengths=self.relevant_lengths,
                                                                      )

        self.remapped_read_positions = {name: info['position_counts']
                                        for name, info in remapped_gene_infos.iteritems()}

        if self.adapter_type == 'polyA':
            gene_infos = positions.get_Transcript_position_counts(self.merged_file_names['unambiguous_bam'],
                                                                  piece_CDSs,
                                                                  relevant_lengths=self.relevant_lengths,
                                                                 )
            self.unambiguous_read_positions = {name: info['position_counts']
                                               for name, info in gene_infos.iteritems()}
            self.write_file('unambiguous_read_positions', self.unambiguous_read_positions)

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

    def get_metagene_positions(self):
        piece_CDSs, max_gene_length = self.get_CDSs()
        read_positions = self.load_read_positions()
        from_starts_and_ends = positions.compute_metagene_positions(read_positions, max_gene_length)
        remapped_read_positions = self.load_read_positions(modifier='remapped')
        from_starts_and_ends_remapped = positions.compute_metagene_positions(remapped_read_positions, max_gene_length)

        self.write_file('from_starts_and_ends', from_starts_and_ends)
        self.write_file('from_starts_and_ends_remapped', from_starts_and_ends_remapped)
        
        if self.adapter_type == 'polyA':
            read_positions = self.load_read_positions(modifier='unambiguous')
            from_starts_and_ends = positions.compute_metagene_positions(read_positions, max_gene_length)

            self.write_file('unambiguous_from_starts', from_starts_and_ends)
        
    def get_error_profile(self):
        piece_CDSs, _ = self.get_CDSs()
        type_counts = ribosomes.error_profile(self.merged_file_names['clean_bam'],
                                              piece_CDSs,
                                              self.relevant_lengths,
                                              self.max_read_length,
                                             )
        self.write_file('mismatches', type_counts)

    def get_recycling_ratios(self):
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
        read_positions = self.load_read_positions()

        buffered_codon_counts = {}
        codon_counts = {}
        codon_counts_stringent = {}
        for name, position_counts in read_positions.iteritems():
            buffered_counts = positions.compute_codon_counts(position_counts, self.offset_type)
            buffered_counts_stringent = positions.compute_codon_counts(position_counts, self.offset_type + '_stringent')
            buffered_codon_counts[name] = {'relaxed': buffered_counts,
                                           'stringent': buffered_counts_stringent,
                                          }
            num_codons = buffered_counts.CDS_length
            # + 1 is to include the stop codon
            codon_counts[name] = buffered_counts['start_codon', :num_codons + 1]
            codon_counts_stringent[name] = buffered_counts_stringent['start_codon', :num_codons + 1]

        self.write_file('buffered_codon_counts', buffered_codon_counts)
        self.write_file('codon_counts', codon_counts)
        self.write_file('codon_counts_stringent', codon_counts_stringent)

    def compute_total_read_counts(self):
        read_positions = self.load_read_positions()
        read_counts = positions.compute_read_counts(read_positions, 0, 0)
        self.write_file('read_counts', read_counts)
        read_counts_exclude_edges = positions.compute_read_counts(read_positions, 30, 4)
        self.write_file('read_counts_exclude_edges', read_counts_exclude_edges)

    def compute_metacodon_counts(self):
        codon_counts = self.read_file('codon_counts')
        metacodon_counts = positions.compute_metacodon_counts(codon_counts,
                                                              self.file_names['genes'],
                                                              self.file_names['genome'],
                                                             )
        self.write_file('metacodon_counts', metacodon_counts)

        #codon_counts_stringent = positions.read_codon_counts_file(self.file_names['codon_counts_stringent'])
        #metacodon_counts_stringent = positions.compute_metacodon_counts(codon_counts_stringent,
        #                                                                self.file_names['genes'],
        #                                                                self.file_names['genome'],
        #                                                               )
        #self.write_file('metacodon_counts_stringent', metacodon_counts_stringent)
        #
        #read_positions = self.load_read_positions()
        #metacodon_counts_nucleotide_resolution = positions.compute_metacodon_counts_nucleotide_resolution(read_positions,
        #                                                                                                  self.file_names['genes'],
        #                                                                                                  self.file_names['genome'],
        #                                                                                                 )
        #self.write_file('metacodon_counts_nucleotide_resolution', metacodon_counts_nucleotide_resolution)

    def compute_mean_densities(self):
        codon_counts = self.read_file('buffered_codon_counts', merged=True)
        mean_densities = positions.compute_averaged_codon_densities(codon_counts)
        self.write_file('mean_densities', mean_densities)

        if self.possibly_misannotated_file_name != None:
            possibly_misannotated_names = {line.strip() for line in open(self.possibly_misannotated_file_name)}
            mean_densities = positions.compute_averaged_codon_densities(codon_counts, possibly_misannotated_names)
            self.write_file('mean_densities_no_misannotated', mean_densities)
        
    def compute_RPKMs(self, exclude_from_start=0, exclude_from_end=0):
        gene_infos = self.read_file('read_counts', merged=True)
        RPKMs = positions.compute_RPKMs(gene_infos, 0, 0)
        self.write_file('RPKMs', RPKMs)
        
        gene_infos = self.read_file('read_counts_exclude_edges', merged=True)
        RPKMs_exclude_edges = positions.compute_RPKMs(gene_infos, 30, 4)
        self.write_file('RPKMs_exclude_edges', RPKMs_exclude_edges)

if __name__ == '__main__':
    script_path = os.path.realpath(__file__)
    map_reduce.controller(RibosomeProfilingExperiment, script_path)
