import matplotlib
matplotlib.use('Agg', warn=False)
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import os
import pysam
from itertools import product
from collections import Counter
import trim
import ribosomes
import positions
import contaminants
import composition
import codons
import visualize
import rna_experiment
import pausing
import examine_specific_codon
from Sequencing import mapping_tools, fastq, sam, utilities, fasta, annotation
import Sequencing.genomes as genomes
import Sequencing.Visualize
from Sequencing import visualize_structure
from Sequencing.Parallel import map_reduce, piece_of_list, split_file
from Sequencing.Serialize import (array_1d,
                                  array_2d,
                                  array_3d,
                                  mismatches,
                                  counts,
                                 )
from Serialize import (read_positions,
                       coverage,
                       codon_counts,
                       expression,
                       RPKMs,
                      )
from find_polyA_cython import predominantly_A

class RibosomeProfilingExperiment(rna_experiment.RNAExperiment):
    num_stages = 2
    
    specific_results_files = [
        ('clean_composition', array_3d, '{name}_clean_composition.npy'),
        ('clean_composition_perfect', array_3d, '{name}_clean_composition_perfect.npy'),
        ('unmapped_composition', array_3d, '{name}_unmapped_composition.npy'),
        ('quality', array_2d, '{name}_quality.txt'),

        ('bowtie_error', '', '{name}_bowtie_error.txt'),
        ('trimmed_reads', 'fastq', '{name}_trimmed.fastq'),

        ('too_short_lengths', array_1d, '{name}_too_short_lengths.txt'),
        ('trimmed_lengths', array_1d, '{name}_trimmed_lengths.txt'),
        ('tRNA_lengths', array_1d, '{name}_tRNA_lengths.txt'),
        ('rRNA_lengths', array_1d, '{name}_rRNA_lengths.txt'),
        ('phiX_lengths', array_1d, '{name}_phiX_lengths.txt'),
        ('synthetic_lengths', array_1d, '{name}_synthetic_lengths.txt'),
        ('other_ncRNA_lengths', array_1d, '{name}_ncRNA_lengths.txt'),
        ('clean_lengths', array_1d, '{name}_clean_lengths.txt'),
        ('clean_trimmed_lengths', array_1d, '{name}_clean_trimmed_lengths.txt'),
        ('remapped_lengths', array_1d, '{name}_remapped_lengths.txt'),
        ('merged_mapping_lengths', array_1d, '{name}_merged_mapping_lengths.txt'),
        ('unmapped_lengths', array_1d, '{name}_unmapped_lengths.txt'),
        ('long_polyA_lengths', array_1d, '{name}_long_polyA_lengths.txt'),
        ('unambiguous_lengths', array_1d, '{name}_unambiguous_lengths.txt'),
        ('oligo_hit_lengths', array_2d, '{name}_oligo_hit_lengths.txt'),

        ('rRNA_bam', 'bam', '{name}_rRNA.bam'),
        ('clean_bam', 'bam', '{name}_clean.bam'),
        ('more_rRNA_bam', 'bam', '{name}_more_rRNA.bam'),
        ('tRNA_bam', 'bam', '{name}_tRNA.bam'),
        ('other_ncRNA_bam', 'bam', '{name}_other_ncRNA.bam'),
        ('unambiguous_bam', 'bam', '{name}_unambiguous.bam'),

        ('codons_to_examine', counts, '{name}_codons_to_examine.txt'),

        ('phiX_bam', 'bam', '{name}_phiX.bam'),

        ('clean_trimmed_bam', 'bam', '{name}_clean_trimmed.bam'),
        ('merged_mappings', 'bam', '{name}_merged_mappings.bam'),

        ('unmapped_structures', None, '{name}_unmapped_structures.txt'),
        ('common_unmapped', counts, '{name}_common_unmapped.txt'),

        ('mismatches', mismatches, '{name}_mismatches.npy'),

        ('read_positions', read_positions, '{name}_read_positions.hdf5'),
        
        ('codon_counts', codon_counts, '{name}_codon_counts.txt'),
        ('codon_counts_anisomycin', codon_counts, '{name}_codon_counts_anisomycin.txt'),
        ('codon_counts_stringent', codon_counts, '{name}_codon_counts_stringent.txt'),
        ('buffered_codon_counts', read_positions, '{name}_buffered_codon_counts.hdf5'),
        ('metacodon_counts', read_positions, '{name}_metacodon_counts.hdf5'),
        ('metanucleotide_counts', read_positions, '{name}_metanucleotide_counts.hdf5'),
        ('mean_densities', read_positions, '{name}_mean_densities.hdf5'),
        ('mean_densities_anisomycin', read_positions, '{name}_mean_densities_anisomycin.hdf5'),
        ('mean_densities_no_misannotated', read_positions, '{name}_mean_densities_no_misannotated.hdf5'),
        ('RPKMs', RPKMs, '{name}_RPKMs.txt'),
        ('RPKMs_exclude_edges', RPKMs, '{name}_RPKMs_exclude_edges.txt'),
        ('read_counts', expression, '{name}_read_counts.txt'),
        ('read_counts_exclude_edges', expression, '{name}_read_counts_exclude_edges.txt'),

        ('rRNA_coverage', coverage, '{name}_rRNA_coverage.hdf5'),
        ('dominant_stretches', 'PH', '{name}_dominant_stretches.txt'),

        ('unmapped_trimmed_fastq', 'fastq', '{name}_unmapped_trimmed.fastq'),
        ('tophat_remapped_polyA_dir', 'dir', 'tophat_remapped_polyA'),
        ('remapped_accepted_hits', 'bam', 'tophat_remapped_polyA/accepted_hits.bam'),
        ('remapped_unmapped_bam', 'bam', 'tophat_remapped_polyA/unmapped.bam'),
        ('remapped_extended', 'bam', '{name}_remapped_extended.bam'),
    
        ('recycling_ratios', 'ratios', '{name}_recycling_ratios.txt'),

        ('reciprocal_rates', 'pickle', '{name}_reciprocal_rates.pkl'),
        ('reciprocal_rates_exclude_50', 'pickle', '{name}_reciprocal_rates_exclude_50.pkl'),

        ('yield', '', '{name}_yield.txt'),
    ]

    specific_figure_files = [
        ('quality', '{name}_quality.png'),
        ('clean_composition', '{name}_clean_composition.pdf'),
        ('clean_composition_perfect', '{name}_clean_composition_perfect.pdf'),
        ('unmapped_composition', '{name}_unmapped_composition.pdf'),
        ('all_lengths', '{name}_all_lengths.pdf'),
        ('clean_lengths', '{name}_clean_lengths.pdf'),
        ('rRNA_coverage', '{name}_rRNA_coverage.pdf'),
        ('oligo_hit_lengths', '{name}_oligo_hit_lengths.pdf'),
        ('dominant_stretch_lengths', '{name}_dominant_stretch_lengths.pdf'),
        ('starts_and_ends', '{name}_starts_and_ends.pdf'),
        ('three_prime_starts_and_ends', '{name}_three_prime_starts_and_ends.pdf'),
        ('ends', '{name}_ends.pdf'),
        ('starts_and_ends_heatmap', '{name}_starts_and_ends_heatmap.pdf'),
        ('starts_and_ends_heatmap_zoomed_out', '{name}_starts_and_ends_heatmap_zoomed_out.pdf'),
        ('mean_densities', '{name}_mean_densities.pdf'),
        ('mean_nucleotide_densities', '{name}_mean_nucleotide_densities.pdf'),
        ('mean_densities_anisomycin', '{name}_mean_densities_anisomycin.pdf'),
        ('mismatches', '{name}_mismatches.pdf'),
        ('frames', '{name}_frames.pdf'),
        ('unambiguous_frames', '{name}_unambiguous_frames.pdf'),
        ('metacodon_counts', '{name}_metacodon_counts.pdf'),
        ('metacodon_counts_stringent', '{name}_metacodon_counts_stringent.pdf'),
        ('metacodon_mean_enrichments', '{name}_metacodon_mean_enrichments.pdf'),
        ('metacodon_counts_nucleotide_resolution', '{name}_metacodon_counts_nucleotide_resolution.pdf'),
        ('metanucleotide_counts', '{name}_metanucleotide_counts.pdf'),
        ('codon_enrichments', '{name}_codon_enrichments.pdf'),
    ]

    specific_outputs = [
        ['clean_composition',
         'clean_composition_perfect',
         'unmapped_composition',
         'too_short_lengths',
         'trimmed_lengths',
         'tRNA_lengths',
         'rRNA_lengths',
         'phiX_lengths',
         'synthetic_lengths',
         'other_ncRNA_lengths',
         'clean_lengths',
         'clean_trimmed_lengths',
         'remapped_lengths',
         'merged_mapping_lengths',
         'unmapped_lengths',
         'long_polyA_lengths',
         'rRNA_coverage',
         'common_unmapped',
         'merged_mappings',
         'rRNA_bam',
         'more_rRNA_bam',
         'tRNA_bam',
         'other_ncRNA_bam',
         'mismatches',
         'codons_to_examine',
        ],
        ['read_positions',
         'metagene_positions',
         'buffered_codon_counts',
         'codon_counts',
         'read_counts',
         'read_counts_exclude_edges',
        ],
    ]
    
    specific_work = [
        ['preprocess',
         'map_full_lengths',
         'process_initially_unmapped',
         'merge_mapping_pathways',
         'process_remapped_unmapped',
         'compute_base_composition',
         'find_unambiguous_lengths',
         'get_rRNA_coverage',
         'examine_locii',
        ],
        ['get_read_positions',
         'get_metagene_positions',
         'compute_total_read_counts',
         'compute_codon_occupancy_counts',
        ],
    ]

    specific_cleanup = [
        ['compute_yield',
         'plot_base_composition',
         'plot_lengths',
         'plot_rRNA_coverage',
         'plot_mismatches',
         'visualize_unmapped',
        ],
        ['compute_RPKMs',
         'compute_mean_densities',
         'plot_starts_and_ends',
         'plot_frames',
         'plot_pausing',
        ],
    ]

    def __init__(self, **kwargs):
        super(RibosomeProfilingExperiment, self).__init__(**kwargs)

        self.adapter_type = kwargs['adapter_type']
        self.possibly_misannotated_file_name = kwargs.get('possibly_misannotated_file_name', None)
        
        self.phiX_bowtie2_index_prefix = kwargs['phiX_bowtie2_index_prefix']
        
        # A fasta file of synthetic sequences (markers, adapters) that need to
        # be filtered out can be provided. 
        self.synthetic_fasta = kwargs.get('synthetic_fasta', None)
        
        # Which mapping of length to A-site offset to use. Varies by experiment
        # because of different digestion and library preparation strategies.
        self.offset_type = kwargs['offset_type']

        # Which codon table to use when checking that claimed coding sequences
        # are valid. E. coli should have this set to 11.
        self.codon_table = kwargs.get('codon_table', 1)
        
        length_range = kwargs.get('relevant_lengths', None)
        if length_range == None:
            self.relevant_lengths = range(27, 32)
        else:
            start, stop = map(int, length_range.split(','))
            self.relevant_lengths = range(start, min(stop, self.max_read_length) + 1)
        
        self.max_interesting_length = int(kwargs.get('max_interesting_length', 51))
        
        #if self.adapter_type == 'polyA':
        #    specific_outputs[0].extend(['unambiguous_lengths',
        #                                'unambiguous_bam',
        #                               ])
        #    specific_outputs[1].extend(['unambiguous_read_positions',
        #                                'unambiguous_from_starts'
        #                               ])

        self.trim_function = trim.bound_trim[self.adapter_type]

        self.codons_to_examine = []

        locii = kwargs.get('codons_to_examine')
        if locii:
            locii = locii.split(';')
            for locus in locii:
                gene_name, codon_number = locus.split(',')
                self.codons_to_examine.append((gene_name, int(codon_number)))

    def compute_base_composition(self):
        seq_info_pairs = composition.get_seq_info_pairs(self.file_names['clean_bam'])
        all_array, perfect_array = composition.length_stratified_composition(seq_info_pairs, self.max_read_length)
        
        self.write_file('clean_composition', all_array)
        self.write_file('clean_composition_perfect', perfect_array)

    def quality_distribution(self):
        q_array, c_array, _ = fastq.quality_and_complexity(self.get_reads(), self.max_read_length)
        self.write_file('quality', q_array)

    def plot_base_composition(self):
        all_array = self.read_file('clean_composition', merged=True)
        composition.length_stratified_plot(all_array,
                                           self.figure_file_names['clean_composition'],
                                           lengths=self.relevant_lengths + ['all', 'all_from_right'],
                                          )
        
        perfect_array = self.read_file('clean_composition_perfect', merged=True)
        composition.length_stratified_plot(perfect_array,
                                           self.figure_file_names['clean_composition_perfect'],
                                           lengths=self.relevant_lengths + ['all', 'all_from_right'],
                                          )
        
        unmapped_array = self.read_file('unmapped_composition', merged=True)
        composition.length_stratified_plot(unmapped_array,
                                           self.figure_file_names['unmapped_composition'],
                                           lengths=self.relevant_lengths + ['all', 'all_from_right'],
                                          )

    def examine_locii(self):
        locii = {}
        
        CDSs, _ = self.get_CDSs(force_all=True)
        CDSs = {c.name: c for c in CDSs}

        for gene_name, codon_number in self.codons_to_examine:
            gene = CDSs[gene_name]
            reads = fastq.reads(self.file_names['preprocessed_reads'])
            triplets = examine_specific_codon.count_triplets(reads, gene, codon_number)

            locii[gene_name, codon_number] = triplets

        self.write_file('codons_to_examine', locii)

    def preprocess(self):
        reads = self.get_reads()
        trimmed_reads = self.trim_reads(reads)

        rRNA_filtered_reads = contaminants.pre_filter(self.file_names['rRNA_index'],
                                                      trimmed_reads,
                                                      self.file_names['rRNA_bam'],
                                                      self.file_names['bowtie_error'],
                                                     )

        with open(self.file_names['preprocessed_reads'], 'w') as preprocessed_fh:
            for read in rRNA_filtered_reads:
                preprocessed_fh.write(str(read))

    def filter_synthetic_sequences(self, reads):
        if self.synthetic_fasta:
            synthetic_sequences = [read.seq for read in fasta.reads(self.synthetic_fasta)]
        else:
            synthetic_sequences = []

        synthetic_lengths = np.zeros(self.max_read_length + 1)
        for read in reads:
            if contaminants.is_synthetic(read, synthetic_sequences):
                synthetic_lengths[len(read.seq)] += 1
            else:
                yield read

        self.write_file('synthetic_lengths', synthetic_lengths)

    def map_full_lengths(self):
        self.map_tophat()
        self.post_filter_contaminants()

    def process_initially_unmapped(self):
        unmapped_reads = sam.bam_to_fastq(self.file_names['unmapped_bam'])
        no_phiX_reads = self.filter_phiX(unmapped_reads)
        self.remap_polyA_trimmed(no_phiX_reads)

    def filter_phiX(self, reads):
        filtered_reads = contaminants.pre_filter(self.phiX_bowtie2_index_prefix,
                                                 reads,
                                                 self.file_names['phiX_bam'],
                                                )
        for read in filtered_reads:
            yield read

        phiX_length_counts = sam.get_length_counts(self.file_names['phiX_bam'])
        phiX_lengths = self.zero_padded_array(phiX_length_counts)
        self.write_file('phiX_lengths', phiX_lengths)

    def map_tophat(self):
        mapping_tools.map_tophat([self.file_names['preprocessed_reads']],
                                 self.file_names['bowtie2_index_prefix'],
                                 self.file_names['genes'],
                                 self.file_names['transcriptome_index'],
                                 self.file_names['tophat_dir'],
                                )
    
    def process_full_length_mappings(self):
        clean_bam = pysam.Samfile(self.file_names['clean_bam'])
        
        type_shape = (self.max_read_length + 1,
                      self.max_read_length,
                      fastq.MAX_EXPECTED_QUAL + 1,
                      6,
                      6,
                     )
        type_counts = np.zeros(type_shape, int)

        # To avoid counting mismatches in non-unique mappings multiple times,
        # a dummy secondary_type_counts array is passed to
        # trim_mismatches_from_start for secondary mappings.
        secondary_type_counts = np.zeros(type_shape, int)
        
        clean_trimmed_length_counts = Counter()
    
        region_fetcher = genomes.build_region_fetcher(self.file_names['genome'],
                                                      load_references=True,
                                                      sam_file=clean_bam,
                                                     )

        for mapping in clean_bam:
            if mapping.is_secondary:
                counts_array = secondary_type_counts
            else:
                counts_array = type_counts

            trimmed_from_start = trim.trim_mismatches_from_start(mapping,
                                                                 region_fetcher,
                                                                 counts_array,
                                                                )
            trimmed_from_end = trim.trim_nongenomic_polyA_from_end(trimmed_from_start,
                                                                   region_fetcher,
                                                                  )
            if not trimmed_from_end.is_unmapped and not trimmed_from_end.is_secondary:
                clean_trimmed_length_counts[trimmed_from_end.qlen] += 1

            yield trimmed_from_end

        self.write_file('mismatches', type_counts)
        
        clean_trimmed_lengths = self.zero_padded_array(clean_trimmed_length_counts)
        self.write_file('clean_trimmed_lengths', clean_trimmed_lengths)

    def remap_polyA_trimmed(self, reads):
        trim.trim_polyA_from_unmapped(reads,
                                      self.file_names['unmapped_trimmed_fastq'],
                                      second_time=True,
                                     )
        if os.path.getsize(self.file_names['unmapped_trimmed_fastq']) == 0:
            # tophat crashes on an empty file, so create empty files of all the
            # output of tophat that we need to exist.
            if not os.path.isdir(self.file_names['tophat_remapped_polyA_dir']):
                os.mkdir(self.file_names['tophat_remapped_polyA_dir'])
            template = pysam.Samfile(self.file_names['accepted_hits'], 'rb')
            empty_accepted_hits = pysam.Samfile(self.file_names['remapped_accepted_hits'], 'wb', template=template)
            empty_accepted_hits.close()
            empty_unmapped = pysam.Samfile(self.file_names['remapped_unmapped_bam'], 'wb', template=template)
            empty_unmapped.close()
        else:
            mapping_tools.map_tophat([self.file_names['unmapped_trimmed_fastq']],
                                     self.file_names['bowtie2_index_prefix'],
                                     self.file_names['genes'],
                                     self.file_names['transcriptome_index'],
                                     self.file_names['tophat_remapped_polyA_dir'],
                                    )

    def process_remapped(self):
        clean_bam = pysam.Samfile(self.file_names['remapped_accepted_hits'])
        
        type_shape = (self.max_read_length + 1,
                      self.max_read_length,
                      fastq.MAX_EXPECTED_QUAL + 1,
                      6,
                      6,
                     )
        type_counts = np.zeros(type_shape, int)
        remapped_length_counts = Counter()
    
        region_fetcher = genomes.build_region_fetcher(self.file_names['genome'],
                                                      load_references=True,
                                                      sam_file=clean_bam,
                                                     )

        for mapping in clean_bam:
            trimmed_from_start = trim.trim_mismatches_from_start(mapping,
                                                                 region_fetcher,
                                                                 type_counts,
                                                                )
            # Add back any genomic A's that were trimmed as part of mappings and
            # any remaining A's from the first non-genomic onward as soft clipped
            # bases for visualization in IGV.
            extended = trim.extend_polyA_end(trimmed_from_start,
                                             region_fetcher,
                                             trimmed_twice=True,
                                            )
            if not extended.is_unmapped and not extended.is_secondary:
                remapped_length_counts[extended.qlen] += 1

            yield extended

        remapped_lengths = self.zero_padded_array(remapped_length_counts)
        self.write_file('remapped_lengths', remapped_lengths)
    
    def process_remapped_unmapped(self):
        unmapped_lengths = np.zeros(self.max_read_length + 1)
        unmapped_seq_counts = Counter()

        long_polyA_lengths = np.zeros(self.max_read_length + 1)
        long_polyA_counts = Counter()
        
        unmapped_reads = sam.bam_to_fastq(self.file_names['remapped_unmapped_bam'])
        unretrimmed_reads = trim.untrim_reads(unmapped_reads, second_time=True)
        synthetic_filtered_reads = self.filter_synthetic_sequences(unretrimmed_reads)

        def record_common(reads):
            for read in reads:
                if predominantly_A(read.seq):
                    long_polyA_lengths[len(read.seq)] += 1
                    long_polyA_counts[read.seq] += 1
                else:
                    unmapped_lengths[len(read.seq)] += 1
                    unmapped_seq_counts[read.seq] += 1
                yield read

        synthetic_filtered_reads = record_common(synthetic_filtered_reads)
        
        seq_info_pairs = ((read.seq, False) for read in synthetic_filtered_reads)
        all_array, _ = composition.length_stratified_composition(seq_info_pairs, self.max_read_length)
        
        self.write_file('unmapped_composition', all_array)
        

        self.write_file('unmapped_lengths', unmapped_lengths)
        self.write_file('long_polyA_lengths', long_polyA_lengths)
        
        non_long_polyA = Counter(dict(unmapped_seq_counts.most_common(100)))
        long_polyA = Counter(dict(long_polyA_counts.most_common(100)))
        
        common_unmapped = {'non_long_polyA': non_long_polyA,
                           'long_polyA': long_polyA,
                           }
        self.write_file('common_unmapped', common_unmapped)

    def visualize_unmapped(self):
        bowtie2_targets = [(self.file_names['genome'], self.file_names['bowtie2_index_prefix'], 'C,20,0'),
                          ]
        sw_genome_dirs = ['/home/jah/genomes/truseq']
        extra_targets = [fasta.Read('smRNA_linker', trim.smRNA_linker)]
        if self.synthetic_fasta:
            extra_targets.extend(list(fasta.reads(self.synthetic_fasta)))
        
        def get_reads():
            for i, (seq, count) in enumerate(self.read_file('common_unmapped')['non_long_polyA'].most_common()):
                read = fastq.Read('{0}_{1}'.format(i, count),
                                  seq,
                                  fastq.encode_sanger([40]*len(seq)),
                                 )
                yield read

        visualize_structure.visualize_unpaired_alignments(get_reads,
                                                          sw_genome_dirs,
                                                          extra_targets,
                                                          bowtie2_targets,
                                                          self.file_names['unmapped_structures'],
                                                         )

    def merge_mapping_pathways(self):
        sam_file = pysam.Samfile(self.file_names['clean_bam'])
        alignment_sorter = sam.AlignmentSorter(sam_file.references,
                                               sam_file.lengths,
                                               self.file_names['merged_mappings'],
                                              )
        with alignment_sorter:
            for trimmed in self.process_full_length_mappings():
                alignment_sorter.write(trimmed)
            for extended in self.process_remapped():
                alignment_sorter.write(extended)

        merged_mapping_lengths = self.read_file('clean_trimmed_lengths') + self.read_file('remapped_lengths')
        self.write_file('merged_mapping_lengths', merged_mapping_lengths)

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
        data = contaminants.produce_rRNA_coverage(self.file_names['rRNA_bam'], self.max_read_length)
        self.write_file('rRNA_coverage', data)
    
    def get_oligo_hit_lengths(self):
        lengths = contaminants.get_oligo_hit_lengths(self.file_names['rRNA_bam'],
                                                     self.file_names['oligos'],
                                                     self.file_names['oligos_sam'],
                                                     self.max_read_length,
                                                    )
        self.write_file('oligo_hit_lengths', lengths)

    def compute_yield(self):
        trimmed_lengths = self.read_file('trimmed_lengths', merged=True)
        phiX_lengths = self.read_file('phiX_lengths', merged=True)
        too_short_lengths = self.read_file('too_short_lengths', merged=True)
        rRNA_lengths = self.read_file('rRNA_lengths', merged=True)
        synthetic_lengths = self.read_file('synthetic_lengths', merged=True)
        tRNA_lengths = self.read_file('tRNA_lengths', merged=True)
        other_ncRNA_lengths = self.read_file('other_ncRNA_lengths', merged=True)
        unmapped_lengths = self.read_file('unmapped_lengths', merged=True)
        long_polyA_lengths = self.read_file('long_polyA_lengths', merged=True)
        clean_lengths = self.read_file('clean_lengths', merged=True)
        remapped_lengths = self.read_file('remapped_lengths', merged=True)

        total_reads = trimmed_lengths.sum() + too_short_lengths.sum()
        phiX_reads = phiX_lengths.sum()
        long_enough_reads = trimmed_lengths.sum()
        rRNA_reads = rRNA_lengths.sum()
        synthetic_reads = synthetic_lengths.sum()
        tRNA_reads = tRNA_lengths.sum()
        other_ncRNA_reads = other_ncRNA_lengths.sum()
        unmapped_reads = unmapped_lengths.sum()
        clean_reads = clean_lengths.sum()
        remapped_reads = remapped_lengths.sum()
        long_polyA_reads = long_polyA_lengths.sum()

        dominant_reads, boundaries = contaminants.identify_dominant_stretches(self.read_file('rRNA_coverage', merged=True),
                                                                              total_reads,
                                                                              self.max_read_length,
                                                                              self.merged_file_names['rRNA_bam'],
                                                                             )
        contaminants.plot_dominant_stretch_lengths(boundaries, self.figure_file_names['dominant_stretch_lengths'])
        other_reads = rRNA_reads - dominant_reads

        with open(self.merged_file_names['dominant_stretches'], 'w') as dominant_stretches_file:
            for rname in sorted(boundaries):
                for start, stop in boundaries[rname]:
                    dominant_stretches_file.write('{0}:{1}-{2}\n'.format(rname, start, stop))

        with open(self.file_names['yield'], 'w') as yield_file:
            yield_file.write('Total reads: {0:,}\n'.format(total_reads))
            for category, count in [('Long enough reads', long_enough_reads),
                                    ('phiX reads', phiX_reads),
                                    ('rRNA reads', rRNA_reads),
                                    ('(rRNA reads from non-dominant stetches)', other_reads),
                                    ('tRNA reads', tRNA_reads),
                                    ('Other ncRNA reads', other_ncRNA_reads),
                                    ('Clean reads', clean_reads),
                                    ('Reads mapped after polyA trimming', remapped_reads),
                                    ('Reads that start wth long polyA', long_polyA_reads),
                                    ('Synthetic reads', synthetic_reads),
                                    ('Unaccounted-for reads', unmapped_reads),
                                   ]:
                fraction = float(count) / total_reads
                line = '{0}: {1:,} ({2:.2%})\n'.format(category,
                                                       count,
                                                       fraction,
                                                      )
                yield_file.write(line)

    def plot_lengths(self):
        lengths = {'too short': (self.read_file('too_short_lengths'), 'purple'),
                   'rRNA': (self.read_file('rRNA_lengths'), 'red'),
                   'phiX': (self.read_file('phiX_lengths'), 'brown'),
                   'synthetic': (self.read_file('synthetic_lengths'), 'orange'),
                   'tRNA': (self.read_file('tRNA_lengths'), 'blue'),
                   'other ncRNA': (self.read_file('other_ncRNA_lengths'), 'black'),
                   'unmapped': (self.read_file('unmapped_lengths'), 'cyan'),
                   'long_polyA': (self.read_file('long_polyA_lengths'), 'yellow'),
                   'clean': (self.read_file('clean_lengths'), 'green'),
                  }

        fig_all, ax_all = plt.subplots(figsize=(12, 8))
    
        for key, (counts, color) in lengths.items():
            if len(counts) != self.max_read_length + 1:
                raise ValueError

            if self.max_read_length > self.max_interesting_length:
                counts[self.max_interesting_length] = counts[self.max_interesting_length:].sum()
                counts[self.max_interesting_length + 1:] = 0
            ax_all.plot(counts, '.-', label=key, color=color)

        ax_all.set_xlim(0, self.max_interesting_length)
        ax_all.set_title('{0}\n{1}\nFragment length distribution by source'.format(self.group, self.name))
        ax_all.set_xlabel('Length of original RNA fragment')
        ax_all.set_ylabel('Number of reads')
        ax_all.legend(loc='upper left', framealpha=0.5)
        Sequencing.Visualize.add_commas_to_yticks(ax_all)
        
        fig_all.savefig(self.figure_file_names['all_lengths'])
        plt.close(fig_all)
        
        fig_clean, ax_clean = plt.subplots(figsize=(12, 8))
        
        lengths['clean_trimmed'] = (self.read_file('clean_trimmed_lengths'), 'blue')
        lengths['remapped'] = (self.read_file('remapped_lengths'), 'red')

        for key in ('clean', 'clean_trimmed', 'remapped'):
            if key in lengths:
                counts, color = lengths[key]
                normalized_counts = np.true_divide(counts, counts.sum())
                ax_clean.plot(normalized_counts, '.-', label=key, color=color)
        
        ax_clean.axvspan(27.5, 28.5, color='green', alpha=0.2)
        ax_clean.set_xlim(0, self.max_interesting_length)
        ax_clean.set_title('{0}\n{1}\nFragment length distribution by source'.format(self.group, self.name))
        ax_clean.set_xlabel('Length of original RNA fragment')
        ax_clean.set_ylabel('Fraction of clean reads')
        ax_clean.legend(loc='upper left', framealpha=0.5)
    
        fig_clean.savefig(self.figure_file_names['clean_lengths'])
        plt.close(fig_clean)

    def get_total_reads(self):
        trimmed_lengths = self.read_file('trimmed_lengths')
        too_short_lengths = self.read_file('too_short_lengths')
        total_reads = trimmed_lengths.sum() + too_short_lengths.sum()
        return total_reads

    def plot_rRNA_coverage(self):
        counts = self.read_file('rRNA_coverage')
        coverage_data = {self.name: (self.get_total_reads(), counts, 'blue')}
        contaminants.plot_rRNA_coverage(coverage_data,
                                        self.file_names['oligos_sam'],
                                        self.figure_file_names['rRNA_coverage'],
                                       )
        
    def plot_oligo_hit_lengths(self):
        lengths = self.read_file('oligo_hit_lengths')
        contaminants.plot_oligo_hit_lengths(self.file_names['oligos'],
                                            lengths,
                                            self.figure_file_names['oligo_hit_lengths'],
                                           )

    def plot_starts_and_ends(self):
        short_lengths = range(20, 24)
        long_lengths = range(27, 32)

        metagene_positions = self.read_file('metagene_positions')

        visualize.plot_metagene_positions(metagene_positions,
                                          self.figure_file_names['starts_and_ends'],
                                          long_lengths + ['all', 'all_nonunique'],
                                          title='{0}\n{1}'.format(self.group, self.name),
                                         )
        
        visualize.plot_metagene_positions(metagene_positions,
                                          self.figure_file_names['three_prime_starts_and_ends'],
                                          ['three_prime_genomic', 'three_prime_nongenomic', 'three_prime_nonunique'],
                                          title='{0}\n{1}'.format(self.group, self.name),
                                         )
        
        visualize.plot_just_ends(metagene_positions['stop_codon'],
                                 self.figure_file_names['ends'],
                                 title='{0}\n{1}'.format(self.group, self.name),
                                )

        visualize.plot_metagene_positions_heatmap(metagene_positions,
                                                  self.figure_file_names['starts_and_ends_heatmap'],
                                                  normalize_to_max_in=(range(30, 34), 'start', slice(1, 1000)),
                                                 )
        visualize.plot_metagene_positions_heatmap(metagene_positions,
                                                  self.figure_file_names['starts_and_ends_heatmap_zoomed_out'],
                                                  zoomed_out=True,
                                                  normalize_to_max_in=(range(30, 34), 'start', slice(1, 1000)),
                                                 )

        visualize.plot_averaged_codon_densities([(self.name, self.read_file('mean_densities'), 0)],
                                                self.figure_file_names['mean_densities'],
                                                past_edge=10,
                                                plot_up_to=1000,
                                                smooth=False,
                                               )
        
        visualize.plot_averaged_nucleotide_densities([(self.name, metagene_positions, 0)],
                                                     self.figure_file_names['mean_nucleotide_densities'],
                                                     past_edge=10,
                                                     plot_up_to=2000,
                                                     smooth=False,
                                                    )

        #visualize.plot_averaged_codon_densities([(self.name, self.read_file('mean_densities_anisomycin'), 0)],
        #                                        self.figure_file_names['mean_densities_anisomycin'],
        #                                        past_edge=10,
        #                                        plot_up_to=100,
        #                                        smooth=False,
        #                                       )

    def plot_frames(self):
        metagene_positions = self.read_file('metagene_positions')

        visualize.plot_frames(metagene_positions['start_codon'],
                              self.figure_file_names['frames'],
                             )

        #if self.adapter_type == 'polyA':
        #    from_starts = self.read_file('unambiguous_from_starts')

        #    visualize.plot_frames(from_starts['from_starts'],
        #                          self.figure_file_names['unambiguous_frames'],
        #                         )

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

    def plot_pausing(self):
        relevant_experiments = [self]

        allowed_at_pause = set(codons.non_stop_codons)
        not_allowed_at_stall = set()

        codon_counts_dict = {}
        around_lists_dict = {}
        stratified_mean_enrichments_dict = {}

        for experiment in relevant_experiments:
            name = experiment.name
            codon_counts_dict[name] = experiment.read_file('buffered_codon_counts',
                                                           specific_keys={'relaxed', 'identities'},
                                                          )
            around_lists_dict[name] = pausing.metacodon_around_pauses(codon_counts_dict[name],
                                                                      allowed_at_pause,
                                                                      not_allowed_at_stall,
                                                                     )
            stratified_mean_enrichments_dict[name] = pausing.compute_stratified_mean_enrichments(around_lists_dict[name])
            del around_lists_dict[name]

        with PdfPages(self.figure_file_names['codon_enrichments']) as pdf:
            for amino_acid in codons.full_back_table:
                if amino_acid == '*':
                    continue

                fig = pausing.plot_codon_enrichments(relevant_experiments,
                                                     stratified_mean_enrichments_dict,
                                                     amino_acid,
                                                     min_x=-60,
                                                     max_x=60,
                                                    )
                pdf.savefig(figure=fig, bbox_inches='tight')
                plt.close(fig)

    def plot_mismatches(self):
        type_counts = self.read_file('mismatches')
        visualize.plot_mismatch_type_by_position(type_counts,
                                                 self.relevant_lengths,
                                                 self.figure_file_names['mismatches'],
                                                )
            
    def zero_padded_array(self, counts):
        array = utilities.counts_to_array(counts)
        if len(array) < self.max_read_length + 1:
            padded_array = np.zeros(self.max_read_length + 1, int)
            padded_array[:len(array)] += array
        else:
            padded_array = array
        return padded_array

    def get_read_positions(self):
        piece_CDSs, max_gene_length = self.get_CDSs()
        gene_infos = positions.get_Transcript_position_counts(self.merged_file_names['merged_mappings'],
                                                              piece_CDSs,
                                                              self.relevant_lengths,
                                                             )

        self.read_positions = {}
        for name, info in gene_infos.iteritems():
            five_prime_counts = info['five_prime_positions']
            three_prime_counts = info['three_prime_positions']

            all_positions = {'three_prime_genomic': three_prime_counts[0],
                             'three_prime_nongenomic': three_prime_counts['all'] - three_prime_counts[0],
                             'three_prime_nonunique': three_prime_counts['all_nonunique'],
                             'sequence': info['sequence'],
                            }
            all_positions.update(five_prime_counts)

            self.read_positions[name] = all_positions

        self.write_file('read_positions', self.read_positions)
        
    def get_metagene_positions(self):
        piece_CDSs, max_gene_length = self.get_CDSs()
        read_positions = self.load_read_positions()
        metagene_positions = positions.compute_metagene_positions(piece_CDSs,
                                                                  read_positions,
                                                                  max_gene_length,
                                                                 )

        self.write_file('metagene_positions', metagene_positions)

    def compute_codon_occupancy_counts(self):
        read_positions = self.load_read_positions()

        buffered_codon_counts = {}
        codon_counts = {}
        codon_counts_stringent = {}
        codon_counts_anisomycin = {}
        for name, position_counts in read_positions.iteritems():
            buffered_counts, identities = positions.compute_codon_counts(position_counts, self.offset_type)
            buffered_counts_stringent, _ = positions.compute_codon_counts(position_counts, self.offset_type + '_stringent')
            buffered_counts_anisomycin, _ = positions.compute_codon_counts(position_counts, self.offset_type + '_anisomycin')
            buffered_codon_counts[name] = {'relaxed': buffered_counts,
                                           'stringent': buffered_counts_stringent,
                                           'anisomycin': buffered_counts_anisomycin,
                                           'identities': identities,
                                          }
            num_codons = buffered_counts.CDS_length

            # + 1 is to include the stop codon
            codon_counts[name] = buffered_counts['start_codon', :num_codons + 1]
            codon_counts_stringent[name] = buffered_counts_stringent['start_codon', :num_codons + 1]
            codon_counts_anisomycin[name] = buffered_counts_anisomycin['start_codon', :num_codons + 1]

        self.write_file('buffered_codon_counts', buffered_codon_counts)
        self.write_file('codon_counts', codon_counts)
        self.write_file('codon_counts_anisomycin', codon_counts_anisomycin) 
        self.write_file('codon_counts_stringent', codon_counts_stringent)

    def compute_total_read_counts(self):
        read_positions = self.load_read_positions()
        read_counts = positions.compute_read_counts(read_positions, 0, 0)
        self.write_file('read_counts', read_counts)
        read_counts_exclude_edges = positions.compute_read_counts(read_positions, 30, 4)
        self.write_file('read_counts_exclude_edges', read_counts_exclude_edges)

    def compute_metacodon_counts(self):
        codon_counts = self.read_file('buffered_codon_counts')
        metacodon_counts = positions.compute_metacodon_counts(codon_counts)
        self.write_file('metacodon_counts', metacodon_counts)

        read_positions = self.load_read_positions()
        metanucleotide_counts = positions.compute_metanucleotide_counts(read_positions)
        self.write_file('metanucleotide_counts', metanucleotide_counts)

    def compute_mean_densities(self):
        codon_counts = self.read_file('buffered_codon_counts', merged=True)
        mean_densities = positions.compute_averaged_codon_densities(codon_counts)
        self.write_file('mean_densities', mean_densities)
        mean_densities_anisomycin = positions.compute_averaged_codon_densities(codon_counts, offset_key='anisomycin')
        self.write_file('mean_densities_anisomycin', mean_densities_anisomycin)

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
