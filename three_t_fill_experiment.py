from Sequencing import adapters, utilities, fastq, mapping_tools, sam
from Sequencing.Parallel import map_reduce
from Sequencing.Serialize import array_1d, counts
from Serialize import read_positions
from Circles import paired_end
import rna_experiment
import os
import pysam
from collections import Counter
from itertools import izip
import contaminants
import trim
import positions
import visualize

class ThreeTFillExperiment(rna_experiment.RNAExperiment):
    num_stages = 2

    specific_results_files = [
        ('R1_preprocessed', 'fastq', '{name}_R1_preprocessed.fastq'),
        ('R2_preprocessed', 'fastq', '{name}_R2_preprocessed.fastq'),

        ('rRNA_bam', 'bam', '{name}_rRNA.bam'),

        ('barcodes', counts, '{name}_barcodes.txt'),

        ('R1_non_rRNA_flipped', 'fastq', '{name}_R1_non_rRNA_flipped.fastq'),

        ('R1_tophat_dir', 'dir', 'tophat_R1'),
        ('R1_accepted_hits', 'bam_by_name', 'tophat_R1/accepted_hits.bam'),
        ('R1_unmapped', 'bam_by_name', 'tophat_R1/unmapped.bam'),

        ('R2_tophat_dir', 'dir', 'tophat_R2'),
        ('R2_accepted_hits', 'bam_by_name', 'tophat_R2/accepted_hits.bam'),
        ('R2_unmapped', 'bam_by_name', 'tophat_R2/unmapped.bam'),
        
        ('combined', 'bam', '{name}_combined.bam'),

        ('trimmed_lengths', array_1d, '{name}_trimmed_lengths.txt'),
        ('tlens', array_1d, '{name}_tlens.txt'),
        
        ('three_prime_read_positions', read_positions, '{name}_three_prime_read_positions.hdf5'),
    ]

    specific_figure_files = []

    specific_outputs = [
        ['barcodes',
         'trimmed_lengths',
         'tlens',
         'combined',
        ],
        ['from_starts_and_ends',
         'three_prime_read_positions',
        ],
    ]

    specific_work = [
        ['preprocess',
         'map_tophat',
         'combine_mappings',
        ],
        ['get_polyA_positions',
         'get_metagene_positions',
        ],
    ]

    specific_cleanup = [
        [],
        ['plot_starts_and_ends',
        ],
    ]

    def __init__(self, **kwargs):
        super(ThreeTFillExperiment, self).__init__(**kwargs)

        self.barcode = kwargs['barcode']
        full_adapter_in_R1 = utilities.reverse_complement(self.barcode) + utilities.reverse_complement(adapters.primers['PE']['R2']) 
        self.adapter_in_R1 = full_adapter_in_R1[:19]
        
    def trim_reads(self, read_pairs):
        total_reads = 0
        long_enough_reads = 0
        trimmed_lengths = Counter()
        barcodes = Counter()
        
        for R1, R2 in read_pairs:
            total_reads += 1
            barcodes[R2.seq[:len(self.barcode)]] += 1

            # R2 isn't expected to have adapters sequence because it will
            # have to get through the A tail first.
            position = adapters.find_adapter(self.adapter_in_R1, 3, R1.seq)
            trimmed_lengths[position] += 1
            if position < 12:
                continue
            long_enough_reads += 1

            R1_slice = slice(None, position)
            # position points to where the barcode starts in R1. The length
            # of the trimmed R2 read should be equal to position.
            R2_slice = slice(len(self.barcode), len(self.barcode) + position)

            processed_R1 = fastq.Read(R1.name, R1.seq[R1_slice], R1.qual[R1_slice])
            processed_R2 = fastq.Read(R2.name, R2.seq[R2_slice], R2.qual[R2_slice])
            
            yield processed_R1, processed_R2

        trimmed_lengths = utilities.counts_to_array(trimmed_lengths)
        self.write_file('trimmed_lengths', trimmed_lengths)
        self.write_file('barcodes', barcodes)
        self.summary.extend(
            [('Total read pairs', total_reads),
             ('Long enough', long_enough_reads),
            ]
        )

    def preprocess(self):
        read_pairs = self.get_read_pairs()
        trimmed_pairs = self.trim_reads(read_pairs)
        filtered_pairs = self.filter_rRNA(trimmed_pairs)
        
        with open(self.file_names['R1_preprocessed'], 'w') as R1_fh, \
             open(self.file_names['R2_preprocessed'], 'w') as R2_fh:

            for R1, R2 in filtered_pairs:
                R1_fh.write(str(R1))
                R2_fh.write(str(R2))

    def filter_rRNA(self, read_pairs):
        filtered_pairs = contaminants.pre_filter_paired(self.file_names['rRNA_index'],
                                                        read_pairs,
                                                        self.file_names['rRNA_bam'],
                                                        '/dev/null',
                                                       )
        non_rRNA_reads = 0
        for pair in filtered_pairs:
            non_rRNA_reads += 1
            yield pair

        self.summary.extend(
            [('Non-rRNA reads', non_rRNA_reads),
            ]
        )

    def map_tophat(self):
        mapping_tools.map_tophat([self.file_names['R1_preprocessed']],
                                 self.file_names['bowtie2_index_prefix'],
                                 self.file_names['genes'],
                                 self.file_names['transcriptome_index'],
                                 self.file_names['R1_tophat_dir'],
                                 no_sort=True,
                                )

        mapping_tools.map_tophat([self.file_names['R2_preprocessed']],
                                 self.file_names['bowtie2_index_prefix'],
                                 self.file_names['genes'],
                                 self.file_names['transcriptome_index'],
                                 self.file_names['R2_tophat_dir'],
                                 no_sort=True,
                                )
    
    def combine_mappings(self):
        num_unmapped = 0
        num_R1_unmapped = 0
        num_R2_unmapped = 0
        num_nonunique = 0
        num_discordant = 0
        num_disoriented = 0
        num_concordant = 0

        tlens = Counter()

        R1_mappings = pysam.Samfile(self.file_names['R1_accepted_hits'])
        R1_unmapped = pysam.Samfile(self.file_names['R1_unmapped'])
        all_R1 = sam.merge_by_name(R1_mappings, R1_unmapped)
        R1_grouped = utilities.group_by(all_R1, lambda m: m.qname)

        R2_mappings = pysam.Samfile(self.file_names['R2_accepted_hits'])
        R2_unmapped = pysam.Samfile(self.file_names['R2_unmapped'])
        all_R2 = sam.merge_by_name(R2_mappings, R2_unmapped)
        R2_grouped = utilities.group_by(all_R2, lambda m: m.qname)

        group_pairs = izip(R1_grouped, R2_grouped)

        alignment_sorter = sam.AlignmentSorter(R1_mappings.references,
                                               R1_mappings.lengths,
                                               self.file_names['combined'],
                                              )

        with alignment_sorter:
            for (R1_qname, R1_group), (R2_qname, R2_group) in group_pairs:
                #print R1_qname, R2_qname
                if fastq.get_pair_name(R1_qname) != fastq.get_pair_name(R2_qname):
                    # Ensure that the iteration through pairs is in sync.
                    print R1_qname, R2_qname
                    raise ValueError
                
                R1_unmapped = any(m.is_unmapped for m in R1_group)
                R2_unmapped = any(m.is_unmapped for m in R2_group)
                if R1_unmapped:
                    num_R1_unmapped += 1
                if R2_unmapped:
                    num_R2_unmapped += 1
                if R1_unmapped or R2_unmapped:
                    num_unmapped += 1
                    continue

                R1_nonunique = len(R1_group) > 1 or any(m.mapq < 40 for m in R1_group)
                R2_nonunique = len(R2_group) > 1 or any(m.mapq < 40 for m in R2_group)
                if R1_nonunique or R2_nonunique:
                    num_nonunique += 1
                    continue
                
                R1_m = R1_group.pop()
                R2_m = R2_group.pop()

                R1_strand = sam.get_strand(R1_m)
                R2_strand = sam.get_strand(R2_m)

                tlen = max(R1_m.aend, R2_m.aend) - min(R1_m.pos, R2_m.pos)
                discordant = (R1_m.tid != R2_m.tid) or (R1_strand) == (R2_strand) or (tlen > 10000)
                if discordant:
                    num_discordant += 1
                    continue
                
                # Reminder: the protocol produces anti-sense reads.
                if R1_strand == '-':
                    if R1_m.pos < R2_m.pos:
                        num_disoriented += 1
                        continue

                elif R1_strand == '+':
                    if R2_m.pos < R1_m.pos:
                        num_disoriented += 1
                        continue
                
                combined_read = paired_end.combine_paired_mappings(R1_m, R2_m)
                
                tlens[tlen] += 1

                if combined_read:
                    # Flip combined_read back to the sense strand.
                    if combined_read.is_reverse:
                        combined_read.is_reverse = False
                    else:
                        combined_read.is_reverse = True

                    trim.set_nongenomic_length(combined_read, 0)
                    
                    alignment_sorter.write(combined_read)

                    num_concordant += 1

        self.summary.extend(
            [('Unmapped', num_unmapped),
             ('R1 unmapped', num_R1_unmapped),
             ('R2 unmapped', num_R2_unmapped),
             ('Nonunique', num_nonunique),
             ('Discordant', num_discordant),
             ('Unexpected orientation', num_disoriented),
             ('Concordant', num_concordant),
            ],
        )

        tlens = utilities.counts_to_array(tlens)
        self.write_file('tlens', tlens)

    def get_polyA_positions(self):
        piece_CDSs, max_gene_length = self.get_CDSs()
        gene_infos = positions.get_Transcript_polyA_position_counts(self.merged_file_names['combined'],
                                                                    piece_CDSs,
                                                                    max_nongenomic_length=-1,
                                                                    left_buffer=500,
                                                                    right_buffer=500,
                                                                   )

        self.three_prime_read_positions = {name: info['position_counts']
                                           for name, info in gene_infos.iteritems()}

        self.write_file('three_prime_read_positions', self.three_prime_read_positions)
    
    def get_metagene_positions(self):
        piece_CDSs, max_gene_length = self.get_CDSs()
        read_positions = self.load_read_positions(modifier='three_prime')
        
        processed_read_positions = {}
        for name in read_positions:
            gene = {'three_prime': read_positions[name]['all']}
            processed_read_positions[name] = gene
    
        from_starts_and_ends = positions.compute_metagene_positions(processed_read_positions, max_gene_length)

        self.write_file('from_starts_and_ends', from_starts_and_ends)
    
    
    def plot_starts_and_ends(self):
        from_starts_and_ends = self.read_file('from_starts_and_ends')

        visualize.plot_metagene_positions(from_starts_and_ends['from_starts'],
                                          from_starts_and_ends['from_ends'],
                                          self.figure_file_names['starts_and_ends_zoomed_out'],
                                          zoomed_out=True,
                                         )
    
    def get_total_eligible_reads(self):
        summary_pairs = self.read_file('summary')
        summary_dict = {name: values[0] for name, values in summary_pairs}
        total_mapped_reads = summary_dict['Nonunique'] + summary_dict['Concordant']
        return total_mapped_reads

if __name__ == '__main__':
    script_path = os.path.realpath(__file__)
    map_reduce.controller(ThreeTFillExperiment, script_path)
