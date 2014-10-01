from Sequencing import adapters, utilities, fastq, mapping_tools, sam
from Sequencing.Parallel import map_reduce
from Sequencing.Serialize import array_1d, counts
from Serialize import read_positions
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

        ('R1_non_rRNA', 'fastq', '{name}_R1_non_rRNA.fastq'),
        ('R2_non_rRNA', 'fastq', '{name}_R2_non_rRNA.fastq'),
        ('rRNA_sam', 'sam', '{name}_rRNA.sam'),
        ('rRNA_bam', 'bam', '{name}_rRNA.bam'),

        ('barcodes', counts, '{name}_barcodes.txt'),

        ('R1_non_rRNA_flipped', 'fastq', '{name}_R1_non_rRNA_flipped.fastq'),

        ('R1_tophat_dir', 'dir', 'tophat_R1'),
        ('R1_accepted_hits', 'bam', 'tophat_R1/accepted_hits.bam'),
        ('R1_unmapped', 'bam', 'tophat_R1/unmapped.bam'),
        ('R1_sorted_by_name', 'bam', '{name}_R1_by_name.bam'),

        ('R2_tophat_dir', 'dir', 'tophat_R2'),
        ('R2_accepted_hits', 'bam', 'tophat_R2/accepted_hits.bam'),
        ('R2_unmapped', 'bam', 'tophat_R2/unmapped.bam'),
        ('R2_sorted_by_name', 'bam', '{name}_R2_by_name.bam'),
        
        ('combined', 'bam', '{name}_combined.bam'),
        ('combined_sorted', 'bam', '{name}_combined_sorted.bam'),

        ('trimmed_lengths', array_1d, '{name}_trimmed_lengths.txt'),
        ('tlens', array_1d, '{name}_tlens.txt'),
        
        ('three_prime_read_positions', read_positions, '{name}_three_prime_read_positions.hdf5'),
    ]

    specific_figure_files = []

    specific_outputs = [
        [#'barcodes',
         'trimmed_lengths',
         'tlens',
         'combined_sorted',
        ],
        ['from_starts_and_ends',
         'three_prime_read_positions',
        ],
    ]

    specific_work = [
        ['preprocess',
         'filter_rRNA',
         'flip_R1',
         'map_tophat',
         'sort_by_name',
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
        full_adapter_in_R1 = utilities.reverse_complement(self.barcode) + adapters.paired_end_R2_rc 
        self.adapter_in_R1 = full_adapter_in_R1[:19]
        
    def preprocess(self):
        total_reads = 0
        long_enough_reads = 0
        trimmed_lengths = Counter()
        barcodes = Counter()
        read_pairs = self.get_read_pairs()
        
        with open(self.file_names['R1_preprocessed'], 'w') as R1_fh, \
             open(self.file_names['R2_preprocessed'], 'w') as R2_fh:

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

                processed_R1_record = fastq.make_record(R1.name,
                                                        R1.seq[R1_slice],
                                                        R1.qual[R1_slice],
                                                       )
                processed_R2_record = fastq.make_record(R2.name,
                                                        R2.seq[R2_slice],
                                                        R2.qual[R2_slice],
                                                       )
                R1_fh.write(processed_R1_record)
                R2_fh.write(processed_R2_record)

        trimmed_lengths = utilities.counts_to_array(trimmed_lengths)
        self.write_file('trimmed_lengths', trimmed_lengths)
        self.write_file('barcodes', barcodes)
        self.log.extend(
            [('Total read pairs', total_reads),
             ('Long enough', long_enough_reads),
            ]
        )

    def filter_rRNA(self):
        contaminants.pre_filter_paired(self.file_names['R1_preprocessed'],
                                       self.file_names['R2_preprocessed'],
                                       self.file_names['R1_non_rRNA'],
                                       self.file_names['rRNA_index'],
                                       self.file_names['rRNA_sam'],
                                       self.file_names['rRNA_bam'],
                                       '/dev/null',
                                      )
        non_rRNA_reads = sum(1 for r in fastq.reads(self.file_names['R1_non_rRNA']))
        self.log.extend(
            [('Non-rRNA reads', non_rRNA_reads),
            ]
        )

    def flip_R1(self):
        with open(self.file_names['R1_non_rRNA_flipped'], 'w') as R1_flipped_fh:
            for read in fastq.reads(self.file_names['R1_non_rRNA']):
                flipped_record = fastq.make_record(read.name,
                                                   utilities.reverse_complement(read.seq),
                                                   read.qual[::-1],
                                                  )
                R1_flipped_fh.write(flipped_record)

    def map_tophat(self):
        mapping_tools.map_tophat([self.file_names['R1_non_rRNA_flipped']],
                                 self.file_names['bowtie2_index_prefix'],
                                 self.file_names['genes'],
                                 self.file_names['transcriptome_index'],
                                 self.file_names['R1_tophat_dir'],
                                )
        mapping_tools.map_tophat([self.file_names['R2_non_rRNA']],
                                 self.file_names['bowtie2_index_prefix'],
                                 self.file_names['genes'],
                                 self.file_names['transcriptome_index'],
                                 self.file_names['R2_tophat_dir'],
                                )
    
    def sort_by_name(self):
        for edge in ['R1', 'R2']:
            accepted_fn = self.file_names['{}_accepted_hits'.format(edge)]
            accepted_prefix, _ = os.path.splitext(accepted_fn)
            accepted_sorted_fn = '{}_sorted.bam'.format(accepted_prefix)

            unmapped_fn = self.file_names['{}_unmapped'.format(edge)]
            unmapped_prefix, _ = os.path.splitext(unmapped_fn)
            unmapped_sorted_fn = '{}_sorted.bam'.format(unmapped_prefix)

            merged_fn = self.file_names['{}_sorted_by_name'.format(edge)]

            sam.sort_bam(accepted_fn, accepted_sorted_fn, by_name=True)
            sam.sort_bam(unmapped_fn, unmapped_sorted_fn, by_name=True)
            pysam.merge('-nf', merged_fn, accepted_sorted_fn, unmapped_sorted_fn)
    
    def combine_mappings(self):
        num_unmapped = 0
        num_R1_unmapped = 0
        num_R2_unmapped = 0
        num_nonunique = 0
        num_discordant = 0
        num_disoriented = 0
        num_concordant = 0

        tlens = Counter()

        R1_mappings = pysam.Samfile(self.file_names['R1_sorted_by_name'])
        R2_mappings = pysam.Samfile(self.file_names['R2_sorted_by_name'])

        R1_grouped = utilities.group_by(R1_mappings, lambda m: m.qname)
        R2_grouped = utilities.group_by(R2_mappings, lambda m: m.qname)
        with pysam.Samfile(self.file_names['combined'], 'wb', template=R1_mappings) as combined_bam_file:
            group_pairs = izip(R1_grouped, R2_grouped)
            for (R1_qname, R1_group), (R2_qname, R2_group) in group_pairs:
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

                R1_strand = '-' if R1_m.is_reverse else '+'
                R2_strand = '-' if R2_m.is_reverse else '+'

                tlen = max(R1_m.aend, R2_m.aend) - min(R1_m.pos, R2_m.pos)
                discordant = (R1_m.tid != R2_m.tid) or (R1_strand) != (R2_strand) or (tlen > 10000)
                if discordant:
                    num_discordant += 1
                    continue

                tlens[tlen] += 1
                
                if R1_strand == '+':
                    first_read = R2_m
                    second_read = R1_m
                elif R1_strand == '-':
                    first_read = R1_m
                    second_read = R2_m
                
                if R1_strand == '+':
                    extra_pairs = [(read, ref) for read, ref in R2_m.aligned_pairs if ref != None and read != None and ref < R1_m.pos]
                    if extra_pairs == []:
                        bases_from_R2 = 0
                        gap_cigar = []
                    else:
                        last_read, last_ref = extra_pairs[-1]
                        bases_from_R2 = last_read + 1
                        if last_ref == None:
                            print extra_pairs
                            print R1_m
                            print R2_m
                        gap = R1_m.pos - last_ref - 1
                        if gap < 0:
                            raise ValueError
                        gap_cigar = [(sam.BAM_CREF_SKIP, gap)]

                    combined_read = pysam.AlignedRead()
                    combined_read.qname = fastq.get_pair_name(R1_m.qname)
                    combined_read.tid = R1_m.tid
                    combined_read.pos = R2_m.pos
                    
                    if R1_m.pos < R2_m.pos:
                        num_disoriented += 1
                        continue

                    combined_read.is_reverse = False
                    combined_read.mapq = min(R1_m.mapq, R2_m.mapq)
                    combined_read.rnext = -1
                    combined_read.pnext = -1
                
                    combined_read.seq = R2_m.seq[:bases_from_R2] + R1_m.seq
                    combined_read.qual = R2_m.qual[:bases_from_R2] + R1_m.qual
                    truncated_R2_cigar = sam.truncate_cigar_blocks(R2_m.cigar, bases_from_R2)
                    combined_read.cigar = truncated_R2_cigar + gap_cigar + R1_m.cigar

                elif R1_strand == '-':
                    extra_pairs = [(read, ref) for read, ref in R2_m.aligned_pairs if read != None and ref >= R1_m.aend]
                    if extra_pairs == []:
                        bases_from_R2 = 0
                        gap_cigar = []
                    else:
                        first_read, first_ref = extra_pairs[0]
                        bases_from_R2 = R2_m.qlen - first_read
                        gap = first_ref - R1_m.aend
                        if gap < 0:
                            raise ValueError
                        gap_cigar = [(3, gap)]

                    combined_read = pysam.AlignedRead()
                    combined_read.qname = fastq.get_pair_name(R1_m.qname)
                    combined_read.tid = R1_m.tid
                    combined_read.pos = R1_m.pos
                    
                    if R2_m.pos < R1_m.pos:
                        num_disoriented += 1
                        continue

                    combined_read.is_reverse = True
                    combined_read.mapq = min(R1_m.mapq, R2_m.mapq)
                    combined_read.rnext = -1
                    combined_read.pnext = -1
                
                    combined_read.seq = R1_m.seq + R2_m.seq[-bases_from_R2:]
                    combined_read.qual = R1_m.qual + R2_m.qual[-bases_from_R2:]
                    truncated_R2_cigar = sam.truncate_cigar_blocks(R2_m.cigar[::-1], bases_from_R2)[::-1]
                    combined_read.cigar = R1_m.cigar + gap_cigar + truncated_R2_cigar

                trim.set_nongenomic_length(combined_read, 0)
                combined_bam_file.write(combined_read)
                num_concordant += 1

        self.log.extend(
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

        sam.sort_bam(self.file_names['combined'], self.file_names['combined_sorted'])
    
    def get_polyA_positions(self):
        piece_CDSs, max_gene_length = self.get_CDSs()
        gene_infos = positions.get_Transcript_polyA_position_counts(self.merged_file_names['combined_sorted'],
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
        log_pairs = self.read_file('log')
        log_dict = {name: values[0] for name, values in log_pairs}
        total_mapped_reads = log_dict['Nonunique'] + log_dict['Concordant']
        return total_mapped_reads

if __name__ == '__main__':
    script_path = os.path.realpath(__file__)
    map_reduce.controller(ThreeTFillExperiment, script_path)
