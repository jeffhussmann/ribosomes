import os.path
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import subprocess
import pysam
from itertools import izip
from collections import defaultdict
from Circles import mapping_tools
import gtf
import recycling
from Circles import sam
from Circles import fastq
from Circles import Serialize
from Circles.utilities import base_order, base_to_index, base_to_complement_index

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'BlueViolet', 'Gold']
    
# To some extent, this is linked to the definition of simple_CDS, which
# excludes genes that have another gene within 50 bp.
edge_buffer = 50

def map_tophat(reads_file_name,
               bowtie2_index,
               gtf_file_name,
               transcriptome_index,
               tophat_dir,
              ):
    tophat_command = ['tophat2',
                      '--GTF', gtf_file_name,
                      '--no-novel-juncs',
                      '--num-threads', '1',
                      #'--bt2-mm',
                      '--output-dir', tophat_dir,
                      '--transcriptome-index', transcriptome_index,
                      bowtie2_index,
                      reads_file_name,
                     ]
    # tophat maintains its own logs of everything that is written to the
    # console, so discard output.
    subprocess.check_output(tophat_command, stderr=subprocess.STDOUT)

A_site_offsets = {'ingolia_cell': {29: 15,
                                   30: 15,
                                   31: 16,
                                   32: 16,
                                   33: 16,
                                   34: 17,
                                   35: 17,
                                  },
                  'guo_nature':   {27: 15,
                                   28: 15,
                                   29: 15,
                                   30: 15,
                                   31: 15,
                                   32: 15,
                                  },
                  'yeast':        {28: 15,
                                   29: 15,
                                   30: 16,
                                  },
                 }

def get_codon_counts(gene_info, offset_type):
    # gene_info['CDS_length'] is the index of the first nucleotide of the
    # stop codon. Ingolia's original model never has the stop codon in the
    # A site, but subsequent data show an accumulation of (typically
    # length 29 or 30) reads that do advance this far. The +1 has been
    # added to include this.
    codon_counts = np.zeros(gene_info['CDS_length'] // 3 + 1, int)
    for length in set(gene_info['position_counts']) & set(A_site_offsets[offset_type]):
        position_counts = gene_info['position_counts'][length]
        start_index = -A_site_offsets[offset_type][length]
        end_index = gene_info['CDS_length'] - A_site_offsets[offset_type][length] + 3
        codon_counts += position_counts[start_index:end_index:3] + \
                        position_counts[start_index - 1:end_index - 1:3] + \
                        position_counts[start_index + 1:end_index + 1:3]

    return codon_counts

def get_total_read_count(gene_info):
    # TODO: this needs to be updated to new offset handling strategy
    A_site_offset = 15
    start_index = -A_site_offset
    end_index = gene_info['CDS_length'] - A_site_offset + 3
    count = gene_info['position_counts']['all'][start_index:end_index].sum()
    return count

def compute_read_counts(gene_infos, stringency):
    for gene_name in gene_infos:
        if stringency == 'everything':
            read_count = get_total_read_count(gene_infos[gene_name])
        else:
            read_count = get_codon_counts(gene_infos[gene_name], stringency).sum()
        expression = np.array([read_count, 0])
        gene_infos[gene_name]['expression'] = expression
    return gene_infos

def compute_RPKMs(gene_infos):
    RPKMs = {}
    total_mapped_reads = 0
    for gene_name in gene_infos:
        read_count = gene_infos[gene_name]['expression'][0]
        length = gene_infos[gene_name]['CDS_length']
        RPKMs[gene_name] = float(read_count) / length
        total_mapped_reads += read_count

    for gene_name in RPKMs:
        RPKMs[gene_name] = 1.e9 * RPKMs[gene_name] / total_mapped_reads

    return RPKMs

def get_TPMs(rpf_positions_dict):
    TPMs = {}
    T = 0
    for gene_name in rpf_positions_dict:
        counts = get_codon_counts(rpf_positions_dict[gene_name])
        length = rpf_positions_dict[gene_name]['CDS_length']
        reads = counts.sum()
        TPMs[gene_name] = float(reads) / length
        T += TPMs[gene_name]

    for gene_name in rpf_positions_dict:
        length = rpf_positions_dict[gene_name]['CDS_length']
        TPMs[gene_name] = 1.e6 * TPMs[gene_name] / T

    return TPMs

def determine_ambiguity_of_positions(extent, genome_index):
    edge_overlap = 50

    seqname, gene_strand, start, end = extent
    gene_length = abs(end - start) + 1  
    
    ref_seq = str(genome_index[seqname].seq)[start - 2 * edge_overlap:end + edge_overlap + 1]
    
    shape = (3, 3 * edge_overlap + gene_length)
    position_ambiguity = np.empty(shape, int)
    for length in [28, 29, 30]:
        for start in range(3 * edge_overlap + gene_length - length):
            # Explanation of encoding:
            # 0 means a fragment of this length starting at this position would
            # end in an A and therefore can't be the result of poly-A trimming.
            # 1 means a fragment of this length starting at this position would
            # end in a non-A followed by an A, so the read produced can't be
            # distinguished from a read of greater length starting at the same
            # position.
            # 2 means a fragment of this length starting at this position would
            # end in a non-A followed by another non-A and is therefore
            # unambiguous.
            if ref_seq[start + length - 1].upper() == 'A':
                position_ambiguity[length - 28, start] = 0
            elif ref_seq[start + length].upper() == 'A':
                position_ambiguity[length - 28, start] = 1
            else:
                position_ambiguity[length - 28, start] = 2
                
    return position_ambiguity

def read_codon_counts_file(fn, reports_P_site=False):
    genes = {}
    for line in open(fn):
        name, values = line.strip().split('\t', 1)
        values = np.fromstring(values, dtype=float, sep='\t')
        # If the P site is reported, shift over by one to convert to A site.
        if reports_P_site:
            values = np.concatenate(([0], values))
        genes[name] = values
            
    return genes

def scatter_positions(rpf_positions_list):
    edge_overlap = 50

    fig, ax = plt.subplots(figsize=(12, 16))

    counts_list = []

    for rpf_positions in rpf_positions_list:
        experiment_counts = []
        for (name, length), counts in zip(*rpf_positions): 
            start_at_codon = -8
            end_at_codon = length / 3
            codons = np.arange(start_at_codon, end_at_codon)
            codon_starts = np.arange(2 * edge_overlap + 3 * start_at_codon,
                                     2 * edge_overlap + 3 * end_at_codon,
                                     3,
                                    )
            counts = [counts[0][c] for c in codon_starts]
            experiment_counts.extend(counts)
        counts_list.append(experiment_counts)
      
    ax.scatter(counts_list[0], counts_list[1], s=1)

def get_ratios(first, second):
    assert set(first) == set(second)
    ratios = {key: np.divide(float(first[key]), second[key]) for key in first}
    return ratios

def plot_frameshifts(rpf_counts_list,
                     position_ambiguity_list,
                     edge_overlap,
                     gene_name,
                     gene_length,
                     exp_name,
                    ):
    ambiguity_to_color = {0: 'red',
                          1: 'green',
                          2: 'black',
                         }

    length_data = zip([28, 29, 30], rpf_counts_list, position_ambiguity_list)
    for fragment_length, rpf_counts, position_ambiguity in length_data:
        start_at_codon = 0 
        # codon_starts[i] is the index into a position array at which codon i
        # starts.
        codon_starts = np.arange(2 * edge_overlap + (start_at_codon * 3),
                                 2 * edge_overlap + gene_length,
                                 3,
                                )
        codon_numbers = start_at_codon + np.arange(len(codon_starts))
        
        # frame_counts_list[i, j] will be the number of RPF's starting at frame i of
        # codon j
        frame_counts_list = np.zeros((3, len(codon_starts)), int)
        # frame_colors_list will be used to visualize the ambiguity of each position
        frame_colors_list = [['']*len(codon_starts) for frame in range(3)]

        for c, codon_start in enumerate(codon_starts):
            for frame in range(3):
                frame_counts_list[frame, c] = rpf_counts[codon_start + frame]
                ambiguity = position_ambiguity[codon_start + frame]
                frame_colors_list[frame][c] = ambiguity_to_color[ambiguity]

        frames_so_far = frame_counts_list.cumsum(axis=1)
        fraction_frames_so_far = np.true_divide(frames_so_far, frames_so_far.sum(axis=0))

        frames_remaining = np.fliplr(np.fliplr(frame_counts_list).cumsum(axis=1))
        fraction_frames_remaining = np.true_divide(frames_remaining, frames_remaining.sum(axis=0))
        
        fig, axs = plt.subplots(4, 1, sharex=True)
        cumulative_ax = axs[0]
        frame_axs = axs[1:]

        for frame, (ax, frame_counts, frame_colors) in enumerate(zip(frame_axs, frame_counts_list, frame_colors_list)):
            ax.scatter(codon_numbers, frame_counts, s=20, c=frame_colors, linewidths=0)
            ax.set_ylim(-1, frame_counts_list.max() + 1)
            ax.set_xlim(codon_numbers[0], codon_numbers[-1])
            ax.set_title('Frame {0}'.format(frame))

        for frame, so_far, remaining, color in zip([0, 1, 2], fraction_frames_so_far, fraction_frames_remaining, colors):
            cumulative_ax.plot(codon_numbers, so_far, color=color, label='{0} so far'.format(frame))
            cumulative_ax.plot(codon_numbers, remaining, color=color, linestyle='--', label='{0} remaining'.format(frame))
            #difference = so_far - remaining
            #cumulative_ax.plot(codon_numbers, difference, color=color, linestyle=':')
            #difference = remaining - so_far
            #cumulative_ax.plot(codon_numbers, difference, color=color, linestyle=':')
        cumulative_ax.set_xlim(codon_numbers[0])
        cumulative_ax.set_ylim(-0.02, 1.02)
        
        cumulative_ax.legend(loc='upper right', framealpha=0.5)
        fig.suptitle('{0} - length {1} fragments\n{2}'.format(gene_name, fragment_length, exp_name))

    return codon_numbers, frame_counts_list, fraction_frames_so_far, fraction_frames_remaining

def plot_RPKMs():
    experiments = [
        ('geranshenko1', '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep1_foot/results/Initial_rep1_foot_rpf_positions.txt'),
        #('geranshenko2', '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep2_foot/results/Initial_rep2_foot_rpf_positions.txt'),
        ('ingolia1', '/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-1/results/Footprints-rich-1_rpf_positions.txt'),
        ('ingolia2', '/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-2/results/Footprints-rich-2_rpf_positions.txt'),
        ('ingolia1_mRNA', '/home/jah/projects/arlen/experiments/ingolia_science/mRNA-rich-1/results/mRNA-rich-1_rpf_positions.txt'),
        #('R98S', '/home/jah/projects/arlen/experiments/belgium_8_6_13/R98S_cDNA_sample/results/R98S_cDNA_sample_rpf_positions.txt'),
        #('suppressed', '/home/jah/projects/arlen/experiments/belgium_8_6_13/Suppressed_R98S_cDNA_sample/results/Suppressed_R98S_cDNA_sample_rpf_positions.txt'),
        #('WT',  '/home/jah/projects/arlen/experiments/belgium_8_6_13/WT_cDNA_sample/results/WT_cDNA_sample_rpf_positions.txt'),
    ]

    names = [name for name, _ in experiments]
    rpf_positions_lists = [Serialize.read_file(fn, 'rpf_positions') for _, fn in experiments]
    RPKMs_list = [get_TPMs(rpf_positions_list) for rpf_positions_list in rpf_positions_lists]
    vals = [[val for name, val in sorted(RPKMs.items())] for RPKMs in RPKMs_list]
    
    fig = plt.figure(figsize=(12, 12))
    for r in range(len(names)):
        for c in range(r + 1, len(names)):
            ax = fig.add_subplot(len(names) - 1, len(names) - 1, r * (len(names) - 1) + (c - 1) + 1)
            first = names[c]
            second = names[r]
            xs = vals[r]
            ys = vals[c]

            print sum(xs)
            print sum(ys)
            print
            
            ax.scatter(xs, ys, s=1)
            ax.set_yscale('log')
            ax.set_xscale('log')
            ax.set_xlim(1e-1, 5e4)
            ax.set_ylim(1e-1, 5e4)
            ax.plot([1e-2, 1e5], [1e-2, 1e5], '-', color='red', alpha=0.2)
            ax.set_xlabel(first)
            ax.set_ylabel(second)

def density_length_correlation():
    experiments = [
        ('ingolia1',
         '/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-1/results/Footprints-rich-1_rpf_positions.txt',
         '/home/jah/projects/arlen/experiments/ingolia_science/mRNA-rich-1/results/mRNA-rich-1_rpf_positions.txt',
        ),
        #('geranshenko1',
        # '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep1_foot/results/Initial_rep1_foot_rpf_positions.txt',
        # '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep1_mRNA/results/Initial_rep1_mRNA_rpf_positions.txt',
        #),
        ('R98S',
         '/home/jah/projects/arlen/experiments/belgium_8_6_13/R98S_cDNA_sample/results/R98S_cDNA_sample_rpf_positions.txt',
         '/home/jah/projects/arlen/experiments/belgium_8_6_13/R98S_cDNA_mRNA/results/R98S_cDNA_mRNA_rpf_positions.txt',
        ),
    ]

    for name, rpf_fn, mRNA_fn in experiments:
        rpfs = Serialize.read_file(rpf_fn, 'rpf_positions')
        rpf_TPMs = get_TPMs(rpfs)
        mRNAs = Serialize.read_file(mRNA_fn, 'rpf_positions')
        mRNA_TPMs = get_TPMs(mRNAs)
        ratios = get_ratios(rpf_TPMs, mRNA_TPMs)
        ratios_vals = [np.log2(val) for _, val in sorted(ratios.items())]
        lengths = [rpfs[gene_name]['CDS_length'] for gene_name in rpfs]
        fig, ax = plt.subplots()
        ax.scatter(lengths, ratios_vals, s=1)
        print name
        print scipy.stats.pearsonr(lengths, ratios_vals)
        print scipy.stats.spearmanr(lengths, ratios_vals)
        ax.set_xlabel('CDS length')
        ax.set_ylabel('Translational efficiency (log_2 ratio of TPMs)')
        ax.set_title(name)
    
def plot_all_starts_ends():
    experiments = [
        ('geranshenko1', 
         '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep1_foot/results/Initial_rep1_foot_from_starts.txt',
         '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep1_foot/results/Initial_rep1_foot_from_ends.txt',
        ),
        ('ingolia1',
         '/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-1/results/Footprints-rich-1_from_starts.txt',
         '/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-1/results/Footprints-rich-1_from_ends.txt',
        ),
        ('ingolia2',
         '/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-2/results/Footprints-rich-2_from_starts.txt',
         '/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-2/results/Footprints-rich-2_from_ends.txt',
        ),
        ('R98S',
         '/home/jah/projects/arlen/experiments/belgium_8_6_13/R98S_cDNA_sample/results/R98S_cDNA_sample_from_starts.txt',
         '/home/jah/projects/arlen/experiments/belgium_8_6_13/R98S_cDNA_sample/results/R98S_cDNA_sample_from_ends.txt',
        ),
        ('suppressed',
         '/home/jah/projects/arlen/experiments/belgium_8_6_13/Suppressed_R98S_cDNA_sample/results/Suppressed_R98S_cDNA_sample_from_starts.txt',
         '/home/jah/projects/arlen/experiments/belgium_8_6_13/Suppressed_R98S_cDNA_sample/results/Suppressed_R98S_cDNA_sample_from_ends.txt',
        ),
        ('WT',
         '/home/jah/projects/arlen/experiments/belgium_8_6_13/WT_cDNA_sample/results/WT_cDNA_sample_from_starts.txt',
         '/home/jah/projects/arlen/experiments/belgium_8_6_13/WT_cDNA_sample/results/WT_cDNA_sample_from_ends.txt',
        ),
    ]

    names = [name for name, _, _ in experiments]
    from_starts_list = [Serialize.read_file(fn, 'array') for _, fn, _ in experiments]
    from_ends_list = [Serialize.read_file(fn, 'array') for _, _, fn in experiments]
    fig_file_name = '/home/jah/projects/arlen/results/compare_starts_and_ends.pdf'
    plot_starts_and_ends_new(from_starts_list, from_ends_list, names, fig_file_name)

def error_profile(bam_file_name, simple_CDSs, relevant_lengths):
    type_shape = (50,
                  fastq.MAX_EXPECTED_QUAL + 1,
                  len(base_order),
                  len(base_order),
                 )
    type_counts = {length: np.zeros(type_shape, int) for length in relevant_lengths}

    bamfile = pysam.Samfile(bam_file_name, 'rb')
    for CDS in simple_CDSs:
        reads = bamfile.fetch(CDS.seqname,
                              CDS.start - edge_buffer,
                              CDS.end + edge_buffer,
                             )
        for read in reads:
            if read.mapq != 50:
                # Non-unique mapping
                continue
            elif read.qlen not in relevant_lengths:
                continue
            elif sam.contains_indel_pysam(read):
                continue
            else:
                strand = '-' if read.is_reverse else '+'
                
                if strand != CDS.strand:
                    continue
                else:
                    alignment = sam.produce_alignment(read, from_pysam=True)

                    if strand == '+':
                        index_lookup = base_to_index
                    else:
                        index_lookup = base_to_complement_index

                    for ref_char, read_char, qual, ref_pos, read_pos in alignment:
                        ref_index = index_lookup[ref_char]
                        read_index = index_lookup[read_char]
                        coords = (read_pos, qual, ref_index, read_index)
                        type_counts[read.qlen][coords] += 1

    return type_counts

def recycling_ratios(rpf_positions_dict, simple_CDSs, genome):
    edge_overlap = 50
        
    ratio_lists = {amino_acid: defaultdict(list) for amino_acid in recycling.back_table}

    for gene in simple_CDSs:
        gene_name = gtf.parse_attribute(gene.attribute)['protein_id']
        codon_counts = get_codon_counts(rpf_positions_dict[gene_name], stringent=False)
        aa_locations = recycling.get_amino_acid_locations(gene, genome)
        if aa_locations == None:
            # MT or internal stop codon
            continue
        for amino_acid in aa_locations:
            location_pairs = izip(aa_locations[amino_acid], aa_locations[amino_acid][1:])
            for (first_position, first_codon), (second_position, second_codon) in location_pairs:
                first_count = codon_counts[first_position]
                second_count = codon_counts[second_position]
                ratio = (first_count, second_count)
                ratio_lists[amino_acid][first_codon, second_codon].append(ratio)
        
    return ratio_lists

def make_codon_counts_file(gene_infos, codon_counts_fn, offset_type):
    with open(codon_counts_fn, 'w') as codon_counts_fh:
        for gene_name in sorted(gene_infos):
            # Temporary hack
            if gene_infos[gene_name]['CDS_length'] % 3:
                print gene_name
                continue
            counts = get_codon_counts(gene_infos[gene_name], offset_type)
            counts_string = '\t'.join(str(count) for count in counts)
            line = '{0}\t{1}\n'.format(gene_name, counts_string)
            codon_counts_fh.write(line)

def make_RPKMs_file(gene_infos, RPKMs_fn):
    RPKMs = compute_RPKMs(gene_infos)
    with open(RPKMs_fn, 'w') as RPKMs_fh:
        for gene_name in sorted(RPKMs):
            line = '{0}\t{1:0.2f}\n'.format(gene_name, RPKMs[gene_name])
            RPKMs_fh.write(line)

def read_RPKMs_file(RPKMs_fn):
    def line_to_gene(line):
        name, value = line.strip().split()
        value = float(value)
        return name, value

    return dict(line_to_gene(line) for line in open(RPKMs_fn))

def read_bartel_file(bartel_file):
    RPF_dict = {}
    mRNA_dict = {}
    fh = open(bartel_file)
    # Skip column labels line
    fh.readline()
    for line in fh:
        name, mRNA_value, RPF_value = line.strip().split()
        mRNA_value = float(mRNA_value)
        RPF_value = float(RPF_value)
        mRNA_dict[name] = mRNA_value
        RPF_dict[name] = RPF_value

    return RPF_dict, mRNA_dict

def plot_metagene_averaged():
    # Generators that yields arrays of counts
    def counts_from_rpf_positions_fn(rpf_positions_fn):
        rpf_positions_dict = Serialize.read_file(rpf_positions_fn, 'rpf_positions')
        for gene_name in rpf_positions_dict:
            counts = get_codon_counts(rpf_positions_dict[gene_name],
                                      stringent=False,
                                     )
            yield counts

    def counts_from_premal_fn(premal_fn):
        counts_dict = read_premal_file(premal_fn)
        for gene_name in counts_dict:
            counts = counts_dict[gene_name]
            yield counts
        
    rpf_experiments = [
        #('Geranshenko1', '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep1_foot/results/Initial_rep1_foot_rpf_positions.txt'),
        #('Geranshenko2', '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep2_foot/results/Initial_rep2_foot_rpf_positions.txt'),
        #('Geranshenko_mrna', '/home/jah/projects/arlen/experiments/gerashchenko_pnas/5min_rep1_mRNA/results/5min_rep1_mRNA_rpf_positions.txt'),

        ('Ingolia1', '/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-1/results/Footprints-rich-1_rpf_positions.txt'),
        #('Ingolia2', '/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-2/results/Footprints-rich-2_rpf_positions.txt'),
        #('Ingolia1_mRNA', '/home/jah/projects/arlen/experiments/ingolia_science/mRNA-rich-1/results/mRNA-rich-1_rpf_positions.txt'),

        ('Brar1' ,'/home/jah/projects/arlen/experiments/brar_science/s_tA-fp_100211_l3_sequence/results/s_tA-fp_100211_l3_sequence_rpf_positions.txt'),
        ('Brar2' ,'/home/jah/projects/arlen/experiments/brar_science/s_gb15exp_veg_-fp_100219_l4_sequence/results/s_gb15exp_veg_-fp_100219_l4_sequence_rpf_positions.txt'),
        ('Brar3' ,'/home/jah/projects/arlen/experiments/brar_science/s_14201exp_veg_-fp_100219_l6_sequence/results/s_14201exp_veg_-fp_100219_l6_sequence_rpf_positions.txt'),
        ('Brar4' ,'/home/jah/projects/arlen/experiments/brar_science/s_t1-fp_090807_l4/results/s_t1-fp_090807_l4_rpf_positions.txt'),

        #('suppressed', '/home/jah/projects/arlen/experiments/belgium_8_6_13/Suppressed_R98S_cDNA_sample/results/Suppressed_R98S_cDNA_sample_rpf_positions.txt'),
        #('R98S', '/home/jah/projects/arlen/experiments/belgium_8_6_13/R98S_cDNA_sample/results/R98S_cDNA_sample_rpf_positions.txt'),
        #('WT',  '/home/jah/projects/arlen/experiments/belgium_8_6_13/WT_cDNA_sample/results/WT_cDNA_sample_rpf_positions.txt'),
    ]

    premal_experiments = [
        ('Bartel', '/home/jah/projects/arlen/experiments/plotkin/genePosReads.txt'),
    ]

    all_experiments = [(name, counts_from_rpf_positions_fn(fn)) for name, fn in rpf_experiments] + \
                      [(name, counts_from_premal_fn(fn)) for name, fn in premal_experiments]

    fig, ax = plt.subplots(figsize=(12, 12))

    plot_up_to = 500

    for name, counts_generator in all_experiments:
        sum_of_normalized = np.zeros(10000)
        long_enough_genes = np.zeros(10000)

        for counts in counts_generator:
            if counts.sum() < 64:
                continue

            num_codons = len(counts)
            density = counts.sum() / float(num_codons)
            normalized = counts / float(density)
            sum_of_normalized[:num_codons] += normalized
            long_enough_genes[:num_codons] += np.ones(num_codons)
        
        mean_densities = sum_of_normalized / long_enough_genes 

        ax.plot(mean_densities[:plot_up_to], '.-', label=name)
    
    ax.legend(loc='upper right', framealpha=0.5)

    ax.set_xlabel('Position (codons)')
    ax.set_ylabel('Normalized mean reads')
    
    ax.plot(np.ones(plot_up_to), color='black', alpha=0.5)
    ax.set_ylim(0, 3)

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ax.set_aspect((xmax - xmin) / (ymax - ymin))

def plot_metagene_from_codon_counts():
    # Generators that yields arrays of counts
    def counts_from_codon_counts_fn(codon_counts_fn, reports_P_site):
        counts_dict = read_codon_counts_file(codon_counts_fn, reports_P_site)
        for gene_name in counts_dict:
            counts = counts_dict[gene_name]
            yield counts
        
    experiments = [
        #('Bartel', '/home/jah/projects/arlen/experiments/plotkin/genePosReads.txt', True),
        ('Weinberg_stringent', '/home/jah/projects/arlen/experiments/weinberg/RPF/results/RPF_stringent_codon_counts.txt', False),
        ('Weinberg', '/home/jah/projects/arlen/experiments/weinberg/RPF/results/RPF_codon_counts.txt', False),
        #('Ingolia', '/home/jah/projects/arlen/results/Ingolia_RPF_codon_counts.txt', False),
        #('Gerashchenko', '/home/jah/projects/arlen/results/Gerashchenko_RPF_codon_counts.txt', False),
        #('Brar', '/home/jah/projects/arlen/results/Brar_RPF_codon_counts.txt', False),
        #('Hussmann_WT', '/home/jah/projects/arlen/results/UT_WT_RPF_codon_counts.txt', False),
        #('Zinshteyn_1', '/home/jah/projects/arlen/results/Zinshteyn_1_RPF_codon_counts.txt', False),
        #('McManus', '/home/jah/projects/arlen/experiments/mcmanus_gr/S._cerevisiae_Ribo-seq_Rep_1/results/S._cerevisiae_Ribo-seq_Rep_1_codon_counts.txt', False),
        #('Zinshteyn_1_mRNA', '/home/jah/projects/arlen/results/Zinshteyn_1_mRNA_codon_counts.txt', False),
        #('Zinshteyn_2', '/home/jah/projects/arlen/results/Zinshteyn_2_RPF_codon_counts.txt', False),
        ('no_aditive', '/home/jah/projects/arlen/experiments/guydosh/wild-type_no_additive/results/wild-type_no_additive_codon_counts.txt', False),
        ('CHX', '/home/jah/projects/arlen/experiments/guydosh/wild-type_CHX/results/wild-type_CHX_codon_counts.txt', False),
        #('no_aditive', '/home/jah/projects/arlen/experiments/guydosh/wild-type_no_additive/results/wild-type_no_additive_codon_counts.txt', False),
        #('no_aditive', '/home/jah/projects/arlen/experiments/guydosh/wild-type_no_additive/results/wild-type_no_additive_codon_counts.txt', False),
    ]

    all_experiments = [(name, counts_from_codon_counts_fn(fn, reports_P_site))
                       for name, fn, reports_P_site in experiments]

    fig, ax = plt.subplots(figsize=(12, 12))

    plot_up_to = 500

    for name, counts_generator in all_experiments:
        sum_of_normalized = np.zeros(10000)
        long_enough_genes = np.zeros(10000)

        for counts in counts_generator:
            if counts.sum() < 64:
                continue

            num_codons = len(counts)
            density = counts.sum() / float(num_codons)
            normalized = counts / float(density)
            sum_of_normalized[:num_codons] += normalized
            long_enough_genes[:num_codons] += np.ones(num_codons)
        
        mean_densities = sum_of_normalized / long_enough_genes 

        ax.plot(mean_densities[:plot_up_to], '.-', label=name)
    
    ax.legend(loc='upper right', framealpha=0.5)

    ax.set_xlabel('Position (codons)')
    ax.set_ylabel('Normalized mean reads')
    
    ax.plot(np.ones(plot_up_to), color='black', alpha=0.5)
    ax.set_ylim(0, 3)

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ax.set_aspect((xmax - xmin) / (ymax - ymin))

def plot_metagene_unaveraged(from_end=False):
    # Generators that yields arrays of counts
    def counts_from_rpf_postiions_fn(rpf_positions_fn):
        rpf_positions_dict = Serialize.read_file(rpf_positions_fn, 'rpf_positions')
        for gene_name in rpf_positions_dict:
            counts = get_codon_counts(rpf_positions_dict[gene_name],
                                      stringent=False,
                                     )
            yield counts

    def counts_from_premal_fn(premal_fn):
        counts_dict = read_premal_file(premal_fn)
        for gene_name in counts_dict:
            counts = counts_dict[gene_name]
            yield counts
        
    rpf_experiments = [
        #('Geranshenko1', '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep1_foot/results/Initial_rep1_foot_rpf_positions.txt'),
        #('Geranshenko2', '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep2_foot/results/Initial_rep2_foot_rpf_positions.txt'),
        #('Geranshenko_mrna', '/home/jah/projects/arlen/experiments/gerashchenko_pnas/5min_rep1_mRNA/results/5min_rep1_mRNA_rpf_positions.txt'),

        #('Ingolia1', '/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-1/results/Footprints-rich-1_rpf_positions.txt'),
        #('Ingolia2', '/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-2/results/Footprints-rich-2_rpf_positions.txt'),
        #('Ingolia1_mRNA', '/home/jah/projects/arlen/experiments/ingolia_science/mRNA-rich-1/results/mRNA-rich-1_rpf_positions.txt'),

        #('Brar1' ,'/home/jah/projects/arlen/experiments/brar_science/s_tA-fp_100211_l3_sequence/results/s_tA-fp_100211_l3_sequence_rpf_positions.txt'),
        #('Brar2' ,'/home/jah/projects/arlen/experiments/brar_science/s_gb15exp_veg_-fp_100219_l4_sequence/results/s_gb15exp_veg_-fp_100219_l4_sequence_rpf_positions.txt'),
        #('Brar3' ,'/home/jah/projects/arlen/experiments/brar_science/s_14201exp_veg_-fp_100219_l6_sequence/results/s_14201exp_veg_-fp_100219_l6_sequence_rpf_positions.txt'),
        #('Brar4' ,'/home/jah/projects/arlen/experiments/brar_science/s_t1-fp_090807_l4/results/s_t1-fp_090807_l4_rpf_positions.txt'),
        
        ('Nagalakshmi_RH_ori', '/home/jah/projects/arlen/experiments/nagalakshmi_science/RH_ori/results/RH_ori_rpf_positions.txt'),
        ('Nagalakshmi_RH_bio', '/home/jah/projects/arlen/experiments/nagalakshmi_science/RH_ori/results/RH_ori_rpf_positions.txt'),
        ('Nagalakshmi_dT_ori', '/home/jah/projects/arlen/experiments/nagalakshmi_science/dT_ori/results/dT_ori_rpf_positions.txt'),
        ('Nagalakshmi_dT_bio', '/home/jah/projects/arlen/experiments/nagalakshmi_science/dT_bio/results/dT_bio_rpf_positions.txt'),

        #('suppressed', '/home/jah/projects/arlen/experiments/belgium_8_6_13/Suppressed_R98S_cDNA_sample/results/Suppressed_R98S_cDNA_sample_rpf_positions.txt'),
        #('R98S', '/home/jah/projects/arlen/experiments/belgium_8_6_13/R98S_cDNA_sample/results/R98S_cDNA_sample_rpf_positions.txt'),
        #('WT',  '/home/jah/projects/arlen/experiments/belgium_8_6_13/WT_cDNA_sample/results/WT_cDNA_sample_rpf_positions.txt'),
    ]

    premal_experiments = [
        #('Bartel', '/home/jah/projects/arlen/experiments/plotkin/genePosReads.txt'),
    ]

    all_experiments = [(name, counts_from_rpf_postiions_fn(fn)) for name, fn in rpf_experiments] + \
                      [(name, counts_from_premal_fn(fn)) for name, fn in premal_experiments]

    plot_to = 500
    fig_cumulative, ax_cumulative = plt.subplots()

    if from_end:
        xs = np.arange(0, -plot_to, -1)
    else:
        xs = np.arange(plot_to)

    for name, counts_generator in all_experiments:
        print name

        expected_counts = np.zeros(5000)
        actual_counts = np.zeros(5000)

        for counts in counts_generator:
            num_codons = len(counts)
            r_g = counts.sum()
            uniform_counts = np.ones(num_codons) * r_g / num_codons
            
            actual_counts[:num_codons] += counts
            expected_counts[:num_codons] += uniform_counts

        normalized_cumulative = np.cumsum(actual_counts - expected_counts) / actual_counts.sum()

        #ax_cumulative.plot(xs, normalized_cumulative[:plot_to], label=name)
        ax_cumulative.plot(xs, ((actual_counts - expected_counts) / expected_counts)[:plot_to], label=name)

    ax_cumulative.plot(xs, np.zeros(plot_to), 'k--')
    ax_cumulative.legend()
    if from_end:
        xlabel = 'Codon position relative to end'
    else:
        xlabel = 'Codon position relative to start'
    ax_cumulative.set_xlabel(xlabel)
    ax_cumulative.set_ylabel('Cumulative sum of (actual - expected)')

if __name__ == '__main__':
    genome_index = mapping_tools.get_genome_index('/home/jah/projects/arlen/data/organisms/saccharomyces_cerevisiae/genome', explicit_path=True)

    #clean_bam_fn = '/home/jah/remote/slate/projects/arlen/experiments/belgium_8_6_13/WT_cDNA_sample/results/WT_cDNA_sample_clean.bam'
    clean_bam_fn = '/home/jah/remote/slate/projects/arlen/experiments/belgium_8_6_13/R98S_cDNA_sample/results/R98S_cDNA_sample_clean.bam'
    #clean_bam_fn = '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep1_foot/results/Initial_rep1_foot_clean.bam'
    #clean_bam_fn = '/home/jah/remote/slate/projects/arlen/experiments/gerashchenko_pnas/Initial_rep2_foot/results/Initial_rep2_foot_clean.bam'
    _, exp_name = os.path.split(clean_bam_fn)
    
    gtf_fn = '/home/jah/projects/arlen/data/organisms/saccharomyces_cerevisiae/transcriptome/genes.gtf'
    
    #gene_name = 'YOR239W'
    #gene_name = 'YPL052W'
    #gene_name = 'YBL030C'
    #gene_name = 'YLR318W'   # EST2
    #gene_name = 'YIL009C-A'
    #gene_name = 'YDL220C' # CDC13
    #gene_name = 'YLR233C' # EST1
    #gene_name = 'YDR082W' # STN1
    #gene_name = 'YER103W' # higher mRNA in Suppressed than R98S
    #gene_name = 'YNL067W' # lower mRNA in R98S than WT
    gene_name = 'YHR096C' # higher mRNA in R98S than WT
    
    extent = gtf.get_extent_by_name(gtf_fn, gene_name)
    gene_length = extent[3] - extent[2] + 1
    position_counts, expression_counts = get_extent_positions(clean_bam_fn, extent, genome_index)
    position_ambiguity = determine_ambiguity_of_positions(extent, genome_index)
    codon_numbers, frames, so_far, remaining = plot_frameshifts(position_counts,
                                                                position_ambiguity,
                                                                50,
                                                                gene_name,
                                                                gene_length,
                                                                exp_name,
                                                               )

#if __name__ == '__main__':
#    density_length_correlation()
