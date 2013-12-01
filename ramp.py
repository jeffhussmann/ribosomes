import Serialize
import ribosomes
import numpy as np
import matplotlib.pyplot as plt

samples = [
    ('Ingolia1' ,'/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-1/results/Footprints-rich-1_rpf_positions.txt'),
    #('R98S' ,'/home/jah/projects/arlen/experiments/belgium_8_6_13/R98S_cDNA_sample/results/R98S_cDNA_sample_rpf_positions.txt'),
    #('Suppressed' ,'/home/jah/projects/arlen/experiments/belgium_8_6_13/Suppressed_R98S_cDNA_sample/results/Suppressed_R98S_cDNA_sample_rpf_positions.txt'),
    ('WT' ,'/home/jah/projects/arlen/experiments/belgium_8_6_13/WT_cDNA_sample/results/WT_cDNA_sample_rpf_positions.txt'),
    ('Gerashchenko1' ,'/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep1_foot/results/Initial_rep1_foot_rpf_positions.txt'),
    #('Gerashchenko2' ,'/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep2_foot/results/Initial_rep2_foot_rpf_positions.txt'),
    ('Brar' ,'/home/jah/projects/arlen/experiments/brar_science/s_tA-fp_100211_l3_sequence/results/s_tA-fp_100211_l3_sequence_rpf_positions.txt'),
    #('Ingolia2' ,'/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-2/results/Footprints-rich-2_rpf_positions.txt'),
    #('R98S_mRNA' ,'/home/jah/projects/arlen/experiments/belgium_8_6_13/R98S_cDNA_mRNA/results/R98S_cDNA_mRNA_rpf_positions.txt'),
    #('Suppressed_mRNA' ,'/home/jah/projects/arlen/experiments/belgium_8_6_13/Suppressed_R98S_cDNA_mRNA/results/Suppressed_R98S_cDNA_mRNA_rpf_positions.txt'),
    #('WT_mRNA' ,'/home/jah/projects/arlen/experiments/belgium_8_6_13/WT_cDNA_mRNA/results/WT_cDNA_mRNA_rpf_positions.txt'),
    #('Gerashchenko_mRNA1' ,'/home/jah/projects/arlen/experiments/gerashchenko_pnas/5min_rep1_mRNA/results/5min_rep1_mRNA_rpf_positions.txt'),
    #('Ingolia_mRNA1' ,'/home/jah/projects/arlen/experiments/ingolia_science/mRNA-rich-1/results/mRNA-rich-1_rpf_positions.txt'),
    #('Ingolia_mRNA2' ,'/home/jah/projects/arlen/experiments/ingolia_science/mRNA-rich-2/results/mRNA-rich-2_rpf_positions.txt'),
]

fig_cumulative, ax_cumulative = plt.subplots()
fig_raw, ax_raw = plt.subplots()

exclude_initiation = True
from_end = True
min_num_codons = 0
max_num_codons = 10000

plot_to = 700

if from_end:
    xs = np.arange(0, -plot_to, -1)
else:
    xs = np.arange(plot_to)

for sample_name, rpf_fn in samples:
    print sample_name

    rpf_positions_list = Serialize.read_file(rpf_fn, 'rpf_positions')

    expected_counts = np.zeros(5000)
    actual_counts = np.zeros(5000)

    for gene_name in rpf_positions_list:
        counts = rpf_positions_list[gene_name]['counts'][28]
        A_site_offset = 15
        start_index = 100 - A_site_offset
        if exclude_initiation:
            start_index += 6
        end_index = -(50 + A_site_offset)
        in_frames = counts[start_index:end_index + 30:3]
        if from_end:
            in_frames = in_frames[::-1]

        num_codons = len(in_frames)
        r_g = in_frames.sum()
        # uniform, by definition, has the same expected counts if reversed
        uniform_counts = np.ones(num_codons) * r_g / num_codons
        
        if min_num_codons <= num_codons <= max_num_codons:
            actual_counts[:num_codons] += in_frames
            expected_counts[:num_codons] += uniform_counts

    normalized_cumulative = np.cumsum(actual_counts - expected_counts) / actual_counts.sum()
    
    normalized_expected = expected_counts# / actual_counts.sum()
    normalized_actual = actual_counts# / actual_counts.sum()

    ax_cumulative.plot(xs, normalized_cumulative[:plot_to], label=sample_name)
    ax_raw.plot(xs, normalized_expected[:plot_to],
                label='{0}, expected if uniform'.format(sample_name),
               )
    ax_raw.plot(xs, normalized_actual[:plot_to],
                '.', label='{0}, actually observed'.format(sample_name),
               )

premal_fn = '/home/jah/projects/arlen/experiments/plotkin/genePosReads.txt'
genes = ribosomes.read_premal_file(premal_fn)

expected_counts = np.zeros(5000)
actual_counts = np.zeros(5000)

for gene_name in genes:
    counts = genes[gene_name]
    if exclude_initiation:
        counts = counts[1:]
    if from_end:
        counts = counts[::-1]
    num_codons = len(counts)
    r_g = counts.sum()
    # uniform, by definition, has the same expected counts if reversed
    uniform_counts = np.ones(num_codons) * r_g / num_codons
    
    if min_num_codons <= num_codons <= max_num_codons:
        actual_counts[:num_codons] += counts
        expected_counts[:num_codons] += uniform_counts

normalized_cumulative = np.cumsum(actual_counts - expected_counts) / actual_counts.sum()

normalized_expected = expected_counts# / actual_counts.sum()
normalized_actual = actual_counts# / actual_counts.sum()

ax_cumulative.plot(xs, normalized_cumulative[:plot_to], label='bartel')
ax_raw.plot(xs, normalized_expected[:plot_to],
            label='{0}, expected if uniform'.format('bartel'),
           )
ax_raw.plot(xs, normalized_actual[:plot_to],
            '.', label='{0}, actually observed'.format('bartel'),
           )

if from_end:
    xlabel = 'Codon position relative to end'
else:
    xlabel = 'Codon position relative to start'

ax_cumulative.plot(xs, np.zeros(plot_to), 'k--')
ax_cumulative.legend()
ax_cumulative.set_xlabel(xlabel)
ax_cumulative.set_ylabel('Cumulative sum of (actual - expected)')

ax_raw.legend()
ax_raw.set_xlabel(xlabel)
ax_raw.set_ylabel('Counts')
#ax_raw.set_ylim(0, 3500)
