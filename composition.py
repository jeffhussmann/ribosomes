from composition_cython import *
import matplotlib.pyplot as plt
import pysam
from Sequencing import utilities
from Sequencing.Visualize import igv_colors

igv_colors = igv_colors.normalized_rgbs

def get_seq_info_pairs(clean_bam_fn):
    clean_bam_file = pysam.Samfile(clean_bam_fn)
    
    for aligned_read in clean_bam_file:
        if aligned_read.is_unmapped or aligned_read.is_secondary:
            continue
        perfect_and_unique = dict(aligned_read.tags)['NM'] == 0 and aligned_read.mapq == 50
        if aligned_read.is_reverse:
            seq = utilities.reverse_complement(aligned_read.seq)
        else:
            seq = aligned_read.seq

        yield seq, perfect_and_unique

def length_stratified_plot(base_counts, fig_file_name, lengths):
    ''' Plot fractions of all base calls that are each base at each cycle.
    '''
    max_length, cycles, _ = base_counts.shape

    bases = fastq.base_order[:5]
    fig, axs = plt.subplots(len(lengths), 1, figsize=(15, len(lengths) * 6))

    from_right_base_counts = np.zeros_like(base_counts)
    for length in range(max_length):
        from_right_base_counts[length, cycles - length:] += base_counts[length, :length]

    for length, ax in zip(lengths, axs):
        if length == 'all':
            xs = np.arange(cycles)
            counts = base_counts.sum(axis=0)
        elif length == 'all_from_right':
            xs = np.arange(-cycles + 1, 1)
            counts = from_right_base_counts.sum(axis=0)
        else:
            xs = np.arange(length)
            counts = base_counts[length]

        denominators = np.maximum(1, counts.sum(axis=1)).astype(float)
        total = int(max(denominators))

        for i, base in enumerate(bases):
            fractions = counts.T[i] / denominators
            ax.plot(xs, fractions[:len(xs)], '.-', color=igv_colors[base], linewidth=0.5, label=base) 

        ax.set_ylim(0, 0.6)
        ax.set_title('Length {1}: {0:,d} total reads'.format(total, length))

        if length == 'all_from_right':
            ax.set_xlim(-len(fractions) + 1, 0.5)
            ax.legend(loc='upper left', framealpha=0.5)
        else:
            ax.set_xlim(-0.5, len(fractions) - 1)
            ax.legend(loc='upper right', framealpha=0.5)
    
    #fig.suptitle('Base composition vs. position in trimmed read')
    fig.savefig(fig_file_name, bbox_inches='tight')
    plt.close()
