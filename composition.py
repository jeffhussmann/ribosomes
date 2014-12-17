from composition_cython import *
import matplotlib.pyplot as plt
import pysam
from Sequencing import utilities
from Sequencing.Visualize import igv_colors

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

def length_stratified_plot(base_counts, fig_file_name, lengths=[28, 29, 30, 'all']):
    ''' Plot fractions of all base calls that are each base at each cycle.
    '''
    _, cycles, _ = base_counts.shape

    bases = fastq.base_order[:5]
    fig, axs = plt.subplots(len(lengths), 1, figsize=(15, len(lengths) * 6))

    for length, ax in zip(lengths, axs):
        if length == 'all':
            counts = base_counts.sum(axis=0)
        else:
            counts = base_counts[length]

        denominators = np.maximum(1, counts.sum(axis=1)).astype(float)
        total = int(denominators[0])

        for i, base in enumerate(bases):
            fractions = counts.T[i] / denominators
            ax.plot(fractions, '.-', color=igv_colors.normalized_rgbs[base], linewidth=0.5, label=base) 

        ax.set_ylim(ymin=0)
        ax.set_xlim(0, len(fractions) - 1)
        ax.legend(loc='upper right', framealpha=0.5)
        ax.set_title('Length {1}: {0:,d} total reads'.format(total, length))
    
    #fig.suptitle('Base composition vs. position in trimmed read')
    fig.savefig(fig_file_name, bbox_inches='tight')
    plt.close()
