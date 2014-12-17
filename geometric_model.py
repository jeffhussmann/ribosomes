import numpy as np
import matplotlib.pyplot as plt
import Serialize
import positions
import brewer2mpl
from scipy.optimize import leastsq
from itertools import cycle

mRNA_experiments = [
    ##('test', '/home/jah/projects/arlen/experiments/weinberg/mRNA/results/test_read_positions.txt'),
    #('Weinberg (RiboZero)', '/home/jah/projects/arlen/experiments/weinberg/mRNA/results/mRNA_read_positions_all.txt'),
    #('Weinberg (dynabeads)', '/home/jah/projects/arlen/experiments/weinberg/Dynabeads/results/Dynabeads_read_positions_all.txt'),
    ('Ingolia', '/home/jah/projects/arlen/experiments/ingolia_science/mRNA-rich-1/results/mRNA-rich-1_read_positions_all.txt'),
    #('McManus', '/home/jah/projects/arlen/experiments/mcmanus_gr/S._cerevisiae_RNA-seq_Rep_1/results/S._cerevisiae_RNA-seq_Rep_1_read_positions_all.txt'),
    #('Guydosh', '/home/jah/projects/arlen/experiments/guydosh_cell/wild-type_mRNA-Seq/results/wild-type_mRNA-Seq_read_positions_all.txt'),
    #('Gerashchenko', '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep1_mRNA/results/Initial_rep1_mRNA_read_positions_all.txt'),
    #('Arlen', '/home/jah/projects/arlen/experiments/belgium_8_6_13/WT_cDNA_mRNA/results/WT_cDNA_mRNA_read_positions_all.txt'),
    #('Zinshteyn', '/home/jah/projects/arlen/experiments/zinshteyn_plos_genetics/WT_RNA-seq_1/results/WT_RNA-seq_1_read_positions_all.txt'),
    #('Weinberg (dynabeads)', '/home/jah/projects/arlen/experiments/weinberg/Dynabeads/results/Dynabeads_read_positions_all.txt'),
    
    ('Weinberg (RiboZero)', '/home/jah/projects/arlen/experiments/weinberg/mRNA/results/mRNA_read_positions_all.txt'),
    ('Weinberg (Unselected)', '/home/jah/projects/arlen/experiments/weinberg/Unselected/results/Unselected_read_positions_all.txt'),


    #('Zinshteyn_mRNA_1', '/home/jah/projects/arlen/experiments/zinshteyn_plos_genetics/WT_RNA-seq_1/results/WT_RNA-seq_1_read_positions_all.txt'),
    #('Weinberg (unselected)', '/home/jah/projects/arlen/experiments/weinberg/Unselected/results/Unselected_read_positions_all.txt'),
    
    #('Nagalakshmi_RH_ori', '/home/jah/projects/arlen/experiments/nagalakshmi_science/RH_ori/results/RH_ori_read_positions.txt'),
    #('Nagalakshmi_RH_bio', '/home/jah/projects/arlen/experiments/nagalakshmi_science/RH_bio/results/RH_bio_read_positions.txt'),
    #('Nagalakshmi_dT_ori', '/home/jah/projects/arlen/experiments/nagalakshmi_science/dT_ori/results/dT_ori_read_positions.txt'),
    #('Nagalakshmi_dT_bio', '/home/jah/projects/arlen/experiments/nagalakshmi_science/dT_bio/results/dT_bio_read_positions.txt'),
]

def smoothed(array, window_size):
    smoothed_array = np.zeros_like(array)
    for i in range(window_size):
        smoothed_array[i] = array[:i + 1].sum() / float(i + 1)
    for i in range(window_size, len(array) - window_size):
        smoothed_array[i] = array[i - window_size:i + window_size + 1].sum() / float(2 * window_size + 1)
    for i in range(len(array) - window_size, len(array)):
        smoothed_array[i] = array[i:].sum() / float(len(array) - i)
    return smoothed_array

def plot_mRNA_metagene_unaveraged(from_end, min_length, max_length):
    # Generators that yields arrays of counts
    def counts_from_read_positions_fn(read_positions_fn, from_end):
        gene_infos = Serialize.read_file(read_positions_fn, 'read_positions')
        for gene_name in gene_infos:
            #if gene_name == 'YLR256W':
            #    continue
            if gene_name in {'YLR249W', 'YPL106C', 'YGL008C'}:
                continue
            if from_end:
                counts = gene_infos[gene_name]['all'].relative_to_end
            else:
                counts = gene_infos[gene_name]['all']
            yield counts

    bmap = brewer2mpl.get_map('Set1', 'qualitative', 9)
    colors = cycle(bmap.mpl_colors[:5] + bmap.mpl_colors[6:])

    experiments = [(name, counts_from_read_positions_fn(fn, from_end)) for name, fn in mRNA_experiments]

    plot_to = 2000
    fig_cumulative, ax_cumulative = plt.subplots()

    xs = np.arange(-49, plot_to)
    if from_end:
        xs = -xs

    for (name, counts_generator), color in zip(experiments, colors):
        print name

        expected_counts = positions.PositionCounts(100000, 100, 100, dtype=float)
        actual_counts = positions.PositionCounts(100000, 100, 100, dtype=float)

        for counts in counts_generator:
            if not min_length <= counts.extent_length <= max_length:
                continue
            num_positions = counts.extent_length
            edge_slice = slice(-49, counts.extent_length)
            r_g = counts[edge_slice].sum()
            uniform_counts = np.ones(num_positions) * r_g / num_positions

            actual_counts[edge_slice] += counts[edge_slice]
            expected_counts[:num_positions] += uniform_counts

        print actual_counts.sum()
        print expected_counts.sum()
        ax_cumulative.plot(xs, expected_counts[-49:plot_to] / expected_counts[0], '--', color=color) 
        #ax_cumulative.plot(xs, actual_counts[-49:plot_to] / expected_counts[0], 'o', color=color, markersize=1, markeredgewidth=0, label='{0}, actual'.format(name)) 
        ax_cumulative.plot(xs, smoothed(actual_counts[-49:plot_to], 15) / expected_counts[0], '-', label='{0}'.format(name), color=color) 
        ax_cumulative.set_ylim(0.8, 1.5)

    #ax_cumulative.plot(xs, np.zeros(plot_to), 'k--')
    ax_cumulative.legend(loc='upper left', framealpha=0.5)
    if from_end:
        xlabel = 'Position relative to end of CDS'
    else:
        xlabel = 'Position relative to start of CDS'
    ax_cumulative.set_xlabel(xlabel)
    ax_cumulative.set_xlim(min(xs), max(xs))
    ax_cumulative.set_ylabel('Mapped read counts, normalized across data sets')
    #ax_cumulative.set_title('Read counts in the final {0} bases of CDSs at least {0} long'.format(min_length))

    fig_cumulative.set_size_inches(18, 12)
    plt.savefig('mRNA_bias_comparison_4.png')
    plt.savefig('mRNA_bias_comparison.pdf')

def end_bias(from_end, min_length, max_length):
    # Generators that yields arrays of counts
    def counts_from_read_positions_fn(read_positions_fn, from_end):
        gene_infos = Serialize.read_file(read_positions_fn, 'read_positions')
        for gene_name in gene_infos:
            #if gene_name == 'YLR256W':
            #    continue
            if from_end:
                counts = gene_infos[gene_name]['all'].relative_to_end
            else:
                counts = gene_infos[gene_name]['all']
            yield gene_name, counts

    experiment = ('Weinberg_mRNA', '/home/jah/projects/arlen/experiments/weinberg/mRNA/results/mRNA_read_positions.txt')
    #experiment = ('test', '/home/jah/projects/arlen/experiments/weinberg/mRNA/results/test_read_positions.txt')

    plot_to = max_length

    xs = np.arange(-49, plot_to)
    if from_end:
        xs = -xs

    name, fn = experiment
    
    excess_list = []
    bigs = []

    for gene_name, counts in counts_from_read_positions_fn(fn, from_end):
        if counts.extent_length < min_length:
            continue
            
        num_positions = counts.extent_length
        edge_slice = slice(0, counts.extent_length)
        
        r_g = counts[:counts.extent_length].sum()
        uniform_counts = np.ones(num_positions) * r_g / num_positions

        excess = (counts[:200] - uniform_counts[:200]).sum()
        excess_list.append(excess)

        if excess > 500:
            print gene_name, excess
            bigs.append(counts[:counts.extent_length])

    return excess_list, bigs

def make_L_distribution(l_g, p):
    l = np.arange(l_g + 1)
    L_distribution = (1 - p)**l * p
    L_distribution[l_g] = (1 - p)**l_g
    return L_distribution

def make_P_distribution(l_g, p):
    l = np.arange(l_g + 1)
    uniform_factor = np.zeros(l_g + 1)
    uniform_factor[1:] = 1. / l[1:]
    L_distribution = make_L_distribution(l_g, p)
    L_distribution = L_distribution / L_distribution[1:].sum()
    sum_terms = uniform_factor * L_distribution
    P_distribution = sum_terms
    for x in range(l_g - 1, 0, -1):
        P_distribution[x] += P_distribution[x + 1]
    return P_distribution

edge_overlap = 50

def counts_from_genes(genes):
    for gene_name in sorted(genes):
        CDS_slice = slice(2 * edge_overlap, 2 * edge_overlap + genes[gene_name]['CDS_length'])
        counts = sum(genes[gene_name]['position_counts'][length][CDS_slice]
                     for length in genes[gene_name]['position_counts'])
        yield counts

def get_actual_counts(genes, min_length):
    total_actual_counts = np.zeros(100000)

    for counts in counts_from_genes(genes):
        l_g = len(counts)
        if l_g < min_length:
            continue
        counts = counts[::-1]
        total_actual_counts[:l_g] += counts

    return total_actual_counts

def get_geometric_counts(genes, p, min_length):
    total_geometric_counts = np.zeros(100000)

    for counts in counts_from_genes(genes):
        l_g = len(counts)
        if l_g < min_length:
            continue
        r_g = counts.sum()

        # Some weirdness about whether to include 0 or not in make_P_distribution
        geometric_counts = make_P_distribution(l_g, p)[1:] * r_g

        total_geometric_counts[:l_g] += geometric_counts

    return total_geometric_counts

def residuals(p, genes, min_length):
    print 'testing', p
    err = get_actual_counts(genes, min_length) - get_geometric_counts(genes, p, min_length)
    return err

def fit_p():
    read_positions_fn = '/home/jah/projects/arlen/experiments/ingolia_science/mRNA-rich-1/results/mRNA-rich-1_read_positions.txt'
    #read_positions_fn = '/home/jah/projects/arlen/experiments/ingolia_science/mRNA-rich-2/results/mRNA-rich-2_read_positions.txt'
    #read_positions_fn = '/home/jah/projects/arlen/experiments/nagalakshmi_science/RH_ori/results/RH_ori_read_positions.txt'
    genes = Serialize.read_file(read_positions_fn, 'read_positions')
    min_length = 5000
    p_lsq = leastsq(residuals, 7.8e-5, args=(genes, min_length), full_output=True)
    print p_lsq
    return p_lsq

def total_counts_given_p():
    mRNA_experiments = [
        ('Ingolia_mRNA_1', '/home/jah/projects/arlen/experiments/ingolia_science/mRNA-rich-1/results/mRNA-rich-1_read_positions.txt'),
        #('Ingolia_mRNA_2', '/home/jah/projects/arlen/experiments/ingolia_science/mRNA-rich-2/results/mRNA-rich-2_read_positions.txt'),
        #('Nagalakshmi_RH_ori', '/home/jah/projects/arlen/experiments/nagalakshmi_science/RH_ori/results/RH_ori_read_positions.txt'),
        #('Nagalakshmi_RH_bio', '/home/jah/projects/arlen/experiments/nagalakshmi_science/RH_bio/results/RH_bio_read_positions.txt'),
        #('Nagalakshmi_dT_ori', '/home/jah/projects/arlen/experiments/nagalakshmi_science/dT_ori/results/dT_ori_read_positions.txt'),
        #('Nagalakshmi_dT_bio', '/home/jah/projects/arlen/experiments/nagalakshmi_science/dT_bio/results/dT_bio_read_positions.txt'),
    ]

    experiments = [(name, counts_from_genes(Serialize.read_file(fn, 'read_positions')))
                    for name, fn in mRNA_experiments]

    min_length = 0000
    max_length = 10000
    plot_to = 2000
    fig_cumulative, ax_cumulative = plt.subplots()

    xs = np.arange(0, -plot_to, -1)

    for name, counts_generator in experiments:
        print name

        total_uniform_counts = np.zeros(100000)
        total_geometric_counts = np.zeros(100000)
        total_actual_counts = np.zeros(100000)

        for counts in counts_generator:
            l_g = len(counts)
            if not (min_length < l_g < max_length):
                continue
            r_g = counts.sum()

            uniform_counts = np.ones(l_g) / l_g * r_g
            # Some weirdness about whether to include 0 or not in make_P_distribution
            #geometric_counts = make_P_distribution(l_g, 4e-4)[1:] * r_g
            geometric_counts = make_P_distribution(l_g, 7.8e-5)[1:] * r_g
            #geometric_counts = make_P_distribution(l_g, 4.6e-5)[1:] * r_g

            counts = counts[::-1]
            
            total_actual_counts[:l_g] += counts
            total_uniform_counts[:l_g] += uniform_counts
            total_geometric_counts[:l_g] += geometric_counts

        print total_actual_counts.sum()
        print total_uniform_counts.sum()
        print total_geometric_counts.sum()
        ax_cumulative.plot(xs, smoothed(total_actual_counts[:plot_to], 5), '-', label=name + '_actual_smoothed') 
        ax_cumulative.plot(xs, total_uniform_counts[:plot_to], '-', linewidth=2, label=name + '_uniform') 
        ax_cumulative.plot(xs, total_geometric_counts[:plot_to], '-', linewidth=2, label=name + '_geometric') 

    ax_cumulative.legend()
    xlabel = 'Position relative to end'
    
    ax_cumulative.set_xlabel(xlabel)
    ax_cumulative.set_ylabel('Mapped read counts')
    #ax_cumulative.set_title('Read counts in the final {0} bases of CDSs at least {0} long'.format(min_length))

#if __name__ == '__main__':
#    read_positions_fn = '/home/jah/projects/arlen/experiments/ingolia_science/mRNA-rich-1/results/mRNA-rich-1_read_positions.txt'
#    genes = Serialize.read_file(read_positions_fn, 'read_positions')
#    tail_counts = {}
#    for n in [100]:
#        tail_counts[n] = [sum(counts[-n:]) for counts in counts_from_genes(genes)
#                          if sum(counts) > 100 and len(counts) > 100 and sum(counts[-n:]) > 0]
#    read_density = [sum(counts) / float(len(counts)) for counts in counts_from_genes(genes)
#                    if sum(counts) > 100 and len(counts) > 100 and sum(counts[-n:]) > 0]
#
#    ratios = [t / r_d for t, r_d in zip(tail_counts[100], read_density)]
#    lengths = [len(counts) for counts in counts_from_genes(genes)
#               if sum(counts) > 100 and len(counts) > 100 and sum(counts[-n:]) > 0]
#
#    log_lengths = np.log10(lengths)
#    log_ratios = np.log10(ratios)
