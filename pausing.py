from __future__ import division
import positions
import codons
import numpy as np
import itertools
from collections import defaultdict, Counter
import brewer2mpl
import matplotlib.pyplot as plt
import matplotlib.cm
import matplotlib.colors
import matplotlib.transforms
from matplotlib.backends.backend_pdf import PdfPages
import Sequencing.utilities
import Sequencing.Visualize
import itertools
import scipy.stats
import os
from pausing_cython import fast_stratified_mean_enrichments, StratifiedMeanEnrichments

igv_colors = Sequencing.Visualize.igv_colors.normalized_rgbs

CHX_bmap = brewer2mpl.get_map('PuOR', 'Diverging', 10)
dark_CHX = CHX_bmap.mpl_colors[0]
light_CHX = CHX_bmap.mpl_colors[2]
dark_noCHX = CHX_bmap.mpl_colors[-2]
light_noCHX = CHX_bmap.mpl_colors[-3]

def means_of_rest(a):
    ''' Constructs an array m whose i'th element is the mean of all the elements in a except the i'th. '''
    m = (np.sum(a) - a) / float(len(a) - 1)
    return m

def compute_pause_scores(codon_counts_dict, special_sets={}):
    experiment_names = codon_counts_dict.keys()
    
    ratios_dict = {name: defaultdict(list) for name in experiment_names}
    raw_counts_dict = {name: defaultdict(list) for name in experiment_names}

    representative_counts = codon_counts_dict.itervalues().next()
    gene_names = representative_counts.keys()

    cds_slice = slice(('start_codon', 2), 'stop_codon')

    is_specials = {name: np.vectorize(special_set.__contains__) for name, special_set in special_sets.items()}

    qualifying_genes = 0

    for gene_name in gene_names:
        all_counts = [codon_counts_dict[experiment_name][gene_name]['relaxed'][cds_slice]
                      for experiment_name in experiment_names]
        medians = np.median(all_counts, axis=1)
        
        qualifies = np.all(medians >= 2)
        if not qualifies:
            continue
        
        qualifying_genes += 1

        identities = representative_counts[gene_name]['identities'][cds_slice]
        locations = {name: is_special(identities) for name, is_special in is_specials.items()}
        
        if not special_sets:
            locations['not_special'] = slice(None)
        else:
            locations['not_special'] = ~(np.any(locations.values(), axis=0))

        for experiment_name, counts in zip(experiment_names, all_counts):          
            count_ratios = counts / means_of_rest(counts)
            
            for set_name, mask in locations.items():
                ratios_dict[experiment_name][set_name].extend(count_ratios[mask])
                raw_counts_dict[experiment_name][set_name].extend(counts[mask])

    return ratios_dict, raw_counts_dict

def order_by_mean_density(codon_counts,
                          count_type='relaxed',
                          num_before=90,
                          num_after=90,
                         ):
    cds_slice = slice(('start_codon', 2), 'stop_codon')

    means = {}
    for gene_name in codon_counts:
        counts = codon_counts[gene_name][count_type][cds_slice]
        
        length = len(counts)
        if length < num_before + num_after + 1:
            mean = -1
        else:
            mean = np.mean(counts[num_before:length - num_after])
            
        means[gene_name] = mean

    sorted_names = sorted(codon_counts.keys(),
                          key=lambda n: (means[n], n),
                          reverse=True,
                         )
    means = [means[name] for name in sorted_names]

    return sorted_names, means

def get_highly_expressed_gene_names(codon_counts_dict,
                                    min_mean=1,
                                    min_median=0,
                                    count_type='relaxed',
                                    num_before=90,
                                    num_after=90,
                                   ):
    all_gene_names = codon_counts_dict.itervalues().next().keys()
    high_gene_names = []
    
    gene_means = defaultdict(list)
    gene_medians = defaultdict(list)
    
    cds_slice = slice(('start_codon', 2), 'stop_codon')
    
    def all_high(gene_name):
        to_return = True

        for exp_name in codon_counts_dict:
            counts = codon_counts_dict[exp_name][gene_name][count_type][cds_slice]

            if len(counts) < num_before + num_after + 1:
                return False
            
            gene_mean = np.mean(counts[num_before:-num_after])
            gene_median = np.median(counts[num_before:-num_after])

            gene_means[exp_name].append(gene_mean)
            gene_medians[exp_name].append(gene_median)

            if gene_mean <= min_mean or gene_median < min_median:
                to_return = False

        return to_return

    high_gene_names = [name for name in all_gene_names if all_high(name)]

    return high_gene_names, gene_means, gene_medians

def metacodon_around_pauses(codon_counts,
                            relevant_at_pause,
                            not_allowed_at_offset,
                            gene_names,
                            keep_count_context=False,
                            TEs=defaultdict(int),
                            count_type='relaxed',
                            num_before=90,
                            num_after=90,
                           ):
    old_settings = np.seterr(all='raise')

    cds_slice = slice(('start_codon', 2), 'stop_codon')
    
    counts_around_list = []
    ratios_around_list = []
    #codons_around_list = []
    nucleotides_around_list = []
    TE_list = []
    
    def is_relevant(position, codons):
        if codons[position] not in relevant_at_pause:
            return False

        for offset, not_allowed in not_allowed_at_offset.items():
            if codons[position + offset] in not_allowed:
                return False
        
        return True
    
    for gene_name in gene_names:
        counts = codon_counts[gene_name][count_type][cds_slice]
        codons = codon_counts[gene_name]['identities'][cds_slice]

        gene_TE = TEs[gene_name]

        median = np.median(counts)
        mean = np.mean(counts)

        if len(counts) < num_before + num_after + 1:
            raise ValueError
            
        #denominators = means_of_rest(counts)
        denominators = np.mean(counts[num_before:-num_after])

        if denominators == 0:
            raise ValueError(gene_name)
            
        ratios = np.true_divide(counts, denominators)
            
        for position in range(num_before, len(counts) - num_after):
            if is_relevant(position, codons):
                around_slice = slice(position - num_before, position + num_after + 1)
                
                if keep_count_context:
                    counts_around = counts[around_slice]
                    ratios_around = ratios[around_slice]
                else:
                    counts_around = counts[position]
                    ratios_around = ratios[position]
                
                codons_around = codons[around_slice]
                nucleotides_around = np.array(list(''.join(codons_around)))
                
                counts_around_list.append(counts_around)
                ratios_around_list.append(ratios_around)
                #codons_around_list.append(codons_around)
                nucleotides_around_list.append(nucleotides_around)
                TE_list.append(gene_TE)
                    
    around_lists = {'counts': np.asarray(counts_around_list),
                    'ratios': np.asarray(ratios_around_list),
                    'TEs': np.asarray(TE_list),
                    #'codons': np.asarray(codons_around_list),
                    'nucleotides': np.asarray(nucleotides_around_list),
                    'num_before': num_before,
                    'num_after': num_after,
                   }
    
    np.seterr(**old_settings)
    
    return around_lists

def compute_binned_base_compositions(around_lists, masks):
    binned_base_compositions = {}
    for bin_name, mask in masks.items():
        fractions = {}
        for offset in range(-90, 93):
            column = around_lists['nucleotides'][:, 90 + offset][mask]
            counts = Counter(column)
            fractions[offset] = {base: counts[base] / float(len(column)) for base in 'TCAG'}
            
        binned_base_compositions[bin_name] = fractions
        
    return binned_base_compositions

def compute_stratified_mean_enrichments(around_lists):
    stratified_mean_enrichments = {}
    
    num_before = around_lists['num_before']
    num_after = around_lists['num_after']
    
    ratios = around_lists['ratios']
    
    masks = make_base_identity_masks(around_lists)
    
    relevant_offsets = range(-(num_before * 3), (num_after + 1) * 3)
    for offset in relevant_offsets:
        stratified_mean_enrichments[offset] = {}
        for base in 'TCAG':
            mask = masks[offset][base]
            masked_ratios = ratios[mask]
            if len(masked_ratios) > 0:
                stratified_mean_enrichments[offset][base] = np.mean(masked_ratios)
            else:
                stratified_mean_enrichments[offset][base] = 0
    
    relevant_offsets = range(-16, 0) + range(0, 3) + range(3, 13)
    for offsets in itertools.combinations(relevant_offsets, 2):
        first_offset, second_offset = offsets
        stratified_mean_enrichments[offsets] = {}
        for bases in itertools.product('TCAG', repeat=2):
            first_base, second_base = bases
            mask = masks[first_offset][first_base] & masks[second_offset][second_base]
            masked_ratios = ratios[mask]
            if len(masked_ratios) > 0:
                stratified_mean_enrichments[offsets][bases] = np.mean(masked_ratios)
            else:
                stratified_mean_enrichments[offsets][bases] = 0
                
    relevant_codons = [(i, i + 1, i + 2) for i in range(-(num_before * 3), (num_after + 1) * 3, 3)]
    for codon_positions in relevant_codons:
        stratified_mean_enrichments[codon_positions] = {}
        for codon_id in codons.non_stop_codons:
            p1, p2, p3 = codon_positions
            b1, b2, b3 = codon_id
            mask = masks[p1][b1] & masks[p2][b2] & masks[p3][b3]
            masked_ratios = ratios[mask]
            if len(masked_ratios) > 0:
                stratified_mean_enrichments[codon_positions][codon_id] = np.mean(masked_ratios)
            else:
                stratified_mean_enrichments[codon_positions][codon_id] = 0

    relevant_P_starts = range(-(num_before * 3), num_after * 3, 3)
    for P_start in relevant_P_starts:
        print 'starting', P_start
        P_p0, P_p1, P_p2 = P_start + 0, P_start + 1, P_start + 2
        A_p0, A_p1, A_p2 = P_start + 3, P_start + 4, P_start + 5
        A_masks = {}
        P_masks = {}
        for codon_id in codons.non_stop_codons:
            b0, b1, b2 = codon_id
            P_masks[codon_id] = masks[P_p0][b0] & masks[P_p1][b1] & masks[P_p2][b2]
            A_masks[codon_id] = masks[A_p0][b0] & masks[A_p1][b1] & masks[A_p2][b2]
        
        means = {}
        for P_codon_id in codons.non_stop_codons:
            for A_codon_id in codons.non_stop_codons:
                mask = A_masks[A_codon_id] & P_masks[P_codon_id]
                masked_ratios = ratios[mask]
                if len(masked_ratios) > 0:
                    means[P_codon_id, A_codon_id] = np.mean(masked_ratios)
                else:
                    means[P_codon_id, A_codon_id] = 0
        stratified_mean_enrichments[(P_p0, P_p1, P_p2), (A_p0, A_p1, A_p2)] = means
                    
    stratified_mean_enrichments['all'] = np.mean(ratios)
    
    return stratified_mean_enrichments

def make_base_identity_masks(around_lists):
    num_before = around_lists['num_before']
    num_after = around_lists['num_after']
    
    masks = {}
    for offset in range(-(num_before * 3), (num_after + 1) * 3):
        masks[offset] = {base: around_lists['nucleotides'][:, num_before * 3 + offset] == base for base in 'TCAG'}
    return masks

def split_into_bins(around_lists, quantize_at, num_quantiles): 
    ratios_around = around_lists['ratios']
    num_before = around_lists['num_before']

    quantiles = scipy.stats.mstats.mquantiles(ratios_around[:, num_before + quantize_at],
                                              prob=np.linspace(0, 1, num_quantiles + 1),
                                             )
    first_nonzero = quantiles.nonzero()[0][0]
    if first_nonzero > 0:
        quantiles = quantiles[first_nonzero - 1:]

    boundaries = zip(quantiles, quantiles[1:])
    bins = [(ratios_around[:, num_before + quantize_at] >= start) & (ratios_around[:, num_before + quantize_at] < end) for start, end in boundaries]
    binned = [ratios_around[mask] for mask in bins]
    
    bins = {i: mask for i, mask in enumerate(bins)}
    bins['all'] = slice(None)
    
    return bins, binned, quantiles

def split_into_TE_bins(around_lists, num_quantiles): 
    ratios_around = around_lists['ratios']

    TEs = around_lists['TEs']

    quantiles = scipy.stats.mstats.mquantiles(TEs,
                                              prob=np.linspace(0, 1, num_quantiles + 1),
                                             )

    boundaries = zip(quantiles, quantiles[1:])
    bins = [(TEs >= start) & (TEs < end) for start, end in boundaries]
    binned = [ratios_around[mask] for mask in bins]
    
    bins = {i: mask for i, mask in enumerate(bins)}
    bins['all'] = slice(None)
    
    return bins, binned, quantiles

def plot_cdfs(ratios_dict, ordinal_colors=False, sorting_key=lambda x: x):
    ''' Plots 1 - cdf of pause scores that have been stratified into different sets. '''
    if ordinal_colors:
        colors = get_jet_colors(len(ratios_dict))
    else:
        bmap = brewer2mpl.get_map('Set1', 'qualitative', 9)
        colors = bmap.mpl_colors[:5] + bmap.mpl_colors[6:] + ['black']
        colors = colors + colors

    fig, ax = plt.subplots(figsize=(16, 12))

    color_iter = iter(colors)
    for name in sorted(ratios_dict, key=sorting_key):
        ratios = ratios_dict[name]
        for key in ratios:
            if ratios[key] == []:
                continue
            sorted_values, cumulative = Sequencing.utilities.empirical_cdf(ratios[key])
            ax.semilogy(sorted_values, 1 - cumulative, color=color_iter.next(), label=name + '_' + key)

    ax.legend(framealpha=0.5)
    ax.set_ylabel('Fraction of codons with enrichment >= x')
    ax.set_xlabel('Ratio of reads at codon to (mean across coding sequence excluding codon)')
    
    return ax
   
def scatter_enrichments(ratios_dict, x_name, y_name, set_name, draw_labels=True, ax=None):
    if ax == None:
        fig, ax = plt.subplots(figsize=(12, 12))

    ratios = {'x': ratios_dict[x_name][set_name],
              'y': ratios_dict[y_name][set_name],
             }

    smallest_nonzeros = {key: min(r for r in rs if r > 0) for key, rs in ratios.items()}  
    sanitized = {key: np.log2(np.maximum(smallest_nonzeros[key], ratios[key])) for key in smallest_nonzeros}

    Sequencing.Visualize.enhanced_scatter(sanitized['x'],
                                          sanitized['y'],
                                          ax,
                                          do_fit=False,
                                          show_p_value=False,
                                          hists_height=0.2,
                                          #color_by_density=False,
                                         )
    smallest = min(min(sanitized['x']), min(sanitized['y']))
    largest = max(max(sanitized['x']), max(sanitized['y']))

    ax.set_xlim(smallest - 0.2, largest + 0.2)
    ax.set_ylim(smallest - 0.2, largest + 0.2)

    Sequencing.Visualize.draw_diagonal(ax)
    
    if draw_labels:
        ax.set_xlabel('log2({0})'.format(x_name))
        ax.set_ylabel('log2({0})'.format(y_name))
 
def scatter_all_codon_id_enrichments(stratified_lists_dict,
                                     x_name,
                                     y_name,
                                     codon_positions,
                                     just_hists=False,
                                    ):
    fig, axs = plt.subplots(16, 4, figsize=(4 * 4, 16 * 4))

    for codon_id, ax in zip(codons.non_stop_codons, axs.flatten()):
        lists = {'x': stratified_lists_dict[x_name][codon_positions][codon_id],
                 'y': stratified_lists_dict[y_name][codon_positions][codon_id],
                }


        if just_hists:
            sanitized = {key: np.log2(np.maximum(2**-7, lists[key])) for key in lists}
            kwargs = {'bins': np.linspace(-7, 7, 50),
                      'histtype': 'step',
                     }
            ax.hist(sanitized['x'], **kwargs)
            ax.hist(sanitized['y'], **kwargs)
        else:
            smallest_nonzeros = {key: min(r for r in rs if r > 0) for key, rs in lists.items()}  
            sanitized = {key: np.log2(np.maximum(smallest_nonzeros[key], lists[key])) for key in smallest_nonzeros}
            
            Sequencing.Visualize.enhanced_scatter(sanitized['x'],
                                                  sanitized['y'],
                                                  ax,
                                                  do_fit=False,
                                                  show_p_value=False,
                                                  #color_by_density=False,
                                                  #hists_height=0.2,
                                                 )

            smallest = min(min(sanitized['x']), min(sanitized['y']))
            largest = max(max(sanitized['x']), max(sanitized['y']))

            ax.set_xlim(smallest - 0.2, largest + 0.2)
            ax.set_ylim(smallest - 0.2, largest + 0.2)

            Sequencing.Visualize.draw_diagonal(ax)
    
        #ax.set_xlabel('log2({0})'.format(x_name))
        #ax.set_ylabel('log2({0})'.format(y_name))

        ax.set_title('{0}: {1:,}'.format(codon_id, len(sanitized['x'])))

def get_jet_colors(n):
    jet_colors = [matplotlib.cm.jet(x) for x in np.linspace(0, 1, n)]
    return jet_colors

def plot_binned_base_compositions(binned_base_compositions, normalized=False, show_A_site=True):
    fig, axs = plt.subplots(4, 1, figsize=(16, 12 * 4))
    
    jet_colors = get_jet_colors(len(binned_base_compositions) - 1)
    
    for base, ax in zip('TCAG', axs):
        xs = np.arange(-25, 26)

        for bin_i in range(len(binned_base_compositions) - 1):
            numerators = [binned_base_compositions[bin_i][x][base] for x in xs]
            denominators = [max(1e-1, binned_base_compositions['all'][x][base]) for x in xs]
            if normalized:
                ys = np.true_divide(numerators, denominators)
            else:
                ys = np.array(numerators)
                
            if not show_A_site:
                masked_xs = list(xs)
                masked_ys = list(ys)
                for i in range(25, 25 + 3):
                    masked_xs[i] = None
                    masked_ys[i] = None
            else:
                masked_xs = xs
                masked_ys = ys

            ax.plot(masked_xs, masked_ys, '.-', label=bin_i, color=jet_colors[bin_i])

        ax.legend(framealpha=0.5)
        ax.set_title(base)
        
        if normalized:
            ax.set_ylim(0.2, 1.8)
        else:
            ax.set_ylim(ymax=0.5)

        five_prime_edge = -15
        three_prime_edge = five_prime_edge + 28 - 1
        
        for x in xs:
            if x in [five_prime_edge, three_prime_edge]:
                alpha = 0.7
            else:
                alpha = 0.1
            ax.axvline(x, color='black', alpha=alpha)
            
        ax.set_xlim(min(xs), max(xs))
        
        tRNA_sites = [('A', 0, 'red'),
                      ('P', -3, 'blue'),
                      ('E', -6, 'green'),
                     ]

        for site, position, color in tRNA_sites: 
            ax.axvspan(position - 0.5, position + 2.5, color=color, alpha=0.1)
            ax.annotate(site,
                        xy=(position + 1, 1),
                        xycoords=('data', 'axes fraction'),
                        xytext=(0, -25),
                        textcoords='offset points',
                        horizontalalignment='center',
                        size=20,
                       )

def plot_nucleotide_enrichments(enrichments,
                                plot_A_site=True,
                                min_x=-30,
                                max_x=32, 
                                ax=None,
                                flip=True,
                                dense_lines=True,
                                marker_size=6,
                                line_width=1,
                                minimal_ticks=False,
                                legend_kwargs={'loc': 'upper right'},
                               ):
    if plot_A_site:
        xs = range(min_x, 0) + range(0, 3) + range(3, max_x + 1)
    else:
        xs = range(min_x, 0) + [None, None, None] + range(3, max_x + 1)

    if ax == None:
        fig, ax = plt.subplots(figsize=(16, 12))
    
    for base in 'TCAG':
        ys = []
        for x in xs:
            if x == None:
                y = None
            else:
                y = enrichments['nucleotide', x, base]
            ys.append(y)
        ax.plot(xs, ys, '.-',
                color=igv_colors[base],
                label=base,
                markersize=marker_size,
                lw=line_width,
               )

    ax.legend(framealpha=0.5, **legend_kwargs)

    five_prime_edge = -15
    three_prime_edge = five_prime_edge + 28 - 1

    for x in range(min_x, max_x + 1):
        lw = 1
        if x in [five_prime_edge, three_prime_edge]:
            alpha = 1
            lw = 1.5
        elif dense_lines and x % 3 == 0:
            alpha = 0.2
        elif dense_lines:
            alpha = 0.05
        else:
            alpha = 0

        ax.axvline(x, color='black', alpha=alpha, lw=lw)
    
    if minimal_ticks:
        xticks = [-15, 0, 12]
    else:
        xticks = range(min_x, max_x + 1, 3)

    ax.set_xticks(xticks)
    ax.set_xticklabels(map(str, xticks))

    ax.set_xlim(min_x, max_x)

    tRNA_sites = [('A', 0, 'red'),
                  ('P', -3, 'blue'),
                  ('E', -6, 'green'),
                 ]

    for site, position, color in tRNA_sites: 
        ax.axvspan(position - 0.5, position + 2.5, color=color, alpha=0.1, lw=0)
        ax.annotate(site,
                    xy=(position + 1, 1),
                    xycoords=('data', 'axes fraction'),
                    xytext=(0, -25),
                    textcoords='offset points',
                    horizontalalignment='center',
                    size=20,
                   )

    ax.set_ylabel('Mean relative enrichment')
    ax.set_xlabel('Offset (nucleotides)')
        
    if flip:
        ax.invert_xaxis()
        flipped_labels = [str(int(-x)) for x in ax.get_xticks()]
        ax.set_xticklabels(flipped_labels)
        
def plot_dinucleotide_effects(stratified_mean_enrichments, relevant_offsets, min_difference, fancy=True, log=True):
    baseline = stratified_mean_enrichments['all']
    actuals = []
    expecteds = []
    sites = []
    for first_offset, second_offset in itertools.combinations(relevant_offsets, 2):
        for first_base, second_base in itertools.product('TCAG', repeat=2):
            expected = stratified_mean_enrichments[first_offset][first_base] * stratified_mean_enrichments[second_offset][second_base] / baseline**2
            actual = stratified_mean_enrichments[first_offset, second_offset][first_base, second_base] / baseline

            expecteds.append(expected)
            actuals.append(actual)
            sites.append('({0}, {1}) x ({2}, {3})'.format(first_offset, first_base, second_offset, second_base))

    fig, ax = plt.subplots(figsize=(16, 12))

    if log:
        xs = np.log2(expecteds)
        ys = np.log2(actuals)
    else:
        xs = np.array(expecteds)
        ys = np.array(actuals)
    
    if fancy:
        Sequencing.Visualize.enhanced_scatter(xs, ys, ax, do_fit=False)
    else:
        ax.scatter(xs, ys, s=2)
    
    ax.set_xlabel('Expected if multiplicative')
    ax.set_ylabel('Actual')
    
    ax.set_aspect(1.)

    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    lower = min(x_min, y_min)
    upper = max(x_max, y_max)

    ax.plot([lower, upper], [lower, upper], color='black', alpha=0.2, scalex=False, scaley=False);

    label_scatter_plot(ax, xs, ys, sites, min_difference)

def plot_dicodon_effects(enrichments,
                         min_difference,
                         special_Ps,
                         special_As,
                         fancy=True,
                         log=True,
                         offset=0,
                        ):
    baseline = 1
    actuals = []
    expecteds = []
    sites = []
    colors = []

    black = matplotlib.colors.colorConverter.to_rgba('black', alpha=0.1)
    red = matplotlib.colors.colorConverter.to_rgba('red', alpha=1)

    for P_codon_id in codons.non_stop_codons:
        for A_codon_id in codons.non_stop_codons:
            P_value = enrichments['codon', offset - 1, P_codon_id]
            A_value = enrichments['codon', offset, A_codon_id]
            expected =  P_value * A_value
            actual = enrichments['dicodon', offset, (P_codon_id, A_codon_id)] / baseline

            expecteds.append(expected)
            actuals.append(actual)
            label = '{0}-{1}, ({2}-{3}), {4}'.format(P_codon_id,
                                                    A_codon_id,
                                                    codons.forward_table[P_codon_id],
                                                    codons.forward_table[A_codon_id],
                                                    enrichments['dicodon_occurences', 0, (P_codon_id, A_codon_id)],
                                                   )
            sites.append(label)
            if A_codon_id in special_As and P_codon_id in special_Ps:
                color = red
            else:
                color = black
            colors.append(color)

    fig, ax = plt.subplots(figsize=(16, 12))

    if log:
        xs = np.log2(expecteds)
        ys = np.log2(actuals)
    else:
        xs = np.array(expecteds)
        ys = np.array(actuals)
    
    if fancy:
        Sequencing.Visualize.enhanced_scatter(xs, ys, ax, do_fit=False)
    else:
        ax.scatter(xs, ys, c=colors, s=8, linewidths=(0,))
    
    ax.set_xlabel('Expected if multiplicative')
    ax.set_ylabel('Actual')
    
    ax.set_aspect(1.)

    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    lower = min(x_min, y_min)
    upper = max(x_max, y_max)

    ax.plot([lower, upper], [lower, upper], color='black', alpha=0.2, scalex=False, scaley=False);

    to_label = np.abs(xs - ys) > min_difference
    to_label = ys > min_difference
    label_scatter_plot(ax, xs, ys, sites, to_label)
    
def plot_codon_effects(stratified_mean_enrichments, relevant_codons, fancy=True, log=True, ax=None):
    baseline = stratified_mean_enrichments['all']
    actuals = []
    expecteds = []
    labels = []
    for codon in relevant_codons:
        for bases in codons.non_stop_codons:
            expected = stratified_mean_enrichments[codon[0]][bases[0]] * stratified_mean_enrichments[codon[1]][bases[1]] * stratified_mean_enrichments[codon[2]][bases[2]] / baseline**3
            actual = stratified_mean_enrichments[codon][bases] / baseline

            expecteds.append(expected)
            actuals.append(actual)

            bases = ''.join(bases)
            labels.append('{0}: {1} ({2})'.format(codon[0] // 3, bases, codons.forward_table[bases]))

    if ax == None:
        fig, ax = plt.subplots(figsize=(16, 12))

    if log:
        xs = np.log2(expecteds)
        ys = np.log2(actuals)
    else:
        xs = expecteds
        ys = actuals

    if fancy:
        Sequencing.Visualize.enhanced_scatter(xs, ys, ax, do_fit=False, show_p_value=False, draw_diagonal=True)
    else:
        ax.scatter(xs, ys, s=2)
    
    ax.set_xlabel('Expected if multiplicative')
    ax.set_ylabel('Actual')
    
    smallest = min(min(xs), min(ys))
    largest = max(max(xs), max(ys))
    
    ax.set_xlim(smallest - 0.1, largest + 0.1)
    ax.set_ylim(smallest - 0.1, largest + 0.1)
    ax.set_aspect(1.)

    to_label = [i for i, d in enumerate(xs - ys) if abs(d) > 2.5]
    label_scatter_plot(ax, xs, ys, labels, to_label)
    
def label_scatter_plot(ax, xs, ys, labels, to_label,
                       vector='orthogonal',
                       initial_distance=50,
                       arrow_alpha=0.2,
                       manual_ratios=None,
                       manual_alignments=None,
                       text_kwargs={'size': 10},
                      ):
    def attempt_text(x, y, site, distance):
        if vector == 'orthogonal':
            x_offset = np.sign(x - y) * distance
            y_offset = -np.sign(x - y) * distance
            ha = 'center'
            va = 'top' if y_offset < 0 else 'bottom'
        elif vector == 'radial':
            norm = np.linalg.norm([x, y])
            x_offset = x * distance / norm
            y_offset = y * distance / norm
            ha, va = manual_alignments
        elif vector == 'sideways':
            x_offset = distance
            y_offset = 0
            ha, va = 'center', 'top'
        elif vector == 'manual':
            x_ratio, y_ratio = manual_ratios
            ha, va = manual_alignments
            x_offset = distance * x_ratio
            y_offset = distance * y_ratio

        text = ax.annotate(site,
                           xy=(x, y),
                           xycoords=('data', 'data'),
                           xytext=(x_offset, y_offset),
                           textcoords='offset points',
                           ha=ha,
                           va=va,
                           **text_kwargs)
        ax.figure.canvas.draw()

        return text, text.get_window_extent(), (x, y, x_offset, y_offset)

    ax.figure.canvas.draw()
    starting_labels = [ax.xaxis.get_label(), ax.yaxis.get_label()] + ax.get_yticklabels() + ax.get_xticklabels()
    bboxes = [label.get_window_extent() for label in starting_labels]

    tuples = itertools.izip(np.asarray(xs)[to_label],
                            np.asarray(ys)[to_label],
                            np.asarray(labels)[to_label],
                            )
    for x, y, label in tuples:
        distance = initial_distance
        text, bbox, coords = attempt_text(x, y, label, distance)
        while any(bbox.fully_overlaps(other_bbox) for other_bbox in bboxes):
            text.remove()
            distance += 10
            text, bbox, coords = attempt_text(x, y, label, distance)
            if distance >= 500:
                break
        
        x, y, x_offset, y_offset = coords
        ax.annotate('',
                    xy=(x, y),
                    xycoords=('data', 'data'),
                    xytext=(x_offset, y_offset),
                    textcoords=('offset points', 'offset points'),
                    arrowprops={'arrowstyle': '->', 'alpha': 0.5},
                   )

        bboxes.append(bbox)

def label_enrichment_across_conditions_plot(ax, start_ys, end_ys, labels, xs, label_offset=30, size=10):
    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.renderer

    def attempt_text(y, label, distance, side, y_sign):
        if side == 'left':
            x = 0
            #x_offset = -distance
            x_offset = -label_offset
            y_offset = distance * y_sign
            ha = 'right'
        elif side == 'right':
            x = max(xs)
            #x_offset = 1 + distance
            x_offset = label_offset
            y_offset = distance * y_sign
            ha = 'left'
            
        text = ax.annotate(label,
                           xy=(x, y),
                           xycoords=('data', 'data'),
                           xytext=(x_offset, y_offset),
                           textcoords=('offset points', 'offset points'),
                           ha=ha,
                           va='center',
                           size=size,
                          )
        ax.figure.canvas.draw()
        return text, text.get_window_extent(renderer), (x, y, x_offset, y_offset)

    starting_labels = []
    #starting_labels = [ax.xaxis.get_label(), ax.yaxis.get_label()] + ax.get_yticklabels()
    bboxes = [label.get_window_extent(renderer) for label in starting_labels]

    left_tuples = sorted(zip(start_ys, labels, itertools.repeat('left')), reverse=True)
    right_tuples = sorted(zip(end_ys, labels, itertools.repeat('right')), reverse=True)
   
    middle_index = len(left_tuples) // 2
    up_tuples = left_tuples[middle_index::-1] + right_tuples[middle_index::-1]
    down_tuples = left_tuples[middle_index + 1:] + right_tuples[middle_index + 1:]

    up_tuples = [(y, label, side, 1) for y, label, side in up_tuples]
    down_tuples = [(y, label, side, -1) for y, label, side in down_tuples]

    for y, label, side, y_sign in up_tuples + down_tuples:
        distance = 0
        text, bbox, coords = attempt_text(y, label, distance, side, y_sign)
        while any(bbox.fully_overlaps(other_bbox) for other_bbox in bboxes):
            text.remove()
            distance += 10
            text, bbox, coords = attempt_text(y, label, distance, side, y_sign)
            if distance >= 500:
                break
        
        x, y, x_offset, y_offset = coords
        ax.annotate('',
                    xy=(x, y),
                    xycoords=('data', 'data'),
                    xytext=(x_offset, y_offset),
                    textcoords=('offset points', 'offset points'),
                    arrowprops={'arrowstyle': '->', 'alpha': 0.5},
                   )

        bboxes.append(bbox)
            
def load_tRNA_copy_numbers(kind='copy number'):
    if kind == 'copy number':
        expected_times = {}
        fn = '{0}/projects/translation_elongation/data/S.cer.code'.format(os.environ['HOME'])
        for line in open(fn):
            codon, _, copy_number, wobble = line.strip().split()
            expected_times[codon] = 1. / int(copy_number)
    elif kind == 'wobble':
        expected_times = {}
        fn = '{0}/projects/translation_elongation/data/S.cer.code'.format(os.environ['HOME'])
        for line in open(fn):
            codon, _, copy_number, wobble = line.strip().split()
            expected_times[codon] = 1. / (int(copy_number) * float(wobble))
    elif kind == 'tAI':
        expected_times = {}
        fn = '{0}/projects/translation_elongation/data/tAI.txt'.format(os.environ['HOME'])
        for line in open(fn):
            codon, anticodon, tAI = line.strip().split()
            expected_times[codon] = 1. / float(tAI)
    elif kind == 'premal tAI':
        expected_times = {}
        fn = '{0}/projects/translation_elongation/data/premal_tAI.txt'.format(os.environ['HOME'])
        for line in open(fn):
            codon, anticodon, _, _, _, tAI = line.strip().split()
            expected_times[codon] = 1. / float(tAI)
    return expected_times

def assign_colors_by_values(codon_to_ranks):
    colors = get_jet_colors(len(codon_to_ranks))
    codon_to_color = {codon: colors[codon_to_ranks[codon][0]] for codon in codon_to_ranks}
    return codon_to_color

def build_codon_to_ranks(stratified_mean_enrichments_dict, position, name_order):
    all_ranks = []
    all_values = []
    for name in name_order:
        enrichments = stratified_mean_enrichments_dict[name]
        values = [enrichments['codon', position, codon_id] for codon_id in codons.non_stop_codons]
        ranks = np.array(values).argsort().argsort()
        all_values.append(values)
        all_ranks.append(ranks)

    all_values = np.array(all_values).T
    all_ranks = np.array(all_ranks).T

    codon_to_values = {codon: values for codon, values in zip(codons.non_stop_codons, all_values)}
    codon_to_ranks = {codon: ranks for codon, ranks in zip(codons.non_stop_codons, all_ranks)}

    return codon_to_ranks, codon_to_values

def build_dicodon_to_ranks(stratified_mean_enrichments_dict,
                           position,
                           name_order,
                           allowed_dicodons=None,
                          ):
    all_ranks = []
    all_values = []
    if allowed_dicodons == None:
        allowed_dicodons = list(itertools.product(codons.non_stop_codons, repeat=2))
    for name in name_order:
        enrichments = stratified_mean_enrichments_dict[name]
        values = [enrichments['dicodon', position, dicodon_id] for dicodon_id in allowed_dicodons]
        ranks = np.array(values).argsort().argsort()
        all_values.append(values)
        all_ranks.append(ranks)

    all_values = np.array(all_values).T
    all_ranks = np.array(all_ranks).T

    dicodon_to_values = {dicodon: values for dicodon, values in zip(allowed_dicodons, all_values)}
    dicodon_to_ranks = {dicodon: ranks for dicodon, ranks in zip(allowed_dicodons, all_ranks)}

    return dicodon_to_ranks, dicodon_to_values

def plot_dicodon_enrichments_across_conditions(stratified_mean_enrichments_dict,
                                               position,
                                               name_order,
                                               allowed_at_P,
                                               allowed_at_A,
                                               by_rank=False,
                                               highlight_movement=True,
                                               rank_or_log='rank',
                                               log_scale=False,
                                              ):
    allowed_dicodons = list(itertools.product(allowed_at_P, allowed_at_A))
    dicodon_to_ranks, dicodon_to_values = build_dicodon_to_ranks(stratified_mean_enrichments_dict, position, name_order, allowed_dicodons)
    dicodon_to_color = assign_colors_by_values(dicodon_to_ranks)

    rank_deltas = {dicodon: dicodon_to_ranks[dicodon][0] - dicodon_to_ranks[dicodon][-1]
                   for dicodon in allowed_dicodons}
    
    log_deltas = {dicodon: np.log2(max(2**-8, dicodon_to_values[dicodon][0])) - np.log2(max(2**-8, dicodon_to_values[dicodon][-1]))
                  for dicodon in allowed_dicodons}
    
    fig, ax = plt.subplots(figsize=(16, 12))

    xs = range(len(name_order))

    for dicodon in allowed_dicodons:
        if by_rank:
            ys = dicodon_to_ranks[dicodon]
        else:
            ys = np.maximum(2**-8, dicodon_to_values[dicodon])

        if highlight_movement:
            if rank_or_log == 'rank':
                alpha = min(1, abs(rank_deltas[dicodon]) / float(5000))
                width = abs(rank_deltas[dicodon]) / float(5000)
            elif rank_or_log == 'log':
                alpha = min(1, abs(log_deltas[dicodon]) / 2)
                width = abs(log_deltas[dicodon] / 4)
        else:
            alpha = 1
            width = 1

        ax.plot(xs, ys, '-', alpha=alpha, lw=width, color=dicodon_to_color[dicodon])
    
    if by_rank:
        ax.set_yticks([])
        ax.set_ylim(-0.5, len(allowed_dicodons) - 0.5)
    else:
        if log_scale:
            ax.set_yscale('log', basey=2)


def plot_enrichments_across_conditions(enrichments,
                                       position,
                                       name_order,
                                       highlight_movement=True,
                                       force_highlight=set(),
                                       by_rank=False,
                                       label_rules=(0, 'rank', abs),
                                       log_scale=True,
                                       force_label=set(),
                                       ax=None,
                                       force_ylims=None,
                                       print_deltas=True,
                                       label_offset=32,
                                       label_size=10,
                                       marker_size=6,
                                       heat_exception=False,
                                       ylabel_size=14,
                                      ):
    if ax == None:
        fig, ax = plt.subplots(figsize=(16, 12))
    else:
        fig = ax.get_figure()

    codon_to_ranks, codon_to_values = build_codon_to_ranks(enrichments, position, name_order)
    codon_to_color = assign_colors_by_values(codon_to_ranks)

    rank_deltas = {codon_id: codon_to_ranks[codon_id][0] - codon_to_ranks[codon_id][-1]
                   for codon_id in codons.non_stop_codons}
    
    log_deltas = {codon_id: np.log2(codon_to_values[codon_id][0] / codon_to_values[codon_id][-1])
                  for codon_id in codons.non_stop_codons}

    num_to_label, rank_or_log, transform = label_rules
    if rank_or_log == 'rank':
        deltas = rank_deltas
    elif rank_or_log == 'log':
        deltas = log_deltas

    xs = range(len(name_order))
    if heat_exception:
        xs = [0, 4]
        buffer_factor = 2
    else:
        buffer_factor = 1
        
    for codon_id in codons.non_stop_codons:
        if highlight_movement:
            if rank_or_log == 'rank':
                alpha = min(1, min(30, abs(rank_deltas[codon_id])) / float(50))
                width = min(30, abs(rank_deltas[codon_id])) / float(15)
            elif rank_or_log == 'log':
                alpha = min(1, abs(log_deltas[codon_id]))
                width = 1
                width = abs(log_deltas[codon_id] / 1.5)
        elif force_highlight:
            if codon_id in force_highlight:
                alpha = 1
                width = 2
            else:
                alpha = 0.1
                width = 0.5
        else:
            alpha = 1
            width = 1
                  
        if by_rank:
            ys = codon_to_ranks[codon_id]
        else:
            ys = codon_to_values[codon_id]

        ax.plot(xs, ys, '.-',
                color=codon_to_color[codon_id],
                alpha=alpha,
                lw=width,
                markersize=marker_size,
               )
             
    sorted_deltas = sorted(deltas, key=lambda c: transform(deltas[c]), reverse=True)

    if num_to_label > 0 and print_deltas:
        if rank_or_log == 'rank':
            print 'Rank changes'
            for codon_id in sorted_deltas[:num_to_label]:
                print '{0}: {1:2>d}'.format(codon_id, deltas[codon_id])
        elif rank_or_log == 'log':
            print 'log_2 changes'
            for codon_id in sorted_deltas[:num_to_label]:
                print '{0}: {1:0.3f}'.format(codon_id, deltas[codon_id])
    
    start_ys  = []
    end_ys = []
    labels = []

    first_enrichments = enrichments[name_order[0]]
    last_enrichments = enrichments[name_order[-1]]

    for codon_id in set(sorted_deltas[:num_to_label]) | set(force_label):
        label = '{0} ({1})'.format(codon_id, codons.full_forward_table[codon_id])
        
        if by_rank:
            start_y = codon_to_ranks[codon_id][0]
            end_y = codon_to_ranks[codon_id][-1]
        else:
            start_y = first_enrichments['codon', position, codon_id]
            end_y = last_enrichments['codon', position, codon_id]
            
        start_ys.append(start_y)
        end_ys.append(end_y)
        labels.append(label)

    if by_rank:
        ax.set_yticks([])
        ax.set_ylim(-0.5, 61.5)
    else:
        if force_ylims:
            ax.set_ylim(force_ylims)
        if log_scale:
            ax.set_yscale('log', basey=2)

        ax.tick_params(labelright=True)

        big_bold = matplotlib.font_manager.FontProperties(size=ylabel_size, weight='bold')
        for label in ax.get_yticklabels():
            label.set_fontproperties(big_bold)

        ax.yaxis.grid(True, which='major', linestyle='-', alpha=0.6)
    
    label_enrichment_across_conditions_plot(ax, start_ys, end_ys, labels, xs,
                                            label_offset=label_offset,
                                            size=label_size,
                                           )

    ax.set_xticks(xs)
    ax.set_xticklabels(name_order, rotation=30, ha='right', size=12)
    x_span = max(xs) - min(xs)
    ax.set_xlim(min(xs) - 0.02 * buffer_factor * x_span, max(xs) + 0.02 * buffer_factor * x_span)

    return fig

def plot_correlations_across_conditions(enrichments,
                                        name_order,
                                        position=0,
                                        plot_p_values=True,
                                        plot_correlations=True,
                                        x_labels=None,
                                        tRNA_value_source='tAI',
                                       ):
    rhos = []
    ps = []

    tAIs = load_tRNA_copy_numbers(tRNA_value_source)
    tAI_values = [tAIs[codon] for codon in codons.non_stop_codons]

    name_length = max(map(len, name_order))
    for name in name_order:
        e_values = enrichments[name]['codon', position, codons.non_stop_codons]
        rho, p = scipy.stats.spearmanr(tAI_values, e_values)
        rhos.append(rho)
        ps.append(p)
        print '{0:<{name_length}}{1:10.4f}\t{2:0.1e}'.format(name, rho, p, name_length=name_length)
        
    fig, ax = plt.subplots(figsize=(16, 12))
    xs = np.arange(len(rhos))
    if plot_correlations:
        ax.plot(xs, rhos, 'o-', color='blue', label=r'$\rho$')
    ax.set_ylabel(r'Spearman $\rho$', color='blue', size=16)

    if x_labels is None:
        x_labels = name_order

    ax.set_xticks(xs)
    ax.set_xticklabels(x_labels, rotation=30, ha='right', size=14)

    ax.axhline(0, color='black', alpha=1)

    offset_string = ' ({0:+d})'.format(position) if position != 0 else ''
    title = 'Rank correlation of codon identity A-site occupancy{1} with 1 / tAI'.format(position, offset_string)
    ax.set_title(title, size=20)


    big_bold = matplotlib.font_manager.FontProperties(size=12, weight='bold')
    for label in ax.get_yticklabels():
        label.set_fontproperties(big_bold)
    
    p_ax = ax.twinx()
    p_ax.set_ylim(1e-8, 1)
    p_ax.set_yscale('log')
    p_ax.set_ylabel('p value', color='red', size=16)
    
    for label in p_ax.get_yticklabels():
        label.set_fontproperties(big_bold)
    
    if plot_p_values:
        p_ax.plot(xs, ps, 'o-', color='red', alpha=0.5, label='p value')
        p_ax.grid(alpha=0.7)

    ax.set_ylim(-0.5, 1)
    ax.set_xlim(min(xs) - 0.1, max(xs) + 0.1)

    return fig
                                       
def plot_codon_enrichments(names,
                           enrichments,
                           codons_to_highlight,
                           min_x=-30,
                           max_x=30,
                           log_scale=False,
                           force_ylims=None,
                           split_by_codon=False,
                           only_show_highlights=False,
                           sample_to_label=None,
                           flip=False,
                           ax=None,
                           mark_stalls=False,
                           legend_kwargs={'loc': 'upper right'},
                           marker_size=6,
                           line_width=1,
                           mark_active_sites=True,
                          ):

    if ax != None and ((split_by_codon and len(codons_to_highlight) > 1) or (not split_by_codon and len(names) > 1)):
        raise ValueError('Can\'t supply ax if more than one frame will be produced')
    
    bmap = brewer2mpl.get_map('Set1', 'qualitative', 9)
    colors = bmap.mpl_colors[:5] + bmap.mpl_colors[6:]
    colors_iter = itertools.cycle(iter(colors))

    all_xs = {}
    all_ys = {}

    for sample in names:
        for codon_id in codons.non_stop_codons:
            xs = np.arange(min_x, max_x + 1)
            ys = enrichments[sample]['codon', min_x:max_x + 1, codon_id]

            all_xs[sample, codon_id] = xs
            all_ys[sample, codon_id] = ys
        
    if sample_to_label == None:
        sample_to_label = {name: name for name in names}
            
    if split_by_codon:
        sample_to_color = {name: colors_iter.next() for name in names}

        if ax == None:
            fig, axs = plt.subplots(len(codons_to_highlight), 1,
                                    figsize=(16, 12 * len(codons_to_highlight)),
                                    squeeze=False,
                                   )
        else:
            fig = ax.get_figure()
            axs = np.array([ax])

        for codon_id, codon_ax in zip(codons_to_highlight, axs.flatten()):
            amino_acid = codons.full_forward_table[codon_id]
            for sample in names:
                kwargs = {'color': sample_to_color[sample],
                          'alpha': 1,
                          'label': sample_to_label[sample],
                          'markersize': marker_size,
                          'linewidth': line_width,
                         }
                xs = all_xs[sample, codon_id]
                ys = all_ys[sample, codon_id]
                codon_ax.plot(xs, ys, '.-', **kwargs)

            codon_ax.set_title('{0} ({1})'.format(codon_id, amino_acid), size=16)
            codon_ax.legend(framealpha=0.5, **legend_kwargs)
            codon_ax.axhline(1, color='black')

    else:
        codon_to_color = {codon_id: colors_iter.next() for codon_id in codons_to_highlight}
        codon_to_handle = {}

        if ax == None:
            fig, axs = plt.subplots(len(names), 1,
                                    figsize=(16, 12 * len(names)),
                                    squeeze=False,
                                   )
        else:
            fig = ax.get_figure()
            axs = np.array([ax])

        for sample, sample_ax in zip(names, axs.flatten()):
            for codon_id in codons.non_stop_codons:
                if codon_id in codons_to_highlight:
                    amino_acid = codons.full_forward_table[codon_id]

                    kwargs = {'color': codon_to_color[codon_id],
                              'alpha': 1,
                              'label': '{0} ({1})'.format(codon_id, amino_acid),
                              'zorder': 3,
                              'markersize': marker_size,
                              'linewidth': line_width,
                              }
                else:
                    if only_show_highlights:
                        continue

                    if not codons_to_highlight:
                        alpha = 0.5
                    else:
                        alpha = 0.1

                    kwargs = {'color': 'black',
                              'alpha': alpha,
                              }

                handle, = sample_ax.plot(all_xs[sample, codon_id], all_ys[sample, codon_id], '.-', **kwargs)
                codon_to_handle[codon_id] = handle

            sample_ax.set_title(sample_to_label[sample])
            sample_ax.axhline(1, color='black')
            handles = [codon_to_handle[codon_id] for codon_id in codons_to_highlight]
            labels = [handle.get_label() for handle in handles]
            if len(codons_to_highlight) > 0:
                sample_ax.legend(handles, labels, framealpha=0.5, **legend_kwargs)

    for ax in axs.flatten():
        ax.set_xlim(min_x, max_x)
        if log_scale:
            ax.set_yscale('log', basey=2)

        if mark_active_sites:
            mark_active_sites_and_borders(ax, alpha=0.7)

        if mark_stalls:
            for p in range(10, max_x, 10):
                ax.axvline(p, ls=':', color='black', alpha=0.2)
    
        if force_ylims:
            ax.set_ylim(force_ylims)

        ax.set_xlabel('Offset (codons)')
        ax.set_ylabel('Mean relative enrichment')

        if flip:
            ax.invert_xaxis()
            flipped_labels = [str(int(-x)) for x in ax.get_xticks()]
            ax.set_xticklabels(flipped_labels)

    return fig

def mark_active_sites_and_borders(ax, alpha=0.4):
    tRNA_sites = [('A', 0, 'red'),
                  ('P', -1, 'blue'),
                  ('E', -2, 'green'),
                 ]

    read_borders = [('left', -5, 'black'),
                    ('right', 4, 'black'),
                   ]

    for site, position, color in tRNA_sites + read_borders: 
        ax.axvline(position, color=color, alpha=alpha)

def plot_dicodon_enrichments(names,
                             stratified_mean_enrichments_dict,
                             allowed_at_P,
                             allowed_at_A,
                             dicodons_to_highlight,
                             min_x=-30,
                             max_x=30,
                             log_scale=False,
                             force_ylims=None,
                             split_by_codon=False,
                             sample_to_label=None,
                             flip=False,
                            ):

    bmap = brewer2mpl.get_map('Set1', 'qualitative', 9)
    colors = bmap.mpl_colors[:5] + bmap.mpl_colors[6:]
    colors_iter = itertools.cycle(iter(colors))

    all_xs = {}
    all_ys = {}
    
    allowed_dicodons = list(itertools.product(allowed_at_P, allowed_at_A))

    for sample in names:
        for dicodon in allowed_dicodons:
            xs = np.arange(min_x, max_x + 1)
            ys = stratified_mean_enrichments_dict[sample]['dicodon', min_x:max_x + 1, dicodon]

            all_xs[sample, dicodon] = xs
            all_ys[sample, dicodon] = ys
            
    dicodon_to_color = {dicodon: colors_iter.next() for dicodon in dicodons_to_highlight}
    dicodon_to_handle = {}

    fig, axs = plt.subplots(len(names), 1,
                            figsize=(16, 12 * len(names)),
                            squeeze=False,
                           )

    for sample, sample_ax in zip(names, axs.flatten()):
        for dicodon in allowed_dicodons:
            if dicodon in dicodons_to_highlight:
                first_codon, second_codon = dicodon
                first_aa = codons.full_forward_table[first_codon]
                second_aa = codons.full_forward_table[second_codon]

                kwargs = {'color': dicodon_to_color[dicodon],
                          'alpha': 1,
                          'label': '{0} ({1}) - {2} ({3})'.format(first_codon, first_aa, second_codon, second_aa),
                          }
            else:
                if not dicodons_to_highlight:
                    alpha = 0.3
                else:
                    alpha = 0.1

                kwargs = {'color': 'black',
                          'alpha': alpha,
                          }

            handle, = sample_ax.plot(all_xs[sample, dicodon], all_ys[sample, dicodon], '.-', **kwargs)
            dicodon_to_handle[dicodon] = handle

        sample_ax.set_title(sample)
        handles = [dicodon_to_handle[dicodon] for dicodon in dicodons_to_highlight]
        labels = [handle.get_label() for handle in handles]
        if len(dicodons_to_highlight) > 0:
            sample_ax.legend(handles, labels, framealpha=0.5)

    for ax in axs.flatten():
        ax.set_xlim(min_x, max_x)
        if log_scale:
            ax.set_yscale('log', basey=2)

        mark_active_sites_and_borders(ax)

        for p in range(10, max_x, 10):
            ax.axvline(p, ls=':', color='black', alpha=0.2)
    
        if force_ylims:
            ax.set_ylim(force_ylims)

        ax.set_xlabel('Offset (codons)')
        ax.set_ylabel('Mean relative enrichment')
        
        if flip:
            ax.invert_xaxis()
            flipped_labels = [str(int(-x)) for x in ax.get_xticks()]
            ax.set_xticklabels(flipped_labels)

    return fig

def plot_codon_enrichments_all_amino_acids(enrichments, figure_file_name):
    with PdfPages(figure_file_name) as pdf:
        for amino_acid in codons.full_back_table:
            if amino_acid == '*':
                continue

            fig = plot_codon_enrichments(sorted(enrichments),
                                         enrichments,
                                         codons.full_back_table[amino_acid],
                                         min_x=-90,
                                         max_x=89,
                                         flip=True,
                                        )
            pdf.savefig(figure=fig, bbox_inches='tight')
            plt.close(fig)

def load_TEs(RPF_experiment, mRNA_experiment):
    def experiment_to_RPKMs(experiment):
        read_counts = experiment.read_file('read_counts')
        counts = {gene_name: read_counts[gene_name]['expression'][0] for gene_name in read_counts}
        total = sum(counts.values())
        RPKMs = {gene_name: max(0.1, (1.e9 / total) * counts[gene_name] / CDS_lengths[gene_name]) for gene_name in counts}
        return RPKMs

    transcripts, _ = RPF_experiment.get_CDSs()
    CDS_lengths = {t.name: t.CDS_length for t in transcripts}

    RPF_rpkms = experiment_to_RPKMs(RPF_experiment)
    mRNA_rpkms = experiment_to_RPKMs(mRNA_experiment)
          
    TEs = {gene_name: RPF_rpkms[gene_name] / mRNA_rpkms[gene_name] for gene_name in RPF_rpkms}
    return TEs

def compute_binned_means(around_lists, quantize_at, num_bins):
    num_before = around_lists['num_before']
    
    if quantize_at == 'TE':
        bins, binned, quantiles = split_into_TE_bins(around_lists, num_bins)
    else:
        bins, binned, quantiles = split_into_bins(around_lists, quantize_at, num_bins)
    
    offsets = np.arange(-40, 41)

    all_binned_means = np.asarray([np.mean(b, axis=0) for b in binned])
    binned_means = {offset: all_binned_means[:, num_before + offset] for offset in offsets}
    
    if quantize_at == 'TE':
        binned_means['TE'] = [np.mean(around_lists['TEs'][bins[i]]) for i in np.arange(num_bins)]
        
    return binned_means, quantiles

def plot_offset_relationships(around_lists,
                              quantize_at,
                              num_bins,
                              offsets_to_highlight=[0, -10, -20, -30, -40, 10],
                             ):
    fig, ax = plt.subplots(figsize=(16, 12))
    
    binned_means, quantiles = compute_binned_means(around_lists, quantize_at, num_bins)
    
    if quantize_at == 'TE':
        xs = binned_means['TE']
    else:
        xs = binned_means[quantize_at]

    for offset in offsets_to_highlight:
        ax.plot(xs, binned_means[offset], 'o-', label='{0:+}'.format(offset))
        
    for offset in binned_means:
        if offset in offsets_to_highlight or offset == 'TE':
            continue
            
        ax.plot(xs, binned_means[offset], '.-', color='black', alpha=0.1)
        
    for q in quantiles[1:-1]:
        ax.axvline(q, color='black', alpha=0.2)

    ax.set_xscale('log', basex=2)
    ax.set_yscale('log', basey=2)

    if quantize_at == 'TE':
        ax.set_xlabel('(arbitrary scaling of) average ribosome density per message') 
    else:
        ax.set_xlabel('Average enrichment in bin at His')
        
    ax.set_ylabel('Average enrichment in bin at offset')
    ax.legend(loc='upper right', framealpha=0.5)
    
def scatter_offset_relationships(around_lists,
                                 x_offset,
                                 y_offset,
                                 num_bins,
                                 log_scales=True,
                                 aspect_equal=True,
                                ):
    fig, ax = plt.subplots(figsize=(14, 14))
    
    binned_means, quantiles = compute_binned_means(around_lists, x_offset, num_bins)
    num_before = around_lists['num_before']
    ratios = {offset: around_lists['ratios'][:, num_before + offset] for offset in binned_means if offset != 'TE'}
    if x_offset == 'TE':
        ratios['TE'] = around_lists['TEs']
        
    smallest_nonzeros = {offset: min(r for r in rs if r > 0) for offset, rs in ratios.items()}
    
    def sanitize(offset):
        values = np.maximum(smallest_nonzeros[offset], ratios[offset])
        if log_scales:
            values = np.log2(values)
        return values
        
    sanitized = {offset: sanitize(offset) for offset in smallest_nonzeros}
    
    Sequencing.Visualize.enhanced_scatter(sanitized[x_offset],
                                          sanitized[y_offset],
                                          ax,
                                          do_fit=False,
                                          hists_height=0.2,
                                         )

    if x_offset == 'TE':
        ax.set_xlabel('(arbitrary scaling of) average ribosome density per message')
    else:
        ax.set_xlabel('log2(enrichment at {0})'.format(x_offset))
    
    ax.set_ylabel('log2(enrichment at {0})'.format(y_offset))
    
    x_means = binned_means[x_offset]
    y_means = binned_means[y_offset]
    if log_scales:
        x_means = np.log2(x_means)
        y_means = np.log2(y_means)
    
    if aspect_equal:
        ax.set_aspect(1.)
    
    ax.plot(x_means, y_means, 'o-', label='0', color='black')

    for q in quantiles[1:-1]:
        if log_scales:
            q = np.log2(q)
        ax.axvline(q, color='black', alpha=0.2)

def area_under_curve_additive(enrichments, CHX_name, no_CHX_name, x_min, x_max, ax):
    tAIs = load_tRNA_copy_numbers()
    tAI_values = [tAIs[codon] for codon in codons.non_stop_codons]
    areas = []

    for codon_id in codons.non_stop_codons:
        ys = enrichments[CHX_name]['codon', x_min:x_max + 1, codon_id]
        area = sum(ys - 1)
        areas.append(area)

    areas = np.array(areas)
    no_CHX_A = enrichments[no_CHX_name]['codon', 0, codons.non_stop_codons]
    no_CHX_P = enrichments[no_CHX_name]['codon', -1, codons.non_stop_codons]
    CHX_A = enrichments[CHX_name]['codon', 0, codons.non_stop_codons]
    CHX_P = enrichments[CHX_name]['codon', -1, codons.non_stop_codons]

    xs = areas
    ys = no_CHX_A + no_CHX_P - CHX_A - CHX_P

    Sequencing.Visualize.enhanced_scatter(xs, ys, ax,
                                          color_by_density=False,
                                          marker_size=10,
                                          #do_fit=False,
                                         )

    labels = ['{0} ({1})'.format(codon_id, codons.forward_table[codon_id]) for codon_id in codons.non_stop_codons]
    to_label = np.array([codon_id in ['CGA', 'CGG', 'CCG', 'CAC', 'CAT'] for codon_id in codons.non_stop_codons])
    
    #label_scatter_plot(ax, xs, ys, labels, to_label, vector='sideways', arrow_alpha=0.5)
    ax.set_ylabel('Net change in enrichment\nat A + P sites', size=14)
    ax.set_xlabel('Area under downstream wave\nwith CHX treatment', size=14)

    ax.axvline(0, color='black', alpha=0.5)
    ax.axhline(0, color='black', alpha=0.5)

def area_under_curve(stratified_mean_enrichments_dict, CHX_name, no_CHX_name, x_min, x_max):
    tAIs = load_tRNA_copy_numbers()
    tAI_values = [tAIs[codon] for codon in codons.non_stop_codons]
    areas = []

    for codon_id in codons.non_stop_codons:
        ys = stratified_mean_enrichments_dict[CHX_name]['codon', x_min:x_max + 1, codon_id]
        area = sum(ys - 1)
        areas.append(area)

    areas = np.array(areas)
    no_CHX_A = np.array([stratified_mean_enrichments_dict[no_CHX_name]['codon', 0, codon_id] for codon_id in codons.non_stop_codons])
    no_CHX_P = np.array([stratified_mean_enrichments_dict[no_CHX_name]['codon', -1, codon_id] for codon_id in codons.non_stop_codons])
    CHX_A = np.array([stratified_mean_enrichments_dict[CHX_name]['codon', 0, codon_id] for codon_id in codons.non_stop_codons])
    CHX_P = np.array([stratified_mean_enrichments_dict[CHX_name]['codon', -1, codon_id] for codon_id in codons.non_stop_codons])

    xs_list = [('areas', areas),
               ('areas', areas),
               ('areas', areas),
               ('areas', areas),
               ('areas', areas),
               ('CHX_A', CHX_A),
               ('no_CHX_A', no_CHX_A),
               ('no_CHX_A', no_CHX_A),
              ]

    ys_list = [('no_CHX_A', no_CHX_A),
               ('no_CHX_A - CHX_A', no_CHX_A - CHX_A),
               ('no_CHX_A + no_CHX_P - CHX_A - CHX_p', no_CHX_A + no_CHX_P - CHX_A - CHX_P),
               ('tAI_values', tAI_values),
               ('CHX_A', CHX_A),
               ('tAI_values', tAI_values),
               ('tAI_values', tAI_values),
               ('no_CHX_P', no_CHX_P),
              ]

    num_rows = int(np.ceil((len(xs_list) + 1) / 2.))
    
    fig, axs = plt.subplots(num_rows, 2, figsize=(16, 8 * num_rows))

    for ax, (x_label, xs), (y_label, ys) in zip(axs.flatten(), xs_list, ys_list):
        Sequencing.Visualize.enhanced_scatter(xs,
                                              ys,
                                              ax,
                                              color_by_density=False,
                                              marker_size=10,
                                             )

        labels = ['{0} ({1})'.format(codon_id, codons.forward_table[codon_id]) for codon_id in codons.non_stop_codons]
        to_label = np.array([codon_id in ['CGA', 'CGG', 'CCG', 'CAC', 'CAT'] for codon_id in codons.non_stop_codons])
        
        #label_scatter_plot(ax, xs, ys, labels, to_label, vector='sideways', arrow_alpha=0.5)
        ax.set_ylabel(y_label)
        ax.set_xlabel(x_label)
    
    A_ax = axs.flatten()[0]
    rho, p = scipy.stats.spearmanr(CHX_A, tAI_values)
    A_ax.annotate('A site occupancy - rho = {0:+0.2f} (p = {1:0.2e})'.format(rho, p),
                  xy=(0, 1),
                  xycoords='axes fraction',
                  xytext=(5, -15),
                  textcoords='offset points',
                  family='monospace',
                 )
    
    rho, p = scipy.stats.spearmanr(areas, tAI_values)
    A_ax.annotate('Area under curve - rho = {0:+0.2f} (p = {1:0.2e})'.format(rho, p),
                  xy=(0, 1),
                  xycoords='axes fraction',
                  xytext=(5, -30),
                  textcoords='offset points',
                  family='monospace',
                 )

    rho, p = scipy.stats.spearmanr(areas + CHX_A, tAI_values)
    A_ax.annotate('Area + A site    - rho = {0:+0.2f} (p = {1:0.2e})'.format(rho, p),
                  xy=(0, 1),
                  xycoords='axes fraction',
                  xytext=(5, -45),
                  textcoords='offset points',
                  family='monospace',
                 )

    alphas = np.linspace(-5, 5, 1000)
    rhos = []
    for alpha in alphas:
        rho, p = scipy.stats.pearsonr(areas + alpha * CHX_A, no_CHX_A)
        rhos.append(rho)
        
    ax = axs.flatten()[-1]
    ax.plot(alphas, rhos)

    max_alpha = alphas[np.argmax(rhos)]
    ax.axvline(max_alpha, color='black', alpha=0.5)
    ax.annotate(r'$\alpha$' + ' = {0:0.3f}'.format(max_alpha), 
                xy=(max_alpha, 0),
                xycoords=('data', 'axes fraction'),
                xytext=(5, 15),
                textcoords='offset points',
                size=20,
                )
    ax.set_xlim(min(alphas), max(alphas))
    
    return fig

def dicodon_area_under_curve(stratified_mean_enrichments_dict, CHX_name, no_CHX_name, x_min, x_max):
    tAIs = load_tRNA_copy_numbers()
    tAI_values = [tAIs[codon] for codon in codons.non_stop_codons]
    areas = []

    all_dicodons = list(itertools.product(codons.non_stop_codons, repeat=2))
    for dicodon in all_dicodons:
        ys = stratified_mean_enrichments_dict[CHX_name]['dicodon', x_min:x_max + 1, dicodon]
        area = sum(ys - 1)
        areas.append(area)

    areas = np.array(areas)
    no_CHX = np.array([stratified_mean_enrichments_dict[no_CHX_name]['dicodon', 0, dicodon] for dicodon in all_dicodons])
    CHX = np.array([stratified_mean_enrichments_dict[CHX_name]['dicodon', 0, dicodon] for dicodon in all_dicodons])

    xs_list = [('areas', areas),
               ('areas', areas),
               ('no_CHX', no_CHX),
              ]

    ys_list = [('no_CHX', no_CHX),
               ('no_CHX - CHX', no_CHX - CHX),
               ('CHX', CHX),
              ]

    num_rows = int(np.ceil((len(xs_list) + 1) / 2.))
    
    fig, axs = plt.subplots(num_rows, 2, figsize=(16, 8 * num_rows))

    for ax, (x_label, xs), (y_label, ys) in zip(axs.flatten(), xs_list, ys_list):
        Sequencing.Visualize.enhanced_scatter(xs,
                                              ys,
                                              ax,
                                              color_by_density=True,
                                              marker_size=10,
                                             )

        #labels = ['{0} ({1})'.format(codon_id, codons.forward_table[codon_id]) for codon_id in codons.non_stop_codons]
        #to_label = np.array([codon_id in ['CGA', 'CGG', 'CCG', 'CAC', 'CAT'] for codon_id in codons.non_stop_codons])
        
        #label_scatter_plot(ax, xs, ys, labels, to_label, vector='sideways', arrow_alpha=0.5)
        ax.set_ylabel(y_label)
        ax.set_xlabel(x_label)
    
    #A_ax = axs.flatten()[0]
    #rho, p = scipy.stats.spearmanr(CHX_A, tAI_values)
    #A_ax.annotate('A site occupancy - rho = {0:+0.2f} (p = {1:0.2e})'.format(rho, p),
    #              xy=(0, 1),
    #              xycoords='axes fraction',
    #              xytext=(5, -15),
    #              textcoords='offset points',
    #              family='monospace',
    #             )
    #
    #rho, p = scipy.stats.spearmanr(areas, tAI_values)
    #A_ax.annotate('Area under curve - rho = {0:+0.2f} (p = {1:0.2e})'.format(rho, p),
    #              xy=(0, 1),
    #              xycoords='axes fraction',
    #              xytext=(5, -30),
    #              textcoords='offset points',
    #              family='monospace',
    #             )

    #rho, p = scipy.stats.spearmanr(areas + CHX_A, tAI_values)
    #A_ax.annotate('Area + A site    - rho = {0:+0.2f} (p = {1:0.2e})'.format(rho, p),
    #              xy=(0, 1),
    #              xycoords='axes fraction',
    #              xytext=(5, -45),
    #              textcoords='offset points',
    #              family='monospace',
    #             )

    alphas = np.linspace(-5, 5, 1000)
    rs = []
    for alpha in alphas:
        r, p = scipy.stats.pearsonr(areas + alpha * CHX, no_CHX)
        rs.append(r)
        
    ax = axs.flatten()[-1]
    ax.plot(alphas, rs)

    max_alpha = alphas[np.argmax(rs)]
    ax.axvline(max_alpha, color='black', alpha=0.5)
    ax.annotate(r'$\alpha$' + ' = {0:0.3f}'.format(max_alpha), 
                xy=(max_alpha, 0),
                xycoords=('data', 'axes fraction'),
                xytext=(5, 15),
                textcoords='offset points',
                size=20,
                )
    ax.set_xlim(min(alphas), max(alphas))
    
    return fig

def offset_difference_correlation(enrichments, names,
                                  plot_lims=(-90, 89),
                                  enrichment_ylims=None,
                                  r_ylims=(-1, 1),
                                  size=1,
                                  min_p=-10,
                                  use_P_sites=False,
                                  show_A_site=False,
                                  variance_explained=False,
                                  p_value_panels=True,
                                  marker_size=6,
                                  line_width=1.5,
                                  text_size=14,
                                  withhold_results=False,
                                  annotate_maximum=True,
                                 ):
    if any('noCHX' in name for name in names):
        noCHX_name = [name for name in names if 'noCHX' in name][0]
        CHX_names = [name for name in names if 'noCHX' not in name]
    else:
        noCHX_name = names[0]
        CHX_names = names[1:]
    
    noCHX_As = enrichments[noCHX_name]['codon', 0, codons.non_stop_codons]
    noCHX_Ps = enrichments[noCHX_name]['codon', -1, codons.non_stop_codons]
    
    if p_value_panels:
        gs_kwargs = dict(hspace=0.07, wspace=0.1, height_ratios=[0.5, 1, 0.5])
        fig, axs = plt.subplots(3, len(CHX_names),
                                figsize=(size * 4 * len(CHX_names), size * 3 * 3),
                                gridspec_kw=gs_kwargs,
                                squeeze=False,
                               )
    else:
        gs_kwargs = dict(hspace=0.07, wspace=0.1, height_ratios=[0.3, 1])
        fig, axs = plt.subplots(2, len(CHX_names),
                                figsize=(size * 4 * len(CHX_names), size * 3 * 2),
                                gridspec_kw=gs_kwargs,
                                squeeze=False,
                               )
    
    for CHX_name, ax_col in zip(CHX_names, axs.T):
        if p_value_panels:
            enrichment_ax, r_ax, p_ax = ax_col
        else:
            enrichment_ax, r_ax = ax_col

        CHX_As = enrichments[CHX_name]['codon', 0, codons.non_stop_codons]
        CHX_Ps = enrichments[CHX_name]['codon', -1, codons.non_stop_codons]

        x_min, x_max = -90, 90
        xs = np.arange(x_min, x_max)

        ys_strings = [
            'noCHX_As - CHX_As',
            'noCHX_As + noCHX_Ps - CHX_As - CHX_Ps',
        ]
        
        labels = [
            'A site',
            'A site + P site',
        ]
        
        if not use_P_sites:
            ys_strings = ys_strings[:1]

        if use_P_sites and not show_A_site:
            ys_strings = ys_strings[1:]

        for ys_string, label in zip(ys_strings, labels):
            ys = eval(ys_string)
            x_rs, x_ps = np.array([scipy.stats.pearsonr(enrichments[CHX_name]['codon', x, codons.non_stop_codons], ys) for x in xs]).T
            
            if use_P_sites and show_A_site and ys_string == ys_strings[0]:
                alpha = 0.5
            else:
                alpha = 1.0
                
            if variance_explained:
                ys = x_rs**2
            else:
                ys = x_rs

            if not p_value_panels and show_A_site and ys_string == ys_strings[0]:
                color = 'green'
            else:
                color = 'blue'

            if not withhold_results:
                r_ax.plot(xs, ys, '.-',
                          alpha=alpha,
                          color=color,
                          label=label,
                          markersize=marker_size,
                          linewidth=line_width,
                         )

            if p_value_panels:
                p_ax.plot(xs, np.maximum(10**min_p, x_ps), '.-', alpha=alpha, color='green', label=label)

        plot_codon_enrichments([CHX_name], enrichments,
                               ['CGA'],
                               min_x=x_min,
                               max_x=x_max,
                               ax=enrichment_ax,
                               legend_kwargs={'loc': 'upper right'},
                               only_show_highlights=True,
                               mark_active_sites=False,
                               line_width=line_width,
                               marker_size=marker_size,
                              )

        if enrichment_ylims:
            enrichment_ax.set_ylim(*enrichment_ylims)
            enrichment_ax.set_yticks(np.arange(0, enrichment_ylims[1] + 0.1))

        r_ax.set_yticks([-1.0, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0])
        r_ax.yaxis.grid(linestyle='-', alpha=0.3)

        if variance_explained:
            r_ax.set_ylim(0, 1)
        else:
            r_ax.set_ylim(*r_ylims)
        
        r_ax.axhline(0, color='black')

        if p_value_panels:
            p_ax.set_yscale('log')
            p_ax.set_ylim(10**min_p, 1)
            p_ax.set_yticks([eval('1e{0}'.format(p)) for p in np.arange(0, min_p - 1, -4)])
            p_ax.yaxis.grid(linestyle='-', alpha=0.3)
            p_ax.set_ylabel('P-value of correlation', size=16)
            if len(ys_strings) > 1:
                p_ax.legend(framealpha=0.5, loc='lower right')
        
        label_kwargs = {'size': 16, 'family': 'serif'}

        ax_col[-1].set_xlabel('Offset (codons)', label_kwargs)

        for ax in ax_col:
            mark_active_sites_and_borders(ax, alpha=0.6)

            ax.set_xlim(*plot_lims)
            ax.invert_xaxis()
            flipped_labels = [str(int(-x)) for x in ax.get_xticks()]
            ax.set_xticklabels(flipped_labels)
            
            for label in ax.get_yticklabels() + ax.get_xticklabels():
                label.set_size(12)

        if len(ys_strings) > 1:
            if variance_explained:
                loc = 'upper right'
            else:
                loc = 'lower right'
            r_ax.legend(framealpha=0.5, loc=loc)
        else:
            r_ax.set_ylabel(ys_strings[0])

        if not use_P_sites:
            if variance_explained:
                r_ax.set_ylabel('$r^2$' + 'between 61 enrichments at offset\nand A-site changes', **label_kwargs)
            else:
                r_ax.set_ylabel('Correlation of 61 enrichments at offset\nwith A-site changes', **label_kwargs)
        else:
            if variance_explained:
                r_ax.set_ylabel('$r^2$' + 'between 61 enrichments at offset\nand active site changes', **label_kwargs)
            else:
                r_ax.set_ylabel('Correlation of 61 enrichments at offset\nwith active site changes', **label_kwargs)
        
        enrichment_ax.set_ylabel('Mean relative\nenrichment', **label_kwargs)
        
        # Annotate the maximum correlation
        i = np.argmax(ys)
        x = xs[i]
        y = ys[i]
        if variance_explained:
            label = r'$r^2 =\ {0:0.2f}$'.format(y)
        else:
            label = r'$r = {0:0.2f}$'.format(y)

        if not withhold_results and annotate_maximum:
            r_ax.annotate(label,
                          xy=(x, y),
                          xycoords='data',
                          xytext=(20, 0.93),
                          textcoords=('offset points', 'axes fraction'),
                          ha='left',
                          va='bottom',
                          arrowprops=dict(arrowstyle='->',
                                          color='black',
                                          linewidth=1.5,
                                          connectionstyle='arc3,rad=0.2',
                                          alpha=0.6,
                                          relpos=(0., 0.4),
                                         ),
                          size=text_size,
                         )
        
    for ax in axs[:, 1:].flatten():
        ax.set_ylabel('')
        ax.set_yticklabels([])
        
    for ax in axs[:-1, :].flatten():
        ax.set_xlabel('')
        ax.set_xticklabels([])

    if p_value_panels and any(x_ps <= 10**min_p):
        fig.canvas.draw()
        labels = [label.get_text() for label in axs[-1, 0].get_yticklabels()]    
        labels[-1] = labels[-1][:1] + '\leq' + labels[-1][1:]
        axs[-1, 0].set_yticklabels(labels)

    return fig

def offset_tAI_correlation(enrichments, CHX_names, plot_lims,
                           enrichment_ylims=None,
                           rho_ylims=(-0.6, 0.6),
                           size=6,
                           min_p=-7,
                           p_value_panels=True,
                           marker_size=6,
                           line_width=1.5,
                           text_size=14,
                           withhold_results=False,
                           p_offsets=None,
                           tRNA_value_source='tAI',
                          ):
    tAIs = load_tRNA_copy_numbers(tRNA_value_source)
    tAI_values = [tAIs[codon] for codon in codons.non_stop_codons]

    if p_offsets == None:
        p_offsets = [20]*len(CHX_names)
    
    #if p_value_panels:
    #    height_ratios = [0.3, 1, 0.3]
    #    gs_kwargs = dict(hspace=0.07, wspace=0.1, height_ratios=height_ratios, left=0, right=1)
    #    fig, axs = plt.subplots(3, len(CHX_names),
    #                            figsize=(0.7 * size * len(CHX_names), size * sum(height_ratios)),
    #                            gridspec_kw=gs_kwargs,
    #                            squeeze=False,
    #                           )
    #else:
    #    height_ratios = [0.3, 1]
    #    gs_kwargs = dict(hspace=0.07, wspace=0.1, height_ratios=height_ratios)
    #    fig, axs = plt.subplots(2, len(CHX_names),
    #                            figsize=(size * len(CHX_names), size * sum(height_ratios)),
    #                            gridspec_kw=gs_kwargs,
    #                            squeeze=False,
    #                           )

    num_exps = len(CHX_names)
    exps_per_row = 4
    full_rows, leftover = divmod(num_exps, exps_per_row)
    if leftover > 0:
        num_rows = full_rows + 1
    else:
        num_rows = full_rows
    height_ratios = [0.3, 1, 0.3]
    fig = plt.figure(figsize=(0.7 * size * exps_per_row, size * num_rows * sum(height_ratios)))

    columns = []

    high_level = matplotlib.gridspec.GridSpec(num_rows, 1, left=0, right=1, bottom=0, top=1)
    for row in high_level:
        low_level = matplotlib.gridspec.GridSpecFromSubplotSpec(3, exps_per_row, subplot_spec=row, hspace=0.07, wspace=0.1, height_ratios=height_ratios)
        row_axes = []
        for p in low_level:
            ax = plt.Subplot(fig, p)
            fig.add_subplot(ax)
            row_axes.append(ax)
            
        for column in np.array(row_axes).reshape(low_level.get_geometry()).T:
            columns.append(column)

    for CHX_name, ax_col, p_offset in itertools.izip_longest(CHX_names, columns, p_offsets):
        if p_value_panels:
            enrichment_ax, rho_ax, p_ax = ax_col
        else:
            enrichment_ax, rho_ax = ax_col

        if CHX_name == None:
            for ax in ax_col:
                fig.delaxes(ax)
            continue

        As = enrichments[CHX_name]['codon', 0, codons.non_stop_codons]
        
        x_rhos = []
        x_ps = []
        
        x_min, x_max = plot_lims
        xs = np.arange(x_min, x_max + 1)
            
        x_rhos, x_ps = np.array([scipy.stats.spearmanr(enrichments[CHX_name]['codon', x, codons.non_stop_codons], tAI_values) for x in xs]).T

        if withhold_results:
            for x, rho in zip(xs, x_rhos):
                if x == 0:
                    rho_ax.scatter([x], [rho],
                                   color='blueviolet',
                                   s=20,
                                  )
        else:
            rho_ax.plot(xs, x_rhos, '.-',
                        alpha=1,
                        color='blueviolet',
                        markersize=marker_size,
                        linewidth=line_width,
                       )

        label_kwargs = {'size': 16, 'family': 'serif'}

        if p_value_panels:
            if withhold_results:
                for x, p in zip(xs, x_ps):
                    if x == 0:
                        p_ax.scatter([x], [p],
                                     color='limegreen',
                                     s=20,
                                    )
            else:
                p_ax.plot(xs, np.maximum(10**min_p, x_ps), '.-',
                          alpha=1,
                          color='limegreen',
                          markersize=marker_size,
                          linewidth=line_width,
                         )
                for x, p in zip(xs, x_ps):
                    if p < 10**min_p:
                        p_ax.scatter([x], [10**min_p], s=20, color='limegreen', clip_on=False)
            
            p_ax.set_yscale('log')
            p_ax.set_ylim(ymax=1, ymin=10**min_p)
            p_ax.set_yticks([eval('1e{0}'.format(p)) for p in np.arange(0, min_p - 1, -2)])
            p_ax.set_ylabel('P-value of\nrank correlation', **label_kwargs)
            p_ax.yaxis.grid(linestyle='-', alpha=0.3)

        plot_codon_enrichments([CHX_name], enrichments,
                               ['CGA'],
                               min_x=x_min,
                               max_x=x_max,
                               ax=enrichment_ax,
                               only_show_highlights=True,
                               mark_active_sites=False,
                               line_width=line_width,
                               marker_size=marker_size,
                              )
        if enrichment_ylims:
            enrichment_ax.set_ylim(*enrichment_ylims)
            enrichment_ax.set_yticks(np.arange(0, enrichment_ylims[1] + 0.1))

        for ax in ax_col:
            mark_active_sites_and_borders(ax, alpha=0.6)

            ax.set_xlim(*plot_lims)
            ax.invert_xaxis()
            flipped_labels = [str(int(-x)) for x in ax.get_xticks()]
            ax.set_xticklabels(flipped_labels)

        rho_ax.set_ylim(*rho_ylims)
        
        rho_ax.set_ylabel('Rank correlation of\n61 enrichments with 1 / tAI', **label_kwargs)
        enrichment_ax.set_ylabel('Mean relative\nenrichment', **label_kwargs)
        
        ax_col[-1].set_xlabel('Offset (codons)', **label_kwargs)
        
        for ax in ax_col:
            for label in ax.get_yticklabels() + ax.get_xticklabels():
                label.set_size(12)
        
        rho_ax.axhline(0, color='black')
        enrichment_ax.set_title(CHX_name)
        
        rho_ax.yaxis.grid(linestyle='-', alpha=0.3)
        
        # Annotate the maximum correlation
        i = np.argmax(x_rhos)
        x = xs[i]
        rho = x_rhos[i]
        p = x_ps[i]
        if not withhold_results:
            rho_ax.annotate(r'$\rho = {0:0.2f}$'.format(rho),
                            xy=(x, rho),
                            xycoords='data',
                            xytext=(5, 1.05),
                            textcoords=('offset points', 'axes fraction'),
                            ha='left',
                            va='bottom',
                            arrowprops=dict(arrowstyle='->',
                                            color='black',
                                            linewidth=1.5,
                                            connectionstyle='arc3,rad=0.2',
                                            alpha=0.6,
                                            relpos=(0., 0.),
                                           ),
                            size=text_size,
                           )

            if p_value_panels:
                mantissa, exponent = '{0:0.1e}'.format(p).split('e')
                if p < 10**min_p:
                    p_string = r'p < 10^{{{0}}}'.format(min_p)
                    p = 10**min_p
                else:
                    p_string = r'p = {0} \times 10^{{{1}}}'.format(mantissa, int(exponent))
                p_ax.annotate(r'${0}$'.format(p_string),
                              xy=(x, p),
                              xycoords='data',
                              xytext=(p_offset, 0.05),
                              textcoords=('offset points', 'axes fraction'),
                              ha=('left' if p_offset > 0 else 'right'),
                              va='bottom',
                              arrowprops=dict(arrowstyle='->',
                                              color='black',
                                              linewidth=1.5,
                                              connectionstyle='arc3,rad=-0.2' if p_offset > 0 else 'arc3,rad=0.2',
                                              alpha=0.6,
                                              relpos=(0., 0.3) if p_offset > 0 else (1., 0.3),
                                             ),
                              size=text_size,
                             )

    for i, ax_col in enumerate(columns):
        if i % exps_per_row != 0:
            for ax in ax_col:
                ax.set_ylabel('')
                ax.set_yticklabels([])

        for ax in ax_col[:-1]:
            ax.set_xlabel('')
            ax.set_xticklabels([])

    ##for ax in axs[:, 1:].flatten():
    ##    if ax.legend_:
    ##        ax.legend_.remove()
    #
    #if p_value_panels and any(x_ps <= 10**min_p):
    #    fig.canvas.draw()
    #    labels = [label.get_text() for label in axs[-1, 0].get_yticklabels()]    
    #    labels[-1] = labels[-1][:1] + '\leq' + labels[-1][1:]
    #    axs[-1, 0].set_yticklabels(labels)

    return fig, columns

def correlation_heatmap(enrichments, names, labels=None, offset=0, ax=None, cmap=matplotlib.cm.RdBu_r):
    if ax == None:
        fig, ax = plt.subplots(figsize=(8, 8))
    else:
        fig = ax.figure

    if labels == None:
        labels = names

    all_values = [enrichments[name]['codon', offset, codons.non_stop_codons] for name in names]
    correlations = np.zeros((len(names), len(names)))
    for row, row_values in enumerate(all_values):
        for col, col_values in enumerate(all_values):
            correlations[row, col], _ = scipy.stats.pearsonr(row_values, col_values)

    im = ax.imshow(correlations,
                   interpolation='nearest',
                   cmap=cmap,
                   vmin=-1, vmax=1,
                  )
    ax.set_yticks(range(len(labels)))
    ax.set_xticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.tick_params(labeltop=True, labelbottom=False)
    ax.set_xticklabels(labels, rotation=45, ha='left')
    for t in ax.yaxis.get_ticklines() + ax.xaxis.get_ticklines():
        t.set_visible(False)

    p = ax.get_position()
    height = p.height / 2.
    start = (p.y1 + p.y0) / 2. - height / 2.
    colorbar_ax = fig.add_axes([p.x1 + 0.05 * p.width, start, p.width * 0.05, height])
    fig.colorbar(im, cax=colorbar_ax, ticks=[-1, 0, 1])
    
    colorbar_ax.annotate('r',
                         xy=(0.5, 1),
                         xycoords='axes fraction',
                         xytext=(0, 5),
                         textcoords='offset points',
                         size=16,
                         weight='bold',
                         ha='center',
                         va='bottom',
                         family='serif',
                        )

    return fig
