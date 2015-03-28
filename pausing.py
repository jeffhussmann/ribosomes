import positions
import codons
import numpy as np
import itertools
from collections import defaultdict, Counter
import brewer2mpl
import matplotlib.pyplot as plt
import matplotlib.cm
from matplotlib.backends.backend_pdf import PdfPages
import Sequencing.utilities
import Sequencing.Visualize
import itertools
import scipy.stats
import os

igv_colors = Sequencing.Visualize.igv_colors.normalized_rgbs

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
    quantize_at = 0
    
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

def plot_nucleotide_enrichments(stratified_mean_enrichments, plot_A_site=True, min_x=-30, max_x=32, ax=None):
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
                y = stratified_mean_enrichments[x][base] / stratified_mean_enrichments['all']
            ys.append(y)
        ax.plot(xs, ys, '.-', color=igv_colors[base], label=base, markersize=7)

    ax.legend(framealpha=0.5)

    five_prime_edge = -15
    three_prime_edge = five_prime_edge + 28 - 1

    for x in range(min_x, max_x + 1):
        if x in [five_prime_edge, three_prime_edge]:
            alpha = 1
        elif x % 3 == 0:
            alpha = 0.2
        else:
            alpha = 0.1
        ax.axvline(x, color='black', alpha=alpha)

    ax.set_xlim(min_x, max_x)

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

    ax.set_ylabel('Relative mean enrichment at His codons with specific base at offset')
    ax.set_xlabel('Offset')
        
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
        xs = expecteds
        ys = actuals
    
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
                      ):
    def attempt_text(x, y, site, distance):
        if vector == 'orthogonal':
            x_offset = np.sign(x - y) * distance
            y_offset = -np.sign(x - y) * distance
        elif vector == 'radial':
            norm = np.linalg.norm([x, y])
            x_offset = x * distance / norm
            y_offset = y * distance / norm
        elif vector == 'sideways':
            x_offset = -distance
            y_offset
            
        text = ax.annotate(site,
                           xy=(x, y),
                           xycoords=('data', 'data'),
                           xytext=(x_offset, y_offset),
                           textcoords='offset points',
                           ha='center',
                           size=10,
                           arrowprops={'arrowstyle': '->',
                                       'alpha': arrow_alpha,
                                      },
                          )
        ax.figure.canvas.draw()
        return text, text.get_window_extent()

    ax.figure.canvas.draw()
    starting_labels = [ax.xaxis.get_label(), ax.yaxis.get_label()] + ax.get_yticklabels() + ax.get_xticklabels()
    bboxes = [label.get_window_extent() for label in starting_labels]

    tuples = itertools.izip(np.asarray(xs)[to_label],
                            np.asarray(ys)[to_label],
                            np.asarray(labels)[to_label],
                            )
    for x, y, label in tuples:
        distance = initial_distance
        text, bbox = attempt_text(x, y, label, distance)
        while any(bbox.fully_overlaps(other_bbox) for other_bbox in bboxes):
            text.remove()
            distance += 10
            text, bbox = attempt_text(x, y, label, distance)
        bboxes.append(bbox)

def label_enrichment_across_conditions_plot(ax, start_ys, end_ys, labels, num_conditions):
    def attempt_text(y, label, distance, side):
        if side == 'left':
            x = 0
            x_offset = -distance
            ha = 'right'
        elif side == 'right':
            x = num_conditions - 1
            x_offset = 1 + distance
            ha = 'left'
            
        text = ax.annotate(label,
                           xy=(x, y),
                           xycoords=('data', 'data'),
                           xytext=(x_offset, 0),
                           textcoords=('axes fraction', 'offset points'),
                           ha=ha,
                           va='center',
                           size=10,
                           arrowprops={'arrowstyle': '->', 'alpha': 0.5},
                          )
        ax.figure.canvas.draw()
        return text, text.get_window_extent()

    ax.figure.canvas.draw()

    starting_labels = [ax.xaxis.get_label(), ax.yaxis.get_label()] + ax.get_yticklabels()
    bboxes = [label.get_window_extent() for label in starting_labels]

    left_tuples = itertools.izip(start_ys, labels, itertools.repeat('left'))
    right_tuples = itertools.izip(end_ys, labels, itertools.repeat('right'))
    tuples = itertools.chain(left_tuples, right_tuples)

    for y, label, side in tuples:
        distance = 0.015
        text, bbox = attempt_text(y, label, distance, side)
        while any(bbox.fully_overlaps(other_bbox) for other_bbox in bboxes):
            text.remove()
            distance += 0.01
            text, bbox = attempt_text(y, label, distance, side)
        bboxes.append(bbox)
            
def load_premal_elongation_times():
    expected_times = {}
    fn = '{0}/projects/translation_elongation/data/S.cer.code'.format(os.environ['HOME'])
    for line in open(fn):
        codon, _, copy_number, wobble = line.strip().split()
        expected_times[codon] = 1. / (int(copy_number) * float(wobble))
    return expected_times

def assign_colors_by_values(codon_to_ranks):
    colors = get_jet_colors(len(codon_to_ranks))
    codon_to_color = {codon: colors[codon_to_ranks[codon][0]] for codon in codon_to_ranks}
    return codon_to_color

def build_codon_to_ranks(stratified_mean_enrichments_dict, position, name_order):
    codon_slice = (position, position + 1, position + 2)
    all_ranks = []
    all_values = []
    for name in name_order:
        values = [stratified_mean_enrichments_dict[name][codon_slice][codon] for codon in codons.non_stop_codons]
        ranks = np.array(values).argsort().argsort()
        all_values.append(values)
        all_ranks.append(ranks)

    all_values = np.array(all_values).T
    all_ranks = np.array(all_ranks).T

    codon_to_values = {codon: values for codon, values in zip(codons.non_stop_codons, all_values)}
    codon_to_ranks = {codon: ranks for codon, ranks in zip(codons.non_stop_codons, all_ranks)}

    return codon_to_ranks, codon_to_values

def plot_enrichments_across_conditions(stratified_mean_enrichments_dict,
                                       position,
                                       name_order,
                                       highlight_movement=True,
                                       force_highlight=set(),
                                       by_rank=False,
                                       label_rules=(0, 'rank', abs),
                                       log_scale=True,
                                       force_label=set(),
                                      ):
    fig, ax = plt.subplots(figsize=(16, 12))

    codon_slice = (position, position + 1, position + 2)

    codon_to_ranks, codon_to_values = build_codon_to_ranks(stratified_mean_enrichments_dict, position, name_order)
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
        
    for codon_id in codons.non_stop_codons:
        if highlight_movement:
            if rank_or_log == 'rank':
                alpha = min(1, abs(rank_deltas[codon_id]) / float(50))
                width = abs(rank_deltas[codon_id]) / float(30)
            elif rank_or_log == 'log':
                alpha = min(1, abs(log_deltas[codon_id]))
                width = abs(log_deltas[codon_id] / 2)**2
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

        ax.plot(xs, ys, '.-', color=codon_to_color[codon_id], alpha=alpha, lw=width)
             
    sorted_deltas = sorted(deltas, key=lambda c: transform(deltas[c]), reverse=True)

    if num_to_label > 0:
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

    for codon_id in set(sorted_deltas[:num_to_label]) | set(force_label):
        label = '{0} ({1})'.format(codon_id, codons.full_forward_table[codon_id])
        
        if by_rank:
            start_y = codon_to_ranks[codon_id][0]
            end_y = codon_to_ranks[codon_id][-1]
        else:
            start_y = stratified_mean_enrichments_dict[name_order[0]][codon_slice][codon_id]
            end_y = stratified_mean_enrichments_dict[name_order[-1]][codon_slice][codon_id]
            
        start_ys.append(start_y)
        end_ys.append(end_y)
        labels.append(label)

    if by_rank:
        ax.set_yticks([])
        ax.set_ylim(-0.5, 61.5)
    else:
        if log_scale:
            ax.set_yscale('log', basey=2)

        ax.tick_params(labelright=True)

        big_bold = matplotlib.font_manager.FontProperties(size=12, weight='bold')
        for label in ax.get_yticklabels():
            label.set_fontproperties(big_bold)

        ax.yaxis.grid(True, which='major', linestyle='-', alpha=0.3)
    
    label_enrichment_across_conditions_plot(ax, start_ys, end_ys, labels, len(name_order))

    ax.set_xticks(xs)
    ax.set_xticklabels(name_order, rotation=45, ha='right', size=12)
    ax.set_xlim(min(xs) - 0.1, max(xs) + 0.1)

def plot_correlations_across_conditions(stratified_mean_enrichments_dict,
                                        name_order,
                                        codon_position=0,
                                       ):
    rhos = []
    ps = []

    tAIs = load_premal_elongation_times()
    tAI_values = [tAIs[codon] for codon in codons.non_stop_codons]

    codon_positions = (3 * codon_position,
                       3 * codon_position + 1,
                       3 * codon_position + 2,
                      )

    for name in name_order:
        e_values = [stratified_mean_enrichments_dict[name][codon_positions][codon] for codon in codons.non_stop_codons]
        rho, p = scipy.stats.spearmanr(tAI_values, e_values)
        rhos.append(rho)
        ps.append(p)
        print '{0:<30}\t{1:10.4f}\t{2:0.4e}'.format(name, rho, p)
        
    fig, ax = plt.subplots(figsize=(16, 12))
    xs = np.arange(len(rhos))
    ax.plot(xs, rhos, 'o-', label=r'$\rho$')
    ax.set_ylabel(r'Spearman $\rho$', color='blue', size=16)

    p_ax = ax.twinx()
    p_ax.plot(xs, ps, 'o-', color='red', alpha=0.5, label='p value')
    p_ax.set_ylim(1e-10, 1)
    p_ax.set_yscale('log')
    p_ax.set_ylabel('p value', color='red', size=16)

    ax.legend(loc='upper left', framealpha=0.5)
    p_ax.legend(loc='upper right', framealpha=0.5)
    
    ax.set_xticks(xs)
    ax.set_xticklabels(name_order, rotation=45, ha='right', size=12)
    ax.set_xlim(min(xs) - 0.1, max(xs) + 0.1)

    ax.axhline(0, color='black', alpha=1)

    ax.set_title('Rank correlation of mean occupancy based on codon id at offset {0} from A-site with estimated tRNA scarcity'.format(codon_position))

    ax.set_ylim(-0.5, 1)

    p_ax.grid(alpha=0.7)
    
    big_bold = matplotlib.font_manager.FontProperties(size=12, weight='bold')
    for label in ax.get_yticklabels() + p_ax.get_yticklabels():
        label.set_fontproperties(big_bold)

    return fig
                                       
def plot_codon_enrichments(names,
                           stratified_mean_enrichments_dict,
                           codons_to_highlight,
                           min_x=-30,
                           max_x=30,
                           log_scale=False,
                           force_ylims=None,
                           split_by_codon=False,
                          ):

    bmap = brewer2mpl.get_map('Set1', 'qualitative', 9)
    colors = bmap.mpl_colors[:5] + bmap.mpl_colors[6:]
    colors_iter = itertools.cycle(iter(colors))

    all_xs = {}
    all_ys = {}

    for sample in names:
        for codon_id in codons.non_stop_codons:
            xs = []
            ys = []
            for i in range(min_x, max_x + 1):
                xs.append(i)
                codon_positions = (i * 3, i * 3 + 1, i * 3 + 2)
                ys.append(stratified_mean_enrichments_dict[sample][codon_positions][codon_id])

            all_xs[sample, codon_id] = xs
            all_ys[sample, codon_id] = ys
            

    if split_by_codon:
        sample_to_color = {name: colors_iter.next() for name in names}
        fig, axs = plt.subplots(len(codons_to_highlight), 1,
                                figsize=(16, 12 * len(codons_to_highlight)),
                                squeeze=False,
                               )
        for codon_id, codon_ax in zip(codons_to_highlight, axs.flatten()):
            amino_acid = codons.full_forward_table[codon_id]
            for sample in names:
                kwargs = {'color': sample_to_color[sample],
                          'alpha': 1,
                          'label': sample,
                         }
                codon_ax.plot(all_xs[sample, codon_id], all_ys[sample, codon_id], '.-', **kwargs)

            codon_ax.set_title('{0} ({1})'.format(codon_id, amino_acid))
            codon_ax.legend(loc='upper left', framealpha=0.5)
            codon_ax.axhline(1, color='black')

    else:
        codon_to_color = {codon_id: colors_iter.next() for codon_id in codons_to_highlight}
        codon_to_handle = {}

        fig, axs = plt.subplots(len(names), 1,
                                figsize=(16, 12 * len(names)),
                                squeeze=False,
                               )
        for sample, sample_ax in zip(names, axs.flatten()):
            for codon_id in codons.non_stop_codons:
                if codon_id in codons_to_highlight:
                    amino_acid = codons.full_forward_table[codon_id]

                    kwargs = {'color': codon_to_color[codon_id],
                              'alpha': 1,
                              'label': '{0} ({1})'.format(codon_id, amino_acid),
                              }
                else:
                    kwargs = {'color': 'black',
                              'alpha': 0.1,
                              }

                handle, = sample_ax.plot(all_xs[sample, codon_id], all_ys[sample, codon_id], '.-', **kwargs)
                codon_to_handle[codon_id] = handle

        #offset = -11
        #codon_positions = (offset * 3, offset * 3 + 1, offset * 3 + 2)
        #xs = [offset]*len(codons.non_stop_codons)
        #ys = [stratified_mean_enrichments_dict[sample][codon_positions][codon_id] for codon_id in codons.non_stop_codons]
        #labels = ['{0} ({1})'.format(codon_id, codons.full_forward_table[codon_id]) for codon_id in codons.non_stop_codons]

        #to_label = sorted(range(len(codons.non_stop_codons)), key=ys.__getitem__)
        #pausing.label_scatter_plot(ax, xs, ys, labels, to_label[-10:], vector='orthogonal', initial_distance=100)
        
            sample_ax.set_title(sample)
            handles = [codon_to_handle[codon_id] for codon_id in codons_to_highlight]
            labels = [handle.get_label() for handle in handles]
            if len(codons_to_highlight) > 0:
                sample_ax.legend(handles, labels, framealpha=0.5)

    for ax in axs.flatten():
        ax.set_xlim(min_x, max_x)
        if log_scale:
            ax.set_yscale('log', basey=2)

        tRNA_sites = [('A', 0, 'red'),
                      ('P', -1, 'blue'),
                      ('E', -2, 'green'),
                     ]

        read_borders = [('left', -5, 'black'),
                        ('right', 4, 'black'),
                       ]

        for site, position, color in tRNA_sites + read_borders: 
            ax.axvline(position, color=color, alpha=0.2)

        for p in range(10, max_x, 10):
            ax.axvline(p, ls=':', color='black', alpha=0.2)
    
        if force_ylims:
            ax.set_ylim(force_ylims)

    return fig

def plot_codon_enrichments_all_amino_acids(relevant_experiments, figure_file_name):
    stratified_mean_enrichments_dict = {exp.name: exp.read_file('stratified_mean_enrichments') for exp in relevant_experiments}
    
    with PdfPages(figure_file_name) as pdf:
        for amino_acid in codons.full_back_table:
            if amino_acid == '*':
                continue

            fig = plot_codon_enrichments(relevant_experiments,
                                         stratified_mean_enrichments_dict,
                                         amino_acid,
                                         min_x=-60,
                                         max_x=60,
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
