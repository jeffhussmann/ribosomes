import positions
import codons
import numpy as np
import itertools
from collections import defaultdict, Counter
import scipy.stats
import brewer2mpl
import matplotlib.pyplot as plt
import matplotlib.cm
import Sequencing.utilities
import Sequencing.Visualize
import itertools

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

def metacodon_around_pauses(codon_counts, allowed_at_pause, not_allowed_at_stall):
    np.seterr(all='raise')
    gene_names = codon_counts.keys()

    cds_slice = slice(('start_codon', 2), 'stop_codon')
    
    counts_around_list = []
    ratios_around_list = []
    codons_around_list = []
    nucleotides_around_list = []
    
    num_before = 90
    num_after = 90
    
    for gene_name in gene_names:
        counts = codon_counts[gene_name]['relaxed'][cds_slice]
        codons = codon_counts[gene_name]['identities'][cds_slice]

        median = np.median(counts)
        mean = np.mean(counts)
        
        if mean < 1:
            continue
        
        if len(counts) < num_before + num_after + 1:
            continue
            
        #denominators = means_of_rest(counts)
        denominators = np.mean(counts[num_before:-num_after])
            
        if denominators == 0:
            continue
            
        ratios = np.true_divide(counts, denominators)
            
        for offset in range(num_before, len(counts) - num_after):
            if codons[offset] in allowed_at_pause and codons[offset - 10] not in not_allowed_at_stall:
                around_slice = slice(offset - num_before, offset + num_after + 1)
                counts_around = counts[around_slice]
                ratios_around = ratios[around_slice]
                codons_around = codons[around_slice]
                nucleotides_around = np.array(list(''.join(codons_around)))
                
                counts_around_list.append(counts_around)
                ratios_around_list.append(ratios_around)
                codons_around_list.append(codons_around)
                nucleotides_around_list.append(nucleotides_around)
                    
    around_lists = {'counts': np.asarray(counts_around_list),
                    'ratios': np.asarray(ratios_around_list),
                    'codons': np.asarray(codons_around_list),
                    'nucleotides': np.asarray(nucleotides_around_list),
                    'num_before': num_before,
                    'num_after': num_after,
                   }
    
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
    
    ratios = around_lists['ratios'][:, num_before]
    
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
    num_quantiles = 10
    
    ratios_around = around_lists['ratios']

    quantiles = scipy.stats.mstats.mquantiles(ratios_around[:, 30 + quantize_at], prob=np.linspace(0, 1, num_quantiles + 1))
    first_nonzero = quantiles.nonzero()[0][0]
    quantiles = quantiles[first_nonzero - 1:]

    boundaries = zip(quantiles, quantiles[1:])
    bins = [(ratios_around[:, 30 + quantize_at] >= start) & (ratios_around[:, 30 + quantize_at] < end) for start, end in boundaries]
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
                                          draw_diagonal=True,
                                          hists_height=0.2,
                                         )

    if draw_labels:
        ax.set_xlabel('log2({0})'.format(x_name))
        ax.set_ylabel('log2({0})'.format(y_name))
        ax.set_title('Consistency of position-specific enrichments ({0} codons)'.format(set_name))

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
    
def label_scatter_plot(ax, xs, ys, labels, to_label, vector='orthogonal', initial_distance=50):
    def attempt_text(x, y, site, distance):
        if vector == 'orthogonal':
            x_offset = np.sign(x - y) * distance
            y_offset = -np.sign(x - y) * distance
        elif vector == 'radial':
            norm = np.linalg.norm([x, y])
            x_offset = x * distance / norm
            y_offset = y * distance / norm
            
        text = ax.annotate(site,
                           xy=(x, y),
                           xycoords=('data', 'data'),
                           xytext=(x_offset, y_offset),
                           textcoords='offset points',
                           ha='center',
                           size=10,
                           arrowprops={'arrowstyle': '->', 'alpha': 0.2},
                          )
        ax.figure.canvas.draw()
        return text, text.get_window_extent()

    ax.figure.canvas.draw()
    bboxes = [ax.xaxis.get_label().get_window_extent(),
              ax.yaxis.get_label().get_window_extent(),
             ]

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
            
def load_premal_elongation_times():
    expected_times = {}
    fn = '/home/jah/projects/translation_elongation/data/S.cer.code'
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
    for name in name_order:
        vals = [stratified_mean_enrichments_dict[name][codon_slice][codon] for codon in codons.non_stop_codons]
        ranks = np.array(vals).argsort().argsort()
        all_ranks.append(ranks)
    all_ranks = np.array(all_ranks).T

    codon_to_ranks = {codon: ranks for codon, ranks in zip(codons.non_stop_codons, all_ranks)}
    return codon_to_ranks

def plot_enrichments_across_conditions(stratified_mean_enrichments_dict,
                                       position,
                                       name_order,
                                       highlight_movement=True,
                                       by_rank=False,
                                       label_largest=0,
                                       log_scale=True,
                                       force_label=set(),
                                      ):
    fig, ax = plt.subplots(figsize=(16, 12))

    codon_slice = (position, position + 1, position + 2)

    codon_to_ranks = build_codon_to_ranks(stratified_mean_enrichments_dict, position, name_order)
    codon_to_color = assign_colors_by_values(codon_to_ranks)

    deltas = {codon_id: codon_to_ranks[codon_id][0] - codon_to_ranks[codon_id][-1] for codon_id in codons.non_stop_codons}
    
    xs = range(len(name_order))
        
    for codon_id in codons.non_stop_codons:
        if highlight_movement:
            alpha = min(1, abs(deltas[codon_id]) / float(50))
            width = abs(deltas[codon_id]) / float(30)
        else:
            alpha = 1
            width = 1
                  
        if by_rank:
            ys = codon_to_ranks[codon_id]
        else:
            ys = [stratified_mean_enrichments_dict[name][codon_slice][codon_id] for name in name_order]
        ax.plot(xs, ys, '.-', color=codon_to_color[codon_id], alpha=alpha, lw=width)
             
    sorted_deltas = sorted(deltas, key=lambda c: abs(deltas[c]), reverse=True)
    for codon_id in set(sorted_deltas[:label_largest]) | set(force_label):
        label = '{0} ({1})'.format(codon_id, codons.full_forward_table[codon_id])
        
        if by_rank:
            start_y = codon_to_ranks[codon_id][0]
            end_y = codon_to_ranks[codon_id][-1]
        else:
            start_y = stratified_mean_enrichments_dict[name_order[0]][codon_slice][codon_id]
            end_y = stratified_mean_enrichments_dict[name_order[-1]][codon_slice][codon_id]
            
        ax.annotate(label,
                    xy=(0, start_y),
                    xycoords=('data', 'data'),
                    xytext=(-15, 0),
                    textcoords='offset points',
                    ha='right',
                    va='center',
                    size=10,
                    arrowprops={'arrowstyle': '->', 'alpha': 0.3},
                   )
        
        ax.annotate(label,
                    xy=(len(name_order) - 1, end_y),
                    xycoords=('data', 'data'),
                    xytext=(15, 0),
                    textcoords='offset points',
                    ha='left',
                    va='center',
                    size=10,
                    arrowprops={'arrowstyle': '->', 'alpha': 0.3},
                   )
        
    if by_rank:
        ax.set_yticks([])
        ax.set_ylim(-0.5, 61.5)
    else:
        if log_scale:
            ax.set_yscale('log', basey=2)
            
    
    ax.set_xticks(xs)
    ax.set_xticklabels(name_order, rotation=45, ha='right')
    ax.set_xlim(min(xs) - 0.1, max(xs) + 0.1)

def plot_codon_enrichments(relevant_experiments,
                           stratified_mean_enrichments_dict,
                           aa_to_highlight,
                          ):
    fig, axs = plt.subplots(len(relevant_experiments), 1,
                            figsize=(16, 12 * len(relevant_experiments)),
                            squeeze=False,
                           )

    bmap = brewer2mpl.get_map('Set1', 'qualitative', 9)
    colors = bmap.mpl_colors[:5] + bmap.mpl_colors[6:]

    for experiment, ax in zip(relevant_experiments, axs.flatten()):
        sample = experiment.name
        colors_iter = cycle(colors)

        for codon_id in codons.non_stop_codons:
            amino_acid = codons.full_forward_table[codon_id]

            if amino_acid == aa_to_highlight:
                kwargs = {'color': colors_iter.next(),
                          'alpha': 1,
                          'label': '{0} ({1})'.format(codon_id, amino_acid),
                          }
            else:
                kwargs = {'color': 'black',
                          'alpha': 0.1,
                          }

            xs = []
            ys = []
            for i in range(-60, 31):
                xs.append(i)
                ys.append(stratified_mean_enrichments_dict[sample]['codon', i][codon_id])

            ax.plot(xs, ys, '.-', **kwargs)

        offset = -11
        xs = [offset]*len(codons.non_stop_codons)
        ys = [stratified_mean_enrichments_dict[sample]['codon', offset][codon_id] for codon_id in codons.non_stop_codons]
        labels = ['{0} ({1})'.format(codon_id, codons.full_forward_table[codon_id]) for codon_id in codons.non_stop_codons]

        to_label = sorted(range(len(codons.non_stop_codons)), key=ys.__getitem__)
        #pausing.label_scatter_plot(ax, xs, ys, labels, to_label[-10:], vector='orthogonal', initial_distance=100)
        
        ax.set_title(sample)
        ax.legend(framealpha=0.5)

    return fig
