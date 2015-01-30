import matplotlib
matplotlib.use('Agg', warn=False)
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import positions
import numpy as np
import brewer2mpl
import Sequencing.Visualize
import Sequencing.utilities as utilities
import Sequencing.genomes as genomes
import gtf
import operator
import Circles.variants as variants
import codons
import scipy.stats
    
bmap = brewer2mpl.get_map('Set1', 'qualitative', 9)
colors = bmap.mpl_colors[:5] + bmap.mpl_colors[6:] + ['black']
colors = colors + colors

igv_colors = Sequencing.Visualize.igv_colors.normalized_rgbs

def smoothed(position_counts, window_size):
    smoothed_array = positions.PositionCounts(position_counts.landmarks,
                                              position_counts.left_buffer,
                                              position_counts.right_buffer,
                                              data=position_counts.data,
                                             )
    for i in range(window_size):
        smoothed_array[i] = position_counts[:i + 1].sum() / float(i + 1)
    for i in range(window_size, smoothed_array.extent_length - window_size):
        smoothed_array[i] = position_counts[i - window_size:i + window_size + 1].sum() / float(2 * window_size + 1)
    for i in range(smoothed_array.extent_length - window_size, smoothed_array.extent_length):
        smoothed_array[i] = position_counts[i:].sum() / float(smoothed_array.extent_length - i)
    return smoothed_array

def plot_metagene_positions(metagene_positions,
                            figure_fn,
                            lengths,
                            title=None,
                            by_base=False,
                           ):

    landmarks = [('start', 'end'),
                 ('start_codon', 'stop_codon'),
                ]

    pdf = matplotlib.backends.backend_pdf.PdfPages(figure_fn)

    for start_landmark, end_landmark in landmarks:
        fig, axs = plt.subplots(2, 2, figsize=(24, 16))

        for zoomed_out, (start_ax, end_ax) in zip([False, True], axs):
            if zoomed_out:
                start_xs = np.arange(-190, 290)
                end_xs = np.arange(-290, 190)
                tick_interval = 30
            else:
                start_xs = np.arange(-21, 19)
                end_xs = np.arange(-33, 7)
                tick_interval = 3
        
            if by_base:
                for base in 'TCAG':
                    #start_key = '{0}_{1}'.format(start_landmark, base)
                    #end_key = '{0}_{1}'.format(end_landmark, base)
                    #start_ys = metagene_positions[start_key]['all'][start_landmark, start_xs]
                    #end_ys = metagene_positions[end_key]['all'][end_landmark, end_xs]

                    #kwargs = {'label': base,
                    #          'color': igv_colors[base],
                    #          'marker': '.',
                    #         }
                    #
                    #start_ax.plot(start_xs, start_ys, **kwargs)
                    #end_ax.plot(end_xs, end_ys, **kwargs)
                    
                    start_key = '{0}_{1}_uniform'.format(start_landmark, base)
                    end_key = '{0}_{1}_uniform'.format(end_landmark, base)
                    start_ys = metagene_positions[start_key]['all'][start_landmark, start_xs]
                    end_ys = metagene_positions[end_key]['all'][end_landmark, end_xs]

                    kwargs = {'label': '{0}_uniform'.format(base),
                              'color': igv_colors[base],
                              'marker': '.',
                              'alpha': 0.3,
                             }
                    
                    start_ax.plot(start_xs, start_ys, **kwargs)
                    end_ax.plot(end_xs, end_ys, **kwargs)
            else:
                start_counts = metagene_positions[start_landmark]
                end_counts = metagene_positions[end_landmark]
                recorded_lengths = sorted(set(lengths) & set(start_counts.keys()))

                for length in recorded_lengths:
                    start_ys = start_counts[length][start_landmark, start_xs]
                    end_ys = end_counts[length][end_landmark, end_xs]
                    
                    kwargs = {'label': length,
                              'marker': '.',
                             }
                    if length == 'all':
                        kwargs['color'] = 'black'
                        kwargs['alpha'] = 0.7

                    start_ax.plot(start_xs, start_ys, **kwargs)
                    end_ax.plot(end_xs, end_ys, **kwargs)

            start_ax.set_xlim(min(start_xs), max(start_xs))
            start_xticks = [x for x in start_xs if x % tick_interval == 0]
            start_ax.set_xticks(start_xticks)
            for x in start_xticks:
                start_ax.axvline(x, color='black', alpha=0.1)
            start_ax.set_xlabel('Position of read relative to {0}'.format(start_landmark))
            start_ax.set_ylabel('Number of uniquely mapped reads of specified length')

            Sequencing.Visualize.add_commas_to_yticks(start_ax)
            
            end_ax.set_xlim(min(end_xs), max(end_xs))
            end_xticks = [x for x in end_xs if x % tick_interval == 0]
            end_ax.set_xticks(end_xticks)
            for x in end_xticks:
                end_ax.axvline(x, color='black', alpha=0.1)
            end_ax.set_xlabel('Position of read relative to {0}'.format(end_landmark))
            end_ax.set_ylabel('Number of uniquely mapped reads of specified length')
            
            Sequencing.Visualize.add_commas_to_yticks(end_ax)

            start_ax.legend(loc='upper right', framealpha=0.5)
            end_ax.legend(loc='upper right', framealpha=0.5)

            ymax = max(start_ax.get_ylim()[1], end_ax.get_ylim()[1])
            start_ax.set_ylim(0, ymax)
            end_ax.set_ylim(0, ymax)

        if title:
            fig.suptitle(title)
        
        pdf.savefig(fig)
        plt.close(fig)

    pdf.close()

def plot_just_ends(from_ends, figure_fn, title):
    fig, (linear_ax, log_ax) = plt.subplots(2, 1, figsize=(24, 16))

    end_xs = np.arange(-63, 97)
    tick_interval = 3
    
    lengths = range(27, 32)

    recorded_lengths = sorted(set(lengths) & set(from_ends.keys()))

    for length in recorded_lengths:
        linear_ax.plot(end_xs, from_ends[length]['stop_codon', end_xs], '.-', label=length)
        log_ax.plot(end_xs, from_ends[length]['stop_codon', end_xs], '.-', label=length)

    for ax in [linear_ax, log_ax]:
        ax.set_xlim(min(end_xs), max(end_xs))
        end_xticks = [x for x in end_xs if x % tick_interval == 0]
        ax.set_xticks(end_xticks)
        for x in end_xticks:
            ax.axvline(x, color='black', alpha=0.1)
        ax.set_xlabel('Position of read relative to stop codon')
        ax.set_ylabel('Number of uniquely mapped reads of specified length')
        
        ax.set_ylim(ymin=1)
        ax.legend(loc='upper right', framealpha=0.5)
        
    log_ax.set_yscale('log')
    Sequencing.Visualize.add_commas_to_yticks(linear_ax)

    fig.suptitle(title)

    fig.savefig(figure_fn)
    plt.close(fig)

def plot_metagene_positions_all_lengths(from_starts, from_ends, figure_fn, zoomed_out=False, title=None):
    lengths = sorted(from_starts['A'].keys())
    
    fig, axs = plt.subplots(len(lengths), 2, figsize=(24, 8 * len(lengths)))

    if zoomed_out:
        start_xs = np.arange(-190, 190)
        end_xs = np.arange(-190, 190)
        tick_interval = 30
    else:
        start_xs = np.arange(3, 91)
        end_xs = np.arange(-91, 0)
        tick_interval = 3
    
    for length, (start_ax, end_ax) in zip(lengths, axs):
        for base in 'TCAG':
            color = igv_colors[base]
            start_ys = from_starts[base][length]['start_codon', start_xs]
            end_ys = from_ends[base][length]['stop_codon', end_xs]
            start_ax.plot(start_xs, start_ys, '.-', label=base, color=color)
            end_ax.plot(end_xs, end_ys, '.-', label=base, color=color)

        start_ax.set_xlim(min(start_xs), max(start_xs))
        start_xticks = [x for x in start_xs if x % tick_interval == 0]
        start_ax.set_xticks(start_xticks)
        for x in start_xticks:
            start_ax.axvline(x, color='black', alpha=0.1)
        start_ax.set_xlabel('Position of read relative to start of CDS')
        start_ax.set_ylabel('Number of uniquely mapped reads of specified length')
        
        end_ax.set_xlim(min(end_xs), max(end_xs))
        end_xticks = [x for x in end_xs if x % tick_interval == 0]
        end_ax.set_xticks(end_xticks)
        for x in end_xticks:
            end_ax.axvline(x, color='black', alpha=0.1)
        end_ax.set_xlabel('Position of read relative to stop codon')
        end_ax.set_ylabel('Number of uniquely mapped reads of specified length')

        start_ax.legend(loc='upper right', framealpha=0.5)
        end_ax.legend(loc='upper right', framealpha=0.5)

        start_ax.set_title('Length {0}'.format(length))
        end_ax.set_title('Length {0}'.format(length))

        ymax = max(start_ax.get_ylim()[1], end_ax.get_ylim()[1])
        start_ax.set_ylim(0, ymax)
        end_ax.set_ylim(0, ymax)

    if title:
        fig.suptitle(title)
    
    fig.savefig(figure_fn, bbox_inches='tight')
    plt.close(fig)

def plot_metacodon_positions(position_counts, figure_fn, key='feature'):
    fig, (long_ax, short_ax) = plt.subplots(2, 1, figsize=(12, 16))

    xs = np.arange(-30, 30)
    
    short_lengths = range(20, 24)
    long_lengths = range(27, 32)

    for length in long_lengths:
        long_ax.plot(xs, position_counts[length][key, xs], '.-', label=length)
    for length in short_lengths:
        short_ax.plot(xs, position_counts[length][key, xs], '.-', label=length)

    mod_3 = [x for x in xs if x % 3 == 0]
    for ax in (long_ax, short_ax):
        ax.set_xlim(min(xs), max(xs))
        ax.set_xticks(mod_3)
        for x in mod_3:
            ax.axvline(x, color='black', alpha=0.1)
        ax.set_xlabel('Position of read relative to codon')
        ax.set_ylabel('Number of uniquely mapped reads of specified length')
        
        ax.legend(loc='upper right', framealpha=0.5)

    fig.savefig(figure_fn)
    plt.close(fig)

def plot_metagene_positions_heatmap(metagene_positions,
                                    figure_fn,
                                    zoomed_out=False,
                                    normalize_to_max_in=None,
                                   ):
    landmarks = [('start', 'end'),
                 ('start_codon', 'stop_codon'),
                ]

    pdf = matplotlib.backends.backend_pdf.PdfPages(figure_fn)

    if zoomed_out:
        start_xs = np.arange(-90, 60)
        end_xs = np.arange(-60, 90)
    else:
        start_xs = np.arange(-21, 31)
        end_xs = np.arange(-33, 19)

    if len(start_xs) != len(end_xs):
        raise ValueError(len(start_xs), len(end_xs))
    
    for start_landmark, end_landmark in landmarks:
        from_starts = metagene_positions[start_landmark]
        from_ends = metagene_positions[end_landmark]

        lengths = sorted([k for k in from_starts if isinstance(k, int)], reverse=True)
        from_starts_5_array = [from_starts[l][start_landmark, start_xs] for l in lengths]
        from_ends_5_array = [from_ends[l][end_landmark, end_xs] for l in lengths]
        from_starts_3_array = [positions.convert_to_three_prime(from_starts[l], l)[start_landmark, start_xs] for l in lengths]
        from_ends_3_array = [positions.convert_to_three_prime(from_ends[l], l)[end_landmark, end_xs] for l in lengths]

        from_starts_5_array = np.asarray(from_starts_5_array)
        from_starts_3_array = np.asarray(from_starts_3_array)
        from_ends_5_array = np.asarray(from_ends_5_array)
        from_ends_3_array = np.asarray(from_ends_3_array)

        if normalize_to_max_in:
            allowed_lengths, which_end, allowed_positions = normalize_to_max_in
            if which_end == 'start':
                counts = from_starts
                landmark = start_landmark
            elif which_end == 'end':
                counts = from_ends
                landmark = end_landmark
            allowed_counts = [counts[l][landmark, allowed_positions] for l in allowed_lengths]
            overall_max = np.asarray(allowed_counts).max()
        else:
            overall_max = max(from_starts_5_array.max(),
                              from_ends_5_array.max(),
                              from_starts_3_array.max(),
                              from_ends_3_array.max(),
                             )

        blues_cdict = {'red': ((0.0, 1.0, 1.0),
                               (1.0, 0.0, 0.0)),
                       'green': ((0.0, 1.0, 1.0),
                                 (1.0, 0.0, 0.0)),
                       'blue': ((0.0, 1.0, 1.0),
                                (1.0, 1.0, 1.0)),
                      }

        blues = matplotlib.colors.LinearSegmentedColormap('blues', blues_cdict, 10000)
        blues.set_over('red')
        
        reds_cdict = {'red': ((0.0, 1.0, 1.0),
                              (1.0, 1.0, 1.0)),
                      'green': ((0.0, 1.0, 1.0),
                                (1.0, 0.0, 0.0)),
                      'blue': ((0.0, 1.0, 1.0),
                               (1.0, 0.0, 0.0)),
                     }

        reds = matplotlib.colors.LinearSegmentedColormap('blues', reds_cdict, 10000)
        reds.set_over('blue')

        kwargs_five = {'interpolation': 'nearest',
                       'cmap': blues,
                       'vmin': 0,
                       'vmax': overall_max,
                      }
        kwargs_three = {'interpolation': 'nearest',
                        'cmap': reds,
                        'vmin': 0,
                        'vmax': overall_max,
                       }
        
        width = len(start_xs) * 2. / 3.
        height = len(lengths) * 2 / 3.
        fig, ((start_5_ax, end_5_ax), (start_3_ax, end_3_ax)) = plt.subplots(2, 2, figsize=(width, height))

        start_5_ax.imshow(from_starts_5_array, **kwargs_five)
        end_5_ax.imshow(from_ends_5_array, **kwargs_five)
        start_3_ax.imshow(from_starts_3_array, **kwargs_three)
        end_3_ax.imshow(from_ends_3_array, **kwargs_three)

        index_to_x = {i: x for i, x in enumerate(start_xs)}
        if zoomed_out:
            modulus = 30
        else:
            modulus = 3
        xticks = [i for i in range(len(start_xs)) if index_to_x[i] % modulus == 0]
        xticklabels = [str(index_to_x[i]) for i in xticks]

        for ax in (start_5_ax, start_3_ax):
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticklabels)

        start_5_ax.set_xlabel('Position of 5\' end of read relative to {0}'.format(start_landmark))
        start_3_ax.set_xlabel('Position of 3\' end of read relative to {0}'.format(start_landmark))
        
        index_to_x = {i: x for i, x in enumerate(end_xs)}
        xticks = [i for i in range(len(end_xs)) if index_to_x[i] % modulus == 0]
        xticklabels = [str(index_to_x[i]) for i in xticks]

        for ax in (end_5_ax, end_3_ax):
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticklabels)

        end_5_ax.set_xlabel('Position of 5\' end of read relative to {0}'.format(end_landmark))
        end_3_ax.set_xlabel('Position of 3\' end of read relative to {0}'.format(end_landmark))
        
        y_ticks = np.arange(len(lengths))
        y_tick_labels = [str(l) for l in lengths]
        
        for ax in (start_5_ax, start_3_ax, end_5_ax, end_3_ax):
            ax.set_yticks(y_ticks)
            ax.set_yticklabels(y_tick_labels)
        
        fig.set_size_inches((len(start_xs) * 2. / 3., len(lengths) * 2 / 3.))
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
    pdf.close()

@Sequencing.Visualize.optional_ax
def plot_joint_positions_scatter(joint_position_counts, transcript, name, ax=None):
    xs = []
    ys = []
    cs = []
    ss = []

    sorted_counts = sorted(joint_position_counts.items(), key=lambda x: x[1])
    for (x, y), c in sorted_counts:
        xs.append(x)
        ys.append(y)
        cs.append(c)
        #s = 3 * max(1, np.log10(c))
        s = 3 * max(1, c / 10.)
        ss.append(s)

    transcript.build_coordinate_maps()
    kwargs = {'color': 'black',
              'alpha': 0.5,
             }
    ax.axvline(0, **kwargs)
    ax.axvline(transcript.CDS_length, **kwargs)
    ax.axhline(0, **kwargs)
    ax.axhline(-transcript.CDS_length, **kwargs)
    x_lims = [-500, transcript.CDS_length + 500]
    y_lims = [-500 - transcript.CDS_length, 500]
    ax.plot(x_lims, y_lims, **kwargs)
    ax.scatter(xs, ys, ss, c=cs, linewidth=(0.,))
    ax.set_xlim(x_lims)
    ax.set_ylim(y_lims)
    ax.set_aspect(1.)
    ax.set_title(name)
    ax.set_xlabel('5\' end of transcript relative to start codon')
    ax.set_ylabel('3\' end of transcript relative to stop codon')

def plot_frames(from_starts, figure_fn):
    ''' Plots bar graphs of the distribution of reading frame that the 5' end of
        reads fell, stratified by read length.
    '''
    def bar_plot(values, name, tick_labels, ax):
        N = len(values)

        starts = np.arange(N)
        gap = 1 # in multiples of bar width

        width = 1. / (1 + gap)

        ax.bar(starts, values, width, label=name)

        ax.set_xlim(xmin=-width / 2, xmax=N - width / 2)

        ax.set_xticks(starts + (width * 1 / 2))
        ax.set_xticklabels(tick_labels)
        ax.xaxis.set_ticks_position('none')
        ax.set_ylim(0, 1)

    relevant_lengths = sorted(from_starts.keys())

    nonzero_lengths = []
    for length in relevant_lengths:
        frame_counts = np.zeros(3, int)
        for p in range(-15, from_starts[length].CDS_length - 15 + 3):
            frame_counts[p % 3] += from_starts[length]['start_codon', p]

        if frame_counts.sum() != 0:
            nonzero_lengths.append((length, frame_counts))

    fig, axs = plt.subplots(len(nonzero_lengths), 1, figsize=(6, 3 * len(nonzero_lengths)))
    for (length, frame_counts), ax in zip(nonzero_lengths, axs):
        frame_fractions = np.true_divide(frame_counts, frame_counts.sum())

        tick_labels = ['0', '1', '2']

        bar_plot(frame_fractions, str(length), tick_labels, ax)
        label = 'Length {0}\n{1:,} total reads'.format(length, frame_counts.sum())
        ax.text(0.99, 0.98, label,
                horizontalalignment='right',
                verticalalignment='top',
                transform=ax.transAxes,
               )
        
    axs[-1].set_xlabel('Frame')
    axs[len(axs) // 2].set_ylabel('Fraction of uniquely mapped reads of specified length')

    fig.savefig(figure_fn, bbox_inches='tight')

def plot_averaged_nucleotide_densities(data_sets,
                                       figure_fn,
                                       show_start=False,
                                       past_edge=10,
                                       smooth=False,
                                       plot_up_to=1000,
                                      ):
    if show_start:
        fig, (start_ax, end_ax) = plt.subplots(1, 2, figsize=(24, 12))
    else:
        fig, end_ax = plt.subplots(figsize=(12, 12))

    start_xs = np.arange(-past_edge, plot_up_to + 1)
    end_xs = np.arange(-plot_up_to, past_edge)

    for name, metagene_positions, color_index in data_sets:
        marker = '' if smooth else '.'
        linewidth = 2 if smooth else 1

        densities = metagene_positions['end_sum_of_enrichments']['three_prime_genomic'] / metagene_positions['end_num_eligible']['three_prime_genomic']
        if smooth:
            densities = smoothed(densities, 5)
        end_densities = densities['end', end_xs]
        end_ax.plot(end_xs,
                    end_densities,
                    '.-',
                    label=name,
                    color=colors[color_index],
                    marker=marker,
                    linewidth=linewidth,
                   )

        if show_start:
            densities = metagene_positions['start_sum_of_enrichments']['all'] / metagene_positions['start_num_eligible']['all']
            if smooth:
                densities = smoothed(densities, 5)
            start_densities = densities['start', start_xs]

            start_ax.plot(start_xs,
                          start_densities,
                          '.-',
                          label=name,
                          color=colors[color_index],
                          marker=marker,
                          linewidth=linewidth,
                         )
        
    
    if len(data_sets) > 1:
        start_ax.legend(loc='upper right', framealpha=0.5)

    if show_start:
        start_ax.set_xlabel('Nucleotides from start')
        start_ax.set_ylabel('Mean normalized read density')
        start_ax.set_xlim(min(start_xs), max(start_xs))
        start_ax.plot(start_xs, [1 for x in start_xs], color='black', alpha=0.5)
   
    end_ax.set_xlabel('Nucleotides from end')
    end_ax.yaxis.tick_right()
    end_ax.set_xlim(min(end_xs), max(end_xs))
    end_ax.plot(end_xs, [1 for x in end_xs], color='black', alpha=0.5)
    
    axs = [end_ax]
    if show_start:
        axs.append(end_ax)

    ymax = max(ax.get_ylim()[1] for ax in axs)
    
    for ax in axs:
        ax.set_ylim(0, ymax + 0.1)
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        ax.set_aspect((xmax - xmin) / (ymax - ymin))
        #ax.set_aspect(1)

    #fig.set_size_inches(12, 12)
    fig.savefig(figure_fn, bbox_inches='tight')
    plt.close(fig)

def plot_averaged_codon_densities(data_sets,
                                  figure_fn,
                                  show_end=True,
                                  past_edge=10,
                                  smooth=False,
                                  plot_up_to=100,
                                 ):
    if show_end:
        fig, (start_ax, end_ax) = plt.subplots(1, 2, figsize=(24, 12))
    else:
        fig, start_ax = plt.subplots(figsize=(12, 12))

    start_xs = np.arange(-past_edge, plot_up_to + 1)
    end_xs = np.arange(-plot_up_to, past_edge)

    for name, mean_densities, color_index in data_sets:
        densities = mean_densities['from_start']['codons']
        if smooth:
            densities = smoothed(densities, 5)
        start_densities = densities['start_codon', start_xs]

        marker = '' if smooth else '.'
        linewidth = 2 if smooth else 1

        start_ax.plot(start_xs,
                      start_densities,
                      '.-',
                      label=name,
                      color=colors[color_index],
                      marker=marker,
                      linewidth=linewidth,
                     )
        if show_end:
            densities = mean_densities['from_end']['codons']
            if smooth:
                densities = smoothed(densities, 5)
            end_densities = densities['stop_codon', end_xs]
            end_ax.plot(end_xs,
                        end_densities,
                        '.-',
                        label=name,
                        color=colors[color_index],
                        marker=marker,
                        linewidth=linewidth,
                       )
    
    if len(data_sets) > 1:
        start_ax.legend(loc='upper right', framealpha=0.5)

    start_ax.set_xlabel('Number of codons from start codon')
    start_ax.set_ylabel('Mean normalized read density')
    start_ax.set_xlim(min(start_xs), max(start_xs))
    start_ax.plot(start_xs, [1 for x in start_xs], color='black', alpha=0.5)
   
    if show_end:
        end_ax.set_xlabel('Number of codons from stop codon')
        end_ax.yaxis.tick_right()
        end_ax.set_xlim(min(end_xs), max(end_xs))
        end_ax.plot(end_xs, [1 for x in end_xs], color='black', alpha=0.5)
    
    axs = [start_ax]
    if show_end:
        axs.append(end_ax)

    ymax = max(ax.get_ylim()[1] for ax in axs)
    
    for ax in axs:
        ax.set_ylim(0, ymax + 0.1)
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        ax.set_aspect((xmax - xmin) / (ymax - ymin))
        #ax.set_aspect(1)

    #fig.set_size_inches(12, 12)
    fig.savefig(figure_fn, bbox_inches='tight')
    plt.close(fig)

def make_metacodon_panel(ax, metacodon_counts, codon_id, style, keys_to_plot):
    id_counts = metacodon_counts[codon_id]
    random_counts = id_counts.itervalues().next()
    xs = np.arange(-random_counts.left_buffer, random_counts.right_buffer)

    ax.set_title(codon_id)

    if style == 'enrichment':
        average_enrichment = id_counts['sum_of_enrichments'] / id_counts['num_eligible']
        ys = average_enrichment['codon', xs]
        ax.plot(xs, ys, '.-')
        
        ones = np.ones_like(xs)
        ax.plot(xs, ones, color='black', alpha=0.5)
    elif style == 'actual':
        for key in keys_to_plot:
            ys = id_counts[key]['codon', xs]
            line, = ax.plot(xs, ys, '.-', label=key)

            if isinstance(key, int):
                # -key + 1 is the position a read starts at if position 0
                # is the last base in it
                ax.axvline(-key + 1, color=line.get_color(), alpha=0.2)

        #ax.plot(xs, metacodon_counts[codon_id]['uniform'].counts, '.-')
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 3))

    ax.set_ylim(ymin=0)
    ax.set_xlim(min(xs), max(xs))
    ax.axvline(0, ls='--', color='black', alpha=0.5)

def plot_metacodon_counts(metacodon_counts,
                          fig_fn,
                          codon_ids='all',
                          style='actual',
                          keys_to_plot=['actual'],
                         ):
    bases = 'TCAG'

    if codon_ids == 'all':
        fig, axs = plt.subplots(16, 4, figsize=(16, 64))

        for first in range(4):
            for second in range(4):
                for third in range(4):
                    i = first * 4 + third
                    j = second
                    codon_id = ''.join(bases[first] + bases[second] + bases[third])
                    make_metacodon_panel(axs[i, j], metacodon_counts, codon_id, style, keys_to_plot)
    else:
        fig, axs = plt.subplots(len(codon_ids), 1, figsize=(8, 8 * len(codon_ids)))
        for codon_id, ax in zip(codon_ids, axs):
            make_metacodon_panel(ax, metacodon_counts, codon_id, style)

    fig.savefig(fig_fn, bbox_inches='tight')

def plot_metacodon_comparison(names,
                              metacodon_counts_list,
                              fig_fn,
                              style='actual',
                              keys_to_plot=['actual'],
                             ):
    codon_ids = codons.all_codons[:30]
    fig, axs = plt.subplots(len(names),
                            len(codon_ids),
                            figsize=(8 * len(codon_ids), 8 * len(names)),
                           )
    for n, name in enumerate(names):
        for c, codon_id in enumerate(codon_ids):
            make_metacodon_panel(axs[n, c], metacodon_counts_list[n], codon_id, style, keys_to_plot)

        axs[n, 0].set_ylabel(name)

    fig.savefig(fig_fn, bbox_inches='tight')

def plot_single_length_metacodon_counts(metacodon_counts,
                                        fig_fn,
                                        codon_ids,
                                        length=28,
                                        edge="5'",
                                       ):
    random_codon_id = metacodon_counts.iterkeys().next()
    random_counts = metacodon_counts[random_codon_id].itervalues().next()
    keys = metacodon_counts[random_codon_id].keys()
    xs = np.arange(-30, 30)
    minor_xticks = [x for x in xs if x % 3 == 0]
    major_xticks = [x for x in xs if x % 15 == 0]
    bases = 'TCAG'

    def make_plot(ax, codon_id):
        ax.set_title(codon_id)

        if codon_id not in metacodon_counts:
            return

        counts = metacodon_counts[codon_id][length]
        if edge == "3'":
            positions = xs - length + 1
        else:
            positions = xs
        line, = ax.plot(xs, counts['feature', positions], '.-')
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 3))

        ax.set_ylim(ymin=0)
        ax.set_xlim(min(xs), max(xs))
        #ax.axvline(0, ls='--', color='black', alpha=0.5)
        ax.set_xticks(major_xticks)
        for x in minor_xticks:
            ax.axvline(x, color='black', alpha=0.1)

    
    if codon_ids == 'all':
        fig, axs = plt.subplots(16, 4, figsize=(16, 64))

        for first in range(4):
            for second in range(4):
                for third in range(4):
                    i = first * 4 + third
                    j = second
                    codon_id = ''.join(bases[first] + bases[second] + bases[third])
                    make_plot(axs[i, j], codon_id)
    else:
        fig, axs = plt.subplots(len(codon_ids), 1, figsize=(8, 8 * len(codon_ids)))
        for codon_id, ax in zip(codon_ids, axs):
            make_plot(ax, codon_id)

    fig.savefig(fig_fn, bbox_inches='tight')

def plot_frameshifts(gtf_fn,
                     bam_fns,
                     gene_name,
                     exp_name,
                     genome_dir,
                     show_fractions=False,
                    ):
    codon_buffer = 10 
    start_codon = 'ATG'
    stop_codons = {'TAA', 'TAG', 'TGA'}
    A_site_offset = 5

    CDSs = {c.name: c for c in gtf.get_CDSs(gtf_fn, genome_dir, '/dev/null')}
    lengths = [28]
    
    transcript = CDSs[gene_name]

    left_buffer = 30 + codon_buffer * 3
    right_buffer = (codon_buffer + 1) * 3
    transcript.build_extent_maps(left_buffer, right_buffer)

    
    experiment_counts = []
    for bam_fn in bam_fns:
        counts = positions.get_Transcript_extent_position_counts(transcript,
                                                                 bam_fn,
                                                                 lengths,
                                                                 left_buffer=left_buffer,
                                                                 right_buffer=right_buffer,
                                                                )
        experiment_counts.append(counts)

    # Get the sequence of the extent.
    extent_sequence = transcript.get_extent_sequence(left_buffer=left_buffer,
                                                     right_buffer=right_buffer,
                                                    )

    for length in lengths:
        A_site_offset = positions.A_site_offsets['yeast'][length]

        length_counts = reduce(operator.add, [counts[length] for counts in experiment_counts])
        codon_numbers = np.arange(-codon_buffer, transcript.extent_length / 3 + codon_buffer)
        
        # frame_counts_list[i, j] will be the number of RPF's starting at frame i of
        # codon j
        frame_counts_list = np.zeros((3, len(codon_numbers)), int)
        start_codon_locations = [[] for _ in range(3)]
        stop_codon_locations = [[] for _ in range(3)]

        for c, codon_number in enumerate(codon_numbers):
            codon_start = 3 * codon_number
            for frame in range(3):
                frame_counts_list[frame, c] = length_counts['start', codon_start + frame - A_site_offset]
                codon = extent_sequence['start', codon_start + frame:codon_start + frame + 3]
                codon = ''.join(codon) 
                
                if codon == start_codon:
                    start_codon_locations[frame].append(codon_number)
                
                if codon in stop_codons:
                    stop_codon_locations[frame].append(codon_number)

        if show_fractions:
            fig, axs = plt.subplots(4, 1, sharex=True)
            cumulative_ax = axs[0]
            frame_axs = axs[1:]
        else:
            fig, frame_axs = plt.subplots(3, 1, sharex=True)

        for frame, (ax, frame_counts) in enumerate(zip(frame_axs, frame_counts_list)):
            nonzero_codon_numbers = [c_n for c_n, f_c in zip(codon_numbers, frame_counts) if f_c != 0]
            nonzero_frame_counts = [f_c for c_n, f_c in zip(codon_numbers, frame_counts) if f_c != 0]
            ax.plot(nonzero_codon_numbers, nonzero_frame_counts, '.')
            ax.set_ylim(0, frame_counts_list.max() + 1)
            ax.set_xlim(codon_numbers[0], codon_numbers[-1])
            ax.set_title('Frame {0}'.format(frame))
            ax.set_ylabel('Read counts')

            for x in start_codon_locations[frame]:
                ax.axvspan(x - 0.5, x + 0.5, facecolor='green', edgecolor='none', alpha=0.2)

            for x in stop_codon_locations[frame]:
                ax.axvspan(x - 0.5, x + 0.5, facecolor='red', edgecolor='none', alpha=0.2)

        frame_axs[-1].set_xlabel('Codons from start codon')

        if show_fractions:
            frames_so_far = frame_counts_list.cumsum(axis=1)
            fraction_frames_so_far = np.true_divide(frames_so_far,
                                                    np.maximum(1, frames_so_far.sum(axis=0)),
                                                   )

            frames_remaining = np.fliplr(np.fliplr(frame_counts_list).cumsum(axis=1))
            fraction_frames_remaining = np.true_divide(frames_remaining,
                                                       np.maximum(1, frames_remaining.sum(axis=0)),
                                                      )
            
            for frame in [0, 1, 2]:
                so_far = fraction_frames_so_far[frame]
                remaining = fraction_frames_remaining[frame]
                color = colors[frame]
                cumulative_ax.plot(codon_numbers, so_far, color=color, label='{0} so far'.format(frame))
                cumulative_ax.plot(codon_numbers, remaining, color=color, linestyle='--', label='{0} remaining'.format(frame))
            
            cumulative_ax.set_xlim(codon_numbers[0])
            cumulative_ax.set_ylim(-0.02, 1.02)
            cumulative_ax.set_ylabel('Fraction of reads in extent')
        
            cumulative_ax.legend(loc='upper right', framealpha=0.5)
        
        fig.suptitle('{2}\n{0}\nlength {1} fragments'.format(gene_name, length, exp_name))

        return fig

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

def plot_mismatch_type_by_position(type_counts, relevant_lengths, figure_fn):
    # To get the total number of reads for each length, get the total number of
    # base calls at read position 0 for that length. 
    length_totals = type_counts[:, 0].sum(axis=(1, 2, 3))
    nonzero_lengths = [length for length, count in enumerate(length_totals) if count > 0]
    
    fig, axs = plt.subplots(len(nonzero_lengths), 2, figsize=(24, 8 * len(nonzero_lengths)))
    
    base_order = utilities.base_order[:4]

    for length, (absolute_ax, zoomed_ax) in zip(nonzero_lengths, axs):
        length_counts = type_counts[length]
        all_rates = variants.compute_fractions_miscalled_as(length_counts, min_q=30)

        xs = np.arange(length)
        for ax in [absolute_ax, zoomed_ax]:
            for i, b in enumerate(base_order):
                ax.plot(xs,
                        all_rates[:length, i],
                        '.-',
                        color=igv_colors[b],
                        label=b,
                        markersize=10,
                        linewidth=0.5,
                       )
            ax.legend(framealpha=0.5)

            ax.set_xlim(-1, max(nonzero_lengths))
            ax.set_title('Length {0} - {1:,} total reads'.format(length, length_totals[length]))
            ax.axvline(0, color='black', alpha=0.5)
            ax.axvline(length - 1, color='black', alpha=0.5)
    
            ax.set_xlabel('Position in read')
        absolute_ax.set_ylim(0, 1)

    axs[len(axs) // 2, 0].set_ylabel('Fraction of bases miscalled as each type')

    fig.savefig(figure_fn, bbox_inches='tight')
    plt.close(fig)

if __name__ == '__main__1':
    from ribosome_profiling_experiment import RibosomeProfilingExperiment
    import subprocess

    job_fns = ['/home/jah/projects/ribosomes/experiments/belgium_2014_12_10/WT_1_FP/job/description.txt',
               '/home/jah/projects/ribosomes/experiments/belgium_2014_12_10/WT_2_FP/job/description.txt',
               '/home/jah/projects/ribosomes/experiments/belgium_2014_12_10/R98S_1_FP/job/description.txt',
               '/home/jah/projects/ribosomes/experiments/belgium_2014_12_10/R98S_2_FP/job/description.txt',
               '/home/jah/projects/ribosomes/experiments/belgium_2013_08_06/WT_cDNA_sample/job/description.txt',
               #'/home/jah/projects/ribosomes/experiments/belgium_2013_08_06/R98S_cDNA_sample/job/description.txt',
               #'/home/jah/projects/ribosomes/experiments/belgium_2013_08_06/Suppressed_R98S_cDNA_sample/job/description.txt',
              ]

    frameshift_genes = ['YIL009C-A',
                        'YOR239W',
                        'YPL052W',
                        #'YAL038W',
                       ]

    gene_fns = {gene: [] for gene in frameshift_genes}

    for job_fn in job_fns:
        experiment = RibosomeProfilingExperiment.from_description_file_name(job_fn)
        for gene in frameshift_genes:
            fig = plot_frameshifts(experiment.file_names['genes'],
                                   [experiment.file_names['merged_mappings']],
                                   gene,
                                   '{0}: {1}'.format(experiment.group, experiment.name),
                                   experiment.file_names['genome'],
                                  )
            fig.set_size_inches(16, 12)

            fn = '{0}/{1}.pdf'.format(experiment.work_results_dir, gene)
            fig.savefig(fn, bbox_inches='tight')
            plt.close(fig)
            gene_fns[gene].append(fn)
    
    for gene in frameshift_genes:
        merged_fn = '/home/jah/projects/ribosomes/experiments/belgium_2014_10_27/{0}.pdf'.format(gene)
        subprocess.check_call(['pdftk'] + gene_fns[gene] + ['cat', 'output', merged_fn])

if __name__ == '__main__1':
    import select_work
    experiments = select_work.build_all_experiments(verbose=False)
    chx_experiments = experiments['gerashchenko_nar']
    names = sorted(chx_experiments, key=select_work.gerashchenko_nar_sorting_key)
    metacodon_counts_list = [chx_experiments[name].read_file('metacodon_counts')
                             for name in names]

    plot_metacodon_comparison(names, metacodon_counts_list, 'test.png', style='enrichment')

if __name__ == '__main__':
    import select_work
    names = ['dom34KO_3-AT',
             #'dom34KO_CHX',
             #'dom34KO_mRNA-Seq',
             #'dom34KO_no_additive',
             'wild-type_3-AT',
             #'wild-type_CHX',
             #'wild-type_mRNA-Seq',
             #'wild-type_no_additive',
            ]

    experiments = select_work.build_all_experiments(verbose=False)
    relevant_experiments = [e for e in experiments['guydosh_cell'].values() if e.name in names]
    sorted_experiments = sorted(relevant_experiments, key=lambda e: e.name)

    #sorted_experiments = select_work.get_gerashchenko_nar_experiments()

    fig, ax = plt.subplots(figsize=(12, 12))

    special_codons = {'H': set(codons.full_back_table['H'])}
    #special_codons = {c : c for c in codons.full_back_table['H']}
    #special_codons = []

    color_iter = iter(colors)
    codon_counts_list = []
    for exp in sorted_experiments:
        print 'Reading {0}'.format(exp.name)
        codon_counts = exp.read_file('buffered_codon_counts', specific_keys={'relaxed', 'identities'})
        codon_counts_list.append(codon_counts)

    ratios_list, raw_counts_list = positions.compute_pause_scores(codon_counts_list, special_codons)

    for exp, ratios in zip(sorted_experiments, ratios_list):
        for key in ratios:
            if ratios[key] == []:
                continue
            sorted_values, cumulative = utilities.empirical_cdf(ratios[key])
            ax.semilogy(sorted_values, 1 - cumulative, color=color_iter.next(), label=exp.name + '_' + key)

    ax.legend(framealpha=0.5)
    ax.set_ylabel('Fraction of codons with enrichment >= x')
    ax.set_xlabel('Ratio of reads at codon to median reads across coding sequence')
    
    fig.savefig('../results/stalling/1_pausing_at_his.png', bbox_inches='tight')
    plt.close(fig)


    counts_around, ratios_around = positions.metacodon_around_pauses(codon_counts_list[0], special_codons=codons.full_back_table['H'], show_progress=True)
    quantiles = scipy.stats.mstats.mquantiles(ratios_around[:, 30], prob=np.linspace(0, 1, 11))
    boundaries = zip(quantiles, quantiles[1:])
    bins = [(ratios_around[:, 30] >= start) & (ratios_around[:, 30] < end) for start, end in boundaries]
    binned = [ratios_around[mask] for mask in bins] 

    fig, ax = plt.subplots(figsize=(12, 12))

    bin_means = np.asarray([np.mean(b, axis=0) for b in binned])[:, 30]
    before_means = np.asarray([np.mean(b, axis=0) for b in binned])[:, 30 - 10]
    two_before_means = np.asarray([np.mean(b, axis=0) for b in binned])[:, 30 - 20]
    control_means = np.asarray([np.mean(b, axis=0) for b in binned])[:, 30 + 10]

    ax.plot(bin_means, bin_means, 'o-', label='0')
    ax.plot(bin_means, before_means, 'o-', label='-10')
    ax.plot(bin_means, two_before_means, 'o-', label='-20')
    ax.plot(bin_means, control_means, 'o-', label='+10')
    for q in quantiles:
        ax.axvline(q, color='black', alpha=0.5)

    ax.set_xscale('log', basex=2)
    ax.set_yscale('log', basey=2)

    ax.set_xlabel('Average enrichment in bin at His')
    ax.set_ylabel('Average enrichment in bin at offset')
    ax.legend(loc='upper right', framealpha=0.5)
    ax.set_title('Stalling behind His codons is independant of read density at the His codon')

    fig.savefig('../results/stalling/5_mean_stalling_vs_mean_pausing.png', bbox_inches='tight')
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(12, 12))

    xs = np.arange(-30, 30)
    ys = ratios_around.mean(axis=0)

    ax.plot(xs, ys, '.-')
    ax.set_ylabel('Average enrichment')
    ax.set_xlabel('Offset relative to His codon')
    fig.savefig('../results/stalling/3_enrichment_around_his.png', bbox_inches='tight')
    plt.close(fig)
    
    fig, ax = plt.subplots(figsize=(12, 12))

    xs = np.arange(-30, 30)
    less_than = ratios_around[ratios_around[:, 30] <= 1].mean(axis=0)
    greater_than = ratios_around[ratios_around[:, 30] > 1].mean(axis=0)

    ax.plot(xs, less_than, '.-', label='His enrichment <= 1')
    ax.plot(xs, greater_than, '.-', label='His enrichment > 1')
    ax.legend(loc='upper right', framealpha=0.5)
    ax.set_ylabel('Average enrichment')
    ax.set_xlabel('Offset relative to His codon')
    fig.savefig('../results/stalling/4_enrichment_around_his_stratified.png', bbox_inches='tight')
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(12, 12))

    Sequencing.Visualize.enhanced_scatter(np.maximum(0.01, ratios_list[0]['H']),
                                          np.maximum(0.01, ratios_list[1]['H']),
                                          relevant_experiments[0].name,
                                          relevant_experiments[1].name,
                                          'Position-specific enrichments are consistent across replicates',
                                          ax,
                                          do_fit=False,
                                          show_p_value=False,
                                         )
    ax.set_yscale('log')
    ax.set_xscale('log')
    
    fig.savefig('../results/stalling/2_pausing_reproducibility.png', bbox_inches='tight')
    plt.close(fig)
