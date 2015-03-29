import ribosomes
import numpy as np
import matplotlib.pyplot as plt
import gtf
import explore_UTRs
import select_work
import Sequencing.Visualize
from collections import Counter

bad_gene_names = {'YER109C', 'YOR031W'}

experiments = select_work.build_all_experiments(verbose=False)

d = {'belgium_2014_12_10': {'WT_1': {'mRNA': 'WT_1_mRNA',
                                     'RPF': 'WT_1_FP',
                                    },
                            'WT_2': {'mRNA': 'WT_2_mRNA',
                                     'RPF': 'WT_2_FP',
                                    },
                            'R98S_1': {'mRNA': 'R98S_1_mRNA',
                                       'RPF': 'R98S_1_FP'
                                      },
                            'R98S_2': {'mRNA': 'R98S_2_mRNA',
                                       'RPF': 'R98S_2_FP'
                                      },
                           },
     'weinberg': {'Unselected': {'mRNA': 'Unselected',
                                 'RPF': 'RPF',
                                },
                  'Dynabeads': {'mRNA': 'Dynabeads',
                                'RPF': 'RPF',
                               },
                 },
    }

random_group = d.keys().pop()
random_condition = d[random_group].keys().pop()
random_experiment = experiments[random_group][d[random_group][random_condition]['RPF']]

transcripts, _ = random_experiment.get_CDSs()
gene_names = [t.name for t in transcripts]
CDS_lengths = [t.CDS_length for t in transcripts]

def experiment_to_read_counts(experiment):
    read_counts = experiment.read_file('read_counts')
    array = np.asarray([read_counts[gene]['expression'][0] for gene in gene_names])
    nonzero_array = np.maximum(array, 0.1)
    return nonzero_array

counts = {}

for group in d:
    counts[group] = {}
    for condition in d[group]:
        counts[group][condition] = {}
        for kind in d[group][condition]:
            sample = d[group][condition][kind]
            experiment = experiments[group][sample]
            counts[group][condition][kind] = experiment_to_read_counts(experiment)

def make_counts_array_file():
    fn = 'all_counts.txt'
    read_counts = {}
    full_experiments = []
    for experiment in sorted(read_counts_fns):
        read_counts[experiment] = {}
        for kind in sorted(read_counts_fns[experiment]):
            full_experiment = '{0}:{1}'.format(experiment, kind)
            full_experiments.append(full_experiment)
            read_counts[experiment][kind] = Serialize.read_file(read_counts_fns[experiment][kind], 'expression')

    gene_names = sorted(read_counts[experiment][kind].keys())
    gene_lengths = [read_counts[experiment][kind][name]['CDS_length'] for name in gene_names]

    full_array = [gene_lengths]

    for full_experiment in full_experiments:
        experiment, kind = full_experiment.split(':')
        counts = [read_counts[experiment][kind][name]['expression'][0] for name in gene_names]
        full_array.append(counts)

    full_array = np.asarray(full_array).T

    with open(fn, 'w') as fh:
        fh.write('name\tCDS_length\t{0}\n'.format('\t'.join(full_experiments)))
        for gene_name, row in zip(gene_names, full_array):
            fh.write('{0}\t'.format(gene_name))
            fh.write('{0}\n'.format('\t'.join(map(str, row))))

def read_counts_array_file():
    fn = 'all_counts.txt'
    fh = open(fn)
    fields = {name: i - 1 for i, name in enumerate(fh.readline().strip().split())}
    full_array = []
    for line in fh:
        full_array.append(map(int, line.strip().split()[1:]))

    return fields, np.asarray(full_array)

def filtered_low_counts(min_count):
    fields, array = read_counts_array_file()
    at_least = np.all(array[:, 1:] > min_count, axis=1)
    filtered = np.asarray([row for row, accept in zip(array, at_least) if accept])
    return fields, filtered

def where_all_at_least(array, min_count):
    where_at_least = np.where(np.min(array, axis=0) > min_count)

    return where_at_least

def compare_TEs():
    names = [('belgium_2014_12_10', 'WT_1'),
             ('belgium_2014_12_10', 'WT_2'),
             ('belgium_2014_12_10', 'R98S_1'),
             ('belgium_2014_12_10', 'R98S_2'),
            ]

    all_mRNA_counts = []
    all_RPF_counts = []
    for group, condition in names:
        mRNA_counts = counts[group][condition]['mRNA']
        RPF_counts = counts[group][condition]['RPF']

        all_mRNA_counts.append(mRNA_counts)
        all_RPF_counts.append(RPF_counts)

    where_at_least = where_all_at_least(np.asarray(all_mRNA_counts + all_RPF_counts), 10)

    all_TEs = []
    for i, (group, condition) in enumerate(names):
        mRNA_normalized = np.true_divide(all_mRNA_counts[i], all_mRNA_counts[i].sum())
        RPF_normalized = np.true_divide(all_RPF_counts[i], all_RPF_counts[i].sum())
        TEs = (np.log2(RPF_normalized) - np.log2(mRNA_normalized))[where_at_least]
        all_TEs.append(TEs)
    
    labels = [condition for group, condition in names]
    full_labels = ['{1} TE (log2 (mRNA RPKM / RPF RPKM))'.format(group, condition) for group, condition in names]
    
    fig = plt.figure(figsize=(9 * len(names), 9 * len(names)))
    
    for r in range(len(names)):
        for c in range(r, len(names)):
            ax = fig.add_subplot(len(names), len(names), r * len(names) + c + 1)

            xs = all_TEs[c]
            ys = all_TEs[r]

            if r == c:
                x_label = full_labels[c]
                y_label= full_labels[r]
            else:
                x_label = labels[c]
                y_label= labels[r]

            Sequencing.Visualize.enhanced_scatter(xs, ys, x_label, y_label, '', ax, lims=(-4, 4))

    fig.savefig('/home/jah/projects/ribosomes/results/TE_correlations.png', bbox_inches='tight')

def TE_vs_length():
    fields, filtered = filtered_low_counts(500)
    log10lengths = np.log10(filtered[:, fields['CDS_length']])
    lengths = np.asarray(filtered[:, fields['CDS_length']], dtype=float)
    print len(lengths)

    TEs = {#'Ingolia': np.log2(filtered[:, fields['ingolia:RPF_1']] + filtered[:, fields['ingolia:RPF_2']]) - np.log2(filtered[:, fields['ingolia:mRNA_1']] + filtered[:, fields['ingolia:mRNA_2']]),
           #'Weinberg_RiboZero': np.log2(filtered[:, fields['weinberg:RPF']]) - np.log2(filtered[:, fields['weinberg:RiboZero']]),
           #'Weinberg_Dynabeads': np.log2(filtered[:, fields['weinberg:RPF']]) - np.log2(filtered[:, fields['weinberg:Dynabeads']]),
           #'Weinberg_Unselected': np.log2(filtered[:, fields['weinberg:RPF']]) - np.log2(filtered[:, fields['weinberg:Unselected']]),
           #'Mcmanus': np.log2(filtered[:, fields['mcmanus:RPF_1']]) - np.log2(filtered[:, fields['mcmanus:mRNA_1']]),
           #'Arlen': np.log2(filtered[:, fields['belgium:RPF']]) - np.log2(filtered[:, fields['belgium:mRNA']]),
           #'Gerashchenko': np.log2(filtered[:, fields['gerashchenko:RPF']]) - np.log2(filtered[:, fields['gerashchenko:mRNA']]),
           #'Guydosh': np.log2(filtered[:, fields['guydosh:RPF']]) - np.log2(filtered[:, fields['guydosh:mRNA']]),
           'Zinshteyn': np.log2(filtered[:, fields['zinshteyn:RPF']]) - np.log2(filtered[:, fields['zinshteyn:mRNA']]),

           #'Weinberg_Unselected_mRNA': np.log2(filtered[:, fields['weinberg:Unselected']] / lengths),
           #'Weinberg_Dynabeads_mRNA': np.log2(filtered[:, fields['weinberg:Dynabeads']]) - np.log2(filtered[:, fields['weinberg:Unselected']]),
           #'Weinberg_RiboZero_mRNA': np.log2(filtered[:, fields['weinberg:RiboZero']]) - np.log2(filtered[:, fields['weinberg:Unselected']]),
           #'Weinberg_RPF': np.log2(filtered[:, fields['weinberg:RPF']] / lengths),
          }

    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    name = 'Zinshteyn'
    colored_scatter(log10lengths,
                    TEs[name] - np.mean(TEs[name]),
                    'log10(CDS length)',
                    'ribosome density - log2(RPF RPKM / mRNA RPKM)',
                    name,
                    ax,
                   )
    return

    #fig, axs = plt.subplots(2, 3, figsize=(24, 16))
    #for name, ax in zip(['Ingolia', 'Mcmanus', 'Arlen', 'Gerashchenko', 'Guydosh', 'Zinshteyn'], axs.flatten()):
    #    colored_scatter(log10lengths,
    #                    TEs[name] - np.mean(TEs[name]),
    #                    'log10(length)',
    #                    'ribosome density - log2(RPF RPKM / mRNA RPKM)',
    #                    name,
    #                    ax,
    #                   )

    #fig.savefig('/home/jah/scratch/joint_group_5_14/six_length_correlations.png')
    #return

    num_rows, num_cols = 1, 3
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(8 * num_cols, 8 * num_rows))
    for name, ax in zip(['Weinberg_Dynabeads', 'Weinberg_RiboZero', 'Weinberg_Unselected'], axs.flatten()):
        colored_scatter(log10lengths,
                        TEs[name] - np.mean(TEs[name]),
                        'log10(length)',
                        'ribosome density - log2(RPF RPKM / mRNA RPKM)',
                        name,
                        ax,
                       )
    fig.savefig('/home/jah/scratch/joint_group_5_14/weinberg_TEs.png')
    return
        
    #fig, axs = plt.subplots(1, 2, figsize=(16, 8))
    #for name, ax in zip(['Weinberg_Dynabeads_mRNA', 'Weinberg_RiboZero_mRNA'], axs.flatten()):
    #    colored_scatter(log10lengths,
    #                    TEs[name] - np.mean(TEs[name]),
    #                    'log10(length)',
    #                    'log2( mRNA RPKM / Weinberg_Unselected RPKM)',
    #                    name,
    #                    ax,
    #                   )
    #    ax.set_ylim(-3, 3)
    #fig.savefig('/home/jah/scratch/joint_group_5_14/weinberg_vs_unselected.png')
    
def new_TE_distribution():
    fields, filtered = filtered_low_counts(0)

    TEs = {'ingolia_1': np.log2(filtered[:, fields['ingolia:RPF_1']]) - np.log2(filtered[:, fields['ingolia:mRNA_1']]),
           'ingolia_2': np.log2(filtered[:, fields['ingolia:RPF_2']]) - np.log2(filtered[:, fields['ingolia:mRNA_2']]),
           'ingolia_both': np.log2(filtered[:, fields['ingolia:RPF_1']] + filtered[:, fields['ingolia:RPF_2']]) - np.log2(filtered[:, fields['ingolia:mRNA_1']] + filtered[:, fields['ingolia:mRNA_2']]),
           'weinberg_RiboZero': np.log2(filtered[:, fields['weinberg:RPF']]) - np.log2(filtered[:, fields['weinberg:RiboZero']]),
           'weinberg_Dynabeads': np.log2(filtered[:, fields['weinberg:RPF']]) - np.log2(filtered[:, fields['weinberg:Dynabeads']]),
           'weinberg_Unselected': np.log2(filtered[:, fields['weinberg:RPF']]) - np.log2(filtered[:, fields['weinberg:Unselected']]),
           'artificial': np.log2(filtered[:, fields['weinberg:RPF']]) - np.log2(filtered[:, fields['weinberg:Dynabeads']] / np.asarray(filtered[:, fields['CDS_length']], dtype=float)),
           'artificial2': np.log2(filtered[:, fields['weinberg:RPF']]) - np.log2(filtered[:, fields['weinberg:Dynabeads']] * np.asarray(filtered[:, fields['CDS_length']], dtype=float)),
          }

    for name in ['weinberg_RiboZero', 'weinberg_Dynabeads', 'weinberg_Unselected', 'ingolia_both', 'artificial', 'artificial2']:
        plt.hist(TEs[name] - np.mean(TEs[name]), histtype='step', bins=100, range=(-4, 4), label=name)

    plt.legend()
    plt.xlabel('log2(RPF RPKM / mRNA RPKM')
    plt.ylabel('Number of genes')

    explore_UTRs.scatter_with_hists_colors(TEs['weinberg_RiboZero'] - np.mean(TEs['weinberg_RiboZero']),
                                           TEs['weinberg_Unselected'] - np.mean(TEs['weinberg_Unselected']),
                                           'weinberg_Ribozero',
                                           'weinberg_Unselected',
                                           'Joint distribution of TEs',
                                          )

    print scipy.stats.pearsonr(TEs['weinberg_RiboZero'], TEs['weinberg_Unselected'])
    print scipy.stats.spearmanr(TEs['weinberg_RiboZero'], TEs['weinberg_Unselected'])

    return TEs

def plot_RPKMs_scatter():
    names = ['Unselected',
             'Dynabeads',
             'RiboMinus',
             'RiboZero',
            ]

    vals = [arrays['weinberg'][name] for name in names]

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

def mRNA_RPKM_length_bias():
    gtf_fn = '/home/jah/projects/arlen/data/organisms/saccharomyces_cerevisiae/EF4/transcriptome/genes.gtf'
    genome_dir = '/home/jah/projects/arlen/data/organisms/saccharomyces_cerevisiae/EF4/genome'
    coding_sequence_fetcher = gtf.make_coding_sequence_fetcher(gtf_fn, genome_dir)
    lengths = np.asarray([len(coding_sequence_fetcher(name)) for name in gene_names])

    explore_UTRs.scatter_with_hists_colors(lengths,
                                           #np.log2(arrays['weinberg']['Dynabeads']) - np.log2(arrays['weinberg']['Unselected']),
                                           np.log2(arrays['weinberg']['RPF']) - np.log2(arrays['ingolia']['RPF']),
                                           'coding sequence length',
                                           'log2(Ingolia mRNA RPKM / Weinberg mRNA RPKM)',
                                           '',
                                          )

    plt.ylim(-7, 7)
    plt.gcf().set_size_inches(12, 8)
    #plt.savefig('mRNA_RPKM_length_bias.pdf')
    #plt.savefig('mRNA_RPKM_length_bias.png')
    #plt.close()

def TE_length_bias():
    gtf_fn = '/home/jah/projects/arlen/data/organisms/saccharomyces_cerevisiae/EF4/transcriptome/genes.gtf'
    genome_dir = '/home/jah/projects/arlen/data/organisms/saccharomyces_cerevisiae/EF4/genome'
    coding_sequence_fetcher = gtf.make_coding_sequence_fetcher(gtf_fn, genome_dir)
    lengths = np.asarray([len(coding_sequence_fetcher(name)) for name in gene_names])
    explore_UTRs.scatter_with_hists_colors(lengths,
                                           #TEs['ingolia']['mRNA'] - TEs['weinberg']['RiboZero'],
                                           #TEs['ingolia2']['mRNA'] - TEs['ingolia']['mRNA'],
                                           #TEs['weinberg']['Dynabeads'] - TEs['weinberg']['Unselected'],
                                           TEs['weinberg']['Dynabeads'] - TEs['weinberg']['RiboZero'],
                                           'coding sequence length',
                                           'log2(Ingolia TE / Weinberg TE)',
                                           '',
                                          )

    plt.ylim(-7, 7)
    plt.gcf().set_size_inches(12, 8)
    #plt.savefig('TE_length_bias.pdf')
    #plt.savefig('TE_length_bias.png')
    #plt.close()

def TE_distribution():
    plt.figure()
    plt.hist(TEs['ingolia']['mRNA'] - np.mean(TEs['ingolia']['mRNA']), histtype='step', range=(-5, 5), bins=100, label='Ingolia 1 TEs')
    plt.hist(TEs['ingolia2']['mRNA'] - np.mean(TEs['ingolia2']['mRNA']), histtype='step', range=(-5, 5), bins=100, label='Ingolia 2 TEs')
    #plt.hist(TEs['ingolia2']['mRNA'] - TEs['ingolia']['mRNA'], histtype='step', range=(-5, 5), bins=100, label='Ingolia 1 / 2')
    plt.hist(TEs['weinberg']['Dynabeads'] - np.mean(TEs['weinberg']['Dynabeads']), histtype='step', range=(-5, 5), bins=100, label='Weinberg Dynabeads TEs')
    plt.hist(TEs['weinberg']['RiboZero'] - np.mean(TEs['weinberg']['RiboZero']), histtype='step', range=(-5, 5), bins=100, label='Weinberg RiboZero TEs')
    #plt.hist(TEs['weinberg']['RiboZero'] - TEs['weinberg']['Dynabeads'], histtype='step', range=(-5, 5), bins=50, label='Weinberg RiboZero TEs')
    #plt.hist(TEs['weinberg']['RiboMinus'] - np.mean(TEs['weinberg']['RiboMinus']), histtype='step', range=(-5, 5), bins=100, label='Weinberg RiboMinus TEs')
    #plt.hist(TEs['weinberg']['Unselected'] - np.mean(TEs['weinberg']['Unselected']), histtype='step', range=(-5, 5), bins=100, label='Weinberg Unselected TEs')
    plt.xlim(-5, 5)
    plt.xlabel('log2 TE (RPF RPKM / mRNA RPKM)')
    plt.ylabel('Number of genes')
    plt.legend()
    plt.gcf().set_size_inches(12, 8)
    #plt.savefig('translation_regulation.pdf')
    #plt.savefig('translation_regulation.png')
