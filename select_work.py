import os
import glob
import ribosome_profiling_experiment
import subprocess
import numpy as np
import visualize
import contaminants
from collections import Counter

def build_all_experiments(verbose=True):
    families = ['zinshteyn_plos_genetics',
                'ingolia_science',
                'weinberg',
                'dunn_elife',
                'gerashchenko_pnas',
                'gerashchenko_nar',
                'guydosh_cell',
                'mcmanus_gr',
                'artieri',
                'lareau_elife',
                'belgium_2014_08_07',
                'belgium_2014_03_05',
                'belgium_2013_08_06',
               ]
    experiments = {}
    for family in families:
        if verbose:
            print family
        experiments[family] = {}
        prefix = '/home/jah/projects/arlen/experiments/{0}/'.format(family)
        dirs = [path for path in glob.glob('{}*'.format(prefix)) if os.path.isdir(path)]
        for d in sorted(dirs):
            _, name = os.path.split(d)
            if verbose:
                print '\t', name
            description_file_name = '{0}/job/description.txt'.format(d)
            experiments[family][name] = ribosome_profiling_experiment.RibosomeProfilingExperiment.from_description_file_name(description_file_name)

    return experiments

def read_counts_and_RPKMS():
    experiments = build_all_experiments()
    for family in experiments:
        print family
        for name in experiments[family]:
            print '\t', name
            experiments[family][name].compute_total_read_counts()
            experiments[family][name].compute_RPKMs()

def package_files(key):
    prefix = '/home/jah/projects/arlen/'
    os.chdir(prefix)
    full_file_names = []
    package_file_name = 'all_{}.tar.gz'.format(key)

    experiments = build_all_experiments()
    for family in experiments:
        if 'belgium' in family:
            continue
        for name in experiments[family]:
            full_file_names.append(experiments[family][name].file_names[key])

    def strip_prefix(fn, prefix):
        if not fn.startswith(prefix):
            raise ValueError(fn)
        return fn[len(prefix):]

    relative_file_names = [strip_prefix(fn, prefix) for fn in full_file_names]
    tar_command = ['tar', '-czf', package_file_name] + relative_file_names
    subprocess.check_call(tar_command)

def make_counts_array_file(exclude_edges=False):
    prefix = '/home/jah/projects/arlen/'
    os.chdir(prefix)
    if exclude_edges:
        fn = 'all_read_counts_exclude_edges.txt'
    else:
        fn = 'all_read_counts.txt'
    read_counts = {}
    full_experiments = []
    experiments = build_all_experiments()
    for family in sorted(experiments):
        read_counts[family] = {}
        for name in sorted(experiments[family]):
            full_experiment = '{0}:{1}'.format(family, name)
            full_experiments.append(full_experiment)
            if exclude_edges:
                read_counts[family][name] = experiments[family][name].read_file('read_counts_exclude_edges')
            else:
                read_counts[family][name] = experiments[family][name].read_file('read_counts')

    gene_names = sorted(read_counts[family][name].keys())
    gene_lengths = [read_counts[family][name][gene_name]['CDS_length'] for gene_name in gene_names]

    full_array = [gene_lengths]

    for full_experiment in full_experiments:
        family, name = full_experiment.split(':')
        counts = [read_counts[family][name][gene_name]['expression'][0] for gene_name in gene_names]
        full_array.append(counts)

    full_array = np.asarray(full_array).T

    with open(fn, 'w') as fh:
        fh.write('name\tlength\t{0}\n'.format('\t'.join(full_experiments)))
        for gene_name, row in zip(gene_names, full_array):
            fh.write('{0}\t'.format(gene_name))
            fh.write('{0}\n'.format('\t'.join(map(str, row))))

def make_restricted_starts_and_ends_plots():
    all_experiments = build_all_experiments(verbose=False)

    relevant_lengths = range(19, 25)
    for name in all_experiments['gerashchenko_pnas']:
        if 'rep1' in name and 'foot' in name and 'Initial' not in name:
            print name
            experiment = all_experiments['gerashchenko_pnas'][name]
            #experiment.plot_starts_and_ends()
            experiment.plot_mismatch_types()
            #position_counts = experiment.read_file('from_starts_and_ends')
            #visualize.plot_metagene_positions(position_counts['from_starts'],
            #                                  position_counts['from_ends'],
            #                                  experiment.figure_file_names['starts_and_ends'],
            #                                  relevant_lengths=relevant_lengths,
            #                                 )
            #visualize.plot_metacodon_positions(position_counts['from_starts'],
            #                                   experiment.figure_file_names['starts_and_ends'],
            #                                   key='start_codon',
            #                                  )

def make_mismatch_position_plots():
    all_experiments = build_all_experiments(verbose=False)

    for group in all_experiments:
        for name in all_experiments[group]:
            if group == 'dunn_elife':
                continue
            print name
            experiment = all_experiments[group][name]
            #experiment.plot_mismatch_positions()
            #experiment.plot_starts_and_ends()
            experiment.plot_lengths()

def make_multipage_pdf(figure_name):
    all_experiments = build_all_experiments(verbose=False)
    all_fn = '/home/jah/projects/arlen/results/guydosh_{0}.pdf'.format(figure_name)
    fns = []
    for group in sorted(all_experiments):
        if 'guydosh' not in group:
            continue
        for name in sorted(all_experiments[group]):
            fns.append(all_experiments[group][name].figure_file_names[figure_name])

    pdftk_command = ['pdftk'] + fns + ['cat', 'output', all_fn]
    subprocess.check_call(pdftk_command)

def get_read_lengths():
    all_experiments = build_all_experiments(verbose=False)
    read_lengths = Counter()
    for group in sorted(all_experiments):
        for name in sorted(all_experiments[group]):
            experiment = all_experiments[group][name]
            read_lengths[experiment.max_read_length] += 1
            if experiment.max_read_length == 76:
                print group, name
    print read_lengths.most_common()

def make_averaged_codon_densities_plot():
    experiments = build_all_experiments(verbose=False)
    
    def transform(experiment):
        _, concentration = experiment.name.split('_', 1)
        concentration, _ = concentration.split('CHX')
        if concentration == 'no':
            concentration = 0
        else:
            concentration = concentration.strip('_x')
            if '_' in concentration:
                num, denom = concentration.split('_')
                concentration = float(num) / float(denom)
            else:
                concentration = int(concentration)

        return concentration

    sorted_experiments = sorted(experiments['gerashchenko_nar'].values(), key=transform)
    data_sets = [(experiment.name, experiment.read_file('mean_densities'), i)
                 for i, experiment in enumerate(sorted_experiments)]

    visualize.plot_averaged_codon_densities(data_sets,
                                            'test.pdf',
                                            past_edge=10,
                                            plot_up_to=100,
                                            show_end=True,
                                           )

def make_rRNA_coverage_plot():
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k'] * 10
    experiments = build_all_experiments(verbose=False)
    all_experiments = [exp for exp in experiments['belgium_2014_08_07'].values() if 'FP' in exp.name]

    coverage_data = {exp.name: (exp.get_total_reads(), exp.read_file('rRNA_coverage'), color)
                     for exp, color in zip(all_experiments, colors)}

    contaminants.plot_rRNA_coverage(coverage_data,
                                    all_experiments[0].file_names['oligos_sam'],
                                    'belgium_2014_08_07_rRNA_coverage_{0}.pdf',
                                   )
