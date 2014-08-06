import TIF_seq_experiment
import three_p_experiment
import three_t_fill_experiment
import ribosome_profiling_experiment
import TL_seq_experiment
import os
import glob
import matplotlib.pyplot as plt
import numpy as np
import h5py
import Serialize.read_positions
import gtf
import visualize
from collections import defaultdict
from itertools import izip

def build_all_experiments(verbose=True):
    families = {'TIF_seq': (TIF_seq_experiment.TIFSeqExperiment, ['pelechano_nature']),
                'three_p_seq': (three_p_experiment.ThreePExperiment, ['three_p_seq']),
                'TL_seq': (TL_seq_experiment.TLSeqExperiment, ['arribere_gr', 'park_nar']),
                'three_t_fill_seq': (three_t_fill_experiment.ThreeTFillExperiment, ['wilkening_nar']),
                'weinberg': (ribosome_profiling_experiment.RibosomeProfilingExperiment, ['weinberg']),
               }

    experiments = {}
    for kind in families:
        experiments[kind] = {}
        Experiment, kind_families = families[kind]
        for family in kind_families:
            if verbose:
                print family
            experiments[kind][family] = {}
            prefix = '/home/jah/projects/arlen/experiments/{0}/'.format(family)
            dirs = [path for path in glob.glob('{}*'.format(prefix)) if os.path.isdir(path)]
            for d in sorted(dirs):
                _, name = os.path.split(d)
                if verbose:
                    print '\t', name
                description_file_name = '{0}/job/description.txt'.format(d)
                experiments[kind][family][name] = Experiment.from_description_file_name(description_file_name)

    return experiments

def test_joint_plot():
    gene_name = 'YBR068C'
    experiments = build_all_experiments(verbose=False)
    CDSs, _ = experiments['TIF_seq']['pelechano_nature']['ypd_bio1_lib1'].get_CDSs()
    transcripts = {t.name: t for t in CDSs}
    names = ['ypd_bio1_lib1', 'nypd_bio2_lib1', 'ypd_bio1_lib4']
    fig, axs = plt.subplots(len(names), 1)
    for name, ax in zip(names, axs):
        experiment = experiments['TIF_seq']['pelechano_nature'][name]
        joint_position_counts = experiment.read_file('joint_positions')
        visualize.plot_joint_positions_scatter(joint_position_counts[gene_name],
                                               transcripts[gene_name],
                                               name,
                                               ax=ax,
                                              )

def find_internal_polyA():
    experiments = build_all_experiments(verbose=False)
    experiment = experiments['TIF_seq']['pelechano_nature']['ypd_bio1_lib1']
    CDSs, _ = experiment.get_CDSs()
    for transcript in CDSs:
        transcript.build_coordinate_maps()
        fn = experiment.file_names['three_prime_read_positions']
        f = h5py.File(fn, 'r')
        gene = Serialize.read_positions.build_gene(f[transcript.name])
        num_internal = gene['all']['stop_codon', -transcript.CDS_length:0].sum()
        if num_internal > 50:
            print transcript.name, gene['all']['stop_codon', -transcript.CDS_length:0].sum()
            print transcript.seqname, transcript.start, transcript.end
            raw_input()


_handles = defaultdict(list)
_labels = defaultdict(list)
def plot_with_line_alpha(ax, xs, ys, **kwargs):
    label = kwargs.pop('label', None)
    line,  = ax.plot(xs, ys, '-', alpha=0.3, **kwargs)
    nonzero_xs = [x for x, y in izip(xs, ys) if y != 0]
    nonzero_ys = [y for x, y in izip(xs, ys) if y != 0]
    circles, = ax.plot(nonzero_xs, nonzero_ys, '.', color=line.get_color(), **kwargs)

    _handles[ax].append((line, circles)) 
    _labels[ax].append(label) 
    ax.legend(_handles[ax], _labels[ax])

if __name__ == '__main__':
    gene_name = 'YAL035W'

    experiments = build_all_experiments(verbose=False)

    five_prime_experiments = [(name, experiment) for name, experiment in sorted(experiments['TL_seq']['arribere_gr'].items()) if 'TLSeq1' in name] + \
                             [(name, experiment) for name, experiment in sorted(experiments['TL_seq']['park_nar'].items()) if name == 'SMORE-seq_WT_TAP+_rep1'] + \
                             [(name, experiment) for name, experiment in sorted(experiments['TIF_seq']['pelechano_nature'].items()) if name == 'ypd_bio1_lib1']
                             
    three_prime_experiments = [(name, experiment) for name, experiment in sorted(experiments['three_p_seq']['three_p_seq'].items())] + \
                              [(name, experiment) for name, experiment in sorted(experiments['TIF_seq']['pelechano_nature'].items()) if name == 'ypd_bio1_lib1'] + \
                              [(name, experiment) for name, experiment in sorted(experiments['three_t_fill_seq']['wilkening_nar'].items()) if '3tfill_ypd_rep1' in name]
    
    rna_seq_experiments = [(name, experiment) for name, experiment in sorted(experiments['weinberg']['weinberg'].items()) if 'RPF' not in name]

    fig, ax = plt.subplots()

    for name, experiment in five_prime_experiments:
        fn = experiment.file_names['five_prime_read_positions']
        f = h5py.File(fn, 'r')
        gene = Serialize.read_positions.build_gene(f[gene_name])

        xs = np.arange(-200, gene['all'].CDS_length)
        counts = gene['all']['start_codon', xs]
        #normalization = float(experiment.get_total_eligible_reads()) / 1.e6
        #normalization = 1.
        normalization = counts.sum()
        normalized = np.true_divide(counts, normalization)
        plot_with_line_alpha(ax, xs, normalized, label=name)
    
    fig, ax = plt.subplots()
    
    for name, experiment in three_prime_experiments:
        fn = experiment.file_names['three_prime_read_positions']
        f = h5py.File(fn, 'r')
        gene = Serialize.read_positions.build_gene(f[gene_name])

        xs = np.arange(-gene['all'].CDS_length, 200)
        total_reads = experiment.get_total_eligible_reads()
        #total_reads = 1.e6
        print name, total_reads

        if 0 in gene:
            genomic_counts = gene[0]['stop_codon', xs]
            nongenomic_counts = gene['all']['stop_codon', xs] - gene[0]['stop_codon', xs]
            #normalization = float(experiment.get_total_eligible_reads()) / 1.e6
            #normalization = 1.
            normalization = (genomic_counts + nongenomic_counts).sum()
            genomic_normalized = np.true_divide(genomic_counts, normalization)
            nongenomic_normalized = np.true_divide(nongenomic_counts, normalization)

            plot_with_line_alpha(ax, xs, genomic_normalized, label='{}_genomic'.format(name))
            plot_with_line_alpha(ax, xs, nongenomic_normalized, label='{}_nongenomic'.format(name))
        else:
            counts = gene['all']['stop_codon', xs]
            #normalization = float(experiment.get_total_eligible_reads()) / 1.e6
            #normalization = 1.
            normalization = counts.sum()
            normalized = np.true_divide(counts, normalization)
            plot_with_line_alpha(ax, xs, normalized, label=name)
