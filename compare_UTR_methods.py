import TIF_seq_experiment
import three_p_experiment
import three_t_fill_experiment
import ribosome_profiling_experiment
import TL_seq_experiment
import os
import matplotlib.pyplot as plt
import numpy as np
import h5py
import Serialize.read_positions
import gff
import glob
import positions
import visualize
import Sequencing.utilities as utilities
import Sequencing.genomes as genomes
from collections import defaultdict, Counter
from itertools import izip

def build_all_experiments(verbose=True):
    families = {'TIF_seq': (TIF_seq_experiment.TIFSeqExperiment, ['pelechano_nature']),
                'three_p_seq': (three_p_experiment.ThreePExperiment, ['three_p_seq']),
                'TL_seq': (TL_seq_experiment.TLSeqExperiment, ['arribere_gr', 'park_nar']),
                'three_t_fill_seq': (three_t_fill_experiment.ThreeTFillExperiment, ['wilkening_nar']),
                'ribosome_profiling': (ribosome_profiling_experiment.RibosomeProfilingExperiment, ['weinberg', 'guydosh_cell']),
               }

    experiments = {}
    for kind in families:
        experiments[kind] = {}
        Experiment, kind_families = families[kind]
        for family in kind_families:
            if verbose:
                print family
            experiments[kind][family] = {}
            prefix = '/home/jah/projects/ribosomes/experiments/{0}/'.format(family)
            dirs = [path for path in glob.glob('{}*'.format(prefix)) if os.path.isdir(path)]
            for d in sorted(dirs):
                _, name = os.path.split(d)
                if verbose:
                    print '\t', name
                description_file_name = '{0}/job/description.txt'.format(d)
                experiments[kind][family][name] = Experiment.from_description_file_name(description_file_name)

    return experiments

def test_joint_plot():
    #gene_name = 'YBR068C'
    gene_name = 'YAL035W'
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
    ax.legend(_handles[ax], _labels[ax], framealpha=0.5)

if __name__ == '__main__1':
    #gene_name = 'YAL035W'
    gene_name = 'YIL127C'
    experiments = build_all_experiments(verbose=False)

    ribosome_profiling_experiments = [(n, e) for n, e in sorted(experiments['ribosome_profiling']['weinberg'].items()) if 'Dynabeads' in n]
    print ribosome_profiling_experiments
    
    for name, experiment in ribosome_profiling_experiments:
        fig, ax = plt.subplots()
        fn = experiment.file_names['three_prime_read_positions']
        f = h5py.File(fn, 'r')
        gene = Serialize.read_positions.build_gene(f[gene_name])

        xs = np.arange(-10, gene['all'].CDS_length + 100)
        counts = gene['all']['start_codon', xs]
        length_so_far = np.arange(len(xs)) + 1
        length_remaining = length_so_far[::-1]
        counts_so_far = counts.cumsum()
        counts_remaining = np.cumsum(counts[::-1])[::-1]
        density_so_far = np.true_divide(counts_so_far, length_so_far)
        density_remaining = np.true_divide(counts_remaining, length_remaining)

        n = 20
        density_over_previous = np.asarray([gene['all']['start_codon', x - 10:x].sum() / 10. for x in xs])

        ax.plot(xs, density_so_far, label='density so far')
        ax.plot(xs, density_remaining, label='density over remaining')
        ax.plot(xs, density_over_previous, label='density over last {0}'.format(n))
        ax.axhline(np.true_divide(counts.sum(), len(xs)), label='overall density', color='black')
        ax.legend()
        ax.set_title(name)
        ax.set_xlim(min(xs), max(xs))

        raw_counts_ax = ax.twinx()
        raw_counts_ax.plot(xs, counts, '.')

def plot_gene(gene_name):
    #gene_name = 'YAL035W'
    #gene_name = 'YMR005W'
    #gene_name = 'YIL127C'
    #gene_name = 'YKL064W'
    #gene_name = 'YCR012W'

    experiments = build_all_experiments(verbose=False)

    five_prime_experiments = [(n, e) for n, e in sorted(experiments['TL_seq']['arribere_gr'].items()) if 'TLSeq1' in n] + \
                             [(n, e) for n, e in sorted(experiments['TL_seq']['park_nar'].items()) if n == 'SMORE-seq_WT_TAP+_rep1'] + \
                             [(n, e) for n, e in sorted(experiments['TIF_seq']['pelechano_nature'].items()) if n == 'ypd_bio1_lib1' or n == 'ypd_bio1_lib4']
                             
    three_prime_experiments = [(n, e) for n, e in sorted(experiments['three_p_seq']['three_p_seq'].items())] + \
                              [(n, e) for n, e in sorted(experiments['TIF_seq']['pelechano_nature'].items()) if n == 'ypd_bio1_lib1' or n == 'ypd_bio1_lib4'] + \
                              [(n, e) for n, e in sorted(experiments['three_t_fill_seq']['wilkening_nar'].items()) if '3tfill_ypd_rep1' in n]

    rna_seq_experiments = [(n, e) for n, e in sorted(experiments['ribosome_profiling']['weinberg'].items()) if 'RPF' not in n]

    composition_fn = '/home/jah/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/transcript_recent_As.hdf5'
    composition_file = h5py.File(composition_fn, 'r')
    gene_composition = Serialize.read_positions.build_gene(composition_file[gene_name])

    fig, ax = plt.subplots()

    for name, experiment in five_prime_experiments:
        print name
        fn = experiment.file_names['read_positions']
        f = h5py.File(fn, 'r')
        gene = Serialize.read_positions.build_gene(f[gene_name])
        print gene.keys()
        xs = np.arange(-200, gene['all'].CDS_length)

        counts = gene['all']['start_codon', xs]
        #normalization = float(experiment.get_total_eligible_reads()) / 1.e6
        #normalization = 1.
        normalization = counts.sum()
        normalized = np.true_divide(counts, normalization)
        plot_with_line_alpha(ax, xs, normalized, label=name)

    composition_ax = ax.twinx()
    xs = np.arange(-200, gene_composition[20].CDS_length)
    plot_with_line_alpha(composition_ax,
                         xs,
                         np.true_divide(gene_composition[20]['start_codon', xs], 20),
                        )
    composition_ax.set_ylim(0, 1)
    
    fig, ax = plt.subplots()
    
    for name, experiment in three_prime_experiments:
        fn = experiment.file_names['read_positions']
        f = h5py.File(fn, 'r')
        gene = Serialize.read_positions.build_gene(f[gene_name])
        xs = np.arange(-gene.itervalues().next().CDS_length, 400)

        total_reads = experiment.get_total_eligible_reads()
        #total_reads = 1.e6
        print name, total_reads

        if 0 in gene:
        #if False:
            genomic_counts = gene[0]['stop_codon', xs]
            nongenomic_counts = gene['all']['stop_codon', xs] - gene[0]['stop_codon', xs]
            #normalization = float(experiment.get_total_eligible_reads()) / 1.e6
            #normalization = 1.
            #normalization = (genomic_counts + nongenomic_counts).sum()
            genomic_normalized = np.true_divide(genomic_counts, genomic_counts.sum())
            nongenomic_normalized = np.true_divide(nongenomic_counts, nongenomic_counts.sum())

            plot_with_line_alpha(ax, xs, genomic_normalized, label='{}_genomic'.format(name))
            plot_with_line_alpha(ax, xs, nongenomic_normalized, label='{}_nongenomic'.format(name))
        else:
            print gene.keys()
            counts = gene['all']['stop_codon', xs]
            #normalization = float(experiment.get_total_eligible_reads()) / 1.e6
            #normalization = 1.
            normalization = counts.sum()
            normalized = np.true_divide(counts, normalization)
            plot_with_line_alpha(ax, xs, normalized, label=name)
    
    #composition_ax = ax.twinx()
    #xs = np.arange(-gene_composition[20].CDS_length, 100)
    #plot_with_line_alpha(composition_ax,
    #                     xs,
    #                     np.true_divide(gene_composition[20]['stop_codon', xs], 20),
    #                    )
    #composition_ax.set_ylim(0, 1)

def produce_transcript_base_compositions():
    gff_fn = '/home/jah/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/transcriptome/genes.gff'
    genome_dir = '/home/jah/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/genome/'
    composition_fn = '/home/jah/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/transcript_recent_As.hdf5'
    CDSs = gff.get_CDSs(gff_fn, genome_dir)

    left_buffer = 500
    right_buffer = 500
    genes = {}

    windows = [5, 10, 20]

    for transcript in utilities.progress_bar(len(CDSs), CDSs):
        genes[transcript.name] = {}
        transcript.build_coordinate_maps()
        landmarks = {'start': 0,
                     'start_codon': transcript.transcript_start_codon,
                     'stop_codon': transcript.transcript_stop_codon,
                     'end': transcript.transcript_length,
                    }
        sequence = transcript.get_transcript_sequence(left_buffer, right_buffer)
        
        A_locations = positions.PositionCounts(landmarks,
                                               left_buffer,
                                               right_buffer,
                                               data=(sequence.data == 'A'),
                                              )
        for window in windows:
            recent_As = positions.PositionCounts(landmarks,
                                                 left_buffer,
                                                 right_buffer,
                                                )
            for left_edge in range(-left_buffer, transcript.CDS_length + right_buffer - window):
                num_As = sum(A_locations['start', left_edge:left_edge + window])
                recent_As['start', left_edge] = num_As

            genes[transcript.name][window] = recent_As

        transcript.delete_coordinate_maps()

    Serialize.read_positions.write_file(genes, composition_fn)

def call_3p_peaks():
    gtf_fn = '/home/jah/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/transcriptome/genes.gtf'
    genome_dir = '/home/jah/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/genome/'
    composition_fn = '/home/jah/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/transcript_recent_As.hdf5'

    output_fn = '/home/jah/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/transcript_3p_lengths.txt'

    region_fetcher = genomes.build_region_fetcher(genome_dir)
    CDSs = gtf.get_CDSs(gtf_fn)
    CDS_dict = {t.name: t for t in CDSs}
    
    experiments = build_all_experiments(verbose=False)
    
    three_prime_experiments = [(n, e) for n, e in sorted(experiments['three_p_seq']['three_p_seq'].items())] + \
                              [(n, e) for n, e in sorted(experiments['three_t_fill_seq']['wilkening_nar'].items()) if '3tfill_ypd_rep1' in n] + \
                              [(n, e) for n, e in sorted(experiments['TIF_seq']['pelechano_nature'].items()) if n == 'ypd_bio1_lib1' or n == 'ypd_bio1_lib4']
    
    argmaxes = {}
    fractions = {}
    joints = {}
    for name, experiment in three_prime_experiments:
        print name
        argmaxes[name] = {}
        fractions[name] = []
        joints[name] = []
        fn = experiment.file_names['three_prime_read_positions']
        f = h5py.File(fn, 'r')
        for transcript in utilities.progress_bar(len(CDSs), CDSs):
            if transcript.name not in f:
                continue
            gene = Serialize.read_positions.build_gene(f[transcript.name])
            xs = np.arange(0, 400)

            argmax = gene['all'].argmax_over_slice('stop_codon', xs)
            argmaxes[name][transcript.name] = argmax
            most = gene['all']['stop_codon', argmax]
            total = gene['all']['stop_codon', xs].sum()
            if total > 9:
                fraction = np.true_divide(most, total)
                fractions[name].append(fraction)
                joints[name].append((argmax, fraction))
    
    with open(output_fn, 'w') as output_fh:
        name_order = sorted(argmaxes['Cerevisiae_3Pseq'], key=argmaxes['Cerevisiae_3Pseq'].get)
        for name in name_order:
            output_fh.write('{0}\t'.format(str(CDS_dict[name])))
            for exp_name, _ in three_prime_experiments:
                output_fh.write('{0}\t'.format(argmaxes[exp_name][name]))
            output_fh.write('\n')

def call_5p_peaks():
    gtf_fn = '/home/jah/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/transcriptome/genes.gtf'
    genome_dir = '/home/jah/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/genome/'
    region_fetcher = genomes.build_region_fetcher(genome_dir)
    CDSs = gtf.get_CDSs(gtf_fn)
    
    experiments = build_all_experiments(verbose=False)
    
    five_prime_experiments = [(n, e) for n, e in sorted(experiments['TL_seq']['arribere_gr'].items()) if 'TLSeq1' in n] + \
                             [(n, e) for n, e in sorted(experiments['TL_seq']['park_nar'].items()) if n == 'SMORE-seq_WT_TAP+_rep1'] + \
                             [(n, e) for n, e in sorted(experiments['TIF_seq']['pelechano_nature'].items()) if n == 'ypd_bio1_lib1' or n == 'ypd_bio1_lib4']

    argmaxes = {}
    fractions = {}
    joints = {}
    for name, experiment in five_prime_experiments:
        print name
        argmaxes[name] = Counter()
        fractions[name] = []
        joints[name] = []
        fn = experiment.file_names['five_prime_read_positions']
        f = h5py.File(fn, 'r')
        for transcript in utilities.progress_bar(len(CDSs), CDSs):
            if transcript.name not in f:
                continue
            gene = Serialize.read_positions.build_gene(f[transcript.name])
            xs = np.arange(-300, 0)

            argmax = gene['all'].argmax_over_slice('start_codon', xs)
            argmaxes[name][argmax] += 1
            most = gene['all']['start_codon', argmax]
            total = gene['all']['start_codon', xs].sum()
            if total == 0:
                print transcript
            if total > 9:
                fraction = np.true_divide(most, total)
                fractions[name].append(fraction)
                joints[name].append((argmax, fraction))

if __name__ == '__main__':
    gff_fn = '/home/jah/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/transcriptome/genes.gff'
    genome_dir = '/home/jah/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/genome/'
    composition_fn = '/home/jah/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/transcript_recent_As.hdf5'
    CDSs = gff.get_CDSs(gff_fn, genome_dir)

    import select_work
    exps = select_work.build_all_experiments(verbose=False)

    reads_fn = exps['belgium_2014_12_10']['WT_1_mRNA'].file_names['three_prime_read_positions']
    reads_fh = h5py.File(reads_fn, 'r')

    meta_counts = positions.PositionCounts({'A':0}, left_buffer=100000, right_buffer=100000)

    f = h5py.File(composition_fn, 'r')
    for t in utilities.progress_bar(len(CDSs), CDSs):
        if t.name not in reads_fh:
            continue
        gene = Serialize.read_positions.build_gene(f[t.name])
        t.build_coordinate_maps()

        if t.transcript_length < 301:
            continue
        end = t.transcript_length - 200
        sl = ('start', np.arange(100, end))
        A_rich_position = gene[10].argmax_over_slice(*sl)
        if gene[10]['start', A_rich_position] > 9:
            counts = Serialize.read_positions.build_gene(reads_fh[t.name])
            before_counts = counts['all']['start', 0:A_rich_position]
            after_counts = counts['all']['start', A_rich_position:A_rich_position + 200]
            meta_counts['A', -len(before_counts):0] += before_counts
            meta_counts['A', :len(after_counts)] += after_counts
