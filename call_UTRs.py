import numpy as np
import Serialize.read_positions
import h5py
import Sequencing.utilities as utilities
import matplotlib.pyplot as plt
import gff
import transcript
import glob
import os

def build_all_experiments(verbose=True):
    import three_p_experiment
    import TL_seq_experiment
    import TIF_seq_experiment
    import ribosome_profiling_experiment
    import three_t_fill_experiment

    families = {'TIF_seq': (TIF_seq_experiment.TIFSeqExperiment, ['pelechano_nature']),
                'three_p_seq': (three_p_experiment.ThreePExperiment, ['three_p_seq']),
                'TL_seq': (TL_seq_experiment.TLSeqExperiment, ['arribere_gr', 'park_nar']),
                'three_t_fill_seq': (three_t_fill_experiment.ThreeTFillExperiment, ['wilkening_nar']),
                'ribosome_profiling': (ribosome_profiling_experiment.RibosomeProfilingExperiment, ['weinberg']),
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

def call_UTR_boundaries(boundaries_fn, diagnostic_fn='/dev/null'):
    experiments = build_all_experiments(verbose=False)
    
    five_prime_exp = experiments['TL_seq']['arribere_gr']['S288C_TLSeq1']
    three_prime_exp = experiments['three_p_seq']['three_p_seq']['Cerevisiae_3Pseq']

    other_five_prime_exps = [experiments['TL_seq']['park_nar']['SMORE-seq_WT_TAP+_rep1'],
                             experiments['TIF_seq']['pelechano_nature']['ypd_bio1_lib1'],
                            ]

    
    other_three_prime_exps = [experiments['three_t_fill_seq']['wilkening_nar']['3tfill_ypd_rep1'],
                              experiments['TIF_seq']['pelechano_nature']['ypd_bio1_lib1'],
                             ]

    five_prime_fh = h5py.File(five_prime_exp.file_names['five_prime_read_positions'], 'r')
    three_prime_fh = h5py.File(three_prime_exp.file_names['three_prime_read_positions'], 'r')
    
    other_five_prime_fhs = [h5py.File(exp.file_names['five_prime_read_positions'], 'r') for exp in other_five_prime_exps]
    other_three_prime_fhs = [h5py.File(exp.file_names['three_prime_read_positions'], 'r') for exp in other_three_prime_exps]

    transcripts, _ = five_prime_exp.get_CDSs()

    UTR_boundaries = {}
    
    with open(diagnostic_fn, 'w') as diagnostic_fh:
        progress = utilities.progress_bar(len(transcripts), sorted(transcripts))
        for transcript in progress:
            name = transcript.name

            transcript.build_coordinate_maps(left_buffer=500, right_buffer=500)

            five_prime_gene = Serialize.read_positions.build_gene(five_prime_fh[name], specific_keys={'all'})
            other_genes = [Serialize.read_positions.build_gene(other_fh[name], specific_keys={'all'}) for other_fh in other_five_prime_fhs]
            five_xs = np.arange(-500, transcript.CDS_length)
            five_slice = ('start_codon', five_xs)
            
            five_counts = five_prime_gene['all']
            five_sum = five_counts[five_slice].sum()
            if five_sum == 0:
                five_offset = 0
            else:
                five_offset = five_counts.argmax_over_slice('start_codon', five_xs)

            n_largest = five_counts.n_largest_over_slice(10, five_slice)
            five_prime_diagnostic = []
            for i in n_largest:
                row = []
                for gene in [five_prime_gene] + other_genes:
                    count = gene['all']['start_codon', i]
                    total = gene['all'][five_slice].sum()
                    if row == []:
                        genomic = transcript.transcript_to_genomic[transcript.transcript_start_codon + i]
                        row.append('{0}\t({1:,})\t'.format(i, genomic))
                    row.append('{0}\t{1:0.2%}'.format(count, count / float(total)))
                five_prime_diagnostic.append('\t'.join(row))
            five_prime_diagnostic = '\n'.join(five_prime_diagnostic)
            
            three_prime_gene = Serialize.read_positions.build_gene(three_prime_fh[name], specific_keys={'all', '0'})
            other_genes = [Serialize.read_positions.build_gene(other_fh[name], specific_keys={'all', '0'}) for other_fh in other_three_prime_fhs]
            three_xs = np.arange(-transcript.CDS_length, 500)
            three_slice = ('stop_codon', three_xs)
            
            three_counts = three_prime_gene['all']# - three_prime_gene[0]
            three_sum = three_counts[three_slice].sum()
            if three_sum == 0:
                three_offset = 3
            else:
                three_offset = three_counts.argmax_over_slice('stop_codon', three_xs)

            n_largest = three_counts.n_largest_over_slice(10, three_slice)
            three_prime_diagnostic = []
            for i in n_largest:
                row = []
                for gene in [three_prime_gene] + other_genes:
                    count = gene['all']['stop_codon', i]
                    total = gene['all'][three_slice].sum()
                    if row == []:
                        genomic = transcript.transcript_to_genomic[transcript.transcript_stop_codon + i]
                        row.append('{0}\t({1:,})\t'.format(i, genomic))
                    row.append('{0}\t{1:0.2%}'.format(count, count / float(total)))
                three_prime_diagnostic.append('\t'.join(row))
            three_prime_diagnostic = '\n'.join(three_prime_diagnostic)

            diagnostic_fh.write('{0}\n'.format(str(transcript)))
            diagnostic_fh.write('{0}\n'.format(five_prime_diagnostic))
            diagnostic_fh.write('\n')
            diagnostic_fh.write('{0}\n'.format(three_prime_diagnostic))
            diagnostic_fh.write('\n')

            five_pos = transcript.transcript_to_genomic[transcript.transcript_start_codon + five_offset]
            three_pos = transcript.transcript_to_genomic[transcript.transcript_stop_codon + three_offset]
            
            transcript.delete_coordinate_maps()

            UTR_boundaries[name] = (transcript.seqname, transcript.strand, five_pos, three_pos)

    write_UTR_file(UTR_boundaries, boundaries_fn)

def look_at_densities():
    import ribosome_profiling_experiment
    description_fn = '/home/jah/projects/ribosomes/experiments/weinberg/RiboZero/job/description.txt'
    exp = ribosome_profiling_experiment.RibosomeProfilingExperiment.from_description_file_name(description_fn)
    
    names = []
    zero_ratios = []

    hdf5_file = h5py.File(exp.file_names['three_prime_read_positions'], 'r')
    progress = utilities.progress_bar(len(hdf5_file), hdf5_file)
    for gene_name in progress:
        gene = Serialize.read_positions.build_gene(hdf5_file[gene_name], specific_keys={'0'})
        zero_counts = gene[0]
        before = zero_counts['polyA', -100:1].sum()
        after = zero_counts['polyA', 1:102].sum()
        names.append(gene_name)
        zero_ratios.append((before, after))

    return names, zero_ratios

def write_UTR_file(UTR_boundaries, UTR_fn):
    def sort_key(name):
        seqname, strand, five_pos, three_pos = UTR_boundaries[name]
        return (seqname, min(five_pos, three_pos), max(five_pos, three_pos), strand)

    with open(UTR_fn, 'w') as UTR_fh:
        for name in sorted(UTR_boundaries, key=sort_key):
            seqname, strand, five_pos, three_pos = UTR_boundaries[name]
            line = '{0}\t{1}\t{2}\t{3}\t{4}\n'.format(name,
                                                      seqname,
                                                      strand,
                                                      five_pos,
                                                      three_pos,
                                                     )
            UTR_fh.write(line)

def read_UTR_file(UTR_fn):
    UTR_boundaries = {}
    for line in open(UTR_fn):
        name, seqname, strand, five_pos, three_pos = line.strip().split()
        UTR_boundaries[name] = (seqname, strand, int(five_pos), int(three_pos))

    return UTR_boundaries
