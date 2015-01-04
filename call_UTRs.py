import numpy as np
import Serialize.read_positions
import h5py
import Sequencing.utilities as utilities
import matplotlib.pyplot as plt
import gff
import transcript

def call_UTR_boundaries(boundaries_fn):
    import three_p_experiment
    import TL_seq_experiment
    three_prime_description_fn = '/home/jah/projects/ribosomes/experiments/three_p_seq/Cerevisiae_3Pseq/job/description.txt'
    three_prime_exp = three_p_experiment.ThreePExperiment.from_description_file_name(three_prime_description_fn)


    three_prime_genes = Serialize.read_positions.read_file(three_prime_exp.file_names['three_prime_read_positions'],
                                                           specific_keys={'all', '0'},
                                                          )
    
    five_prime_description_fn = '/home/jah/projects/ribosomes/experiments/arribere_gr/S288C_TLSeq1/job/description.txt'
    five_prime_exp = TL_seq_experiment.TLSeqExperiment.from_description_file_name(five_prime_description_fn)

    five_prime_genes = Serialize.read_positions.read_file(five_prime_exp.file_names['five_prime_read_positions'],
                                                          specific_keys={'all'},
                                                         )

    three_xs = np.arange(3, 500)
    five_xs = np.arange(-200, 0)

    transcripts = gff.get_CDSs(three_prime_exp.file_names['genes'],
                               three_prime_exp.file_names['genome'],
                              )
    transcripts = {t.name: t for t in transcripts}

    UTR_boundaries = {}
    
    for name in three_prime_genes:
        t = transcripts[name]
        t.build_coordinate_maps(left_buffer=500, right_buffer=500)

        five_counts = five_prime_genes[name]['all']
        five_count = five_counts['start_codon', five_xs].sum()
        if five_count == 0:
            five_offset = 0
        else:
            five_offset = five_counts.argmax_over_slice('start_codon', five_xs)

        three_counts = three_prime_genes[name]['all'] - three_prime_genes[name][0]
        three_count = three_counts['stop_codon', three_xs].sum()
        if three_count == 0:
            three_offset = 3
        else:
            three_offset = three_counts.argmax_over_slice('stop_codon', three_xs)

        five_pos = t.transcript_to_genomic[t.transcript_start_codon + five_offset]
        three_pos = t.transcript_to_genomic[t.transcript_stop_codon + three_offset]
        
        t.delete_coordinate_maps()

        UTR_boundaries[name] = (t.seqname, t.strand, five_pos, three_pos)

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
