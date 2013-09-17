#import matplotlib
#matplotlib.use('Agg', warn=False)
import trim
from collections import Counter
import glob
import os.path
import mutations
import numpy as np
import time
import pysam
import ribosomes
import matplotlib.pyplot as plt
import subprocess

def make_file_names(group, sample):
    data_dir = '/home/jah/projects/arlen/data/{0}/{1}/'.format(group, sample)
    results_dir = '/home/jah/projects/arlen/results/{0}/{1}/'.format(group, sample)
    igenomes_dir = '/home/jah/projects/arlen/data/igenomes/Saccharomyces_cerevisiae/Ensembl/EF2/'
    tophat_dir = results_dir + '/tophat/'
    
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)

    file_names = {
        'results': results_dir,
        'data': data_dir,
        'tophat': tophat_dir,
        'igenomes': igenomes_dir,
        'genome': igenomes_dir + '/Sequence/Chromosomes',
        'R1_reads': glob.glob(data_dir + '*.fastq'),
        'gtf': igenomes_dir + 'Annotation/Genes/genes.gtf',
        'bowtie2_index': igenomes_dir + 'Sequence/Bowtie2Index/genome',
        'trimmed_reads':     results_dir + 'trimmed.fastq',
        'trimmed_lengths':   results_dir + 'trimmed_lengths.txt',
        'adapter_distances': results_dir + 'adapter_distances.txt',
        'no_rRNA_reads':     results_dir + 'no_rRNA.fastq',
        'rRNA_mapping_log':  results_dir + 'rRNA_mapping_log.txt',
        'rRNA_sam':          results_dir + 'rRNA.sam',
        'no_rRNA_lengths':   results_dir + 'no_rRNA_lengths.txt',
        'filtered_bam':      results_dir + 'filtered.bam',
        'unambiguous_bam':   results_dir + 'unambiguous.bam',
        'tRNA_lengths':      results_dir + 'tRNA_lengths.txt',
        'more_rRNA_lengths': results_dir + 'more_rRNA_lengths.txt',
        'from_ends':         results_dir + 'from_ends.txt',
        'from_starts':       results_dir + 'from_starts.txt',
        'from_ends_unambiguous':         results_dir + 'from_ends_unambiguous.txt',
        'from_starts_unambiguous':       results_dir + 'from_starts_unambiguous.txt',
        'fractions':         results_dir + 'fraction.txt',
        'summary':           results_dir + 'summary.txt',
        'lengths_fig':       results_dir + 'lengths.png',
        'accepted_hits':     tophat_dir + 'accepted_hits.bam',
        'unmapped':          tophat_dir + 'unmapped.bam',
    }

    return file_names

def get_all_lengths(bam_fn):
    bamfile = pysam.Samfile(bam_fn, 'rb')
    qlens = Counter(ar.qlen for ar in bamfile)
    return mutations.counts_to_array(qlens)

def plot_lengths(file_names):
    trimmed_lengths = np.loadtxt(file_names['trimmed_lengths'])
    no_rRNA_lengths = np.loadtxt(file_names['no_rRNA_lengths'])
    unmapped_lengths = get_all_lengths(file_names['unmapped'])
    tRNA_lengths = np.loadtxt(file_names['tRNA_lengths'])
    more_rRNA_lengths = np.loadtxt(file_names['more_rRNA_lengths'])

    rRNA_lengths = trimmed_lengths - no_rRNA_lengths + more_rRNA_lengths
    rpf_lengths = no_rRNA_lengths - more_rRNA_lengths - tRNA_lengths - unmapped_lengths

    fig, (ax_all, ax_rpf) = plt.subplots(2, 1, figsize=(12, 16))
    
    ax_all.plot(rRNA_lengths, '.-', label='rRNA', color='green')
    ax_all.plot(tRNA_lengths, '.-', label='tRNA', color='blue')
    ax_all.plot(rpf_lengths, '.-', label='RPF', color='red')
    ax_all.legend()

    ax_rpf.plot(rpf_lengths, '.-', label='RFP', color='red')
    ax_rpf.legend()

    fig.savefig(file_names['lengths_fig'])
    plt.close()

samples = [
    #('belgium_8_6_13', 'R98S_cDNA_mRNA', trim.trim_adapters),
    #('belgium_8_6_13', 'Suppressed_R98S_cDNA_mRNA', trim.trim_adapters),
    #('belgium_8_6_13', 'WT_cDNA_mRNA', trim.trim_adapters),
    ('belgium_8_6_13', 'R98S_cDNA_sample', trim.trim_adapters),
    ('belgium_8_6_13', 'Suppressed_R98S_cDNA_sample', trim.trim_adapters),
    ('belgium_8_6_13', 'WT_cDNA_sample', trim.trim_adapters),
    #('ingolia_science', 'mRNA-rich-1', trim.trim_poly_A),
    #('ingolia_science', 'mRNA-rich-2', trim.trim_poly_A),
    ('ingolia_science', 'Footprints-rich-1', trim.trim_poly_A),
    ('ingolia_science', 'Footprints-rich-2', trim.trim_poly_A),
    #('nagalakshmi_science', 'dT_ori', trim.trim_nothing),
    #('nagalakshmi_science', 'dT_tech', trim.trim_nothing),
    #('nagalakshmi_science', 'RH_ori', trim.trim_nothing),
]

def make_plots():
    file_names_list = [make_file_names(group, sample) for group, sample, _ in samples]

    from_starts_fns = [fns['from_starts'] for fns in file_names_list]
    from_ends_fns = [fns['from_ends'] for fns in file_names_list]
    fractions_fns = [fns['fractions'] for fns in file_names_list]
    read_lengths = [trim.get_read_length(fns['R1_reads']) for fns in file_names_list]
    names = [name for group, name, _ in samples]
    summary_fns = [fns['summary'] for fns in file_names_list]

    #from_starts_fns.extend([fns['from_starts_unambiguous'] for fns in file_names_list])
    #from_ends_fns.extend([fns['from_ends_unambiguous'] for fns in file_names_list])
    #fractions_fns.extend([fns['fractions'] for fns in file_names_list])
    #read_lengths.extend([trim.get_read_length(fns['R1_reads']) for fns in file_names_list])
    #names.extend(['{0}_unambiguous'.format(name) for group, name, _ in samples])

    #ribosomes.plot_RPF(from_starts_fns,
    #                   from_ends_fns,
    #                   names,
    #                   read_lengths,
    #                  )
    
    #fragment_lengths = slice(20, None)
    #ribosomes.plot_mRNA(from_starts_fns,
    #                    from_ends_fns,
    #                    fractions_fns,
    #                    names,
    #                    read_lengths,
    #                    fragment_lengths,
    #                   )

    ribosomes.process_summary(summary_fns, names)

def make_data():
    for group, sample, trim_function in samples:
        print group, sample
        file_names = make_file_names(group, sample)
        
        read_length = trim.get_read_length(file_names['R1_reads'])

        start_time = time.time()
        
        #trim_function(file_names['R1_reads'],
        #              file_names['trimmed_reads'],
        #              file_names['trimmed_lengths'],
        #             )

        #ribosomes.filter_contaminant('yeast_rRNA',
        #                             file_names['trimmed_reads'],
        #                             file_names['no_rRNA_reads'],
        #                             file_names['rRNA_sam'],
        #                             file_names['rRNA_mapping_log'],
        #                             file_names['no_rRNA_lengths'],
        #                            )

        #ribosomes.map_tophat(file_names)

        #ribosomes.further_filtering(file_names['accepted_hits'],
        #                            file_names['gtf'],
        #                            file_names['filtered_bam'],
        #                            file_names['tRNA_lengths'],
        #                            file_names['more_rRNA_lengths'],
        #                            read_length,
        #                           )

        trim.unambiguously_trimmed(file_names['filtered_bam'],
                                   file_names['unambiguous_bam'],
                                   file_names['genome'],
                                  )

        #plot_lengths(file_names)

        #from_starts, from_ends, fractions = ribosomes.get_aggregate_positions_fetch(file_names['gtf'],
        #                                                                            file_names['filtered_bam'],
        #                                                                            file_names['summary'],
        #                                                                            read_length,
        #                                                                           )
        #np.savetxt(file_names['from_starts'], from_starts, fmt='%d')
        #np.savetxt(file_names['from_ends'], from_ends, fmt='%d')
        #np.savetxt(file_names['fractions'], fractions, fmt='%d')

        from_starts, from_ends, fractions = ribosomes.get_aggregate_positions_fetch(file_names['gtf'],
                                                                                    file_names['unambiguous_bam'],
                                                                                    file_names['summary'],
                                                                                    read_length,
                                                                                   )
        np.savetxt(file_names['from_starts_unambiguous'], from_starts, fmt='%d')
        np.savetxt(file_names['from_ends_unambiguous'], from_ends, fmt='%d')

        end_time = time.time()
        print 'Processed in {0:0.2f} s'.format(end_time - start_time)

if __name__ == '__main__':
    #make_data()
    make_plots()
