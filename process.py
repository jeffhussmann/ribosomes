import matplotlib
matplotlib.use('Agg', warn=False)
import trim
from collections import Counter
import glob
import os.path
import mutations
import numpy as np
import time
import pysam
import sam
import ribosomes
import matplotlib.pyplot as plt
import subprocess

def make_file_names(group, sample):
    data_dir = '/home/jah/projects/arlen/data/{0}/{1}/'.format(group, sample)
    results_dir = '/home/jah/projects/arlen/results/{0}/{1}/'.format(group, sample)
    organism_dir = '/home/jah/projects/arlen/data/organisms/saccharomyces_cerevisiae/'
    tophat_dir = results_dir + '/tophat/'
    oligos_sam_fn = '/home/jah/projects/arlen/data/organisms/saccharomyces_cerevisiae/contaminant/subtraction_oligos.sam'
    
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)

    file_names = {
        'results': results_dir,
        'data': data_dir,
        'tophat': tophat_dir,
        'genome': organism_dir + '/genome',
        'R1_reads': glob.glob(data_dir + '*.fastq'),
        'gtf': organism_dir + 'Annotation/Genes/genes.gtf',
        'bowtie2_index': organism_dir + 'genome',
        'trimmed_reads':     results_dir + 'trimmed.fastq',

        'oligos_sam': oligos_sam_fn,

        'trimmed_lengths':   results_dir + 'trimmed_lengths.txt',
        'filtered_lengths':  results_dir + 'filtered_lengths.txt',
        'tRNA_lengths':      results_dir + 'tRNA_lengths.txt',
        'rRNA_lengths': results_dir + 'rRNA_lengths.txt',
        'clean_lengths':     results_dir + 'clean_lengths.txt',
        'unambiguous_lengths':     results_dir + 'unambiguous_lengths.txt',

        'filtered_reads':     results_dir + 'filtered.fastq',
        'rRNA_mapping_log':  results_dir + 'rRNA_mapping_log.txt',
        'rRNA_sam':          results_dir + 'rRNA.sam',
        'rRNA_bam':          results_dir + 'rRNA.bam',
        'rRNA_coverage':     results_dir + 'rRNA_coverage.txt',
        'clean_bam':         results_dir + 'clean.bam',
        'unambiguous_bam':   results_dir + 'unambiguous.bam',
        'from_ends':         results_dir + 'from_ends.txt',
        'from_starts':       results_dir + 'from_starts.txt',
        'from_ends_unambiguous':         results_dir + 'from_ends_unambiguous.txt',
        'from_starts_unambiguous':       results_dir + 'from_starts_unambiguous.txt',
        'fractions':         results_dir + 'fraction.txt',
        'summary':           results_dir + 'summary.txt',

        'all_lengths_fig':       results_dir + 'all_lengths.png',
        'clean_lengths_fig':       results_dir + 'clean_lengths.png',
        'accepted_hits_bam':     tophat_dir + 'accepted_hits.bam',
        'unmapped_bam':          tophat_dir + 'unmapped.bam',
        'oligo_lengths_fig': results_dir + 'oligo_lengths.png',

        'log':              results_dir + 'log.txt',
    }

    return file_names

def get_all_lengths(bam_fn):
    bamfile = pysam.Samfile(bam_fn, 'rb')
    qlens = Counter(ar.qlen for ar in bamfile)
    return mutations.counts_to_array(qlens)

def plot_lengths(file_names, unambiguous):
    trimmed_lengths = np.loadtxt(file_names['trimmed_lengths'])
    filtered_lengths = np.loadtxt(file_names['filtered_lengths'])
    unmapped_lengths = get_all_lengths(file_names['unmapped_bam'])
    tRNA_lengths = np.loadtxt(file_names['tRNA_lengths'])
    rRNA_lengths = np.loadtxt(file_names['rRNA_lengths'])
    clean_lengths = get_all_lengths(file_names['clean_bam'])
    np.savetxt(file_names['clean_lengths'], clean_lengths, fmt='%d')

    fig_all, ax_all = plt.subplots(figsize=(12, 8))
    
    ax_all.plot(rRNA_lengths, '.-', label='rRNA', color='red')
    ax_all.plot(tRNA_lengths, '.-', label='tRNA', color='blue')
    ax_all.plot(clean_lengths, '.-', label='clean', color='green')
    ax_all.plot(unmapped_lengths, '.-', label='unmapped', color='c')
    ax_all.set_xlim(0, 50)
    ax_all.set_title('Fragment length distribution by source')
    ax_all.set_xlabel('Length of original RNA fragment')
    ax_all.set_ylabel('Number of reads')

    fig_clean, ax_clean = plt.subplots(figsize=(12, 8))
    
    ax_clean.plot(clean_lengths, '.-', label='clean', color='green')
    ax_clean.axvspan(27.5, 28.5, color='green', alpha=0.2)
    ax_clean.set_xlim(0, 50)
    ax_clean.legend()
    ax_clean.set_title('Fragment length distribution by source')
    ax_clean.set_xlabel('Length of original RNA fragment')
    ax_clean.set_ylabel('Number of reads')
    
    if unambiguous:
        unambiguous_lengths = get_all_lengths(file_names['unambiguous_bam'])
        np.savetxt(file_names['unambiguous_lengths'], unambiguous_lengths, fmt='%d')
        ax_all.plot(unambiguous_lengths, '.-', label='clean, unambiguous', color='purple')
        ax_clean.plot(unambiguous_lengths, '.-', label='clean, unambiguous', color='purple')

    ax_all.legend()
    fig_all.savefig(file_names['all_lengths_fig'])
    fig_clean.savefig(file_names['clean_lengths_fig'])
    
    total_reads = int(trimmed_lengths.sum())
    rRNA_reads = int(rRNA_lengths.sum())
    tRNA_reads = int(tRNA_lengths.sum())
    unmapped_reads = int(unmapped_lengths.sum())
    clean_reads = int(clean_lengths.sum())

    with open(file_names['log'], 'w') as log_fh:
        log_fh.write('Total reads: {0:,}\n'.format(total_reads))
        log_fh.write('rRNA reads: {0:,} ({1:0.2f}%)\n'.format(rRNA_reads, 100 * float(rRNA_reads) / total_reads))
        log_fh.write('tRNA reads: {0:,} ({1:0.2f}%)\n'.format(tRNA_reads, 100 * float(tRNA_reads) / total_reads))
        log_fh.write('Unmapped reads: {0:,} ({1:0.2f}%)\n'.format(unmapped_reads, 100 * float(unmapped_reads) / total_reads))
        log_fh.write('Clean reads: {0:,} ({1:0.2f}%)\n'.format(clean_reads, 100 * float(clean_reads) / total_reads))

samples = [
    #('belgium_8_6_13', 'R98S_cDNA_mRNA', trim.trim_adapters, False),
    #('belgium_8_6_13', 'Suppressed_R98S_cDNA_mRNA', trim.trim_adapters, False),
    #('belgium_8_6_13', 'WT_cDNA_mRNA', trim.trim_adapters, False),
    ('belgium_8_6_13', 'R98S_cDNA_sample', trim.trim_adapters, False),
    ('belgium_8_6_13', 'Suppressed_R98S_cDNA_sample', trim.trim_adapters, False),
    ('belgium_8_6_13', 'WT_cDNA_sample', trim.trim_adapters, False),
    #('ingolia_science', 'mRNA-rich-1', trim.trim_poly_A, True),
    #('ingolia_science', 'mRNA-rich-2', trim.trim_poly_A, True),
    ('ingolia_science', 'Footprints-rich-1', trim.trim_poly_A, True),
    #('ingolia_science', 'Footprints-rich-2', trim.trim_poly_A, True),
    #('nagalakshmi_science', 'dT_ori', trim.trim_nothing, False),
    #('nagalakshmi_science', 'dT_tech', trim.trim_nothing, False),
    #('nagalakshmi_science', 'RH_ori', trim.trim_nothing, False),
    ('brar_science', 's_tA-fp_100211_l3_sequence', trim.trim_poly_A, True),
]

def make_plots():
    file_names_list = [make_file_names(group, sample) for group, sample, _, _ in samples]
    unambs = [u for _, _, _, u in samples]

    from_starts_fns = [fns['from_starts'] for fns in file_names_list]
    from_ends_fns = [fns['from_ends'] for fns in file_names_list]
    fractions_fns = [fns['fractions'] for fns in file_names_list]
    trimmed_lengths_fns = [fns['trimmed_lengths'] for fns in file_names_list]
    clean_length_fns = [fns['clean_lengths'] if not unamb else fns['unambiguous_lengths']
                        for fns, unamb in zip(file_names_list, unambs)]
    read_lengths = [trim.get_read_length(fns['R1_reads']) for fns in file_names_list]
    names = [name for group, name, _, _ in samples]
    summary_fns = [fns['summary'] for fns in file_names_list]

    #ribosomes.plot_trimmed_lengths(trimmed_lengths_fns, names)

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

    #ribosomes.plot_frames(summary_fns, names)

    #ribosomes.plot_clean_lengths(clean_length_fns, names)

    ribosomes.compare(summary_fns, clean_length_fns, names)

def make_data():
    for group, sample, trim_function, unambiguous in samples:
        start_time = time.time()
        
        print group, sample
        file_names = make_file_names(group, sample)
        read_length = trim.get_read_length(file_names['R1_reads'])
        
        #trim_function(file_names['R1_reads'],
        #              file_names['trimmed_reads'],
        #              file_names['trimmed_lengths'],
        #             )

        ribosomes.pre_filter('yeast_rRNA',
                             file_names['trimmed_reads'],
                             file_names['filtered_reads'],
                             file_names['rRNA_sam'],
                             file_names['rRNA_mapping_log'],
                             file_names['filtered_lengths'],
                            )

        sam.sam_to_sorted_indexed_bam(file_names['rRNA_sam'],
                                      file_names['rRNA_bam'],
                                     )

        ribosomes.produce_rRNA_coverage(file_names['rRNA_bam'],
                                        file_names['oligos_sam'],
                                        file_names['rRNA_coverage'],
                                       )
        
        #trimmed_lengths = np.loadtxt(file_names['trimmed_lengths'])
        #filtered_lengths = np.loadtxt(file_names['filtered_lengths'])
        #pre_rRNA_lengths = trimmed_lengths - filtered_lengths

        #ribosomes.map_tophat(file_names)

        #tRNA_lengths, rRNA_lengths = ribosomes.filter(file_names['accepted_hits_bam'],
        #                                              file_names['gtf'],
        #                                              file_names['clean_bam'],
        #                                              read_length,
        #                                             )

        #rRNA_lengths += pre_rRNA_lengths
        #np.savetxt(file_names['rRNA_lengths'], rRNA_lengths, fmt='%d')
        #np.savetxt(file_names['tRNA_lengths'], tRNA_lengths, fmt='%d')
        #
        #if unambiguous:
        #    trim.unambiguously_trimmed(file_names['clean_bam'],
        #                               file_names['unambiguous_bam'],
        #                               file_names['genome'],
        #                              )

        #plot_lengths(file_names, unambiguous)

        #if not unambiguous:
        #    from_starts, from_ends, fractions = ribosomes.get_aggregate_positions_fetch(file_names['gtf'],
        #                                                                                file_names['clean_bam'],
        #                                                                                file_names['summary'],
        #                                                                                read_length,
        #                                                                               )
        #    np.savetxt(file_names['from_starts'], from_starts, fmt='%d')
        #    np.savetxt(file_names['from_ends'], from_ends, fmt='%d')
        #    np.savetxt(file_names['fractions'], fractions, fmt='%d')
        #else:
        #    from_starts, from_ends, fractions = ribosomes.get_aggregate_positions_fetch(file_names['gtf'],
        #                                                                                file_names['unambiguous_bam'],
        #                                                                                file_names['summary'],
        #                                                                                read_length,
        #                                                                               )
        #    np.savetxt(file_names['from_starts'], from_starts, fmt='%d')
        #    np.savetxt(file_names['from_ends'], from_ends, fmt='%d')

        end_time = time.time()
        print 'Processed in {0:0.2f} s'.format(end_time - start_time)


def plot_coverage():
    file_names_list = [make_file_names(group, sample) for group, sample, _, _ in samples]
    #coverage_fns = [fns['rRNA_coverage'] for fns in file_names_list]
    names = [sample for _, sample, _, _ in samples]
    #ribosomes.plot_rRNA_coverage(names, coverage_fns, file_names_list[0]['oligos_sam'])

    for file_names, name in zip(file_names_list, names):
        ribosomes.plot_oligo_hit_lengths(file_names['rRNA_bam'],
                                         file_names['oligos_sam'],
                                         file_names['oligo_lengths_fig'],
                                         name,
                                        )

if __name__ == '__main__':
    #make_data()
    #make_plots()
    plot_coverage()
