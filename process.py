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

def make_file_names(sample):
    data_dir = '/home/jah/projects/arlen/data'
    results_dir = '/home/jah/projects/arlen/results'
    igenomes_dir = data_dir + '/igenomes/Saccharomyces_cerevisiae/Ensembl/EF2'
    
    if not os.path.isdir('{0}/{1}/'.format(results_dir, sample)):
        os.makedirs('{0}/{1}/'.format(results_dir, sample))

    data_base = '{0}/{1}/{1}'.format(data_dir, sample)
    results_base = '{0}/{1}/{1}'.format(results_dir, sample)
    tophat_base = '{0}/{1}/tophat/'.format(results_dir, sample)

    file_names = {
        'results': results_dir,
        'data': data_dir,
        'tophat': tophat_base,
        'igenomes': igenomes_dir,
        'R1_reads': glob.glob(data_base + '*.fastq')[0],
        'gtf': igenomes_dir + '/Annotation/Genes/genes.gtf',
        'bowtie2_index': igenomes_dir + '/Sequence/Bowtie2Index/genome',
    }

    tails = {
        'trimmed_reads':     '_trimmed.fastq',
        'trimmed_lengths':   '_trimmed_lengths.txt',
        'adapter_distances': '_adapter_distances.txt',
        'no_rRNA_reads':     '_no_rRNA.fastq',
        'rRNA_mapping_log':  '_rRNA_mapping_log.txt',
        'rRNA_sam':          '_rRNA.sam',
        'no_rRNA_lengths':   '_no_rRNA_lengths.txt',
        'filtered_bam':      '_filtered.bam',
        'tRNA_lengths':      '_tRNA_lengths.txt',
        'more_rRNA_lengths': '_more_rRNA_lengths.txt',
        'from_ends':         '_from_ends.txt',
        'from_starts':       '_from_starts.txt',
        'lengths_fig':       '_lengths.png',
    }

    tophat_tails = {
        'accepted_hits': 'accepted_hits.bam',
        'unmapped': 'unmapped.bam',
    }

    results = {name: results_base + tail for name, tail in tails.items()}
    tophat_results = {name: tophat_base + tail for name, tail in tophat_tails.items()}
    file_names.update(results)
    file_names.update(tophat_results)
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


if __name__ == '__main__':
    samples = [
        'R98S_cDNA_mRNA',
        'R98S_cDNA_sample',
        #'Suppressed_R98S_cDNA_mRNA',
        #'Suppressed_R98S_cDNA_sample',
        'WT_cDNA_mRNA',
        #'WT_cDNA_sample',
    ]

    for sample in samples:
        print sample
        file_names = make_file_names(sample)

        start_time = time.time()
        
        trim.trim_adapters(file_names['R1_reads'],
                           file_names['trimmed_reads'],
                           file_names['trimmed_lengths'],
                          )

        ribosomes.filter_contaminant('yeast_rRNA',
                                     file_names['trimmed_reads'],
                                     file_names['no_rRNA_reads'],
                                     file_names['rRNA_sam'],
                                     file_names['rRNA_mapping_log'],
                                     file_names['no_rRNA_lengths'],
                                    )


        ribosomes.map_tophat(file_names)

        ribosomes.further_filtering(file_names['accepted_hits'],
                                    file_names['gtf'],
                                    file_names['filtered_bam'],
                                    file_names['tRNA_lengths'],
                                    file_names['more_rRNA_lengths'],
                                   )

        plot_lengths(file_names)

        from_starts, from_ends = ribosomes.get_aggregate_positions(file_names['gtf'],
                                                                   file_names['filtered_bam'],
                                                                  )
        np.savetxt(file_names['from_starts'], from_starts, fmt='%d')
        np.savetxt(file_names['from_ends'], from_ends, fmt='%d')

        end_time = time.time()
        print 'Processed in {0:0.2f} s'.format(end_time - start_time)
