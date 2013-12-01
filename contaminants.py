import matplotlib
matplotlib.use('Agg', warn=False)
import mapping as mapping_tools
import fasta
import sam
import gtf
import pysam
import ribosomes
import random
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from itertools import chain, izip

def pre_filter(contaminant_index,
               trimmed_reads_fn,
               filtered_reads_fn,
               sam_fn,
               bam_fn,
               error_fn,
              ):
    ''' Maps reads in trimmed_reads_fn to contaminant_index. Records reads that
        don't map in filtered_reads_fn. 
    '''
    mapping_tools.map_bowtie2(trimmed_reads_fn,
                              contaminant_index,
                              sam_fn,
                              unaligned_reads_file_name=filtered_reads_fn,
                              threads=1,
                              report_all=True,
                              explicit_path=True,
                              error_file_name=error_fn,
                             )
    sam.make_sorted_indexed_bam(sam_fn, bam_fn)

def post_filter(input_bam_fn,
                gtf_fn,
                clean_bam_fn,
                more_rRNA_bam_fn,
                tRNA_bam_fn,
               ):
    ''' Removes any remaining mappings to tRNA or rRNA genes.
        Removes secondary mappings so that each read is represented at most
        once, but this DOES NOT mean that only unique mappings are retained.
    '''
    qnames_already_used = set()

    rRNA_genes = gtf.get_rRNA_genes(gtf_fn)
    tRNA_genes = gtf.get_tRNA_genes(gtf_fn)

    input_bam_file = pysam.Samfile(input_bam_fn, 'rb')
   
    # Find reads with any mappings that overlap rRNA or tRNA genes and write one
    # mapping for each such read to a contaminant bam file.
    for genes, bam_fn in [(rRNA_genes, more_rRNA_bam_fn),
                          (tRNA_genes, tRNA_bam_fn)]:
        contaminant_reads = defaultdict(list)
        for gene in genes:
            overlapping_mappings = input_bam_file.fetch(gene.seqname,
                                                        gene.start,
                                                        gene.end,
                                                       )
            for mapping in overlapping_mappings:
                contaminant_reads[mapping.qname].append(mapping)

        with pysam.Samfile(bam_fn, 'wb', header=input_bam_file.header) as bam_file:
            for qname in contaminant_reads:
                if qname not in qnames_already_used:
                    mapping = random.choice(contaminant_reads[qname])
                    bam_file.write(mapping)
                    qnames_already_used.add(qname)

    input_bam_file.close()
         
    # Create a new clean bam file consisting of the primary mapping of each
    # read that wasn't flagged as a contaminant.
    input_bam_file = pysam.Samfile(input_bam_fn, 'rb')
    with pysam.Samfile(clean_bam_fn, 'wb', header=input_bam_file.header) as clean_bam_file:
        for mapping in input_bam_file:
            if mapping.qname not in qnames_already_used and not mapping.is_secondary:
                clean_bam_file.write(mapping)

    pysam.index(clean_bam_fn)

def produce_rRNA_coverage(bam_file_names):
    ''' Counts the number of mappings that overlap each position in the
        reference sequences that bam_file_names were mapped to.
    
        total_reads: number of distinct QNAME's present in bam_file_names
        counts: dict (keyed by RNAME) of arrays representing counts for each
                position in RNAME
    '''
    bam_files = [pysam.Samfile(fn, 'rb') for fn in bam_file_names]
    all_reads = chain.from_iterable(bam_files)

    rnames = bam_files[0].references
    lengths = bam_files[0].lengths
    counts = {name: np.zeros(length, int) for name, length in zip(rnames, lengths)}
    
    read_names = set()

    for mapping in all_reads:
        read_names.add(mapping.qname)
        array = counts[rnames[mapping.tid]]
        for position in mapping.positions:
            array[position] += 1

    total_reads = len(read_names)

    return total_reads, counts

def plot_rRNA_coverage(coverage_data, oligos_sam_fn, fig_fn_template):
    ''' Plots the number of mappings that overlap each position in the reference
        sequences mapped to. Highlights the regions targeted by oligos.
    '''
    oligos_sam_file = pysam.Samfile(oligos_sam_fn, 'r')
    rnames = oligos_sam_file.references
    lengths = oligos_sam_file.lengths
    oligo_mappings = load_oligo_mappings(oligos_sam_fn)
    
    figs = {}
    axs = {}
    for i, (rname, length) in enumerate(zip(rnames, lengths)):
        figs[rname], axs[rname] = plt.subplots(figsize=(18, 12))
        axs[rname].set_title('rRNA identity - ' + rname)
        axs[rname].set_xlim(0, length)

    for experiment_name in coverage_data:
        total_reads, counts = coverage_data[experiment_name]
        for rname in counts:
            normalized_counts = np.true_divide(counts[rname], total_reads)
            axs[rname].plot(normalized_counts, label=experiment_name)

    for oligo, color in izip(oligo_mappings, ribosomes.colors):
        for rname, start, end in oligo_mappings[oligo]:
            axs[rname].axvspan(start, end, color=color, alpha=0.12, linewidth=0)
            _, y_max = axs[rname].get_ylim()
            axs[rname].text(float(start + end) / 2,
                            y_max,
                            oligo, 
                            horizontalalignment='center',
                            verticalalignment='top',
                           )
    for rname in rnames:
        leg = axs[rname].legend(loc='upper right', fancybox=True)
        leg.get_frame().set_alpha(0.5)
        axs[rname].set_xlabel('Position in rRNA')
        axs[rname].set_ylabel('Fraction of all reads mapping to position')
        figs[rname].savefig(fig_fn_template.format(rname))
        plt.close(figs[rname])

def load_oligo_mappings(oligos_sam_fn):
    oligos_sam_file = pysam.Samfile(oligos_sam_fn, 'r')
    oligo_mappings = defaultdict(list)
    for aligned_read in oligos_sam_file:
        positions = aligned_read.positions
        rname = oligos_sam_file.getrname(aligned_read.tid)
        extent = (rname, min(positions), max(positions))
        oligo_mappings[aligned_read.qname].append(extent)
    return oligo_mappings

def get_oligo_hit_lengths(bam_fn,
                          oligos_fasta_fn,
                          oligos_sam_fn,
                          max_read_length):
    oligo_mappings = load_oligo_mappings(oligos_sam_fn)
    bam_file = pysam.Samfile(bam_fn, 'rb')

    oligo_names = [read.name for read in fasta.reads(oligos_fasta_fn)]
    lengths = np.zeros((len(oligo_names), max_read_length + 1), int)

    for oligo_number, oligo_name in enumerate(oligo_names):
        for rname, start, end in oligo_mappings[oligo_name]:
            reads = bam_file.fetch(rname, start, end)
            for aligned_read in reads:
                if not aligned_read.is_secondary:
                    lengths[oligo_number][aligned_read.qlen] += 1
    
    return lengths

def plot_oligo_hit_lengths(oligos_fasta_fn, lengths, fig_fn):
    oligo_names = [read.name for read in fasta.reads(oligos_fasta_fn)]
    if len(oligo_names) == 0:
        # If no oligos have been defined, there is no picture to make.
        return None
    
    fig, ax = plt.subplots(figsize=(18, 12))
    for oligo_name, oligo_lengths, color in zip(oligo_names, lengths, ribosomes.colors):
        normalized_lengths = np.true_divide(oligo_lengths, oligo_lengths.sum())
        ax.plot(normalized_lengths, 'o-', color=color, label=oligo_name)
    
    leg = ax.legend(loc='upper right', fancybox=True)
    leg.get_frame().set_alpha(0.5)
    
    ax.set_xlim(0, lengths.shape[1] - 1)

    ax.set_xlabel('Length of original RNA fragment')
    ax.set_ylabel('Number of fragments')
    ax.set_title('Distribution of fragment lengths overlapping each oligo')
    
    fig.savefig(fig_fn)
    plt.close(fig)

def extract_rRNA_sequences(genome, rRNA_genes, rRNA_sequences_fn):
    with open(rRNA_sequences_fn, 'w') as rRNA_sequences_fh:
        for gene in rRNA_genes:
            name = gtf.parse_attribute(gene.attribute)['gene_name']
            seq = genome[gene.seqname][gene.start:gene.end + 1]
            record = fasta.make_record(name, seq)
            rRNA_sequences_fh.write(record)
