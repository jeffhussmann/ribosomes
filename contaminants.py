import matplotlib
matplotlib.use('Agg', warn=False)
import random
import os
import numpy as np
import matplotlib.pyplot as plt
import pysam
from collections import defaultdict
from itertools import chain, izip, cycle
from Circles import mapping_tools
from Circles import fasta
from Circles import sam
import gtf

colors = cycle(['b', 'g', 'r', 'c', 'm', 'y', 'BlueViolet', 'Gold'])

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
                              suppress_unaligned_SAM=True,
                              seed_mismatches=1,
                              seed_interval_function='C,1,0',
                              error_file_name=error_fn,
                             )
    sam.make_sorted_indexed_bam(sam_fn, bam_fn)
    os.remove(sam_fn)

def post_filter(input_bam_fn,
                gtf_fn,
                clean_bam_fn,
                more_rRNA_bam_fn,
                tRNA_bam_fn,
                tRNA_bam_sorted_fn,
                other_ncRNA_bam_fn,
                other_ncRNA_bam_sorted_fn,
               ):
    ''' Removes any remaining mappings to rRNA transcripts and any mappings
        to tRNA or other noncoding RNA transcripts.
        If a read has any mapping to an rRNA transcript, write all such mappings
        to more_rRNA_bam_fn with exactly one flagged primary.
        If a read has any mapping to a tRNA transcript, write all such mappings
        to tRNA_bam_fn, with exactly one flagged primary only if there were no
        rRNA mappings.
        If a read has any mapping to any other noncoding RNA transcript, write
        all such mappings to other_ncRNA_bam_fn, with exactly one flagged
        only if there were no rRNA or tRNA mappings.
        Write all reads with no mappings to any noncoding RNA to clean_bam_fn.
    '''
    contaminant_qnames = set()

    rRNA_transcripts, tRNA_transcripts, other_ncRNA_transcripts = gtf.get_noncoding_RNA_transcripts(gtf_fn)

    input_bam_file = pysam.Samfile(input_bam_fn)
   
    # Find reads with any mappings that overlap rRNA or tRNA transcripts and write any
    # such mappings to a contaminant bam file.
    for transcripts, bam_fn in [(rRNA_transcripts, more_rRNA_bam_fn),
                                (tRNA_transcripts, tRNA_bam_fn),
                                (other_ncRNA_transcripts, other_ncRNA_bam_fn),
                               ]:
        with pysam.Samfile(bam_fn, 'wb', template=input_bam_file) as bam_file:
            for transcript in transcripts:
                transcript.build_coordinate_maps()
                overlapping_mappings = input_bam_file.fetch(transcript.seqname,
                                                            transcript.start,
                                                            transcript.end,
                                                           )
                for mapping in overlapping_mappings:
                    # Confirm that there is at least one base from the read
                    # mapped to a position in the transcript (i.e. it isn't just
                    # a spliced read whose junction contains the transcript).
                    if any(p in transcript.genomic_to_transcript and
                           0 <= transcript.genomic_to_transcript[p] < transcript.transcript_length
                           for p in mapping.positions):
                        if mapping.qname not in contaminant_qnames:
                            # This is the first time seeing this qname, so flag
                            # it as primary.
                            mapping.is_secondary = False
                            contaminant_qnames.add(mapping.qname)
                        else:
                            # This qname has already been seen, so flag it as
                            # secondary.
                            mapping.is_secondary = True
                        bam_file.write(mapping)

    input_bam_file.close()

    sam.sort_bam(tRNA_bam_fn, tRNA_bam_sorted_fn)
    sam.sort_bam(other_ncRNA_bam_fn, other_ncRNA_bam_sorted_fn)
         
    # Create a new clean bam file consisting of all mappings of each
    # read that wasn't flagged as a contaminant.
    input_bam_file = pysam.Samfile(input_bam_fn, 'rb')
    with pysam.Samfile(clean_bam_fn, 'wb', template=input_bam_file) as clean_bam_file:
        for mapping in input_bam_file:
            if mapping.qname not in contaminant_qnames:
                clean_bam_file.write(mapping)

def produce_rRNA_coverage(bam_file_name, specific_length=None):
    ''' Counts the number of mappings that overlap each position in the
        reference sequences that bam_file_names were mapped to.
    
        counts: dict (keyed by RNAME) of arrays representing counts for each
                position in RNAME
    '''
    bam_file = pysam.Samfile(bam_file_name)

    if specific_length:
        eligible_reads = (read for read in bam_file if read.qlen == specific_length)
    else:
        eligible_reads = bam_file

    rnames = bam_file.references
    lengths = bam_file.lengths
    counts = {name: np.zeros(length, int) for name, length in zip(rnames, lengths)}
    
    read_names = set()

    for mapping in eligible_reads:
        read_names.add(mapping.qname)
        array = counts[rnames[mapping.tid]]
        for position in mapping.positions:
            array[position] += 1

    return counts

threshold = 0.02
def identify_dominant_stretches(counts, total_reads, bam_fn):
    ''' Identify connected stretches of positions where the fraction of total
        reads mapping is greater than a threshold.
    '''
    boundaries = {}

    for rname in counts:
        # Zero added to the beginning and end so that if a dominant stretch
        # starts at beginning or end, there will be a transition for np.diff
        # to find.
        augmented_counts = np.concatenate(([0], counts[rname], [0]))
        #augemented_counts = np.concatenate(([counts[rname]]))
        normalized_counts = np.true_divide(augmented_counts, total_reads)
        above_threshold = normalized_counts >= threshold
        if not np.any(above_threshold):
            continue
        # np.diff(above_threshold) is True wherever above_threshold changes
        # between True and False. Adapted from a stackoverflow answer.
        # + 1 is so this will be the first thing over the threshold or first
        # thing under it, - 1 is to under the shift of adding a zero at the
        # front.
        above_threshold_boundaries = np.where(np.diff(above_threshold))[0] + 1 - 1

        first_start = np.min(np.where(above_threshold))

        iter_boundaries = iter(above_threshold_boundaries)
        pairs = list(izip(iter_boundaries, iter_boundaries))
        boundaries[rname] = pairs

    # Count the number of reads that overlap any of the dominant stretches.
    overlapping_qnames = set()
    all_qnames = set()
    bam_file = pysam.Samfile(bam_fn)
    for rname in boundaries:
        for start, end in boundaries[rname]:
            reads = bam_file.fetch(rname, start, end)
            for aligned_read in reads:
                overlapping_qnames.add(aligned_read.qname)

    dominant_reads = len(overlapping_qnames)

    return dominant_reads, boundaries

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
    legends = {}
    for i, (rname, length) in enumerate(zip(rnames, lengths)):
        figs[rname], axs[rname] = plt.subplots(figsize=(0.003 * length, 12))
        axs[rname].set_title('rRNA identity: {0}'.format(rname))
        axs[rname].set_xlim(0, length)

    for experiment_name in coverage_data:
        total_reads, counts, color = coverage_data[experiment_name]
        for rname in counts:
            normalized_counts = np.true_divide(counts[rname], total_reads)
            axs[rname].plot(normalized_counts, color=color, label=experiment_name)
            axs[rname].axhline(threshold, linestyle='--', color='black', alpha=0.5)
            legends[rname] = axs[rname].legend(loc='upper right', framealpha=0.5)
            axs[rname].figure.canvas.draw()

    bboxes = {rname: [legends[rname].get_window_extent()] for rname in rnames}
    
    for oligo_name, color in izip(sorted(oligo_mappings), colors):
        for rname, start, end in oligo_mappings[oligo_name]:
            axs[rname].axvspan(start, end, color=color, alpha=0.12, linewidth=0)

            # Annotate the coloring of oligos with their names, avoiding
            # overlapping any other annotations or the legend.
            def attempt_text(y):
                text = axs[rname].annotate(oligo_name,
                                           xy=(float(start + end) / 2, y),
                                           xycoords=('data', 'axes fraction'),
                                           ha='center',
                                           va='top',
                                          )
                axs[rname].figure.canvas.draw()
                return text, text.get_window_extent()
            
            y = 0.995
            text, this_bbox = attempt_text(y)
            while any(this_bbox.fully_overlaps(other_bbox) for other_bbox in bboxes[rname]):
                text.remove()
                y -= 0.01
                text, this_bbox = attempt_text(y)
            bboxes[rname].append(this_bbox)

    for rname in rnames:
        axs[rname].set_xlabel('Position in rRNA')
        axs[rname].set_ylabel('Fraction of all reads mapping to position')
        figs[rname].savefig(fig_fn_template.format(rname), bbox_inches='tight')
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
    for oligo_name, oligo_lengths, color in zip(oligo_names, lengths, colors):
        denominator = np.maximum(oligo_lengths.sum(), 1)
        normalized_lengths = np.true_divide(oligo_lengths, denominator)
        ax.plot(normalized_lengths, 'o-', color=color, label=oligo_name)
    
    ax.legend(loc='upper right', framealpha=0.5)
    
    ax.set_xlim(0, lengths.shape[1] - 1)

    ax.set_xlabel('Length of original RNA fragment')
    ax.set_ylabel('Number of fragments')
    ax.set_title('Distribution of fragment lengths overlapping each oligo')
    
    fig.savefig(fig_fn)
    plt.close(fig)

def extract_rRNA_sequences(genome, rRNA_genes, rRNA_sequences_fn):
    with open(rRNA_sequences_fn, 'w') as rRNA_sequences_fh:
        for gene in rRNA_genes:
            name = gene.attribute['gene_name']
            seq = genome[gene.seqname][gene.start:gene.end + 1]
            record = fasta.make_record(name, seq)
            rRNA_sequences_fh.write(record)
