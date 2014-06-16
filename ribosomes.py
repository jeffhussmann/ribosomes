import numpy as np
import matplotlib.pyplot as plt
import subprocess
import pysam
from itertools import izip
from collections import defaultdict
import gtf
import recycling
import trim
from Sequencing import sam, fastq, mapping_tools, Serialize
from Sequencing.utilities import base_order, base_to_index, base_to_complement_index

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'BlueViolet', 'Gold']
    
# To some extent, this is linked to the definition of simple_CDS, which
# excludes genes that have another gene within 50 bp.
edge_buffer = 50

def map_tophat(reads_file_names,
               bowtie2_index,
               gtf_file_name,
               transcriptome_index,
               tophat_dir,
               num_threads=1,
              ):
    tophat_command = ['tophat2',
                      '--GTF', gtf_file_name,
                      '--no-novel-juncs',
                      '--num-threads', str(num_threads),
                      #'--bt2-mm',
                      '--output-dir', tophat_dir,
                      '--transcriptome-index', transcriptome_index,
                      bowtie2_index,
                      ','.join(reads_file_names),
                     ]
    # tophat maintains its own logs of everything that is written to the
    # console, so discard output.
    subprocess.check_output(tophat_command, stderr=subprocess.STDOUT)

def determine_ambiguity_of_positions(extent, genome_index):
    edge_overlap = 50

    seqname, gene_strand, start, end = extent
    gene_length = abs(end - start) + 1  
    
    ref_seq = str(genome_index[seqname].seq)[start - 2 * edge_overlap:end + edge_overlap + 1]
    
    shape = (3, 3 * edge_overlap + gene_length)
    position_ambiguity = np.empty(shape, int)
    for length in [28, 29, 30]:
        for start in range(3 * edge_overlap + gene_length - length):
            # Explanation of encoding:
            # 0 means a fragment of this length starting at this position would
            # end in an A and therefore can't be the result of poly-A trimming.
            # 1 means a fragment of this length starting at this position would
            # end in a non-A followed by an A, so the read produced can't be
            # distinguished from a read of greater length starting at the same
            # position.
            # 2 means a fragment of this length starting at this position would
            # end in a non-A followed by another non-A and is therefore
            # unambiguous.
            if ref_seq[start + length - 1].upper() == 'A':
                position_ambiguity[length - 28, start] = 0
            elif ref_seq[start + length].upper() == 'A':
                position_ambiguity[length - 28, start] = 1
            else:
                position_ambiguity[length - 28, start] = 2
                
    return position_ambiguity

def get_ratios(first, second):
    assert set(first) == set(second)
    ratios = {key: np.divide(float(first[key]), second[key]) for key in first}
    return ratios

def error_profile(bam_file_name, simple_CDSs, relevant_lengths, max_read_length):
    type_shape = (len(relevant_lengths),
                  max_read_length,
                  fastq.MAX_EXPECTED_QUAL + 1,
                  6,
                  6,
                 )
    type_counts = np.zeros(type_shape, int)
    length_to_index = {length: i for i, length in enumerate(relevant_lengths)}

    bamfile = pysam.Samfile(bam_file_name, 'rb')
    for CDS in simple_CDSs:
        reads = bamfile.fetch(CDS.seqname,
                              CDS.start - edge_buffer,
                              CDS.end + edge_buffer,
                             )
        for read in reads:
            if read.mapq != 50:
                # Non-unique mapping
                continue
            elif read.qlen not in relevant_lengths:
                continue
            elif sam.contains_indel_pysam(read):
                continue
            else:
                strand = '-' if read.is_reverse else '+'
                
                if strand != CDS.strand:
                    continue
                else:
                    alignment = sam.produce_alignment(read, from_pysam=True)

                    if strand == '+':
                        index_lookup = base_to_index
                    else:
                        index_lookup = base_to_complement_index

                    for ref_char, read_char, qual, ref_pos, read_pos in alignment:
                        ref_index = index_lookup[ref_char]
                        read_index = index_lookup[read_char]
                        coords = (length_to_index[read.qlen],
                                  read_pos,
                                  qual,
                                  ref_index,
                                  read_index,
                                 )
                        type_counts[coords] += 1

    return type_counts

def recycling_ratios(rpf_positions_dict, simple_CDSs, genome):
    edge_overlap = 50
        
    ratio_lists = {amino_acid: defaultdict(list) for amino_acid in recycling.back_table}

    for gene in simple_CDSs:
        gene_name = gtf.parse_attribute(gene.attribute)['protein_id']
        codon_counts = get_codon_counts(rpf_positions_dict[gene_name], stringent=False)
        aa_locations = recycling.get_amino_acid_locations(gene, genome)
        if aa_locations == None:
            # MT or internal stop codon
            continue
        for amino_acid in aa_locations:
            location_pairs = izip(aa_locations[amino_acid], aa_locations[amino_acid][1:])
            for (first_position, first_codon), (second_position, second_codon) in location_pairs:
                first_count = codon_counts[first_position]
                second_count = codon_counts[second_position]
                ratio = (first_count, second_count)
                ratio_lists[amino_acid][first_codon, second_codon].append(ratio)
        
    return ratio_lists

def trim_polyA_from_unmapped(unmapped_bam_file_name,
                             trimmed_fastq_file_name,
                             min_length,
                             max_read_length,
                             second_time=False,
                            ):
    reads = (fastq.Read(read.qname, read.seq, read.qual) for read in pysam.Samfile(unmapped_bam_file_name))
    trim.trim(reads,
              trimmed_fastq_file_name,
              min_length,
              max_read_length,
              lambda seq: 0,
              trim.find_poly_A,
              second_time=second_time,
             )
