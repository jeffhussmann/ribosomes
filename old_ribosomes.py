import numpy as np
import matplotlib.pyplot as plt
import subprocess
import pysam
from itertools import izip
from collections import defaultdict
import gtf
import recycling
import trim
from Sequencing import sam, fastq, mapping_tools, Serialize, utilities

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'BlueViolet', 'Gold']
    
# To some extent, this is linked to the definition of simple_CDS, which
# excludes genes that have another gene within 50 bp.
edge_buffer = 50

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
