from __future__ import division
cimport cython
import numpy as np
from collections import defaultdict
import codons

@cython.boundscheck(False)
def fast_stratified_mean_enrichments(codon_counts, gene_names, long num_before, long num_after):
    cdef int position, codon_offset, nucleotide_offset, length
    cdef unsigned int absolute_index, last_absolute_index, absolute_position, codon_index, last_codon_index, nuc_index, i, j
    cdef double ratio, numerator, denominator
    
    cdef long [:, :] occurences = np.zeros((num_before + num_after + 1, 64), int)   
    cdef double [:, :] total_enrichment = np.zeros((num_before + num_after + 1, 64), float)
    
    cdef long [:, :, :] dicodon_occurences = np.zeros((num_before + num_after + 1, 64, 64), int)
    cdef double [:, :, :] dicodon_total_enrichment = np.zeros((num_before + num_after + 1, 64, 64), float)
    
    cdef long [:, :] nuc_occurences = np.zeros((3 * (num_before + num_after + 1), 4), int)   
    cdef double [:, :] nuc_total_enrichment = np.zeros((3 * (num_before + num_after + 1), 4), float)
    
    cdef double [:] ratios
    cdef long [:] codon_indices
    cdef long [:] nucleotide_indices
    
    codon_to_index = {codon: i for i, codon in enumerate(codons.all_codons)}
    index_to_codon = codons.all_codons
    nucleotide_to_index = {n: i for i, n in enumerate('TCAG')}
    index_to_nucleotide = 'TCAG'
    cds_slice = slice(('start_codon', 2), 'stop_codon')

    for gene_name in gene_names:
        counts = codon_counts[gene_name]['relaxed'][cds_slice]
        codon_ids = codon_counts[gene_name]['identities'][cds_slice]
        codon_indices = np.array([codon_to_index[codon] for codon in codon_ids])
        nucleotides = np.array(''.join(codon_ids), dtype='c')
        nucleotide_indices = np.array([nucleotide_to_index[n] for n in nucleotides])
        
        if len(counts) < num_before + num_after + 1:
            raise ValueError(gene_name)

        mean = np.mean(counts[num_before:-num_after])
        if mean == 0:
            raise ValueError(gene_name)
        
        ratios = counts / mean

        length = len(counts)
        
        for position in range(num_before, length - num_after):
            ratio = ratios[position]
            for nucleotide_offset in range(-num_before * 3, num_after * 3):
                absolute_index = position * 3 + nucleotide_offset
                absolute_position = num_before * 3 + nucleotide_offset
                nuc_index = nucleotide_indices[absolute_index]
                nuc_occurences[absolute_position, nuc_index] += 1
                nuc_total_enrichment[absolute_position, nuc_index] += ratio
            
            for codon_offset in range(-num_before, num_after):
                absolute_index = position + codon_offset
                absolute_position = num_before + codon_offset
                codon_index = codon_indices[absolute_index]
                occurences[absolute_position, codon_index] += 1
                total_enrichment[absolute_position, codon_index] += ratio
                
                if codon_offset > -num_before:
                    last_absolute_index = absolute_index - 1
                    last_codon_index = codon_indices[last_absolute_index]
                    dicodon_occurences[absolute_position, codon_index, last_codon_index] += 1
                    dicodon_total_enrichment[absolute_position, codon_index, last_codon_index] += ratio
                
    stratified_mean_enrichments = {}
    for position in range(-num_before, num_after + 1):
        codon_positions = (3 * position, 3 * position + 1, 3 * position + 2)
        dicodon_positions = ((3 * position, 3 * position + 1, 3 * position + 2),
                             (3 * position + 3, 3 * position + 4, 3 * position + 5),
                            )
        stratified_mean_enrichments[codon_positions] = {}
        stratified_mean_enrichments[dicodon_positions] = {}
        absolute_position = num_before + position
        for i in range(64):
            codon_id = index_to_codon[i]
            numerator = total_enrichment[absolute_position, i]
            denominator = max(occurences[absolute_position, i], 1)
            stratified_mean_enrichments[codon_positions][codon_id] = numerator / denominator
            if position < num_after:
                for j in range(64):
                    next_codon_id = index_to_codon[j]
                    numerator = dicodon_total_enrichment[absolute_position, i, j]
                    denominator = max(dicodon_occurences[absolute_position, i, j], 1)
                    stratified_mean_enrichments[dicodon_positions][codon_id, next_codon_id] = numerator / denominator
                    
    for position in range(-num_before * 3, num_after * 3):
        stratified_mean_enrichments[position] = {}
        absolute_position = num_before * 3 + position
        for i in range(4):
            nuc_id = index_to_nucleotide[i]
            numerator = nuc_total_enrichment[absolute_position, i]
            denominator = max(nuc_occurences[absolute_position, i], 1)
            stratified_mean_enrichments[position][nuc_id] = numerator / denominator
                
    return stratified_mean_enrichments
