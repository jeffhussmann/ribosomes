from __future__ import division
cimport cython
import numpy as np
import codons

def make_arrays(num_before, num_after, dtype=float):
    arrays = {'nucleotide': np.zeros((3 * (num_before + num_after + 1), 4), dtype),
              'codon': np.zeros((num_before + num_after + 1, 64), dtype),
              'dicodon': np.zeros((num_before + num_after + 1, 64, 64), dtype),
             }
    return arrays

@cython.boundscheck(False)
@cython.wraparound(False)
def fast_stratified_mean_enrichments(codon_counts,
                                     sorted_gene_names,
                                     breakpoints,
                                     long num_before,
                                     long num_after,
                                     count_type='relaxed',
                                    ):
    cdef int position, codon_offset, nucleotide_offset, length
    cdef int absolute_index, last_absolute_index, absolute_position, codon_index, last_codon_index, nuc_index, i, j
    cdef double ratio, numerator, denominator, mean

    occurence_arrays = make_arrays(num_before, num_after, int)
    total_enrichment_arrays = make_arrays(num_before, num_after, float)

    cdef long [:, ::1] occurences = occurence_arrays['codon']
    cdef double [:, ::1] total_enrichment = total_enrichment_arrays['codon']
    
    cdef long [:, :, ::1] dicodon_occurences = occurence_arrays['dicodon']
    cdef double [:, :, ::1] dicodon_total_enrichment = total_enrichment_arrays['dicodon']
    
    cdef long [:, ::1] nuc_occurences = occurence_arrays['nucleotide']
    cdef double [:, ::1] nuc_total_enrichment = total_enrichment_arrays['nucleotide']
    
    cdef double [::1] ratios
    cdef long [::1] codon_indices
    cdef long [::1] nucleotide_indices
    
    cds_slice = slice(('start_codon', 2), 'stop_codon')

    enrichment_arrays = {}

    total_relevant_counts = 0

    for gene_name in sorted_gene_names:
        counts = codon_counts[gene_name][count_type][cds_slice]
        codon_ids = codon_counts[gene_name]['identities'][cds_slice]
        codon_indices = np.array([codons.codon_to_index[codon] for codon in codon_ids])
        nucleotides = ''.join(codon_ids)
        nucleotide_indices = np.array([codons.nucleotide_to_index[n] for n in nucleotides])
        
        length = len(counts)
        
        if length < num_before + num_after + 1:
            mean = 0.
        else:
            total_relevant_counts += counts[num_before:length - num_after].sum()
            mean = np.mean(counts[num_before:length - num_after])

        if mean != 0.:
            ratios = counts / mean

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
                        dicodon_occurences[absolute_position, last_codon_index, codon_index] += 1
                        dicodon_total_enrichment[absolute_position, last_codon_index, codon_index] += ratio

        if gene_name in breakpoints:
            label = breakpoints[gene_name]
            enrichment_arrays[label] = {}
            for key in occurence_arrays:
                enrichment_arrays[label][key] = total_enrichment_arrays[key] / np.maximum(1, occurence_arrays[key])
                enrichment_arrays[label][key + '_occurences'] = np.copy(occurence_arrays[key])
                enrichment_arrays[label]['total_relevant_counts'] = total_relevant_counts

    stratified_mean_enrichments = StratifiedMeanEnrichments(num_before, num_after, enrichment_arrays)

    return stratified_mean_enrichments

class StratifiedMeanEnrichments(object):
    def __init__(self, num_before, num_after, arrays):
        self.num_before = num_before
        self.num_after = num_after
        self.arrays = arrays
        
    def __getitem__(self, slice_):
        if len(slice_) == 3:
            # Backwards compatibility with calls before cutoff was introduced
            feature, position_slice, label = slice_
            cutoff = '0.10'
        else:
            feature, cutoff, position_slice, label = slice_
        
        if feature.startswith('nucleotide'):
            multiple = 3
        else:
            multiple = 1

        if isinstance(position_slice, slice):
            absolute_start = position_slice.start + self.num_before * multiple
            absolute_stop = position_slice.stop + self.num_before * multiple 
            absolute_slice = slice(absolute_start, absolute_stop)
        elif isinstance(position_slice, (int, long)):
            absolute_slice = position_slice + self.num_before * multiple 
        
        if feature.startswith('nucleotide'):
            if len(label) > 1:
                index = [codons.nucleotide_to_index[l] for l in label]
            else:
                index = codons.nucleotide_to_index[label]
            full_slice = (absolute_slice, index)
        elif feature.startswith('codon'):
            if isinstance(label, list):
                index = [codons.codon_to_index[l] for l in label]
            else:
                index = codons.codon_to_index[label]
            full_slice = (absolute_slice, index)
        elif feature.startswith('dicodon'):
            first_codon, second_codon = label
            first_index = codons.codon_to_index[first_codon]
            second_index = codons.codon_to_index[second_codon]
            full_slice = (absolute_slice, first_index, second_index)
        
        return self.arrays[cutoff][feature][full_slice]
