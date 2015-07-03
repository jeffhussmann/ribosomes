from __future__ import division
cimport cython
import numpy as np
import codons

def find_minimum_density_breakpoints(codon_counts,
                                     min_means,
                                     exclude_from_start,
                                     exclude_from_end,
                                     count_type,
                                    ):
    ''' min_means: list of minimum mean density breakpoints
        exclude_from_edges: (exclude_from_start, exclude_from_end) tuples
    '''
    cds_slice = slice(('start_codon', 2), 'stop_codon')

    breakpoints = {}
    means = {}
    for gene_name in codon_counts:
        counts = codon_counts[gene_name][count_type][cds_slice]
        
        length = len(counts)
        if length <= exclude_from_start + exclude_from_end:
            mean = -1
        else:
            mean = np.mean(counts[exclude_from_start:length - exclude_from_end])
            
        means[gene_name] = mean

    sorted_names = sorted(codon_counts.keys(),
                          key=lambda n: (means[n], n),
                          reverse=True,
                         )
    for min_mean in min_means:
        name = [name for name in sorted_names if means[name] > min_mean][-1]
        breakpoints[name] = min_mean

    return sorted_names, breakpoints

def make_arrays(num_around, dtype=float):
    arrays = {'nucleotide': np.zeros((3 * (2 * num_around + 1), 4), dtype),
              'codon': np.zeros((2 * num_around + 1, 64), dtype),
              'dicodon': np.zeros((2 * num_around + 1, 64, 64), dtype),
             }
    return arrays

def make_hdf5_key(min_mean, exclude_from_start, exclude_from_end):
    key = '{0:0.2f},{1},{2}'.format(min_mean, exclude_from_start, exclude_from_end)
    return key

@cython.boundscheck(False)
@cython.wraparound(False)
def fast_stratified_mean_enrichments(codon_counts,
                                     exclude_from_edges,
                                     min_means,
                                     long num_around,
                                     count_type='relaxed',
                                     keys=['codon', 'nucleotide'],
                                    ):
    cdef int position, codon_offset, nucleotide_offset, length
    cdef int absolute_index, last_absolute_index, absolute_position, codon_index, last_codon_index, nuc_index, i, j, exclude_from_start, exclude_from_end, offset_start, offset_end
    cdef double ratio, numerator, denominator, mean

    cdef long [:, ::1] occurences
    cdef double [:, ::1] total_enrichment
    
    cdef long [:, :, ::1] dicodon_occurences
    cdef double [:, :, ::1] dicodon_total_enrichment
    
    cdef long [:, ::1] nuc_occurences
    cdef double [:, ::1] nuc_total_enrichment
    
    cdef double [::1] ratios
    cdef long [::1] codon_indices
    cdef long [::1] nucleotide_indices
    
    enrichment_arrays = {}
   
    for exclude_from_start, exclude_from_end in exclude_from_edges: 
        occurence_arrays = make_arrays(num_around, int)
        total_enrichment_arrays = make_arrays(num_around, float)
        
        occurences = occurence_arrays['codon']
        total_enrichment = total_enrichment_arrays['codon']
        
        dicodon_occurences = occurence_arrays['dicodon']
        dicodon_total_enrichment = total_enrichment_arrays['dicodon']
        
        nuc_occurences = occurence_arrays['nucleotide']
        nuc_total_enrichment = total_enrichment_arrays['nucleotide']
        
        cds_slice = slice(('start_codon', 2), 'stop_codon')

        total_relevant_counts = 0

        sorted_gene_names, breakpoints = find_minimum_density_breakpoints(codon_counts,
                                                                          min_means,
                                                                          exclude_from_start,
                                                                          exclude_from_end,
                                                                          count_type,
                                                                         )

        for gene_name in sorted_gene_names:
            counts = codon_counts[gene_name][count_type][cds_slice]
            codon_ids = codon_counts[gene_name]['identities'][cds_slice]
            codon_indices = np.array([codons.codon_to_index[codon] for codon in codon_ids])
            nucleotides = ''.join(codon_ids)
            nucleotide_indices = np.array([codons.nucleotide_to_index[n] for n in nucleotides])
            
            length = len(counts)
            
            if length <= exclude_from_start + exclude_from_end:
                mean = 0.
            else:
                total_relevant_counts += counts[exclude_from_start:length - exclude_from_end].sum()
                mean = np.mean(counts[exclude_from_start:length - exclude_from_end])

            if mean != 0.:
                ratios = counts / mean

                for position in range(exclude_from_start, length - exclude_from_end):
                    ratio = ratios[position]
                    offset_start = max(-position, -num_around)
                    offset_end = min(length - position, num_around + 1)
                    for nucleotide_offset in range(offset_start * 3, offset_end * 3):
                        absolute_index = position * 3 + nucleotide_offset
                        absolute_position = num_around * 3 + nucleotide_offset
                        nuc_index = nucleotide_indices[absolute_index]
                        nuc_occurences[absolute_position, nuc_index] += 1
                        nuc_total_enrichment[absolute_position, nuc_index] += ratio
                    
                    for codon_offset in range(offset_start, offset_end):
                        absolute_index = position + codon_offset
                        absolute_position = num_around + codon_offset
                        codon_index = codon_indices[absolute_index]
                        occurences[absolute_position, codon_index] += 1
                        total_enrichment[absolute_position, codon_index] += ratio
                        
                        if codon_offset > offset_start:
                            last_absolute_index = absolute_index - 1
                            last_codon_index = codon_indices[last_absolute_index]
                            dicodon_occurences[absolute_position, last_codon_index, codon_index] += 1
                            dicodon_total_enrichment[absolute_position, last_codon_index, codon_index] += ratio

            if gene_name in breakpoints:
                label = make_hdf5_key(breakpoints[gene_name], exclude_from_start, exclude_from_end)
                enrichment_arrays[label] = {}
                for key in keys:
                    enrichment_arrays[label][key] = total_enrichment_arrays[key] / np.maximum(1, occurence_arrays[key])
                    enrichment_arrays[label][key + '_occurences'] = np.copy(occurence_arrays[key])
                    enrichment_arrays[label]['total_relevant_counts'] = total_relevant_counts

    stratified_mean_enrichments = StratifiedMeanEnrichments(num_around, enrichment_arrays)

    return stratified_mean_enrichments

class StratifiedMeanEnrichments(object):
    def __init__(self, num_around, arrays):
        self.num_around = num_around
        self.arrays = arrays
        
    def __getitem__(self, slice_):
        if len(slice_) == 3:
            # Backwards compatibility with calls before condition was introduced
            feature, position_slice, label = slice_
            condition = make_hdf5_key(0.10, 90, 90)
        else:
            condition, feature, position_slice, label = slice_
            condition = make_hdf5_key(*condition)
        
        if feature.startswith('nucleotide'):
            multiple = 3
        else:
            multiple = 1

        if isinstance(position_slice, slice):
            absolute_start = position_slice.start + self.num_around * multiple
            absolute_stop = position_slice.stop + self.num_around * multiple 
            absolute_slice = slice(absolute_start, absolute_stop, position_slice.step)
        elif isinstance(position_slice, (int, long)):
            absolute_slice = position_slice + self.num_around * multiple 
        
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
        
        array = self.arrays[condition][feature][()]
        return array[full_slice]
