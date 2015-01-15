''' Utilities for counting reads mapped to each position in a gene. '''

import numpy as np
import numbers
import pysam
from collections import Counter
from itertools import cycle, product
import codons
import gtf
import Sequencing.Serialize
import trim

class PositionCounts(object):
    ''' Wrapper around an array of counts with annotated landmark positions '''
    def __init__(self, landmarks, left_buffer, right_buffer, data=None, dtype=np.int32):
        self.landmarks = landmarks
        self.left_buffer = left_buffer
        self.right_buffer = right_buffer

        leftmost_name = min(landmarks, key=landmarks.get)
        leftmost_value = landmarks[leftmost_name]
        rightmost_name = max(landmarks, key=landmarks.get)
        rightmost_value = landmarks[rightmost_name]

        # For historical reasons, rightmost_value is assumed to point to one
        # past the 'interesting' region, so that it is the first index of the
        # left buffer. This means that there is NOT a '+ 1' at the end of this
        # length expression as might be expected.
        length = (rightmost_value + right_buffer) - (leftmost_value - left_buffer)
        self.landmark_to_index = {name: left_buffer + (landmarks[name] - leftmost_value)
                                  for name in landmarks}
       
        if data == None:
            self.data = np.zeros(length, dtype)
        else:
            if len(data) != length:
                raise ValueError(len(data), length)

            self.data = data

    def transform_relative_slice(self, relative_slice, end_allowed=False):
        if isinstance(relative_slice, str):
            # relative_slice is just a landmark
            absolute_slice = self.landmark_to_index[relative_slice]

        elif isinstance(relative_slice, tuple):
            # relative_slice is a (landmark, slice) tuple 
            landmark, key = relative_slice

            if key == None:
                absolute_slice = self.transform_relative_slice(landmark)

            elif isinstance(key, (int, long)):
                absolute_slice = self.landmark_to_index[landmark] + key
                if absolute_slice < 0 or absolute_slice > len(self.data):
                    raise IndexError('Length of data - {0}, attempted to access {1} from key of {2}'.format(len(self.data), absolute_slice, key))
                # When it is the stop of a slice, adjusted_key must be allowed to be
                # equal to len(self.data) to allow a slice to contain the end point.
                if absolute_slice == len(self.data) and not end_allowed:
                    raise IndexError('Attempted to access end point')

            elif isinstance(key, slice):
                if key.stop == None:
                    raise ValueError('None not allowed as stop in slice')

                if key.step != None and key.step < 0:
                    raise ValueError('Negative step not allowed')

                start = self.transform_relative_slice((landmark, key.start))
                stop = self.transform_relative_slice((landmark, key.stop), end_allowed=True)
                absolute_slice = slice(start, stop, key.step)

            elif isinstance(key, (list, np.ndarray)):
                absolute_slice = [self.transform_relative_slice((landmark, k)) for k in key]

            else:
                raise TypeError(type(key))

        elif isinstance(relative_slice, slice):
            # relative_slice.start and relative_slice.stop are each one of the
            # two previous categories
            start = self.transform_relative_slice(relative_slice.start)
            stop = self.transform_relative_slice(relative_slice.stop)
            absolute_slice = slice(start, stop, relative_slice.step)

        else:
            raise NotImplementedError(relative_slice)

        return absolute_slice

    def __getitem__(self, relative_slice):
        absolute_slice = self.transform_relative_slice(relative_slice)
        return self.data[absolute_slice]

    def __setitem__(self, relative_slice, value):
        absolute_slice = self.transform_relative_slice(relative_slice)
        self.data[absolute_slice] = value

    def __iadd__(self, other):
        self.check_compatibility(other)
        self.data += other.data

        return self

    def sum(self):
        return self.data.sum()

    def argmax_over_slice(self, landmark, key):
        if isinstance(key, (list, np.ndarray)):
            counts = self[landmark, key]
            return key[counts.argmax()]
        else:
            raise NotImplementedError

    def n_largest_over_slice(self, n, landmark_and_key):
        landmark, key = landmark_and_key
        
        if isinstance(key, (list, np.ndarray)):
            counts = self[landmark, key]
            n_largest = counts.argsort()[:-(n + 1):-1]
            return key[n_largest]
        else:
            raise NotImplementedError

    @property
    def CDS_length(self):
        return self.landmarks['stop_codon'] - self.landmarks['start_codon']

    def check_compatibility(self, other):
        ''' Check if self and other can be used for element-wise operations. '''
        if self.landmarks != other.landmarks:
            raise ValueError('Unequal landmarks:', self.landmarks, other.landmarks)
        if self.left_buffer != other.left_buffer or self.right_buffer != other.right_buffer:
            raise ValueError('Unequal buffer lengths')

    def __div__(self, other):
        if isinstance(other, PositionCounts):
            self.check_compatibility(other)
            return PositionCounts(self.landmarks,
                                  self.left_buffer,
                                  self.right_buffer,
                                  data=np.true_divide(self.data, other.data),
                                 )
        elif isinstance(other, numbers.Number):
            return PositionCounts(self.landmarks,
                                  self.left_buffer,
                                  self.right_buffer,
                                  data=np.true_divide(self.data, other),
                                 )
    
    def __add__(self, other):
        if isinstance(other, PositionCounts):
            self.check_compatibility(other)
            return PositionCounts(self.landmarks,
                                  self.left_buffer,
                                  self.right_buffer,
                                  data=(self.data + other.data),
                                 )
        else:
            raise ValueError('bad types in PositionCounts addition')
    
    def __sub__(self, other):
        if isinstance(other, PositionCounts):
            self.check_compatibility(other)
            return PositionCounts(self.landmarks,
                                  self.left_buffer,
                                  self.right_buffer,
                                  data=(self.data - other.data),
                                 )
        else:
            raise ValueError('bad types in PositionCounts subtraction')

def convert_to_three_prime(position_counts, length):
    ''' Shift position counts that represent 5' edges of fragments of given
        length to represetnt 3' edges.
    '''
    if length >= len(position_counts.data):
        # All counts will be destroyed, so just return zeros.
        three_prime_counts = PositionCounts(position_counts.landmarks,
                                            position_counts.left_buffer,
                                            position_counts.right_buffer,
                                           )
    else:
        converted_data = np.concatenate((np.zeros(length - 1, position_counts.data.dtype),
                                         position_counts.data[:-(length - 1)],
                                        ),
                                       )
        three_prime_counts = PositionCounts(position_counts.landmarks,
                                            position_counts.left_buffer,
                                            position_counts.right_buffer,
                                            data=converted_data,
                                           )
    return three_prime_counts

# Number of nucleotide positions to include on the left and right side of
# counts of read positions.
edge_buffer = 50
left_buffer = 4 * edge_buffer
right_buffer = 4 * edge_buffer

# Number of codons to include on either side of counts of codon positions.
codon_buffer = 10

def get_Transcript_extent_position_counts(transcript,
                                          clean_bam_fn,
                                          relevant_lengths,
                                          left_buffer=0,
                                          right_buffer=0,
                                         ):
    bam_file = pysam.Samfile(clean_bam_fn)
    landmarks = {'start': 0,
                 'stop': transcript.extent_length,
                }
    position_counts = {l: PositionCounts(landmarks, left_buffer, right_buffer)
                       for l in relevant_lengths + ['all']}
    
    # fetch raises a ValueError if given a negative start, but it doesn't 
    # care if the end is valid.
    left_edge = max(0, min(transcript.genomic_to_extent) - left_buffer)
    right_edge = max(transcript.genomic_to_extent) + right_buffer
    overlapping_reads = bam_file.fetch(transcript.seqname, left_edge, right_edge)
    for read in overlapping_reads:
        if read.mapq != 50:
            continue

        read_strand = '-' if read.is_reverse else '+'
        if read_strand != transcript.strand:
            continue

        if read_strand == '+':
            five_prime_position = read.pos
        elif read_strand == '-':
            five_prime_position = read.aend - 1

        if five_prime_position in transcript.genomic_to_extent:
            extent_coord = transcript.genomic_to_extent[five_prime_position]

            position_counts['all']['start', extent_coord] += 1
            if read.qlen in relevant_lengths:
                position_counts[read.qlen]['start', extent_coord] += 1

    return position_counts

def get_Transcript_position_counts(clean_bam_fn,
                                   transcripts,
                                   relevant_lengths,
                                   left_buffer=left_buffer,
                                   right_buffer=right_buffer,
                                  ):
    gene_infos = {}
    bam_file = pysam.Samfile(clean_bam_fn)
    
    max_nongenomic_length = 5

    for transcript in transcripts:
        transcript.build_coordinate_maps(left_buffer, right_buffer)

        nonunique = 0
        alternatively_spliced = 0
        
        landmarks = {'start': 0,
                     'start_codon': transcript.transcript_start_codon,
                     'stop_codon': transcript.transcript_stop_codon,
                     'end': transcript.transcript_length,
                    }
        five_prime_positions = {l: PositionCounts(landmarks, left_buffer, right_buffer)
                                for l in relevant_lengths + ['all']}
        
        three_prime_positions = {l: PositionCounts(landmarks, left_buffer, right_buffer)
                                 for l in range(max_nongenomic_length + 1) + ['all']}

        transcript_sequence = transcript.get_transcript_sequence(left_buffer, right_buffer)
        five_prime_positions['sequence'] = transcript_sequence
        three_prime_positions['sequence'] = transcript_sequence
        
        # fetch raises a ValueError if given a negative start, but it doesn't 
        # care if the end is valid.
        left_edge = max(0, transcript.start - left_buffer)
        right_edge = transcript.end + right_buffer
        overlapping_reads = bam_file.fetch(transcript.seqname, left_edge, right_edge)
        for read in overlapping_reads:
            if any(transcript.is_spliced_out(position) for position in read.positions):
                alternatively_spliced += 1
                continue
            
            if read.mapq != 50:
                nonunique += 1
                continue

            read_strand = '-' if read.is_reverse else '+'
            if read_strand != transcript.strand:
                continue
            
            left_edge = read.pos
            right_edge = read.aend - 1
            
            if read_strand == '+':
                five_prime_position = left_edge
                three_prime_position = right_edge
            elif read_strand == '-':
                five_prime_position = right_edge
                three_prime_position = left_edge

            if five_prime_position in transcript.genomic_to_transcript:
                transcript_coord = transcript.genomic_to_transcript[five_prime_position]
                five_prime_positions['all']['start', transcript_coord] += 1
                
                if read.qlen in relevant_lengths:
                    five_prime_positions[read.qlen]['start', transcript_coord] += 1
            
            if three_prime_position in transcript.genomic_to_transcript:
                transcript_coord = transcript.genomic_to_transcript[three_prime_position]
                three_prime_positions['all']['start', transcript_coord] += 1

                nongenomic_length = trim.get_nongenomic_length(read)
                if nongenomic_length <= max_nongenomic_length:
                    three_prime_positions[nongenomic_length]['start', transcript_coord] += 1

        gene_infos[transcript.name] = {'CDS_length': transcript.CDS_length,
                                       'five_prime_positions': five_prime_positions,
                                       'three_prime_positions': three_prime_positions,
                                       'nonunique': nonunique,
                                       'alternatively_spliced': alternatively_spliced,
                                      }
        transcript.delete_coordinate_maps()

    return gene_infos

def get_joint_position_counts_sparse(extended_bam_fn, transcript, left_buffer=left_buffer, right_buffer=right_buffer):
    bam_file = pysam.Samfile(extended_bam_fn)
    transcript.build_coordinate_maps(left_buffer, right_buffer)

    landmarks = {'start': 0,
                 'start_codon': transcript.transcript_start_codon,
                 'stop_codon': transcript.transcript_stop_codon,
                 'end': transcript.transcript_length,
                }
    joint_position_counts = Counter()
    
    # fetch raises a ValueError if given a negative start, but it doesn't 
    # care if the end is valid.
    left_edge = max(0, transcript.start - left_buffer)
    right_edge = transcript.end + right_buffer
    overlapping_reads = bam_file.fetch(transcript.seqname, left_edge, right_edge)
    for read in overlapping_reads:
        if any(transcript.is_spliced_out(position) for position in read.positions):
            continue
        
        if read.mapq != 50:
            continue

        read_strand = '-' if read.is_reverse else '+'
        if read_strand != transcript.strand:
            continue

        if read_strand == '+':
            five_prime_position = read.pos
            three_prime_position = read.aend - 1
        elif read_strand == '-':
            five_prime_position = read.aend - 1
            three_prime_position = read.pos

        if five_prime_position in transcript.genomic_to_transcript and \
           three_prime_position in transcript.genomic_to_transcript:

            five_prime_transcript_coord = transcript.genomic_to_transcript[five_prime_position]
            three_prime_transcript_coord = transcript.genomic_to_transcript[three_prime_position] - transcript.transcript_stop_codon
            joint_position_counts[five_prime_transcript_coord, three_prime_transcript_coord] += 1

    transcript.delete_coordinate_maps()

    return joint_position_counts

A_site_offsets = {'ingolia_cell': {29: 15,
                                   30: 15,
                                   31: 16,
                                   32: 16,
                                   33: 16,
                                   34: 17,
                                   35: 17,
                                  },
                  'guo_nature':   {27: 15,
                                   28: 15,
                                   29: 15,
                                   30: 15,
                                   31: 15,
                                   32: 15,
                                  },
                  'yeast':        {28: 15,
                                   29: 15,
                                   30: 16,
                                  },
                  'yeast_stringent': {28: 15,
                                     },
                  'yeast_anisomycin': {20: 15,
                                       21: 15,
                                       22: 16,
                                       23: 16,
                                      }
                 }

def compute_codon_counts(position_counts, offset_type):
    CDS_length = position_counts.values()[0].CDS_length

    if CDS_length % 3 != 0:
        raise ValueError('CDS length not divisible by 3')

    # Note: CDS_length is the index of the first nucleotide of the stop codon.
    # Ingolia's original model never has the stop codon in the A site, but
    # subsequent data show an accumulation of (typically length 29 or 30) reads
    # that do advance this far.
    num_codons = CDS_length // 3
    landmarks = {'start_codon': 0,
                 'stop_codon': num_codons,
                }
    codon_counts = PositionCounts(landmarks, codon_buffer, codon_buffer)

    recorded_lengths = set(position_counts.keys())
    known_A_site_lengths = set(A_site_offsets[offset_type].keys())

    for length in recorded_lengths & known_A_site_lengths:
        A_site_offset = A_site_offsets[offset_type][length]
        start_index = -A_site_offset - (codon_buffer * 3)
        end_index = CDS_length - A_site_offset + (codon_buffer * 3)
        in_frame = slice(start_index, end_index, 3)
        one_behind = slice(start_index - 1, end_index - 1, 3)
        one_ahead = slice(start_index + 1, end_index + 1, 3)
        codon_counts.data += position_counts[length]['start_codon', in_frame] + \
                             position_counts[length]['start_codon', one_behind] + \
                             position_counts[length]['start_codon', one_ahead]

    sequence_slice = slice(('start_codon', -codon_buffer * 3), ('stop_codon', codon_buffer * 3))
    sequence = ''.join(position_counts['sequence'][sequence_slice])
    codon_identities_list = list(codons.codons_from_seq(sequence))
    codon_identities = PositionCounts(landmarks,
                                      codon_buffer,
                                      codon_buffer,
                                      data = np.asarray(codon_identities_list),
                                     )
    return codon_counts, codon_identities

def compute_metagene_positions(CDSs, position_counts, max_CDS_length):
    ''' max_CDS_length needs to be passed in because it may reflect the max of
        more than just the CDS being considered here.
    '''
    relevant_lengths, left_buffer, right_buffer = extract_lengths_and_buffers(position_counts)

    landmarks = {'start': 0,
                 'start_codon': 0,
                 'stop_codon': max_CDS_length,
                 'end': max_CDS_length,
                }

    def make_PositionCounts_dictionary(**kwargs):
        d = {length: PositionCounts(landmarks, left_buffer, right_buffer, **kwargs)
             for length in relevant_lengths}
        return d

    metagene_positions = {}
    for landmark in landmarks:
        metagene_positions[landmark] = make_PositionCounts_dictionary()
        uniform_key = '{0}_uniform'.format(landmark)
        metagene_positions[uniform_key] = make_PositionCounts_dictionary(dtype=float)
        for b in 'TCAG':
            key = '{0}_{1}'.format(landmark, b)
            metagene_positions[key] = make_PositionCounts_dictionary()
            
            uniform_key = '{0}_{1}_uniform'.format(landmark, b)
            metagene_positions[uniform_key] = make_PositionCounts_dictionary(dtype=float)

    for CDS in CDSs:
        # Skip CDSs that overlap other qualifying features or that have another
        # too close.
        if CDS.num_overlapping > 0:
            continue

        CDS.build_coordinate_maps(left_buffer, right_buffer)
        downstream = CDS.genomic_to_transcript.get(CDS.downstream_transcript)
        if downstream != None:
            if downstream - CDS.transcript_stop_codon < 100:
                continue

        counts = position_counts[CDS.name]
        CDS_length = CDS.CDS_length
        
        landmark_slices = [('start', slice(-left_buffer, CDS_length)),
                           ('start_codon', slice(-left_buffer, CDS_length)),
                           ('stop_codon', slice(-CDS_length, right_buffer)),
                           ('end', slice(-CDS_length, right_buffer)),
                          ]

        for length in relevant_lengths:
            for landmark_slice in landmark_slices:
                landmark, _ = landmark_slice
                sliced_counts = counts[length][landmark_slice]
                metagene_positions[landmark][length][landmark_slice] += sliced_counts

                density = sliced_counts.sum() / float(len(sliced_counts))
                uniform_counts = np.ones_like(sliced_counts) * density
                    
                uniform_key = '{0}_uniform'.format(landmark)
                metagene_positions[uniform_key][length][landmark_slice] += uniform_counts
                
                for b in 'TCAG':
                    base_mask = counts['sequence'][landmark_slice] == b
                    base_counts = np.multiply(sliced_counts, base_mask)
                    key = '{0}_{1}'.format(landmark, b)
                    metagene_positions[key][length][landmark_slice] += base_counts 

                    # Metagene base composition could be skewed simply because
                    # highly expressed genes happen to have a particular base
                    # at a given offset, rather than because that base is really
                    # occuring more often than its neighbors.
                    # To control for expression-weighted composition, compute
                    # the average read density for each gene and sum up
                    # base-masked arrays of it.
                    uniform_key = '{0}_{1}_uniform'.format(landmark, b)
                    uniform_base_counts = np.multiply(uniform_counts, base_mask)
                    metagene_positions[uniform_key][length][landmark_slice] += uniform_base_counts
                    
        CDS.delete_coordinate_maps()

    return metagene_positions

def extract_lengths_and_buffers(position_counts):
    arbitrary_gene = position_counts.itervalues().next()
    relevant_lengths = [l for l in arbitrary_gene if l != 'sequence']
    arbitrary_length = relevant_lengths[0]
    arbitrary_counts = arbitrary_gene[arbitrary_length]
    left_buffer = arbitrary_counts.left_buffer
    right_buffer = arbitrary_counts.right_buffer
    return relevant_lengths, left_buffer, right_buffer

def compute_averaged_codon_densities(codon_counts, offset_key='relaxed', names_to_skip=set()): 
    # To reduce noise, genes with less than min_counts total counts are ignored.
    min_counts = 64
    try:
        max_length = max(counts[offset_key].CDS_length
                         for counts in codon_counts.itervalues()
                         if counts[offset_key].sum() >= min_counts
                        )
    except ValueError:
        # max() arg is an empty sequence
        max_length = 10000

    max_length = max(max_length, 2000)

    landmarks = {'start_codon': 0,
                 'stop_codon': max_length,
                }
    sum_of_normalized_from_start = PositionCounts(landmarks, codon_buffer, codon_buffer, dtype=float)
    sum_of_normalized_from_end = PositionCounts(landmarks, codon_buffer, codon_buffer, dtype=float)
    long_enough_genes_from_start = PositionCounts(landmarks, codon_buffer, codon_buffer)
    long_enough_genes_from_end = PositionCounts(landmarks, codon_buffer, codon_buffer)

    uniform = PositionCounts(landmarks, codon_buffer, codon_buffer)
    uniform.data[...] = 1

    for name, counts_offset_groups in codon_counts.iteritems():
        if name in names_to_skip:
            print 'skipping', name
            continue

        counts = counts_offset_groups[offset_key]
        if counts.sum() < min_counts:
            continue

        num_codons = counts.CDS_length
        density = counts.sum() / float(num_codons)
        normalized = counts / float(density)
        
        start_slice = ('start_codon', slice(-codon_buffer, counts.CDS_length + codon_buffer))
        sum_of_normalized_from_start[start_slice] += normalized[start_slice]
        long_enough_genes_from_start[start_slice] += uniform[start_slice]
        
        end_slice = ('stop_codon', slice(-(counts.CDS_length + codon_buffer), codon_buffer))
        sum_of_normalized_from_end[end_slice] += normalized[end_slice]
        long_enough_genes_from_end[end_slice] += uniform[end_slice]
        
    mean_densities = {'from_start': {'codons': sum_of_normalized_from_start / long_enough_genes_from_start},
                      'from_end': {'codons': sum_of_normalized_from_end / long_enough_genes_from_end},
                     }

    return mean_densities

def normalized_codon_density_distribution(codon_counts, offset_key='relaxed'):
    all_normalized_densities = []
    for name, counts_offset_groups in codon_counts.iteritems():
        counts = counts_offset_groups[offset_key]
        num_codons = counts.CDS_length
        normalization = float(counts['start_codon', 0:num_codons].sum()) / num_codons
        if normalization != 0:
            all_normalized_densities.append(np.true_divide(counts['start_codon', 0:10], normalization))

    return np.asarray(all_normalized_densities)

def compute_metacodon_counts(read_positions, gtf_fn, genome_dir, codon_table=1):
    left_buffer = 50
    right_buffer = 50

    # Figure out what lengths were recorded.
    random_gene = read_positions.iterkeys().next()
    length_keys = read_positions[random_gene].keys()

    dinucleotides =  list(''.join(pair) for pair in product('TCAG', repeat=2))
    features_keys = codons.non_stop_codons + ['T', 'C', 'A', 'G'] + dinucleotides

    metacodon_counts = {features_key: {length: PositionCounts({'feature': 0}, left_buffer, right_buffer)
                                       for length in length_keys
                                      }
                        for features_key in features_keys}

    coding_sequence_fetcher = gtf.make_coding_sequence_fetcher(gtf_fn, genome_dir, codon_table)

    for name, read_counts in read_positions.iteritems():
        coding_sequence = coding_sequence_fetcher(name)
        if coding_sequence == None:
            continue

        for c, codon_id in enumerate(codons.codons_from_seq(coding_sequence)):
            p = 3 * c
            if p >= left_buffer and p <= len(coding_sequence) - right_buffer:
                p_slice = ('start_codon', slice(p - left_buffer, p + left_buffer))
                for length in length_keys:
                    actual_counts = read_positions[name][length][p_slice]
                    
                    feature_slice = ('feature', slice(-left_buffer, right_buffer))
                    metacodon_counts[codon_id][length][feature_slice] += actual_counts
                    nucleotide = codon_id[0]
                    metacodon_counts[nucleotide][length][feature_slice] += actual_counts
                    dinucleotide = codon_id[:2]
                    metacodon_counts[dinucleotide][length][feature_slice] += actual_counts

    return metacodon_counts

def get_total_read_count(position_counts, exclude_from_start, exclude_from_end):
    # TODO: this needs to be updated to new offset handling strategy
    A_site_offset = 15
    start_index = -A_site_offset + exclude_from_start
    end_index = position_counts['all'].CDS_length - A_site_offset + 3 - exclude_from_end
    count = position_counts['all']['start_codon', start_index:end_index].sum()
    return count

def compute_read_counts(read_positions, exclude_from_start, exclude_from_end):
    read_counts = {}
    for name, position_counts in read_positions.iteritems():
        read_count = get_total_read_count(position_counts, exclude_from_start, exclude_from_end)
        
        expression = np.array([read_count, 0])
        read_counts[name] = {'CDS_length': position_counts['all'].CDS_length - exclude_from_start - exclude_from_end,
                             'expression': expression,
                            }
    return read_counts

def compute_RPKMs(gene_infos, exclude_from_start, exclude_from_end):
    RPKMs = {}
    total_mapped_reads = 0
    for gene_name in gene_infos:
        read_count = gene_infos[gene_name]['expression'][0]
        length = gene_infos[gene_name]['CDS_length'] - exclude_from_start - exclude_from_end
        RPKMs[gene_name] = float(read_count) / length
        total_mapped_reads += read_count

    for gene_name in RPKMs:
        RPKMs[gene_name] = 1.e9 * RPKMs[gene_name] / total_mapped_reads

    return RPKMs
