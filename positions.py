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

        try:
            leftmost_name = min(landmarks, key=landmarks.get)
            leftmost_value = landmarks[leftmost_name]
            rightmost_name = max(landmarks, key=landmarks.get)
            rightmost_value = landmarks[rightmost_name]
        except AttributeError:
            print landmarks, left_buffer, right_buffer
            raise

        # For historical reasons, rightmost_value is assumed to point to one
        # past the 'interesting' region, so that it is the first index of the
        # left buffer. This means that there is NOT a '+ 1' at the end of this
        # length expression as might be expected.
        length = (rightmost_value + right_buffer) - (leftmost_value - left_buffer)
        self.landmark_to_index = {name: left_buffer + (landmarks[name] - leftmost_value)
                                  for name in landmarks}
       
        if isinstance(dtype, tuple):
            # 2D array
            first_dtype, second_dtype = dtype
            self.data = [PositionCounts(landmarks, left_buffer, right_buffer) for i in range(length)]
        elif data == None:
            self.data = np.zeros(length, dtype)
        else:
            if len(data) != length:
                raise ValueError(len(data), length)
            self.data = data

    def adjust_relative_to_landmark(self, landmark, key):
        if isinstance(key, (int, long)):
            adjusted_key = self.landmark_to_index[landmark] + key
            # Note that adjusted_key must be allowed to be equal to
            # len(self.data) to allow a slice to contain the end point.
            if adjusted_key < 0 or adjusted_key > len(self.data):
                raise IndexError(len(self.data), adjusted_key, key, self.landmark_to_index, self.left_buffer, self.right_buffer)

        elif isinstance(key, slice):
            if key.start == None:
                start = 0
            else:
                start = key.start

            if key.stop == None:
                raise ValueError('None not allowed as stop in slice')
            else:
                stop = key.stop

            if key.step == None:
                adjusted_step = 1
            elif key.step < 0:
                raise ValueError('Negative step not allowed')
            else:
                adjusted_step = key.step

            adjusted_start = self.adjust_relative_to_landmark(landmark, start)
            adjusted_stop = self.adjust_relative_to_landmark(landmark, stop)
            adjusted_key = slice(adjusted_start, adjusted_stop, adjusted_step)

        elif isinstance(key, (list, np.ndarray)):
            adjusted_key = [self.adjust_relative_to_landmark(landmark, k) for k in key]

        else:
            raise TypeError(type(key))

        return adjusted_key
    
    def __getitem__(self, landmarks_and_keys):
        if len(landmarks_and_keys) == 2:
            # 1D
            landmark, key = landmarks_and_keys
            adjusted_key = self.adjust_relative_to_landmark(landmark, key)
            return self.data[adjusted_key]
        elif len(landmarks_and_keys) == 4:
            # 2D
            first_landark, first_key, second_landmark, second_key = landmarks_and_keys
            first_adjusted_key = self.adjust_relative_to_landmark(first_landark, first_key)
            return self.data[first_adjusted_key][second_landmark, second_key]

    def __setitem__(self, landmarks_and_keys, value):
        if len(landmarks_and_keys) == 2:
            # 1D
            landmark, key = landmarks_and_keys
            adjusted_key = self.adjust_relative_to_landmark(landmark, key)
            self.data[adjusted_key] = value
        elif len(landmarks_and_keys) == 4:
            # 2D
            first_landark, first_key, second_landmark, second_key = landmarks_and_keys
            first_adjusted_key = self.adjust_relative_to_landmark(first_landark, first_key)
            self.data[first_adjusted_key][second_landmark, second_key] = value

    def __iadd__(self, other):
        if self.landmarks != other.landmarks:
            raise ValueError('Unequal landmarks:', self.landmarks, other.landmarks)
        if self.left_buffer != other.left_buffer or self.right_buffer != other.right_buffer:
            raise ValueError('Unequal buffer lengths')
        else:
            self.data += other.data

        return self

    def sum(self):
        return self.data.sum()

    @property
    def CDS_length(self):
        return self.landmarks['stop_codon'] - self.landmarks['start_codon']
    
    def __div__(self, other):
        if isinstance(other, PositionCounts):
            if self.landmarks != other.landmarks:
                raise ValueError('Unequal landmarks:', self.landmarks, other.landmarks)
            if self.left_buffer != other.left_buffer or self.right_buffer != other.right_buffer:
                raise ValueError('Unequal buffer lengths')
            else:
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
    
    def __sub__(self, other):
        if isinstance(other, PositionCounts):
            if self.landmarks != other.landmarks:
                raise ValueError('Unequal landmarks:', self.landmarks, other.landmarks)
            if self.left_buffer != other.left_buffer or self.right_buffer != other.right_buffer:
                raise ValueError('Unequal buffer lengths')
            else:
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
        converted_data = np.concatenate((np.zeros(length - 1), position_counts.data[:-(length - 1)]))
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

def get_Transcript_position_counts(clean_bam_fn, transcripts, relevant_lengths, left_buffer=left_buffer, right_buffer=right_buffer):
    gene_infos = {}
    bam_file = pysam.Samfile(clean_bam_fn)
    for transcript in transcripts:
        transcript.build_coordinate_maps(left_buffer, right_buffer)
        if transcript.CDS_length < 0:
            print transcript.name, transcript.CDS_length

        expression = np.zeros(2, int)
        nonunique = 0
        alternatively_spliced = 0

        landmarks = {'start': 0,
                     'start_codon': transcript.transcript_start_codon,
                     'stop_codon': transcript.transcript_stop_codon,
                     'end': transcript.transcript_length,
                    }
        transcript_position_counts = {l: PositionCounts(landmarks, left_buffer, right_buffer)
                                      for l in relevant_lengths + ['all']}
        
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
                expression[1] += 1
                continue

            expression[0] += 1

            if read_strand == '+':
                five_prime_position = read.pos
            elif read_strand == '-':
                five_prime_position = read.aend - 1

            if five_prime_position in transcript.genomic_to_transcript:
                transcript_coord = transcript.genomic_to_transcript[five_prime_position]
                transcript_position_counts['all']['start', transcript_coord] += 1
                if read.qlen in relevant_lengths:
                    transcript_position_counts[read.qlen]['start', transcript_coord] += 1

        gene_infos[transcript.name] = {'CDS_length': transcript.CDS_length,
                                       'position_counts': transcript_position_counts,
                                       'expression': expression,
                                       'nonunique': nonunique,
                                       'alternatively_spliced': alternatively_spliced,
                                      }
        transcript.delete_coordinate_maps()

    return gene_infos

def get_Transcript_polyA_position_counts(extended_bam_fn, transcripts, max_relevant_length, left_buffer=left_buffer, right_buffer=right_buffer):
    gene_infos = {}
    bam_file = pysam.Samfile(extended_bam_fn)
    for transcript in transcripts:
        transcript.build_coordinate_maps(left_buffer, right_buffer)

        landmarks = {'start': 0,
                     'start_codon': transcript.transcript_start_codon,
                     'stop_codon': transcript.transcript_stop_codon,
                     'end': transcript.transcript_length,
                    }
        transcript_position_counts = {l: PositionCounts(landmarks, left_buffer, right_buffer)
                                      for l in range(max_relevant_length + 1) + ['all']}
        
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
                three_prime_position = read.aend - 1
            elif read_strand == '-':
                three_prime_position = read.pos

            if three_prime_position in transcript.genomic_to_transcript:
                transcript_coord = transcript.genomic_to_transcript[three_prime_position]
                transcript_position_counts['all']['start', transcript_coord] += 1
                nongenomic_length = trim.get_nongenomic_length(read)
                if nongenomic_length <= max_relevant_length:
                    transcript_position_counts[nongenomic_length]['start', transcript_coord] += 1

        gene_infos[transcript.name] = {'CDS_length': transcript.CDS_length,
                                       'position_counts': transcript_position_counts,
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
    codon_counts = PositionCounts({'start_codon': 0, 'stop_codon': num_codons},
                                  codon_buffer,
                                  codon_buffer,
                                 )

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

    return codon_counts

def compute_metagene_positions(position_counts, max_CDS_length):
    ''' max_CDS_length needs to be passed in because it may reflect the max of
        more than just the CDS being considered here.
    '''
    random_gene = position_counts.itervalues().next()
    relevant_lengths = random_gene.keys()
    random_length = relevant_lengths[0]
    random_position_counts = random_gene[random_length]
    left_buffer = random_position_counts.left_buffer
    right_buffer = random_position_counts.right_buffer

    landmarks = {'start_codon': 0,
                 'stop_codon': max_CDS_length,
                }
    from_starts = {length: PositionCounts(landmarks, left_buffer, right_buffer)
                   for length in relevant_lengths}
    from_ends = {length: PositionCounts(landmarks, left_buffer, right_buffer)
                 for length in relevant_lengths}
    
    for name, counts in position_counts.iteritems():
        CDS_length = counts[random_length].CDS_length
        start_slice = ('start_codon', slice(-left_buffer, CDS_length))
        end_slice = ('stop_codon', slice(-CDS_length, right_buffer))
        for length in relevant_lengths:
            from_starts[length][start_slice] += counts[length][start_slice]
            from_ends[length][end_slice] += counts[length][end_slice]

    from_starts_and_ends = {'from_starts': from_starts,
                            'from_ends': from_ends,
                           }
    return from_starts_and_ends

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

    for name, counts in codon_counts.iteritems():
        if name in names_to_skip:
            print 'skipping', name
            continue

        counts = counts[offset_key]
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

def compute_metacodon_counts(read_positions, gtf_fn, genome_dir):
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

    coding_sequence_fetcher = gtf.make_coding_sequence_fetcher(gtf_fn, genome_dir)

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
