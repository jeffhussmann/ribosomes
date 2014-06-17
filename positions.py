''' Utilities for counting reads mapped to each position in a gene. '''

import numpy as np
import numbers
import pysam
from collections import Counter
from itertools import cycle, product
import codons
import gtf
import Sequencing.Serialize

class OldPositionCounts(object):
    ''' Wrapper around an array of counts of positions for an extent and for a
        buffer on either edge that allows for indexing relative to the extent's
        start and end.
    '''
    def __init__(self, extent_length, left_buffer, right_buffer, counts=None, dtype=int):
        self.left_buffer = left_buffer
        self.right_buffer = right_buffer
        self.extent_length = extent_length
        if counts == None:
            self.counts = np.zeros(left_buffer + right_buffer + extent_length, dtype=dtype)
        else:
            assert len(counts) == left_buffer + right_buffer + extent_length
            self.counts = counts

    @classmethod
    def from_string(cls, string, left_buffer, right_buffer):
        try:
            counts = np.array(map(int, string.strip().split()))
            dtype = int
        except ValueError:
            counts = np.array(map(float, string.strip().split()))
            dtype = float

        extent_length = len(counts) - left_buffer - right_buffer
        return cls(extent_length, left_buffer, right_buffer, counts=counts, dtype=dtype)

    def __str__(self):
        string = ' '.join(map(str, self.counts)) + '\n'
        return string

    def sum(self):
        return self.counts.sum()

    def adjust_relative_to_start(self, key):
        if isinstance(key, (int, long)):
            adjusted_key = self.left_buffer + key
            if adjusted_key < 0:
                raise IndexError(adjusted_key, key, self.left_buffer)

        elif isinstance(key, slice):
            if key.start == None:
                start = 0
            else:
                start = key.start

            if key.stop == None:
                stop = self.extent_length
            else:
                stop = key.stop

            if key.step == None:
                adjusted_step = 1
            elif key.step < 0:
                raise ValueError('Negative step not allowed')
            else:
                adjusted_step = key.step

            adjusted_start = self.adjust_relative_to_start(start)
            adjusted_stop = self.adjust_relative_to_start(stop)
            adjusted_key = slice(adjusted_start, adjusted_stop, adjusted_step)

        elif isinstance(key, (list, np.ndarray)):
            adjusted_key = map(self.adjust_relative_to_start, key)

        else:
            raise TypeError(type(key))

        return adjusted_key

    def __getitem__(self, key):
        adjusted_key = self.adjust_relative_to_start(key)
        return self.counts[adjusted_key]

    def __setitem__(self, key, value):
        adjusted_key = self.adjust_relative_to_start(key)
        self.counts[adjusted_key] = value
    
    def __iadd__(self, other):
        assert self.extent_length == other.extent_length
        assert self.left_buffer == other.left_buffer
        assert self.right_buffer == other.right_buffer

        self.counts += other.counts

        return self

    def __div__(self, other):
        if isinstance(other, PositionCounts):
            if self.extent_length != other.extent_length:
                raise ValueError
            if self.left_buffer != other.left_buffer:
                raise ValueError
            if self.right_buffer != other.right_buffer:
                raise ValueError

            return PositionCounts(self.extent_length, self.left_buffer, self.right_buffer, counts=np.true_divide(self.counts, other.counts))
        elif isinstance(other, numbers.Number):
            return PositionCounts(self.extent_length, self.left_buffer, self.right_buffer, counts=np.true_divide(self.counts, other))

    @property
    def relative_to_end(self):
        return PositionCounts.RelativeToEndCounts(self)

    class RelativeToEndCounts(object):
        ''' Convoluted hack to allow PositionCounts.relative_to_end to be an object
            that has __getitem__ and __setitem__ methods so that
            position_counts.relative_to_end[a:b:c] works.
        '''
        def __init__(self, position_counts):
            self.right_buffer = position_counts.right_buffer
            self.extent_length = position_counts.extent_length
            self.counts = position_counts.counts

        def adjust_relative_to_end(self, key):
            ''' Note: relative_to_end[0] will point to one after the end of the
                extent.
            '''
            if isinstance(key, (int, long)):
                adjusted_key = len(self.counts) - self.right_buffer - key
                if adjusted_key == -1:
                    # Slicing quirk to include the last element
                    adjusted_key = None
                elif adjusted_key < -1:
                    raise IndexError

            elif isinstance(key, slice):
                if key.start == None:
                    start = 0
                else:
                    start = key.start

                if key.stop == None:
                    stop = self.extent_length
                else:
                    stop = key.stop

                if key.step == None:
                    adjusted_step = -1
                elif key.step < 0:
                    raise ValueError('Negative step not allowed')
                else:
                    adjusted_step = -key.step

                adjusted_start = self.adjust_relative_to_end(start)
                adjusted_stop = self.adjust_relative_to_end(stop)
                adjusted_key = slice(adjusted_start, adjusted_stop, adjusted_step)

            elif isinstance(key, (list, np.ndarray)):
                adjusted_key = map(self.adjust_relative_to_end, key)

            else:
                raise TypeError(type(key))

            return adjusted_key

        def __getitem__(self, key):
            adjusted_key = self.adjust_relative_to_end(key)
            return self.counts[adjusted_key]

        def __setitem__(self, key, value):
            adjusted_key = self.adjust_relative_to_end(key)
            self.counts[adjusted_key] = value


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
        # left buffer. This means that this is NOT a '+ 1' at the end of this
        # length expression as might be expected.
        length = (rightmost_value + right_buffer) - (leftmost_value - left_buffer)
        self.landmark_to_index = {name: left_buffer + (landmarks[name] - leftmost_value)
                                  for name in landmarks}
        
        if data == None:
            self.data = np.zeros(length, dtype)
        else:
            if len(data) != length:
                raise ValueError(len(data), length)
            else:
                self.data = data

    def adjust_relative_to_landmark(self, landmark, key):
        if isinstance(key, (int, long)):
            adjusted_key = self.landmark_to_index[landmark] + key
            # Note that adjusted_key must be allowed to be equal to
            # len(self.data) to allow a slice to contain the end point.
            if adjusted_key < 0 or adjusted_key > len(self.data):
                raise IndexError(adjusted_key, key, self.landmark_to_index, self.left_buffer, self.right_buffer)

        elif isinstance(key, str):
            if key not in self.landmark_to_index:
                raise ValueError(key)

            adjusted_key = self.landmark_to_index[key]

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
    
    def __getitem__(self, landmark_and_key):
        landmark, key = landmark_and_key
        adjusted_key = self.adjust_relative_to_landmark(landmark, key)
        return self.data[adjusted_key]

    def __setitem__(self, landmark_and_key, value):
        landmark, key = landmark_and_key
        adjusted_key = self.adjust_relative_to_landmark(landmark, key)
        self.data[adjusted_key] = value
    
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

# Number of nucleotide positions to include on the left and right side of
# counts of read positions.
edge_buffer = 50
left_buffer = 4 * edge_buffer
right_buffer = 4 * edge_buffer

# Number of codons to include on either side of counts of codon positions.
codon_buffer = 10

def get_Transcript_position_counts(clean_bam_fn, transcripts, relevant_lengths):
    gene_infos = {}
    bam_file = pysam.Samfile(clean_bam_fn)
    for transcript in transcripts:
        transcript.build_coordinate_maps()
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
            if any(position not in transcript.genomic_to_transcript for position in read.positions):
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
        landmarks = counts[random_length].landmarks
        CDS_length = landmarks['stop_codon'] - landmarks['start_codon']
        start_slice = ('start_codon', slice(-left_buffer, CDS_length))
        end_slice = ('stop_codon', slice(-CDS_length, right_buffer))
        for length in relevant_lengths:
            from_starts[length][start_slice] += counts[length][start_slice]
            from_ends[length][end_slice] += counts[length][end_slice]

    from_starts_and_ends = {'from_starts': from_starts,
                            'from_ends': from_ends,
                           }
    return from_starts_and_ends

def compute_metacodon_counts(codon_counts, gtf_fn, genome_dir):
    window = 30

    landmarks = {'codon': 0}
    metacodon_counts = {codon_id: {'actual': PositionCounts(landmarks, window, window),
                                   'uniform': PositionCounts(landmarks, window, window, dtype=float),
                                   'sum_of_enrichments': PositionCounts(landmarks, window, window, dtype=float),
                                   'num_eligible': PositionCounts(landmarks, window, window),
                                  }
                        for codon_id in codons.non_stop_codons}

    coding_sequence_fetcher = gtf.make_coding_sequence_fetcher(gtf_fn, genome_dir)

    for name, read_counts in codon_counts.iteritems():
        coding_sequence = coding_sequence_fetcher(name)
        if coding_sequence == None:
            continue

        total_counts = codon_counts[name].sum()
        num_codons = len(codon_counts[name])
        density = float(total_counts) / num_codons
        uniform_counts = np.full(2 * window, density)
        num_eligible = np.ones(2 * window)

        for p, codon_id in enumerate(codons.codons_from_seq(coding_sequence)):
            if p >= window and p <= num_codons - window:
                actual_counts = read_counts[p - window:p + window]
                metacodon_counts[codon_id]['actual']['codon', -window:window] += actual_counts
                metacodon_counts[codon_id]['uniform']['codon', -window:window] += uniform_counts
                
                if density > 0:
                    enrichments = actual_counts / density
                    metacodon_counts[codon_id]['sum_of_enrichments']['codon', -window:window] += enrichments
                    metacodon_counts[codon_id]['num_eligible']['codon', -window:window] += num_eligible

    return metacodon_counts

def compute_metacodon_counts_nucleotide_resolution(read_positions, gtf_fn, genome_dir):
    minus_window = 40
    plus_window = 40

    random_gene = read_positions.iterkeys().next()
    length_keys = read_positions[random_gene].keys()

    dinucleotides =  list(''.join(pair) for pair in itertools.product('TCAG', repeat=2))
    keys = codons.non_stop_codons + ['T', 'C', 'A', 'G'] + dinucleotides

    metacodon_counts = {key: {length: PositionCounts(plus_window, minus_window, 0)
                              for length in length_keys
                             }
                        for key in keys}

    coding_sequence_fetcher = gtf.make_coding_sequence_fetcher(gtf_fn, genome_dir)

    for name, read_counts in read_positions.iteritems():
        coding_sequence = coding_sequence_fetcher(name)
        if coding_sequence == None:
            continue

        for c, codon_id in enumerate(codons.codons_from_seq(coding_sequence)):
            p = 3 * c
            if p >= minus_window and p <= len(coding_sequence) - plus_window:
                for length in length_keys:
                    actual_counts = read_positions[name][length][p - minus_window:p + plus_window]
                    
                    metacodon_counts[codon_id][length][-minus_window:plus_window] += actual_counts

                    metacodon_counts[codon_id[0]][length][-minus_window:plus_window] += actual_counts
                    
                    #dinucleotide = codon_id[:2]
                    dinucleotide = coding_sequence[p - 1:p + 1]
                    metacodon_counts[dinucleotide][length][-minus_window:plus_window] += actual_counts

    return metacodon_counts

def compute_averaged_codon_densities(codon_counts, names_to_skip=set()): 
    # To reduce noise, genes with less than min_counts total counts are ignored.
    min_counts = 64
    max_length = max(counts['relaxed'].CDS_length
                     for counts in codon_counts.itervalues()
                     if counts['relaxed'].sum() >= min_counts
                    )

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

        counts = counts['relaxed']
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
