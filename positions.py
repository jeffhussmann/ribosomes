''' Utilities for counting reads mapped to each position in a gene. '''

import numpy as np
import numbers
import pysam
from collections import Counter
from itertools import cycle, product
import codons
import gtf
import Sequencing.Serialize

class PositionCounts(object):
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

        transcript_position_counts = {l: PositionCounts(transcript.transcript_length, left_buffer, right_buffer)
                                      for l in relevant_lengths + ['all']}
        
        # fetch raises a ValueError if given a negative start, but it doesn't 
        # care if the end is valid.
        left_edge = max(0, transcript.start - left_buffer)
        right_edge = transcript.end + right_buffer
        overlapping_reads = bam_file.fetch(transcript.seqname, left_edge, right_edge)
        for read in overlapping_reads:
            if any(position not in transcript.genomic_to_transcript for position in read.positions):
                # Alternative splicing
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
            transcript_position_counts['all'][transcript_coord] += 1
            if read.qlen in relevant_lengths:
                transcript_position_counts[read.qlen][transcript_coord] += 1

        CDS_position_counts = {l: PositionCounts(transcript.CDS_length, left_buffer, right_buffer)
                               for l in relevant_lengths + ['all']}

        for key in CDS_position_counts:
            CDS_slice = slice(transcript.transcript_start_codon - left_buffer, transcript.transcript_stop_codon + right_buffer)
            # [:] is to cause an error if the lengths aren't the same.
            CDS_position_counts[key].counts[:] = transcript_position_counts[key][CDS_slice]
        
        gene_infos[transcript.name] = {'CDS_length': transcript.CDS_length,
                                       'position_counts': CDS_position_counts,
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
    arbitrary_length_counts = position_counts.values()[0]
    CDS_length = arbitrary_length_counts.extent_length
    if CDS_length % 3 != 0:
        raise ValueError('CDS length not divisible by 3')

    # Note: CDS_length is the index of the first nucleotide of the stop codon.
    # Ingolia's original model never has the stop codon in the A site, but
    # subsequent data show an accumulation of (typically length 29 or 30) reads
    # that do advance this far.
    num_codons = CDS_length // 3
    codon_counts = PositionCounts(num_codons, codon_buffer, codon_buffer)

    recorded_lengths = set(position_counts.keys())
    known_A_site_lengths = set(A_site_offsets[offset_type].keys())

    for length in recorded_lengths & known_A_site_lengths:
        A_site_offset = A_site_offsets[offset_type][length]
        start_index = -A_site_offset - (codon_buffer * 3)
        end_index = CDS_length - A_site_offset + (codon_buffer * 3)
        in_frame = slice(start_index, end_index, 3)
        one_behind = slice(start_index - 1, end_index - 1, 3)
        one_ahead = slice(start_index + 1, end_index + 1, 3)
        codon_counts.counts += position_counts[length][in_frame] + \
                               position_counts[length][one_behind] + \
                               position_counts[length][one_ahead]

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

    from_starts = {length: PositionCounts(max_CDS_length, left_buffer, right_buffer)
                   for length in relevant_lengths}
    from_ends = {length: PositionCounts(max_CDS_length, left_buffer, right_buffer)
                 for length in relevant_lengths}
    
    for name, counts in position_counts.iteritems():
        CDS_length = counts[random_length].extent_length
        start_slice = slice(-left_buffer, CDS_length)
        # This is harmlessly wrong, should be in some sense -right_buffer + 1
        end_slice = slice(-right_buffer, CDS_length)
        for length in relevant_lengths:
            from_starts[length][start_slice] += counts[length][start_slice]
            from_ends[length].relative_to_end[end_slice] += counts[length].relative_to_end[end_slice]

    from_starts_and_ends = {'from_starts': from_starts,
                            'from_ends': from_ends,
                           }
    return from_starts_and_ends

def compute_metacodon_counts(codon_counts, gtf_fn, genome_dir):
    window = 30

    metacodon_counts = {codon_id: {'actual': PositionCounts(window, window, 0),
                                   'uniform': PositionCounts(window, window, 0, dtype=float),
                                   'sum_of_enrichments': PositionCounts(window, window, 0, dtype=float),
                                   'num_eligible': PositionCounts(window, window, 0),
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
                metacodon_counts[codon_id]['actual'][-window:window] += actual_counts
                metacodon_counts[codon_id]['uniform'][-window:window] += uniform_counts
                
                if density > 0:
                    enrichments = actual_counts / density
                    metacodon_counts[codon_id]['sum_of_enrichments'][-window:window] += enrichments
                    metacodon_counts[codon_id]['num_eligible'][-window:window] += num_eligible

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
    max_length = max(counts['relaxed'].extent_length
                     for counts in codon_counts.itervalues()
                     if counts['relaxed'].sum() >= min_counts
                    )

    sum_of_normalized_from_start = PositionCounts(max_length, codon_buffer, codon_buffer, dtype=float)
    sum_of_normalized_from_end = PositionCounts(max_length, codon_buffer, codon_buffer, dtype=float)
    long_enough_genes_from_start = PositionCounts(max_length, codon_buffer, codon_buffer)
    long_enough_genes_from_end = PositionCounts(max_length, codon_buffer, codon_buffer)

    uniform = PositionCounts(max_length, codon_buffer, codon_buffer, counts=np.ones(max_length + 2 * codon_buffer))

    for name, counts in codon_counts.iteritems():
        if name in names_to_skip:
            print 'skipping', name
            continue

        counts = counts['relaxed']
        if counts.sum() < min_counts:
            continue

        num_codons = counts.extent_length + 2 * codon_buffer
        density = counts.sum() / float(num_codons)
        normalized = counts / float(density)
        
        start_slice = slice(-codon_buffer, counts.extent_length + codon_buffer)
        sum_of_normalized_from_start[start_slice] += normalized[start_slice]
        long_enough_genes_from_start[start_slice] += uniform[start_slice]
        
        end_slice = slice(-codon_buffer + 1, counts.extent_length + 1 + codon_buffer)
        sum_of_normalized_from_end.relative_to_end[end_slice] += normalized.relative_to_end[end_slice]
        long_enough_genes_from_end.relative_to_end[end_slice] += uniform.relative_to_end[end_slice]
        
    mean_densities = {'from_start': {'codons': sum_of_normalized_from_start / long_enough_genes_from_start},
                      'from_end': {'codons': sum_of_normalized_from_end / long_enough_genes_from_end},
                     }

    return mean_densities

def get_total_read_count(position_counts, exclude_from_start, exclude_from_end):
    # TODO: this needs to be updated to new offset handling strategy
    CDS_length = position_counts['all'].extent_length
    A_site_offset = 15
    start_index = -A_site_offset + exclude_from_start
    end_index = CDS_length - A_site_offset + 3 - exclude_from_end
    count = position_counts['all'][start_index:end_index].sum()
    return count

def compute_read_counts(read_positions, exclude_from_start, exclude_from_end):
    read_counts = {}
    for name, position_counts in read_positions.iteritems():
        read_count = get_total_read_count(position_counts, exclude_from_start, exclude_from_end)
        
        expression = np.array([read_count, 0])
        read_counts[name] = {'CDS_length': position_counts['all'].extent_length,
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
