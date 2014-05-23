''' Utilities for counting reads mapped to each position in a gene. '''

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg', warn=False)
import numpy as np
import numbers
import itertools
import pysam
from collections import Counter
import codons
import gtf
import Circles.Serialize
import brewer2mpl

def smoothed(array, window_size):
    smoothed_array = np.zeros_like(array)
    for i in range(window_size):
        smoothed_array[i] = array[:i + 1].sum() / float(i + 1)
    for i in range(window_size, len(array) - window_size):
        smoothed_array[i] = array[i - window_size:i + window_size + 1].sum() / float(2 * window_size + 1)
    for i in range(len(array) - window_size, len(array)):
        smoothed_array[i] = array[i:].sum() / float(len(array) - i)
    return smoothed_array

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
left_buffer = 2 * edge_buffer
right_buffer = edge_buffer

# Number of codons to include on either side of counts of codon positions.
codon_buffer = 10

def get_CDS_position_counts(clean_bam_fn, simple_CDSs, relevant_lengths):
    gene_infos = {}
    bam_file = pysam.Samfile(clean_bam_fn)
    for CDS in simple_CDSs:
        gene_name = CDS.attribute['transcript_id']
        extent = (CDS.seqname, CDS.strand, CDS.start, CDS.end)
        CDS_length = abs(CDS.end - CDS.start) + 1  
        gene_infos[gene_name] = {'CDS_length': CDS_length}

        try:
            extent_info = get_extent_position_counts(bam_file, extent, relevant_lengths)
            gene_infos[gene_name].update(extent_info)
        except IndexError:
            # For now, this means there were spliced reads in the extent
            print gene_name
            del gene_infos[gene_name]

    return gene_infos

def get_extent_position_counts(bam_file, extent, relevant_lengths):
    ''' Returns counts of reads whose 5' edge falls at each position in extent,
        stratified by length for each length in relevant_lengths, as well as
        cumulative counts for all lengths.
    '''
    seqname, extent_strand, start, end = extent
    extent_length = abs(end - start) + 1
    position_counts = {l: PositionCounts(extent_length, left_buffer, right_buffer)
                       for l in relevant_lengths + ['all']}
    expression = np.zeros(2, int)
    nonunique = 0
    
    region_start = start - edge_buffer
    region_end = end + edge_buffer
    overlapping_reads = bam_file.fetch(seqname, region_start, region_end)
    for read in overlapping_reads:
        if read.mapq != 50:
            nonunique += 1
        else:
            read_strand = '-' if read.is_reverse else '+'

            if read_strand != extent_strand:
                expression[1] += 1
                continue
            else:
                expression[0] += 1
                if extent_strand == '+':
                    from_start = read.pos - start
                elif extent_strand == '-':
                    from_start = end - (read.aend - 1)

            try:
                position_counts['all'][from_start] += 1
            except IndexError:
                raise IndexError
            if read.qlen in relevant_lengths:
                position_counts[read.qlen][from_start] += 1

    extent_info = {'position_counts': position_counts,
                   'expression': expression,
                   'nonunique': nonunique,
                  }

    return extent_info

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
        
        overlapping_reads = bam_file.fetch(transcript.seqname, transcript.start, transcript.end)
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

    return from_starts, from_ends

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

def make_codon_counts_file(codon_counts, file_name, stringency='relaxed'):
    with open(file_name, 'w') as fh:
        for name in sorted(codon_counts):
            buffered_counts = codon_counts[name][stringency]
            num_codons = buffered_counts.extent_length
            # + 1 is to include the stop codon
            counts = buffered_counts[:num_codons + 1]
            counts_string = '\t'.join(str(count) for count in counts)
            line = '{0}\t{1}\n'.format(name, counts_string)
            fh.write(line)

def read_codon_counts_file(fn, reports_P_site=False):
    genes = {}
    for line in open(fn):
        name, values = line.strip().split('\t', 1)
        values = np.fromstring(values, dtype=float, sep='\t')
        # Weinberg's file reports counts in P-sites. We are interested in counts
        # in A-sites, for which no reliable information can be obtained for the
        # initiation codon or first codon after this.
        # Shift everything over by 1 here to report A-sites.
        # Responsibility of analysis to ignore first two values and last value.
        if reports_P_site:
            values = np.concatenate(([0], values))
        genes[name] = values
            
    return genes

def plot_metagene_positions(from_starts, from_ends, figure_fn):
    relevant_lengths = sorted(from_starts.keys())
    fig, (start_ax, end_ax) = plt.subplots(2, 1, figsize=(12, 16))

    start_xs = np.arange(-21, 19)
    end_xs = np.arange(-6, 34)
    for length in relevant_lengths:
        start_ax.plot(start_xs, from_starts[length][start_xs], '.-', label=length)
        end_ax.plot(-end_xs, from_ends[length].relative_to_end[end_xs], '.-', label=length)

    start_ax.set_xlim(min(start_xs), max(start_xs))
    mod_3 = [x for x in start_xs if x % 3 == 0]
    start_ax.set_xticks(mod_3)
    for x in mod_3:
        start_ax.axvline(x, color='black', alpha=0.1)
    start_ax.set_xlabel('Position of read relative to start of CDS')
    start_ax.set_ylabel('Number of uniquely mapped reads of specified length')
    
    end_ax.set_xlim(min(-end_xs), max(-end_xs))
    mod_3 = [x for x in -end_xs if x % 3 == 0]
    end_ax.set_xticks(mod_3)
    for x in mod_3:
        end_ax.axvline(x, color='black', alpha=0.1)
    end_ax.set_xlabel('Position of read relative to stop codon')
    end_ax.set_ylabel('Number of uniquely mapped reads of specified length')

    start_ax.legend(loc='upper right', framealpha=0.5)
    end_ax.legend(loc='upper right', framealpha=0.5)
    
    fig.savefig(figure_fn)

def plot_frames(from_starts, figure_fn):
    ''' Plots bar graphs of the distribution of reading frame that the 5' end of
        reads fell, stratified by read length.
    '''
    def bar_plot(values, name, tick_labels, ax):
        N = len(values)

        starts = np.arange(N)
        gap = 1 # in multiples of bar width

        width = 1. / (1 + gap)

        ax.bar(starts, values, width, label=name)

        ax.set_xlim(xmin=-width / 2, xmax=N - width / 2)

        ax.set_xticks(starts + (width * 1 / 2))
        ax.set_xticklabels(tick_labels)
        ax.xaxis.set_ticks_position('none')
        ax.set_ylim(0, 1)

        leg = ax.legend(loc='upper right', fancybox=True)
        leg.get_frame().set_alpha(0.5)

    relevant_lengths = sorted(from_starts.keys())
    extent_length = from_starts[relevant_lengths[0]].extent_length
    fig, axs = plt.subplots(len(relevant_lengths), 1, figsize=(6, 3 * len(relevant_lengths)))

    for length, ax in zip(relevant_lengths, axs):
        frames = np.zeros(3, int)
        for p in range(-15, extent_length):
            frames[p % 3] += from_starts[length][p]

        frames = np.true_divide(frames, frames.sum())

        tick_labels = ['0', '1', '2']

        bar_plot(frames, str(length), tick_labels, ax)
        
    axs[-1].set_xlabel('Frame')
    axs[len(axs) // 2].set_ylabel('Fraction of uniquely mapped reads of specified length')

    fig.savefig(figure_fn)

def plot_averaged_codon_densities(data_sets, figure_fn, show_end=True, past_edge=codon_buffer, smooth=False):
    bmap = brewer2mpl.get_map('Set1', 'qualitative', 9)
    colors = bmap.mpl_colors[:5] + bmap.mpl_colors[6:] + ['black']
    colors = colors + colors

    plot_up_to = 500
    
    if show_end:
        fig, (start_ax, end_ax) = plt.subplots(1, 2, figsize=(12, 8))
    else:
        fig, start_ax = plt.subplots(figsize=(12, 8))

    start_xs = np.arange(-past_edge, plot_up_to + 20 + 1)
    end_xs = np.arange(-past_edge + 1, plot_up_to + 20 + 1)

    for name, mean_densities, color_index in data_sets:
        densities = mean_densities['from_start']['codons'][0:]
        if smooth:
            densities = smoothed(densities, 5)
        start_densities = densities[start_xs]

        steady_state = mean_densities['from_start']['codons'][400:500].sum() / 100
        #steady_state = 1.

        marker = '' if smooth else '.'
        linewidth = 2 if smooth else 1

        start_ax.plot(start_xs, start_densities / steady_state, '.-', label=name, color=colors[color_index], marker=marker, linewidth=linewidth)
        if show_end:
            end_ax.plot(-end_xs, mean_densities['from_end']['codons'].relative_to_end[end_xs], '.-', label=name, color=colors[color_index])
    
    start_ax.legend(loc='upper right', framealpha=0.5)

    start_ax.set_xlabel('Number of codons from start codon')
    start_ax.set_ylabel('Mean normalized read density')
    start_ax.set_xlim(min(start_xs), max(start_xs))
    start_ax.plot(start_xs, [1 for x in start_xs], color='black', alpha=0.5)
   
    if show_end:
        end_ax.set_xlabel('Number of codons from stop codon')
        end_ax.yaxis.tick_right()
        end_ax.set_xlim(min(-end_xs), max(-end_xs))
        end_ax.plot(-end_xs, [1 for x in end_xs], color='black', alpha=0.5)
    
    axs = [start_ax]
    if show_end:
        axs.append(end_ax)

    ymax = max(ax.get_ylim()[1] for ax in axs)
    
    for ax in axs:
        #ax.set_ylim(0, ymax + 0.1)
        ax.set_ylim(0, 6)
        ax.set_xlim(xmax=plot_up_to)
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        ax.set_aspect((xmax - xmin) / (ymax - ymin))
        #ax.set_aspect(1)

    fig.set_size_inches(12, 12)
    fig.savefig(figure_fn, bbox_inches='tight')

def plot_metacodon_counts(metacodon_counts, fig_fn, codon_ids='all', enrichment=False, keys_to_plot=['actual']):
    random_counts = metacodon_counts['TTT'].itervalues().next()
    xs = np.arange(-random_counts.left_buffer, random_counts.extent_length)

    bases = 'TCAG'

    def make_plot(ax, codon_id):
        ax.set_title(codon_id)

        if codon_id not in metacodon_counts:
            return

        if enrichment:
            average_enrichment = metacodon_counts[codon_id]['sum_of_enrichments'] / metacodon_counts[codon_id]['num_eligible']
            ones = np.ones(len(average_enrichment.counts))
            ax.plot(xs, average_enrichment.counts, 'o-')
            ax.plot(xs, ones, color='black', alpha=0.5)
        else:
            for key in keys_to_plot:
                counts = metacodon_counts[codon_id][key].counts
                line, = ax.plot(xs, counts, '.-', label=key)
                if isinstance(key, int):
                    # -key + 1 is the position a read starts at if position 0
                    # is the last base in it
                    ax.axvline(-key + 1, color=line.get_color(), alpha=0.2)
            #ax.plot(xs, metacodon_counts[codon_id]['uniform'].counts, '.-')
            ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 3))
            if codon_ids != 'all':
                ax.legend(loc='upper right', framealpha=0.5)

        ax.set_ylim(ymin=0)
        ax.set_xlim(min(xs), max(xs))
        ax.axvline(0, ls='--', color='black', alpha=0.5)
    
    if codon_ids == 'all':
        fig, axs = plt.subplots(16, 4, figsize=(16, 64))

        for first in range(4):
            for second in range(4):
                for third in range(4):
                    i = first * 4 + third
                    j = second
                    codon_id = ''.join(bases[first] + bases[second] + bases[third])
                    make_plot(axs[i, j], codon_id)
    else:
        fig, axs = plt.subplots(len(codon_ids), 1, figsize=(8, 8 * len(codon_ids)))
        for codon_id, ax in zip(codon_ids, axs):
            make_plot(ax, codon_id)

    fig.savefig(fig_fn, bbox_inches='tight')

def plot_single_length_metacodon_counts(metacodon_counts, fig_fn, codon_ids, length=28, edge="5'"):
    random_codon_id = metacodon_counts.iterkeys().next()
    random_counts = metacodon_counts[random_codon_id].itervalues().next()
    keys = metacodon_counts[random_codon_id].keys()
    xs = np.arange(-20, 20)

    bases = 'TCAG'

    def make_plot(ax, codon_id):
        ax.set_title(codon_id)

        if codon_id not in metacodon_counts:
            return

        counts = metacodon_counts[codon_id][length]
        if edge == "3'":
            positions = xs - length + 1
        else:
            positions = xs
        line, = ax.plot(xs, counts[positions], '.-')
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 3))

        ax.set_ylim(ymin=0)
        ax.set_xlim(min(xs), max(xs))
        ax.axvline(0, ls='--', color='black', alpha=0.5)
    
    if codon_ids == 'all':
        fig, axs = plt.subplots(16, 4, figsize=(16, 64))

        for first in range(4):
            for second in range(4):
                for third in range(4):
                    i = first * 4 + third
                    j = second
                    codon_id = ''.join(bases[first] + bases[second] + bases[third])
                    make_plot(axs[i, j], codon_id)
    else:
        fig, axs = plt.subplots(len(codon_ids), 1, figsize=(8, 8 * len(codon_ids)))
        for codon_id, ax in zip(codon_ids, axs):
            make_plot(ax, codon_id)

    fig.suptitle('Read counts around every occurence, {0}'.format(edge))
    fig.savefig(fig_fn, bbox_inches='tight')

#if __name__ == '__main__':
#    data_sets = [
#        #('belgium_3_5_14', '/home/jah/projects/arlen/experiments/belgium_3_5_14/wt/results/wt_codon_counts.txt'),
#        #('belgium_8_6_13', '/home/jah/projects/arlen/experiments/belgium_8_6_13/R98S_cDNA_sample/results/R98S_cDNA_sample_codon_counts.txt'),
#        ('guo_nature', '/home/jah/projects/arlen/experiments/guo_nature/Footprint_wild-type_runs1-2/results/guo_nature_codon_counts.txt'),
#        ('ingolia_cell_nothing', '/home/jah/projects/arlen/experiments/ingolia_cell/ES_cell_feeder-free__w__LIF_none_ribo_mesc_nochx_Illumina_GAII/results/ingolia_cell_codon_counts.txt'),
#        ('ingolia_cell_emetine', '/home/jah/projects/arlen/experiments/ingolia_cell/ES_cell_feeder-free__w__LIF_60_s_emetine_ribo_mesc_emet_Illumina_GAII/results/emetine_codon_counts.txt'),
#    ]
#    plot_metagene_averaged(data_sets, from_end=False)

if __name__ == '__main__':
    data_sets = [
        ('Ingolia 1 (CHX)', '/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-1/results/Footprints-rich-1_mean_densities.txt', 3),
        ('Ingolia 2 (CHX)', '/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-2/results/Footprints-rich-2_mean_densities.txt', 0),
        ('McManus (CHX)', '/home/jah/projects/arlen/experiments/mcmanus_gr/S._cerevisiae_Ribo-seq_Rep_1/results/S._cerevisiae_Ribo-seq_Rep_1_mean_densities.txt', 4),
        ('Gerashchenko (CHX)', '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep1_foot/results/Initial_rep1_foot_mean_densities.txt', 5),
        ('Zinshteyn (CHX)', '/home/jah/projects/arlen/experiments/zinshteyn_plos_genetics/WT_Ribosome_Footprint_1/results/WT_Ribosome_Footprint_1_mean_densities.txt', 6),
        ('Dunn (CHX)', '/home/jah/projects/arlen/experiments/dunn_elife/dunn_elife/results/dunn_elife_mean_densities.txt', 1),
        ('Weinberg (flash freezing)', '/home/jah/projects/arlen/experiments/weinberg/RPF/results/RPF_mean_densities.txt', -1),
        ('Arlen 1 (CHX)', '/home/jah/projects/arlen/experiments/belgium_3_5_14/wt/results/wt_mean_densities.txt', 3),
        ('Arlen 2 (CHX)', '/home/jah/projects/arlen/experiments/belgium_8_6_13/WT_cDNA_sample/results/WT_cDNA_sample_mean_densities.txt', 2),
        ##('weinberg_mRNA', '/home/jah/projects/arlen/experiments/weinberg/mRNA/results/mRNA_mean_densities.txt'),
        
        #('Artieri (CHX)', '/home/jah/projects/arlen/experiments/artieri/Mixed_parental_ribosome_protected_fragments_replicate_1/results/artieri_RPF_replicate_1_mean_densities.txt', 2),
        #('Guydosh (flash freezing then CHX)', '/home/jah/projects/arlen/experiments/guydosh_cell/wild-type_CHX/results/wild-type_CHX_mean_densities.txt'),
        
        #('Ingolia mRNA', '/home/jah/projects/arlen/experiments/ingolia_science/mRNA-rich-1/results/mRNA-rich-1_mean_densities.txt', 0),
        
        #('R98S', '/home/jah/projects/arlen/experiments/belgium_8_6_13/R98S_cDNA_sample/results/R98S_cDNA_sample_mean_densities.txt', 2),
        #('Supressed', '/home/jah/projects/arlen/experiments/belgium_8_6_13/Suppressed_R98S_cDNA_sample/results/Suppressed_R98S_cDNA_sample_mean_densities.txt', 3),
        #('weinberg_no_misannotated', '/home/jah/projects/arlen/experiments/weinberg/RPF/results/RPF_mean_densities_no_misannotated.txt'),
        #('weinberg_mRNA', '/home/jah/projects/arlen/experiments/weinberg/mRNA/results/mRNA_mean_densities.txt'),
        #('ingolia_mRNA', '/home/jah/projects/arlen/experiments/ingolia_science/mRNA-rich-1/results/mRNA-rich-1_mean_densities.txt'),
        #('ingolia_science_no_misannotated', '/home/jah/projects/arlen/experiments/ingolia_science/Footprints-rich-1/results/Footprints-rich-1_mean_densities_no_misannotated.txt'),
        #('guydosh_3-AT', '/home/jah/projects/arlen/experiments/guydosh_cell/wild-type_3-AT/results/wild-type_3-AT_mean_densities.txt'),
        #('zinshteyn_WT_2', '/home/jah/projects/arlen/experiments/zinshteyn_plos_genetics/WT_Ribosome_Footprint_2/results/WT_Ribosome_Footprint_2_mean_densities.txt'),
        #('zinshteyn_delta_elp3', '/home/jah/projects/arlen/experiments/zinshteyn_plos_genetics/delta_elp3_Ribosome_Footprint/results/delta_elp3_Ribosome_Footprint_mean_densities.txt'),
        #('zinshteyn_delta_ncs2', '/home/jah/projects/arlen/experiments/zinshteyn_plos_genetics/delta_ncs2_Ribosome_Footprint/results/delta_ncs2_Ribosome_Footprint_mean_densities.txt'),
        #('zinshteyn_delta_ncs6_1', '/home/jah/projects/arlen/experiments/zinshteyn_plos_genetics/delta_ncs6_Ribosome_Footprint_1/results/delta_ncs6_Ribosome_Footprint_1_mean_densities.txt'),
        #('zinshteyn_delta_ncs6_2', '/home/jah/projects/arlen/experiments/zinshteyn_plos_genetics/delta_ncs6_Ribosome_Footprint_2/results/delta_ncs6_Ribosome_Footprint_2_mean_densities.txt'),
        #('belgium_8_6_13', '/home/jah/projects/arlen/experiments/belgium_8_6_13/R98S_cDNA_sample/results/R98S_cDNA_sample_read_counts.txt', 'yeast'),
        #('ingolia_cell_nothing', '/home/jah/projects/arlen/experiments/ingolia_cell/ES_cell_feeder-free__w__LIF_none_ribo_mesc_nochx_Illumina_GAII/results/ingolia_cell_codon_counts.txt'),
        #('ingolia_cell_emetine', '/home/jah/projects/arlen/experiments/ingolia_cell/ES_cell_feeder-free__w__LIF_60_s_emetine_ribo_mesc_emet_Illumina_GAII/results/emetine_read_positions.txt', 'ingolia_cell'),
        #('guo_nature', '/home/jah/projects/arlen/experiments/guo_nature/Footprint_wild-type_runs1-2/results/guo_nature_read_positions.txt', 'guo_nature'),
        #('artieri', '/home/jah/projects/arlen/experiments/artieri/Mixed_parental_ribosome_protected_fragments_replicate_1/results/artieri_RPF_replicate_1_mean_densities.txt'),
        #('artieri', '/home/jah/projects/arlen/experiments/artieri/Mixed_parental_ribosome_protected_fragments_replicate_1/results/artieri_RPF_replicate_1_mean_densities.txt'),
    ]

    data_sets = [(name, Circles.Serialize.read_file(fn, 'read_positions'), color_index)
                 for name, fn, color_index in data_sets]

    plot_averaged_codon_densities(data_sets, 'RPF_5_prime_ramps_4.png', show_end=False, past_edge=0, smooth=True)
