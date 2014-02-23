''' Utilities for counting reads mapped to each position in a gene. '''

import numpy as np
import pysam
import gtf
from collections import Counter
import matplotlib.pyplot as plt

class PositionCounts(object):
    ''' Wrapper around an array of counts of positions for an extent and for a
        buffer on either edge that allows for indexing relative to the extent's
        start and end.
    '''
    def __init__(self, extent_length, left_buffer, right_buffer, counts=None):
        self.left_buffer = left_buffer
        self.right_buffer = right_buffer
        self.extent_length = extent_length
        if counts == None:
            self.counts = np.zeros(left_buffer + right_buffer + extent_length, int)
        else:
            assert len(counts) == left_buffer + right_buffer + extent_length
            self.counts = counts

    @classmethod
    def from_string(cls, string, left_buffer, right_buffer):
        counts = np.array(map(int, string.strip().split()))
        extent_length = len(counts) - left_buffer - right_buffer
        return cls(extent_length, left_buffer, right_buffer, counts=counts)

    def __str__(self):
        string = ' '.join(map(str, self.counts)) + '\n'
        return string

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
                if adjusted_key < 0:
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

edge_buffer = 50
left_buffer = 2 * edge_buffer
right_buffer = edge_buffer

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
                #print read
                raise IndexError
            if read.qlen in relevant_lengths:
                position_counts[read.qlen][from_start] += 1

    extent_info = {'position_counts': position_counts,
                   'expression': expression,
                   'nonunique': nonunique,
                  }

    return extent_info

def dump_CDS_to_sam(bam_file_name, sam_file_name, CDS):
    bam_file = pysam.Samfile(bam_file_name)
    sam_file = pysam.Samfile(sam_file_name, 'wh', template=bam_file)
        
    region_start = CDS.start - edge_buffer
    region_end = CDS.end + edge_buffer
    overlapping_reads = bam_file.fetch(CDS.seqname, region_start, region_end)
    mapqs = Counter()
    for read in overlapping_reads:
        mapqs[read.mapq] += 1
        if read.mapq != 50:
            pass
        else:
            read_strand = '-' if read.is_reverse else '+'

            if read_strand != CDS.strand:
                pass
            else:
                sam_file.write(read)

    print mapqs

def compute_metagene_positions(gene_infos, max_gene_length):
    random_gene = gene_infos.itervalues().next()
    relevant_lengths = sorted(random_gene['position_counts'].keys())
    random_position_counts = random_gene['position_counts'][relevant_lengths[0]]
    left_buffer = random_position_counts.left_buffer
    right_buffer = random_position_counts.right_buffer
    
    from_starts = {length: PositionCounts(max_gene_length, left_buffer, right_buffer)
                   for length in relevant_lengths}
    from_ends = {length: PositionCounts(max_gene_length, left_buffer, right_buffer)
                 for length in relevant_lengths}
    
    for name in gene_infos:
        CDS_length = gene_infos[name]['CDS_length']
        position_counts = gene_infos[name]['position_counts']
        start_slice = slice(-left_buffer, CDS_length)
        end_slice = slice(-right_buffer, CDS_length)
        for length in relevant_lengths:
            from_starts[length][start_slice] += position_counts[length][start_slice]
            from_ends[length].relative_to_end[end_slice] += position_counts[length].relative_to_end[end_slice]

    return from_starts, from_ends

def plot_metagene_positions(from_starts, from_ends, figure_fn):
    relevant_lengths = sorted(from_starts.keys())
    fig, (start_ax, end_ax) = plt.subplots(2, 1, figsize=(12, 16))

    start_xs = np.arange(-21, 19)
    end_xs = np.arange(-6, 34)
    for length in relevant_lengths:
        start_ax.plot(start_xs, from_starts[length][start_xs], '.-', label=length)
        end_ax.plot(-end_xs, from_ends[length].relative_to_end[end_xs], '.-', label=length)

    start_ax.set_xlim(min(start_xs), max(start_xs))
    start_ax.set_xticks([x for x in start_xs if x % 3 == 0])
    start_ax.set_xlabel('Position of read relative to start of CDS')
    start_ax.set_ylabel('Number of uniquely mapped reads of specified length')
    
    end_ax.set_xlim(min(-end_xs), max(-end_xs))
    end_ax.set_xticks([x for x in -end_xs if x % 3 == 0])
    end_ax.set_xlabel('Position of read relative to stop codon')
    end_ax.set_ylabel('Number of uniquely mapped reads of specified length')

    leg = start_ax.legend(loc='upper right', fancybox=True)
    leg.get_frame().set_alpha(0.5)
    leg = end_ax.legend(loc='upper right', fancybox=True)
    leg.get_frame().set_alpha(0.5)
    
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
