''' Utilities for counting reads mapped to each position in a gene. '''

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg', warn=False)
import numpy as np
import numbers
import pysam
from collections import Counter
import ribosomes
import gtf
import Circles.Serialize

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
    representative_position_counts = random_gene['position_counts'][relevant_lengths[0]]
    left_buffer = representative_position_counts.left_buffer
    right_buffer = representative_position_counts.right_buffer
    
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

def plot_metagene_averaged(data_sets, from_end=False):
    fig, ax = plt.subplots(figsize=(12, 12))

    plot_up_to = 40
    min_counts = 64
    
    if from_end:
        xs = np.arange(0, -plot_up_to, -1)
    else:
        xs = np.arange(plot_up_to)

    for name, codon_counts_fn in data_sets:
        codon_counts = ribosomes.read_codon_counts_file(codon_counts_fn)
        max_length = max(len(counts) for counts in codon_counts.itervalues() if counts.sum() >= min_counts)

        sum_of_normalized = np.zeros(max_length)
        long_enough_genes = np.zeros(max_length)

        for gene_name, counts in codon_counts.iteritems():
            if counts.sum() < min_counts:
                continue

            num_codons = len(counts)
            density = counts.sum() / float(num_codons)
            normalized = counts / float(density)
            
            if from_end:
                normalized = normalized[::-1]

            sum_of_normalized[:num_codons] += normalized
            long_enough_genes[:num_codons] += np.ones(num_codons)
        
        mean_densities = sum_of_normalized / long_enough_genes 

        ax.plot(xs, mean_densities[:plot_up_to], '.-', label=name)
    
    ax.legend(loc='upper right', framealpha=0.5)

    ax.set_xlabel('Position (codons)')
    ax.set_ylabel('Normalized mean reads')
    
    ax.plot(np.ones(plot_up_to), color='black', alpha=0.5)
    ax.set_ylim(0, 8)

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ax.set_aspect((xmax - xmin) / (ymax - ymin))

def get_codon_counts_generalized(gene_info, offset_type, common_buffer):
    codon_counts = PositionCounts(gene_info['CDS_length'] // 3, common_buffer, common_buffer)
    for length in set(gene_info['position_counts']) & set(ribosomes.A_site_offsets[offset_type]):
        position_counts = gene_info['position_counts'][length]
        start_index = -ribosomes.A_site_offsets[offset_type][length] - (common_buffer * 3)
        end_index = gene_info['CDS_length'] - ribosomes.A_site_offsets[offset_type][length] + (common_buffer * 3)
        codon_counts[-common_buffer:codon_counts.extent_length + common_buffer] += position_counts[start_index:end_index:3] + \
                                                              position_counts[start_index - 1:end_index - 1:3] + \
                                                              position_counts[start_index + 1:end_index + 1:3]

    return codon_counts

def plot_metagene_averaged_generalized(data_sets, from_end=False):
    fig, ax = plt.subplots(figsize=(12, 12))

    plot_up_to = 40
    min_counts = 64
    common_buffer = 10
    
    xs = np.arange(-common_buffer, plot_up_to + 1)

    for name, read_positions_fn, offset_type in data_sets:
        gene_infos = Circles.Serialize.read_file(read_positions_fn, 'read_positions')
        
        sum_of_normalized = PositionCounts(50000, common_buffer, common_buffer)
        long_enough_genes = PositionCounts(50000, common_buffer, common_buffer)

        for gene_name in gene_infos:
            codon_counts = get_codon_counts_generalized(gene_infos[gene_name], offset_type, common_buffer)

            if codon_counts.counts.sum() < min_counts:
                continue

            num_codons = len(codon_counts.counts)
            density = codon_counts.counts.sum() / float(num_codons)
            normalized = codon_counts / float(density)
            
            start_slice = slice(-common_buffer, codon_counts.extent_length + common_buffer)
            sum_of_normalized[start_slice] += normalized[start_slice]
            long_enough_genes[start_slice] += np.ones(num_codons)
        
        mean_densities = sum_of_normalized / long_enough_genes 

        ax.plot(xs, mean_densities[xs], '.-', label=name)
    
    ax.legend(loc='upper right', framealpha=0.5)

    ax.set_xlabel('Position (codons)')
    ax.set_ylabel('Normalized mean reads')

    ax.set_xlim(min(xs), max(xs))
    
    ax.plot(xs, [1 for x in xs], color='black', alpha=0.5)
    ax.set_ylim(0, 8)

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ax.set_aspect((xmax - xmin) / (ymax - ymin))

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
        ('belgium_3_5_14', '/home/jah/projects/arlen/experiments/belgium_3_5_14/wt/results/wt_read_positions.txt', 'yeast'),
        #('belgium_8_6_13', '/home/jah/projects/arlen/experiments/belgium_8_6_13/R98S_cDNA_sample/results/R98S_cDNA_sample_read_counts.txt', 'yeast'),
        #('ingolia_cell_nothing', '/home/jah/projects/arlen/experiments/ingolia_cell/ES_cell_feeder-free__w__LIF_none_ribo_mesc_nochx_Illumina_GAII/results/ingolia_cell_codon_counts.txt'),
        #('ingolia_cell_emetine', '/home/jah/projects/arlen/experiments/ingolia_cell/ES_cell_feeder-free__w__LIF_60_s_emetine_ribo_mesc_emet_Illumina_GAII/results/emetine_read_positions.txt', 'ingolia_cell'),
        #('guo_nature', '/home/jah/projects/arlen/experiments/guo_nature/Footprint_wild-type_runs1-2/results/guo_nature_read_positions.txt', 'guo_nature'),
    ]
    plot_metagene_averaged_generalized(data_sets, from_end=False)
