import numpy as np
import pysam
from functools import partial
from itertools import chain
from collections import Counter
from Circles import mapping_tools
from Circles import utilities
from Circles import fastq
from Circles.annotation import Annotation_factory
import trim_cython

payload_annotation_fields = [
    ('original_name', 's'),
    ('barcode', 's'),
]
PayloadAnnotation = Annotation_factory(payload_annotation_fields)

def trim(reads,
         trimmed_reads_fn,
         min_length,
         max_read_length,
         find_start,
         find_end,
        ):
    ''' Wrapper that handles the logistics of trimming reads given functions
        find_start and find_end that take a sequence and
        returns a positions that trimming should occur at.
    '''
    trimmed_lengths = np.zeros(max_read_length + 1, int)
    too_short_lengths = np.zeros(max_read_length + 1, int)
    barcode_counts = Counter()
    
    with open(trimmed_reads_fn, 'w') as trimmed_reads_fh:
        for read in reads:
            start = find_start(read.seq)
            end = find_end(read.seq) 
            length = end - start

            if length < min_length:
                too_short_lengths[length] += 1
            else:
                trimmed_lengths[length] += 1
                barcode = read.seq[:start]
                barcode_counts[barcode] += 1
                annotation = PayloadAnnotation(original_name=read.name,
                                               barcode=barcode,
                                              )
                trimmed_record = fastq.make_record(annotation.identifier,
                                                   read.seq[start:end],
                                                   read.qual[start:end],
                                                  )
                trimmed_reads_fh.write(trimmed_record)

    return trimmed_lengths, too_short_lengths, barcode_counts

truseq_R2_rc = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
linker = 'CTGTAGGCACCATCAAT'
weinberg_linker = 'TCGTATGCCGTCTTCTGCTTG'

adapter_prefix_length = 10
max_distance = 1

find_truseq = partial(
    trim_cython.find_adapter,
    truseq_R2_rc[:adapter_prefix_length],
    max_distance,
)
trim_truseq = partial(
    trim,
    find_start=lambda seq: 0,
    find_end=find_truseq,
)

find_linker = partial(
    trim_cython.find_adapter,
    linker[:adapter_prefix_length],
    max_distance,
)
trim_linker = partial(
    trim,
    find_start=lambda seq: 0,
    find_end=find_linker,
)

find_linker_short = partial(
    trim_cython.find_short_adapter,
    linker,
)
trim_linker_short = partial(
    trim,
    find_start=lambda seq: 0,
    find_end=find_linker_short,
)

trim_polyA = partial(
    trim,
    find_start=lambda seq: 0,
    find_end=trim_cython.find_poly_A,
)

find_weinberg_linker = partial(
    trim_cython.find_short_adapter,
    weinberg_linker,
)
trim_weinberg = partial(
    trim,
    find_start=lambda seq: 8,
    find_end=find_weinberg_linker,
)

find_bartel_linker_medium = partial(
    trim_cython.find_medium_adapter,
    weinberg_linker,
)
trim_bartel_linker_medium = partial(
    trim,
    find_start=lambda seq: 0,
    find_end=find_bartel_linker_medium,
)

trim_nothing = partial(
    trim,
    find_start=lambda seq: 0,
    find_end=len,
)

def unambiguously_trimmed(filtered_bam_fn, unambiguous_bam_fn, genome_dir):
    ''' Reads that have had poly-As trimmed may have had some real RPF A's
        trimmed as well. Retains only mapped reads for which the next base in
        the reference is a non-A.
    '''
    def first_non_A(seq):
        for i, c in enumerate(seq):
            if c != 'A':
                return i
        return len(seq)

    loaded_genome = mapping_tools.load_genome('saccharomyces_cerevisiae',
                                              explicit_path=genome_dir,
                                             )
    
    bamfile = pysam.Samfile(filtered_bam_fn, 'rb')
    with pysam.Samfile(unambiguous_bam_fn, 'wb', header=bamfile.header) as unambiguous_bam_fh:
        for read in bamfile:
            seqname = bamfile.getrname(read.tid)
            if not read.is_reverse:
                start = read.pos + len(read.seq)
                end = start + 10
                after = loaded_genome[seqname][start:end]
            else:
                start = read.pos - 10
                end = read.pos
                after = loaded_genome[seqname][start:end]
                after = utilities.reverse_complement(after)

            if first_non_A(after) == 0:
                unambiguous_bam_fh.write(read)

    pysam.index(unambiguous_bam_fn)
