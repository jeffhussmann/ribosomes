import numpy as np
import pysam
from functools import partial
from itertools import chain
from collections import Counter
from Sequencing import genomes, utilities, fastq
from Circles.annotation import Annotation_factory
from trim_cython import *

payload_annotation_fields = [
    ('original_name', 's'),
    ('barcode', 's'),
    ('trimmed', 's'),
]
PayloadAnnotation = Annotation_factory(payload_annotation_fields)

trimmed_twice_annotation_fields = payload_annotation_fields + [('retrimmed', 's')]
TrimmedTwiceAnnotation = Annotation_factory(trimmed_twice_annotation_fields)

def trim(reads, trimmed_fn, min_length, max_read_length, find_start, find_end, second_time=False):
    ''' Wrapper that handles the logistics of trimming reads given functions
        find_start and find_end that take a sequence and
        returns a positions that trimming should occur at.
    '''
    trimmed_lengths = np.zeros(max_read_length + 1, int)
    too_short_lengths = np.zeros(max_read_length + 1, int)
    barcode_counts = Counter()
    
    with open(trimmed_fn, 'w') as trimmed_fh:
        for read in reads:
            start = find_start(read.seq)
            end = find_end(read.seq) 
            length = end - start

            if length < min_length:
                too_short_lengths[length] += 1
            else:
                trimmed_lengths[length] += 1
                barcode = read.seq[:start]
                trimmed = read.seq[end:]
                barcode_counts[barcode] += 1
                if second_time:
                    payload_annotation = PayloadAnnotation.from_identifier(read.name)
                    annotation = TrimmedTwiceAnnotation(retrimmed=trimmed,
                                                        **payload_annotation)
                else:
                    annotation = PayloadAnnotation(original_name=read.name,
                                                   barcode=barcode,
                                                   trimmed=trimmed,
                                                  )
                trimmed_record = fastq.make_record(annotation.identifier,
                                                   read.seq[start:end],
                                                   read.qual[start:end],
                                                  )
                trimmed_fh.write(trimmed_record)

    return trimmed_lengths, too_short_lengths, barcode_counts

truseq_R2_rc = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
smRNA_linker = 'CTGTAGGCACCATCAAT'
bartel_linker = 'TCGTATGCCGTCTTCTGCTTG'

adapter_prefix_length = 10
max_distance = 1

finders = {'truseq':        (lambda seq: 0,
                             partial(find_adapter, truseq_R2_rc[:adapter_prefix_length], max_distance),
                            ),
           'linker':        (lambda seq: 0,
                             partial(find_adapter, smRNA_linker[:adapter_prefix_length], max_distance),
                            ),
           'linker_short':  (lambda seq: 0,
                             partial(find_short_adapter, smRNA_linker),
                            ),
           'polyA':         (lambda seq: 0,
                             find_poly_A,
                            ),
           'weinberg':      (lambda seq: 8,
                             partial(find_short_adapter, bartel_linker),
                            ),
           'bartel_medium': (lambda seq: 0,
                             partial(find_medium_adapter, bartel_linker),
                            ),
           'nothing':       (lambda seq: 0,
                             len,
                            ),
           'jeff':          (find_jeff_start,
                             partial(find_short_adapter, smRNA_linker),
                            ),
          }

max_barcode_length = {'weinberg': 8,
                      'jeff': 18,
                     }

bound_trim = {key: partial(trim, find_start=find_start, find_end=find_end)
              for key, (find_start, find_end) in finders.items()}

def unambiguously_trimmed(bam_fn, unambiguous_bam_fn, genome_dir):
    ''' Reads that have had poly-As trimmed may have had some real RPF A's
        trimmed as well. Retains only mapped reads for which the last aligned
        base and the following base in the reference are both non-A.
    '''
    genome = genomes.load_entire_genome(genome_dir)
    
    bamfile = pysam.Samfile(bam_fn)
    with pysam.Samfile(unambiguous_bam_fn, 'wb', header=bamfile.header) as unambiguous_bam_fh:
        for read in bamfile:
            rname = bamfile.getrname(read.tid)

            if not read.is_reverse:
                if read.positions[-1] == bamfile.lengths[read.tid] - 1:
                    # There is no next base to get
                    continue
                last_position = read.positions[-1]
                last_base, next_base = genome[rname][last_position:last_position + 2]
            else:
                if read.positions[0] == 0:
                    # There is no next base to get
                    continue
                last_position = read.positions[0]
                last_base, next_base = utilities.reverse_complement(genome[rname][last_position - 1:last_position + 1])

            if last_base.upper() != 'A' and next_base.upper() != 'A':
                unambiguous_bam_fh.write(read)

    pysam.index(unambiguous_bam_fn)
