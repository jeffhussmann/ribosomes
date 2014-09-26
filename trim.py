import numpy as np
import string
import pysam
from functools import partial
from itertools import chain
from collections import Counter
from Sequencing import genomes, utilities, fastq, sam
from Sequencing.annotation import Annotation_factory
from trim_cython import *
from Sequencing.adapters_cython import *

payload_annotation_fields = [
    ('original_name', 's'),
    ('barcode', 's'),
    ('trimmed_seq', 's'),
    ('trimmed_qual', 's'),
]
PayloadAnnotation = Annotation_factory(payload_annotation_fields)

retrimmed_fields = [('retrimmed_seq', 's'),
                    ('retrimmed_qual', 's'),
                   ]
trimmed_twice_annotation_fields = payload_annotation_fields + retrimmed_fields
TrimmedTwiceAnnotation = Annotation_factory(trimmed_twice_annotation_fields)

_sanitize_table = string.maketrans('/', chr(ord('/') - 1))
def sanitize_qual(qual):
    ''' If '/' is in a qname, bowtie/tophat truncate the rest of the qname.
        If qual strings of trimmed portions of reads are to be put in qnames,
        '/' needs to be downgraded to chr(ord('/') - 1) = '.'
    '''
    sanitized = qual.translate(_sanitize_table)
    return sanitized

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
                trimmed_seq = read.seq[end:]
                trimmed_qual = sanitize_qual(read.qual[end:])
                barcode_counts[barcode] += 1
                if second_time:
                    payload_annotation = PayloadAnnotation.from_identifier(read.name)
                    annotation = TrimmedTwiceAnnotation(retrimmed_seq=trimmed_seq,
                                                        retrimmed_qual=trimmed_qual,
                                                        **payload_annotation)
                else:
                    annotation = PayloadAnnotation(original_name=read.name,
                                                   barcode=barcode,
                                                   trimmed_seq=trimmed_seq,
                                                   trimmed_qual=trimmed_qual,
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
smRNA_15_linker = 'TGGAATTCTCGGGTGCCAAGG'

adapter_prefix_length = 10
max_distance = 1

finders = {'truseq':        (lambda seq: 0,
                             partial(find_adapter, truseq_R2_rc[:adapter_prefix_length], max_distance),
                            ),
           'linker':        (lambda seq: 0,
                             partial(find_adapter, smRNA_linker[:adapter_prefix_length], max_distance),
                            ),
           'linker_15':     (lambda seq: 0,
                             partial(find_adapter, smRNA_15_linker[:adapter_prefix_length], max_distance),
                            ),
           'polyA':         (lambda seq: 0,
                             find_poly_A,
                            ),
           'weinberg':      (lambda seq: 8,
                             partial(find_adapter, bartel_linker[:adapter_prefix_length], max_distance),
                            ),
           'nothing':       (lambda seq: 0,
                             len,
                            ),
           'jeff':          (find_jeff_start,
                             partial(find_adapter, smRNA_linker[:adapter_prefix_length], max_distance),
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

def trim_polyA_from_unmapped(unmapped_bam_file_name,
                             trimmed_fastq_file_name,
                             min_length,
                             max_read_length,
                             second_time=False,
                            ):
    reads = (fastq.Read(read.qname, read.seq, read.qual) for read in pysam.Samfile(unmapped_bam_file_name))
    trim(reads,
         trimmed_fastq_file_name,
         min_length,
         max_read_length,
         lambda seq: 0,
         find_poly_A,
         second_time=second_time,
        )

def extend_polyA_ends(bam_fn, extended_bam_fn, genome_dir, trimmed_twice=False):
    bam_file = pysam.Samfile(bam_fn)
    region_fetcher = genomes.build_region_fetcher(genome_dir,
                                                  load_references=True,
                                                  sam_file=bam_file,
                                                 )
    with pysam.Samfile(extended_bam_fn, 'wb', template=bam_file) as extended_bam_file:
        for mapping in bam_file:
            if mapping.is_unmapped:
                extended_bam_file.write(mapping)
                continue

            if trimmed_twice:
                # Trailing poly-As were removed by the second trimming step. 
                annotation = TrimmedTwiceAnnotation.from_identifier(mapping.qname)
                polyA_seq = annotation['retrimmed_seq']
                polyA_qual = annotation['retrimmed_qual']
                new_qname = PayloadAnnotation.from_prefix_identifier(mapping.qname).identifier
            else:
                annotation = PayloadAnnotation.from_identifier(mapping.qname)
                polyA_seq = annotation['trimmed_seq']
                polyA_qual = annotation['polyA_qual']
                new_qname = '{0}_{1}'.format(annotation['original_name'],
                                             annotation['barcode'],
                                            )

            num_trimmed = len(polyA_seq)
            
            if mapping.is_reverse:
                after = region_fetcher(mapping.tid, mapping.pos - num_trimmed, mapping.pos)
                after = utilities.reverse_complement(after)
            else:
                after = region_fetcher(mapping.tid, mapping.aend, mapping.aend + num_trimmed)

            extra_genomic_As = 0
            for b in after:
                if b == 'A':
                    extra_genomic_As += 1
                else:
                    break

            nongenomic_length = num_trimmed - extra_genomic_As

            if mapping.is_reverse:
                nongenomic_start = mapping.pos - 1 - extra_genomic_As
            else:
                # Note: 'aend points to one past the last aligned residue'
                nongenomic_start = mapping.aend + extra_genomic_As

            extra_genomic_seq = polyA_seq[:extra_genomic_As]
            soft_clipped_seq = polyA_seq[extra_genomic_As:]
            extra_genomic_qual = polyA_qual[:extra_genomic_As]
            soft_clipped_qual = polyA_qual[extra_genomic_As:]

            extra_seq = extra_genomic_seq + soft_clipped_seq
            extra_qual = extra_genomic_qual + soft_clipped_qual

            if mapping.is_reverse:
                final_cigar_block_index = 0
                extended_seq = utilities.reverse_complement(extra_seq) + mapping.seq
                extended_qual = extra_qual[::-1] + mapping.qual
                mapping.pos = mapping.pos - extra_genomic_As
            else:
                final_cigar_block_index = -1
                extended_seq = mapping.seq + extra_seq
                extended_qual = mapping.qual + extra_qual

            # Note: writing to mapping.seq destroys mapping.qual, so
            # mapping.qual needs to be retrieved above
            mapping.seq = extended_seq
            mapping.qual = extended_qual

            op, length = mapping.cigar[final_cigar_block_index]
            if op != 0:
                raise ValueError
            length += extra_genomic_As
            
            updated_cigar = mapping.cigar
            updated_cigar[final_cigar_block_index] = (op, length)
            if len(soft_clipped_seq) > 0:
                soft_clipped_block = [(4, len(soft_clipped_seq))]
                if final_cigar_block_index == 0:
                    updated_cigar = soft_clipped_block + updated_cigar
                elif final_cigar_block_index == -1:
                    updated_cigar = updated_cigar + soft_clipped_block

            mapping.cigar = updated_cigar

            if mapping.tags:
                # Clear the MD tag since the possible addition of bases to the
                # alignment may have made it inaccurate. 
                filtered_tags = filter(lambda t: t[0] != 'MD', mapping.tags)
                mapping.tags = filtered_tags

            set_nongenomic_length(mapping, nongenomic_length)

            mapping.qname = new_qname
            extended_bam_file.write(mapping)

def get_nongenomic_length(mapping):
    tags = {name: value for name, value in mapping.tags}
    nongenomic_length = tags['ZN']
    return nongenomic_length

def set_nongenomic_length(mapping, nongenomic_length):
    mapping.setTag('ZN', nongenomic_length)

def trim_mismatches_from_start(bam_fn, trimmed_bam_fn, genome_dir, relevant_lengths, max_read_length):
    ''' Remove all consecutive Q30+ mismatches from the beginning of alignments,
        under the assumption that these represent untemplated additions during
        reverse transcription.
        Characterize the mismatches.
    '''
    type_shape = (len(relevant_lengths),
                  max_read_length,
                  fastq.MAX_EXPECTED_QUAL + 1,
                  6,
                  6,
                 )
    type_counts = np.zeros(type_shape, int)
    length_to_index = {length: i for i, length in enumerate(relevant_lengths)}
    
    bam_file = pysam.Samfile(bam_fn)
    region_fetcher = genomes.build_region_fetcher(genome_dir,
                                                  load_references=True,
                                                  sam_file=bam_file,
                                                 )

    with pysam.Samfile(trimmed_bam_fn, 'wb', template=bam_file) as trimmed_bam_file:
        for mapping in bam_file:
            if sam.contains_indel_pysam(mapping):
                continue
            if mapping.is_unmapped:
                trimmed_bam_file.write(mapping)
                continue

            if mapping.is_reverse:
                aligned_pairs = mapping.aligned_pairs[::-1]
                index_lookup = utilities.base_to_complement_index
            else:
                aligned_pairs = mapping.aligned_pairs
                index_lookup = utilities.base_to_index

            decoded_qual = fastq.decode_sanger(mapping.qual)
            length_index = length_to_index.get(mapping.qlen, None)
            
            bases_to_trim = 0
            found_trim_point = False
            for read_index, ref_index in aligned_pairs:
                if read_index == None:
                    # This shouldn't be able to be triggered since alignments
                    # containing indels are ruled out above.
                    continue

                if mapping.is_reverse:
                    corrected_read_index = mapping.qlen - 1 - read_index
                else:
                    corrected_read_index = read_index

                ref_base = region_fetcher(mapping.tid, ref_index, ref_index + 1)
                read_base = mapping.seq[read_index]
                read_qual = decoded_qual[read_index]
                if length_index != None:
                    coords = (length_index,
                              corrected_read_index,
                              read_qual,
                              index_lookup[ref_base],
                              index_lookup[read_base],
                             )
                    type_counts[coords] += 1

                if not found_trim_point:
                    if read_base != ref_base and read_qual >= 30:
                        bases_to_trim += 1
                    else:
                        first_ref_index = ref_index
                        found_trim_point = True

            if bases_to_trim == 0:
                trimmed_mapping = mapping
            else:
                trimmed_mapping = pysam.AlignedRead()
                trimmed_mapping.qname = mapping.qname
                trimmed_mapping.tid = mapping.tid
                
                # first_ref_index has been set above to the be index of the
                # reference base aligned to the first non-trimmed base in the
                # read. If the mapping is forward, this will be the new pos.
                # If the mapping is reverse, the pos won't change.
                if mapping.is_reverse:
                    first_ref_index = mapping.pos
                trimmed_mapping.pos = first_ref_index

                trimmed_mapping.is_reverse = mapping.is_reverse
                trimmed_mapping.mapq = mapping.mapq

                if mapping.is_reverse:
                    # bases_to_trim is never zero here, so there is no danger
                    # of minus zero
                    trimmed_slice = slice(None, -bases_to_trim)
                else:
                    trimmed_slice = slice(bases_to_trim, None)

                trimmed_mapping.seq = mapping.seq[trimmed_slice]
                trimmed_mapping.qual = mapping.qual[trimmed_slice]
                trimmed_mapping.rnext = -1
                trimmed_mapping.pnext = -1

                trimmed_length = len(mapping.seq) - bases_to_trim
                if mapping.is_reverse:
                    # truncate_cigar_blocks removes bases from the end, which is
                    # what we want for reverse reads
                    trimmed_cigar = sam.truncate_cigar_blocks(mapping.cigar, trimmed_length)
                else:
                    # We want to remove from the beginning, so flip the blocks,
                    # remove from the end, then flip back.
                    flipped_cigar = mapping.cigar[::-1]
                    trimmed_cigar = sam.truncate_cigar_blocks(flipped_cigar, trimmed_length)
                    trimmed_cigar = trimmed_cigar[::-1]
                
                trimmed_mapping.cigar = trimmed_cigar

            trimmed_bam_file.write(trimmed_mapping)

    return type_counts
