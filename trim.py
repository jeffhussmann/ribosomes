import numpy as np
import string
import pysam
from functools import partial
from itertools import chain
from collections import Counter
import Sequencing.fastq as fastq
import Sequencing.genomes as genomes
import Sequencing.sam as sam
import Sequencing.utilities as utilities
import Sequencing.sw as sw
from Sequencing.annotation import Annotation_factory
from trim_cython import *
from Sequencing.adapters_cython import *
from Sequencing import sw

payload_annotation_fields = [
    ('original_name', 's'),
    ('left_seq', 's'),
    ('left_qual', 's'),
    ('right_seq', 's'),
    ('right_qual', 's'),
]
PayloadAnnotation = Annotation_factory(payload_annotation_fields)

retrimmed_fields = [('retrimmed_left_seq', 's'),
                    ('retrimmed_left_qual', 's'),
                    ('retrimmed_right_seq', 's'),
                    ('retrimmed_right_qual', 's'),
                   ]
trimmed_twice_annotation_fields = payload_annotation_fields + retrimmed_fields
TrimmedTwiceAnnotation = Annotation_factory(trimmed_twice_annotation_fields)

def trim(reads, find_start=None, find_end=None, second_time=False):
    ''' Wrapper that handles the logistics of trimming reads given functions
        find_start and find_end that take a sequence and returns positions
        that trimming should occur at.
    '''
    if find_start == None:
        find_start = lambda seq: 0
    if find_end == None:
        find_end = len

    for read in reads:
        start = find_start(read.seq)
        end = find_end(read.seq) 

        left_seq = read.seq[:start]
        left_qual = fastq.sanitize_qual(read.qual[:start])
        right_seq = read.seq[end:]
        right_qual = fastq.sanitize_qual(read.qual[end:])
        if second_time:
            payload_annotation = PayloadAnnotation.from_identifier(read.name)
            annotation = TrimmedTwiceAnnotation(retrimmed_left_seq=left_seq,
                                                retrimmed_left_qual=left_qual,
                                                retrimmed_right_seq=right_seq,
                                                retrimmed_right_qual=right_qual,
                                                **payload_annotation)
        else:
            annotation = PayloadAnnotation(original_name=read.name,
                                           left_seq=left_seq,
                                           left_qual=left_qual,
                                           right_seq=right_seq,
                                           right_qual=right_qual,
                                          )
        trimmed_read = fastq.Read(annotation.identifier,
                                  read.seq[start:end],
                                  read.qual[start:end],
                                 )
        yield trimmed_read

truseq_R2_rc = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
smRNA_linker = 'CTGTAGGCACCATCAAT'
bartel_linker = 'TCGTATGCCGTCTTCTGCTTG'
smRNA_15_linker = 'TGGAATTCTCGGGTGCCAAGG'

full_linker = smRNA_linker + truseq_R2_rc

adapter_prefix_length = 15
max_distance = 1

def trim_by_local_alignment(adapter, seq):
    ''' Try to find a near-exact match. If this fails, do a local alignment. '''
    trim_at = find_adapter(adapter[:adapter_prefix_length], 1, seq)

    if trim_at > len(seq) - adapter_prefix_length:
        alignment, = sw.generate_alignments(adapter,
                                            seq,
                                            'unpaired_adapter',
                                            max_alignments=1,
                                           )

        score_diff = 2 * len(alignment['path']) - alignment['score']
        adapter_start_in_seq = sw.first_target_index(alignment['path'])

        if alignment['path'] and score_diff <= 10. / 22 * len(alignment['path']):
            trim_at = adapter_start_in_seq

    return trim_at

finders = {'truseq':        (None,
                             partial(find_adapter, truseq_R2_rc[:adapter_prefix_length], max_distance),
                            ),
           'truseq_local':  (None,
                             partial(trim_by_local_alignment, truseq_R2_rc),
                            ),
           'linker':        (None,
                             partial(find_adapter, smRNA_linker[:adapter_prefix_length], max_distance),
                            ),
           'linker_local':  (None,
                             partial(trim_by_local_alignment, full_linker),
                            ),
           'linker_15':     (None,
                             partial(find_adapter, smRNA_15_linker[:adapter_prefix_length], max_distance),
                            ),
           'polyA':         (None,
                             find_poly_A,
                            ),
           'weinberg':      (lambda seq: 8,
                             partial(find_adapter, bartel_linker[:adapter_prefix_length], max_distance),
                            ),
           'nothing':       (None,
                             None,
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

def trim_polyA_from_unmapped(unmapped_reads,
                             trimmed_fastq_file_name,
                             second_time=False,
                            ):
    trimmed_reads = trim(unmapped_reads,
                         find_start=lambda x: 0,
                         find_end=find_poly_A,
                         second_time=second_time,
                        )
    with open(trimmed_fastq_file_name, 'w') as trimmed_fh:
        for trimmed_read in trimmed_reads:
            trimmed_fh.write(str(trimmed_read))

def untrim_reads(trimmed_reads, second_time=False):
    if second_time:
        Annotation = TrimmedTwiceAnnotation
        left_seq_key = 'retrimmed_left_seq'
        left_qual_key = 'retrimmed_left_qual'
        right_seq_key = 'retrimmed_right_seq'
        right_qual_key = 'retrimmed_right_qual'
    else:
        Annotation = PayloadAnnotation
        left_seq_key = 'left_seq'
        left_qual_key = 'left_qual'
        right_seq_key = 'right_seq'
        right_qual_key = 'right_qual'

    for trimmed_read in trimmed_reads:
        annotation = Annotation.from_identifier(trimmed_read.name)
        name = trimmed_read.name
        seq = annotation[left_seq_key] + trimmed_read.seq + annotation[right_seq_key]
        qual = annotation[left_qual_key] + trimmed_read.qual + annotation[right_qual_key]
        read = fastq.Read(name, seq, qual)
        yield read

def extend_polyA_ends(bam_fn, extended_bam_fn, genome_dir, trimmed_twice=False):
    bam_file = pysam.Samfile(bam_fn)
    region_fetcher = genomes.build_region_fetcher(genome_dir,
                                                  load_references=True,
                                                  sam_file=bam_file,
                                                 )

    # Adding bases to the end of minus strand mappings produces a file
    # that is not necessarily sorted, so re-sort. 
    alignment_sorter = sam.AlignmentSorter(bam_file.references,
                                           bam_file.lengths,
                                           extended_bam_fn,
                                          )

    with alignment_sorter:
        for mapping in bam_file:
            extended_mapping = extend_polyA_end(mapping, region_fetcher, trimmed_twice)
            alignment_sorter.write(extended_mapping)

def extend_polyA_end(mapping, region_fetcher, trimmed_twice=False):
    if mapping.is_unmapped:
        return mapping

    if trimmed_twice:
        # Trailing poly-As were removed by the second trimming step. 
        annotation = TrimmedTwiceAnnotation.from_identifier(mapping.qname)
        polyA_seq = annotation['retrimmed_right_seq']
        polyA_qual = annotation['retrimmed_right_qual']
        new_qname = PayloadAnnotation.from_prefix_identifier(mapping.qname).identifier
    else:
        annotation = PayloadAnnotation.from_identifier(mapping.qname)
        polyA_seq = annotation['right_seq']
        polyA_qual = annotation['right_qual']
        new_qname = '{0}_{1}_{2}'.format(annotation['original_name'],
                                         annotation['left_seq'],
                                         annotation['left_qual'],
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
        soft_clipped_block = [(sam.BAM_CSOFT_CLIP, len(soft_clipped_seq))]
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

    return mapping

def get_nongenomic_length(mapping):
    tags = {name: value for name, value in mapping.tags}
    # If ZN wasn't set, we want to return 0.
    nongenomic_length = tags.get('ZN', 0)
    return nongenomic_length

def set_nongenomic_length(mapping, nongenomic_length):
    # setTag throws a fit if it is is given a long, so coerce to int
    mapping.setTag('ZN', int(nongenomic_length))

def trim_mismatches_from_start(mapping, region_fetcher, type_counts):
    ''' Remove all consecutive Q30+ mismatches from the beginning of alignments,
        under the assumption that these represent untemplated additions during
        reverse transcription.
        Characterize the mismatches into type_counts.
    '''
    if sam.contains_indel_pysam(mapping) or mapping.is_unmapped:
        set_nongenomic_length(mapping, 0)
        return mapping

    if mapping.is_reverse:
        aligned_pairs = mapping.aligned_pairs[::-1]
        index_lookup = utilities.base_to_complement_index
    else:
        aligned_pairs = mapping.aligned_pairs
        index_lookup = utilities.base_to_index

    decoded_qual = fastq.decode_sanger(mapping.qual)
    
    bases_to_trim = 0
    found_trim_point = False
    first_ref_index = None
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
        coords = (mapping.qlen,
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

    if first_ref_index == None:
        raise ValueError('first_ref_index not set')

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
        trimmed_mapping.is_secondary = mapping.is_secondary
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
            # Remove blocks from the end
            trimmed_cigar = sam.truncate_cigar_blocks_up_to(mapping.cigar, trimmed_length)
        else:
            # Remove blocks from the beginning
            trimmed_cigar = sam.truncate_cigar_blocks_from_beginning(mapping.cigar, trimmed_length)
        
        trimmed_mapping.cigar = trimmed_cigar

    return trimmed_mapping

def trim_nongenomic_polyA_from_end(mapping, region_fetcher):
    ''' If a mapping ends in a polyA stretch, soft clip from the first
    nongenomic A onward.
    '''
    if sam.contains_indel_pysam(mapping) or mapping.is_unmapped:
        return mapping

    first_ref_index = None
    if mapping.is_reverse:
        bases_to_trim = 0
        poly_T_end = find_poly_T(mapping.seq)
        for read_index, ref_index in mapping.aligned_pairs[::-1]:
            if read_index == None:
                # indels are filtered out above, so this can only be 
                # a skip from splicing
                continue
            if read_index > poly_T_end:
                first_ref_index = ref_index
                continue
            ref_base = region_fetcher(mapping.tid, ref_index, ref_index + 1)
            if ref_base != 'T':
                bases_to_trim = read_index + 1
                break
            else:
                # first_ref_index needs to be set to the last position
                # that passed that 'are you genomic?' test
                first_ref_index = ref_index
    else:
        first_ref_index = mapping.pos
        bases_to_trim = 0
        poly_A_start = find_poly_A(mapping.seq)
        for read_index, ref_index in mapping.aligned_pairs:
            if read_index < poly_A_start:
                continue
            ref_base = region_fetcher(mapping.tid, ref_index, ref_index + 1)
            if ref_base != 'A':
                bases_to_trim = len(mapping.seq) - read_index
                break
    
    if first_ref_index == None:
        print mapping
        raise ValueError('first_ref_index not set')

    if bases_to_trim > 0:
        mapping.pos = first_ref_index

        trimmed_length = len(mapping.seq) - bases_to_trim
        soft_clipped_block = [(sam.BAM_CSOFT_CLIP, bases_to_trim)]
        if mapping.is_reverse:
            # Remove blocks from the beginning.
            trimmed_cigar = sam.truncate_cigar_blocks_from_beginning(mapping.cigar, trimmed_length)
            updated_cigar = soft_clipped_block + trimmed_cigar
        else:
            # Remove blocks from the end.
            trimmed_cigar = sam.truncate_cigar_blocks_up_to(mapping.cigar, trimmed_length)
            updated_cigar = trimmed_cigar + soft_clipped_block
        
        mapping.cigar = updated_cigar
    
    if mapping.tags:
        # Clear the MD tag since the possible removal of bases to the
        # alignment may have made it inaccurate. 
        # TODO: now have machinery to make it accurate.
        filtered_tags = filter(lambda t: t[0] != 'MD', mapping.tags)
        mapping.tags = filtered_tags

    set_nongenomic_length(mapping, bases_to_trim)

    return mapping

if __name__ == '__main__':
    fastq_fn = '/home/jah/projects/ribosomes/experiments/guydosh_cell/dom34KO_CHX/data/SRR1042854.fastq'
    seqs = [r.seq for _, r in zip(xrange(100000), fastq.reads(fastq_fn))]
    seqs = utilities.progress_bar(len(seqs), seqs)
    adapter = full_linker
    count = 0
    counts = Counter()
    for seq in seqs:
        counts[trim_by_local_alignment(adapter, seq)] += 1
