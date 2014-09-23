import matplotlib
matplotlib.use('Agg', warn=False)
import matplotlib.pyplot as plt
import glob
import os
import pysam
import numpy as np
from collections import Counter

from Sequencing import fastq
from Sequencing import annotation
from Sequencing import utilities
from Sequencing import sam
from Sequencing import Visualize
from Sequencing import mapping_tools
from Sequencing import genomes
from Sequencing.annotation import Annotation_factory
from Sequencing.Parallel import map_reduce
from Sequencing.Parallel import split_file

import trim
import ribosomes
import contaminants
import gtf
from ribosome_profiling_experiment import RibosomeProfilingExperiment

mapping_annotation_fields = [
    ('RNAME', '-<{RNAME_length}s'),
    ('POS', '0{POS_length}d'),
    ('cigar', '-<0{max_cigar_length}.{max_cigar_length}s'),
    ('strand', 's'),
    ('barcode', '-<{max_barcode_length}s'),
    ('original_name', 's'),
]

position_annotation_fields = mapping_annotation_fields[:4]
fragment_annotation_fields = mapping_annotation_fields[:5]

collapsed_annotation_fields = fragment_annotation_fields + [('count', '04d')]

class BarcodedRibosomeProfilingExperiment(RibosomeProfilingExperiment):
    num_stages = 3

    def __init__(self, **kwargs):
        RibosomeProfilingExperiment.__init__(self, **kwargs)

        self.max_barcode_length = trim.max_barcode_length[self.adapter_type]
        self.max_cigar_length = 10

        specific_results_files = [
            ('annotated_clean_sam', 'sam', '{name}_annotated_clean.sam'),
            ('sorted_clean_sam', 'sam_sorted', '{name}_sorted_clean.sam'),

            ('barcode_counts', 'counts', '{name}_barcode_counts.txt'),
            ('collapsed_sam', 'sam', '{name}_collapsed.sam'),
            ('collapsed_bam', 'bam', '{name}_collapsed.bam'),
            ('amplification_counts', 'counts', '{name}_amplification_counts.txt'),
            ('collisions', 'concatenate', '{name}_collisions.txt'),
        ]

        specific_figure_files = []

        specific_outputs = [
            ['sorted_clean_sam'],
            ['amplification_counts',
             'collisions',
             'collapsed_bam',
            ],
            [],
        ]

        specific_work = [
            [(self.annotate_mappings, 'Annotating mappings'),
             (self.sort_mappings, 'Sorting mappings'),
            ],
            [(self.collapse_fragments, 'Collapsing fragments'),
            ],
            [#(self.get_read_positions, 'Counting mapping positions'),
             #(self.get_error_profile, 'Getting error profile'),
             #(self.get_recycling_ratios, 'Getting recycling ratios'),
            ],
        ]

        specific_cleanup = [
            [],
            [self.index_collapsed_bam,
            ],
            [],
        ]

        self.results_files.extend(specific_results_files)
        self.figure_files.extend(specific_figure_files)
        
        # Splice in the collapsing stage in the middle
        self.outputs = [self.outputs[0] + specific_outputs[0],
                        specific_outputs[1],
                        self.outputs[1] + specific_outputs[2],
                       ]
        
        self.work = [self.work[0] + specific_work[0],
                     specific_work[1],
                     self.work[1] + specific_work[2],
                    ]
        
        self.cleanup = [self.cleanup[0] + specific_cleanup[0],
                        specific_cleanup[1],
                        self.cleanup[1] + specific_cleanup[2],
                       ]

        self.make_file_names()

        self.genome_index = genomes.get_genome_index(self.file_names['genome'])
        spec_kwargs = {'RNAME_length': genomes.max_RNAME_length(self.genome_index),
                       'POS_length': genomes.max_POS_length(self.genome_index),
                       'max_barcode_length': self.max_barcode_length,
                       'max_cigar_length': self.max_cigar_length,
                      }
        self.MappingAnnotation = Annotation_factory(mapping_annotation_fields, **spec_kwargs)
        self.PositionAnnotation = Annotation_factory(position_annotation_fields, **spec_kwargs)
        self.FragmentAnnotation = Annotation_factory(fragment_annotation_fields, **spec_kwargs)
        self.CollapsedAnnotation = Annotation_factory(collapsed_annotation_fields, **spec_kwargs)
        
    def annotate_mappings(self):
        clean_bam_file = pysam.Samfile(self.file_names['clean_bam'])
        annotated_sam_file = pysam.Samfile(self.file_names['annotated_clean_sam'],
                                           'wh',
                                           template=clean_bam_file,
                                          )
        for aligned_read in clean_bam_file:
            if aligned_read.mapq == 50:
                payload_annotation = trim.PayloadAnnotation.from_identifier(aligned_read.qname)
                rname = clean_bam_file.getrname(aligned_read.tid)
                strand = '-' if aligned_read.is_reverse else '+'
                mapping_annotation = self.MappingAnnotation(RNAME=rname,
                                                            POS=aligned_read.pos,
                                                            cigar=aligned_read.cigarstring,
                                                            strand=strand,
                                                            **payload_annotation)
                aligned_read.qname = mapping_annotation.identifier
                annotated_sam_file.write(aligned_read)
    
    def sort_mappings(self):
        sam.sort(self.file_names['annotated_clean_sam'],
                 self.file_names['sorted_clean_sam'],
                )
                
    def get_sorted_sam_lines(self):
        ''' Get this piece's set of SAM lines from the merged SAM file, ensuring
            that the piece boundary doesn't break up a set of mappings to the
            same position.
        '''
        get_position = annotation.make_convertor(self.MappingAnnotation,
                                                 self.PositionAnnotation,
                                                )

        sam_lines = split_file.piece(self.merged_file_names['sorted_clean_sam'],
                                     self.num_pieces,
                                     self.which_piece,
                                     'sam',
                                     key=get_position,
                                    )
        return sam_lines

    def collapse_fragments(self):
        get_position = annotation.make_convertor(self.MappingAnnotation,
                                                 self.PositionAnnotation,
                                                )
        get_fragment = annotation.make_convertor(self.MappingAnnotation,
                                                 self.FragmentAnnotation,
                                                )
        amplification_counts = Counter()

        sq_lines = sam.get_sq_lines(self.merged_file_names['sorted_clean_sam'])
        sam_lines = self.get_sorted_sam_lines()

        with open(self.file_names['collapsed_sam'], 'w') as collapsed_fh, \
             open(self.file_names['collisions'], 'w') as collisions_fh:
            
            for sq_line in sq_lines:
                collapsed_fh.write(sq_line)

            position_groups = utilities.group_by(sam_lines, get_position)
            for position_annotation, position_lines in position_groups:
                fragment_counts = Counter()
                position_count = len(position_lines)
                fragment_groups = utilities.group_by(position_lines, get_fragment)
                for fragment_annotation, fragment_lines in fragment_groups:
                    fragment_count = len(fragment_lines)
                    fragment_counts[fragment_count] += 1
                    amplification_counts['{},{}'.format(position_count, fragment_count)] += 1
                    collapsed_annotation = self.CollapsedAnnotation(count=fragment_count, **fragment_annotation)
                    new_line = sam.splice_in_name(fragment_lines[0], collapsed_annotation.identifier)
                    collapsed_fh.write(new_line)
                fragment_counts = utilities.counts_to_array(fragment_counts)
                if position_count > 100:
                    collisions_fh.write(position_annotation.identifier + '\n')
                    collisions_fh.write(','.join(map(str, fragment_counts)) + '\n')

        sam.make_sorted_bam(self.file_names['collapsed_sam'],
                            self.file_names['collapsed_bam'],
                           )

        self.write_file('amplification_counts', amplification_counts)
                
    def index_collapsed_bam(self):
        sam.index_bam(self.merged_file_names['collapsed_bam'])

if __name__ == '__main__':
    script_path = os.path.realpath(__file__)
    map_reduce.controller(BarcodedRibosomeProfilingExperiment, script_path)
