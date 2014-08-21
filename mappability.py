import gtf
import trim
from Sequencing import fastq, genomes
from Sequencing.annotation import Annotation_factory
from Sequencing.Parallel import map_reduce
import os
import ribosome_profiling_experiment
import positions
import pysam
from Serialize import read_positions
from collections import defaultdict

artifical_annotation_fields = [
    ('transcript_name', 's'),
    ('position', '06d'),
]

artifical_annotation = Annotation_factory(artifical_annotation_fields)

def make_artificial_reads(transcript,
                          fragment_length,
                          read_length,
                          adapter_sequence,
                          region_fetcher,
                          common_buffer,
                         ):
    transcript_sequence = transcript.retrieve_sequence(region_fetcher,
                                                       left_buffer=common_buffer,
                                                       right_buffer=common_buffer + fragment_length,
                                                      )
    # Needs to include one non-Solexa value for automatic encoding recognition.
    high_quals = fastq.encode_sanger([25] + [30]*(read_length - 1))
    for i, transcript_position in enumerate(range(-common_buffer, transcript.CDS_length + common_buffer)):
        annotation = artifical_annotation(transcript_name=transcript.name,
                                          position=transcript_position,
                                         )
        fragment_sequence = transcript_sequence[i:i + fragment_length]
        if '-' in fragment_sequence:
            # skip fragments that run off the edge of a reference sequence
            continue

        full_sequence = fragment_sequence + adapter_sequence
        read = fastq.Read(annotation.identifier, full_sequence[:read_length], high_quals)
        yield read

adapter_sequences = {'linker': trim.smRNA_linker + trim.truseq_R2_rc}

class MappabilityExperiment(ribosome_profiling_experiment.RibosomeProfilingExperiment):
    num_stages = 1

    specific_results_files = [
        ('uniqueness', read_positions, '{name}_uniqueness.hdf5'),
    ]

    specific_outputs = [
        ['uniqueness'],
    ]

    def __init__(self, **kwargs):
        super(MappabilityExperiment, self).__init__(**kwargs)

        self.adapter_sequence = adapter_sequences[self.adapter_type]
        self.fragment_length = int(kwargs['fragment_length'])
        self.common_buffer = 100

        self.work = [['record_uniqueness']]
        self.outputs = [['uniqueness']]
        self.cleanup = [[]]

    def get_reads(self):
        CDSs, _ = self.get_CDSs()
        region_fetcher = genomes.build_region_fetcher(self.file_names['genome'],
                                                      load_references=True,
                                                     )
        for transcript in CDSs:
            reads = make_artificial_reads(transcript,
                                          self.fragment_length,
                                          self.max_read_length,
                                          self.adapter_sequence,
                                          region_fetcher,
                                          self.common_buffer,
                                         )
            for read in reads:
                yield read

    def record_uniqueness(self):
        CDSs, _ = self.get_CDSs()
        uniqueness = {}
        transcripts = {}
        
        # For any genomic position that participates in a transcript, this will
        # contain a mapping to a set of all transcripts it participates in.
        genomic_to_all_transcripts = defaultdict(set)

        for transcript in CDSs:
            landmarks = {'start': 0,
                         'start_codon': transcript.transcript_start_codon,
                         'stop_codon': transcript.transcript_stop_codon,
                         'end': transcript.transcript_length,
                        }
            uniqueness[transcript.name] = {self.fragment_length: positions.PositionCounts(landmarks, self.common_buffer, self.common_buffer)}
            transcript.build_coordinate_maps(left_buffer=self.common_buffer, right_buffer=self.common_buffer)
            transcripts[transcript.name] = transcript

            for genomic_position, transcript_position in transcript.genomic_to_transcript.iteritems():
                full_position = (transcript.seqname, transcript.strand, genomic_position)
                genomic_to_all_transcripts[full_position].add((transcript.name, transcript_position))

        bam_file = pysam.Samfile(self.file_names['accepted_hits'])

        for read in bam_file:
            # If this read was incorrectly trimmed, don't record it.
            if read.qlen != self.fragment_length:
                continue

            annotation = artifical_annotation.from_prefix_identifier(read.qname)
            true_transcript = transcripts[annotation['transcript_name']]
            true_position = annotation['position']
            strand = '-' if read.is_reverse else '+'
            if strand == '+':
                five_prime = read.pos
            else:
                five_prime = read.aend - 1
            full_mapped_position = (bam_file.getrname(read.tid), strand, five_prime)

            if read.mapq < 50:
                # Flag the true source of the read as nonunique.
                uniqueness[true_transcript.name][self.fragment_length]['start_codon', true_position] = 2
                
                # Hopefully redundantly, flag the position actually mapped to as
                # nonunqiue.
                for transcript_name, transcript_position in genomic_to_all_transcripts[full_mapped_position]:
                    uniqueness[transcript_name][self.fragment_length]['start_codon', transcript_position] = 2
            else:
                # Check that any read with a MAPQ of 50 is to the expected position.
                full_true_position = (true_transcript.seqname,
                                      true_transcript.strand,
                                      true_transcript.transcript_to_genomic[true_position],
                                     )

                if read.mapq == 50 and (full_mapped_position != full_true_position):
                    raise ValueError(full_mapped_position, full_true_position)
                
                # As long as this hasn't been mapped to by some other fragment,
                # mark it as unique.
                if uniqueness[true_transcript.name][self.fragment_length]['start_codon', true_position] == 0: 
                    uniqueness[true_transcript.name][self.fragment_length]['start_codon', true_position] = 1

        self.write_file('uniqueness', uniqueness)


if __name__ == '__main__':
    script_path = os.path.realpath(__file__)
    map_reduce.controller(MappabilityExperiment, script_path)
