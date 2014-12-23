import gtf
import numpy as np
import Sequencing.genomes as genomes
import Sequencing.utilities as utilities
from collections import defaultdict
import positions
import interval_tree
import Bio.Seq

class Transcript(object):
    def __init__(self,
                 name,
                 features,
                 overlap_finder,
                 region_fetcher,
                 corrections=None,
                 codon_table=1,
                ):
        self.name = name
        self.region_fetcher = region_fetcher
        self.codon_table = codon_table
        self.overlap_finder = overlap_finder.overlapping

        # Further processing assumes that features is sorted by (start, end).
        features = sorted(features)

        strands = {feature.strand for feature in features}
        if len(strands) > 1:
            raise ValueError(self.name)
        self.strand = strands.pop()

        seqnames = {feature.seqname for feature in features}
        if len(seqnames) > 1:
            raise ValueError(self.name)
        self.seqname = seqnames.pop()

        self.exons = [feature for feature in features if feature.feature == 'exon']

        # Apply corrections
        if corrections:
            seqname, strand, five_pos, three_pos = corrections
            if self.seqname != seqname or self.strand != strand:
                raise ValueError

            if self.strand == '+':
                self.exons[0].start = five_pos
                self.exons[-1].end = three_pos
            else:
                self.exons[0].start = three_pos
                self.exons[-1].end = five_pos
                
        # A transcript can have more than one start_codon or stop_codon feature
        # if the codon is split across multiple exons.
        self.start_codons = [feature for feature in features if feature.feature == 'start_codon']
        self.stop_codons = [feature for feature in features if feature.feature == 'stop_codon']
        
        self.CDSs = [None for exon in self.exons]
        CDSs = [feature for feature in features if feature.feature == 'CDS']
        for CDS in CDSs:
            for e, exon in enumerate(self.exons):
                if CDS.is_contained_in(exon):
                    self.CDSs[e] = CDS
                    break
            else:
                raise ValueError(CDS, self.exons)

        if self.start_codons:
            if self.strand == '+':
                self.first_start_codon_position = min(sc.start for sc in self.start_codons)
            elif self.strand == '-':
                self.first_start_codon_position = max(sc.end for sc in self.start_codons)
        else:
            self.first_start_codon_position = None
        
        if self.stop_codons:
            if self.strand == '+':
                self.first_stop_codon_position = min(sc.start for sc in self.stop_codons)
            elif self.strand == '-':
                self.first_stop_codon_position = max(sc.end for sc in self.stop_codons)
        else:
            self.first_stop_codon_position = None

        self.start = min(exon.start for exon in self.exons)
        self.end = max(exon.end for exon in self.exons)

    @property
    def comparison_key(self):
        return self.seqname, self.start, self.end, self.strand

    def __lt__(self, other):
        return self.comparison_key < other.comparison_key

    def build_coordinate_maps(self, left_buffer=0, right_buffer=0):
        ''' Make dictionaries mapping from genomic coordinates to transcript
            coordinates and vice-versa.
        '''
        if self.strand == '+':
            exon_position_lists = [np.arange(exon.start, exon.end + 1) for exon in self.exons]
        elif self.strand == '-':
            exon_position_lists = [np.arange(exon.end, exon.start - 1, -1) for exon in self.exons[::-1]]
        
        exon_positions = np.concatenate(exon_position_lists)

        self.transcript_length = len(exon_positions)
        
        # Add some upstream and downstream bases.
        upstream_transcript = np.arange(-left_buffer, 0)
        downstream_transcript = np.arange(self.transcript_length, self.transcript_length + right_buffer)
        if self.strand == '+':
            upstream_positions = np.arange(self.start - left_buffer, self.start)
            downstream_positions = np.arange(self.end + 1, self.end + 1 + right_buffer)
        elif self.strand == '-':
            upstream_positions = np.arange(self.end + left_buffer, self.end, -1)
            downstream_positions = np.arange(self.start - 1, self.start - 1 - right_buffer, -1)

        self.transcript_to_genomic = dict(enumerate(exon_positions))
        self.transcript_to_genomic.update(zip(upstream_transcript, upstream_positions)) 
        self.transcript_to_genomic.update(zip(downstream_transcript, downstream_positions)) 

        self.genomic_to_transcript = {g: t for t, g in self.transcript_to_genomic.iteritems()}

        if self.first_stop_codon_position != None:
            if self.first_start_codon_position != None:
                genomic_start_codon = self.first_start_codon_position
            else:
                # E. coli genes that aren't initiated with AUG don't have a start
                # codon listed in the gtf file.
                if self.strand == '+':
                    genomic_start_codon = self.start
                elif self.strand == '-':
                    genomic_start_codon = self.end

            self.transcript_start_codon = self.genomic_to_transcript[genomic_start_codon]
            self.transcript_stop_codon = self.genomic_to_transcript[self.first_stop_codon_position]
            # By convention, CDS_length includes no bases of the stop codon.
            self.CDS_length = self.transcript_stop_codon - self.transcript_start_codon
    
    def build_extent_maps(self, left_buffer=0, right_buffer=0):
        ''' Make dictionaries mapping from genomic coordinates to transcript
            coordinates and vice-versa.
        '''
        start = self.first_start_codon_position
        stop = self.first_stop_codon_position
        if self.strand == '+':
            extent_positions = np.arange(start, stop)
        elif self.strand == '-':
            extent_positions = np.arange(start, stop, -1)
        
        self.extent_length = abs(stop - start)
        
        # Add some upstream and downstream bases.
        upstream_extent = np.arange(-left_buffer, 0)
        downstream_extent = np.arange(self.extent_length, self.extent_length + right_buffer)
        if self.strand == '+':
            upstream_positions = np.arange(start - left_buffer, start)
            downstream_positions = np.arange(stop, stop + right_buffer)
        elif self.strand == '-':
            upstream_positions = np.arange(start + left_buffer, start, -1)
            downstream_positions = np.arange(stop, stop - right_buffer, -1)

        self.extent_to_genomic = dict(enumerate(extent_positions))
        self.extent_to_genomic.update(zip(upstream_extent, upstream_positions)) 
        self.extent_to_genomic.update(zip(downstream_extent, downstream_positions)) 
        
        self.genomic_to_extent = {g: e for e, g in self.extent_to_genomic.iteritems()}

    def get_extent_sequence(self, left_buffer=0, right_buffer=0):
        ''' Get the sequence of the extent. Useful for looking at gene with
        annotated frameshifts.
        '''
        sequence = self.region_fetcher(self.seqname,
                                       min(self.genomic_to_extent),
                                       max(self.genomic_to_extent) + 1,
                                      )
        if self.strand == '-':
            sequence = utilities.reverse_complement(sequence)

        extent_landmarks = {'start': 0,
                            'end': self.extent_length,
                           }
        return positions.PositionCounts(extent_landmarks,
                                        left_buffer,
                                        right_buffer,
                                        data=np.asarray(sequence, 'c'),
                                       )

    def get_transcript_sequence(self, left_buffer=0, right_buffer=0):
        ''' Get the sequence of the mature transcript.
        '''
        # Remake coordinate maps to guarantee buffer sizes
        self.build_coordinate_maps(left_buffer, right_buffer)

        transcript_positions = range(-left_buffer,
                                     self.transcript_length + right_buffer,
                                    )
        genomic_positions = [self.transcript_to_genomic[t] for t in transcript_positions]

        bases = [self.region_fetcher(self.seqname, p, p + 1) for p in genomic_positions]
        sequence = ''.join(bases).upper()
        if self.strand == '-':
            sequence = utilities.complement(sequence)
        
        landmarks = {'start': 0,
                     'start_codon': self.transcript_start_codon,
                     'stop_codon': self.transcript_stop_codon,
                     'end': self.transcript_length,
                    }

        return positions.PositionCounts(landmarks,
                                        left_buffer,
                                        right_buffer,
                                        data=np.asarray(sequence, 'c'),
                                       )

    def get_coding_sequence(self):
        transcript_sequence = self.get_transcript_sequence()
        coding_sequence = ''.join(transcript_sequence['start_codon', 0:self.CDS_length + 3])
        
        # Ensure that the coding sequence is well-formed.
        try:
            Bio.Seq.translate(coding_sequence, cds=True, table=self.codon_table)  
        except Bio.Seq.CodonTable.TranslationError as error:
            coding_sequence = None

        return coding_sequence

    def is_spliced_out(self, position):
        ''' Returns True if the genomic position is between the start and end of
            this transcript but not part of it.
        '''
        is_within = self.start < position < self.end
        not_part_of = position not in self.genomic_to_transcript
        return is_within and not_part_of

    def delete_coordinate_maps(self):
        del self.transcript_to_genomic
        del self.genomic_to_transcript

    def __str__(self):
        return '{0} {1}:{2}-{3}'.format(self.name, self.seqname, self.start, self.end)

def write_UTR_file(UTR_boundaries, UTR_fn):
    def sort_key(name):
        seqname, strand, five_pos, three_pos = UTR_boundaries[name]
        return (seqname, min(five_pos, three_pos), max(five_pos, three_pos), strand)

    with open(UTR_fn, 'w') as UTR_fh:
        for name in sorted(UTR_boundaries, key=sort_key):
            seqname, strand, five_pos, three_pos = UTR_boundaries[name]
            line = '{0}\t{1}\t{2}\t{3}\t{4}\n'.format(name,
                                                      seqname,
                                                      strand,
                                                      five_pos,
                                                      three_pos,
                                                     )
            UTR_fh.write(line)

def read_UTR_file(UTR_fn):
    UTR_boundaries = {}
    for line in open(UTR_fn):
        name, seqname, strand, five_pos, three_pos = line.strip().split()
        UTR_boundaries[name] = (seqname, strand, int(five_pos), int(three_pos))

    return UTR_boundaries

def get_transcripts(all_features, genome_dir, utr_fn):
    import gff
    region_fetcher = genomes.build_region_fetcher(genome_dir, load_references=True)

    gff_features = gff.get_all_features('/home/jah/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/saccharomyces_cerevisiae.gff')
    overlap_finder = interval_tree.NamedOverlapFinder(gff_features)

    feature_lists = defaultdict(list)
    for feature in all_features:
        transcript_name = feature.attribute['transcript_id']
        feature_lists[transcript_name].append(feature)

    UTR_boundaries = read_UTR_file(utr_fn)

    transcripts = [Transcript(name, features, overlap_finder, region_fetcher, UTR_boundaries.get(name))
                   for name, features in feature_lists.iteritems()]

    return transcripts

if __name__ == '__main__':
    import gff
    import pprint

    gtf_fn = '/home/jah/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/transcriptome/genes.gtf'
    gff_fn = '/home/jah/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/saccharomyces_cerevisiae.gff'
    utr_fn = '/home/jah/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/transcriptome/inferred_UTR_lengths.txt'
    genome_dir = '/home/jah/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/genome'

    CDSs = gtf.get_CDSs(gtf_fn, genome_dir, utr_fn)
    gff_features = gff.get_all_features(gff_fn)
    overlap_finder = interval_tree.NamedOverlapFinder(gff_features).overlapping
    name_to_object = {f.attribute['Name']: f for f in gff_features if 'Name' in f.attribute}

    def is_interesting(possible, CDS):
        if possible.feature == 'chromosome':
            return False

        same_gene = name_to_object[CDS.name]
        if same_gene.is_ancestor_of(possible):
            return False

        if possible.parent:
            return False

        return True

    for CDS in CDSs:
        overlapping = overlap_finder(CDS.seqname, CDS.start, CDS.end)
        same_gene = name_to_object[CDS.name]
        overlapping = [f for f in overlapping if is_interesting(f, CDS)]
        print CDS
        print
        for f in overlapping:
            print f
            pprint.pprint(f.attribute)

        print
        raw_input()
