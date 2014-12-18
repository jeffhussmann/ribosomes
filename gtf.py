from collections import namedtuple, defaultdict, Counter
from Sequencing import genomes, utilities
import numpy as np
import Bio.Seq
import positions

gtf_fields = ['seqname',
              'source',
              'feature',
              'start',
              'end',
              'score',
              'strand',
              'frame',
              'attribute_string',
              'attribute',
             ]
Feature = namedtuple('Feature', gtf_fields)

# The attribute dictionary can't be hashed, but it is enough to hash the string
# that it was parsed out of.
Feature.__hash__ = lambda self: hash(self[:-1])

def parse_gtf_attribute(attribute):
    fields = attribute.strip(';').split('; ')
    pairs = [field.split() for field in fields]
    parsed = {name: value.strip('"') for name, value in pairs}
    return parsed

def parse_gtf_line(line, attribute_parser=parse_gtf_attribute):
    fields = line.strip().split('\t')
    fields.append(attribute_parser(fields[-1]))
    feature = Feature._make(fields)
    start = int(feature.start) - 1
    # Convert from 1-based indexing to 0-based
    end = int(feature.end) - 1
    if feature.frame != '.':
        frame = int(feature.frame)
    else:
        frame = feature.frame
    feature = feature._replace(start=start, end=end, frame=frame)
    return feature

def feature_to_line(feature):
    feature = feature._replace(start=str(feature.start + 1),
                               end=str(feature.end + 1),
                               frame=str(feature.frame),
                              )
    return '\t'.join(feature[:-1]) + '\n'

def get_all_features(gtf_fn):
    all_features = [parse_gtf_line(line) for line in open(gtf_fn)]
    return all_features

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

def get_noncoding_RNA_transcripts(gtf_fn, genome_dir):
    all_features = get_all_features(gtf_fn)
    transcripts = get_transcripts(all_features, genome_dir)
    rRNA_transcripts = []
    tRNA_transcripts = []
    other_ncRNA_transcripts = []
    for transcript in transcripts:
        if any(exon.source == 'rRNA' for exon in transcript.exons):
            rRNA_transcripts.append(transcript)
        elif any(exon.source == 'tRNA' for exon in transcript.exons):
            tRNA_transcripts.append(transcript)
        elif any('RNA' in exon.source for exon in transcript.exons):
            other_ncRNA_transcripts.append(transcript)
    return rRNA_transcripts, tRNA_transcripts, other_ncRNA_transcripts

def get_all_sources(gtf_fn):
    all_features = get_all_features(gtf_fn)
    sources = Counter(feature.source for feature in all_features)
    return sources

def get_all_exons(all_features):
    exons = [feature for feature in all_features if feature.feature == 'exon']
    return exons

def comparison_key(feature):
    return feature.seqname, feature.start, feature.end, feature.feature, feature.strand

def feature_lt(first, second):
    return comparison_key(first) < comparison_key(second)

Feature.__lt__ = feature_lt

def sort_transcripts(transcripts, by_end=False):
    if by_end:
        def key(transcript):
            return transcript.seqname, transcript.end, transcript.start, transcript.strand
    else:
        def key(transcript):
            return transcript.seqname, transcript.start, transcript.end, transcript.strand

    return sorted(transcripts, key=key)

def get_transcripts(all_features, genome_dir, utr_fn):
    region_fetcher = genomes.build_region_fetcher(genome_dir)

    feature_lists = defaultdict(list)
    for feature in all_features:
        transcript_name = feature.attribute['transcript_id']
        feature_lists[transcript_name].append(feature)

    UTR_boundaries = read_UTR_file(utr_fn)

    transcripts = [Transcript(name, features, region_fetcher, UTR_boundaries.get(name))
                   for name, features in feature_lists.iteritems()]

    return transcripts

def contained_in(exon, CDS):
    return exon.seqname == CDS.seqname and \
           CDS.start >= exon.start and \
           CDS.end <= exon.end

class Transcript(object):
    def __init__(self, name, features, region_fetcher, corrections):
        self.name = name
        self.region_fetcher = region_fetcher

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
                self.exons[0] = self.exons[0]._replace(start=five_pos)
                self.exons[-1] = self.exons[-1]._replace(end=three_pos)
            else:
                self.exons[0] = self.exons[0]._replace(start=three_pos)
                self.exons[-1] = self.exons[-1]._replace(end=five_pos)
                
        # A transcript can have more than one start_codon or stop_codon feature
        # if the codon is split across multiple exons.
        self.start_codons = [feature for feature in features if feature.feature == 'start_codon']
        self.stop_codons = [feature for feature in features if feature.feature == 'stop_codon']
        
        self.CDSs = [None for exon in self.exons]
        CDSs = [feature for feature in features if feature.feature == 'CDS']
        for CDS in CDSs:
            for e, exon in enumerate(self.exons):
                if contained_in(exon, CDS):
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

    def get_extent_sequence(self, genome_dir, left_buffer=0, right_buffer=0):
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
                                        data=sequence,
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
        if self.strand == '-':
            bases = [utilities.complement(b) for b in bases]
        sequence = ''.join(bases).upper()
        
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

    def retrieve_sequence(self, region_fetcher, left_buffer=0, right_buffer=0):
        # Remake coordinate maps to guarantee buffer sizes
        self.build_coordinate_maps(left_buffer, right_buffer)

        transcript_positions = range(self.transcript_start_codon - left_buffer,
                                     self.transcript_stop_codon + right_buffer,
                                    )
        positions = [self.transcript_to_genomic[t] for t in transcript_positions]

        bases = [region_fetcher(self.seqname, p, p + 1) for p in positions]
        if self.strand == '-':
            bases = [utilities.complement(b) for b in bases]
        sequence = ''.join(bases).upper()

        return sequence

def make_coding_sequence_fetcher(gtf_fn, genome_dir, codon_table=1):
    CDSs = get_CDSs(gtf_fn)
    gtf_dict = {t.name: t for t in CDSs}

    region_fetcher = genomes.build_region_fetcher(genome_dir)

    def coding_sequence_fetcher(name):
        ''' Returns the coding sequence for name if name is in the GTF file and
            is a well formed coding sequence, otherwise returns None.
        '''
        if name not in gtf_dict:
            raise ValueError(name)
        
        transcript = gtf_dict[name]
        if transcript.strand == '+':
            seqs = [region_fetcher(transcript.seqname, CDS.start, CDS.end + 1)
                    for CDS in transcript.CDSs if CDS]
            seqs += [region_fetcher(transcript.seqname, stop_codon.start, stop_codon.end + 1)
                     for stop_codon in transcript.stop_codons]
            seq = ''.join(seqs)
        elif transcript.strand == '-':
            seqs = [utilities.reverse_complement(region_fetcher(transcript.seqname, CDS.start, CDS.end + 1))
                    for CDS in transcript.CDSs[::-1] if CDS]
            seqs += [utilities.reverse_complement(region_fetcher(transcript.seqname, stop_codon.start, stop_codon.end + 1))
                     for stop_codon in transcript.stop_codons[::-1]]
            seq = ''.join(seqs)
        try:
            Bio.Seq.translate(seq, cds=True, table=codon_table)
        except Bio.Seq.CodonTable.TranslationError as error:
            return None

        return seq.upper()

    coding_sequence_fetcher.names = set(gtf_dict.keys())

    return coding_sequence_fetcher

def get_nonoverlapping_transcripts(transcripts, edge_buffer=0, require_same_strand=False):
    def overlaps(first, second):
        first_left = first.start - edge_buffer
        first_right = first.end + edge_buffer
        second_left = second.start - edge_buffer
        second_right = second.end + edge_buffer
        
        return first.seqname == second.seqname and \
               second_left <= first_right and \
               second_right >= first_left and \
               (not require_same_strand or first.strand == second.strand)

    overlapping = set()
    nonoverlapping = set()

    overlap_connections = [[] for t in transcripts]

    # Sort by starting position, so that if sorted_transcripts[i] does not overlap
    # sorted_transcripts[j] for i < j, then sorted_transcripts[i] can't overlap
    # sorted_transcripts[k] for k > j.
    sorted_transcripts = sort_transcripts(transcripts) 
    for i in xrange(len(sorted_transcripts)):
        first_feature = sorted_transcripts[i]
        for j in xrange(i + 1, len(sorted_transcripts)):
            second_feature = sorted_transcripts[j]
            if overlaps(first_feature, second_feature):
                overlap_connections[i].append(second_feature)
                overlap_connections[j].append(first_feature)
            else:
                break

    nonoverlapping = [g for i, g in enumerate(sorted_transcripts) if not overlap_connections[i]]

    return nonoverlapping, overlap_connections

def get_translated_transcripts(transcripts):
    ''' Returns a list of all elements in transcripts that are translated.
    '''
    translated_transcripts = [transcript for transcript in transcripts if any(transcript.CDSs)]
    return translated_transcripts

def get_CDSs(gtf_fn, genome_dir, utr_fn):
    all_features = get_all_features(gtf_fn)
    transcripts = get_transcripts(all_features, genome_dir, utr_fn)
    translated_transcripts = get_translated_transcripts(transcripts)

    return sort_transcripts(translated_transcripts)

def feature_list_to_dict(features):
    feature_dict = defaultdict(list)
    for feature in features:
        feature_dict[feature.attribute['gene_id']].append(feature)
    return feature_dict

def make_yeast_list():
    weinberg_fn = '/home/jah/projects/ribosomes/experiments/weinberg/most_weinberg_transcripts.txt'
    yeast_fn = '/home/jah/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/transcript_list.txt'
    gtf_fn = '/home/jah/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/transcriptome/genes.gtf'
    genome_dir = '/home/jah/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/genome'
    
    all_features = get_all_features(gtf_fn)
    transcripts = get_transcripts(all_features, genome_dir)
    translated = get_translated_transcripts(transcripts)
    nonoverlapping, _ = get_nonoverlapping_transcripts(transcripts, require_same_strand=True)
    weinberg_names = {line.strip() for line in open(weinberg_fn)}
    final_set = ({t.name for t in translated} & {n.name for n in nonoverlapping}) | weinberg_names

    with open(yeast_fn, 'w') as yeast_fh:
        for name in sorted(final_set):
            yeast_fh.write('{0}\n'.format(name))
