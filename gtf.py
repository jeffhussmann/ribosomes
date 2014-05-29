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

def get_all_features(gtf_fn):
    all_features = [parse_gtf_line(line) for line in open(gtf_fn)]
    return all_features

def get_noncoding_RNA_transcripts(gtf_fn):
    all_features = get_all_features(gtf_fn)
    transcripts = get_transcripts(all_features)
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

def sort_features(features, by_end=False):
    if by_end:
        def key(feature):
            return feature.seqname, feature.end, feature.start, feature.feature
    else:
        def key(feature):
            return feature.seqname, feature.start, feature.end, feature.feature

    return sorted(features, key=key)

def sort_transcripts(transcripts, by_end=False):
    if by_end:
        def key(transcript):
            return transcript.seqname, transcript.end, transcript.start, transcript.strand
    else:
        def key(transcript):
            return transcript.seqname, transcript.start, transcript.end, transcript.strand

    return sorted(transcripts, key=key)

def get_transcripts(all_features):
    feature_lists = defaultdict(list)
    for feature in all_features:
        transcript_name = feature.attribute['transcript_id']
        feature_lists[transcript_name].append(feature)

    transcripts = [Transcript(name, features) for name, features in feature_lists.iteritems()]

    return transcripts

def contained_in(exon, CDS):
    return exon.seqname == CDS.seqname and \
           CDS.start >= exon.start and \
           CDS.end <= exon.end

class Transcript(object):
    def __init__(self, name, features):
        self.name = name

        # Further processing assumes that features is sorted by (start, end).
        features = sort_features(features)

        self.exons = [feature for feature in features if feature.feature == 'exon']

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

        strands = {feature.strand for feature in features}
        if len(strands) > 1:
            raise ValueError(self.name)
        self.strand = strands.pop()

        seqnames = {feature.seqname for feature in features}
        if len(seqnames) > 1:
            raise ValueError(self.name)
        self.seqname = seqnames.pop()

        self.start = min(exon.start for exon in self.exons)
        self.end = max(exon.end for exon in self.exons)

    def build_coordinate_maps(self):
        ''' Make dictionaries mapping from genomic coordinates to transcript
            coordinates and vice-versa.
        '''
        if self.strand == '+':
            exon_position_lists = [np.arange(exon.start, exon.end + 1) for exon in self.exons]
        elif self.strand == '-':
            exon_position_lists = [np.arange(exon.end, exon.start - 1, -1) for exon in self.exons[::-1]]
        
        exon_positions = np.concatenate(exon_position_lists)

        self.transcript_length = len(exon_positions)
        # Add some upstream and downstream bases to support yeast's lack of
        # annotated UTRs.
        # This could cause problems in other organisms if one isoform is a
        # truncation of another.
        buffer_length = 50
        upstream_transcript = np.arange(-buffer_length, 0)
        downstream_transcript = np.arange(self.transcript_length, self.transcript_length + buffer_length)
        if self.strand == '+':
            upstream_positions = np.arange(self.start - buffer_length, self.start)
            downstream_positions = np.arange(self.end + 1, self.end + 1 + buffer_length)
        elif self.strand == '-':
            upstream_positions = np.arange(self.end + buffer_length, self.end, -1)
            downstream_positions = np.arange(self.start - 1, self.start - 1 - buffer_length, -1)

        self.transcript_to_genomic = dict(enumerate(exon_positions))
        self.transcript_to_genomic.update(zip(upstream_transcript, upstream_positions)) 
        self.transcript_to_genomic.update(zip(downstream_transcript, downstream_positions)) 

        self.genomic_to_transcript = {g: t for t, g in self.transcript_to_genomic.iteritems()}

        if self.start_codons and self.stop_codons:
            if self.strand == '+':
                genomic_start_codon = min(sc.start for sc in self.start_codons)
                genomic_stop_codon = min(sc.start for sc in self.stop_codons)
            elif self.strand == '-':
                genomic_start_codon = max(sc.end for sc in self.start_codons)
                genomic_stop_codon = max(sc.end for sc in self.stop_codons)
            self.transcript_start_codon = self.genomic_to_transcript[genomic_start_codon]
            self.transcript_stop_codon = self.genomic_to_transcript[genomic_stop_codon]
            # By convention, CDS_lengths includes no bases of the stop codon.
            self.CDS_length = self.transcript_stop_codon - self.transcript_start_codon

    def delete_coordinate_maps(self):
        del self.transcript_to_genomic
        del self.genomic_to_transcript

    def __str__(self):
        return '{0} {1}:{2}-{3}'.format(self.name, self.seqname, self.start, self.end)

def make_coding_sequence_fetcher(gtf_fn, genome_dir):
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
            Bio.Seq.translate(seq, cds=True)
        except Bio.Seq.CodonTable.TranslationError as error:
            print name
            return None

        return seq.upper()

    return coding_sequence_fetcher

def get_extent_by_name(gtf_fn, name):
    all_features = get_all_features(gtf_fn)
    entries = [feature for feature in all_features if feature.attribute['gene_id'] == name]
    start_codon = [entry for entry in entries if entry.feature == 'start_codon'][0]
    stop_codon = [entry for entry in entries if entry.feature == 'stop_codon'][0]
    if any(entry.strand == '-' for entry in entries):
        start = stop_codon.end
        end = start_codon.start + 2
    else:
        start = start_codon.start
        # Haven't decided what the convention should be for which base is the end
        end = stop_codon.end
    seqname = start_codon.seqname
    strand = start_codon.strand
    return seqname, strand, start, end

def get_nonoverlapping(features, edge_buffer=0):
    ''' Returns all elements of features that do not overlap any other element of
        features.
    '''
    def overlaps(first, second):
        first_left = first.start - edge_buffer
        first_right = first.end + edge_buffer
        second_left = second.start - edge_buffer
        second_right = second.end + edge_buffer
        return first.seqname == second.seqname and \
               second_left <= first_right and \
               second_right >= first_left

    overlapping = set()
    nonoverlapping = set()

    overlap_connections = [[] for f in features]

    # Sort by starting position, so that if sorted_features[i] does not overlap
    # sorted_features[j] for i < j, then sorted_features[i] can't overlap
    # sorted_features[k] for k > j.
    sorted_features = sort_features(features) 
    for i in xrange(len(sorted_features)):
        first_feature = sorted_features[i]
        for j in xrange(i + 1, len(sorted_features)):
            second_feature = sorted_features[j]
            if overlaps(first_feature, second_feature):
                overlap_connections[i].append(second_feature)
                overlap_connections[j].append(first_feature)
            else:
                break

    nonoverlapping = [g for i, g in enumerate(sorted_features) if not overlap_connections[i]]

    return nonoverlapping, overlap_connections

def exclude_dubious_ORFs(CDSs, dubious_ORFs_fn):
    dubious_ORF_names = {line.strip() for line in open(dubious_ORFs_fn)}
    non_dubious_ORFs = [CDS for CDS in CDSs if CDS.attribute['gene_id'] not in dubious_ORF_names]
    return non_dubious_ORFs

def get_single_exon_transcripts(transcripts):
    ''' Returns a list of all elements in transcripts that have a single exon.
    '''
    single_exon_transcripts = [transcript for transcript in transcripts if len(transcript.exons) == 1]
    return single_exon_transcripts

def get_translated_transcripts(transcripts):
    ''' Returns a list of all elements in transcripts that are translated.
    '''
    translated_transcripts = [transcript for transcript in transcripts if any(transcript.CDSs)]
    return translated_transcripts

def get_simple_CDSs(gtf_fn, exclude_from=None):
    ''' Returns all single exon CDSs that do not overlap any other CDS.
        Any gene IDs listed in any file name in exclude_from are not eligible
        but also do not disqualify a CDS by overlapping it - think dubious ORFs
        or uncharacterized genes.
    '''
    all_features = get_all_features(gtf_fn)
    exons = get_all_exons(all_features)

    relevant_exons = exons
    if exclude_from:
        for name_list_fn in exclude_from:
            relevant_exons = exclude_dubious_ORFs(relevant_exons, name_list_fn)

    nonoverlapping, _ = get_nonoverlapping(relevant_exons, edge_buffer=20)
    nonoverlapping = set(nonoverlapping)

    transcripts = get_transcripts(all_features)
    translated_transcripts = get_translated_transcripts(transcripts)
    single_exon_transcripts = get_single_exon_transcripts(translated_transcripts)
    simple_CDSs = [transcript.CDSs[0] for transcript in single_exon_transcripts
                   if transcript.exons[0] in nonoverlapping]
    simple_CDSs = sort_features(list(simple_CDSs))

    return simple_CDSs

def get_CDSs(gtf_fn):
    all_features = get_all_features(gtf_fn)
    transcripts = get_transcripts(all_features)
    translated_transcripts = get_translated_transcripts(transcripts)

    return sort_transcripts(translated_transcripts)

def feature_list_to_dict(features):
    feature_dict = defaultdict(list)
    for feature in features:
        feature_dict[feature.attribute['gene_id']].append(feature)
    return feature_dict
