from collections import defaultdict, Counter
from Sequencing import genomes, utilities
import positions
import transcript

class Feature(object):
    def __init__(self, line):
        fields = line.strip().split('\t')
        
        self.seqname = fields[0].strip('chr')
        self.source = fields[1]
        self.feature = fields[2]
        self.start = int(fields[3]) - 1
        self.end = int(fields[4]) - 1
        self.score = fields[5]
        self.strand = fields[6]
        self.frame = fields[7]
        if self.frame != '.':
            self.frame = int(self.frame)
        self.attribute_string = fields[8]
        
        self.attribute = self.parse_attribute_string()
    
    def parse_attribute_string(self):
        fields = self.attribute_string.strip(';').split('; ')
        pairs = [field.split() for field in fields]
        parsed = {name: value.strip('"') for name, value in pairs}
        return parsed

    def __str__(self):
        fields = (self.seqname,
                  self.source,
                  self.feature,
                  str(self.start + 1),
                  str(self.end + 1),
                  self.score,
                  self.strand,
                  str(self.frame),
                  #self.attribute_string,
                 )
        line = '\t'.join(fields)
        return line
    
    def __hash__(self):
        return hash(str(self))

    @property
    def comparison_key(self):
        key = (self.seqname,
               self.start,
               self.end,
               self.feature,
               self.strand,
              )
        return key

    def __lt__(self, other):
        return self.comparison_key < other.comparison_key

    def is_contained_in(self, other):
        return self.seqname == other.seqname and \
               self.start >= other.start and \
               self.end <= other.end

def get_all_features(gtf_fn):
    all_features = [Feature(line) for line in open(gtf_fn)]
    return all_features

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
    transcripts = transcript.get_transcripts(all_features, genome_dir, utr_fn)
    translated_transcripts = get_translated_transcripts(transcripts)

    return sorted(translated_transcripts)

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
