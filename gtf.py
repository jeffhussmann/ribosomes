from collections import defaultdict, Counter
from Sequencing import genomes, utilities
import positions
import transcript as transcript_utils

class Feature(object):
    def __init__(self, line=None):
        if line == None:
            # Allow __init__ to be called with no arguments to allow the
            # @classmethod constructor below.
            return

        fields = line.strip().split('\t')
        
        self.seqname = fields[0]
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
        
        self.parse_attribute_string()
    
    @classmethod
    def from_fields(cls, seqname, source, feature, start, end, score, strand, frame, attribute_string):
        obj = cls()
        obj.seqname = seqname
        obj.source = source
        obj.feature = feature
        obj.start = start
        obj.end = end
        obj.score = score
        obj.strand = strand
        obj.frame = frame
        obj.attribute_string = attribute_string
        obj.parse_attribute_string()
        return obj

    @classmethod
    def sequence_edge(cls, seqname, position):
        obj = cls()
        obj.seqname = seqname
        obj.feature = 'edge'
        obj.start = position
        obj.end = position
        obj.strand = '.'
        obj.source = '.'
        obj.score = '.'
        obj.frame = '.'
        obj.attribute_string = '.'
        obj.parse_attribute_string()
        return obj
    
    def parse_attribute_string(self):
        if self.attribute_string == '.':
            parsed = {}
        else:
            fields = self.attribute_string.strip(';').split('; ')
            pairs = [field.split() for field in fields]
            parsed = {name: value.strip('"') for name, value in pairs}

        self.attribute = parsed

    def __str__(self):
        fields = (self.seqname,
                  self.source,
                  self.feature,
                  str(self.start + 1),
                  str(self.end + 1),
                  self.score,
                  self.strand,
                  str(self.frame),
                  self.attribute_string,
                 )
        line = '\t'.join(fields)
        return line
    
    @property
    def pasteable(self):
        return '{0}:{1}-{2}'.format(self.seqname, self.start, self.end)
    
    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return str(self) == str(other)

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

def get_noncoding_RNA_transcripts(gtf_fn):
    all_features = get_all_features(gtf_fn)
    transcripts = transcript_utils.get_transcripts(all_features, '/dev/null', '/dev/null')
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
