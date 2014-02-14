from collections import namedtuple, defaultdict, Counter
from itertools import islice

gtf_fields = ['seqname',
              'source',
              'feature',
              'start',
              'end',
              'score',
              'strand',
              'frame',
              'attribute',
             ]
Gene = namedtuple('Gene', gtf_fields)

# Work-around (temporary?) to reconcile automatic parsing of attribute string
# into a dictionary with the need to hash Genes to check for set membership
Gene.__hash__ = lambda self: hash(self[:-1])

def parse_attribute(attribute):
    fields = attribute.strip(';').split('; ')
    pairs = [field.split() for field in fields]
    parsed = {name: value.strip('"') for name, value in pairs}
    return parsed

def parse_gtf_line(line):
    fields = line.strip().split('\t')
    fields[-1] = parse_attribute(fields[-1])
    gene = Gene._make(fields)
    start = int(gene.start) - 1
    # Convert from 1-based indexing to 0-based
    end = int(gene.end) - 1
    if gene.frame != '.':
        frame = int(gene.frame)
    else:
        frame = gene.frame
    gene = gene._replace(start=start, end=end, frame=frame)
    return gene

def get_all_features(gtf_fn):
    all_features = [parse_gtf_line(line) for line in open(gtf_fn)]
    return all_features

def get_noncoding_RNA_genes(gtf_fn):
    all_features = get_all_features(gtf_fn)
    rRNA_genes = []
    tRNA_genes = []
    other_ncRNA_genes = []
    for feature in all_features:
        if feature.source == 'rRNA':
            rRNA_genes.append(feature)
        elif feature.source == 'tRNA':
            tRNA_genes.append(feature)
        elif 'RNA' in feature.source:
            other_ncRNA_genes.append(feature)
    return rRNA_genes, tRNA_genes, other_ncRNA_genes

def get_all_sources(gtf_fn):
    all_features = get_all_features(gtf_fn)
    sources = Counter(feature.source for feature in all_features)
    return sources

def get_all_CDSs(all_features):
    # This used to also check if gene.source == 'protein_coding', but as I
    # generalize to other organisms, I am encountering GTF files that don't fill
    # in meaningful source values.
    CDSs = [feature for feature in all_features if feature.feature == 'CDS']

    return CDSs

def get_all_exons(all_features):
    exons = [feature for feature in all_features if feature.feature == 'exon']
    return exons

def sort_features(features):
    def key(feature):
        return feature.seqname, feature.start, feature.feature

    return sorted(features, key=key)

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
        self.exons = [feature for feature in features if feature.feature == 'exon']
        self.CDSs = [None for exon in self.exons]
        
        CDSs = [feature for feature in features if feature.feature == 'CDS']
        for CDS in CDSs:
            for e, exon in enumerate(self.exons):
                if contained_in(exon, CDS):
                    self.CDSs[e] = CDS
                    break
            else:
                raise ValueError(CDS, self.exons)

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
    ''' Returns all transcripts that have a single exon and a CDS.
    '''
    single_exon_transcripts = [transcript for transcript in transcripts
                               if len(transcript.exons) == 1 and transcript.CDSs[0]]
    return single_exon_transcripts

def get_simple_CDSs(gtf_fn, exclude_from=None):
    ''' Returns all single exon CDSs that do not overlap any other CDS.
        Any gene IDs listed in any file name in exclude_from are not eligible
        but also do not disqualify a CDS by overlapping it - think dubious ORFs
        or uncharacterized genes.
    '''
    all_features = get_all_features(gtf_fn)
    #CDSs = get_all_CDSs(all_features)
    exons = get_all_exons(all_features)
    transcripts = get_transcripts(all_features)

    relevant_exons = exons
    if exclude_from:
        for name_list_fn in exclude_from:
            relevant_exons = exclude_dubious_ORFs(relevant_exons, name_list_fn)

    nonoverlapping, _ = get_nonoverlapping(relevant_exons, edge_buffer=20)
    nonoverlapping = set(nonoverlapping)
    single_exon_transcripts = get_single_exon_transcripts(transcripts)
    simple_CDSs = [transcript.CDSs[0] for transcript in single_exon_transcripts
                   if transcript.exons[0] in nonoverlapping]
    simple_CDSs = sort_features(list(simple_CDSs))

    return simple_CDSs

def feature_list_to_dict(features):
    feature_dict = defaultdict(list)
    for feature in features:
        feature_dict[feature.attribute['gene_id']].append(feature)
    return feature_dict
