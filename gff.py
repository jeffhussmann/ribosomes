import gtf
import urllib
import pprint
import gff_transcript
from collections import Counter

class Feature(gtf.Feature):
    def __init__(self, line=None):
        super(Feature, self).__init__(line)
        self.parent = None
        self.children = set()

    def parse_attribute_string(self):
        fields = self.attribute_string.split(';')
        pairs = [field.split('=') for field in fields]
        parsed = {name: urllib.unquote(value).strip('"') for name, value in pairs}
        return parsed

    def populate_connections(self, id_to_object):
        parent_id = self.attribute.get('Parent')
        if parent_id:
            self.parent = id_to_object[parent_id]
            self.parent.children.add(self)

    def get_descendants(self):
        descendants = set(self.children)
        for child in self.children:
            descendants.update(child.get_descendants())
        return descendants

    def print_family(self, level=0):
        print '{0}{1}'.format('\t'*min(level, 1), self)
        if level == 0:
            pprint.pprint(self.attribute)
        for child in self.children:
            child.print_family(level=level + 1)

    def is_ancestor_of(self, possible_descendant):
        if self == possible_descendant:
            return True
        else:
            return any(child.is_ancestor_of(possible_descendant) for child in self.children)

new_comments = '''\
# seqnames replaced.
# Gene's with a five_prime_UTR_intron have had their mRNA extended upstream 50 bps of the intron.
'''

def convert_SGD_file(old_fn, new_fn):
    ''' Replaces seqnames in a file from SGD to agree with those from ENSEMBL - 
    removes 'chr' from beginning and replaces mt with Mito.
    '''
    def standardize_seqname(seqname):
        seqname = seqname[3:]
        if seqname == 'mt':
            seqname = 'Mito'
        return seqname

    def attribute_to_string(attribute):
        entries = []
        for key, value in sorted(attribute.items()):
            key = urllib.quote(key, safe='')
            value = urllib.quote(value, safe='')
            entry = '{0}={1}'.format(key, value)
            entries.append(entry)

        return ';'.join(entries)

    def process_features_exons(features):
        new_features = []

        num_previous_exons = Counter()
        for feature in sorted(features):
            if feature.feature == 'CDS':
                ancestor = feature
                while ancestor.parent:
                    ancestor = ancestor.parent

                gene_id = ancestor.attribute['ID']
                parent_id = feature.parent.attribute['ID']

                exon_id = '{0}_exon_{1}'.format(gene_id, num_previous_exons[gene_id] + 1)
                num_previous_exons[gene_id] += 1

                exon_attribute = {'ID': exon_id,
                                  'Parent': parent_id,
                                 }

                exon = Feature.from_fields(feature.seqname,
                                           feature.source,
                                           'exon',
                                           feature.start,
                                           feature.end,
                                           feature.score,
                                           feature.strand,
                                           '.',
                                           attribute_to_string(exon_attribute),
                                          )
                
                feature.attribute['Parent'] = exon_id
                feature.attribute_string = attribute_to_string(feature.attribute)
                
                new_features.append(exon)

        with_exons = sorted(features + new_features)

        populate_connections(with_exons)

        upstream_exons = []
            
        for feature in with_exons:
            if feature.feature == 'five_prime_UTR_intron':
                feature.feature = 'intron'
                
                ancestor = feature
                while ancestor.parent:
                    ancestor = ancestor.parent

                gene_id = ancestor.attribute['ID']
                parent_id = feature.parent.attribute['ID']

                exon_id = '{0}_exon_{1}'.format(gene_id, num_previous_exons[gene_id] + 1)
                num_previous_exons[gene_id] += 1

                # Find the exon of this mRNA that is closest to the intron
                # and extend it to the edge of the intron if necessary.
                # Temporarily, add a fake exon upstream of the intron.

                mRNA = feature.parent
                if mRNA.feature != 'mRNA':
                    raise ValueError
                
                exons = [c for c in mRNA.children if c.feature == 'exon']

                if feature.strand == '+':
                    distance = lambda e: e.start
                elif feature.strand == '-':
                    distance = lambda e: e.end

                closest = min(exons, key=distance)
                
                exon_attribute = {'ID': exon_id,
                                  'Parent': parent_id,
                                 }

                if feature.strand == '+':
                    closest.start = feature.end + 1

                    exon_start = feature.start - 50
                    exon_end = feature.start - 1

                    ancestor.start = exon_start
                    mRNA.start = exon_start

                elif feature.strand == '-':
                    closest.end = feature.start - 1
                    
                    exon_start = feature.end + 1
                    exon_end = feature.end + 50

                    ancestor.end = exon_end
                    mRNA.end = exon_end

                upstream_exon = Feature.from_fields(feature.seqname,
                                                    feature.source,
                                                    'exon',
                                                    exon_start,
                                                    exon_end,
                                                    feature.score,
                                                    feature.strand,
                                                    '.',
                                                    attribute_to_string(exon_attribute),
                                                   )

                upstream_exons.append(upstream_exon)

        processed_features = sorted(with_exons + upstream_exons)
        populate_connections(processed_features)

        return processed_features

    features = get_all_features(old_fn)
    for feature in features:
        feature.seqname = standardize_seqname(feature.seqname)
    
    with open(new_fn, 'w') as new_fh:
        original_lines = open(old_fn)

        for line in original_lines:
            if line.startswith('#'):
                new_fh.write(line)
            else:
                break
        
        new_fh.write(new_comments)
        
        for line in original_lines:
            if line.startswith('#'):
                break

        features = process_features_exons(features)

        for feature in features:
            new_fh.write(str(feature) + '\n')
        
        ## The first line to start with a comment after the features still needs
        ## to be written.
        #new_fh.write(line)

        #for line in original_lines:
        #    new_fh.write(line)
                   
def populate_connections(features):
    for f in features:
        f.children = set()
        f.parent = None

    id_to_object = {f.attribute['ID']: f for f in features if 'ID' in f.attribute}

    for f in features:
        f.populate_connections(id_to_object)

def get_all_features(gff_fn):
    ''' Ignore any line starting with '#' and all lines after any lines startin with '##FASTA'
    '''
    def relevant_lines(gff_fn):
        for line in open(gff_fn):
            if line.startswith('##FASTA'):
                break
            elif line.startswith('#'):
                continue
            else:
                yield line

    all_features = [Feature(line) for line in relevant_lines(gff_fn)]
    populate_connections(all_features)    

    return all_features

def get_CDSs(gff_fn, genome_dir, utr_fn):
    all_features = get_all_features(gff_fn)
    genes = gff_transcript.get_genes(gff_fn, genome_dir)
    translated_genes = [g for g in genes if g.CDSs]
    return translated_genes
