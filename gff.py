import gtf
import urllib
import pprint
import call_UTRs
import transcript
import interval_tree
from collections import Counter

class Feature(gtf.Feature):
    def __init__(self, line=None):
        super(Feature, self).__init__(line)
        self.parent = None
        self.children = set()

    def parse_attribute_string(self):
        if self.attribute_string == '.':
            parsed = {}
        else:
            fields = self.attribute_string.split(';')
            pairs = [field.split('=') for field in fields]
            parsed = {name: urllib.unquote(value).strip('"') for name, value in pairs}

        self.attribute = parsed
    
    def unparse_attribute_string(self):
        entries = []
        for key, value in sorted(self.attribute.items()):
            key = urllib.quote(str(key), safe='')
            value = urllib.quote(str(value), safe='')
            entry = '{0}={1}'.format(key, value)
            entries.append(entry)

        self.attribute_string = ';'.join(entries)

    def populate_connections(self, id_to_object):
        parent_id = self.attribute.get('Parent')
        if parent_id:
            self.parent = id_to_object[parent_id]
            self.parent.children.add(self)

    @property
    def descendants(self):
        descendants = set(self.children)
        for child in self.children:
            descendants.update(child.descendants)
        return descendants

    def print_family(self, level=0):
        print '{0}{1}'.format('\t'*min(level, 1), self)
        if level == 0:
            pprint.pprint(self.attribute)
        for child in self.children:
            child.print_family(level=level + 1)

def convert_SGD_file(old_fn, new_fn):
    ''' Replaces seqnames in a file from SGD to agree with those from ENSEMBL - 
    removes 'chr' from beginning and replaces mt with Mito.
    Adds exon features to genes.
    '''
    def standardize_seqname(seqname):
        seqname = seqname[3:]
        if seqname == 'mt':
            seqname = 'Mito'
        return seqname

    def add_exons(features):
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

                exon = Feature.from_fields(feature.seqname,
                                           feature.source,
                                           'exon',
                                           feature.start,
                                           feature.end,
                                           feature.score,
                                           feature.strand,
                                           '.',
                                           '.',
                                          )
                
                exon.attribute = {'ID': exon_id,
                                  'Parent': parent_id,
                                 }
                exon.unparse_attribute_string()
                
                feature.attribute['Parent'] = exon_id
                feature.unparse_attribute_string()
                
                new_features.append(exon)

        with_exons = sorted(features + new_features)

        populate_all_connections(with_exons)

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
                                                    '.',
                                                   )
                upstream_exon.attribute = {'ID': exon_id,
                                           'Parent': parent_id,
                                          }
                upstream_exon.unparse_attribute_string()

                upstream_exons.append(upstream_exon)

        processed_features = sorted(with_exons + upstream_exons)
        populate_all_connections(processed_features)

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
        
        new_fh.write('''\
# seqnames replaced.
# CDS's have been wrapped in exons.
# Gene's with a five_prime_UTR_intron have had their mRNA extended upstream of the intron.
''')

        features = add_exons(features)

        for feature in features:
            new_fh.write(str(feature) + '\n')
        
        #for line in original_lines:
        #    if line.startswith('#'):
        #        break

        ## The first line to start with a comment after the features still needs
        ## to be written.
        #new_fh.write(line)

        #for line in original_lines:
        #    new_fh.write(line)

def extend_UTRs(old_fn, new_fn, UTR_fn, genome_dir):
    UTR_boundaries = call_UTRs.read_UTR_file(UTR_fn)
    all_features = get_all_features(old_fn)

    genes = transcript.get_gff_transcripts(all_features, '/dev/null')

    genes = {g.name: g for g in genes}

    smallest_start = lambda e: e.start
    largest_end = lambda e: e.end

    for name in UTR_boundaries:
        gene = genes[name]
        if len(gene.mRNAs) != 1:
            raise ValueError('not exactly one mRNA')

        _, _, five_pos, three_pos = UTR_boundaries[name]

        if gene.strand == '+':
            leftmost_exon = min(gene.exons, key=smallest_start)
            leftmost_exon.start = five_pos
            gene.mRNAs[0].start = five_pos
            gene.top_level_feature.start = five_pos

            rightmost_exon = max(gene.exons, key=largest_end)
            rightmost_exon.end = three_pos
            gene.mRNAs[0].end = three_pos
            gene.top_level_feature.end = three_pos
        elif gene.strand == '-':
            leftmost_exon = max(gene.exons, key=largest_end)
            leftmost_exon.end = five_pos
            gene.mRNAs[0].end = five_pos
            gene.top_level_feature.end = five_pos

            rightmost_exon = min(gene.exons, key=smallest_start)
            rightmost_exon.start = three_pos
            gene.mRNAs[0].start = three_pos
            gene.top_level_feature.start = three_pos
    
    with open(new_fn, 'w') as new_fh:
        original_lines = open(old_fn)

        for line in original_lines:
            if line.startswith('#'):
                new_fh.write(line)
            else:
                break
        
        new_fh.write('''\
# UTRs have been extended.
# Top-level features have had distances to the closest other top-level features annotated. 
''')

        mark_nearby(all_features, genome_dir)

        for feature in all_features:
            new_fh.write(str(feature) + '\n')

def populate_all_connections(features):
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
    populate_all_connections(all_features)    

    return all_features

def get_top_level_features(features):
    top_level_features = [f for f in features if f.parent == None]
    return top_level_features

def get_CDSs(gff_fn, genome_dir, annotate_nearby=False):
    all_features = get_all_features(gff_fn)
    if annotate_nearby:
        mark_nearby(all_features, genome_dir)
    genes = transcript.get_gff_transcripts(all_features, genome_dir)
    translated_genes = [g for g in genes if g.CDSs]
    return translated_genes

def get_noncoding_RNA_transcripts(gff_fn):
    all_features = get_all_features(gff_fn)
    genes = transcript.get_gff_transcripts(all_features, '/dev/null')
    rRNA_transcripts = []
    tRNA_transcripts = []
    other_ncRNA_transcripts = []
    for gene in genes:
        if gene.top_level_feature.feature == 'rRNA':
            rRNA_transcripts.append(gene)
        elif gene.top_level_feature.feature == 'tRNA':
            tRNA_transcripts.append(gene)
        elif 'RNA' in gene.top_level_feature.feature:
            other_ncRNA_transcripts.append(gene)
    return rRNA_transcripts, tRNA_transcripts, other_ncRNA_transcripts

def mark_nearby(all_features, genome_dir):
    def is_nontrivial(possible):
        if possible.feature in ['chromosome', 'landmark', 'ARS', 'region']:
            return False
        elif possible.attribute.get('orf_classification') == 'Dubious':
            return False
        else:
            return True
    
    top_level_features = get_top_level_features(all_features)
    nontrivial_features = filter(is_nontrivial, top_level_features)
    overlap_finder = interval_tree.NamedOverlapFinder(nontrivial_features, genome_dir)

    def is_relevant_to(possible, main):
        if possible == main:
            return False
        elif main.strand != '.' and possible.strand != '.' and main.strand != possible.strand:
            return False
        else:
            return True

    for top_level in top_level_features:
        overlapping = overlap_finder.overlapping(top_level.seqname,
                                                 top_level.start,
                                                 top_level.end,
                                                )
        overlapping = [f for f in overlapping if is_relevant_to(f, top_level)]

        before = overlap_finder.find_closest_before(top_level.seqname,
                                                    top_level.strand,
                                                    top_level.start,
                                                   )
        after = overlap_finder.find_closest_after(top_level.seqname,
                                                  top_level.strand,
                                                  top_level.end,
                                                 )

        top_level.attribute['closest_left'] = before[0].end
        top_level.attribute['closest_right'] = after[0].start
        top_level.attribute['overlapping'] = len(overlapping)

        top_level.unparse_attribute_string()

if __name__ == '__main__':
    boundaries_fn = '/home/jah/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/transcriptome/inferred_UTR_lengths.txt'
    genome_dir = '/home/jah/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/genome/'
    original_gff_fn = '/home/jah/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/saccharomyces_cerevisiae.gff'
    with_exons_fn = '/home/jah/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/transcriptome/genes_with_exons.gff'
    final_fn = '/home/jah/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/transcriptome/genes.gff'

    ##convert_SGD_file(original_gff_fn, with_exons_fn)

    ## This assumes that the experiments used in call_UTR_boundaries have been
    ## run using with_exons_fn as the source of gene models.
    #call_UTRs.call_UTR_boundaries(boundaries_fn)

    extend_UTRs(with_exons_fn, final_fn, boundaries_fn, genome_dir)
