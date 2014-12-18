import gtf
import bisect

class GFF_line(object):
    def __init__(self, line):
        self.feature = gtf.parse_gtf_line(line, parse_gff_attribute)
        self.parent = None
        self.children = set()

    def populate_connections(self, name_to_object):
        parent_name = self.feature.attribute.get('Parent')
        if parent_name:
            self.parent = name_to_object[parent_name]
            self.parent.children.add(self)

def parse_gff_attribute(attribute):
    fields = attribute.split(';')
    pairs = [field.split('=') for field in fields]
    parsed = {name: value.strip('"') for name, value in pairs}
    return parsed

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

    all_features = [GFF_line(line) for line in relevant_lines(gff_fn)]
    name_to_object = {}
    for f in all_features:
        name = f.feature.attribute.get('Name')
        if name:
            name_to_object[name] = f

    for f in all_features:
        f.populate_connections(name_to_object)

    return all_features
