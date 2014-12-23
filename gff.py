import gtf
import re
import pprint

ascii_hex_pattern = re.compile(r'%([0-9A-F]{2})')

def decoder(match):
    ascii_hex = match.group(1)
    return ascii_hex.decode('hex')

def sanitize(value):
    ''' Decodes escaped hex into ASCII and removes surrounding quotes. '''
    hex_replaced = ascii_hex_pattern.sub(decoder, value)
    return hex_replaced.strip('"')

class Feature(gtf.Feature):
    def __init__(self, line):
        super(Feature, self).__init__(line)
        self.parent = None
        self.children = set()

    def parse_attribute_string(self):
        fields = self.attribute_string.split(';')
        pairs = [field.split('=') for field in fields]
        parsed = {name: sanitize(value) for name, value in pairs}
        return parsed
    
    def populate_connections(self, name_to_object):
        parent_name = self.attribute.get('Parent')
        if parent_name:
            self.parent = name_to_object[parent_name]
            self.parent.children.add(self)

    def print_family(self, level=0):
        print '{0}{1}'.format('\t'*level, self)
        if level == 0:
            pprint.pprint(self.attribute)
        for child in self.children:
            child.print_family(level=level + 1)

    def is_ancestor_of(self, possible_descendant):
        if self == possible_descendant:
            return True
        else:
            return any(child.is_ancestor_of(possible_descendant) for child in self.children)

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
    
    name_to_object = {f.attribute['Name']: f for f in all_features
                      if 'Name' in f.attribute}

    for f in all_features:
        f.populate_connections(name_to_object)

    return all_features
