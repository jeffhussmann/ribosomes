import gtf
import bisect

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

    all_features = [gtf.parse_gtf_line(line, parse_gff_attribute) for line in relevant_lines(gff_fn)]
    return all_features

def make_starts_in_region_finder(gff_fn):
    def get_complete_start(feature):
        return (feature.seqname, feature.start)

    all_features = get_all_features(gff_fn)
    sorted_features = gtf.sort_features(all_features)
    starts = [get_complete_start(f) for f in sorted_features]

    def starts_in_region_finder(region_start, region_end):
        ''' Finds all features that start betwen region_start and region_end.
        '''
        in_region = []
        # Find the index of the first feature that starts after region start.
        i = bisect.bisect_right(starts, region_start)
        while i < len(starts) and get_complete_start(sorted_features[i]) < region_end:
            in_region.append(sorted_features[i])
            i += 1
        return in_region

    return starts_in_region_finder
