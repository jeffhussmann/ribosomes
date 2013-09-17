from collections import namedtuple, defaultdict

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

def parse_attribute(attribute):
    fields = attribute.strip(';').split('; ')
    pairs = [field.split() for field in fields]
    parsed = {name: value.strip('"') for name, value in pairs}
    return parsed

def parse_gtf_line(line):
    gene = Gene._make(line.strip().split('\t'))
    start = int(gene.start) - 1
    end = int(gene.end) - 1
    if gene.frame != '.':
        frame = int(gene.frame)
    else:
        frame = gene.frame
    gene = gene._replace(start=start, end=end, frame=frame)
    return gene

def get_all_genes(gtf_fn):
    all_genes = [parse_gtf_line(line) for line in open(gtf_fn)]
    return all_genes

def get_contaminant_genes(gtf_fn):
    def is_contaminant(gene):
        return gene.source == 'rRNA' or gene.source == 'tRNA'
    
    all_genes = get_all_genes(gtf_fn)
    contaminant_genes = [gene for gene in all_genes if is_contaminant(gene)]

    return contaminant_genes

def get_all_CDSs(gtf_fn):
    def is_CDS(gene):
        return gene.source == 'protein_coding' and gene.feature == 'CDS'# and gene.end - gene.start > 1000
    
    all_genes = get_all_genes(gtf_fn)
    CDSs = [gene for gene in all_genes if is_CDS(gene)]

    return CDSs

def get_simple_CDSs(gtf_fn):
    ''' Returns all single exon CDSs that do not overlap any other CDS. '''
    CDSs = get_all_CDSs(gtf_fn)
    nonoverlapping = get_nonoverlapping_50(CDSs)
    single_exons = get_single_exons(CDSs)
    simple_CDSs = set(nonoverlapping) & set(single_exons)
    simple_CDSs = sort_genes(list(simple_CDSs))
    return simple_CDSs

def get_nonoverlapping(genes):
    ''' Returns all elements of genes that do not overlap any other element of
        genes. Requires genes to be sorted by start.
    '''
    def overlaps(left, right):
        return (right.seqname, right.start) <= (left.seqname, left.end)

    overlapping = set()
    nonoverlapping = []
    for i, gene in enumerate(genes):
        j = i + 1
        while j < len(genes) and overlaps(gene, genes[j]):
            overlapping.add(gene)
            overlapping.add(genes[j])
            j += 1
        if gene not in overlapping:
            nonoverlapping.append(gene)

    return nonoverlapping

def get_nonoverlapping_50(genes):
    ''' Returns all elements of genes that do not overlap any other element of
        genes. Requires genes to be sorted by start.
    '''
    def overlaps(left, right):
        return (right.seqname, right.start - 50) <= (left.seqname, left.end + 50)

    overlapping = set()
    nonoverlapping = []
    for i, gene in enumerate(genes):
        j = i + 1
        while j < len(genes) and overlaps(gene, genes[j]):
            overlapping.add(gene)
            overlapping.add(genes[j])
            j += 1
        if gene not in overlapping:
            nonoverlapping.append(gene)

    return nonoverlapping

def sort_genes(genes):
    def key(gene):
        return gene.seqname, gene.start, gene.feature

    return sorted(genes, key=key)

def get_single_exons(CDSs):
    ''' Returns all CDSs that only have a single exon. '''
    proteins = defaultdict(list)
    for gene in CDSs:
        attribute = parse_attribute(gene.attribute)
        proteins[attribute['protein_id']].append(gene)
    
    single_exons = [exons[0] for protein, exons in proteins.items() if len(exons) == 1]
    single_exons = sort_genes(single_exons)
    return single_exons
