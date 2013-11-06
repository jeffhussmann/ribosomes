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
    # Convert from 1-based indexing to 0-based
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

def get_tRNA_genes(gtf_fn):
    all_genes = get_all_genes(gtf_fn)
    tRNA_genes = [gene for gene in all_genes if gene.source == 'tRNA']
    return tRNA_genes

def get_rRNA_genes(gtf_fn):
    all_genes = get_all_genes(gtf_fn)
    rRNA_genes = [gene for gene in all_genes if gene.source == 'rRNA']
    return rRNA_genes

def get_all_CDSs(gtf_fn):
    def is_CDS(gene):
        return gene.source == 'protein_coding' and gene.feature == 'CDS'
    
    all_genes = get_all_genes(gtf_fn)
    CDSs = [gene for gene in all_genes if is_CDS(gene)]

    return CDSs

def sort_genes(genes):
    def key(gene):
        return gene.seqname, gene.start, gene.feature

    return sorted(genes, key=key)

def get_extent_by_name(gtf_fn, name):
    all_genes = get_all_genes(gtf_fn)
    entries = [gene for gene in all_genes
               if parse_attribute(gene.attribute)['gene_id'] == name]
    start_codon = [entry for entry in entries if entry.feature == 'start_codon'][0]
    stop_codon = [entry for entry in entries if entry.feature == 'stop_codon'][0]
    if any(entry.strand == '-' for entry in entries):
        raise RuntimeError, 'minus strand NYI'
    start = start_codon.start
    # Haven't decided what the convetion should be for which base is the end
    end = stop_codon.end
    seqname = start_codon.seqname
    strand = start_codon.strand
    return seqname, strand, start, end

def get_nonoverlapping(genes, edge_buffer=0):
    ''' Returns all elements of genes that do not overlap any other element of
        genes..
    '''
    def overlaps(l, r):
        ''' r is guaranteed to start after l because it comes after l
            in a sorted list of genes, so r and l overlap iff r
            hasn't ended by the time l starts.
        '''
        return (r.seqname, r.start - edge_buffer) <= (l.seqname, l.end + edge_buffer)

    overlapping = set()
    nonoverlapping = []

    genes = sort_genes(genes)
    for i, gene in enumerate(genes):
        j = i + 1
        while j < len(genes) and overlaps(gene, genes[j]):
            overlapping.add(gene)
            overlapping.add(genes[j])
            j += 1
        if gene not in overlapping:
            nonoverlapping.append(gene)

    return nonoverlapping

def get_single_exons(CDSs):
    ''' Returns all CDSs that only have a single exon. '''
    proteins = defaultdict(list)
    for gene in CDSs:
        attribute = parse_attribute(gene.attribute)
        proteins[attribute['protein_id']].append(gene)
    
    single_exons = [exons[0] for protein, exons in proteins.items() if len(exons) == 1]
    single_exons = sort_genes(single_exons)
    return single_exons

def get_simple_CDSs(gtf_fn):
    ''' Returns all single exon CDSs that do not overlap any other CDS. '''
    CDSs = get_all_CDSs(gtf_fn)
    nonoverlapping = get_nonoverlapping(CDSs, edge_buffer=50)
    single_exons = get_single_exons(CDSs)
    simple_CDSs = set(nonoverlapping) & set(single_exons)
    simple_CDSs = sort_genes(list(simple_CDSs))
    return simple_CDSs
