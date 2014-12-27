import gff
import Sequencing.genomes as genomes
import transcript

class Gene(transcript.Transcript):
    def __init__(self,
                 feature,
                 region_fetcher,
                ):
        self.name = feature.attribute['ID']
        self.region_fetcher = region_fetcher

        self.strand = feature.strand
        self.seqname = feature.seqname
        self.feature = feature

        descendants = feature.get_descendants()

        self.exons = sorted(c for c in descendants if c.feature == 'exon')

        self.CDSs = sorted(c for c in descendants if c.feature == 'CDS')

        if self.CDSs:
            if self.strand == '+':
                self.first_start_codon_position = min(cds.start for cds in self.CDSs)
                self.first_stop_codon_position = max(cds.end for cds in self.CDSs) - 2
            elif self.strand == '-':
                self.first_start_codon_position = max(cds.end for cds in self.CDSs)
                self.first_stop_codon_position = min(cds.start for cds in self.CDSs) + 2
        else:
            print self.name
            self.first_start_codon_position = None
            self.first_stop_codon_position = None

        self.start = min(exon.start for exon in self.exons)
        self.end = max(exon.end for exon in self.exons)

def get_genes(gff_fn, genome_dir):
    all_features = gff.get_all_features(gff_fn)
    region_fetcher = genomes.build_region_fetcher(genome_dir)
    genes = []
    for feature in all_features:
        if feature.feature in ['gene', 'transposable_element_gene']:
            gene = Gene(feature, region_fetcher)
            genes.append(gene)

    return genes
