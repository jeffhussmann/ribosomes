import Bio.SeqIO
import gtf
import mapping
import mutations

gtf_fn = '/home/jah/projects/arlen/data/organisms/saccharomyces_cerevisiae/transcriptome/genes.gtf'
genes = gtf.get_simple_CDSs(gtf_fn)
genome = mapping.load_genome('/home/jah/projects/arlen/data/organisms/saccharomyces_cerevisiae/genome/', explicit_path=True)

for gene in genes:
    if gene.seqname == 'MT':
        continue
    print gene.source, gene.feature, gene.seqname
    print gene.attribute
    if gene.strand == '+':
        # gene.end is the last base before the stop codon
        seq = genome[gene.seqname][gene.start:gene.end + 4]
        Bio.Seq.translate(seq, cds=True)
    elif gene.strand == '-':
        # gene.start is the first base after the stop codon
        rc_seq = genome[gene.seqname][gene.start - 3:gene.end + 1]
        seq = mutations.reverse_complement(rc_seq)
        Bio.Seq.translate(seq, cds=True)

