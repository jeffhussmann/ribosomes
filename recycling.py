import Bio.Seq
import Bio.Data.CodonTable
import gtf
import mapping
import mutations
from collections import defaultdict
from itertools import izip

forward_table = Bio.Data.CodonTable.standard_dna_table.forward_table
back_table = defaultdict(list)
for codon in sorted(forward_table):
    amino_acid = forward_table[codon]
    back_table[amino_acid].append(codon)

def seq_to_codons(seq):
    assert len(seq) % 3 == 0
    nucs = iter(seq)
    triplets = izip(*[nucs]*3)
    for triplet in triplets:
        yield ''.join(triplet)

def get_amino_acid_locations(gene, genome):
    amino_acid_locations = defaultdict(list)
    if gene.seqname == 'MT':
        # Ignore these for now - diffent genetic code and tRNAs presumably means
        # different translation dynamics
        return None
    try:
        if gene.strand == '+':
            # gene.end is the last base before the stop codon
            seq = genome[gene.seqname][gene.start:gene.end + 4]
            translation = Bio.Seq.translate(seq, cds=True)
        elif gene.strand == '-':
            # gene.start is the first base after the (rc of the) stop codon
            # gene.end is the last base of the (rc of the) start codon
            rc_seq = genome[gene.seqname][gene.start - 3:gene.end + 1]
            seq = mutations.reverse_complement(rc_seq)
            translation = Bio.Seq.translate(seq, cds=True)
    except Bio.Seq.CodonTable.TranslationError, err:
        print err
        print gene.source, gene.feature, gene.seqname
        print gene.attribute
        return None
            
    codons = seq_to_codons(seq)
    for c, (codon, amino_acid) in enumerate(zip(codons, translation)):
        location = (c, codon)
        amino_acid_locations[amino_acid].append(location)

    return amino_acid_locations
