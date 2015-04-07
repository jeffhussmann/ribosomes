import Bio.Data.CodonTable
import copy

nucleotide_order = 'TCAG'
nucleotide_to_index = {b: i for i, b in enumerate(nucleotide_order)}
def codon_sorting_key(codon):
    return [nucleotide_to_index[b] for b in codon]

forward_table = Bio.Data.CodonTable.standard_dna_table.forward_table

non_stop_codons = sorted(forward_table.keys(), key=codon_sorting_key)
stop_codons = sorted(Bio.Data.CodonTable.standard_dna_table.stop_codons, key=codon_sorting_key)

full_forward_table = copy.copy(forward_table)
for stop_codon in stop_codons:
    full_forward_table[stop_codon] = '*'

all_codons = sorted(non_stop_codons + stop_codons, key=codon_sorting_key)
codon_to_index = {codon: i for i, codon in enumerate(all_codons)}

full_back_table = {}
for codon in all_codons:
    amino_acid = full_forward_table[codon]
    if amino_acid not in full_back_table:
        full_back_table[amino_acid] = [codon]
    else:
        full_back_table[amino_acid].append(codon)

amino_acids = sorted(full_back_table.keys())

codon_to_amino_acid_and_index = {}
for amino_acid in full_back_table:
    for i, codon in enumerate(full_back_table[amino_acid]):
        codon_to_amino_acid_and_index[codon] = (amino_acid, i)

degeneracy = {amino_acid: len(full_back_table[amino_acid]) for amino_acid in amino_acids}

def codons_from_seq(seq):
    for i in range(0, len(seq), 3):
        yield seq[i:i + 3]
