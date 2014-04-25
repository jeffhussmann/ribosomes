import Bio.Data.CodonTable

forward_table = Bio.Data.CodonTable.standard_dna_table.forward_table
non_stop_codons = sorted(forward_table.keys())

def codons_from_seq(seq):
    for i in range(0, len(seq), 3):
        yield seq[i:i + 3]
