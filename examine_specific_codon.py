import Sequencing.sw as sw
import Sequencing.fastq as fastq
from collections import Counter

def count_triplets(reads, gene, codon_number):
    ''' Counts triplets that occur at codon_number in the coding sequence of
    gene in reads by local alignment of each read sequence to the context
    around the codon with the codon itself replaced with 'NNN'.
    '''
    sequence = gene.get_coding_sequence()

    start = (codon_number - 1) * 3
    around = 28
    context = sequence[start - around:start] + 'NNN' + sequence[start + 3: start + 3 + around]

    def relevant_alignment(alignment, context, seq):
        path = dict(alignment['path'])
        relevant = True
        triplet = ''
        for position in range(around, around + 3):
            if position not in path:
                relevant = False
            elif min(path[position], len(seq) - 1 - path[position]) < 3:
                relevant = False
            else:
                triplet += seq[path[position]]

        return relevant, triplet

    triplets = Counter()

    for read in reads:
        alignment = sw.generate_alignments(context, read.seq, 'overlap')[0]
        if len(alignment['path']) >= 24 and len(alignment['mismatches']) <= 6 and alignment['XO'] == 0:
            relevant, triplet = relevant_alignment(alignment, context, read.seq)
            if relevant:
                triplets[triplet] += 1

    return triplets
