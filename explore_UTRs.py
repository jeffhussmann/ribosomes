from collections import Counter, defaultdict
from Circles import utilities
import gtf
import pyliftover
from glob import glob
import os.path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm
import scipy.stats

def maybe_int(string):
    try:
        value = int(string)
    except ValueError:
        value = string
    return value
    
def read_weinberg_file():
    fn = 'Cerevisiae_refFlat.txt'

    fh = open(fn)
    for i in range(3):
        fh.readline()

    labels = fh.readline().strip().split()

    genes = {}
    for line in fh:
        fields = line.strip().split('\t')
        pairs = zip(labels, map(maybe_int, fields))
        genes[fields[0]] = dict(pairs[2:])

    return genes

def read_nagalakshmi_file(fn):
    fh = open(fn)
    for i in range(2):
        fh.readline()

    labels = fh.readline().strip().split()

    genes = {}
    for line in fh:
        fields = line.strip('\n').split('\t')
        pairs = zip(labels, map(maybe_int, fields))
        genes[fields[0]] = dict(pairs[1:])

    return genes

def sanity_check_nagalakshmi_file():
    gtf_fn = '/home/jah/projects/arlen/data/organisms/saccharomyces_cerevisiae/SGD1.01/transcriptome/genes.gtf'
    CDSs = gtf.get_CDSs(gtf_fn)
    gtf_dict = {t.name: t for t in CDSs}
    genes = read_nagalakshmi_file('nagalakshmi_annotations.txt')
    discrepancies_after = {}
    for name in genes:
        if name not in gtf_dict:
            print name, 'not in gtf_dict'
            continue

        start, end = genes[name]['SGD_Start'] - 1, genes[name]['SGD_End'] - 1
        chrom = genes[name]['Chrom']
        if start > end:
            start, end = end, start

        if start != gtf_dict[name].start or end != gtf_dict[name].end:
            #print name, chrom
            #print start, end
            #print gtf_dict[name].start, gtf_dict[name].end
            #print
            if chrom not in discrepancies_after:
                discrepancies_after[chrom] = start
            else:
                discrepancies_after[chrom] = min(start, discrepancies_after[chrom])

    return discrepancies_after

def sanity_check_lifted_nagalakshmi_file(description):
    gtf_fn = '/home/jah/projects/arlen/data/organisms/saccharomyces_cerevisiae/EF4/transcriptome/genes.gtf'
    CDSs = gtf.get_CDSs(gtf_fn)
    gtf_dict = {t.name: t for t in CDSs}
    genes = read_nagalakshmi_file('nagalakshmi_annotations_lifted_{0}.txt'.format(description))
    discrepancies_after = defaultdict(list)
    for name in genes:
        if name not in gtf_dict:
            print name, 'not in gtf_dict'
            continue

        start, end = genes[name]['SGD_Start'], genes[name]['SGD_End']
        chrom = genes[name]['Chrom']
        if start > end:
            start, end = end, start

        end -= 1

        if start != gtf_dict[name].start or end != gtf_dict[name].end:
            #print name, chrom
            #print start, end
            #print gtf_dict[name].start, gtf_dict[name].end
            #raw_input()
            discrepancies_after[chrom].append(name)

    return discrepancies_after

def sanity_check_weinberg_file():
    gtf_fn = '/home/jah/projects/arlen/data/organisms/saccharomyces_cerevisiae/EF4/transcriptome/genes.gtf'
    CDSs = gtf.get_CDSs(gtf_fn)
    gtf_dict = {t.name: t for t in CDSs}
    genes = read_weinberg_file()
    discrepancies_after = defaultdict(list)
    for name in genes:
        if name not in gtf_dict:
            print name, 'not in gtf_dict'
            continue

        start, end = genes[name]['CdsStart'], genes[name]['CdsEnd']
        chrom = genes[name]['Chromosome']
        
        end -= 1

        if start != gtf_dict[name].start or end != gtf_dict[name].end:
            print name, chrom
            print start, end
            print gtf_dict[name].start, gtf_dict[name].end
            raw_input()
            discrepancies_after[chrom].append(name)

    return discrepancies_after

def liftover_new(chain_fn, description):
    original_fn = 'nagalakshmi_annotations.txt'
    lifted_fn = 'nagalakshmi_annotations_lifted_{0}.txt'.format(description)
    
    original_fh = open(original_fn)
    for i in range(2):
        original_fh.readline()

    labels = original_fh.readline().strip().split()

    keys_to_convert = ['SGD_Start', 'SGD_End', '5\'-UTR_Start', '3\'-UTR_End']

    lo = pyliftover.LiftOver(chain_fn)
    
    with open(lifted_fn, 'w') as lifted_fh:
        original_fh = open(original_fn)
        for i in range(2):
            lifted_fh.write(original_fh.readline())
        
        labels_line = original_fh.readline()
        lifted_fh.write(labels_line)
        labels = labels_line.strip().split()

        for line in original_fh:
            fields = line.strip('\n').split('\t')
            name = fields[0]
            if name == 'YBR013C':
                # This gets its 5' UTR deleted by liftover. Ignore it for now
                continue
            #if name == 'YJR122W':
            #    # This has its coding sequence misannotated in nagalakshmi.
            #    continue

            pairs = zip(labels, map(maybe_int, fields))
            gene = dict(pairs[1:])

            if gene['Chrom'] == 'chrMito':
                # Renamed in EF4, and not included in weinberg anyways.
                continue

            bad_lift = False
            for key in keys_to_convert:
                if gene[key] != '':
                    lift = lo.convert_coordinate(gene['Chrom'], gene[key] - 1)
                    if lift == []:
                        print gene, 'empty list'
                        bad_lift = True
                        break
                    seqname, coord, _, _ = lift[0]
                    gene[key] = coord
            if bad_lift:
                continue

            if gene['SGD_Start'] < gene['SGD_End']:
                # plus strand
                gene['SGD_End'] = gene['SGD_End'] + 1
            elif gene['SGD_Start'] > gene['SGD_End']:
                # minus strand
                gene['SGD_Start'] = gene['SGD_Start'] + 1
            else:
                raise ValueError(name)

            lifted_line = '\t'.join([name] + [str(gene[key]) for key in labels[1:]]) + '\n'
            lifted_fh.write(lifted_line)

def lift_all():
    chrs = [
        'chrI', 
        'chrII', 
        'chrIII', 
        'chrIV', 
        'chrIX', 
        'chrV', 
        'chrVI', 
        'chrVII', 
        'chrVIII', 
        'chrX', 
        'chrXI', 
        'chrXII', 
        'chrXIII', 
        'chrXIV', 
        'chrXV', 
        'chrXVI', 
    ]
    num_disagreeing = {key: [] for key in chrs}
    chain_fns = glob('/home/jah/projects/arlen/data/organisms/saccharomyces_cerevisiae/chains/*sanitized.over.chain')
    for chain_fn in sorted(chain_fns):
        head, tail = os.path.split(chain_fn)
        print tail[:3]
        if '57' not in tail[:3]:
            continue
        liftover_new(chain_fn, tail[:3])
        discs = sanity_check_lifted_nagalakshmi_file(tail[:3])
        for key in chrs:
            num_disagreeing[key].append(len(discs[key]))

    for key in chrs:
        print key, num_disagreeing[key]

    array = np.asarray([num_disagreeing[key] for key in chrs])
    return array

if __name__ == '__main__':
    weinberg_genes = read_weinberg_file()
    nagalakshmi_genes = read_nagalakshmi_file('nagalakshmi_annotations_lifted_V57.txt')

    # Check that the liftover worked.
    bad_genes = defaultdict(list)
    UTRs = {}
    weinberg_fives = []
    weinberg_threes = []
    nagalakshmi_fives = []
    nagalakshmi_threes = []

    for name in set(weinberg_genes) & set(nagalakshmi_genes):
        UTR_dict = {}
        weinberg_gene = weinberg_genes[name]
        nagalakshmi_gene = nagalakshmi_genes[name]


        if weinberg_genes[name]['Strand'] == '+':
            disagree = (weinberg_gene['CdsStart'] != nagalakshmi_gene['SGD_Start'] or
                        weinberg_gene['CdsEnd'] != nagalakshmi_gene['SGD_End'])
        else:
            disagree = (weinberg_gene['CdsStart'] != nagalakshmi_gene['SGD_End'] or
                        weinberg_gene['CdsEnd'] != nagalakshmi_gene['SGD_Start'])
        if disagree:
            bad_genes[weinberg_gene['Chromosome']].append(name)
        else:
            before = weinberg_gene['CdsStart'] - weinberg_gene['TxStart']
            after = weinberg_gene['TxEnd'] - weinberg_gene['CdsEnd']
            if weinberg_gene['Strand'] == '+':
                weinberg_five = before
                weinberg_three = after
            elif weinberg_gene['Strand'] == '-':
                weinberg_five = after
                weinberg_three = before

        nagalakshmi_five = nagalakshmi_gene['5\'-UTR_length']
        if isinstance(nagalakshmi_five, str):
            nagalakshmi_five = 0
        nagalakshmi_three = nagalakshmi_gene['3\'-UTR_length']
        if isinstance(nagalakshmi_three, str):
            nagalakshmi_three = 0

        UTRs[name] = {5: {'weinberg': weinberg_five,
                          'nagalakshmi': nagalakshmi_five,
                         },
                      3: {'weinberg': weinberg_five,
                          'nagalakshmi': nagalakshmi_five,
                         },
                     }
    
        weinberg_fives.append(weinberg_five)
        weinberg_threes.append(weinberg_three)
        nagalakshmi_fives.append(nagalakshmi_five)
        nagalakshmi_threes.append(nagalakshmi_three)

    for xs, ys in [(weinberg_fives, nagalakshmi_fives), (weinberg_threes, nagalakshmi_threes)]:
        fig, ax = plt.subplots()
        points = np.vstack([xs, ys])
        kernel = scipy.stats.gaussian_kde(points)
        colors = kernel(points)
        ax.scatter(xs, ys, c=colors, cmap=matplotlib.cm.jet, s=2, linewidths=(0,))
        ax.set_xlim(xmin=-20)
        ax.set_ylim(ymin=-20)

#fives = Counter()
#threes = Counter()
#
#five_length_to_name = defaultdict(list)
#three_length_to_name = defaultdict(list)
#
#for name in genes:
#    before = int(genes[name]['CdsStart']) - int(genes[name]['TxStart'])
#    after = int(genes[name]['TxEnd']) - int(genes[name]['CdsEnd'])
#    if genes[name]['Strand'] == '+':
#        fives[before] += 1
#        threes[after] += 1
#        five_length_to_name[before].append(name)
#        three_length_to_name[after].append(name)
#    elif genes[name]['Strand'] == '-':
#        fives[after] += 1
#        threes[before] += 1
#        five_length_to_name[after].append(name)
#        three_length_to_name[before].append(name)
#
#fives = utilities.counts_to_array(fives)
#threes = utilities.counts_to_array(threes)
