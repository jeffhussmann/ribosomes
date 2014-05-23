from collections import Counter, defaultdict
from Circles import utilities
import gtf
import gff
import pyliftover
from glob import glob
import os.path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm
import scipy.stats
from mpl_toolkits.axes_grid1 import make_axes_locatable

def maybe_int(string):
    try:
        value = int(string)
    except ValueError:
        value = string
    return value
    
def read_weinberg_file():
    fn = '/home/jah/projects/arlen/data/organisms/saccharomyces_cerevisiae/EF4/transcriptome/Cerevisiae_refFlat.txt'

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

def sanity_check_lifted_nagalakshmi_file(fn):
    gtf_fn = '/home/jah/projects/arlen/data/organisms/saccharomyces_cerevisiae/EF4/transcriptome/genes.gtf'
    CDSs = gtf.get_CDSs(gtf_fn)
    gtf_dict = {t.name: t for t in CDSs}
    genes = read_nagalakshmi_file(fn)
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

def do_liftover(chain_fn, description):
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
        do_liftover(chain_fn, tail[:3])
        discs = sanity_check_lifted_nagalakshmi_file(tail[:3])
        for key in chrs:
            num_disagreeing[key].append(len(discs[key]))

    for key in chrs:
        print key, num_disagreeing[key]

    array = np.asarray([num_disagreeing[key] for key in chrs])
    return array

def scatter_with_hists_colors(x_list, y_list, x_label, y_label, title):
    sampled_points = np.vstack([x_list[:10000], y_list[:10000]])
    points = np.vstack([x_list, y_list])
    kernel = scipy.stats.gaussian_kde(sampled_points)
    colors = kernel(points)

    fig, ax_scatter = plt.subplots()
    ax_scatter.scatter(x_list, y_list, c=colors, cmap=matplotlib.cm.jet, s=4, linewidths=(0.1,))
    ax_scatter.set_aspect(1.)

    divider = make_axes_locatable(ax_scatter)
    ax_hist_x = divider.append_axes('top', 1.2, pad=0.1, sharex=ax_scatter)
    ax_hist_y = divider.append_axes('right', 1.2, pad=0.1, sharey=ax_scatter)

    plt.setp(ax_hist_x.get_xticklabels() + ax_hist_x.get_yticklabels(),
             visible=False)
    plt.setp(ax_hist_x.get_xticklines() + ax_hist_x.get_yticklines(),
             visible=False)
    plt.setp(ax_hist_y.get_xticklabels() + ax_hist_y.get_yticklabels(),
             visible=False)
    plt.setp(ax_hist_y.get_xticklines() + ax_hist_y.get_yticklines(),
             visible=False)
    
    ax_scatter.set_xlabel(x_label)
    ax_scatter.set_ylabel(y_label)

    ax_hist_x.hist(x_list, histtype='step', bins=100)
    ax_hist_y.hist(y_list, histtype='step', bins=100, orientation='horizontal')
    
    #ax_scatter.set_xlim(left=-20)
    #ax_scatter.set_ylim(bottom=-20)

    fig.suptitle(title)

def check_for_overlap():
    weinberg_genes = read_weinberg_file()

    gff_fn = '/home/jah/projects/arlen/data/organisms/saccharomyces_cerevisiae/EF4/saccharomyces_cerevisiae.gff'
    starts_in_region_finder = gff.make_starts_in_region_finder(gff_fn)

    print '- strand'

    for name, gene in weinberg_genes.iteritems():
        if gene['Strand'] != '-':
            continue

        region_start = (gene['Chromosome'], gene['CdsEnd'])
        region_end = (gene['Chromosome'], gene['TxEnd'])

        starts_in_region = starts_in_region_finder(region_start, region_end)
        if starts_in_region != []:
            print name
            print '{0}:{1}-{2}'.format(gene['Chromosome'], gene['CdsEnd'], gene['TxEnd'])
            for feature in starts_in_region:
                print feature.attribute['Name'], feature.start
            raw_input()

    print '+ strand'
    
    #sorted_transcripts = gtf.sort_transcripts(all_transcripts, by_end=True)
    #name_to_sorted_index = {t.name: i for i, t in enumerate(sorted_transcripts)}

    #for name, gene in weinberg_genes.iteritems():
    #    if name not in name_to_sorted_index:
    #        continue
    #    i = name_to_sorted_index[name]
    #    if (sorted_transcripts[i - 1].end <= gene['CdsStart'] and
    #        gene['TxStart'] <= sorted_transcripts[i - 1].end and
    #        gene['Strand'] == '+' and sorted_transcripts[i - 1].strand == '+') or \
    #        name == 'YNL054W':

    #        print name
    #        print gene['TxStart']
    #        print sorted_transcripts[i]
    #        print sorted_transcripts[i - 1]

if __name__ == '__main__':
    nagalakshmi_lifted_fn = '/home/jah/projects/arlen/data/organisms/saccharomyces_cerevisiae/EF4/transcriptome/nagalakshmi_annotations_lifted_V57.txt'
    weinberg_genes = read_weinberg_file()
    nagalakshmi_genes = read_nagalakshmi_file(nagalakshmi_lifted_fn)

    # Check that the liftover worked.
    bad_genes = defaultdict(list)
    UTRs_dict = {}
    weinberg_fives = []
    weinberg_threes = []
    nagalakshmi_fives = []
    nagalakshmi_threes = []
    
    discrepancies = sanity_check_lifted_nagalakshmi_file(nagalakshmi_lifted_fn) 
    changed_names = set(reduce(lambda x, y: x + y, discrepancies.values(), []))

    valid_genes = set(weinberg_genes) & set(nagalakshmi_genes) - changed_names
    for name in set(weinberg_genes) & set(nagalakshmi_genes):
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

        UTRs_dict[name] = {5: {'weinberg': weinberg_five,
                               'nagalakshmi': nagalakshmi_five,
                              },
                           3: {'weinberg': weinberg_three,
                               'nagalakshmi': nagalakshmi_three,
                              },
                          }

        weinberg_fives.append(weinberg_five)
        weinberg_threes.append(weinberg_three)
        nagalakshmi_fives.append(nagalakshmi_five)
        nagalakshmi_threes.append(nagalakshmi_three)

    #for title, xs, ys in [('5\' UTR lengths', weinberg_fives, nagalakshmi_fives),
    #                      ('3\' UTR lengths', weinberg_threes, nagalakshmi_threes),
    #                     ]:
    #    scatter_with_hists_colors(xs, ys, 'Weinberg/Arribere annotation', 'Nagalakshi annotation', title)

    #for name, UTR in UTRs_dict.items():
    #    if UTR[5]['nagalakshmi'] > 0 and UTR[5]['weinberg'] - UTR[5]['nagalakshmi'] > 200:
    #        print name, 'changed' if name in changed_names else ''
    #        print weinberg_genes[name]['Chromosome'], weinberg_genes[name]['TxStart']
    #        print UTR[5]

    sorted_UTRs = sorted(UTRs_dict.items(), key=lambda (name, UTR): UTR[5]['weinberg'], reverse=True)
    for name, UTR in sorted_UTRs[:20]:
        print name
        print '{0}:{1}-{2}'.format(weinberg_genes[name]['Chromosome'], weinberg_genes[name]['CdsEnd'], weinberg_genes[name]['TxEnd'])
        print UTR[5]

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
