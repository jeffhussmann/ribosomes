import pysam
import os
import subprocess
from collections import Counter
from Circles import mapping_tools

#fastq_fn = '/home/jah/projects/arlen/experiments/guo_nature/Footprint_wild-type_runs1-2/data/SRR065774.fastq'
#index_prefix = '/home/jah/projects/arlen/data/organisms/mus_musculus/mm10/genome/genome'
#fastq_fn = '/home/jah/projects/arlen/experiments/belgium_3_5_14/wt/data/wt_cDNA.140219.HiSeq2500.FCB.lane1.R1.fastq'
fastq_fn = '/home/jah/projects/arlen/experiments/guydosh_cell/wild-type_short_footprints/data/SRR1042876.fastq'
index_prefix = '/home/jah/projects/arlen/data/organisms/saccharomyces_cerevisiae/EF4/genome/genome'

root, ext = os.path.splitext(fastq_fn)
small_fastq_fn = '{0}_small.fastq'.format(root)
small_sam_fn = '{0}_small.sam'.format(root)

head_command = ['head', '-n', '100000', fastq_fn]
subprocess.check_call(head_command, stdout=open(small_fastq_fn, 'w'))

mapping_tools.map_bowtie2(small_fastq_fn, index_prefix, small_sam_fn, local=True)

positions = [Counter() for i in range(40)]

randoms = Counter()

for read in pysam.Samfile(small_sam_fn): 
    if read.is_unmapped:
        continue
    randoms[read.qstart] += 1
    trimmed = read.seq[read.qend:]
    for p, b in zip(positions, trimmed):
        p[b] += 1

for p in positions:
    if not p:
        continue
    b, c = p.most_common(1)[0]
    fraction = float(c) / sum(p.values())
    print b, '{0:0.3f}'.format(fraction), sum(p.values())

inferred_adapter = ''.join(p.most_common(1)[0][0] for p in positions if p)
print inferred_adapter

os.remove(small_fastq_fn)
os.remove(small_sam_fn)
