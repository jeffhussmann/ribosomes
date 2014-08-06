import pysam
import os
import subprocess
from collections import Counter
from Sequencing import mapping_tools

#fastq_fn = '/home/jah/projects/arlen/experiments/guo_nature/Footprint_wild-type_runs1-2/data/SRR065774.fastq'
#index_prefix = '/home/jah/projects/arlen/data/organisms/mus_musculus/mm10/genome/genome'
#fastq_fn = '/home/jah/projects/arlen/experiments/belgium_3_5_14/wt/data/wt_cDNA.140219.HiSeq2500.FCB.lane1.R1.fastq'
#fastq_fn = '/home/jah/projects/arlen/experiments/dunn_elife/YCF182_110222_HiSeq.fq'
#fastq_fn = '/home/jah/projects/arlen/experiments/lareau_elife/Cycloheximide_replicate_1/data/SRR1363415.fastq'
#fastq_fn = '/home/jah/projects/arlen/experiments/arribere_gr/S288C_TLSeq2/data/SRR825166.fastq'
#fastq_fn = '/home/jah/projects/arlen/experiments/baudin-baillieu_cell_reports/traductome_PSI-_rep_2/data/SRR594901.fastq'
fastq_fn = '/home/jah/projects/arlen/experiments/baudin-baillieu_cell_reports/Ribo-seq_[PSI+]_rep1/data/SRR1190356.fastq'
index_prefix = '/home/jah/projects/arlen/data/organisms/saccharomyces_cerevisiae/EF4/genome/genome'

root, ext = os.path.splitext(fastq_fn)
small_fastq_fn = '{0}_small.fastq'.format(root)
small_sam_fn = '{0}_small.sam'.format(root)

head_command = ['head', '-n', '100000', fastq_fn]
subprocess.check_call(head_command, stdout=open(small_fastq_fn, 'w'))

mapping_tools.map_bowtie2(small_fastq_fn, index_prefix, small_sam_fn, seed_length=12, local=True)

positions = [Counter() for i in range(40)]

qlens = Counter()

for read in pysam.Samfile(small_sam_fn): 
    if read.is_unmapped:
        continue
    qlens[read.qlen] += 1
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

for i in range(max(qlens)):
    print '{0}:\t{1}'.format(i, qlens[i])

os.remove(small_fastq_fn)
#os.remove(small_sam_fn)
