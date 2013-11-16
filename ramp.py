import Serialize
import numpy as np

#rpf_fn = '/home/jah/projects/arlen/experiments/belgium_8_6_13/R98S_cDNA_sample/results/R98S_cDNA_sample_rpf_positions.txt'
rpf_fn = '/home/jah/projects/arlen/experiments/gerashchenko_pnas/Initial_rep2_foot/results/Initial_rep2_foot_rpf_positions.txt'

rpf_positions_list = Serialize.read_file(rpf_fn, 'rpf_positions')

expected_counts = np.zeros(5000)
actual_counts = np.zeros(5000)

for (name, length), counts in zip(*rpf_positions_list):
    if length < 600:
        continue
    A_site_offset = 15
    length_28s = counts[0][100 - A_site_offset:-(50 + A_site_offset):3]
    r_g = max(1, length_28s.sum())
    
    actual_counts[:length / 3] += np.true_divide(length_28s, r_g)

    num_codons = len(length_28s)
    uniform_counts = np.ones(num_codons) * 1 / num_codons

    expected_counts[:length / 3] += uniform_counts
