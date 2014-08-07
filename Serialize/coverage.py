import h5py
import numpy as np
from itertools import izip

def read_file(file_name):
    counts = {}
    with h5py.File(file_name, 'r') as hdf5_file:
        for name in hdf5_file:
            # str is to convert from unicode
            counts[str(name)] = hdf5_file[name][...]
    return counts

def write_file(counts, file_name):
    with h5py.File(file_name, 'w') as hdf5_file:
        for name in counts:
            hdf5_file[name] = counts[name]

def combine_data(first_counts, second_counts):
    assert first_counts.viewkeys() == second_counts.viewkeys()

    counts = {name: first_counts[name] + second_counts[name] for name in first_counts}

    return counts
