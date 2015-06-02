import h5py
import numpy as np
import pausing

def read_file(file_name, specific_cutoffs=['0.10'], specific_keys=['codon']):
    with h5py.File(file_name, 'r') as hdf5_file:
        if specific_cutoffs == None:
            specific_cutoffs = hdf5_file

        arrays = {}
        for cutoff in specific_cutoffs:
            if specific_keys == None:
                specific_keys = hdf5_file[cutoff].keys()
            arrays[cutoff] = {key: hdf5_file[cutoff][key][()] for key in specific_keys}

        num_before = hdf5_file.attrs['num_before']
        num_after = hdf5_file.attrs['num_after']
        stratified_mean_enrichments = pausing.StratifiedMeanEnrichments(num_before,
                                                                        num_after,
                                                                        arrays,
                                                                       )

    return stratified_mean_enrichments

def write_file(stratified_mean_enrichments, file_name):
    with h5py.File(file_name, 'w') as hdf5_file:
        hdf5_file.attrs['num_before'] = stratified_mean_enrichments.num_before
        hdf5_file.attrs['num_after'] = stratified_mean_enrichments.num_after

        for cutoff in stratified_mean_enrichments.arrays:
            hdf5_file.create_group(cutoff)
            for key in stratified_mean_enrichments.arrays[cutoff]:
                hdf5_file[cutoff][key] = stratified_mean_enrichments.arrays[cutoff][key]
