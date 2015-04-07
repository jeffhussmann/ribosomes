import h5py
import numpy as np
import pausing

def read_file(file_name, specific_keys=None):
    with h5py.File(file_name, 'r') as hdf5_file:
        if specific_keys == None:
            specific_keys = hdf5_file
        arrays = {key: hdf5_file[key][()] for key in specific_keys}

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

        hdf5_file['nucleotide'] = stratified_mean_enrichments.arrays['nucleotide']
        hdf5_file['codon'] = stratified_mean_enrichments.arrays['codon']
        hdf5_file['dicodon'] = stratified_mean_enrichments.arrays['dicodon']
