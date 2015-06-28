import h5py
import numpy as np
import pausing

def read_file(file_name):
    hdf5_file = h5py.File(file_name, 'r')

    num_around = hdf5_file.attrs['num_around']
    stratified_mean_enrichments = pausing.StratifiedMeanEnrichments(num_around, hdf5_file)

    return stratified_mean_enrichments

def write_file(stratified_mean_enrichments, file_name):
    with h5py.File(file_name, 'w') as hdf5_file:
        hdf5_file.attrs['num_around'] = stratified_mean_enrichments.num_around

        for condition in stratified_mean_enrichments.arrays:
            hdf5_file.create_group(condition)
            for key in stratified_mean_enrichments.arrays[condition]:
                hdf5_file[condition][key] = stratified_mean_enrichments.arrays[condition][key]
