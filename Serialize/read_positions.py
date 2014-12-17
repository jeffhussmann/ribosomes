import h5py
import numpy as np
import positions

def is_an_int(string):
    try:
        int(string)
    except ValueError:
        return False
    return True

def build_gene(group, specific_keys=None):
    gene = {} 
    
    if specific_keys == None:
        specific_keys = group

    for key in specific_keys:
        dataset = group[key]
        data = dataset[...]

        attrs = dict(dataset.attrs.items())
        left_buffer = attrs.pop('left_buffer')
        right_buffer = attrs.pop('right_buffer')
        landmarks = attrs

        if is_an_int(key):
            key = int(key)

        gene[key] = positions.PositionCounts(landmarks,
                                             left_buffer,
                                             right_buffer,
                                             data=data,
                                            )
    return gene

def read_file(file_name, specific_keys=None):
    genes = {}
    with h5py.File(file_name, 'r') as hdf5_file:
        for gene_name in hdf5_file:
            genes[gene_name] = build_gene(hdf5_file[gene_name], specific_keys)
    return genes

def write_file(genes, file_name):
    with h5py.File(file_name, 'w') as hdf5_file:
        for gene_name in genes:
            gene_group = hdf5_file.create_group(gene_name)
            for key in genes[gene_name]:
                position_counts = genes[gene_name][key]
                # HDF5 names must be strings
                key = str(key)
                gene_group[key] = np.asarray(position_counts.data)
                gene_group[key].attrs['left_buffer'] = position_counts.left_buffer
                gene_group[key].attrs['right_buffer'] = position_counts.right_buffer
                for name, value in position_counts.landmarks.items():
                    gene_group[key].attrs[name] = value

def combine_data(first_genes, second_genes):
    for gene_name in second_genes:
        if gene_name not in first_genes:
            first_genes[gene_name] = second_genes[gene_name]
        else:
            for key in second_genes[gene_name]:
                if key not in first_genes[gene_name]:
                    first_genes[gene_name][key] = second_genes[gene_name][key]
                else:
                    first_genes[gene_name][key] += second_genes[gene_name][key] 
    return first_genes
