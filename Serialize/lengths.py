import h5py

def read_file(file_name):
    lengths = {}
    with h5py.File(file_name, 'a') as hdf5_file:
        for name in hdf5_file:
            lengths[name] = hdf5_file[name][...]

    return lengths

def write_file(lengths, file_name):
    with h5py.File(file_name, 'a') as hdf5_file:
        for name in lengths:
            if name in hdf5_file:
                hdf5_file[name][...] = lengths[name]
            else:
                hdf5_file[name] = lengths[name]

def combine_data(first_lengths, second_lengths):
    for name in second_lengths:
        if name not in first_lengths:
            first_lengths[name] = second_lengths[name]
        else:
            first_lengths[name] += second_lengths[name]

    return first_lengths
