import os
import glob

template = '''name {name}
work_prefix /home/jah/
scratch_prefix /home/jah/scratch/
data_dir /home/jah/projects/arlen/experiments/{family}/{name}/data/
organism_dir /home/jah/projects/arlen/data/organisms/saccharomyces_cerevisiae/EF4/
transcripts_file_name /home/jah/projects/arlen/experiments/weinberg/most_weinberg_transcripts.txt
relative_results_dir projects/arlen/experiments/{family}/{name}/results
adapter_type {adapter_type}
relevant_lengths 27,31
offset_type yeast
'''

def make_descriptions(family, adapter_type, max_read_length=None, num_pieces=16):
    prefix = '/home/jah/projects/arlen/experiments/{0}/'.format(family)
    bash_fn = '/home/jah/projects/arlen/code/all_{0}.sh'.format(family)
    bash_fh = open(bash_fn, 'w')

    dirs = [path for path in glob.glob('{}*'.format(prefix)) if os.path.isdir(path)]
    for d in sorted(dirs):
        _, name = os.path.split(d)
        job_dir = '{0}/job'.format(d)
        description_fn = '{0}/description.txt'.format(job_dir)
        if not os.path.isdir(job_dir):
            os.mkdir(job_dir)

        with open(description_fn, 'w') as description_fh:
            contents = template.format(family=family,
                                       name=name,
                                       adapter_type=adapter_type,
                                      )
            description_fh.write(contents)
            if max_read_length != None:
                description_fh.write('max_read_length {0}\n'.format(max_read_length))

        bash_fh.write('echo {name}\n'.format(name=name))
        bash_fh.write('python ribosome_profiling_experiment.py --job_dir {0} launch --num_pieces {1}\n'.format(job_dir, num_pieces))
