import os
import glob

template = '''name {name}
group {group}
work_prefix /home/jah/
scratch_prefix /home/jah/scratch/
data_dir /home/jah/projects/ribosomes/experiments/{group}/{name}/data/
organism_dir /home/jah/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/
transcripts_file_name /home/jah/projects/ribosomes/experiments/weinberg/most_weinberg_transcripts.txt
relative_results_dir projects/ribosomes/experiments/{group}/{name}/results
adapter_type {adapter_type}
relevant_lengths {min_relevant_length},{max_relevant_length}
offset_type yeast
phiX_bowtie2_index_prefix /home/jah/bowtie2/phiX_doubled
'''

def make_descriptions(group,
                      adapter_type,
                      max_read_length=None,
                      min_relevant_length=14,
                      max_relevant_length=44,
                      num_pieces=12,
                      bartel_markers=False,
                      stephanie_markers=False,
                     ):
    prefix = '/home/jah/projects/ribosomes/experiments/{0}/'.format(group)
    bash_fn = '/home/jah/projects/ribosomes/code/all_{0}.sh'.format(group)
    bash_fh = open(bash_fn, 'w')

    dirs = [path for path in glob.glob('{}*'.format(prefix)) if os.path.isdir(path)]
    print group
    for d in sorted(dirs):
        _, name = os.path.split(d)
        print '\t', name
        job_dir = '{0}/job'.format(d)
        description_fn = '{0}/description.txt'.format(job_dir)
        if not os.path.isdir(job_dir):
            os.mkdir(job_dir)

        with open(description_fn, 'w') as description_fh:
            contents = template.format(group=group,
                                       name=name,
                                       adapter_type=adapter_type,
                                       min_relevant_length=min_relevant_length,
                                       max_relevant_length=max_relevant_length,
                                      )
            description_fh.write(contents)
            if max_read_length != None:
                description_fh.write('max_read_length {0}\n'.format(max_read_length))
            if bartel_markers:
                description_fh.write('synthetic_index_prefix /home/jah/projects/ribosomes/data/bartel_markers/bartel_markers\n')
            elif stephanie_markers:
                description_fh.write('synthetic_fasta /home/jah/projects/ribosomes/data/stephanie_markers/stephanie_markers.fa\n')

        bash_fh.write('echo {name}\n'.format(name=name))
        bash_fh.write('python ribosome_profiling_experiment.py --job_dir {0} launch --num_pieces {1}\n'.format(job_dir, num_pieces))

    return bash_fn

if __name__ == '__main__':
    kwargs = {'num_pieces': 12,
             }
    all_bash_fn = '/home/jah/projects/ribosomes/code/everything.sh'
    fns = []
    fns.append(make_descriptions('artieri', 'polyA', **kwargs))
    fns.append(make_descriptions('belgium_2013_08_06', 'truseq', **kwargs))
    fns.append(make_descriptions('belgium_2014_03_05', 'linker', **kwargs))
    fns.append(make_descriptions('belgium_2014_08_07', 'linker', stephanie_markers=True, **kwargs))
    fns.append(make_descriptions('belgium_2014_10_27', 'linker_local', stephanie_markers=True, **kwargs))
    fns.append(make_descriptions('belgium_2014_12_10', 'linker_local', stephanie_markers=True, **kwargs))
    fns.append(make_descriptions('dunn_elife', 'linker', **kwargs))
    fns.append(make_descriptions('gerashchenko_pnas', 'polyA', max_read_length=44, **kwargs))
    fns.append(make_descriptions('gerashchenko_nar', 'nothing', max_read_length=50, **kwargs))
    fns.append(make_descriptions('guydosh_cell', 'linker', **kwargs))
    fns.append(make_descriptions('ingolia_science', 'polyA', **kwargs))
    fns.append(make_descriptions('lareau_elife', 'linker', **kwargs))
    fns.append(make_descriptions('mcmanus_gr', 'linker', **kwargs))
    fns.append(make_descriptions('weinberg', 'weinberg', bartel_markers=True, **kwargs))
    fns.append(make_descriptions('zinshteyn_plos_genetics', 'polyA', **kwargs))
    fns.append(make_descriptions('pop_msb', 'linker', **kwargs))
    with open(all_bash_fn, 'w') as all_bash_fh:
        for fn in fns:
            all_bash_fh.write('echo {}\n'.format(fn))
            all_bash_fh.write('bash {}\n'.format(fn))
