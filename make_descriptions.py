import os
import glob
import numpy as np

template = '''\
name {name}
group {group}
work_prefix {home}
scratch_prefix {home}/scratch/
data_dir {home}/projects/ribosomes/experiments/{group}/{name}/data/
organism_dir {home}/projects/ribosomes/data/organisms/saccharomyces_cerevisiae/EF4/
transcripts_file_name {home}/projects/ribosomes/experiments/weinberg/most_weinberg_transcripts.txt
relative_results_dir projects/ribosomes/experiments/{group}/{name}/results
adapter_type {adapter_type}
relevant_lengths {min_relevant_length},{max_relevant_length}
offset_type {offset_type}
phiX_bowtie2_index_prefix {home}/bowtie2/phiX_doubled
'''

def make_descriptions(group,
                      adapter_type,
                      max_read_length=None,
                      min_relevant_length=14,
                      max_relevant_length=44,
                      name_to_offset_type=lambda name: 'yeast',
                      name_to_adapter_type=None,
                      num_pieces=12,
                      markers_fn=None,
                      codons_to_examine=[],
                     ):
    prefix = '{home}/projects/ribosomes/experiments/{0}/'.format(group, home=os.environ['HOME'])
    bash_fn = '{home}/projects/ribosomes/code/all_{0}.sh'.format(group, home=os.environ['HOME'])
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
            if name_to_adapter_type:
                adapter_type = name_to_adapter_type(name)

            contents = template.format(group=group,
                                       name=name,
                                       adapter_type=adapter_type,
                                       min_relevant_length=min_relevant_length,
                                       max_relevant_length=max_relevant_length,
                                       home=os.environ['HOME'],
                                       offset_type=name_to_offset_type(name),
                                      )

            if codons_to_examine:
                locii_string = ';'.join('{0},{1}'.format(g_n, c_n) for g_n, c_n in codons_to_examine)
                line = 'codons_to_examine {0}\n'.format(locii_string)
                contents += line

            description_fh.write(contents)
            
            if max_read_length != None:
                description_fh.write('max_read_length {0}\n'.format(max_read_length))
            
            if markers_fn:
                description_fh.write('synthetic_fasta {0}\n'.format(markers_fn))

        bash_fh.write('echo {name}\n'.format(name=name))
        bash_fh.write('python ribosome_profiling_experiment.py --job_dir {0} launch --num_pieces {1}\n'.format(job_dir, num_pieces))

    return bash_fn

simulation_template = '''\
name {name}
work_prefix {home}/
scratch_prefix {home}/scratch/
relative_results_dir projects/ribosomes/experiments/simulation/{name}/results
template_description_fn {home}/projects/ribosomes/experiments/weinberg/RPF/job/description.txt
initiation_mean_numerator {initiation_mean_numerator}
CHX_mean {CHX_mean}
method {method}
perturbation_model {perturbation_model}
'''

TE_lines = '''\
RPF_description_fn {home}/projects/ribosomes/experiments/guydosh_cell/dom34KO_CHX/job/description.txt
mRNA_description_fn {home}/projects/ribosomes/experiments/guydosh_cell/dom34KO_mRNA-Seq/job/description.txt
'''

new_rates_line = '''\
new_rates_description_fn {home}/projects/ribosomes/experiments/belgium_2014_12_10/WT_1_FP/job/description.txt
'''

def make_noCHX_simulation_descriptions(num_pieces=12):
    initiation_rate_numerators = np.array([0, 1, 5, 10, 20, 30, 50]) * 15
    methods = ['mechanistic'] * len(initiation_rate_numerators)
    methods[0] = 'analytical'
    bash_fn = '{home}/projects/ribosomes/code/all_noCHX_simulation.sh'.format(home=os.environ['HOME'])
    with open(bash_fn, 'w') as bash_fh:
        for initiation_mean_numerator, method in zip(initiation_rate_numerators, methods):
            name = 'noCHX_{0}'.format(initiation_mean_numerator)
            
            job_dir = '{home}/projects/ribosomes/experiments/simulation/{0}/job'.format(name, home=os.environ['HOME'])
            if not os.path.isdir(job_dir):
                os.makedirs(job_dir)
            
            description_fn = '{0}/description.txt'.format(job_dir)
            with open(description_fn, 'w') as description_fh:
                contents = simulation_template.format(name=name,
                                                      initiation_mean_numerator=initiation_mean_numerator,
                                                      method=method,
                                                      CHX_mean=0,
                                                      home=os.environ['HOME'],
                                                     )
                contents += TE_lines.format(home=os.environ['HOME'])
                description_fh.write(contents)
            
            bash_fh.write('echo {name}\n'.format(name=name))
            bash_fh.write('python simulate.py --job_dir {0} launch --num_pieces {1}\n'.format(job_dir, num_pieces))

def make_CHX_simulation_descriptions(variable_TEs, num_pieces=12):
    CHX_means = np.array([0, 10, 50, 100, 200])
    initiation_rate_numerators = np.array([20, 50, 100, 200])

    bash_fn = '{home}/projects/ribosomes/code/all_CHX_simulation.sh'.format(home=os.environ['HOME'])
    with open(bash_fn, 'w') as bash_fh:
        for CHX_mean in CHX_means:
            for initiation_mean_numerator in initiation_rate_numerators:
                name = 'CHX_{0:0>3d}_{1:0>3d}_{2}'.format(CHX_mean,
                                                          initiation_mean_numerator,
                                                          'variable' if variable_TEs else 'fixed',
                                                         )
                
                job_dir = '{home}/projects/ribosomes/experiments/simulation/{0}/job'.format(name, home=os.environ['HOME'])
                if not os.path.isdir(job_dir):
                    os.makedirs(job_dir)
                    
                description_fn = '{0}/description.txt'.format(job_dir)
                with open(description_fn, 'w') as description_fh:
                    contents = simulation_template.format(name=name,
                                                          initiation_mean_numerator=initiation_mean_numerator,
                                                          CHX_mean=CHX_mean,
                                                          method='mechanistic',
                                                          home=os.environ['HOME'],
                                                         )
                    if variable_TEs:
                        contents += TE_lines.format(home=os.environ['HOME'])
                    description_fh.write(contents)

                bash_fh.write('echo {name}\n'.format(name=name))
                bash_fh.write('python simulate.py --job_dir {0} launch --num_pieces {1}\n'.format(job_dir, num_pieces))

def make_perturbation_model_descriptions(num_pieces=12):
    CHX_means = np.array([5, 15])
    initiation_rate_numerators = np.array([100])
    perturbation_models = ['same', 'uniform', 'reciprocal', 'shuffle', 'change_one', 'change_all']

    bash_fn = '{home}/projects/ribosomes/code/all_perturbation_model_simulation.sh'.format(home=os.environ['HOME'])
    with open(bash_fn, 'w') as bash_fh:
        for CHX_mean in CHX_means:
            for initiation_mean_numerator in initiation_rate_numerators:
                for perturbation_model in perturbation_models:
                    name = 'CHX_{0:0>3d}_{1:0>3d}_{2}'.format(CHX_mean,
                                                              initiation_mean_numerator,
                                                              perturbation_model,
                                                             )
                    
                    job_dir = '{home}/projects/ribosomes/experiments/simulation/{0}/job'.format(name, home=os.environ['HOME'])
                    if not os.path.isdir(job_dir):
                        os.makedirs(job_dir)
                        
                    description_fn = '{0}/description.txt'.format(job_dir)
                    with open(description_fn, 'w') as description_fh:
                        contents = simulation_template.format(name=name,
                                                              initiation_mean_numerator=initiation_mean_numerator,
                                                              CHX_mean=CHX_mean,
                                                              method='mechanistic',
                                                              perturbation_model=perturbation_model,
                                                              home=os.environ['HOME'],
                                                             )

                        if perturbation_model == 'change_all':
                            contents += new_rates_line.format(home=os.environ['HOME'])

                        description_fh.write(contents)

                    bash_fh.write('echo {name}\n'.format(name=name))
                    bash_fh.write('python simulate.py --job_dir {0} launch --num_pieces {1}\n'.format(job_dir, num_pieces))

arlen_locii = [('YLR075W', 98),
               ('YLR075W', 104),
               ('YLR075W', 106),
               ('YHR170W', 379),
              ]

markers_fns = {}
for name in ['bartel', 'stephanie', 'guydosh', 'pop', 'jan']:
    markers_fns[name] = '{0}/projects/ribosomes/data/{1}.fa'.format(os.environ['HOME'], name)

if __name__ == '__main__':
    kwargs = {}
    all_bash_fn = '{home}/projects/ribosomes/code/everything.sh'.format(home=os.environ['HOME'])
    fns = []

    fns.append(make_descriptions('artieri', 'polyA', **kwargs))
    fns.append(make_descriptions('artieri_gr_2', 'polyA', **kwargs))
    fns.append(make_descriptions('brar_science', 'polyA', **kwargs))
    fns.append(make_descriptions('belgium_2013_08_06', 'truseq', codons_to_examine=arlen_locii, **kwargs))
    fns.append(make_descriptions('belgium_2014_03_05', 'linker', **kwargs))

    fns.append(make_descriptions('belgium_2014_08_07', 'linker', markers_fn=markers_fns['stephanie'], **kwargs))
    fns.append(make_descriptions('belgium_2014_10_27', 'linker_local', markers_fn=markers_fns['stephanie'], **kwargs))
    fns.append(make_descriptions('belgium_2014_12_10', 'linker_local', codons_to_examine=arlen_locii, markers_fn=markers_fns['stephanie'], **kwargs))
    fns.append(make_descriptions('belgium_2015_03_16', 'linker_local', codons_to_examine=arlen_locii, markers_fn=markers_fns['stephanie'], **kwargs))
    fns.append(make_descriptions('dunn_elife', 'linker_local', **kwargs))

    def name_to_offset_type(name):
        if name in ['5min_rep1_foot', '30min_rep1_foot']:
            offset_type = 'gerashchenko_pnas'
        else:
            offset_type = 'yeast'
        return offset_type

    fns.append(make_descriptions('gerashchenko_pnas', 'polyA', max_read_length=44, name_to_offset_type=name_to_offset_type, **kwargs))
    
    fns.append(make_descriptions('gerashchenko_nar', 'nothing', max_read_length=50, **kwargs))
    
    fns.append(make_descriptions('guydosh_cell', 'linker_local', markers_fn=markers_fns['guydosh'], **kwargs))
    
    fns.append(make_descriptions('ingolia_science', 'polyA', **kwargs))
    fns.append(make_descriptions('lareau_elife', 'linker_local', **kwargs))
    fns.append(make_descriptions('mcmanus_gr', 'linker', **kwargs))
    
    fns.append(make_descriptions('weinberg', 'weinberg', markers_fn=markers_fns, **kwargs))
    
    fns.append(make_descriptions('zinshteyn_plos_genetics', 'polyA', **kwargs))

    def name_to_adapter_type(name):
        if name == 'WT_footprint':
            adapter_type = 'nothing'
        else:
            adapter_type = 'linker'
        return adapter_type
    fns.append(make_descriptions('pop_msb', 'linker', markers_fn=markers_fns['jan'], max_read_length=51, name_to_adapter_type=name_to_adapter_type, **kwargs))

    fns.append(make_descriptions('gardin_elife', 'truseq', **kwargs))
    fns.append(make_descriptions('nedialkova_cell', 'nothing', max_read_length=50, **kwargs))
    
    fns.append(make_descriptions('jan_science', 'linker', markers_fn=markers_fns['jan'], **kwargs))
    fns.append(make_descriptions('williams_science', 'linker', markers_fn=markers_fns['jan'], **kwargs))

    with open(all_bash_fn, 'w') as all_bash_fh:
        for fn in fns:
            all_bash_fh.write('echo {}\n'.format(fn))
            all_bash_fh.write('bash {}\n'.format(fn))
