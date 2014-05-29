import os
import glob
import ribosome_profiling_experiment

def build_all_experiments():
    families = ['zinshteyn_plos_genetics',
                'ingolia_science',
                'weinberg',
                'dunn_elife',
                'gerashchenko_pnas',
                'guydosh_cell',
                #'mcmanus_gr',
                #'artieri',
               ]
    experiments = {}
    for family in families:
        print family
        experiments[family] = {}
        prefix = '/home/jah/projects/arlen/experiments/{0}/'.format(family)
        dirs = [path for path in glob.glob('{}*'.format(prefix)) if os.path.isdir(path)]
        for d in sorted(dirs):
            _, name = os.path.split(d)
            print '\t', name
            description_file_name = '{0}/job/description.txt'.format(d)
            experiments[family][name] = ribosome_profiling_experiment.build_from_description_file_name(description_file_name)

    return experiments

def read_counts_and_RPKMS():
    experiments = build_all_experiments()
    for family in experiments:
        for name in experiments[family]:
            experiments[family][name].compute_total_read_counts()
            experiments[family][name].compute_RPKMs()
