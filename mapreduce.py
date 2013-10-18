class MapReduceExperiment(object):
    def __init__(self, **kwargs):
        self.name = kwargs['name']
        self.work_prefix = kwargs['work_prefix'].rstrip('/')
        self.scratch_prefix = kwargs['scratch_prefix'].rstrip('/')
        self.relative_results_dir = kwargs['relative_results_dir'].strip('/')
        self.num_pieces = kwargs['num_pieces']
        self.which_piece = kwargs['which_piece']
        
        suffix = Parallel.split_file.generate_suffix(self.num_pieces, self.which_piece)

        self.scratch_results_dir = '{0}/{1}{2}'.format(self.scratch_prefix,
                                                       self.relative_results_dir,
                                                       suffix,
                                                      )
        self.work_results_dir = '{0}/{1}'.format(self.work_prefix,
                                                 self.relative_results_dir,
                                                )
        
        if self.which_piece == -1:
            # Sentinel value that indicates this is the merged experiment.
            # Only create the directory in this instance to avoid race
            # conditions.
            if not os.path.isdir(self.work_results_dir):
                os.makedirs(self.work_results_dir)
        else:
            # Only need to create the scratch directory if this ISN't the merged
            # experiment.
            if not os.path.isdir(self.scratch_results_dir):
                os.makedirs(self.scratch_results_dir)
        
        num_stages = len(self.outputs)
        timing_names = [('timing_{0}'.format(i), 'log', '{{name}}_timing_{0}.txt'.format(i))
                        for i in range(num_stages)]
        self.results_files.extend(timing_names)
        
        self.file_names = {}
        self.merged_file_names = {}
        self.file_types = {}
        for key, serialize_type, tail_template in self.results_files:
            file_tail = tail_template.format(name=self.name)
            self.file_names[key] = '{0}/{1}'.format(self.scratch_results_dir, file_tail)
            self.merged_file_names[key] = '{0}/{1}'.format(self.work_results_dir, file_tail)
            self.file_types[key] = serialize_type

        if self.which_piece == -1:
            self.file_names = self.merged_file_names

    def write_file(self, key, data):
        Serialize.write_file(data, self.file_names[key], self.file_types[key])

    def read_file(self, key):
        data = Serialize.read_file(self.file_names[key], self.file_types[key])
        return data

    def do_work(self, stage):
        times = []
        for function, description in self.work[stage]:
            start_time = time.time()
            function()
            end_time = time.time()
            times.append((description, end_time - start_time))

        self.write_file('timing_{0}'.format(stage), times)

    def do_cleanup(self, stage):
        for function in self.cleanup[stage]:
            function()

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--job_dir',
                        required=True,
                       )

    subparsers = parser.add_subparsers(dest='subparser_name')

    parser_process = subparsers.add_parser('process')
    parser_process.add_argument('--num_pieces', type=int,
                                help='how many pieces',
                                default=1,
                               )
    parser_process.add_argument('--which_piece', type=int,
                                help='which piece this is',
                                default=0,
                               )
    parser_process.add_argument('--stage', type=int,
                                help='stage',
                                default=0,
                               )
    
    parser_launch = subparsers.add_parser('launch')
    parser_launch.add_argument('--num_pieces', type=int,
                               default=1,
                              )

    parser_finish = subparsers.add_parser('finish')
    parser_finish.add_argument('--num_pieces', type=int,
                                help='how many pieces',
                                default=1,
                               )
    parser_finish.add_argument('--stage', type=int,
                                help='stage',
                                default=0,
                               )

    args = parser.parse_args()
    return args

def parse_experiment_description(description_fn):
    description = dict(line.strip().split() for line in open(description_fn))
    return description

def launch(args):
    description_file_name = '{0}/description.txt'.format(args.job_dir)
    description = parse_experiment_description(description_file_name)

    script_path = os.path.realpath(__file__)

    def make_process_command(args, which_piece, stage):
        command = ['python', script_path,
                   '--job_dir', args.job_dir,
                   'process',
                   '--num_pieces', str(args.num_pieces),
                   '--which_piece', str(which_piece),
                   '--stage', str(stage),
                  ]
        command_string = ' '.join(command) + '\n'
        return command_string
    
    def make_finish_command(args, stage):
        command = ['python', script_path,
                   '--job_dir', args.job_dir,
                   'finish',
                   '--num_pieces', str(args.num_pieces),
                   '--stage', str(stage),
                  ]
        command_string = ' '.join(command) + '\n'
        return command_string

    for stage in range(1):
        process_file_name = '{0}/process_{1}_stage_{2}'.format(args.job_dir,
                                                               args.num_pieces,
                                                               stage,
                                                              )
        finish_file_name = '{0}/finish_{1}_stage_{2}'.format(args.job_dir,
                                                             args.num_pieces,
                                                             stage,
                                                            )
                              
        with open(process_file_name, 'w') as process_file:
            for which_piece in range(args.num_pieces):
                line = make_process_command(args, which_piece, stage)
                process_file.write(line)

        with open(finish_file_name, 'w') as finish_file:
            line = make_finish_command(args, stage)
            finish_file.write(line)
    
        print 'Launched stage {0} with parallel'.format(stage)
        subprocess.check_call('parallel < {0}'.format(process_file_name), shell=True)
        subprocess.check_call('bash {0}'.format(finish_file_name), shell=True)

def process(args):
    description_file_name = '{0}/description.txt'.format(args.job_dir)
    description = parse_experiment_description(description_file_name)

    experiment = ExperimentClass(num_pieces=args.num_pieces,
                                             which_piece=args.which_piece,
                                             **description)
    experiment.do_work(args.stage)

def finish(args):
    description_file_name = '{0}/description.txt'.format(args.job_dir)
    description = parse_experiment_description(description_file_name)
    
    merged = ExperimentClass(num_pieces=args.num_pieces,
                                         which_piece=-1,
                                         **description)
    pieces = [ExperimentClass(num_pieces=args.num_pieces,
                                          which_piece=which_piece,
                                          **description)
              for which_piece in range(args.num_pieces)]

    for key in merged.outputs[args.stage]:
        piece_file_names = [piece.file_names[key] for piece in pieces]
        merged_file_name = merged.merged_file_names[key]
        Serialize.merge_files(piece_file_names,
                             merged_file_name,
                             merged.file_types[key],
                            )

    merged.do_cleanup(args.stage)

def controller(ExperimentClass, script_path<Mouse>C!!1<Mouse>C!!4):
    args = parse_arguments()
    if args.subparser_name == 'launch':
        launch(args)
    elif args.subparser_name == 'process':
        process(args)
    elif args.subparser_name == 'finish':
        finish(args)
