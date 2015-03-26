import heapq
import numpy as np
import logging
import codons
import positions
import visualize
import pausing
import os
from collections import Counter, defaultdict
import Sequencing.Parallel
import ribosome_profiling_experiment
import Serialize.read_positions as read_positions

experiment_from_fn = ribosome_profiling_experiment.RibosomeProfilingExperiment.from_description_file_name

exponential = np.random.exponential

class Message(object):
    def __init__(self, codon_sequence, initiation_mean, codon_means, CHX_mean):
        self.codon_mean_sequence = [codon_means[codon_id] for codon_id in codon_sequence]
        self.initiation_mean = initiation_mean
        self.events = []
        self.left_edges = set()
        self.leftmost_ribosome = None
        self.ribosomes = {}
        self.current_id_number = 0
        self.current_event_number = 0
        self.current_time = 0
        self.first_runoff_event_number = None

        self.CHX_introduction_time = np.inf
        self.CHX_mean = CHX_mean
        
        self.initiate(0)

    def initiate(self, time):
        if self.leftmost_ribosome and self.leftmost_ribosome.position - 5 <= 4:
            # Occluded from starting
            if self.leftmost_ribosome.arrested:
                # The occlusion will never clear
                pass
            else:
                self.register_next_initiation(time)
        else:
            ribosome = Ribosome(self)
            self.leftmost_ribosome = ribosome
            ribosome.register_next_advance_time(time)
            if time > self.CHX_introduction_time:
                ribosome.register_CHX_arrival_time(time)

            self.register_next_initiation(time)
    
    def register_next_initiation(self, time):
        next_initiation_time = time + exponential(self.initiation_mean)
        heapq.heappush(self.events, (next_initiation_time, 'initiate', None)) 

    def process_next_event(self):
        if not self.events:
            return 'empty'
        else:
            time, event, ribosome = heapq.heappop(self.events)
            self.current_event_number += 1
            self.current_time = time
            
            if event == 'steady_state':
                pass
            elif event == 'initiate':
                self.initiate(time)
            elif event == 'advance':
                was_runoff = ribosome.advance(time)
            
                if was_runoff and self.first_runoff_event_number == None:
                    self.first_runoff_event_number = self.current_event_number
            elif event == 'CHX_arrival':
                ribosome.arrested = True
                ribosome.arrested_at = time

            return event

    def evolve_to_steady_state(self):
        while self.first_runoff_event_number == None:
            event = self.process_next_event()
        
        steady_state_time = np.random.uniform(self.current_time, 2 * self.current_time)
        heapq.heappush(self.events, (steady_state_time, 'steady_state', None))

        while event != 'steady_state':
            event = self.process_next_event()

        if self.CHX_mean != 0:
            self.introduce_CHX(self.current_time)

        while event != 'empty':
            event = self.process_next_event()

        return event

    def introduce_CHX(self, introduction_time):
        self.CHX_introduction_time = introduction_time

        for ribosome in self.ribosomes.values():
            ribosome.register_CHX_arrival_time(introduction_time)

    def collect_measurements(self):
        return Counter(r.position for r in self.ribosomes.itervalues())
    
    def __str__(self):
        description = 'Ribosomes:\n'
        for i in self.ribosomes:
            description += '\t' + str(self.ribosomes[i]) + '\n'

        description += 'Events:\n'
        for time, event, ribosome in sorted(self.events):
            if ribosome == None:
                id_number = None
            else:
                id_number = ribosome.id_number
            event_string = str((time, event, id_number))
            description += '\t' + event_string + '\n'

        return description

class Ribosome(object):
    def __init__(self, message):
        self.message = message
        self.id_number = message.current_id_number
        message.ribosomes[self.id_number] = self
        message.current_id_number += 1
        message.left_edges.add(-5)
        self.position = 0
        self.arrested = False
        self.arrested_at = np.inf

    def __str__(self):
        description = 'id_number: {0}, position: {1}, arrested: {2} ({3})'.format(self.id_number,
                                                                                  self.position,
                                                                                  self.arrested,
                                                                                  self.arrested_at,
                                                                                 )
        return description

    def advance(self, time):
        was_runoff = False

        if self.arrested:
            pass
        elif self.position + 5 in self.message.left_edges:
            # Occluded from advancing
            self.register_next_advance_time(time)
        else:
            self.message.left_edges.remove(self.position - 5)
            
            if self.position == len(self.message.codon_mean_sequence) - 1:
                self.message.ribosomes.pop(self.id_number)
                was_runoff = True
            else:
                self.position += 1
                self.message.left_edges.add(self.position - 5)

                self.register_next_advance_time(time)
            
        return was_runoff

    def register_next_advance_time(self, time):
        mean = self.message.codon_mean_sequence[self.position]
        next_advance_time = time + exponential(mean)
        heapq.heappush(self.message.events, (next_advance_time, 'advance', self))

    def register_CHX_arrival_time(self, time):
        CHX_arrival_time = time + exponential(self.message.CHX_mean)
        heapq.heappush(self.message.events, (CHX_arrival_time, 'CHX_arrival', self)) 

class SimulationExperiment(Sequencing.Parallel.map_reduce.MapReduceExperiment):
    num_stages = 1

    specific_results_files = [
        ('simulated_codon_counts', read_positions, '{name}_simulated_codon_counts.hdf5'),
        ('stratified_mean_enrichments', 'pickle', '{name}_stratified_mean_enrichments.pkl'),
        ('mean_densities', read_positions, '{name}_mean_densities.hdf5'),
    ]

    specific_figure_files = [
        ('mean_densities', '{name}_mean_densities.pdf'),
    ]

    specific_outputs = [
        ['simulated_codon_counts',
        ],
    ]

    specific_work = [
        ['produce_counts',
        ],
    ]

    specific_cleanup = [
        ['compute_stratified_mean_enrichments',
         'compute_mean_densities',
         'plot_mean_densities',
        ]
    ]

    def __init__(self, **kwargs):
        super(SimulationExperiment, self).__init__(**kwargs)

        self.template_experiment = experiment_from_fn(kwargs['template_description_fn'])
        
        if 'RPF_description_fn' in kwargs:
            self.RPF_experiment = experiment_from_fn(kwargs['RPF_description_fn'])
        else:
            self.RPF_experiment = None

        if 'mRNA_description_fn' in kwargs:
            self.mRNA_experiment = experiment_from_fn(kwargs['mRNA_description_fn'])
        else:
            self.mRNA_experiment = None

        self.initiation_mean_numerator = int(kwargs['initiation_mean_numerator'])
        self.CHX_mean = int(kwargs['CHX_mean'])

        self.method = kwargs['method']

    def load_TEs(self):
        if self.RPF_experiment and self.mRNA_experiment:
            TEs = pausing.load_TEs(self.RPF_experiment, self.mRNA_experiment)
        else:
            TEs = defaultdict(lambda: 1)

        return TEs

    def produce_counts(self):
        if self.method == 'mechanistic':
            self.simulate()
        elif self.method == 'analytical':
            self.distribute_analytically()

    def simulate(self):
        buffered_codon_counts = self.template_experiment.read_file('buffered_codon_counts')
        
        stratified_mean_enrichments = self.template_experiment.read_file('stratified_mean_enrichments')
        codon_means = stratified_mean_enrichments[0, 1, 2]
        for codon in codons.stop_codons:
            codon_means[codon] = 1
        
        TEs = self.load_TEs()
        initiation_means = {gene_name: self.initiation_mean_numerator / TEs[gene_name] for gene_name in buffered_codon_counts} 

        all_gene_names = sorted(buffered_codon_counts)
        piece_gene_names = Sequencing.Parallel.piece_of_list(all_gene_names,
                                                             self.num_pieces,
                                                             self.which_piece,
                                                            )
        
        simulated_codon_counts = {}
        cds_slice = slice('start_codon', ('stop_codon', 1))
        for i, gene_name in enumerate(piece_gene_names):
            logging.info('Starting {0} ({1:,} / {2:,})'.format(gene_name, i, len(piece_gene_names) - 1))
            identities = buffered_codon_counts[gene_name]['identities']
            codon_sequence = identities[cds_slice]

            real_counts = buffered_codon_counts[gene_name]['relaxed'][cds_slice]
            total_real_counts = sum(real_counts)

            all_measurements = Counter()
            num_messages = 0
            while sum(all_measurements.values()) < total_real_counts:
                message = Message(codon_sequence, initiation_means[gene_name], codon_means, self.CHX_mean)
                message.evolve_to_steady_state()
                all_measurements.update(message.collect_measurements())
                num_messages += 1

            simulated_counts = positions.PositionCounts(identities.landmarks,
                                                        identities.left_buffer,
                                                        identities.right_buffer,
                                                       )

            for key, value in all_measurements.items():
                simulated_counts['start_codon', key] = value
            
            simulated_codon_counts[gene_name] = {'identities': identities,
                                                 'relaxed': simulated_counts,
                                                }
            logging.info('{0:,} counts generated for {1} from {2:,} messages'.format(sum(all_measurements.values()), gene_name, num_messages))

        self.write_file('simulated_codon_counts', simulated_codon_counts)
    
    def load_codon_means(self):
        stratified_mean_enrichments = self.template_experiment.read_file('stratified_mean_enrichments')
        codon_means = stratified_mean_enrichments[0, 1, 2]
        for codon in codons.stop_codons:
            codon_means[codon] = 1

        return codon_means
    
    def distribute_analytically(self):
        buffered_codon_counts = self.template_experiment.read_file('buffered_codon_counts')
        
        all_gene_names = sorted(buffered_codon_counts)
        piece_gene_names = Sequencing.Parallel.piece_of_list(all_gene_names,
                                                             self.num_pieces,
                                                             self.which_piece,
                                                            )
        
        simulated_codon_counts = {}
        cds_slice = slice('start_codon', ('stop_codon', 1))
        for i, gene_name in enumerate(piece_gene_names):
            identities = buffered_codon_counts[gene_name]['identities']
            codon_sequence = identities[cds_slice]

            real_counts = buffered_codon_counts[gene_name]['relaxed'][cds_slice]
            total_real_counts = sum(real_counts)

            rates_array = np.array([codon_rates[codon_id] for codon_id in codon_sequence])
            fractions_array = rates_array / sum(rates_array)
            
            simulated_counts = positions.PositionCounts(identities.landmarks,
                                                         identities.left_buffer,
                                                         identities.right_buffer,
                                                        )

            for position, fraction in enumerate(fractions_array):
                simulated_counts['start_codon', position] = np.random.binomial(total_real_counts, fraction)
            
            simulated_codon_counts[gene_name] = {'identities': identities,
                                                 'relaxed': simulated_counts,
                                                }

        self.write_file('simulated_codon_counts', simulated_codon_counts)
    
    def compute_stratified_mean_enrichments(self):
        allowed_at_pause = set(codons.non_stop_codons)
        not_allowed_at_stall = {}

        codon_counts = self.read_file('simulated_codon_counts',
                                      specific_keys={'relaxed', 'identities'},
                                     )
        gene_names, _, _ = pausing.get_highly_expressed_gene_names({'self': codon_counts}, min_mean=0)

        around_lists = pausing.metacodon_around_pauses(codon_counts,
                                                       allowed_at_pause,
                                                       not_allowed_at_stall,
                                                       gene_names,
                                                      )
        stratified_mean_enrichments = pausing.compute_stratified_mean_enrichments(around_lists)
        self.write_file('stratified_mean_enrichments', stratified_mean_enrichments)
    
    def compute_mean_densities(self):
        codon_counts = self.read_file('simulated_codon_counts')
        mean_densities = positions.compute_averaged_codon_densities(codon_counts)
        self.write_file('mean_densities', mean_densities)

    def plot_mean_densities(self):
        visualize.plot_averaged_codon_densities([(self.name, self.read_file('mean_densities'), 0)],
                                                self.figure_file_names['mean_densities'],
                                                past_edge=10,
                                                plot_up_to=1000,
                                                smooth=False,
                                               )

if __name__ == '__main__':
    script_path = os.path.realpath(__file__)
    Sequencing.Parallel.map_reduce.controller(SimulationExperiment, script_path)
