import heapq
import numpy as np
import logging
import codons
import positions
import os
from collections import Counter
import Sequencing.Parallel
import ribosome_profiling_experiment
import Serialize.read_positions as read_positions

class Message(object):
    def __init__(self, codon_sequence, initiation_rate, codon_rates):
        self.codon_sequence = codon_sequence
        self.codon_rates = codon_rates
        self.initiation_rate = initiation_rate
        self.events = []
        self.left_edges = set()
        self.leftmost_ribosome = None
        self.ribosomes = {}
        self.current_id_number = 0
        self.current_event_number = 0
        self.first_runoff_event_number = None
        
        self.initiate(0)

    def initiate(self, time):
        if self.leftmost_ribosome and self.leftmost_ribosome.position - 5 <= 4:
            # Occluded from starting
            pass
        else:
            r = Ribosome(self)
            self.leftmost_ribosome = r
            r.register_advance_time(time)
        
        next_time = time + np.random.exponential(self.initiation_rate)
        heapq.heappush(self.events, (next_time, 'initiate'))

    def process_next_event(self):
        next_time, next_event = heapq.heappop(self.events)
        if next_event == 'measure':
            return 'measure'
        elif next_event == 'initiate':
            self.initiate(next_time)
        else:
            ribosome = next_event
            was_runoff = ribosome.advance(next_time)
        
            if was_runoff and self.first_runoff_event_number == None:
                self.first_runoff_event_number = self.current_event_number
        
        self.current_event_number += 1

    def __str__(self):
        description = 'Ribosomes:\n'
        for i in self.ribosomes:
            description += '\t' + str(self.ribosomes[i]) + '\n'

        description += 'Events:\n'
        for event in sorted(self.events):
            description += '\t' + str(event) + '\n'

        return description

    def evolve_to_steady_state(self):
        while self.first_runoff_event_number == None:
            self.process_next_event()
        
        additional_events = np.random.randint(0, self.first_runoff_event_number)

        for i in xrange(additional_events):
            self.process_next_event()

    def collect_measurements(self):
        return Counter(r.position for r in self.ribosomes.itervalues())

class Ribosome(object):
    def __init__(self, message):
        self.message = message
        self.id_number = message.current_id_number
        message.ribosomes[self.id_number] = self
        message.current_id_number += 1
        message.left_edges.add(-5)
        self.position = 0

    def __str__(self):
        return 'id_number: {0}, position: {1}'.format(self.id_number, self.position)

    def advance(self, time):
        was_runoff = False

        if self.position + 5 in self.message.left_edges:
            # Occluded from advancing
            pass
        else:
            self.message.left_edges.remove(self.position - 5)
            if self.position == len(self.message.codon_sequence) - 1:
                self.message.ribosomes.pop(self.id_number)
                was_runoff = True
            else:
                self.position += 1
                self.message.left_edges.add(self.position - 5)
            
        if not was_runoff:
            self.register_advance_time(time)

        return was_runoff

    def register_advance_time(self, time):
        codon_id = self.message.codon_sequence[self.position]
        rate = self.message.codon_rates[codon_id]
        next_time = time + np.random.exponential(rate)
        heapq.heappush(self.message.events, (next_time, self))

class SimulationExperiment(Sequencing.Parallel.map_reduce.MapReduceExperiment):
    num_stages = 1

    specific_results_files = [
        ('simulated_codon_counts', read_positions, '{name}_simulated_codon_counts.hdf5'),
        ('stratified_mean_enrichments', 'pickle', '{name}_stratified_mean_enrichments.pkl'),
    ]

    specific_figure_files = []

    specific_outputs = [
        ['simulated_codon_counts',
        ],
    ]

    specific_work = [
        ['simulate'],
    ]

    specific_cleanup = [
        ['compute_stratified_mean_enrichments',
        ]
    ]

    def __init__(self, **kwargs):
        super(SimulationExperiment, self).__init__(**kwargs)

        self.template_description_fn = kwargs['template_description_fn']
        self.template_experiment = ribosome_profiling_experiment.RibosomeProfilingExperiment.from_description_file_name(self.template_description_fn)

    def simulate(self, **kwargs):
        buffered_codon_counts = self.template_experiment.read_file('buffered_codon_counts')
        
        stratified_mean_enrichments = self.template_experiment.read_file('stratified_mean_enrichments')
        codon_rates = stratified_mean_enrichments[0, 1, 2]
        for codon in codons.stop_codons:
            codon_rates[codon] = 1
        
        initiation_rates = {gene_name: 50 for gene_name in buffered_codon_counts} 

        all_gene_names = sorted(buffered_codon_counts)
        piece_gene_names = Sequencing.Parallel.piece_of_list(all_gene_names,
                                                             self.num_pieces,
                                                             self.which_piece,
                                                            )
        
        simulated_codon_counts = {}
        cds_slice = slice('start_codon', ('stop_codon', 1))
        for i, gene_name in enumerate(piece_gene_names):
            logging.info('Starting {0} ({1:,} / {2:,})'.format(gene_name, i, len(piece_gene_names)))
            identities = buffered_codon_counts[gene_name]['identities']
            codon_sequence = identities[cds_slice]

            real_counts = buffered_codon_counts[gene_name]['relaxed'][cds_slice]
            total_real_counts = sum(real_counts)

            all_measurements = Counter()
            while sum(all_measurements.values()) < total_real_counts:
                message = Message(codon_sequence, initiation_rates[gene_name], codon_rates)
                message.evolve_to_steady_state()
                all_measurements.update(message.collect_measurements())

            simulated_counts = positions.PositionCounts(identities.landmarks,
                                                        identities.left_buffer,
                                                        identities.right_buffer,
                                                       )

            for key, value in all_measurements.items():
                simulated_counts['start_codon', key] = value
            
            simulated_codon_counts[gene_name] = {'identities': identities,
                                                 'relaxed': simulated_counts,
                                                }
            logging.info('{0:,} counts generated for {1}'.format(sum(all_measurements.values()), gene_name))

        self.write_file('simulated_codon_counts', simulated_codon_counts)
    
    def compute_stratified_mean_enrichments(self):
        allowed_at_pause = set(codons.non_stop_codons)
        not_allowed_at_stall = set()

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

if __name__ == '__main__':
    script_path = os.path.realpath(__file__)
    Sequencing.Parallel.map_reduce.controller(SimulationExperiment, script_path)
