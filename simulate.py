import heapq
import numpy as np
import codons
import positions
from collections import Counter

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

def produce_measurements(codon_counts, codon_rates, initiation_rate, num_copies):
    identities = codon_counts['identities']
    codon_sequence = identities['start_codon':('stop_codon', 1)]

    all_measurements = Counter()
    for i in xrange(num_copies):
        message = Message(codon_sequence, initiation_rate, codon_rates)
        message.evolve_to_steady_state()
        all_measurements.update(message.collect_measurements())

    codon_counts['relaxed'] = positions.PositionCounts(identities.landmarks,
                                                         identities.left_buffer,
                                                         identities.right_buffer,
                                                        )
    for key, value in all_measurements.items():
        codon_counts['relaxed']['start_codon', key] = value

    return codon_counts
