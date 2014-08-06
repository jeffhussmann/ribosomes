import three_t_fill_experiment
import os
from Sequencing import utilities, fastq, adapters
from Sequencing.Parallel import map_reduce
from collections import Counter

class WilkeningRNASeqExperiment(three_t_fill_experiment.ThreeTFillExperiment):
    num_stages = 2

    specific_results_files = []
    specific_figure_files = []
    specific_outputs = []
    specific_work = []
    specific_cleanup = []

    def __init__(self, **kwargs):
        super(WilkeningRNASeqExperiment, self).__init__(**kwargs)
        
        full_adapter_in_R2 = utilities.reverse_complement(self.barcode) + adapters.paired_end_R1_rc 
        self.adapter_in_R2 = full_adapter_in_R2[:19]

    def preprocess(self):
        total_reads = 0
        long_enough_reads = 0
        trimmed_lengths = Counter()
        read_pairs = self.get_read_pairs()

        truncated_in_R1 = self.adapter_in_R1[1:]
        truncated_in_R2 = self.adapter_in_R2[1:]
        
        with open(self.file_names['R1_preprocessed'], 'w') as R1_fh, \
             open(self.file_names['R2_preprocessed'], 'w') as R2_fh:

            for R1, R2 in read_pairs:
                total_reads += 1

                # Check for weird thing where expected overhang base doesn't
                # exist in primer dimers.
                R1_dimer_distance = adapters.adapter_hamming_distance(R1.seq,
                                                                      truncated_in_R1,
                                                                      len(R1.seq),
                                                                      len(truncated_in_R1),
                                                                      len(self.barcode),
                                                                     )
                R2_dimer_distance = adapters.adapter_hamming_distance(R2.seq,
                                                                      truncated_in_R2,
                                                                      len(R2.seq),
                                                                      len(truncated_in_R2),
                                                                      len(self.barcode),
                                                                     )
                if R1_dimer_distance <= 3 and R2_dimer_distance <= 3:
                    position = len(self.barcode)
                else:
                    position = adapters.consistent_paired_position(R1.seq,
                                                                   R2.seq,
                                                                   self.adapter_in_R1,
                                                                   self.adapter_in_R2,
                                                                   19,
                                                                   3,
                                                                  )
                if position != None:
                    trimmed_lengths[position] += 1
                    if position - len(self.barcode) < 12:
                        continue
                else:
                    position = len(R1.seq)

                long_enough_reads += 1

                payload_slice = slice(len(self.barcode), position)

                processed_R1_record = fastq.make_record(R1.name,
                                                        R1.seq[payload_slice],
                                                        R1.qual[payload_slice],
                                                       )
                processed_R2_record = fastq.make_record(R2.name,
                                                        R2.seq[payload_slice],
                                                        R2.qual[payload_slice],
                                                       )
                R1_fh.write(processed_R1_record)
                R2_fh.write(processed_R2_record)

        trimmed_lengths = utilities.counts_to_array(trimmed_lengths)
        self.write_file('trimmed_lengths', trimmed_lengths)
        self.log.extend(
            [('Total read pairs', total_reads),
             ('Long enough', long_enough_reads),
            ]
        )

if __name__ == '__main__':
    script_path = os.path.realpath(__file__)
    map_reduce.controller(WilkeningRNASeqExperiment, script_path)
