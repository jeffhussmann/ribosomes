''' Utilities for counting reads mapped to each position in a gene. '''

import numpy as np

class PositionCounts(object):
    ''' Wrapper around an array of counts of positions for an extent and for a
        buffer on either edge that allows for indexing relative to the extent's
        start and end.
    '''
    def __init__(self, extent_length, left_buffer, right_buffer, counts=None):
        self.left_buffer = left_buffer
        self.right_buffer = right_buffer
        self.extent_length = extent_length
        if counts == None:
            self.counts = np.zeros(left_buffer + right_buffer + extent_length, int)
        else:
            assert len(counts) == left_buffer + right_buffer + extent_length
            self.counts = counts

    @classmethod
    def from_string(cls, left_buffer, right_buffer, string):
        counts = np.array(map(int, string.strip().split()))
        extent_length = len(counts) - left_buffer - right_buffer
        return cls(extent_length, left_buffer, right_buffer, counts=counts)

    def __str__(self):
        string = ' '.join(map(str, self.counts)) + '\n'
        return string

    def adjust_relative_to_start(self, key):
        if isinstance(key, int):
            adjusted_key = self.left_buffer + key
            if adjusted_key < 0:
                raise IndexError
            else:
                return adjusted_key
        elif isinstance(key, slice):
            if key.start == None:
                adjusted_start = self.left_buffer
            else:
                adjusted_start = self.left_buffer + key.start
            if adjusted_start < 0:
                raise IndexError

            if key.stop == None:
                adjusted_stop = self.left_buffer + self.extent_length
            else:
                adjusted_stop = self.left_buffer + key.stop
            if adjusted_stop < 0:
                raise IndexError
            
            if key.step == None:
                adjusted_step = 1
            elif key.step < 0:
                raise ValueError('Negative step not allowed')
            else:
                adjusted_step = key.step

            return slice(adjusted_start, adjusted_stop, adjusted_step)

    def __getitem__(self, key):
        adjusted_key = self.adjust_relative_to_start(key)
        return self.counts[adjusted_key]

    def __setitem__(self, key, value):
        adjusted_key = self.adjust_relative_to_start(key)
        self.counts[adjusted_key] = value

    @property
    def relative_to_end(self):
        return PositionCounts.RelativeToEndCounts(self)

    class RelativeToEndCounts(object):
        ''' Convoluted hack to allow PositionCounts.relative_to_end to be an object
            that has __getitem__ and __setitem__ methods so that
            position_counts.relative_to_end[a:b:c] works.
        '''
        def __init__(self, position_counts):
            self.right_buffer = position_counts.right_buffer
            self.extent_length = position_counts.extent_length
            self.counts = position_counts.counts

        def adjust_relative_to_end(self, key):
            if isinstance(key, int):
                adjusted_key = len(self.counts) - self.right_buffer - key
                if adjusted_key < 0:
                    raise IndexError
                else:
                    return adjusted_key
            elif isinstance(key, slice):
                if key.start == None:
                    adjusted_start = len(self.counts) - self.right_buffer
                else:
                    adjusted_start = len(self.counts) - self.right_buffer - key.start
                if adjusted_start < 0:
                    raise IndexError

                if key.stop == None:
                    adjusted_stop = len(self.counts) - self.right_buffer - self.extent_length
                else:
                    adjusted_stop = len(self.counts) - self.right_buffer - key.stop
                if adjusted_stop < 0:
                    raise IndexError
                
                if key.step == None:
                    adjusted_step = -1
                elif key.step < 0:
                    raise ValueError('Negative step not allowed')
                else:
                    adjusted_step = -key.step

                return slice(adjusted_start, adjusted_stop, adjusted_step)

        def __getitem__(self, key):
            adjusted_key = self.adjust_relative_to_end(key)
            return self.counts[adjusted_key]

        def __setitem__(self, key, value):
            adjusted_key = self.adjust_relative_to_end(key)
            self.counts[adjusted_key] = value