import bisect
import Sequencing.genomes as genomes
import gtf
from collections import defaultdict

class IntervalTree(object):
    def __init__(self, intervals):
        sorted_intervals = sorted(intervals)
        
        # Pick a balancing center point by finding the middle of the middle
        # interval in the sorted list.
        middle_interval = sorted_intervals[len(sorted_intervals) / 2]
        self.center_point = (middle_interval.start + middle_interval.end) / 2

        to_left = []
        overlapping = []
        to_right = []
        for interval in sorted_intervals:
            if interval.end < self.center_point:
                to_left.append(interval)
            elif interval.start > self.center_point:
                to_right.append(interval)
            else:
                overlapping.append(interval)

        if to_left:
            self.left_subtree = IntervalTree(to_left)
        else:
            self.left_subtree = None

        if to_right:
            self.right_subtree = IntervalTree(to_right)
        else:
            self.right_subtree = None

        self.overlapping_by_start = sorted(overlapping, key=lambda interval: interval.start)
        self.overlapping_starts = [interval.start for interval in self.overlapping_by_start]
        self.overlapping_by_end = sorted(overlapping, key=lambda interval: interval.end)
        self.overlapping_ends = [interval.end for interval in self.overlapping_by_end]

    def containing_point(self, point):
        if point == self.center_point:
            overlapping = self.overlapping_by_start
        elif point < self.center_point:
            index = bisect.bisect(self.overlapping_starts, point)
            overlapping = self.overlapping_by_start[:index]
            if self.left_subtree:
                overlapping += self.left_subtree.containing_point(point)
        elif point > self.center_point:
            index = bisect.bisect_left(self.overlapping_ends, point)
            overlapping = self.overlapping_by_end[index:]
            if self.right_subtree:
                overlapping += self.right_subtree.containing_point(point)
        
        return overlapping

class IntervalEndpointList(object):
    def __init__(self, intervals):
        starts = [(interval.start, interval) for interval in intervals]
        ends = [(interval.end, interval) for interval in intervals]
        self.endpoint_list = sorted(starts + ends)
        self.endpoints = [endpoint for endpoint, interval in self.endpoint_list]

    def endpoint_contained_in(self, start, end):
        first_index = bisect.bisect_left(self.endpoints, start)
        past_last_index = bisect.bisect_right(self.endpoints, end)

        overlapping_endpoints = self.endpoint_list[first_index:past_last_index]
        # Note: intervals may contain duplicates.
        intervals = [interval for endpoint, interval in overlapping_endpoints]

        return intervals

class OverlapFinder(object):
    def __init__(self, intervals):
        self.interval_tree = IntervalTree(intervals)
        self.endpoint_list = IntervalEndpointList(intervals)

    def overlapping(self, start, end):
        from_tree = self.interval_tree.containing_point(start)
        from_list = self.endpoint_list.endpoint_contained_in(start, end)
        return sorted(set(from_tree + from_list))

class NamedOverlapFinder(object):
    def __init__(self, named_intervals, genome_dir):
        by_name = defaultdict(list)
        for interval in named_intervals:
            by_name[interval.seqname].append(interval)

        genome_index = genomes.get_genome_index(genome_dir)
        for name in by_name:
            start = gtf.Feature.sequence_edge(name, 0)
            end = gtf.Feature.sequence_edge(name, genome_index[name].length)
            by_name[name].append(start)
            by_name[name].append(end)

        self.overlap_finders = {}

        for name in by_name:
            self.overlap_finders[name] = OverlapFinder(by_name[name])

    def overlapping(self, seqname, start, end):
        return self.overlap_finders[seqname].overlapping(start, end)

    def find_closest_before(self, seqname, strand, start):
        window = 100
        before = []
        while before == []:
            overlapping = self.overlapping(seqname, start - window, start)
            before = [f for f in overlapping
                      if f.end <= start and f.strand in ['.', strand]
                     ]
            window *= 10

        before = sorted(before, key=lambda f: f.end, reverse=True)
        return before

    def find_closest_after(self, seqname, strand, end):
        window = 100
        after = []
        while after == []:
            overlapping = self.overlapping(seqname, end, end + window)
            after = [f for f in overlapping
                     if f.start >= end and f.strand in ['.', strand]
                    ]
            window *= 10

        after = sorted(after, key=lambda f: f.start)
        return after
