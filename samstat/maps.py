""" Region Map
This file contains the code for reading GFF3 files into region maps
"""
import os
import warnings
import pysam
from collections import namedtuple
from functools import lru_cache

class AlignmentMap(dict):
    SamIn = namedtuple('InLine',
                       ['alignment_number',
                        'cigar',
                        'rname_positions'])

    def __init__(self, path):
        self.update(self.read_alignment_map(path))

    @classmethod
    def read_alignment_map(cls, path):
        """Reads Alignment map SAM/BAM file into dictionary"""
        amap = {}
        samfile = pysam.AlignmentFile(path, 'r')
        for count, seq_line in enumerate(samfile):
            amap = seq_line.query_name
            qnames.setdefault(qname, cls.SamIn([0], seq_line.cigar, []))
            qnames[qname].alignment_number[0] += 1
            try:
                amap[qname].rname_positions.append((seq_line.reference_name,
                                                      seq_line.reference_start))
            except ValueError:
                warnings.warn('Reference Name is -1, Line #: {}'.format(count))

        return amap

def split_gen(s, delims):
    """iterates a delimited line"""
    start = 0
    for i, c in enumerate(s):
        if c in delims:
            yield s[start:i]
            start = i+1
    yield s[start:]

class RegionMap(object):
    """Reads creates a feature location map from a gff file that can be
    used to determine gene attribute types from sequence location ranges

    Region Map Structure:
        {'RNAME': [Features: [location, strand],], Length}
    """
    Region = namedtuple('Region', ['features', 'length'])
    Feature = namedtuple('Feature', ['location', 'strand'])
    def __init__(self, gff_path, accepted_features='exon'):
        if isinstance(accepted_features, str):
            accepted_features = tuple([accepted_features])
        self.feature_types = accepted_features
        print(accepted_features)
        self.rmap = self.read_gff(gff_path, feature_types=self.feature_types)

    @classmethod
    def read_gff(cls, gff_path, feature_types):
        """Reads a gff file into a region map hash"""
        region_map = {}
        feature_temp = {key: [] for key in feature_types}
        with open(gff_path, 'r') as gff:
            for count, line in enumerate(gff):
                try:
                    if line.startswith('#'):
                        continue
                    line_gen = split_gen(line, '\t ')
                    region = next(line_gen)
                    next(line_gen)
                    feature = next(line_gen)
                    location = [int(next(line_gen)), int(next(line_gen))]
                    next(line_gen)
                    strand = next(line_gen)
                    tmp_feature = cls.Feature(location, strand)
                    try:
                        region_map[region].features[feature].append(tmp_feature)
                    except KeyError:
                        if feature == 'region':
                            region_map.setdefault(region,
                                                  cls.Region(feature_temp,
                                                             location[1]))
                        else:
                            region_map.setdefault(region,
                                                  cls.Region(feature_temp,
                                                             None))
                            try:
                                region_map[region].features[feature]\
                                                  .append(tmp_feature)
                            except KeyError:
                                pass

                except StopIteration:
                    warnings.warn('Invalid line: {} ... skipped'.format(count))
        return cls.calc_missing_region_lengths(region_map)

    @classmethod
    def calc_missing_region_lengths(cls, region_map):
        """Calculates region lengths for regions without a region line"""
        for key, region in region_map.items():
            if region.length is None:
                features = iter(region.features)
                largest = next(features)
                for feat in features:
                    if largest < feat.location[1]: largest = feat.location[1]
                region_map[key] = cls.Region(region.features, largest)
        return region_map

    @lru_cache(maxsize=None)
    def get_location_clasification(self,
                                   region_name,
                                   location_start,
                                   location_stop):
        """Gets location classification from region_map"""
        for key, features in self.rmap[region_name].features.items():
            for feat_range, _ in features:
                start_flag = location_start in range(feat_range[0],
                                                     feat_range[1])
                stop_flag = location_stop in range(feat_range[0],
                                                   feat_range[1])

                if stop_flag is True and start_flag is True:
                    return key
                elif stop_flag != start_flag:
                    return 'combo'
                else:
                    return 'intron'


