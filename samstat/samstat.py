"""Main file for SamStat
Authors: Will Patterson, Amie Romney

Description: TODO
"""
import pysam
import warnings
import timeit
from collections import namedtuple
from functools import lru_cache
from operator import itemgetter


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
    """
    Region = namedtuple('Region', ['features', 'length'])
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
                    try:
                        region_map[region].features[feature].append(location)
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
                                                  .append(location)
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
                largest = next(features)[1] 
                for feat in features:
                    if largest < feat[1]: largest = feat[1]
                region_map[key] = cls.Region(region.features, largest)
        return region_map

    @lru_cache(maxsize=None)
    def get_location_clasification(self,
                                   region_name,
                                   location_start,
                                   location_stop):
        """Gets location classification from region_map"""
        for key, ranges in self.region_map[region_name].features.items():
            for feat_range in ranges:
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

SamIn = namedtuple('InLine',
                   ['alignment_number',
                    'zeros',
                    'sixteens',
                    'cigar',
                    'rname_positions'])

def read_alignment_map(path):
    """Reads Alignment map SAM/BAM file into dictionary"""
    qnames = {}
    samfile = pysam.AlignmentFile(path, 'r')
    for count, seq_line in enumerate(samfile):
        zeros = 0
        sixteens = 0
        if seq_line.flag == 16:
            zeros = 0
            sixteens = 1
        elif seq_line.flag == 0:
            zeros = 1
            sixteens = 0
        qname = seq_line.query_name
        qnames.setdefault(qname, SamIn([0], [0], [0], seq_line.cigar, []))
        qnames[qname].alignment_number[0] += 1
        qnames[qname].zeros[0] += zeros
        qnames[qname].sixteens[0] += sixteens
        try:
            qnames[qname].rname_positions.append((seq_line.reference_name,
                                                  seq_line.reference_start))
        except ValueError:
            warnings.warn('Reference Name is -1, Line #: {}'.format(count))

    return qnames

OutLine = namedtuple('OutLine',
                     ['QNAME',
                      'alignment_number',
                      'zeros',
                      'sixteens',
                      'unique_rnames_low',
                      'unique_rnames_high',
                      'unique_rnames_number',
                      'gff_classification'])

def calculate_statistics(qname_data, region_map):
    """Calculates statistics using SAM and GFF data"""
    out_lines = []
    for qname, qdata in qname_data.items():
        alignment_number = qdata.alignment_number[0]
        zeros = qdata.zeros[0]
        sixteens = qdata.sixteens[0]
        unique_rnames = sorted({(x[1], qdata.rname_positions.count(x))\
                                for x in qdata.rname_positions},
                               key=itemgetter(1))
        unique_rnames_low = unique_rnames[0]
        unique_rnames_high = unique_rnames[-1]
        unique_rnames_number = len(unique_rnames)

        #TODO gff_classification
        out_lines.append(OutLine(qname,
                                 alignment_number,
                                 zeros,
                                 sixteens,
                                 unique_rnames_low,
                                 unique_rnames_high,
                                 unique_rnames_number,
                                 None))

    return out_lines



IN_GFF = '/disk/bioscratch/Will/Drop_Box/GCF_001266775.1_Austrofundulus_limnaeus-1.0_genomic_andMITO.gff'
IN_SAM = '/disk/bioscratch/Will/Drop_Box/HPF_small_RNA_022216.sam'

if __name__ == '__main__':
    start = timeit.default_timer()

    # Current:

    sam_data = read_alignment_map(IN_SAM)
    region_map = RegionMap(IN_GFF)
    calculate_statistics(sam_data, region_map)

    """ OLD:
    [print(x) for x in read_alignment_map(IN_SAM).items()]
    read_alignment_map(IN_SAM)
    print(len(RegionMap.read_gff(IN_GFF, ['exon']).keys()))
    print(read_alignment_map(IN_SAM).keys())
    """

    end = timeit.default_timer()
    print('Total Program Time: {}'.format(end-start))
