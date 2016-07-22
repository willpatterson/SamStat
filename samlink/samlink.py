""" """
import warnings
import timeit
from collections import namedtuple, defaultdict


def split_gen(s, delims):
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
    Region = namedtuple('Region', ['name', 'subregions', 'length'])
    def __init__(self, gff_path, accepted_features=('exon')):
        self.subregion_types = subregion_types
        self.rmap = self.read_gff(gff_path)

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
                        region_map[region].subregions[feature].append(location)
                    except KeyError as e:
                        if feature == 'region':
                            #print('create regoin' + region)
                            region_map.setdefault(region,
                                                  cls.Region(region,
                                                             feature_temp,
                                                             location[1]))
                        else:
                            region_map.setdefault(region,
                                                  cls.Region(region,
                                                             feature_temp,
                                                             None))
                            try:
                                region_map[region].subregions[feature]\
                                                  .append(location)
                            except KeyError:
                                pass

                except StopIteration:
                    warnings.warn('Invalid line: {} ... skipped'.format(count))
            return region_map

IN_GFF = '/disk/bioscratch/Will/Drop_Box/GCF_001266775.1_Austrofundulus_limnaeus-1.0_genomic_andMITO.gff'
IN_SAM = '/disk/bioscratch/Will/Drop_Box/HPF_small_RNA_022216.sam'

if __name__ == '__main__':
    start = timeit.default_timer()
    print(len(RegionMap.read_gff(IN_GFF, ['exon']).keys()))
    end = timeit.default_timer()
    print(end-start)

