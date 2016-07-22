""" """
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
    Region = namedtuple('Region', ['name', 'subregions', 'start', 'end'])
    def __init__(self, gff_path, accepted_features=('exon')):
        self.subregion_types = subregion_types
        self.rmap = self.read_gff(gff_path)

    @classmethod
    def read_gff(cls, gff_path, feature_types):
        """Reads a gff file into a region map hash"""
        region_map = {}
        feature_temp = {key: [] for key in feature_types}
        print(feature_temp)
        with open(gff_path, 'r') as gff:
            for count, line in enumerate(gff):
                try:
                    if line.startswith('#'):
                        print('# line: {}'.format(count))
                        continue
                    line_gen = split_gen(line, '\t').next
                    #print([x for x in line_gen])
                    region_name = line_gen()
                    line_gen()
                    feature = line_gen()
                    if feature == 'region':
                        region_map.setdefault(region_name,
                                              cls.Region(region_name,
                                                         feature_temp,
                                                         int(line_gen()),
                                                         int(line_gen())))
                    elif feature in feature_types:
                        try:
                            region_map[region_name].subregions[feature]\
                                                   .append([int(line_gen()),
                                                            int(line_gen())])
                        except Exception as e: #TODO
                            print('inner try: line {} count {}'.format(line,
                                                                       count))
                            print(e.message)
                except Exception as e:
                    print('outer try: line {} count {}'.format(line, count))
                    print(e.message)
            return region_map

IN_GFF = '/disk/bioscratch/Will/Drop_Box/GCF_001266775.1_Austrofundulus_limnaeus-1.0_genomic_andMITO.gff'
IN_SAM = '/disk/bioscratch/Will/Drop_Box/HPF_small_RNA_022216.sam'

if __name__ == '__main__':
    RegionMap.read_gff(IN_GFF, ['exon'])

