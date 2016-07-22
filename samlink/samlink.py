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

    @staticmethod
    def read_gff(gff_path, feature_types):
        """Reads a gff file into a region map hash"""
        region_map = {}
        feature_template = {key: [] for key in feature_types}
        print(feature_template)
        with open(gff_path, 'r') as gff:
            for line in gff:
                if line.startswith('#'):
                    continue
                line_gen = split_gen(line, '\t')
                #print([x for x in line_gen])
                region_name = line_gen.__next__()
                line_gen.__next__()
                feature = line_gen.__next__()
                if feature == 'region':
                    region_map.setdefault(region_name, RegionMap.Region(region_name, feature_template, line_gen.__next__(), line_gen.__next__()))
                elif feature in feature_types:
                    #try:
                    region_map[region_name].subregions[feature].append([int(line_gen.__next__()), int(line_gen.__next__())]) #except: #TODO
                    #print('warning') #TODO
        return region_map

IN_GFF = '/disk/bioscratch/Will/Drop_Box/GCF_001266775.1_Austrofundulus_limnaeus-1.0_genomic_andMITO.gff'
IN_SAM = '/disk/bioscratch/Will/Drop_Box/HPF_small_RNA_022216.sam'

if __name__ == '__main__':
    print(RegionMap.read_gff(IN_GFF, ['exon']))

