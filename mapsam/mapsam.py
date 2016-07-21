""" """
from collections import namedtuple, defaultdict
from BCBio import GFF

in_gff = '/disk/bioscratch/Will/Drop_Box/GCF_001266775.1_Austrofundulus_limnaeus-1.0_genomic_andMITO.gff'
in_sam = '/disk/bioscratch/Will/Drop_Box/HPF_small_RNA_022216.sam'

class ViviDict(dict):
    """Class for nice nested Hashes"""
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value

def split_gen(s, delims):
    start = 0
    for i, c in enumerate(s):
        if c in delims:
            yield s[start:i]
            start = i+1
    yield s[start:]

class RegionMap(object):

    Region = namedtuple('Region', ['name', 'subregions', 'start', 'end'])
    def __init__(self, gff_path, accepted_features=('exon')):
        self.subregion_types = subregion_types
        self.rmap = self.read_gff(gff_path)

    @staticmethod
    def make_subregion_template(subregion_types):
        subregion_template = {}
        for subregion in subregion_types:
            subregion_template.update({subregion: []})
        return subregion_template

    @staticmethod
    def read_gff(gff_path, feature_types):
        region_map = {}
        with open(in_gff, 'r') as gff:
            for line in gff:
                if line.startswith('#'):
                    pass
                line_gen = split_gen(line, '\t')
                region_name = line_gen.__next__()
                line_gen.__next__()
                feature = line_gen.__next__()
                if feature == 'region':
                    region_map.setdefault(region_name, Region(region_name, {key: [] for key in feature_types}, line_gen.__next__(), line_gen.__next__()))
                    region_map[region_name].feature_types.append((line_gen.__next__(), line_gen.__next__()))
                if feature in feature_types:
                    try:
                        region_map[region_name].feature_types.append((line_gen.__next__(), line_gen.__next__()))
                    except: #TODO
                        print('warning') #TODO
