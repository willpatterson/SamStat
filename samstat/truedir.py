import sys
sys.path.append('..')

from maps import RegionMap
from maps import AlignmentMap

OutLine = namedtuple('OutLine', ['qname', 'rname', 'forward', 'reverse'])

def calculate_statistics(alignment_map, region_map):
    """ """
    #for qname, qdata in alignment_map.items():
    #    for qdata.


def run(in_sam, in_gff, outpath):
    """Runs the rname stats script"""
    amap = AlignmentMap(in_sam)
    rmap = RegionMap(in_gff)

def main(args=None):
    """Comamnd line argument for rname stats"""
    pass

if __name__ == '__main__':
    pass
