import timeit
import sys
import warnings
sys.path.append('..')

from collections import namedtuple
from maps import RegionMap
from maps import AlignmentMap

OutLine = namedtuple('OutLine', ['qname', 'rname', 'forward', 'reverse'])

def calculate_statistics(alignment_map, region_map):
    """ """
    for qname, qdata in alignment_map.items():
        rnames = dict()
        for rname, location, direction in qdata.reference_names:
            rnames.setdefault(rname, [0, 0])
            try:
                direction_count = region_map.get_true_directions(rname, (location, location+qdata.cigar[0][1]-1), direction)
            except TypeError:
                warnings.warn('Invalid direction type')
                continue
            rnames[rname][0] += direction_count.forwards
            rnames[rname][1] += direction_count.reverses
        for rname, directions in rnames.items():
            yield OutLine(qname, rname, directions[0], directions[1])

def format_line_obj(line_obj, ordered_attributes, delimiter):
    """Formats an object into delimited string

    If attribute not found, a ':(' will be inserted instead
    """
    return str(delimiter).join((str(line_obj.__dict__.get(attr, ':(')) for attr in ordered_attributes))


def run(in_sam, in_gff, outpath):
    """Runs SamStat functions"""
    sam_data = AlignmentMap(in_sam)
    region_map = RegionMap(in_gff)
    in_attributes = ['qname', 'rname', 'forward', 'reverse']
    with open(outpath, 'w') as ofile:
        ofile.write('\n'.join([format_line_obj(oline, in_attributes, '\t') for oline in calculate_statistics(sam_data, region_map)]))

def main(args=None):
    """Comamnd line argument for rname stats"""
    pass

IN_GFF = '/disk/bioscratch/Will/Drop_Box/GCF_001266775.1_Austrofundulus_limnaeus-1.0_genomic_andMITO.gff'
IN_SAM = '/disk/bioscratch/Will/Drop_Box/HPF_small_RNA_022216.sam'
OUT_CSV = '/disk/bioscratch/Will/Drop_Box/truedir_output.v1.csv'


if __name__ == '__main__':
    start = timeit.default_timer()
    run(IN_SAM, IN_GFF, OUT_CSV)
    end = timeit.default_timer()
    print('Total Program Time: {}'.format(end-start))
