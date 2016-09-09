import argparse
import os
import pysam
import warnings
import timeit
from collections import namedtuple
from operator import itemgetter

import sys
sys.path.append('..')

from maps import RegionMap
from maps import AlignmentMap


true_dir_out_values = ['qname', 'rname', 'forward', 'reverse']

def calculate_truedirs(alignment_map, region_map):
    """ """
    TrueDirOutValues = namedtuple('TrueDirOutValues', out_values)
    for qname, qdata in alignment_map.items():
        rnames = dict()
        for rname, location, direction in qdata.reference_names:
            rnames.setdefault(rname, [0, 0])
            try:
                direction_count = region_map.get_true_directions(rname, (location, location+qdata.cigar[0][1]-1), direction)
            except TypeError as e:
                warnings.warn('Invalid direction type')
                print(e)
                continue
            rnames[rname][0] += direction_count.forwards
            rnames[rname][1] += direction_count.reverses
        for rname, directions in rnames.items():
            yield TrueDirOutValues(qname, rname, directions[0], directions[1])


qstat_out_values = ['qname',
                    'alignment_number',
                    'unique_rnames_low',
                    'unique_rnames_high',
                    'unique_rnames_number',
                    'exons',
                    'introns',
                    'intergenes',
                    'combos']

def calculate_qstats(qname_data, region_map):
    """Calculates statistics using SAM and GFF data"""
    for qname, qdata in qname_data.items():
        alignment_number = qdata.alignment_number[0]
        try:
            raw_rnames = [x[0] for x in qdata.reference_names]
            unique_rnames = sorted({(x, raw_rnames.count(x)) for x in raw_rnames},
                                   key=itemgetter(1))
            unique_rnames_low = {x[0] for x in unique_rnames if x[1] == unique_rnames[0][1]}
            unique_rnames_high = {x[0] for x in unique_rnames[::-1] if x[1] == unique_rnames[-1][1]}
            unique_rnames_number = len(unique_rnames)
            if unique_rnames_low == unique_rnames_high:
                unique_rnames_low = unique_rnames_high = 'All{}'.format(unique_rnames[0][1])
            else:
                unique_rnames_low = len(unique_rnames_low)
                unique_rnames_high = len(unique_rnames_high)

            gff_classes = [region_map.get_location_clasification(rname, location, location+qdata.cigar[0][1]-1) for rname, location, _ in qdata.reference_names]

            total_exons, total_introns, total_combos, total_intergenes = 0, 0, 0, 0
            for classification in gff_classes:
                total_exons += classification.exons
                total_introns += classification.introns
                total_intergenes += classification.intergene
                total_combos += classification.combos

            yield QstatOutValues(qname,
                          alignment_number,
                          unique_rnames_low,
                          unique_rnames_high,
                          unique_rnames_number,
                          total_exons,
                          total_introns,
                          total_intergenes,
                          total_combos)

        except IndexError:
            warnings.warn('QNAME data cannot be read, Skipping: {}'.format(qname))

def format_line_obj(line_obj, ordered_attributes, delimiter):
    """Formats an object into delimited string

    If attribute not found, a ':(' will be inserted instead
    """
    return str(delimiter).join((str(line_obj.__dict__.get(attr, ':(')) for attr in ordered_attributes))


def run(in_sam, in_gff, outpath, out_values, run_function):
    """Runs SamStat functions"""
    sam_data = AlignmentMap(in_sam)
    region_map = RegionMap(in_gff)
    with open(outpath, 'w') as ofile:
        ofile.write('\n'.join([format_line_obj(oline, out_values, '\t') for oline in run_function(sam_data, region_map)]))

def main():
    """Command line interface for samstat"""
    parser = argparse.ArgumentParser(description="TODO")
    parser.add_argument('operation', type=str, help='Operation to preform on data (qstat, truedir)')
    parser.add_argument('sam_file', type=str, help='Path To SAM input file')
    parser.add_argument('gff_file', type=str, help='Path to GFF input file')
    parser.add_argument('out_path', type=str, help='Path of output file')
    args = parser.parse_args()

    if not os.path.isfile(args.sam_file):
        raise argparse.ArgumentTypeError('sam_file is not a valid file path')
    if not os.path.isfile(args.gff_file):
        raise argparse.ArgumentTypeError('gff_file is not a valid file path')
    if os.path.exists(args.out_path):
        warnings.warn(('Warning output path already exists, data will be '
                       'overwritten. Path {}').format(args.out_path))

    if args.operation is 'qstat':
        run(args.sam_file, args.gff_file, args.out_path, qstat_out_values, calculate_qstats)
    elif args.operation is 'truedir':
        run(args.sam_file, args.gff_file, args.out_path, true_dir_out_values, calculate_truedirs)

IN_GFF = '/disk/bioscratch/Will/Drop_Box/GCF_001266775.1_Austrofundulus_limnaeus-1.0_genomic_andMITO.gff'
IN_SAM = '/disk/bioscratch/Will/Drop_Box/HPF_small_RNA_022216.sam'
OUT_TRUEDIR = '/disk/bioscratch/Will/Drop_Box/truedir_output.v1-0.csv'
OUT_QSTAT = '/disk/bioscratch/Will/Drop_Box/qstat_output.v1-0.csv'


if __name__ == '__main__':
    start = timeit.default_timer()
    run(IN_SAM, IN_GFF, OUT_CSV)
    end = timeit.default_timer()
    print('Total Program Time: {}'.format(end-start))
