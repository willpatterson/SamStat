"""Main file for SamStat
Authors: Will Patterson, Amie Romney

Description: TODO
"""
import argparse
import os
import pysam
import warnings
import timeit
from collections import namedtuple
from operator import itemgetter

import sys
sys.path.append('..')

from region_map import RegionMap

SamIn = namedtuple('InLine',
                   ['alignment_number',
                    'cigar',
                    'rname_positions'])

def read_alignment_map(path):
    """Reads Alignment map SAM/BAM file into dictionary"""
    qnames = {}
    samfile = pysam.AlignmentFile(path, 'r')
    for count, seq_line in enumerate(samfile):
        qname = seq_line.query_name
        qnames.setdefault(qname, SamIn([0], [0], [0], seq_line.cigar, []))
        qnames[qname].alignment_number[0] += 1
        try:
            qnames[qname].rname_positions.append((seq_line.reference_name,
                                                  seq_line.reference_start))
        except ValueError:
            warnings.warn('Reference Name is -1, Line #: {}'.format(count))

    return qnames

OutLine = namedtuple('OutLine',
                     ['qname',
                      'alignment_number',
                      'unique_rnames_low',
                      'unique_rnames_high',
                      'unique_rnames_number',
                      'exons',
                      'introns',
                      'combos'])


def calculate_statistics(qname_data, region_map):
    """Calculates statistics using SAM and GFF data"""
    for qname, qdata in qname_data.items():
        alignment_number = qdata.alignment_number[0]
        try:
            raw_rnames = [x[0] for x in qdata.rname_positions]
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

            gff_classes = [region_map.get_location_clasification(rname, location, location+qdata.cigar[0][1]-1) for rname, location in qdata.rname_positions]
            gff_classes = {x: gff_classes.count(x) for x in gff_classes}

            #TODO gff_classification
            yield OutLine(qname,
                          alignment_number,
                          unique_rnames_low,
                          unique_rnames_high,
                          unique_rnames_number,
                          gff_classes.get('exon', 0),
                          gff_classes.get('intron', 0),
                          gff_classes.get('combo', 0))

        except IndexError:
            warnings.warn('QNAME data cannot be read, Skipping: {}'.format(qname))

def format_line_obj(line_obj, ordered_attributes, delimiter):
    """Formats an object into delimited string

    If attribute not found, a ':(' will be inserted instead
    """
    return str(delimiter).join((str(line_obj.__dict__.get(attr, ':(')) for attr in ordered_attributes))

def run(in_sam, in_gff, outpath):
    """Runs SamStat functions"""
    sam_data = read_alignment_map(in_sam)
    region_map = RegionMap(in_gff)
    in_attributes = ['qname',
                     'alignment_number',
                     'unique_rnames_low',
                     'unique_rnames_high',
                     'unique_rnames_number',
                     'exons',
                     'introns',
                     'combos']
    with open(outpath, 'w') as ofile:
        ofile.write('\n'.join([format_line_obj(oline, in_attributes, '\t') for oline in calculate_statistics(sam_data, region_map)]))

def main():
    """Command line interface for samstat"""
    parser = argparse.ArgumentParser(description="TODO")
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
    run(args.sam_file, args.gff_file, args.out_path)

IN_GFF = '/disk/bioscratch/Will/Drop_Box/GCF_001266775.1_Austrofundulus_limnaeus-1.0_genomic_andMITO.gff'
IN_SAM = '/disk/bioscratch/Will/Drop_Box/HPF_small_RNA_022216.sam'
OUT_CSV = '/disk/bioscratch/Will/Drop_Box/SamStat_output.v2.csv'

if __name__ == '__main__':
    start = timeit.default_timer()
    run(IN_SAM, IN_GFF, OUT_CSV)
    end = timeit.default_timer()
    print('Total Program Time: {}'.format(end-start))
