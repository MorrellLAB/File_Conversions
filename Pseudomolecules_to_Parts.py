#!/usr/bin/env python
"""Convert from IPK full pseudomolecules to parts. Requires argparse."""

import sys
import argparse


# define the length of the part1 pieces as a dictionary constant.
PARTS_SIZES = {
    'chr1H': 312837513,
    'chr2H': 393532674,
    'chr3H': 394310633,
    'chr4H': 355061206,
    'chr5H': 380865482,
    'chr6H': 294822070,
    'chr7H': 325797516,
    'chrUn': 249774706,
    'Pt': 999999999999
    }


def parse_args():
    """Set up an argument parser, and actually parse the args."""
    parser = argparse.ArgumentParser(
        description='Convert pseudomolecule coordinates to parts coordinates',
        add_help=True)
    # Add a mutually exclusive group for the file format. This is so that only
    # one format can be specified.
    fmt = parser.add_mutually_exclusive_group(required=True)
    fmt.add_argument(
        '--vcf',
        '-v',
        action='store_true',
        help='Input file is VCF')
    fmt.add_argument(
        '--bed',
        '-b',
        action='store_true',
        help='Input file is BED')
    fmt.add_argument(
        '--gff',
        '-g',
        action='store_true',
        help='Input file is GFF')
    fmt.add_argument(
        '--reg',
        '-r',
        action='store_true',
        help='Input file is SAMTools/ANGSD regions')
    parser.add_argument(
        'intervals',
        metavar='file',
        help='File with intervals to convert')
    a = parser.parse_args()
    return a


def set_format(args):
    """Set the format that is being used for conversion."""
    arg_dict = vars(args)
    if arg_dict['vcf']:
        return 'VCF'
    elif arg_dict['bed']:
        return 'BED'
    elif arg_dict['gff']:
        return 'GFF'
    elif arg_dict['reg']:
        return 'REG'
    else:
        return None


def vcf_conv(intervals):
    """Convert VCF."""
    # If lines start with '#', print them
    with open(intervals, 'r') as f:
        for line in f:
            if line.startswith('#'):
                print line.strip()
            else:
                tmp = line.strip().split()
                chrom = tmp[0]
                pos = int(tmp[1])
                # Check that the chromosomes are named as we are expecting
                if chrom not in PARTS_SIZES:
                    sys.stderr.write(chrom + ' not reconized. The chromosomes must be named like \'chr1H.\'\n')
                    return
                # Check the parts lengths. If the position is greater than the
                # part1 length, then subtract the part1 length, and change the
                # name to part2
                # Also passing 'chrUn' unchanged, as there is only 1 part
                limit = PARTS_SIZES[chrom]
                if chrom == 'chrUn' or chrom == 'Pt':
                    newchrom = chrom
                    newpos = str(pos)
                elif pos > limit:
                    newchrom = chrom + '_part2'
                    newpos = str(pos - limit)
                else:
                    newchrom = chrom + '_part1'
                    newpos = str(pos)
                print '\t'.join([newchrom, newpos] + tmp[2:])
    return


def bed_conv(intervals):
    """Convert BED."""
    with open(intervals, 'r') as f:
        for line in f:
            if line.startswith('#'):
                print line.strip()
                continue
            tmp = line.strip().split()
            chrom = tmp[0]
            # Check that the chromosomes are named as we are expecting
            if chrom not in PARTS_SIZES:
                sys.stderr.write(chrom + ' not reconized. The chromosomes must be named like \'chr1H.\'\n')
                return
            start = int(tmp[1])
            end = int(tmp[2])
            # Check the coordinates of the start and stop positions with respect
            # to the part1 boundary
            limit = PARTS_SIZES[chrom]
            if chrom == 'ChrUn' or chrom == 'Pt':
                newchrom = chrom
                newstart = str(start)
                newend = str(end)
            elif start+1 > limit and end+1 > limit:
                newchrom = chrom + '_part2'
                newstart = str(start - limit)
                newend = str(end - limit)
            elif (start + 1) <= limit and (end + 1) > limit:
                sys.stderr.write('Interval ' + line.strip() + ' spans parts, omitting.\n')
                continue
            else:
                newchrom = chrom + '_part1'
                newstart = str(start)
                newend = str(end)
            print '\t'.join([newchrom, newstart, newend])
    return


def gff_conv(intervals):
    """Convert GFF."""
    with open(intervals, 'r') as f:
        for line in f:
            if line.startswith('#'):
                print line.strip()
                continue
            tmp = line.strip().split()
            chrom = tmp[0]
            start = int(tmp[3])
            end = int(tmp[4])
            if chrom not in PARTS_SIZES:
                sys.stderr.write(chrom + ' not reconized. The chromosomes must be named like \'chr1H.\'\n')
                return
            limit = PARTS_SIZES[chrom]
            if chrom == 'ChrUn' or chrom == 'Pt':
                newchrom = chrom
                newstart = str(start)
                newend = str(end)
            elif start > limit and end > limit:
                newchrom = chrom + '_part2'
                newstart = str(start - limit)
                newend = str(end - limit)
            elif start+1 <= limit and end+1 > limit:
                sys.stderr.write('Interval ' + line.strip() + ' spans parts, omitting.\n')
                continue
            else:
                newchrom = chrom + '_part1'
                newstart = str(start)
                newend = str(end)
            print '\t'.join([newchrom, tmp[1], tmp[2], newstart, newend] + tmp[5:])
    return


def reg_conv(intervals):
    """Convert SAM region."""
    with open(intervals, 'r') as f:
        for line in f:
            if ':' not in line or '-' not in line:
                sys.stderr.write(line.strip() + ' is not formatted properly. This script expects chr:start-stop format.\n')
                return
            tmp = line.strip().split(':')
            chrom = tmp[0]
            start, end = [int(i) for i in tmp[1].split('-')]
            limit = PARTS_SIZES[chrom]
            if chrom == 'ChrUn' or chrom == 'Pt':
                newchrom = chrom
                newstart = str(start)
                newend = str(end)
            elif start > limit and end > limit:
                newchrom = chrom + '_part2'
                newstart = str(start - limit)
                newend = str(end - limit)
            elif start+1 <= limit and end+1 > limit:
                sys.stderr.write('Interval ' + line.strip() + ' spans parts, omitting.\n')
                continue
            else:
                newchrom = chrom + '_part1'
                newstart = str(start)
                newend = str(end)
            print newchrom + ':' + newstart + '-' + newend
    return


def main():
    """Main function."""
    args = parse_args()
    f = set_format(args)
    # Then, use the format and the supplied file to convert the coordinates
    if f == 'VCF':
        vcf_conv(args.intervals)
    elif f == 'BED':
        bed_conv(args.intervals)
    elif f == 'GFF':
        gff_conv(args.intervals)
    elif f == 'REG':
        reg_conv(args.intervals)
    else:
        print 'Unrecognized format, weird.'
        exit(1)
    return


main()
