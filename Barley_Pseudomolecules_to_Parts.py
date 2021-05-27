#!/usr/bin/env python3
"""Convert from IPK full pseudomolecules to parts. Requires argparse."""

import sys
import argparse


# Define the length of the part1 pieces as a dictionary constant.
# Morex v1 sizes
PARTS_SIZES_V1 = {
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
# Morex v2 sizes including plastids
PARTS_SIZES_V2 = {
    'chr1H': 205502676,
    'chr2H': 305853815,
    'chr3H': 271947776,
    'chr4H': 282386439,
    'chr5H': 205989812,
    'chr6H': 260041240,
    'chr7H': 328847863,
    'chrUn': 85026395,
    'EF115541.1': 136462,
    'AP017301.1': 525599,
    'Pt': 999999999999
    }
# Morex v3 sizes
PARTS_SIZES_V3 = {
    'chr1H': 206486643,
    'chr2H': 301293086,
    'chr3H': 267852507,
    'chr4H': 276149121,
    'chr5H': 204878572,
    'chr6H': 256319444,
    'chr7H': 328847192,
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
        '--gtf',
        '-t',
        action='store_true',
        help='Input file is GTF')
    fmt.add_argument(
        '--reg',
        '-r',
        action='store_true',
        help='Input file is SAMTools/ANGSD regions')
    parser.add_argument(
        'intervals',
        metavar='file',
        help='File with intervals to convert')
    parser.add_argument(
        'ref_version',
        metavar='reference_version',
        help='Input is one of the following: morex_v1, morex_v2, morex_v3')
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
    elif arg_dict['gtf']:
        return 'GTF'
    elif arg_dict['reg']:
        return 'REG'
    else:
        return None


def vcf_conv(intervals, parts_sizes):
    """Convert VCF."""
    # If lines start with '#', print them
    with open(intervals, 'r') as f:
        for line in f:
            if line.startswith('#'):
                print(line.strip())
            else:
                tmp = line.strip().split()
                chrom = tmp[0]
                pos = int(tmp[1])
                # Check that the chromosomes are named as we are expecting
                if chrom not in parts_sizes:
                    sys.stderr.write(chrom + ' not reconized. The chromosomes must be named like \'chr1H.\'\n')
                    return
                # Check the parts lengths. If the position is greater than the
                # part1 length, then subtract the part1 length, and change the
                # name to part2
                # Also passing 'chrUn' unchanged, as there is only 1 part
                limit = parts_sizes[chrom]
                if chrom == 'chrUn' or chrom == 'Pt':
                    newchrom = chrom
                    newpos = str(pos)
                elif pos > limit:
                    newchrom = chrom + '_part2'
                    newpos = str(pos - limit)
                else:
                    newchrom = chrom + '_part1'
                    newpos = str(pos)
                print('\t'.join([newchrom, newpos] + tmp[2:]))
    return


def bed_conv(intervals, parts_sizes):
    """Convert BED."""
    with open(intervals, 'r') as f:
        for line in f:
            if line.startswith('#'):
                print(line.strip())
                continue
            tmp = line.strip().split()
            chrom = tmp[0]
            # Check that the chromosomes are named as we are expecting
            if chrom not in parts_sizes:
                sys.stderr.write(chrom + ' not reconized. The chromosomes must be named like \'chr1H.\'\n')
                return
            start = int(tmp[1])
            end = int(tmp[2])
            # Check the coordinates of the start and stop positions with respect
            # to the part1 boundary
            limit = parts_sizes[chrom]
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
            print('\t'.join([newchrom, newstart, newend]))
    return


def gff_conv(intervals, parts_sizes):
    """Convert GFF."""
    with open(intervals, 'r') as f:
        for line in f:
            if line.startswith('#'):
                if line.startswith('##sequence-region'):
                    tmp_header = line.strip().split()
                    curr_chrom = tmp_header[1]
                    start_pos = tmp_header[2]
                    end_pos = tmp_header[3]
                    if curr_chrom != "chrUn":
                        # Generate ##sequence-region header for split parts
                        curr_part1_size = parts_sizes[curr_chrom]
                        # Print part1
                        print(tmp_header[0], "   ", curr_chrom + '_part1 ', str(start_pos), ' ', str(curr_part1_size))
                        # Print part2
                        print(tmp_header[0], "   ", curr_chrom + '_part2 ', str(curr_part1_size + 1), ' ', str(end_pos))
                    else:
                        # chrUn usually isn't split
                        print(line.strip())
                else:
                    print(line.strip())
                continue
            tmp = line.strip().split()
            chrom = tmp[0]
            start = int(tmp[3])
            end = int(tmp[4])
            if chrom not in parts_sizes:
                sys.stderr.write(chrom + ' not reconized. The chromosomes must be named like \'chr1H.\'\n')
                return
            limit = parts_sizes[chrom]
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
            print('\t'.join([newchrom, tmp[1], tmp[2], newstart, newend] + tmp[5:]))
    return

def gtf_conv(intervals, parts_sizes):
    """Convert GFF."""
    with open(intervals, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            chrom = tmp[0]
            start = int(tmp[3])
            end = int(tmp[4])
            if chrom not in parts_sizes:
                sys.stderr.write(chrom + ' not reconized. The chromosomes must be named like \'chr1H.\'\n')
                return
            limit = parts_sizes[chrom]
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
            print('\t'.join([newchrom, tmp[1], tmp[2], newstart, newend] + tmp[5:]))
    return

def reg_conv(intervals, parts_sizes):
    """Convert SAM region."""
    with open(intervals, 'r') as f:
        for line in f:
            if ':' not in line or '-' not in line:
                sys.stderr.write(line.strip() + ' is not formatted properly. This script expects chr:start-stop format.\n')
                return
            tmp = line.strip().split(':')
            chrom = tmp[0]
            start, end = [int(i) for i in tmp[1].split('-')]
            limit = parts_sizes[chrom]
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
            print(newchrom + ':' + newstart + '-' + newend)
    return


def main():
    """Main function."""
    args = parse_args()
    f = set_format(args)
    # Check which version of the Morex reference we are using
    if args.ref_version == "morex_v1":
        parts_sizes = PARTS_SIZES_V1
    elif args.ref_version == "morex_v2":
        parts_sizes = PARTS_SIZES_V2
    elif args.ref_version == "morex_v3":
        parts_sizes = PARTS_SIZES_V3
    else:
        print('Invalid reference version provided, valid options are: morex_v1, morex_v2, morex_v3')
        exit(1)
    
    # Then, use the format and the supplied file to convert the coordinates
    if f == 'VCF':
        vcf_conv(args.intervals, parts_sizes)
    elif f == 'BED':
        bed_conv(args.intervals, parts_sizes)
    elif f == 'GFF':
        gff_conv(args.intervals, parts_sizes)
    elif f == 'GTF':
        gtf_conv(args.intervals, parts_sizes)
    elif f == 'REG':
        reg_conv(args.intervals, parts_sizes)
    else:
        print('Unrecognized format, weird.')
        exit(1)
    return


main()
