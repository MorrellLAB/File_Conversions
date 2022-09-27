#!/usr/bin/env python3
"""A simple script to convert the pseudomolecule parts into full
pseudomolecules."""

import argparse
import gzip

#   Store the lengths of the barley pseudomolecule parts in a dictionary. We only
#   really need to store the first parts, since we will just add the part2
#   positions to it.
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
    'chrUn': 29110253,
    'EF115541.1': 136462,
    'AP017301.1': 525599,
    'Pt': 999999999999
}

def parse_args():
    """Set up an argument parser, and actually parse the args."""
    parser = argparse.ArgumentParser(
        description='Convert parts coordinates to pseudomolecule coordinates.',
        add_help=True
    )
    # Add mutually exclusive group for the file format. This is so that only
    # one format can be specified.
    fmt = parser.add_mutually_exclusive_group(required=True)
    fmt.add_argument(
        '--vcf',
        '-v',
        action='store_true',
        help='Input file is VCF or VCF.gz'
    )
    fmt.add_argument(
        '--bed',
        '-b',
        action='store_true',
        help='Input file is BED'
    )
    parser.add_argument(
        'intervals',
        metavar='file',
        help='File with intervals to convert'
    )
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
    else:
        return None


def vcf_conv(intervals, parts_sizes):
    """Convert VCF. Remember that VCF is 1-based, which makes the math a lot nicer.
    Handles both VCF and gzipped VCFs."""
    if '.gz' in intervals:
        # If lines start with '#', print them
        with gzip.open(intervals, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    print(line.strip())
                else:
                    tmp = line.strip().split()
                    # modify the chromosome to not have the part
                    chrom = tmp[0].split('_')[0]
                    # And check if we have to modify the position. If the chromosome
                    # name has '_part2' in it, then we have to modify it.
                    if '_part2' in tmp[0]:
                        #offset = parts_sizes[chrom + '_part1']
                        offset = parts_sizes[chrom]
                        newpos = str(int(tmp[1]) + offset)
                    else:
                        newpos = tmp[1]
                    # then print out the VCF line with the new position and chromosome
                    # name
                    toprint = [chrom, newpos] + tmp[2:]
                    print('\t'.join(toprint))
    else:
        # If lines start with '#', print them
        with open(intervals, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    print(line.strip())
                else:
                    tmp = line.strip().split()
                    # modify the chromosome to not have the part
                    chrom = tmp[0].split('_')[0]
                    # And check if we have to modify the position. If the chromosome
                    # name has '_part2' in it, then we have to modify it.
                    if '_part2' in tmp[0]:
                        #offset = parts_sizes[chrom + '_part1']
                        offset = parts_sizes[chrom]
                        newpos = str(int(tmp[1]) + offset)
                    else:
                        newpos = tmp[1]
                    # then print out the VCF line with the new position and chromosome
                    # name
                    toprint = [chrom, newpos] + tmp[2:]
                    print('\t'.join(toprint))


def bed_conv(intervals, parts_sizes):
    """Convert BED. Remember that BED is 0-based (i.e., the first base is 0)."""
    # If lines start with '#', print them
    with open(intervals, 'r') as f:
        for line in f:
            if line.startswith('#'):
                print(line.strip())
                continue
            else:
                tmp = line.strip().split()
                # modify the chromosome to not have th part
                chrom = tmp[0].split('_')[0]
                # And check if we have to modify the position. If the chromosome
                # name has '_part2' in it, then we have to modify it.
                if '_part2' in tmp[0]:
                    #offset = parts_sizes[chrom + '_part1']
                    offset = parts_sizes[chrom]
                    newstartpos = str(int(tmp[1]) + offset)
                    newendpos = str(int(tmp[2]) + offset)
                else:
                    newstartpos = tmp[1]
                    newendpos = tmp[2]
                # Check if standard BED file format with 3 columns
                if len(tmp) == 3:
                    # Then print out BED line with chromosome name and new start and end positions
                    toprint = [chrom, newstartpos, newendpos]
                    print('\t'.join(toprint))
                elif len(tmp) > 3:
                    # Print out BED line with chromosome name and new start and end positions
                    # along with remaining columns
                    toprint = [chrom, newstartpos, newendpos] + tmp[3:]
                    print('\t'.join(toprint))
                else:
                    print('Too few columns, please check input file is in BED format.')
                    exit(1)


def main():
    """Driver function."""
    args = parse_args()
    f = set_format(args)
    # Check which version of the Morex reference we are using
    if args.ref_version == "morex_v1":
        PARTS_SIZES = PARTS_SIZES_V1
    elif args.ref_version == "morex_v2":
        PARTS_SIZES = PARTS_SIZES_V2
    elif args.ref_version == "morex_v3":
        PARTS_SIZES = PARTS_SIZES_V3
    else:
        print('Invalid reference version provided, valid options are: morex_v1, morex_v2, morex_v3')
        exit(1)
    
    # Then, use the format and the supplied file to convert the coordinates
    if f == 'VCF':
        vcf_conv(args.intervals, PARTS_SIZES)
    elif f == 'BED':
        bed_conv(args.intervals, PARTS_SIZES)
    else:
        print('Unrecognized format, weird.')
        exit(1)


main() # Run the program
