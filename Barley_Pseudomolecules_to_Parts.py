#!/usr/bin/env python3
"""Convert from IPK full pseudomolecules to parts. Requires argparse."""

import sys
import argparse
import cyvcf2
from cyvcf2 import VCF


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
    """Convert VCF. Works for both GATK vcfs and larger SVs with "INFO/END" field
    from callers like Sniffles2, cuteSV, Longranger."""
    # If lines start with '#', print them
    vcf = VCF(intervals)
    # Print header lines to stdout
    # Update header with new chromosome part names and lengths
    for l in vcf.raw_header.strip('\n').split('\n'):
        if l.startswith('##contig'):
            if 'chrUn' in l or 'EF115541.1' in l or 'AP017301.1' in l:
                # Print as is, no modifications needed
                print(l)
            else:
                # Split into chromosome parts
                contig = l.split('=')
                chrom_len = contig[3].strip('>')
                chrom_list = contig[2].split(',')
                # Prepare new chromosome parts names
                chrom_part1_name = chrom_list[0] + '_part1'
                chrom_part2_name = chrom_list[0] + '_part2'
                # Get new parts chromosome lengths
                chrom_part1_len = parts_sizes[chrom_list[0]]
                chrom_part2_len = int(chrom_len) - int(chrom_part1_len)
                # Check if chromosome parts lengths add up correctly
                if str(chrom_part1_len + chrom_part2_len) == str(chrom_len):
                    # Prepare updated contig header lines
                    # Each chromosome header line should get split into two parts
                    print('='.join([contig[0], contig[1], chrom_part1_name + ',length', str(chrom_part1_len) + '>']))
                    print('='.join([contig[0], contig[1], chrom_part2_name + ',length', str(chrom_part2_len) + '>']))
        else:
            print(l)
    # Check if we have INFO/END field present in the VCF
    # If so, we'll need to pull the end position from the END field
    # This is a larger SV with the end position in the INFO/END field
    # (e.g., output from callers like Sniffles2, cuteSV, and Longranger)
    if vcf.contains("END"):
        for variant in vcf:
            chrom = variant.CHROM
            start_pos = variant.POS
            # Position in INFO/END field
            end_pos = variant.INFO["END"]
            # Check that the chromosomes are named as we are expecting
            if chrom not in parts_sizes:
                sys.stderr.write(chrom + ' not recognized. The chromosomes must be named like \'chr1H.\'\n')
                return
            # Check the parts lengths. If the position is greater than the
            # part1 length, then subtract the part1 length, and change the
            # name to part2
            # Also passing 'chrUn' unchanged, as there is only 1 part
            limit = parts_sizes[chrom]
            if chrom == 'chrUn' or chrom == 'Pt':
                # No changes made, print as is
                print(str(variant).strip('\n'))
            elif start_pos > limit and end_pos > limit:
                # This is important, otherwise conversion for records that span
                # both chrom parts will be incorrect
                # Modify in current record
                # New chromosome name
                variant.CHROM = newchrom = chrom + '_part2'
                # New start position
                variant.set_pos(start_pos - limit)
                # New end position
                variant.INFO["END"] = end_pos - limit
                # Print modified record to stdout
                print(str(variant).strip('\n'))
            elif start_pos < limit and end_pos <= limit:
                # Modify in current record
                # New chromosome name
                variant.CHROM = chrom + '_part1'
                # Print modified record to stdout
                print(str(variant).strip('\n'))
            elif start_pos <= limit and end_pos > limit:
                # Record spans both chrom parts, exit with error
                sys.stderr.write(chrom + ':' + str(start_pos) + '-' + str(end_pos) + ' spans chrom parts, please check manually before proceeding.\n')
                return
            else:
                sys.stderr.write(chrom + ':' + str(start_pos) + '-' + str(end_pos) + ' is an edge case that is not dealt with properly, please check manually and modify code if necessary before proceeding.\n')
                return
    else:
        # no INFO/END field present (e.g., vcf output from GATK)
        for variant in vcf:
            chrom = variant.CHROM
            pos = variant.POS
            # Check that the chromosomes are named as we are expecting
            if chrom not in parts_sizes:
                sys.stderr.write(chrom + ' not recognized. The chromosomes must be named like \'chr1H.\'\n')
                return
            # Check the parts lengths. If the position is greater than the
            # part1 length, then subtract the part1 length, and change the
            # name to part2
            # Also passing 'chrUn' unchanged, as there is only 1 part
            limit = parts_sizes[chrom]
            if chrom == 'chrUn' or chrom == 'Pt':
                # No changes made, print as is
                print(str(variant).strip('\n'))
            elif pos > limit:
                # Modify in current record
                # New chromosome name
                variant.CHROM = newchrom = chrom + '_part2'
                # New start position
                variant.set_pos(pos - limit)
                # Print modified record to stdout
                print(str(variant).strip('\n'))
            elif pos <= limit:
                # Modify in current record
                # New chromosome name
                variant.CHROM = chrom + '_part1'
                # Print modified record to stdout
                print(str(variant).strip('\n'))
            else:
                sys.stderr.write(chrom + ':' + str(start_pos) + '-' + str(end_pos) + ' is an edge case that is not dealt with properly, please check manually and modify code if necessary before proceeding.\n')
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
            if chrom == 'chrUn' or chrom == 'Pt':
                newchrom = chrom
                newstart = str(start)
                newend = str(end)
            elif (start + 1) > limit and (end + 1) > limit:
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
            tmp = line.strip().split('\t')
            chrom = tmp[0]
            start = int(tmp[3])
            end = int(tmp[4])
            if chrom not in parts_sizes:
                sys.stderr.write(chrom + ' not reconized. The chromosomes must be named like \'chr1H.\'\n')
                return
            limit = parts_sizes[chrom]
            if chrom == 'chrUn' or chrom == 'Pt':
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
            if chrom == 'chrUn' or chrom == 'Pt':
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
            if chrom == 'chrUn' or chrom == 'Pt':
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
