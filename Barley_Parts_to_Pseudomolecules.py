#!/usr/bin/env python3
"""A simple script to convert the pseudomolecule parts into full
pseudomolecules. Remember that VCF is 1-based, which makes the math
a lot nicer."""

import sys
import argparse

#   Store the lengths of the barley pseudomolecule parts in a dictionary. We only
#   really need to store the first parts, since we will just add the part2
#   positions to it.
#   Sizes listed below are for Morex v2
PARTS_SIZES = {
    'chr1H_part1': 205502676,
    'chr2H_part1': 305853815,
    'chr3H_part1': 271947776,
    'chr4H_part1': 282386439,
    'chr5H_part1': 205989812,
    'chr6H_part1': 260041240,
    'chr7H_part1': 328847863
    }

#   Then iterate through the VCF and print out the modified lines
with open(sys.argv[1], 'r') as f:
    for line in f:
        if line.startswith('#'):
            print line.strip()
        else:
            tmp = line.strip().split('\t')
            #   modify the chromosome to not have the part
            chrom = tmp[0].split('_')[0]
            #   And check if we have to modify the position. If the chromosome
            #   name has '_part2' in it, then we have to modify it.
            if '_part2' in tmp[0]:
                offset = part_lengths[chrom + '_part1']
                newpos = str(int(tmp[1]) + offset)
            else:
                newpos = tmp[1]
            #   then print out the VCF line with the new position and chromosome
            #   name
            toprint = [chrom, newpos] + tmp[2:]
            print '\t'.join(toprint)
