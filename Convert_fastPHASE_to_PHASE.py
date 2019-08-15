#!/usr/bin/env python
"""Convert from the fastPHASE input to PHASE. It assumes that all loci are
SNPs (S). This is designed to work with PLINK 1.9 --recode fastphase files.
Prints the PHASE input file to stdout. Takes 1 argument:
    1) fastphase.inp file
"""

import sys


def main(fp):
    """Main function."""
    with open(fp, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                # First line is the number of individuals - print it out
                print line.strip()
            elif index == 1:
                # Second line is the number of loci. Also print it out.
                nloci = int(line.strip())
                print line.strip()
            elif index == 2:
                # Third line is the positions - print it out
                print line.strip()
                # Also print out the locus type!
                print nloci * 'S'
            else:
                # For the rest of the lines, print them out unmodified, except
                # replace the # ID ... with just #...
                if line.startswith('#'):
                    print line.strip().replace(' ID ', '')
                else:
                    print line.strip()
    return


if len(sys.argv) != 2:
    print """Convert from the fastPHASE input to PHASE. It assumes that all loci are
SNPs (S). This is designed to work with PLINK 1.9 --recode fastphase files.
Prints the PHASE input file to stdout. Takes 1 argument:
    1) fastphase.inp file"""
    exit(1)
else:
    main(sys.argv[1])
