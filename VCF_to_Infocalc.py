#!/usr/bin/env python
"""Convert a VCF to the Infocalc (Rosenberg 2003) input file format. Takes two
arguments:
    1) VCF (gzipped)
    2) Population-sample map, two-column whitespace delimited text file with
       sample names in column 1 and population in column 2. No header.
"""

import sys
import gzip

try:
    gzvcf = sys.argv[1]
    popmap = sys.argv[2]
except IndexError:
    sys.stderr.write(__doc__)
    exit(1)


def parse_map(m):
    """Return a dictonary of sample:pop values."""
    pmap = {}
    with open(m, 'r') as f:
        for line in f:
            s, p = line.strip().split()
            pmap[s] = p
    return pmap


def parse_vcf(v):
    """Parse a VCF and return a list of locus names and a dictionary of variant
    data. The variant data is a dictionary of tuple of lists:
    {
        sample1: ([L1 A1, L2 A1, L3 A1], [L1 A2, L2 A2, L3 A2]),
        sample2: ([ . . .], [ . . .]),
        . . .
    }
    Phase is not preserved for infocalc."""
    loci = []
    variants = {}
    with gzip.open(v, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                continue
            elif line.startswith('#CHROM'):
                header = line.strip().split('\t')
                samples = header[9:]
                for s in samples:
                    variants[s] = ([], [])
            else:
                tmp = line.strip().split('\t')
                loci.append(tmp[2])
                for s, g in zip(samples, tmp[9:]):
                    geno = g.split(':')[0]
                    # Any integer 0 or less is considered missing data, but we
                    # use -9 to be consistent with STRUCTURE convention
                    if geno == '.':
                        a1 = '-9'
                        a2 = '-9'
                    else:
                        # 0 is not a valid allele for Infocalc, so we replace
                        # 1 with 2 and 0 with 1.
                        valid_alleles = geno.replace('1', '2').replace('0', '1')
                        a1, a2 = valid_alleles.split('/')
                    variants[s][0].append(a1)
                    variants[s][1].append(a2)
    return(loci, variants)


def validate_map(m, v):
    """Make sure that the sample names in the VCF and the sample names in the
    population map are identical."""
    vcf_samples = list(v.keys())
    pmap_samples = list(m.keys())
    # A bit ugly, but if the sorted samples lists are not exactly equal, then
    # we throw an error
    if sorted(vcf_samples) == sorted(pmap_samples):
        return True
    else:
        return False


def main(v, m):
    """Main function."""
    # Parse the population map
    pmap = parse_map(m)
    # And the VCF
    loci, variants = parse_vcf(v)
    # Check that the map and vcf have the same samples
    if not validate_map(pmap, variants):
        sys.stderr.write(
            'Error: The population map and VCF do not have the same samples!\n'
            )
        exit(2)
    # Start printing the output file. The first line is the list of loci
    sys.stdout.write(' '.join(loci) + '\n')
    # Then write all of the samples, with two lines per sample
    for i, s in enumerate(sorted(variants)):
        # Assign a dummy numeric sample ID
        iid = str(i+1)
        # We are being extra catuious with this check:
        if s not in pmap:
            sys.stderr.write(
                'Error: sample ' + s + ' not found in pop map!\n')
            exit(2)
        samp_pop = pmap[s]
        # Then for each locus:
        for loc in variants[s]:
            # Build the string to print. The first five columns are
            # non-genotype data, and the third column is the population
            # identifier, which is the default column used by infocalc.
            toprint = [s, iid, samp_pop, 'NA', 'NA'] + loc
            sys.stdout.write(' '.join(toprint) + '\n')


# Run the function
main(gzvcf, popmap)
