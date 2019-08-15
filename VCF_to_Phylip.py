#!/usr/bin/env python
"""Convert from a VCF to a PHYLIP input format. Treats heterozygous calls as
missing data. Takes one argument:
    1) VCF"""

import sys

sequences = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        if line.startswith('##'):
            continue
        elif line.startswith('#CHROM'):
            tmp = line.strip().split()
            samples = tmp[9:]
            for s in samples:
                sequences[s] = ''
        else:
            tmp = line.strip().split()
            ref = tmp[3]
            alt = tmp[4]
            genotypes = [g.split(':')[0] for g in tmp[9:]]
            # Start building the genotype strings
            for s, g in zip(samples, genotypes):
                if g == '0/0':
                    sequences[s] += ref
                elif g == '1/1':
                    sequences[s] += alt
                else:
                    sequences[s] += 'N'
            sys.stderr.write('Added ' + tmp[2] + '\n')


# Then print out the phylip input
n_samp = str(len(sequences))
n_sites = str(len(sequences[sequences.keys()[0]]))

print '\t' + n_samp + '\t' + n_sites
for seqname in sorted(sequences):
    # Calculate how many blanks we need
    nblank = 10 - len(seqname)
    if nblank < 0:
        sname = seqname[:10]
    else:
        sname = seqname + ' '*nblank
    print sname + sequences[seqname]
