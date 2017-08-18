#!/usr/bin/evn python
#   Script to convert a Hudson-like polymorphism table to a FASTA file for
#   input into 'compute' (libsequence)
#   The input file format does not have a typical Htable header, and is a 
#   matrix of genotyping data, so we also have to filter on monomorphic
#   markers.

import sys

samples = []
genotypes = []
with open(sys.argv[1], 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        else:
            tmp = line.strip().split('\t')
            samples.append(tmp[0])
            #   Replace the diploid genotypes with haploid ones
            hap_gt = [list(set(a))[0] if len(set(a))==1 else 'N' for a in tmp[1:]]
            #   Then, replaces the 'B' with 'C'
            hap_gt_nuc = ['C' if a=='B' else a for a in hap_gt]
            genotypes.append(hap_gt_nuc)

#   next we transpose the matrix so we can iterate over markers instead of
#   individuals.
genotypes = zip(*genotypes)

newmarkers = []
for marker in genotypes:
    #   cast it to a set, which gives just unique values
    alleles = set(marker)
    #   if there are missing calls, ignore them for now
    #   we use set.discard() because it won't throw a KeyError if there is no N
    alleles.discard('N')
    #   If the remaining set is only one item long, then the marker is 
    #   monomorphic in the population and we remove it
    if len(alleles) == 1:
        continue
    else:
        newmarkers.append(marker)

#   transpose it back so we can tack all the individuals together
newmarkers = zip(*newmarkers)

#   next, we print out the FASTA file
for individual in zip(samples, newmarkers):
    print '>' + individual[0]
    print ''.join(list(individual[1]))
