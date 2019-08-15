#!/usr/bin/env python
"""A script to take a genotyping matrix with population assignment and produce
FASTA files for each population. This was written with K. Thornton's
libsequence tools in mind. This script will also remove monomorphic sites.

Assumes that samples are rows and markers are columns. The first column has a
population of origin and the second column is the individual ID."""

import sys
genotype_matrix = sys.argv[1]
missing = 'NA'

#   This dictionary will contain the population genetic data that we are
#   reading from the input file. It will be of the form
#   {
#       'Pop1_ID': {
#                       'PIs': [PIs],
#                       'Genotypes': [genotypes],
#                       ...
#                  }
#       'Pop2_ID': {
#                       ...
#                  }
#       ...
#   }
popdata = {}

with open(genotype_matrix, 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        else:
            tmp = line.strip().split()
            #   Get the population identifier and PI number
            popid = tmp[0]
            pi_no = tmp[1]
            #   Assign some population IDs. This is to do the combined
            #   population summary stats.
            # if popid == '1' or popid == '2':
            #    popid = '12'
            #   Stick it into the dictionary
            if popid not in popdata:
                popdata[popid] = {
                    'PIs': [],
                    'Genotypes': []
                }
            popdata[popid]['PIs'].append(pi_no)
            popdata[popid]['Genotypes'].append(tmp[2:])


#   Next, we  iterate through the dictionary, removing monomorphic markers
for pop in popdata:
    gen_mat = popdata[pop]['Genotypes']
    #   Transpose the genotype matrix so we can iterate through the markers
    gen_mat = zip(*gen_mat)
    #   And start a new matrix for filtered data
    filtered = []
    for marker in gen_mat:
        #   what are the states?
        alleles = set(marker)
        #   Discard any missing data
        alleles.discard(missing)
        #   If they are monomorphic after discarding missing data, then toss it
        if len(alleles) == 1:
            continue
        else:
            filtered.append(marker)
    #   Then, we transpose it back
    filtered = zip(*filtered)
    #   And write it out!
    handle = open('Pop_' + pop + '.fasta', 'w')
    for pi, genotypes in zip(popdata[pop]['PIs'], filtered):
        #   We need to convert the NA into N
        new_geno = ['N' if x == 'NA' else x for x in genotypes]
        towrite = '>' + pi + '\n' + ''.join(new_geno) + '\n'
        handle.write(towrite)
        handle.flush()
    handle.close()
