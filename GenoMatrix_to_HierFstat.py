#!/usr/bin/env python

#   A script to convert ATCG calls to 0/1 calls for Hierfstat
#   Assumes that columns are markers and rows are individuals

import sys

missing = 'U'

#   Transpose the matrix so that we can iterate over markers instead of
#   individuals. This lets us remove missing data and convert markers into
#   numeric genotypes by marker.
individuals = []
transposed_matrix = []
with open(sys.argv[1], 'r') as f:
    for index, line in enumerate(f):
        if index==0:
            header=line.strip()
        else:
            tmp = line.strip().split('\t')
            individuals.append(tmp[1])
            #   Get the population ID
            popid = tmp[0]
            #   Start at the third column, since the first two columns don't
            #   actually have any data. (0=pop, 1=PI number)
            transposed_matrix.append(tmp[2:])

transposed_matrix = zip(*transposed_matrix)
to_output = []
for index, line in enumerate(transposed_matrix):
    alleles = set(line)
    if missing in alleles:
        alleles.remove(missing)
    alleles = list(alleles)
    numerical= []
    for a in line:
        if a == missing:
            numerical.append('NA')
        else:
            numerical.append(str(alleles.index(a)))
    to_output.append(numerical)

#   print the header
print header
#   Transpose it
to_output = zip(*to_output)
for index, l in enumerate(to_output):
    print '\t'.join([popid, individuals[index]] + list(l))
