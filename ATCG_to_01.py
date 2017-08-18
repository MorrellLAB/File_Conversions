#!/usr/bin/env python

#   A script to convert ATCG calls to 0/1 calls for Structure
#   Assumes that markers are rows and individuals are columns

import sys

missing = '-999'

individuals = []
transposed_matrix = []
with open(sys.argv[1], 'r') as f:
    for index, line in enumerate(f):
        if index % 1000 == 0:
            sys.stderr.write('Read ' + str(index) + ' lines\n')
        if index==0:
            header=line.strip()
        else:
            tmp = line.strip().split('\t')
            individuals.append(tmp[0])
            transposed_matrix.append(tmp[1:])

transposed_matrix = zip(*transposed_matrix)
to_output = []
for index, line in enumerate(transposed_matrix):
    if index % 1000 == 0:
        sys.stderr.write('Processed ' + str(index) + ' lines\n')
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
    if index % 1000 == 0:
        sys.stderr.write('Wrote ' + str(index) + ' lines\n')
    print '\t'.join([individuals[index]] + list(l))
