#!/usr/bin/env python
"""A script to trim down a Hudson SNP Table to the desired samples. Also has 
the option to reorder samples for input into Thornton's libsequence programs.
Takes two files as input: the table, and a file with sample names, one name 
per line. Comments (#) are interpreted as blank lines.

"""

import sys
import os
import argparse
import re

usage = """{0} TABLE LIST [-r]

Will produce a new Hudson table from TABLE, containing only samples in LIST.
SNPs that are monomorphic in the new table will be removed. TABLE is a 
Hudson SNP table, and LIST is a file containing a list of sample names, 
with one sample per line. Comments in the list file are allowed, and samples
listed in LIST that are not present in TABLE will be ignored.

If the -r flag is given, then instead, {0} reorders the 
Hudson table, listing the samples in LIST first in the table. The other 
samples appear in an arbitrary order. This is for generating paritions of the 
data, which go on to programs such as K Thornton's sharedPoly. If -r is given
and LIST is empty (or all comments) then the new table will contain the same
information, just in an arbitrary order.
""".format(os.path.basename(sys.argv[0]))

#   Give our usage/description, and automatically add -h switch
Arguments = argparse.ArgumentParser(usage=usage, add_help=True)
#   An argument for the -r switch
Arguments.add_argument('-r',
                '--reorder',
                action='store_true',
                default=False,
                help='Reorder samples?')
#   An argument for the Hudson table to modify
Arguments.add_argument('table',
                metavar='TABLE',
                type=file,
                help='The Hudson table to process')
#   An argument for the listfile to filter on
Arguments.add_argument('listfile',
                metavar='LIST',
                type=file,
                help='The list of sample names')
#   Parse the arguments, and output help and usage info if necessary
ParsedArgs = Arguments.parse_args()

htable = sys.argv[1]
lfile = sys.argv[2]

#   try/except to see if our files are readable or not. open() raises IOError
#   if the file cannot be read for some reason
try:
    open(htable, 'r')
except IOError:
    print htable + " does not exist, or is not readable!"
    exit(1)

try:
    open(lfile, 'r')
except IOError:
    print lfile + " does not exist, or is not readable!"
    exit(1)

#   A function to read and process the list of sample names
def read_list_file(fname, r):
    #   Returns a list of non-commented sample names
    samples = []
    with open(fname, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                samples.append(line.strip())
    #   It doesn't make sense to have less than 2 samples...
    #   but only if we are not reordering
    if len(samples) < 2 and not r:
        sys.stderr.write('You need to have more than one sample in your list file!\n')
        exit(1)
    return(samples)

#   A function to read and parse the Hudson table supplied as an argument
def read_hudson_table(fname):
    #   Returns a dictionary of the form
    #   {   'sample': [marker states ... ],
    #       'sample': [marker states ... ],
    #       ...
    #   }
    #   Returns a list of marker names
    g_data = {}
    with open(fname, 'r') as f:
        #   We enumerate() over the file so we can get line numbers
        for index, line in enumerate(f):
            #   Just skip the first line.
            if index == 0:
                continue
            #   The second line has marker information
            elif index == 1:
                markers = line.strip().split('\t')
            #   And then, we read the rest of the file
            else:
                #   Strip whitespace, split on tabs
                tmp = line.strip().split('\t')
                sample_name = tmp[0]
                genotypes = tmp[1:]
                #   Check for bad characters with a regex, case-insensitive,
                #   but only if the sample is not the ancestral string
                g_string = ''.join(genotypes)
                if sample_name != 'anc':
                    r = re.search('[^ATCGN]', g_string, re.IGNORECASE)
                    #   If r is defined, then we have a match!
                    if r:
                        sys.stderr.write('Illegal character '+r.group()[0]+' in sample '+sample_name+' at marker number '+str(r.start()+1)+'!\n')
                        exit(1)
                g_data[sample_name] = genotypes
    return(g_data, markers)

#    A function to trim the Hudson table
def trim_hudson(ht, flt, markers):
    #   args: Hudson table data, list of desired samples
    #   Returns a completed Hudson table
    #   A list for the sample names
    samples = []
    #   A list for the genotype matrix
    genotypes = []
    #   We will add the markers to the top of the list
    genotypes.append(markers)
    #   We will put the ancestral line in first, since we know it exists,
    #   and we need it to be at the top
    samples.append('anc')
    genotypes.append(ht['anc'])
    #   iterate through the old Hudson table data, saving only samples we are
    #   interested in
    for ind, geno in ht.iteritems():
        if ind in flt:
            samples.append(ind)
            genotypes.append(geno)
    #   New list for trimmed markers
    trimmed_markers = []
    #   Now, we transpose the genotype matrix, and remove monomorphic SNPs
    t_genotypes = zip(*genotypes)
    for m in t_genotypes:
        #   The genotype data excludes the first two items, since they are the
        #   marker names and the ancestral states
        states = set(m[2:])
        #   We remove missing data, if it is present
        states.discard('N')
        #   If the length of states is less than 2, then we have a monomorphic
        #   marker, so we drop it
        if len(states) < 2:
            continue
        else:
            trimmed_markers.append(m)
    #   Transpose our trimmed markers, and use it as the new matrix
    genotypes = zip(*trimmed_markers)
    #   The number of samples is len(samples)-1, since we remove ancestral
    nsam = len(samples)-1
    #   The number of markers is len(genotypes[0])
    nmarkers = len(genotypes[0])
    #   start building our completed Hudson table!
    completed_table = []
    #   The number of samples, and the number of markers
    completed_table.append(str(nsam)+'\t'+str(nmarkers))
    #   The marker names/positions
    #   We have to do it this way because the Hudson table has positions, and not names
    completed_table.append('\t'+'\t'.join([str(i) for i in range(1,nmarkers+1)]))
    #   Iterate through the zipped samples and marker genotypes
    for s, g in  zip(samples, genotypes[1:]):
        completed_table.append('\t'.join([s] + list(g)))
    return(completed_table)

#   A function to reorder the Hudson table
def reorder_hudson(ht, flt, markers):
    #   parsed table, filter list, markers
    #   Returns a completed Hudson table
    #   A list for the sample names
    samples = []
    #   A list for the genotype matrix
    genotypes = []
    #   We will add the markers to the top of the list
    genotypes.append(markers)
    #   We will put the ancestral line in first, since we know it exists,
    #   and we need it to be at the top
    samples.append('anc')
    genotypes.append(ht['anc'])
    #   Iterate through the list of filter names, appending those we want
    for s in flt:
        if s in ht:
            samples.append(s)
            genotypes.append(ht[s])
    #   iterate through the old Hudson table data, saving the samples
    #   that aren't in the filter list
    #   We also exclude the ancestral state here
    for ind, geno in ht.iteritems():
        if ind not in flt and ind != 'anc':
            samples.append(ind)
            genotypes.append(geno)
    #   The number of samples is len(samples)-1, since we remove ancestral
    nsam = len(samples)-1
    #   The number of markers is len(genotypes[0])
    nmarkers = len(genotypes[0])
    #   start building our completed Hudson table!
    completed_table = []
    #   The number of samples, and the number of markers
    completed_table.append(str(nsam)+'\t'+str(nmarkers))
    #   The marker names/positions
    completed_table.append('\t'+'\t'.join([str(i) for i in range(1,nmarkers+1)]))
    #   Iterate through the zipped samples and marker genotypes
    for s, g in  zip(samples, genotypes[1:]):
        completed_table.append('\t'.join([s] + g))
    return(completed_table)

#   Now we do the work!
parsed_table = read_hudson_table(htable)
filter_list = read_list_file(lfile, ParsedArgs.reorder)
if ParsedArgs.reorder:
    new_table = reorder_hudson(parsed_table[0], filter_list, parsed_table[1])
else:
    new_table = trim_hudson(parsed_table[0], filter_list, parsed_table[1])
sys.stdout.write('\n'.join(new_table)+'\n')
