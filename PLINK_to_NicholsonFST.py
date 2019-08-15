#!/usr/bin/env python
"""A script to generate input files for the Nicholson 2002 FST estimator, as
implemented in the R `popgen' package. This assumes PLINK input file formats.
It acceps a .PED and a .CLST file for genotyping information and cluster
membership. """

import argparse


def parse_args():
    """A function to parse arguments."""
    parser = argparse.ArgumentParser(
        description=('Generate input files for Nicholson 2002 FST estimation '
                     'from PLINK input files. Accepts .PED and .CLST files '
                     'and outputs .NUMA, .N, and .X files.'),
        add_help=True)
    parser.add_argument(
        '--ped-file',
        '-p',
        required=True,
        help='PED file')
    parser.add_argument(
        '--clst-file',
        '-c',
        required=True,
        help='CLST file')
    parser.add_argument(
        '--missing-value',
        '-m',
        required=False,
        help='Missing data value. Defaults to 0',
        default='0')
    parser.add_argument(
        '--output-prefix',
        '-o',
        required=False,
        help='Prefix for output files. Defaults to FST',
        default='FST')
    args = parser.parse_args()
    return args


def summarize_data(clst, ped, miss):
    """A function to summarize the data in the given .CLST and .PED files.
    Will return the number of loci, and the number of populations."""
    #   Empty dictionary to hold the data to summarize
    assignment = {}
    #   First start to read through the cluster file
    with open(clst, 'r') as fhandle:
        for line in fhandle:
            #   Split it on whitespace
            clst_data = line.strip().split()
            #   The two we are interested in are the individual ID, which is
            #   the first element, and the cluster ID, which is the last
            #   element.
            ind_id = clst_data[0]
            clst_id = clst_data[-1]
            assignment[ind_id] = clst_id
    #   What are the cluster IDs?
    clst_ids = sorted(list(set(assignment.values())))
    #   Now, we count how many loci and the alleles at each locus. We will also
    gen_mat = []
    with open(ped, 'r') as fhandle:
        for index, line in enumerate(fhandle):
            #   split up the line on tabs
            linedata = line.strip().split()
            if index == 0:
                #   how many loci are there? We put this after an if, since
                #   we only really need to calculate this once.
                #   PLINK .PED files have six fields of non-genotype data, and
                #   also have diploid genotypes.
                nloci = len(linedata[6:])/2
            gen_mat.append(linedata[6:])
    #   We zip the matrix to transpose it
    gen_mat = zip(*gen_mat)
    #   Then, iterate through the transposed matrix to get the alleles at each
    #   locus. Remeber it is diploid. Instead of iterating through the matrix
    #   itself, we iterate through the range of integers from 0 to nloci.
    alleles = []
    for index in xrange(0, nloci):
        #   Get the diploid calls
        strand1 = ''.join(gen_mat[2*index])
        strand2 = ''.join(gen_mat[(2*index) + 1])
        #   Create a set out of them to count the number of alleles
        all_alleles = set(strand1 + strand2)
        #   Then remove the missing data value, if it's there
        all_alleles.discard(miss)
        #   Then sort the alleles (so we can make sure they are always in the
        #   same order) and return them
        alleles.append(sorted(list(all_alleles)))

    #   Return all the data
    return(assignment, clst_ids, nloci, alleles)


def separate_clusters(pops, ped):
    """A function to separate the genotyping data into separate matrices by
    cluster ID as listed in the cluster file. This function takes two lists
    and returns a dictionary of lists."""
    #   Build the data structure to return
    #   This will be of the form:
    #   {
    #       'pop1_ID': [Genotype matrix for Pop 1],
    #       'pop2_ID': [Genotype matrix for Pop 2],
    #        ...
    #       'popN_ID': [Genotype matrix for Pop 3]
    #   }
    pop_data = {}
    for popid in pops.values():
        pop_data[popid] = []
    #   Next, we iterate through the ped file again, and separate the calls
    #   based on population
    with open(ped, 'r') as fhandle:
        for line in fhandle:
            #   We split on tabs
            ind_data = line.strip().split()
            #   What is the individual ID? It is the second column
            ind_id = ind_data[1]
            #   Based on that, what is the population of origin?
            pop_orig = pops[ind_id]
            #   Then, tack that data onto the appropriate list in pop_data
            #   We just care about the genotypes.
            pop_data[pop_orig].append(ind_data[6:])
    return pop_data


def do_counts(nloci, alleles, miss, population_genotypes):
    """A function to iterate through the population-separated dictionary and
    return the counts for the Nicholson estimator. This will count the number
    of genotypes at each locus for each population, and the number of each
    allele at each locus in each population."""
    #   Start a new dictionary for our counts. We can count alleles right away.
    #   Basically, it will look like this:
    #   {
    #       'NumAlleles': [2, 2, 2, 2, ... N_Loci],
    #       'NumGenotypes': {
    #                           'Pop1_ID': [200, 200, 200, ..., N_Loci],
    #                           'Pop2_ID': [200, 200, 200, ..., N_Loci],
    #                           ...
    #                           'PopN_ID': [200, 200, 200, ..., N_Loci]
    #                       }
    #       'PopAlleles': {
    #                           'Pop1_ID': {
    #                                           [100, ..., N_Loci],
    #                                           [100, ..., N_Loci',
    #                                           ...
    #                                           [100, ..., N_Loci]
    #                                      },
    #                           'Pop2_ID': {
    #                                           ...
    #                                      },
    #                           ...
    #                           'PopN_ID': {
    #                                           ...
    #                                      }
    #                      }
    #   }
    counts = {
        'NumAlleles': [str(len(a)) for a in alleles],
        'NumGenotypes': {},
        'PopAlleles': {}
    }
    #   Start iterating through the dictionary and performing the counts
    for popid, genotypes in population_genotypes.iteritems():
        #   Start sub-lists for our data, since they will ultimately need to
        #   be matrices.
        counts['NumGenotypes'][popid] = []
        counts['PopAlleles'][popid] = []
        #   Then, iterate through the transposed genotypes and count up the
        #   relevant information
        g_trans = zip(*genotypes)
        for index in xrange(0, nloci):
            #   Start counting the number of genotypes
            n_geno = 0
            #   Start a dict of the counts of inidividual alleles
            p_allele = {}
            for a in alleles[index]:
                p_allele[a] = 0
            strand1 = g_trans[2*index]
            strand2 = g_trans[(2*index) + 1]
            #   Zip them and start counting
            for call in zip(strand1, strand2):
                if miss in call:
                    continue
                else:
                    #   No missing data -> add one to n_geno
                    n_geno += 1
                    #   Count up the number of alleles in the pop.
                    for a in alleles[index]:
                        if a in call:
                            p_allele[a] += 1
            #   Then, tack the appropriate information onto the population
            #   specific list
            counts['NumGenotypes'][popid].append(str(n_geno))
            popallele = []
            for a in p_allele:
                popallele.append(str(p_allele[a]))
            counts['PopAlleles'][popid].append(popallele)
    #   Return the big dictionary
    return counts


def write_counts(clusters, numloci, nicholson_data, prefix):
    """This function just writes the matrices into the correct files."""
    #   Start the handles for writing data
    numa_handle = open(prefix + '_NUMA.txt', 'w')
    N_handle = open(prefix + '_N.txt', 'w')
    X_handle = open(prefix + '_X.txt', 'w')
    #   Get the NUMA out of the way
    numa_handle.write('\n'.join(nicholson_data['NumAlleles']))
    numa_handle.flush()
    numa_handle.close()
    #   Then iterate through the data and start writing things!
    for index in xrange(0, numloci):
        N_towrite = []
        X_towrite = []
        for c in clusters:
            N_towrite.append(nicholson_data['NumGenotypes'][c][index])
            X_singlepop = []
            for a in xrange(0, int(nicholson_data['NumAlleles'][index])):
                X_singlepop.append(nicholson_data['PopAlleles'][c][index][a])
            X_towrite += X_singlepop
        N_handle.write('\t'.join(N_towrite) + '\n')
        X_handle.write('\t'.join(X_towrite) + '\n')
    #   Flush those handles and close them
    N_handle.flush()
    N_handle.close()
    X_handle.flush()
    X_handle.close()
    return


def main():
    """Main function."""
    arguments = vars(parse_args())
    ped_file = arguments['ped_file']
    clst_file = arguments['clst_file']
    missing_value = arguments['missing_value']
    output_prefix = arguments['output_prefix']
    #   What are the populations and the number of loci?
    assignment, clusters, n_loci, alleles = summarize_data(
        clst_file,
        ped_file,
        missing_value)
    #   Next, we separate the data matrix by population
    population_genotypes = separate_clusters(assignment, ped_file)
    count_data = do_counts(
        n_loci,
        alleles,
        missing_value,
        population_genotypes)
    write_counts(clusters, n_loci, count_data, output_prefix)
    return


main()
