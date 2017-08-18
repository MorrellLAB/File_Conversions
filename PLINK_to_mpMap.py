#!/usr/bin/env python
"""Converts a parental .ped, offspring .ped, and a .map file into the necessary
input files for mpMap. Pedigree has to be created separately. Note that mpMap
only takes 2, 4, or 8 founder lines. Therefore, the ped files that are fed to
this script should be split by family."""

import sys
import os


def atcg_to_numeric(ped, snpmap, missingval='0'):
    """Creates a dictionary that associates the major allele of a SNP with 1,
    and the minor allele with 0. Major and minor allele states are determined
    by frequency in the PED file. Maps missing and heterozygous calls to NA."""

    def major_minor(calls, miss=missingval):
        """Helper function to return major state and minor state. They are
        arbitrary if the allele frequencies are equal. Missing calls are not
        counted for frequency calcuations."""
        alleles = list(set(calls))
        #   Remove missing calls
        nomiss = [c for c in calls if c != miss]
        #   Count up how many of each allele we see
        counts = {}
        for a in alleles:
            counts[a] = nomiss.count(a)
        #   Get the minimum
        minor_allele = min(counts, key=counts.get)
        #   The major allele is not the minor allele.
        major_allele = [a for a in counts.keys() if a != minor_allele][0]
        return (major_allele, minor_allele)

    snp_ids = []
    with open(snpmap, 'r') as f:
        for line in f:
            tmp = line.strip().split('\t')
            snp_ids.append(tmp[1])
    #   Then, step through the genotype matrix to find the major and minor
    #   alleles for each SNP. They are in the same order as they appear in the
    #   map file.
    g_mat = []
    with open(ped, 'r') as f:
        for line in f:
            #   The first five columns are not genotype data.
            g_mat.append(line.strip().split('\t')[6:])
    #   Then zip it to transpose the matrix, so we can iterate over columns
    g_mat_t = zip(*g_mat)
    #   Have to iterate in steps of 2, since the PED has diploid calls
    snp_numeric = {}
    for snpid, col in zip(snp_ids, xrange(0, len(g_mat_t), 2)):
        all_calls = g_mat_t[col] + g_mat_t[col+1]
        minor, major = major_minor(all_calls)
        #   The ATCG to 01 association will be held in a dictionary so we can
        #   do simple lookups
        snp_numeric[snpid] = {
            minor + minor: '0',
            major + major: '1'
            }
    return snp_numeric


def create_mpmap_mat(ped, snpmap, num, cycle, par=None):
    """Creates a genotyping matrix out of a PED file and the SNP map file.
    Heterozygous sites are treated as missing data. 'par' is an optional
    argument that denotes the parents. If it is supplied, then only the
    specified samples will be printed."""
    snp_order = []
    with open(snpmap, 'r') as f:
        for line in f:
            snp_order.append(line.strip().split('\t')[1])
    g_mat = []
    #   Start the samples list off with an empty string, since in the mpMap docs
    #   the first entry in this matrix should be empty.
    samples = ['']
    with open(ped, 'r') as f:
        for line in f:
            tmp = line.strip().split('\t')
            if not par:
                samples.append(tmp[1])
                g_mat.append(tmp[6:])
            else:
                if tmp[1] in par:
                    samples.append(tmp[1])
                    g_mat.append(tmp[6:])
                else:
                    continue
    g_mat_t = zip(*g_mat)
    founder_mat = []
    for snpid, col in zip(snp_order, xrange(0, len(g_mat_t), 2)):
        calls = []
        for c in zip(g_mat_t[col], g_mat_t[col+1]):
            calls.append(num[snpid].get(c[0]+c[1], 'NA'))
        founder_mat.append([snpid] + calls)
    #   Then put the samples list onto the front of the founder matrix
    founder_mat = [samples] + founder_mat
    return founder_mat


def create_map(snpmap):
    """Re-organizes the PLINK map format into the mpMap format."""
    m_map = [['SNPID', 'Chromosome', 'GeneticPosition']]
    with open(snpmap, 'r') as f:
        for line in f:
            tmp = line.strip().split('\t')
            #   The order of the columns for mpMap are
            #       SNPID, Chromosome, GeneticPosition
            m_map.append([tmp[1], tmp[0], tmp[2]])
    return m_map


def get_parents(offspring):
    """Gets the parental IDs from the offspring PED. These will be written to
    the "founder" matrix. We just read the first line since the PED files
    should all be split by family already."""
    with open(offspring, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            return tmp[2:4]


def main(parents, offspring, snpmap):
    """Main function."""
    snp_num = atcg_to_numeric(parents, snpmap)
    par = get_parents(offspring)
    founders = create_mpmap_mat(parents, snpmap, snp_num, 0, par=par)
    final = create_mpmap_mat(offspring, snpmap, snp_num, 1)
    m_map = create_map(snpmap)
    #   Write the output to tab-separated files
    founder_fname = os.path.basename(parents).replace('.ped', '_mpMap.txt')
    handle = open(founder_fname, 'w')
    for row in founders:
        handle.write('\t'.join(row) + '\n')
    handle.close()
    final_fname = os.path.basename(offspring).replace('.ped', '_mpMap.txt')
    handle = open(final_fname, 'w')
    for row in final:
        handle.write('\t'.join(row) + '\n')
    handle.close()
    map_fname = os.path.basename(snpmap).replace('.map', '_mpMap.map')
    handle = open(map_fname, 'w')
    for row in m_map:
        handle.write('\t'.join(row) + '\n')
    handle.close()
    return

if len(sys.argv) != 4:
    print("""
Usage:
    PLINK_to_mpMap.py [Parents.ped] [Progeny.ped] [SNPs.map]

Will create 'founders' and 'final' genotype matrices and a SNP map for input
into mpMap. The pedigree will have to be created separately.""")
    exit(1)
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3])
