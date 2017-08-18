#!/usr/bin/env python
"""Convert from the HMP format from T3 into PLINK PED files. Will order markers
according to a supplied PLINK Map file. Takes two agruments:
    1) T3 hmp file
    2) PLINK MAP
 """

import sys


def store_alleles(vcf):
    """Reads a VCF and stores a dictionary of the reference and alternate
    alleles for a SNP. SNP IDs are keys, and tuples of (ref, alt) are values."""
    states = {}
    with open(vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split()
                states[tmp[2]] = (tmp[3], tmp[4])
    return states


def parse_hapmap(hmp):
    """Parse the HapMap file and store it as a nested dictionary. The first key
    is the sample name and the second key is the SNP name."""
    hap_dat = {}
    with open(hmp, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                samples = line.strip().split('\t')[11:]
                for s in samples:
                    hap_dat[s] = {}
            else:
                tmp = line.strip().split('\t')
                snpid = tmp[0]
                for samp, geno in zip(samples, tmp[11:]):
                    if geno == 'NN':
                        call = '00'
                    else:
                        call = geno
                    hap_dat[samp].update({snpid: (call[0], call[1])})
    return hap_dat


def generate_ped(ped_data, maporder):
    """Generate the PED in the proper order."""
    ped = []
    for sample in sorted(ped_data):
        fid = '0'
        pid = '0'
        mid = '0'
        sex = '0'
        pheno = '-9'
        geno = []
        for snp in maporder:
            geno += list(ped_data[sample][snp])
        ped.append([fid, sample, pid, mid, sex, pheno] + geno)
    return ped


def main(hmp, plinkmap):
    """Main function."""
    hap = parse_hapmap(hmp)
    # Get the map order
    order = []
    with open(plinkmap, 'r') as f:
        for line in f:
            order.append(line.strip().split()[1])
    pedfile = generate_ped(hap, order)
    for row in pedfile:
        sys.stdout.write('\t'.join(row) + '\n')
    return


main(sys.argv[1], sys.argv[2])
