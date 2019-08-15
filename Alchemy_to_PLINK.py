#!/usr/bin/env python
"""Convert from ALCHEMY report to PLINK PED files. Translate the names from the
old names to the BOPAC names, and put them in physical map order. Check the
alleles and fix as necessary. Takes four arguments:
    1) Alchemy report
    2) BOPA Physical positions VCF
    3) Genetic map CSV
    4) SNP Name Translation CSV
    5) SNP State Translation CSV"""

import sys

# define a dictionary for reverse complement lookups
REVCOMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

def parse_vcf(vcf):
    """Read the VCF, store the SNPs in physical order."""
    vcf_data = {}
    bopa_alleles = {}
    with open(vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split('\t')
                snpid = tmp[2]
                ref = tmp[3]
                alt = tmp[4]
                pos = int(tmp[1])
                chrom = tmp[0]
                # Store just the alleles for checking flip errors
                bopa_alleles[snpid] = (ref, alt)
                if chrom not in vcf_data:
                    vcf_data[chrom] = [(snpid, pos, ref, alt)]
                else:
                    vcf_data[chrom].append((snpid, pos, ref, alt))
    # Sort it by position
    for c in vcf_data:
        vcf_data[c] = sorted(vcf_data[c], key=lambda x: x[1])
    return (bopa_alleles, vcf_data)


def parse_translation(name_trans):
    """Store the translation table as a dictionary."""
    trans = {}
    with open(name_trans, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split(',')
                trans[tmp[3]] = tmp[2]
    return trans


def parse_snpstate(snp_trans):
    """Store the translation table as a dictionary, and nucleotide state of A & B as a tuple."""
    s_trans = {}
    with open(snp_trans, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split(',')
                s_trans[tmp[0]] = {'A': tmp[1], 'B': tmp[2], '0': '0'}
    return s_trans


def parse_alchemy(trans, alchemy):
    """Parse the ALCHEMY report, and store the genotyping data in a dictionary
    of dictionaries. The first key is the individual ID, and the second key is
    the SNP ID, translated to BOPAC."""
    missing = '0'
    thresh = 0.7
    alc_data = {}
    with open(alchemy, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split('\t')
                if tmp[0] in trans:
                    snpid = trans[tmp[0]]
                else:
                    snpid = tmp[0]
                sample = tmp[1]
                ab_call = tmp[2]
                nuc_call = tmp[3]
                prob = float(tmp[4])
                # Check the sample ID and the probability. If the sample ID is
                # 'Blank', we skip it
                if sample == 'Blank':
                    continue
                else:
                    # Fix the sample names so that they match up with other
                    # data sources. Cycle 1 names are MS10S30XXX, Cycle 2 names
                    # are MS11S2XXXX, and Cycle 3 names are MS12_XXXX.
                    if sample.startswith('G'):
                        sample = sample.replace('G10W', 'MS10S3')
                        sample = sample.replace('-', '-0')
                    elif sample.startswith('MS11'):
                        sample = sample
                    elif sample.startswith('MS12'):
                        sample = sample
                    # Then check the probability. If it is less than our
                    # threshold (0.7) we append missing. Else, we append the
                    # call
                    if prob < thresh:
                        if sample in alc_data:
                            alc_data[sample].update({snpid: (missing, missing)})
                        else:
                            alc_data[sample] = {snpid: (missing, missing)}
                    else:
                        #calls = (nuc_call[0], nuc_call[1])
                        calls = (ab_call[0], ab_call[1])
                        if sample in alc_data:
                            alc_data[sample].update({snpid: calls})
                        else:
                            alc_data[sample] = {snpid: calls}
    return alc_data


def fix_alleles(ped_data, bopa_data):
    """Check the alleles in the PED data against the BOPA alleles."""
    fixed_ped = {}
    # First, iterate through the samples
    for samp in ped_data:
        fixed_ped[samp] = {}
        # Then through the SNPs
        for snp in ped_data[samp]:
            alc_1 = ped_data[samp][snp][0]
            alc_2 = ped_data[samp][snp][1]
            # If the BOPA SNP does not have a physical position (and thus, no
            # reference and alternate allele calls), then we can't do anything
            # and just append it...
            if snp not in bopa_data:
                fixed_ped[samp].update({snp: (alc_1, alc_2)})
            elif alc_1 == '0' and alc_2 == '0':
                # If tehy are both missing, then we just append missing calls
                fixed_ped[samp].update({snp: ('0', '0')})
            elif alc_1 in bopa_data[snp] and alc_2 in bopa_data[snp]:
                # If the alleles are in the same orientation as the VCF data,
                # then we're all good.
                fixed_ped[samp].update({snp: (alc_1, alc_2)})
            elif REVCOMP[alc_1] == alc_2:
                # If the two alleles are reverse complements, i.e., A/T or C/G,
                # then we just keep them as-is, since those are hard to detect
                # as strand-flip errors.
                fixed_ped[samp].update({snp: (alc_1, alc_2)})
            elif REVCOMP[alc_1] in bopa_data[snp] and REVCOMP[alc_2] in bopa_data[snp]:
                # Else, if the ALCHEMY calls are RCed with respect to the
                # reference genome, we flip them
                fixed_ped[samp].update({snp: (REVCOMP[alc_1], REVCOMP[alc_2])})
            else:
                # If they cannot be resolved by any of these, we set them to
                # missing.
                fixed_ped[samp].update({snp: ('0', '0')})
    return fixed_ped


def order_snps(ped_data, p_map):
    """Generate the order that the SNPs should go in by iterating through the
    chromosomes of the physical map."""
    snp_order = []
    for chrom in sorted(p_map):
        for snp in p_map[chrom]:
            if snp[0] in ped_data[ped_data.keys()[0]]:
                snp_order.append(snp[0])
    return snp_order


def generate_ped(g_mat, map_order, translation):
    """Generate the PED data with the genotype matrix and the map order of
    the SNPs."""
    ped_file = []
    for sample in sorted(g_mat):
        # Use missing codes for the family ID and the maternal and paternal IDs
        # these will be filled in later.
        fid = '-9'
        mid = '-9'
        pid = '-9'
        sex = '0'
        pheno = '-9'
        geno = []
        for marker in map_order:
            atcg_calls = [
                translation[marker][c]
                for c
                in g_mat[sample][marker]]
            geno += atcg_calls
        ped_file.append([fid, sample, pid, mid, sex, pheno] + geno)
    return ped_file


def generate_map(map_order, genetic_map, physical_map):
    """Generate the PLINK map file with the SNPs in genetic map order."""
    gmap = {}
    with open(genetic_map, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split(',')
                gmap[tmp[0]] = (tmp[1], tmp[2])
    map_file = []
    for snp in map_order:
        if snp not in gmap:
            g_pos = '-9'
            # If the SNP isn't genetically mapped, maybe it's physically mapped.
            # this is not ideal, but it it's something that works ...
            for chrom in physical_map:
                s = [x[0] for x in physical_map[chrom]]
                if snp in s:
                    break
        else:
            chrom = 'chr' + gmap[snp][0]
            g_pos = gmap[snp][1]
        for chromsnp in physical_map[chrom]:
            if snp in chromsnp:
                p_pos = str(chromsnp[1])
                break
        else:
            p_pos = '-9'
        map_file.append([chrom, snp, g_pos, p_pos])
    return map_file


def main(alchemy, bopa, genetic_map, name_trans, snp_trans):
    """Main function."""
    bopa_alleles, phys_map = parse_vcf(bopa)
    name_key = parse_translation(name_trans)
    geno_ab_states = parse_snpstate(snp_trans)
    alchemy_calls = parse_alchemy(name_key, alchemy)\
    #fixed_calls = fix_alleles(alchemy_calls, bopa_alleles)
    #ordered = order_snps(fixed_calls, phys_map)
    ordered = order_snps(alchemy_calls, phys_map)
    #plinkped = generate_ped(fixed_calls, ordered)
    plinkped = generate_ped(alchemy_calls, ordered, geno_ab_states)
    plinkmap = generate_map(ordered, genetic_map, phys_map)
    # Then write out the data
    for row in plinkped:
        sys.stdout.write('\t'.join(row) + '\n')
    for row in plinkmap:
        sys.stderr.write('\t'.join(row) + '\n')
    return


main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
