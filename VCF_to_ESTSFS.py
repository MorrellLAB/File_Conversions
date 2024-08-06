#!/usr/bin/python3 

"""Convert a VCF and BED files with genotype information for outgroup samples 
into input for est-sfs. The primary dependency is the `cyvcf2` library for VCF 
parsing.
"""  

# Script by Jacob Pacheco to convert the output of est-sfs
# [doi:10.1534/genetics.118.301120] into inferred ancestral state.
# Should be suitable for annotation of a VCF or creating a list of
# ancestral state for a list of SNPs.

import gzip
import sys
from cyvcf2 import VCF


def est_vcf_convert(bed_file1, bed_file2, vcf_file, output_file, bed_file3=None):
    '''This script takes three tab-delimited BED files with the following fields:
    1 -chromosome
    2 -zero-based variant position
    3 -one-based variant position
    4 -SNP name (or ".")
    5 -reference allele
    6 -alternate allele
    7 -allele present at the position in an outgroup sample
    8 -dependencies VCF, gzip, sys 
    The BED files need to be provided in order so that outgroup samples most closely
    related to the focal species appear earlier (farther left in file order).
    Outputs are designed for input for ancestral state inference in EST-SFS
    and a VCF compatible file'''

    # GZIPPED VCF FILES CURRENTLY DO NOT WORK
    if vcf_file.endswith('.gz'):
        vcf_file = gzip.open(vcf_file, 'rt')
    else:
        vcf_file = open(vcf_file, 'r')

    if bed_file1.endswith('.gz'):
        bed_file1 = gzip.open(bed_file1, 'rt')
    else:
        bed_file1 = open(bed_file1, 'r')

    if bed_file2.endswith('.gz'):
        bed_file2 = gzip.open(bed_file2, 'rt')
    else:
        bed_file2 = open(bed_file2, 'r')

    if bed_file3:
        if bed_file3.endswith('.gz'):
            bed_file3 = gzip.open(bed_file3, 'rt')
        else:
            bed_file3 = open(bed_file3, 'r')

    # Read BED files
    bed_file1_lines = bed_file1.read().splitlines()
    bed_file2_lines = bed_file2.read().splitlines()
    if bed_file3:
        bed_file3_lines = bed_file3.read().splitlines()

    # Read VCF file
    vcf_reader = VCF(vcf_file)
    VCF_metadata = []

    for record in vcf_reader:
        chrom = record.CHROM if record.CHROM is not None else '.'
        pos = record.POS - 1 if record.POS is not None else '.'
        record_id = record.ID if record.ID is not None else '.'
        ref = record.REF if record.REF is not None else '.'
        alt = record.ALT[0] if record.ALT and record.ALT[0] is not None else '.'
        ns = record.INFO.get('NS', '.')
        ac = record.INFO.get('AC', '.')
        VCF_metadata.append([chrom, pos, record_id, ref, alt, ns, ac])

    # Close files
    vcf_file.close()
    bed_file1.close()
    bed_file2.close()
    if bed_file3:
        bed_file3.close()

    # Process BED and VCF data
    nucleotide_Base_binary = {"A": "1,0,0,0", "C": "0,1,0,0", "G": "0,0,1,0", "T": "0,0,0,1"}
    nucleotide_ref_binary = {"AT": [1,0,0,1], "AG": [1,0,1,0], "AC": [1,1,0,0], "CT": [0,1,0,1], "CG": [0,1,1,0], "CA": [1,1,0,0], "TG": [0,0,1,1], "TC": [0,0,1,1], "TA": [1,0,0,1], "GT": [0,0,1,1], "GC": [0,0,1,1], "GA": [1,0,1,0], 'A': [1,0,0,0], "C": [0,1,0,0], "G": [0,0,1,0], "T": [0,0,0,1]}

    with open(output_file, "w") as outfile:
        for i in range(len(bed_file1_lines)):
            bed_file1_split = bed_file1_lines[i].split()
            bed_file2_split = bed_file2_lines[i].split()
            if bed_file3:
                bed_file3_split = bed_file3_lines[i].split()

            AC = VCF_metadata[i][-1]
            if isinstance(AC, int):
                OutputAC = f'AC={int(AC)}'
                calcAC = AC
            else:
                calcAC = 0
                OutputAC = '.'

            NS = VCF_metadata[i][5]
            if isinstance(NS, int):
                OutputNS = f'NS={int(NS)}'
                calcNS = NS
            else:
                calcNS = 0
                OutputNS = '.'

            ReferencealleleFreq = ((int(calcNS)*2) - (int(calcAC))) // 2
            Alternateallelefreq = int(calcAC) // 2

            reference_allele = str(VCF_metadata[i][4])
            alternate_allele = str(VCF_metadata[i][3])
            genotypesReference = nucleotide_ref_binary.get(reference_allele, [0,0,0,0])
            genotypesReference = list(map(lambda x: x * Alternateallelefreq, genotypesReference))

            genotypesAlternate = nucleotide_ref_binary.get(alternate_allele, [0,0,0,0])
            genotypesAlternate = list(map(lambda x: x * ReferencealleleFreq, genotypesAlternate))

            genotypes = [genotypesReference[j] + genotypesAlternate[j] for j in range(len(genotypesAlternate))]
            result_str = ','.join(str(k) for k in genotypes)

            BED_file1_binary_conversion = nucleotide_Base_binary.get(bed_file1_split[-1], "0,0,0,0")
            BED_file2_binary_conversion = nucleotide_Base_binary.get(bed_file2_split[-1], "0,0,0,0")

            if bed_file3:
                BED_file3_binary_conversion = nucleotide_Base_binary.get(bed_file3_split[-1], "0,0,0,0")
                converted_file = (
                    f"{VCF_metadata[i][1]}\t{VCF_metadata[i][2]}\t{VCF_metadata[i][3]}\t"
                    f"{VCF_metadata[i][4]}\t{OutputAC}\t{OutputNS}\t"
                    f"{bed_file1_split[-1]}\t{bed_file2_split[-1]}\t{bed_file3_split[-1]}\t"
                )
                EST_SFS_VCF_compatible_file = (result_str + "\t" +
                                                BED_file1_binary_conversion + "\t" +
                                                BED_file2_binary_conversion + "\t" +
                                                BED_file3_binary_conversion + "\n")
                outfile.write(EST_SFS_VCF_compatible_file)
            else:
                converted_file = (
                    f"{VCF_metadata[i][2]}\t{VCF_metadata[i][3]}\t{VCF_metadata[i][4]}\t"
                    f"{OutputAC}\t{OutputNS}\t"
                    f"{bed_file1_split[-1]}\t{bed_file2_split[-1]}\t"
                )
                EST_SFS_VCF_compatible_file = (result_str + "\t" +
                                                BED_file1_binary_conversion + "\t" +
                                                BED_file2_binary_conversion + "\n")
                outfile.write(EST_SFS_VCF_compatible_file)


if len(sys.argv) == 6:
    bed_file1 = sys.argv[1]
    bed_file2 = sys.argv[2]
    bed_file3 = sys.argv[3]
    vcf_file = sys.argv[4]
    output_file = sys.argv[5]
    est_vcf_convert(bed_file1, bed_file2, vcf_file, output_file, bed_file3)

elif len(sys.argv) == 5:
    bed_file1 = sys.argv[1]
    bed_file2 = sys.argv[2]
    vcf_file = sys.argv[3]
    output_file = sys.argv[4]
    est_vcf_convert(bed_file1, bed_file2, vcf_file, output_file)
