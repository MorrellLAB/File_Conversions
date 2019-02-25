#!/usr/bin/env python

#by Li Lei, 2018/02/22, in St.Paul
#   A script to convert the vcf file into the input file to XP_CLR. There are three output files: genotype files, and map files:
#usage: python3 vcf2xpclr.py yourSNP.vcf yourSNP.geno yourSNP.map
import sys

#   Start iterating through the file
with open(sys.argv[1], 'r') as f:
    for line in f:
        #   ignore header lines
        if line.startswith('#'):
            continue
        #   This defines how many samples in the VCF
        else:
            tmp = line.strip().split('\t')
            #   Parse out the relevant information
            chromosome = tmp[0]
            bp_pos = tmp[1]
            SNP_id = tmp[2]
            ref_allele = tmp[3]
            alt_alleles = tmp[4]
            genotypes = tmp[9:]
            #print(SNP_id, chromosome,0,bp_pos,ref_allele,alt_alleles,file=open(sys.argv[2], "a"))
            print('\t'.join([SNP_id, chromosome,'0',bp_pos,ref_allele,alt_alleles]),file=open(sys.argv[3], "a"))
            #INFO = tmp[7].split(';')
            #AN = INFO[2].split('=')
            #chr_nb = int(AN[1])/2
            g_column = []
            for g in genotypes:
                if g == '.':
                    g_column.append('9')
                    g_column.append('9')
                else:
                    call = g.split(':')[0]#   In the genotype string, the first element (separated by :) is the actual genotype call
                    if call == '.':
                        g_column.append('9')
                        g_column.append('9')
                    else:
                	#   These are diploid calls, and we are assuming they are unphased
                	#   the are listed in the form allele1/allele2
                	#   with 0 = ref, 1 = alt1, 2 = alt2, and so on...
                        alleles = call.split('/')
                #individual_call = '' #define a dictionary for all of the alleles
                        for x in alleles:
                             if x == '.': #ignor the missing data
                        #print (9,file=open(sys.argv[3], "a"))
                    	         g_column.append('9')
                        		#individual_call += 'N'
                             else:
                                 g_column.append(x)
                        #print (x,file=open(sys.argv[3], "a"))
            print (' '.join(g_column),file=open(sys.argv[2], "a"))

            #print ('\t'.join([chromosome, bp_pos, ref_allele, alt_alleles, chr_nb, maf]))
            #print ('\t'.join([SNP_id, chromosome, "0", bp_pos,ref_allele,alt_alleles],file=open(sys.argv[2], "a"))