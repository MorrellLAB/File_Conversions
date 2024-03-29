#!/usr/bin/python3

#   Script by Jacob Pacheco to convert the output of est-sfs
#   [doi:10.1534/genetics.118.301120] into inferred ancestral state.
#   Should be suitable for anotation of a VCF or creating a list of
#   ancestral state for a list of SNPs.

import vcf
import gzip
import sys


def SecondParse(vcf_file, ESTSFSFILE, output_file, Probabilitypercentagethreshold):
    #  Ancestral allele is original state of nucleotide positon
    #  Alternate allele: also the derived allele, is a variant: 
    #  Reference allele: the state that is found in the reference genome
    #  GZIPPED VCF FILES CURRENTLY DO NOT WORK

    vcf_file = gzip.open(vcf_file, 'rb') if vcf_file.endswith('.gz') else open(vcf_file, 'r') #parses VCF file
    vcf_reader = vcf.Reader(vcf_file)
    VCF_metadata = []

    for record in vcf_reader: # these lines parse VCF data
        chrom = record.CHROM if record.CHROM is not None else '.'
        pos = record.POS - 1 if record.POS is not None else '.'
        record_id = record.ID if record.ID is not None else '.'
        ref = record.REF if record.REF is not None else '.'
        alt = record.ALT[0] if record.ALT and record.ALT[0] is not None else '.'
        ns = record.INFO.get('NS', '.') if record.INFO.get('NS') is not None else '.'
        ac = record.INFO.get('AC', '.') if record.INFO.get('AC') is not None else '.'
        VCF_metadata.append([chrom, pos, record_id, ref, alt, ns, ac])

    

    with open(ESTSFSFILE, "r") as fp1: #opens EST-sfs file and parses data
        with open(output_file, "w") as outfile:
            testfile = fp1.read().splitlines()  
            Minor_or_Major_allele = "0" #zero in case Missing data.
            i = 0  
            for j in range(len(VCF_metadata)):  
                AC = 0 #initializes AC value
                NS = 0 #initializes NS value
                AC = VCF_metadata[j][6][0]
                NS = VCF_metadata[j][5]
                record_id = VCF_metadata[i][2]
    
                
                line = testfile[j]
                Parsed_ESTSFS = line.split() #puts it into a list
        
                if Parsed_ESTSFS[0][0].isdigit() and Parsed_ESTSFS[0][0] != 0 and j> 7:  #Checks if AC and NS values are integer values, otherwises makes them 0
                    if AC == int:
                        AC = AC
                    else:
                        AC = 0

                    if NS == int:
                        NS = NS
                    else:
                        NS = 0

                    if (AC) >= NS/2:  #if AC is larger than NS/2 ancestral is AC
                        Minor_or_Major_allele = VCF_metadata[i][4]
                    else: (AC) <= NS/2 # if NS/2 is larger than AC ancestral is NS
                    Minor_or_Major_allele = VCF_metadata[i][3] #Major allele is the most common, Minor is least common.

         
                    Probability_Percentage = abs((float(Parsed_ESTSFS[2]) *100)) # calculates the probabilty of the ancestral  
                   


                    if Probability_Percentage > Probabilitypercentagethreshold:
                        Final_parsed_string = (VCF_metadata[i][0] + " " + str(VCF_metadata[i][1] + 1) + " "  +str(VCF_metadata[i][3]) + " " +str(VCF_metadata[i][4]) + " " + 'AA=' + "'" +str(Minor_or_Major_allele) + "'" + " " + str(Probability_Percentage/100) + '\n') # output string
                        j += 1
                        outfile.write(Final_parsed_string) #writes final_pasrsed_string to a new file
                    elif Probability_Percentage < Probabilitypercentagethreshold:
                        Probability_Percentage =  abs(100 -Probability_Percentage) 
                        Final_parsed_string = (VCF_metadata[i][0] + " " + str(VCF_metadata[i][1] + 1) + " "  +str(VCF_metadata[i][3]) + " " +str(VCF_metadata[i][4]) + " " + 'AA=' + "'" +str(Minor_or_Major_allele) + "'" + " " + str(Probability_Percentage/100) + '\n') # output string
                        j += 1
                        outfile.write(Final_parsed_string) #writes final_pasrsed_string to a new file if major allele is not most prevalent 
                    i += 1
                else:
                    j += 1


#GZIPPED VCF FILES CURRENTLY DO NOT WORK
vcf_file = sys.argv[1]
ESTSFSFILE = sys.argv[2]
output_file = sys.argv[3]
Probabilitypercentagethreshold = 50

SecondParse(vcf_file, ESTSFSFILE, output_file, Probabilitypercentagethreshold)
