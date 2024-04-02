#!/usr/bin/python3 

"""Convert a VCF and BED files with genotype information for outgroup samples 
into input for est-sfs. The primary dependency is the  PyVCF library for VCF 
parsing.
"""  

import vcf
import gzip
import sys


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

    bed_file1 = gzip.open(bed_file1, 'rt') if bed_file1.endswith('.gz') else open(bed_file1, 'r')  # checks if gzipped then opens the file
    bed_file2 = gzip.open(bed_file2, 'rt') if bed_file2.endswith('.gz') else open(bed_file2, 'r')  # checks if gzipped then opens the file
    if bed_file3:
        bed_file3 = gzip.open(bed_file3, 'rt') if bed_file3.endswith('.gz') else open(bed_file3, 'r')  # checks if gzipped then opens the file

    #  GZIPPED VCF FILES CURRENTLY DO NOT WORK        
    vcf_file = gzip.open(vcf_file, 'r') if vcf_file.endswith('.gz') else open(vcf_file, 'r')  # checks if gzipped then opens the file

    with bed_file1 as f1, bed_file2 as f2:  # splits the bedfiles by white space
        bed_file1 = f1.read().splitlines()
        bed_file2 = f2.read().splitlines()
        if bed_file3:
            with bed_file3 as f3:
                bed_file3 = f3.read().splitlines()
    
    vcf_reader = vcf.Reader(vcf_file)
    VCF_metadata = []
    
    for record in vcf_reader:
        chrom = record.CHROM if record.CHROM is not None else '.'  # Extract chromosome information, replace with '.' if None
        pos = record.POS - 1 if record.POS is not None else '.'  # Extract variant position, subtract 1 to convert to zero-based, replace with '.' if None
        record_id = record.ID if record.ID is not None else '.'  # Extract SNP name, replace with '.' if None
        ref = record.REF if record.REF is not None else '.'  # Extract reference allele, replace with '.' if None
        alt = record.ALT[0] if record.ALT and record.ALT[0] is not None else '.'  # Extract alternate allele, replace with '.' if None
        ns = record.INFO.get('NS', '.') if record.INFO.get('NS') is not None else '.'  # Extract NS (Number of Samples With Data) from INFO, replace with '.' if None
        ac = record.INFO.get('AC', '.') if record.INFO.get('AC') is not None else '.'  # Extract AC (Alternate Allele Count) from INFO, replace with '.' if None
        VCF_metadata.append([chrom, pos, record_id, ref, alt, ns, ac])  # appends all dat into the list
        
    with open(output_file, "w") as outfile:
        for i in range(len(bed_file1)): 
            bed_file1_split = bed_file1[i].split()  # Split the elements at each whitespace
            bed_file2_split = bed_file2[i].split()  # Split the elements at each whitespace
            if bed_file3:
                bed_file3_split = bed_file3[i].split()  # Split the elements at each whitespace
            nucleotide_Base_binary = {"A": "1,0,0,0", "C": "0,1,0,0", "G": "0,0,1,0", "T": "0,0,0,1"}  # Dictionary to map nucleotide bases to binary values
            nucleotide_ref_binary = {"AT": [1,0,0,1], "AG": [1,0,1,0], "AC": [1,1,0,0], "CT": [0,1,0,1], "CG": [0,1,1,0], "CA": [1,1,0,0], "TG": [0,0,1,1], "TC": [0,0,1,1], "TA": [1,0,0,1], "GT": [0,0,1,1], "GC": [0,0,1,1], "GA": [1,0,1,0], 'A':[1,0,0,0], "C": [0,1,0,0], "G": [0,0,1,0], "T": [0,0,0,1]} # Dictionary to map reference alleles to binary values
            
            AC = (VCF_metadata[i][-1][0])    # Extract the AC value from the VCF metadata
            if isinstance(AC, int):  # Check if AC is an integer
                OutputAC = f'AC={int(AC)}'  # If AC is an integer, set the output string with the AC value
                calcAC = AC
            else:
                calcAC = 0  # If AC is not an integer, set calcAC to 0. Calc is used for calculations to find Alternate allele freq, while output is a string representation for output
                OutputAC = '.'
            
            NS = (VCF_metadata[i][5]) 
            if isinstance(NS, int):  # Ensures NS is an int
                OutputNS = f'NS={int(NS)}'  # if it is in int output NS represents the string of NS
                calcNS = NS  # calc is used for NS calculations to calculate the reference allele freq.
            else:
                calcNS = 0
                OutputNS = '.'
            

            ReferencealleleFreq =  ((int(calcNS)*2) +  (int(-calcAC)))//2  # calculates Reference allele frequency 
            Alternateallelefreq =  int(calcAC)//2  # calculates alternate allele frequency 

            reference_allele = str(VCF_metadata[i][4])  # extracts ref allele from input VCF file
            alternate_allele = str(VCF_metadata[i][3])  # extracts ref allele from input VCF file
            genotypesReference = nucleotide_ref_binary.get(reference_allele, [0,0,0,0])
            genotypesReference = list(map(lambda x: x * Alternateallelefreq, genotypesReference))  # assigns the values to a binary representation 

            genotypesAlternate = nucleotide_ref_binary.get(alternate_allele, [0,0,0,0]) 
            genotypesAlternate= list(map(lambda x: x * ReferencealleleFreq, genotypesAlternate))  # assigns the values to a binary representation using dictionary above

            genotypes = [] 
            for j in range(len(genotypesAlternate)):   # puts together genotypes reference and alleles use as a string in the output.
                genotypes.append(genotypesReference[j] + genotypesAlternate[j]) 
            result_str = ','.join(str(k) for k in genotypes)

            bed_file1_split = bed_file1[i].split() 
            BED_file1_binary_conversion = nucleotide_Base_binary.get(bed_file1_split[-1], "0,0,0,0")  # turns binary lists into binary strings. (when inputted before they were casted to list)
            bed_file2_split = bed_file2[i].split()
            BED_file2_binary_conversion = nucleotide_Base_binary.get(bed_file2_split[-1], "0,0,0,0")  # turns binary lists into binary strings.

            if len(sys.argv) == 6:  # Uses this output if there were three bed files in input
                bed_file3_split = bed_file3[i].split()
                BED_file3_binary_conversion = nucleotide_Base_binary.get(bed_file3_split[-1], "0,0,0,0")  # turns binary lists into binary strings.
                converted_file = ((VCF_metadata[i][2]) + "\t" + (VCF_metadata[i][3]) + "\t" + str(VCF_metadata[i][4]) + "\t" + OutputAC + "\t"+ OutputNS + "\t" + bed_file1_split[-1] + "\t" + bed_file2_split[-1] + "\t" + bed_file3_split[-1] +  "\t") #String representation of VCF and Bed info
                EST_SFS_VCF_compatible_file = (result_str + "\t" + BED_file1_binary_conversion + "\t" + BED_file2_binary_conversion + "\t" + BED_file3_binary_conversion + "\n") # adds binary info
                outfile.write(converted_file +  EST_SFS_VCF_compatible_file)

            if len(sys.argv) == 5:  # Uses this output if only two bedfiles were in input
                converted_file = ((VCF_metadata[i][2]) + "\t" + (VCF_metadata[i][3]) + "\t" + str(VCF_metadata[i][4]) + "\t" + OutputAC + "\t"+ OutputNS + "\t" + bed_file1_split[-1] + "\t" + bed_file2_split[-1] +  "\t") #String representation of VCF and Bed info
                EST_SFS_VCF_compatible_file = (result_str + "\t" + BED_file1_binary_conversion + "\t" + BED_file2_binary_conversion + "\n")  # add binary info
                outfile.write(converted_file +  EST_SFS_VCF_compatible_file)

if len(sys.argv) == 6:    # assigns values to be used in command line arguments example command below
    bed_file1 = sys.argv[1] 
    bed_file2 = sys.argv[2]
    bed_file3 = sys.argv[3]
    vcf_file = sys.argv[4]
    output_file = sys.argv[5] 
    est_vcf_convert(sys.argv[1], sys.argv[2], sys.argv[4], sys.argv[5], sys.argv[3])  # example command you type for 3 bedfiles python3 VCF_to_ESTSFS.py bedFile1 bedFile2 bedFile3 vcfFile nameOfFile


elif len(sys.argv) == 5:
    bed_file1 = sys.argv[1]
    bed_file2 = sys.argv[2]
    vcf_file = sys.argv[3]
    output_file = sys.argv[4]
    est_vcf_convert(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])  # example command you type for 2 bedfiles python3 VCF_to_ESTSFS.py bed1File bed2File vcfFile nameOfFile
         
