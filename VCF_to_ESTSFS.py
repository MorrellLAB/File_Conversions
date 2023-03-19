import vcf
import gzip
import sys

def VCF_to_ESTSFS(bed_file1, bed_file2, vcf_file, output_file, bed_file3=None):
    '''This script takes three tab-delimited BED files with the following fields:
    1 -chromosome
    2 -zero-based variant position
    3 -one-based variant position
    4 -SNP name (or ".")
    5 -reference allele
    6 -alternate allele
    7 -allele present at the position in and outgroup sample
    8 -dependencies VCF, gzip, sys 
    The BED files need to provided in order so that outgroup samples most closely
    related to the focal species appear earlier (farther left in file order).
    Outputs are designed for input for ancestral state inference in EST-SFS
    and a VCF compatible file'''

    bed_file1 = gzip.open(bed_file1, 'r') if bed_file1.endswith('.gz') else open(bed_file1, 'r') #checks if gziped and if two or three files are inputted
    bed_file2 = gzip.open(bed_file2, 'r') if bed_file2.endswith('.gz') else open(bed_file2, 'r')
    if bed_file3:
        bed_file3 = gzip.open(bed_file3, 'r') if bed_file3.endswith('.gz') else open(bed_file3, 'r')

    vcf_file = gzip.open(vcf_file, 'r') if vcf_file.endswith('.gz') else open(vcf_file, 'r')

    with bed_file1 as f1, bed_file2 as f2:  #bed allele parsed
        bed_file1 = f1.read().splitlines()
        bed_file2 = f2.read().splitlines()
        if bed_file3:
            with bed_file3 as f3:
                bed_file3 = f3.read().splitlines()
    
    vcf_reader = vcf.Reader(vcf_file) #VCF info parsed
    VCF_metadata = []
    for record in vcf_reader:
        chrom = record.CHROM if record.CHROM is not None else '.'
        pos = record.POS - 1 if record.POS is not None else '.'
        record_id = record.ID if record.ID is not None else '.'
        ref = record.REF if record.REF is not None else '.'
        alt = record.ALT[0] if record.ALT and record.ALT[0] is not None else '.'
        ns = record.INFO.get('NS', '.') if record.INFO.get('NS') is not None else '.'
        ac = record.INFO.get('AC', '.') if record.INFO.get('AC') is not None else '.'
        VCF_metadata.append([chrom, pos, record_id, ref, alt, ns, ac])

    with open(output_file, "w") as outfile:
        for i in range(len(bed_file1)): 
            bed_file1_split = bed_file1[i].split()
            bed_file2_split = bed_file2[i].split()
            if bed_file3:
                bed_file3_split = bed_file3[i].split()
            nucleotide_Base_binary = {"A": "1,0,0,0", "C": "0,1,0,0", "G": "0,0,1,0", "T": "0,0,0,1"} 
            nucleotide_ref_binary = {"AT": [1,0,0,1], "AG": [1,0,1,0], "AC": [1,1,0,0], "CT": [0,1,0,1], "CG": [0,1,1,0], "CA": [1,1,0,0], "TG": [0,0,1,1], "TC": [0,0,1,1], "TA": [1,0,0,1], "GT": [0,0,1,1], "GC": [0,0,1,1], "GA": [1,0,1,0], 'A':[1,0,0,0], "C": [0,1,0,0], "G": [0,0,1,0], "T": [0,0,0,1]}
            
            AC = (VCF_metadata[i][-1][0])  
            if isinstance(AC, int):
                OutputAC = f'AC={int(AC)}'
                calcAC = AC
            else:
                calcAC = 0
                OutputAC = '.'
            
            NS = (VCF_metadata[i][5]) 
            if isinstance(NS, int):
                OutputNS = f'NS={int(NS)}'
                calcNS = NS
            else:
                calcNS = 0
                OutputNS = '.'
            

            ReferencealleleFreq =  ((int(calcNS)*2) +  (int(-calcAC)))//2 
            Alternateallelefreq =  int(calcAC)//2

            reference_allele = str(VCF_metadata[i][4])
            alternate_allele = str(VCF_metadata[i][3])
            genotypesReference = nucleotide_ref_binary.get(reference_allele, [0,0,0,0])
            genotypesReference = list(map(lambda x: x * Alternateallelefreq, genotypesReference))

            genotypesAlternate = nucleotide_ref_binary.get(alternate_allele, [0,0,0,0])
            genotypesAlternate= list(map(lambda x: x * ReferencealleleFreq, genotypesAlternate))

            genotypes = []
            for j in range(len(genotypesAlternate)):
                genotypes.append(genotypesReference[j] + genotypesAlternate[j])
            result_str = ','.join(str(k) for k in genotypes)

            
            bed_file1_split = bed_file1[i].split()
            BED_file1_binary_conversion = nucleotide_Base_binary.get(bed_file1_split[-1], "0,0,0,0")
            bed_file2_split = bed_file2[i].split()
            BED_file2_binary_conversion = nucleotide_Base_binary.get(bed_file2_split[-1], "0,0,0,0")

            if len(sys.argv) == 6: #three bedfiles
                bed_file3_split = bed_file3[i].split()
                BED_file3_binary_conversion = nucleotide_Base_binary.get(bed_file3_split[-1], "0,0,0,0")
                converted_file = ((VCF_metadata[i][2]) + "\t" + (VCF_metadata[i][3]) + "\t" + str(VCF_metadata[i][4]) + "\t" + OutputAC + "\t"+ OutputNS + "\t" + bed_file1_split[-1] + "\t" + bed_file2_split[-1] + "\t" + bed_file3_split[-1] +  "\t")
                EST_SFS_VCF_compatible_file = (result_str + "\t" + BED_file1_binary_conversion + "\t" + BED_file2_binary_conversion + "\t" + BED_file3_binary_conversion + "\n")
                outfile.write(converted_file +  EST_SFS_VCF_compatible_file)

            if len(sys.argv) == 5: #two bedfiles
                converted_file = ((VCF_metadata[i][2]) + "\t" + (VCF_metadata[i][3]) + "\t" + str(VCF_metadata[i][4]) + "\t" + OutputAC + "\t"+ OutputNS + "\t" + bed_file1_split[-1] + "\t" + bed_file2_split[-1] +  "\t")
                EST_SFS_VCF_compatible_file = (result_str + "\t" + BED_file1_binary_conversion + "\t" + BED_file2_binary_conversion + "\n")
                outfile.write(converted_file +  EST_SFS_VCF_compatible_file)

if len(sys.argv) == 6:    #takes in command line arguments example command below
    bed_file1 = sys.argv[1]
    bed_file2 = sys.argv[2]
    bed_file3 = sys.argv[3]
    vcf_file = sys.argv[4]
    output_file = sys.argv[5] 
    est_vcf_convert(sys.argv[1], sys.argv[2], sys.argv[4], sys.argv[5], sys.argv[3]) #example command you type for 3 bedfiles python3 VCF_to_ESTSFS.py bedFile1 bedFile2 bedFile3 vcfFile nameOfFile

elif len(sys.argv) == 5:
    bed_file1 = sys.argv[1]
    bed_file2 = sys.argv[2]
    vcf_file = sys.argv[3]
    output_file = sys.argv[4]
    est_vcf_convert(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]) # example command you type for 2 bedfiles python3 VCF_to_ESTSFS.py bed1File bed2File vcfFile nameOfFile
         