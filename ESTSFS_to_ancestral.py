import vcf
import gzip
import sys


def SecondParse(vcf_file, testfile, output_file):
    
    vcf_file = gzip.open(vcf_file, 'rb') if vcf_file.endswith('.gz') else open(vcf_file, 'r')
    vcf_reader = vcf.Reader(vcf_file)
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

    
    # vcf_reader = vcf.Reader(open('vcf_file', 'r'))
    # info_lines = []
    # for line in vcf_reader.metadata['INFO']:
    #     info_lines.append(str(line))
    # info_lines = info_lines[-5:]  
    # for line in info_lines:
    #     print(line)
            

    with open(testfile, "r") as fp1:
        with open(output_file, "w") as outfile:
            testfile = fp1.read().splitlines()
            Ancestrial_allele = "0"
            i = 0  
            for j in range(len(VCF_metadata)):
                AC = 0
                NS = 0
                AC = VCF_metadata[j][6][0]
                NS = VCF_metadata[j][5]
    
                
                line = testfile[j]
                Parsed_ESTSFS = line.split()
        
                if Parsed_ESTSFS[0][0].isdigit() and Parsed_ESTSFS[0][0] != 0 and j> 7:
                    if AC == int:
                        AC = AC
                    else:
                        AC = 0

                    if NS == int:
                        NS = NS
                    else:
                        NS = 0

                    if (AC) >= NS/2:
                        Ancestrial_allele = VCF_metadata[i][4]
                    else: (AC) <= NS/2
                    Ancestrial_allele = VCF_metadata[i][3]

         
                    Probability_Percentage = abs((float(Parsed_ESTSFS[2]) *100))
                    Final_parsed_string = (VCF_metadata[i][0] + " " + str(VCF_metadata[i][1] + 1) + " "  +str(VCF_metadata[i][3]) + " " +str(VCF_metadata[i][4]) + " " + 'AA=' + "'" +str(Ancestrial_allele) + "'" + " " + str(Probability_Percentage) + '\n')
                    j += 1
                    outfile.write(Final_parsed_string)
                    i += 1
                else:
                    j += 1

vcf_file = sys.argv[1]
testfile_path = sys.argv[2]
output_file = sys.argv[3]
SecondParse(vcf_file, testfile_path, output_file)
