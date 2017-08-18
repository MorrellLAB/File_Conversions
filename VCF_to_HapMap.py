#!/usr/bin/env python
"""Dumb script to convert a VCF to a hapmap format for use with TASSEL. Most of
the information is the same (and it is in the same direction!) so it should be
pretty simple. The HapMap fields are:
    1) SNP ID
    2) Alleles
    3) Chromosome
    4) Position
    5) Strand
    6) Assembly ver.
    7) Center (NA for us)
    8) protLSID (NA again)
    9) assayLSID (NA)
    10) panelLSID
    11) QCode (NA)
    12 ... end) Sample IDs"""

import sys

with open(sys.argv[1], 'r') as f:
    for line in f:
        if line.startswith('##'):
            continue
        elif line.startswith('#CHROM'):
            samplenames = line.strip().split('\t')[9:]
            #   Print a header
            print 'rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanel\tQCode\t' + '\t'.join(samplenames)
        else:
            tmp = line.strip().split('\t')
            rs = tmp[2]
            alleles = tmp[3] + '/' + tmp[4]
            chrom = tmp[0]
            pos = tmp[1]
            strand = '+'
            a_ver = 'RefSeq1.0'
            center = 'NA'
            protLSID = 'NA'
            assayLSID = 'NA'
            panelLSID = 'Genomic Prediction'
            qcode = 'NA'
            genotypes = []
            for g in tmp[9:]:
                gt = g.split(':')[0]
                if gt == './.' or gt == '.':
                    genotypes.append('NN')
                elif gt == '0/0':
                    genotypes.append(tmp[3] + tmp[3])
                elif gt == '0/1':
                    genotypes.append(tmp[3] + tmp[4])
                elif gt == '1/0':
                    genotypes.append(tmp[3] + tmp[4])
                elif gt == '1/1':
                    genotypes.append(tmp[4] + tmp[4])
            print '\t'.join(
                [rs, alleles, chrom, pos, strand, a_ver, center, protLSID, assayLSID, panelLSID, qcode] + genotypes
                )
