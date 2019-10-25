# File_Conversions
Scripts to convert between file formats for various analyses

| Script | Author | Project Nickname | Description |
|--------|--------|------------------|-------------|
| ANNOVAR_to_Effects.py| Tom Kono | Genomic Prediction | ANNOVAR output to SNP effects table |
| ATCG_to_01.py | Tom Kono | Soja Env. Association | Convert ATCG genotypes to 0/1 for STRUCTURE |
| Alchemy_to_PLINK.py | Tom Kono | Genomic Prediction | Convert ALCHEMY output report to PLINK PED |
| GenoMatrix_to_FASTA.py | Tom Kono | Soja Env. Association | Convert a genotyping matrix to FASTA for libsequence summaries. Removes monomorphic markers. |
| GenoMatrix_to_HierFstat.py | Tom Kono | Soja Env. Association | Convert a genotyping matrix to haploid hierfstat input |
| HapMap_to_PLINK.py | Tom Kono | Genomic Prediction | Convert a TASSEL HapMap file to PLINK PED |
| PLINK_to_NicholsonFST.py | Tom Kono | Soja Env. Association | Convert a PLINK PED+CLST combination to inputs for Nicholson FST estimation in the R `popgen` package |
| PLINK_to_mpMap.py | Tom Kono | Genomic Prediction | Convert PLINK PED/MAP to inputs for the R `mpMap` package |
| PolyTable_to_Fasta.py | Tom Kono | North America FST | Convert a polytable-like format to FASTA, removing monomorphic positions |
| VCF_to_Htable.py | Tom Kono | North America FST | Convert from a VCF to a polytable-like format. Codes heterozygous genotypes as missing |
| VCF_to_HapMap.py | Tom Kono | Genomic Prediction | Convert from VCF to a TASSEL HapMap format |
| VCF_to_Infocalc.py | Tom Kono | Fly GBS | Convert from VCF to Infocalc (Rosenberg 2003; Rosenberg 2005) input format |
| VCF_to_Phylip.py | Tom Kono | NA | Convert from a VCF to an input for PHYLIP programs |
| Format_HTable.py | Tom Kono | North America FST | Subset or reorder samples in a Hudson table for input into libsequence tools that require partitioning |
| Barley_Parts_to_Pseudomolecules.py | Tom Kono | Genomic Prediction | Convert VCF or BED coordinates from the IPK psuedomolecule parts assembly to the non-split coordinates. |
| VCF_to_XPCLR.py | Li Lei | Selective Sweep | Convert VCF to XPCLR geno and map file format |
