This repository contains the codes used to produce the main and supplementary figured of the manuscript  "Association mapping of colour variations in a butterfly provides evidences that a supergene locks together a cluster of adaptive loci"

The collection of datasets used in the manuscript can be found at Figshare (doi: 10.6084/m9.figshare.19706320)

### Datasets that contain the result of association analyses
AllHn0_Antpost.PC_1_PC_2_PC_3_PC_4.Assoc.mqfam.total.reformate.CHRPos.WithEmpP # Hn0 multivariate association; Only with 4 components (Figure 3)
AllSample.SupergeneAndFlank.snps.Hn0.CPM_Antpost_Hn0.AllPCs.Assoc.mqfam.total.reformate.CHRPos.WithEmpP # Hn0 multivariate associations ; Wwith all components (Figure S15)
AllSample.SupergeneAndFlank.snps.Hn123.CPM_Antpost_Hn123.PC_1_PC_2_PC_3.Assoc.mqfam.total.reformate.CHRPos.WithEmpP # Hn123 multivariate associations ; Only with 4 components (Figure 3)
AllSample.SupergeneAndFlank.snps.Hn123.CPM_Antpost_Hn123.AllPCs.Assoc.mqfam.total.reformate.CHRPos.WithEmpP # Hn123 multivariate associations ; With all components (Figure S14)
AllSample.SupergeneAndFlank.snps.Hn123.CPM_Post_Hn123.AllPartition.Genotypic # Hn123 Univariate associations (Figure 3 and S16-17)
AllSample.SupergeneAndFlank.snps.Hn0.AllPartition.Genotypic # Hn0 Univariate associations (Figure 3 and S18-19)

### Genomic and phenotypic PCA analyses
AllSubSamp_AllRegion.PCA.txt # Contain the result of genomic PCA (for Fig S1-S5)
CPM_AllPartition_AnnotedSample.txt # Contain the result of phenotypic PCA computed on specific part of the wing (for Fig S9)
CPM_Antpost_sPee_WithGeno_WithPheno # Contain the result of phenotypic PCA (for Fig 1)
Expression24HtableWithPos.txt #RNAseq expression 24h post pupae formation (for Fig 20)

### RNAseq analyses
ExpressionPPtableWithPos.txt #RNAseq expression at pre-pupae stage (for Fig 20)
Hmel2_AnnotationGoodName_CHRPos.gff #Gene annotation of the Hmel2 genome

### Hmel 2 genome information
Hmel2_chromosomes.agp # Scaffold mapping of the Hmel2 genome
Weight_With_Pos_Supergene_Twistt4.Header.sort.SW25000.txt #Result of Twistt analyses

### Code to produce figures from the datasets:
GWAS_PaperPlot_Submit.R
