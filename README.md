# Transcriptomic profiling of PD brains

## Background

This repository contains analysis code for the publication found [here](add link). In this pubblication we performed RNA-sequencing after ribosomal RNA removal of 84 frontal cortex samples. Samples were collected postmortem from 23 non-neurological individuals and 61 individuals with varying degree of α-synuclein-related pathology. The latter were samples collected from donors with iLDB, PD or PDD and were split into three neuropathological groups based on their Braak stage: one group consisting of individuals at Braak stages 1, 2, 3 and 4 (n = 19), one group consisting of individuals at Braak stage 5 (n = 19) and one group consisting of individuals at Braak stage 6 (n = 23). Cell composition was assessed for all the samples. Differential expression analysis was performed both across neuropathological stages and in a pairwise manner.

## Overview

R_scripts folder contains:
- 1.Transcript_to_Gene_&_Pre-filtering.R
- 2.Demographics.R
- 3.Quality_&_Seq_Stats.R
- 4.Cell_Composition.R
- 5.Covariate_Selection.R
- 6.DE_LRT_&_Patterns.R
- 7.Heatmap.R
- 8.DE_Wald.R
- 9.Functional_Enrichment.R
- 10.eQTLs.R

Files folder contains:
- Coldata_Final_84.csv
- ClinicopathVUmc_forChiara.csv
- Sample_ID.csv
- Seq_stats.csv
- degradation.matrix.84.txt
- MReads_84.cs

## Citation

If you use any of the code or data from this repository, please cite our paper.



