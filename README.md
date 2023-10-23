# thy_tau22_scRNAseq

## Name
Single cell transcriptome analysis of the THY-Tau22 mouse model of Alzheimer's disease reveals sex-dependent dysregulations.

## Description
In this project, we analyze single cell RNAseq data from mouse cortex to gain insight into AD-associated molecular changes and their sex dependency for tau-associated pathologies in the cortex as one of the most severely affected brain regions. By examining cell type specific and cell type agnostic AD-related changes and their sex dimorphism for individual genes, pathways and cellular subnetworks, we aimed to both identify statistically significant changes and interpret the upstream mechanisms controlling them.

## Directory structure

- **scripts** contains all scripts including pre-processing & QC, differential expression analysis, pathway enrichment analysis, and gene regulatory network analysis.
- **data** contains input data; intermediate and output data files are generated by the scripts in this directory.
- **figures** contains all figures generated for the manuscript by the scripts; the precomputed figures are available from https://webdav-r3lab.uni.lu/public/data/sbpj-x776/figures/
- **GRN** contains the JAR files needed to run the GRN algorithm and output files for various GRN analyses. The .jar files and the generated output files can be downloaded from https://webdav-r3lab.uni.lu/public/data/sbpj-x776/GRN.zip.
