# thy_tau22_scRNAseq

## Name
Disentangling sex-specific transcriptomic alterations in the THY-Tau22 mouse model of Alzheimer disease (AD).

## Description
In this project we analyze single cell RNAseq data from mouse cortex to gain insights into AD-associated molecular changes and their sex-dependence for tau-associated pathologies in the cortex as one of the most severely affected brain regions. By studying cell-type specific and cell-type agnostic AD-related changes and their sex-dimorphism for single genes, pathways and cellular sub-networks, we aimed to identify both statistically significant alterations and interpret the upstream mechanisms that control them.

## Project status
We performed the planned computational analysis and drafter the first version of manuscript. At the moment we are actively reviewing and polishing the manuscript draft for submission

## Directory structure

- **scripts** contains all the script including *pre-processing & QC*, *differential expression analysis*, *pathway enrichment analysis*, and *gene regulatory network analysis*
- **data** must contain input data; intermediate and output data files will be created in the directory by the scripts
- **figures** will contain all the figures generated for the manuscript by the scripts; the precomputed figures are available from https://webdav-r3lab.uni.lu/public/data/sbpj-x776/figures/
- **GRN** must contain the JAR files required for running GRN algorithm; there are also output files for different GRN analyses. You can download the .jar files and the generated output files from https://webdav-r3lab.uni.lu/public/data/sbpj-x776/GRN.zip
