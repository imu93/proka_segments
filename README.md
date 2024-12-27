# proka_segments
 
## Introduction

Due to the compactness of prokaryotic genomes, some genomic features can occupy the same genomic space. This can be problematic for downstream analyses, such as RNA annotation. The scripts in this repository are designed to manipulate genomic ranges (GFF files) to produce a non-overlapping annotation, where all bases are categorized as:

	1.	Protein-coding sequences

	2.	Non-coding RNAs

	3.	Repetitive elements

	4.	Intergenic regions


In this first release, the repository contains individual (but modular) scripts.

## Installation

Since all scripts are written in R and dependencies are managed via the pacman package, they can be executed using the following conda environment:

```
conda create -n r_env -c conda-forge -c bioconda r-base r-pacman r-dplyr bioconductor-rtracklayer bioconductor-rsamtools
```

## Usage
After activating the conda environment, you can edit each script to suit your specific genome and then run it using Rscript.

Note: This usage section will be updated in future versions with more specific instructions.

## Script Descriptions
1.	01.disjoin_protein_coding_genes.R

This script assigns overlapping bases between genes (operons). Overlapping bases are assigned based on the length of the overlapping elements. After the initial assignment, the ranges are updated in case multiple overlaps occur within a single range.

	2.	02.reduce_ncrna_genes.R

Although bacterial genomes typically have a low proportion of overlapping bases among non-coding RNAs (ncRNAs), this script uses the reduce function to collapse overlapping elements of the same category (e.g., tRNA-tRNA or rRNA-rRNA). The metadata of the longer element is preserved.

	3.	03.reduce_repeats.R

Similar to ncRNAs, bacterial genomes generally exhibit a low proportion of transposable elements (TEs). In cases of overlap between TEs, the reduce function is used to collapse them.

	4.	04.segmented_prokka.R (strand-aware)

This script produces the final annotation, resolving overlaps among genes, ncRNAs, and TEs based on a hierarchical classification. Currently, ncRNAs (e.g., rRNA, tRNA) are placed at the top of the hierarchy, followed by genes, and then TEs (which “kill” genes). The order is arbitrary but can be adjusted according to the user’s preferences.


Note: Contiguous ranges can cause issues with setdiff, so consider using desired_shrink with moderation.

Release Description

Beta release of the proka_segments package scripts.
