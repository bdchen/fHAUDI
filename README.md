# fHAUDI

## Introduction
fHAUDI extends [HAUDI] (https://github.com/frankp-0/HAUDI) to incorporate functional annotations into a local ancestry-informed PRS. fHAUDI uses annotation-specific penalties to prioritize more annotated variants, which are more likely to be causal.

## Installation

```{r}
remotes::install_github("frankp-0/HAUDI")
```

## File Preparation

fHAUDI uses the same file preparation and functions as HAUDI. For clarity, the file preparation steps are repeated here. 

### 1. Prepare input VCF file

fHAUDI requires as input a tabix-indexed vcf file, with AN1 and AN2 subfields
for haplotype local ancestry (as is produced by [flare](https://github.com/browning-lab/flare)).
If you use flare to estimate local ancestry, you can annotate your original vcf
file using a command like: `bcftools annotate -c FORMAT -a flare.anc.vcf.gz
target.vcf.gz -Oz -o target.anc.vcf.gz`.

This will produce a VCF file where the AN1, AN2 fields are missing for variants
not included in the reference panel. To interpolate local ancestry in
the following step, ensure that the VCF file is sorted.

### 2. Creating file-backed matrices

fHAUDI uses file-backed matrices, implemented in the `bigstatsr`
package to store genotype/ancestry data. To create this file-backed
matrix, use the function `make_fbm`, which takes the VCF file prepared
in the previous step as input. An example may be:

```{r}
fbm_result <- make_fbm(
  vcf_file = "target.anc.vcf.gz",
  fbm_pref = "target",
  chunk_size = 400,
  min_ac = 10,
  geno_format = "GT",
  anc_names = c("Pop_01", "Pop_02", "Pop_03") 
)
```

This command would return a list containing:

1) an object of class FBM.code256 containing genotype/ancestry data
2) a data frame containing SNP information (chromosome, position, etc.)

Users may save these objects with: `saveRDS(fbm_result$FBM, file="target.rds")`
and `write.table(fbm_result$info, file="target_info.txt")`


## Run fHAUDI

To run fHAUDI, simply use the function `fHAUDI`. An example is provided below
