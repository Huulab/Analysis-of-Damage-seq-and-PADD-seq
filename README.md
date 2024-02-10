# Codes for Damage-seq and PADD-seq
This repository stores original codes used in ***Coordination of transcription-coupled repair and repair-independent release of stalled RNA polymerase II in response to transcription-blocking lesions***.

## Software Requiremnents

### Linux Tools
+ snakemake >= 5.5.4
+ FastQC >= 0.11.9
+ bgzip >= 1.9
+ Cutadapt >= 1.12
+ trim_galore >= 0.6.7
+ samtools >= 1.9
+ bwa >= 0.7.17
+ sambamba >= 0.8.1
+ bedtools >= 2.30.0
+ deeptools >= 3.5.1
+ MACS2 >= 2.2.7.1

### Python
+ Pyhton version >= 3.7.x
+ numpy >= 1.19.2
+ pandas >= 1.1.2
+ py2bit >= 0.3.0
+ re >= 2.2.1
+ pyfaidx >= 0.5.9.1
+ pysam >= 0.16.0.1

### R
+ R version >= 3.6.1
+ ggplot2 >= 3.3.0
+ stringr >= 1.4.0

### Windows
+ IGV >=  2.9.2


---

## Workflow Example

### 1. Analysis of PADD-seq

You can download the raw sequencing data from [NCBI Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra/) under PRJNA1074553 to `PADD-seq/projects/project1/fastq/`.
Here we use the data example file `PADD-seq/projects/project1/fastq/Sample_DataExample` to demonstrate the workflow.

run `PADD-seq/projects/project1/fastq/Check.sh` and get `PADD-seq/projects/project1/fastq/check.log` to verify the integrity of downloaded files.

Download the reference genome and build a reference genome index for the alignment program (bwa). Edit the `PADD-seq/projects/project1/conf/config`, and replace the path of the variables (`hg38_FA`, `hg38_INDEX_PATH`, etc.) with your own path. 

We provide an example file including reference-point sites, `PADD-seq/annotation/BJ_hg38_FPKM-gt1_GeneDistance-gt2k_GeneLength-gt50k.TSSdat.txt`, which is used to perform meta-gene analysis. You can make your own files following the data structure and replace it. The data structure of the file is shown below:
| Gene ID | Chromasome | Reference-point site | Strand | Notes (optional) |
| :---: | :---: | :---: | :---: | :---: |
| ENSG00000173193 | chr3 | 122680839 | + | 50001 |
| ENSG00000112335 | chr6 | 108261246 | - | 50024 |
| ... | ... | ... | ... | ... |

Then replace the `UPSTREAM` and `DOWNSTREAM` variables in `PADD-seq/projects/project1/conf/config` with distance (bp) upstream and downstream of the reference-point (such as TSS or TES) selected for meta-gene analysis, and replace the `REGION` variable with the length (bp) between `UPSTREAM` and `DOWNSTREAM`.

Run `PADD-seq/projects/project1/start.sh`, and it will automatically perform adaptor cutting, Alignment, meta-gene analysis, nucleotide analysis, etc. 

After completion, you will get `PADD-seq/projects/project1/DataExample.TSS(TES).rawValue` files of meta gene analysis which can be used to plot meta-gene profile with `PADD-seq/projects/project1/plotProfileOfMeta-gene.R` in R, files in `PADD-seq/projects/project1/BaseCount/` of the nucleotide analysis, files with the .bw suffix which can be used to plot screenshots in IGV, and quantified files `PADD-seq/projects/project1/DataExample.GeneByGeneRPKM.TS(NTS).txt` which can be used to plot violin plots.

### 2. Analysis of Damage-seq

refer to the "Analysis of PADD-seq".

## Authors
Codes were uploaded by Yongchang Zhu, Institutes of Biomedical Sciences, Fudan University, Shanghai 200032, China. Contact information: yongchangzhu21@m.fudan.edu.cn
