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

### 1. Quality Check, Alignment and meta-gene analysis
You can download the raw sequencing data from [NCBI Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra/) under PRJNA1074553 to `projects/project1/fastq/`.
Here we use the data example file `projects/project1/fastq/Sample_DataExample` to demonstrate the workflow.

run `projects/project1/fastq/Check.sh` and get `projects/project1/fastq/check.log` to verify the integrity of downloaded files.

Download the reference genome and build a reference genome index for the alignment program (bwa). Edit the `projects/project1/conf/config`, and replace the `hg38_FA`, `hg38_INDEX_PATH` variables with the path of your reference genome and reference genome index, respectively. 

We provide an example file including reference-point sites, `projects/conf/config/Reference-point_site_Example.txt`, which is used to perform meta-gene analysis. You can make your own files following the data structure and replace it. The data structure of the file is shown below:
| Gene ID | Chromasome | Reference-point site | Strand | Notes (optional) |
| :---: | :---: | :---: | :---: | :---: |
| ENSG00000173193 | chr3 | 122680839 | + | 50001 |
| ENSG00000112335 | chr6 | 108261246 | - | 50024 |
| ... | ... | ... | ... | ... |

Then replace the `UPSTREAM` and `DOWNSTREAM` variables in `projects/conf/config` with distance (bp) upstream and downstream of the reference-point (such as TSS or TES) selected for meta-gene analysis, and replace the `REGION` variable with the length (bp) between `UPSTREAM` and `DOWNSTREAM`.

Run `projects/start.sh`, and it will automatically perform adaptor cutting, Alignment, meta-gene analysis and nucleotide analysis. 

After completion, you will get `projects/DataExample.rawValue` file of meta gene analysis which can be used to plot meta-gene profile with `projects/plotProfileOfMeta-gene.R` in R, and files in `projects/BaseCount/` of the nucleotide analysis.

### 2. Plot strand-specific heatmap and screenshot
Run `projects/generateStrandSpecificBW.sh` to get the strand-specific bw files:
```
bash generateStrandSpecificBW.sh DataExample.bed
```
After that, you will get the files of `.forward.bw` and `.reverse.bw`, which can be used to plot screenshot in IGV.

Run `projects/plotHeatmap.sh` to get the heatmap of template strand signal, `projects/DataExample.bed.TS.pdf`, and non-template strand signal, `projects/DataExample.bed.NTS.pdf`:
```
bash plotHeatmap.sh DataExample.bed
```

### 3. Analysis of sample correlations
Run `projects/plotCorrelationForGenes.sh` to get the analysis result of sample correlations among active genes, `projects/correlationScatter.pdf`:
```
bash plotCorrelationForGenes.sh DataExample1.bed DataExample2.bed
```

Or run `projects/bed2bam.sh` and `projects/plotCorrelationForBins.sh` to get the analysis result of sample correlations for the entire genome, `correlationHeatmapPerBin.pdf`:
```
bash bed2bam.sh DataExample3.bed
bash bed2bam.sh DataExample4.bed
bash bed2bam.sh DataExample5.bed
bash bed2bam.sh DataExample6.bed

bash plotCorrelationForBins.sh DataExample3.bed DataExample4.bed DataExample5.bed DataExample6.bed
```

### 4. Peak calling
Run `generate_bed.py` to get files in bed format for peak calling and perform peak calling using MACS2 with no model, 0 bp shift, and 0 bp extension option :
```
macs2 callpeak -t ChIP1.bed -c Input1.bed -g hs -p 0.01 -n ChIP1 --nomodel -f BED -B --SPMR
```
---


## Authors
Codes were generated by Maoxiang Qian, Jiayong Yin and Yongchang Zhu, Institutes of Biomedical Sciences, Fudan University, Shanghai 200032, China. Contact information: yongchangzhu21@m.fudan.edu.cn
