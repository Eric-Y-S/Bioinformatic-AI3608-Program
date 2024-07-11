# Bioinformatic-AI3608-Project
## Description
This repository deposits the codes and results of 2024 Bioinformatic course (AI3608) project.

Our project reanalysed the RNA-seq data and reproduced partial results from the paper:

- Lin, X., Swedlund, B., Ton, ML.N. et al. Mesp1 controls the chromatin and enhancer landscapes essential for spatiotemporal patterning of early cardiovascular progenitors. Nat Cell Biol 24, 1114–1128 (2022). https://doi.org/10.1038/s41556-022-00947-3

## Contributors
Zikun Yang(杨子坤), Zhou Peng(彭周)

## Results
#### 1. Implemented a transcripome data processing pipeline
We implemented a transcripome data processing pipeline based on Snakemake. The pipeline can automatically process raw fastq data from SRA and go on quality control, mapping, count reads and normalization. The detailed procedure is shown below:
![image](https://github.com/Eric-Y-S/Bioinformatic-AI3608-Program/assets/109678718/4f362a09-a2cb-4c88-96a3-7624cbfcfa35)
The pipeline requires:
- a csv file : link SRR ID to the readable and understandable name you set
- a json file : record software path, data type(single end/pair end) and the normalization method
- raw fastq files : SRRxxxxxx_1/2.fastq.gz for pair end data and SRRxxxxxx.fastq.gz for single end data

After running the snakemake pipeline, we get the expression martix for downstream analysis.

#### 2. Identified and classified differentially expressing genes
Then We use R package DESeq2 to identify the differentially expressing genes and separately classify upregulated/downregulated genes into 3 categories: Early, Constant, Late. Please refer to the paper for definitions and methods.

We also analysed the samples' correlation. The results are shown below:

![image](https://github.com/Eric-Y-S/Bioinformatic-AI3608-Program/assets/109678718/fd86e389-48a7-4a09-a031-97e338c358ae)

#### 3. Analysed calling results and successfully reproduced key findings
We plot the expression patterns of genes in different stages and conditions. Our result closely mirrors the paper's, although they appear slightly more disordered. This is probably due to our more lenient quality control standard.(We do not filter "genes with a s.d. higher than 50% of the mean expression" as we did not understand the meaning of this requirement)
![image](https://github.com/Eric-Y-S/Bioinformatic-AI3608-Program/assets/109678718/f6e487b8-05ec-43c1-8ce2-650414a46bb6)

The number of genes we identified is shown below. We broadly replicate the result shown in paper but get more genes for almost all categories.
![image](https://github.com/Eric-Y-S/Bioinformatic-AI3608-Program/assets/109678718/6f0d1604-c52e-4a6a-bb29-e7e84ed64423)

## Acknowledgement
Thanks to Prof. Ya Guo(郭亚) and T.A. Zhiyu Zhang(张之宇) for their guidance.

## File declaration
- `./RNApipeline` contain codes and configuration files for the pipeline.
- `./analysis` contain R scripts and data for out results.
