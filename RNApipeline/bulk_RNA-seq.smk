#!/usr/bin/env python

import re
import subprocess
import pandas as pd

# necessary configuration and information files #
configfile: "config/bulk_RNA-seq.json"
manifest = pd.read_csv("config/bulk_RNA-seq_data.csv",header=0,comment="#",dtype={"original":str,"new":str})
#gene_id2name_list = pd.read_csv("config/gene_id2name.csv",header=0,comment="#")

# get parameters #
FASTP_PATH = config.get("fastp_path", "fastp")
STAR_PATH = config.get("star_path", "STAR")
BEDTOOLS_PATH = config.get("bedtools_path", "bedtools")
SAMTOOLS_PATH = config.get("samtools_path", "samtools")
WD = config.get("wd")
SE_PE = config.get("SE/PE") # select one from SE/PE
GTF_PATH = config.get("gtf_path")
GENOME_PATH = config.get("genome_path")
NOR_METHOD = config.get("normalization_method") # select one from RPKM/FPKM/TPM/NONE

# directory explanation #
# downloaded raw data under ./raw
# renamed and cleaned data under ./clean
# alignment data under ./align



localrules: all,rename_file

rule all:
    input:
        done = expand("raw/{name}.done",name = manifest["new"].tolist()),
        html = expand("clean/{name}.html",name = manifest["new"].tolist()),
        bam = expand("align/{name}/{name}_Aligned.sortedByCoord.out.bam",name = manifest["new"].tolist()),
        final = "results/samples_{}.csv".format(NOR_METHOD)


# link quality_control map #
dic = {}
for i in range(manifest.shape[0]):
    dic[manifest.iloc[i]["original"]] = manifest.iloc[i]["new"]
    dic[manifest.iloc[i]["new"]] = manifest.iloc[i]["original"]

def get_ori_name(wildcards,typef):
    return "{}/raw/{}_{}.fastq.gz".format(WD,dic[wildcards.name],typef)


rule rename_file:
    input:
        ori1 = lambda wildcards: get_ori_name(wildcards,"1"),
        ori2 = lambda wildcards: get_ori_name(wildcards,"2")
    output:
        done = "raw/{name}.done"
    params:
        renamed1 = "raw/{name}_R1.fastq.gz",
        renamed2 = "raw/{name}_R2.fastq.gz"
    threads:
        1
    shell:"""
        if [ ! -f {params.renamed1} ]; then
            ln -f -s {input.ori1} {params.renamed1}
        fi
        if [ ! -f {params.renamed2} ]; then
            ln -f -s {input.ori2} {params.renamed2}
        fi
        touch {output.done}
    """


rule qc:
    input:
        done = rules.rename_file.output.done
    params:
        r1 = "raw/{name}_R1.fastq.gz",
        r2 = "raw/{name}_R2.fastq.gz"
    output:
        clean_r1 = "clean/{name}_R1_clean.fastq.gz",
        clean_r2 = "clean/{name}_R2_clean.fastq.gz",
        html = "clean/{name}.html"
    shell:"""
        if [ ! -d "clean" ];then
             mkdir clean
        fi
        if [ ! -f {output.clean_r1} ];then
            {FASTP_PATH} -i {params.r1} \\
            -I {params.r2} \\
            -o {output.clean_r1} \\
            -O {output.clean_r2} \\
            -q 20 \\
            --length_required=50 \\
            -h {output.html}
        else
            touch {output.clean_r1}
            touch {output.clean_r2}
            touch {output.html}
        fi        
    """


rule star:
    input:
        clean_r1 = "clean/{name}_R1_clean.fastq.gz",
        clean_r2 = "clean/{name}_R2_clean.fastq.gz"
    output:
        bam = "align/{name}/{name}_Aligned.sortedByCoord.out.bam"
    params:
        gtf = GTF_PATH,
        genome_path = GENOME_PATH
    threads:
        8
    shell:"""
    if [ ! -d "align" ];then
            mkdir align
    fi
    if [ ! -d "align/{wildcards.name}" ];then
            mkdir align/{wildcards.name}
    fi
    if [ ! -f {output.bam} ];then
        STAR \\
        --genomeDir {params.genome_path} \\
        --sjdbGTFfile {params.gtf} \\
        --runThreadN {threads} \\
        --readFilesCommand zcat \\
        --alignIntronMin 30 \\
        --alignIntronMax 50000 \\
        --outSJfilterReads All \\
        --outFilterMultimapNmax 1 \\
        --outFileNamePrefix align/{wildcards.name}/{wildcards.name}_ \\
        --readFilesIn {input.clean_r1} {input.clean_r2} \\
        --outSAMtype BAM SortedByCoordinate \\
        --outReadsUnmapped Fastx \\
        --outSAMmode Full \\
        --outTmpDir align/{wildcards.name}/_STARtmp
    fi
    """


# reads count and normalize #
BAM = ""
sample_list = manifest["new"].tolist()
sample_list.sort()
for i in sample_list:
    BAM += "align/{}/{}_Aligned.sortedByCoord.out.bam ".format(i,i)
BAM = BAM.strip()

rule count_reads:
    input:
        bam = expand("align/{sample}/{sample}_Aligned.sortedByCoord.out.bam",sample = manifest["new"].tolist())
    output:
        counts = "results/samples_count.tsv"
    params:
        gtf = GTF_PATH,
        bam = BAM,
        datatype = SE_PE
    shell:"""
        if [ ! -d "results" ];then
            mkdir results
        fi
        if [ ! -f "results/samples_count.tsv" ];then
            if [ params.datatype = SE ];then
                featureCounts \\
                -T 10 \\
                -g gene_name \\
                -a {params.gtf} \\
                -o {output.counts} \\
                {params.bam}
            else
                featureCounts \\
                -p \\
                --countReadPairs \\
                -T 10 \\
                -g gene_name \\
                -a {params.gtf} \\
                -o {output.counts} \\
                {params.bam}
            fi
        fi
    """

rule normalize:
    input:
        tsv = "results/samples_count.tsv"
    output:
        processed = "results/samples_{}.csv".format(NOR_METHOD)
    run:
        # column:sample   row:gene
        # read data
        ori = pd.read_csv("results/samples_count.tsv",sep='\t',comment="#")
        df = pd.DataFrame()
        df['gene'] = ori['Geneid'].copy()
        df['gene_length(kb)'] = ori['Length'].astype(int)/1000.00
        
        
        for i in range(len(sample_list)):
            df[sample_list[i]] = ori.iloc[:,6+i].astype(int).copy()

        # dont normalize
        df.to_csv("results/samples_NONE.csv")

        # normalization
        if(NOR_METHOD == 'TPM'):
            for i in sample_list:
                df[i] = df[i]/df['gene_length(kb)']
            for sample in sample_list:
                total = df[sample].sum()
                df[sample] = df[sample]/total * 1000000
            df.to_csv("results/samples_{}.csv".format(NOR_METHOD))

        if(NOR_METHOD == 'RPKM'):
            print("RPKM developing")

        if(NOR_METHOD == 'FPKM'):
            print("FPKM developing")

        if(NOR_METHOD == 'RPM'):
            for sample in sample_list:
                total = df[sample].sum()
                df[sample] = df[sample]/total * 1000000
            df.to_csv("results/samples_{}.csv".format(NOR_METHOD))
                                                            
