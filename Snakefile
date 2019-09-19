"""
Author: Miles Smith
Affiliation: OMRF
Aim: Snakemake workflow to process HyperPrep RNA-Seq data
Date: 2019/07/22
"""
version: 2.2

from pathlib import Path
from itertools import chain, combinations
from os.path import join
from os import getcwd, environ
import glob
import re

#configfile: "config.yaml"
# this is entirely because it does not seem to be possible to concatenate values
# in either YAML or JSON
BASE_DIR = config["BASE_DIR"]
SOURCE_DIR = BASE_DIR + config["SOURCE_DIR"]
PROJECT_DIR = SOURCE_DIR + config["PROJECT_DIR"]
RAW_DATA_DIR = PROJECT_DIR + config["RAW_DATA_DIR"]
RESULTS_DIR = PROJECT_DIR + config["RESULTS_DIR"]
LOG_DIR = PROJECT_DIR + config["LOG_DIR"]
REF_DIR = config["REF_DIR"]
SEQUENCES_DIR = REF_DIR + config["SEQUENCES_DIR"]
GTF = SEQUENCES_DIR + config["GTF"]
FASTA = SEQUENCES_DIR + config["FASTA"]
STAR_INDEX = SEQUENCES_DIR + config["STAR_INDEX"]
KALLISTO_INDEX = SEQUENCES_DIR + config["KALLISTO_INDEX"]
RESOURCE_DIR = REF_DIR + config["RESOURCE_DIR"]
GENOME_BUILD = config["GENOME_BUILD"]
POLY_A = RESOURCE_DIR + config["POLY_A"]
TRUSEQ_RNA = RESOURCE_DIR + config["TRUSEQ_RNA"]
TRUSEQ = RESOURCE_DIR + config["TRUSEQ"]
RRNAREF = RESOURCE_DIR + config["RRNAREF"]

USER = environ.get("USER")

THREADS = 4

print(f"Checking for samples in {RAW_DATA_DIR}")
# The list of samples to be processed
SAMPLES = glob.glob(f"{RAW_DATA_DIR}**/*.fastq.gz",
                    recursive=False)
#SAMPLES = GS.glob_wildcards(RAW_DATA_DIR + "{samplename}.fastq.gz")
SAMPLES = [sample.replace(f"{RAW_DATA_DIR}/","").replace(".fastq.gz","") 
           for sample
           in SAMPLES]
SAMPLES = [("_").join(sample.split("_")[:-2])
           for sample
           in SAMPLES]
SAMPLES = [_.split("/")[-1]
           for _
           in SAMPLES]

print(f"Found {len(SAMPLES)} samples")

rule initial_qc:
    """Use Fastqc to examine the quality of the fastqs from the CGC."""
    input:
        R1=RAW_DATA_DIR+"/{sample}_R1_001.fastq.gz",
    params:
        threads=f"--threads {THREADS}",
        #outdir=RESULTS_DIR+"/qc/{sample}"
    output:
        html=RESULTS_DIR+"/qc/{sample}/{sample}_fastqc.html",
        zip=RESULTS_DIR+"/qc/{sample}/{sample}_fastqc.zip"
    # singularity:
    #     "docker://milescsmith/fastqc"
    log:
        LOG_DIR+"/fastqc/fastqc_{sample}.log"
    wrapper:
        "0.35.2/bio/fastqc"

rule initial_qc_all:
    """Target rule to run just the inital Fastqc"""
    input:
        expand(RESULTS_DIR+"/qc/{sample}/{sample}_fastqc.html", 
               sample=SAMPLES)
    version: 2.0

rule perfom_trimming:
    """Use BBmap to trim known adaptors, low quality reads, 
    and polyadenylated sequences and filter out ribosomal reads"""
    input:
        R1=RAW_DATA_DIR+"/{sample}_R1_001.fastq.gz",
        wait=RESULTS_DIR+"/qc/{sample}/{sample}_fastqc.html"
    params:
        out_dir="trimmed",
        phred_cutoff=5,
        polyA_ref=POLY_A,
        truseq_rna_adapter_ref=TRUSEQ_RNA,
        truseq_adapter_ref=TRUSEQ,
        rRNA_ref=RRNAREF
    output:
        filteredR1=RESULTS_DIR+"/trimmed/{sample}.R1.fq.gz",
        wasteR1=RESULTS_DIR+"/trimmed/removed_{sample}.R1.fq.gz",
        # to collect metrics on how many ribosomal reads were eliminated
        contam=LOG_DIR+"/trimmed/contam_{sample}.csv"
    #singularity:
    #    "docker://milescsmith/bbmap"
    version: 2.1
    shell:
        """
        bbduk.sh \
                in={input.R1} \
                outu={output.filteredR1} \
                outm={output.wasteR1} \
                ref={params.polyA_ref},{params.truseq_adapter_ref},{params.truseq_rna_adapter_ref},{params.rRNA_ref} \
                stats={output.contam} \
                statscolumns=3 \
                k=13 \
                ktrim=r \
                useshortkmers=t \
                mink=5 \
                qtrim=r \
                trimq=10 \
                minlength=20
        """

rule expand_trimming:
    input:
        R1=expand(RESULTS_DIR+"/trimmed/{sample}.R1.fq.gz",
                  sample = SAMPLES),
        wasteR1=expand(RESULTS_DIR+"/trimmed/removed_{sample}.R1.fq.gz",
                       sample = SAMPLES),
        contam=expand(LOG_DIR+"/trimmed/contam_{sample}.csv",
                      sample = SAMPLES)

# rule star_align:
#     input:
#         fq1=RESULTS_DIR+"/trimmed/{sample}.R1.fq.gz",
#         annotation=GTF
#     output:
#         # see STAR manual for additional output files
#         RESULTS_DIR+"/star/{sample}/Aligned.out.bam"
#         # "star/{sample}/Log.final.out"
#         # "star/{sample}/ReadsPerGene.out.tab"
#     log:
#         LOG_DIR+"/star/{sample}/Log.final.out"
#     params:
#         # path to STAR reference genome index
#         index=STAR_INDEX
#     threads: 8
#     wrapper:
#         "0.38.0/bio/star/align"

rule star_se:
    input:
        RESULTS_DIR+"/trimmed/{sample}.R1.fq.gz"
    output:
        RESULTS_DIR+"/star/{sample}/Aligned.out.bam"
        # "star/{sample}/Log.final.out",
        # "star/{sample}/ReadsPerGene.out.tab"
    log:
        LOG_DIR+"/star/{sample}.log"
    params:
        threads=THREADS,
        # path to STAR reference genome index
        index=STAR_INDEX,
        outdir=RESULTS_DIR+"/star/{sample}"
        # optional parameters
    shell:
        """
        STAR \
            --runThreadN {params.threads} \
            --genomeDir {params.index} \
            --readFilesIn {input} \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --outFilterMismatchNoverLmax 0.6 \
            --alignIntronMin 20 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --outFilterType BySJout \
            --alignSJoverhangMin 8 \
            --outSAMattrIHstart 0 \
            --outSAMattributes NH HI AS nM XS NM MD \
            --outReadsUnmapped Fastx \
            --readFilesCommand gunzip -c \
            --outFileNamePrefix {params.outdir}/
        """

rule align_all:
    input:
        expand(RESULTS_DIR+"/star/{sample}/Aligned.out.bam", sample=SAMPLES)

rule sort:
    input: 
        RESULTS_DIR+"/star/{sample}/Aligned.out.bam"
    output: 
        RESULTS_DIR+"/sorted/{sample}.sorted.bam"
    params:
        ""
    threads: THREADS
    wrapper:
        "0.36.0/bio/samtools/sort"

rule stringtie_quant:
    input:
        # merged_gtf="stringtie/merged.gtf",
        genome_gtf=GTF,
        sample_bam=RESULTS_DIR+"/sorted/{sample}.sorted.bam"
    output:
        gtf=RESULTS_DIR+"/stringtie/{sample}/{sample}.gtf",
        ctabs=expand(
            RESULTS_DIR+"/stringtie/{{sample}}/{name}.ctab",
            name=["i2t", "e2t", "i_data", "e_data", "t_data"]
        )
    threads: THREADS
    shell:
        "stringtie -e -B -p {threads} -G {input.genome_gtf} -o {output.gtf} {input.sample_bam}"

rule stringtie_all_samples:
    input: expand(RESULTS_DIR+"/stringtie/{sample}/{sample}.gtf", sample=SAMPLES)

rule generate_count_matrices:
    input: expand(RESULTS_DIR+"/stringtie/{sample}/{sample}.gtf", sample=SAMPLES)
    output: 
        gcm=RESULTS_DIR+"/gene_count_matrix.csv",
        tcm=RESULTS_DIR+"/transcript_count_matrix.csv"
    params: 
        inputdir = RESULTS_DIR+"/stringtie"
    shell:
        "prepDE.py -i {params.inputdir} -g {output.gcm} -t {output.tcm}"

rule quant_all_samples:
    input: RESULTS_DIR+"/transcript_count_matrix.csv"

rule run_post_star_multiqc:
    input:
        stringtie_results = RESULTS_DIR+"/transcript_count_matrix.csv",
        log_files = LOG_DIR
    output:
        name=LOG_DIR+"/multiqc_star_align_report.html"
    params:
        proj_dir=PROJECT_DIR
    log:
        LOG_DIR+"/multiqc.html"
    version: 1.2
    #singularity:
    #    "docker://ewels/multiqc"
    shell:
        "multiqc --force {params.proj_dir} -n {output} {input.stringtie_results} {input.log_files}"

rule star_with_qc:
    input: LOG_DIR+"/multiqc_star_align_report.html"
    version: 1.1