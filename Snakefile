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
SALMON_INDEX = SEQUENCES_DIR + config["SALMON_INDEX"]
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
        "0.38.0/bio/fastqc"

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

rule kallisto:
    """Psuedoalign sequences using Kallisto. MUCH faster than STAR and 
    I"m not convinced that STAR is any better at the alignment."""
    input:
        fq1=RESULTS_DIR+"/trimmed/{sample}.R1.fq.gz",
        GTF=GTF
    output:
        RESULTS_DIR+"/kallisto/{sample}/abundance.h5",
        RESULTS_DIR+"/kallisto/{sample}/abundance.tsv",
        RESULTS_DIR+"/kallisto/{sample}/run_info.json"
    params:
        index=KALLISTO_INDEX,
        threads=THREADS,
        out_dir=RESULTS_DIR+"/kallisto/{sample}/"
    log:
        LOG_DIR+"/kallisto/kallisto_{sample}.log"
    #singularity:
    #    "docker://milescsmith/kallisto"
    version: 1.2
    shell:
        """
        kallisto quant \
            --threads={params.threads} \
            --output-dir={params.out_dir} \
            --index={params.index} \
            --single \
            --fragment-length=76 \
            --sd=5 \
            --bootstrap-samples=12 \
            --gtf {input.GTF} \
            --bias {input.fq1}
        """
rule kallisto_quant_all:
    """Target rule to force alignement of all the samples. If aligning 
    with Kallisto, use this as the target run since Kallisto typically does 
    not make the bam files needed below."""
    input: expand(RESULTS_DIR+"/kallisto/{sample}/abundance.h5", sample=SAMPLES)

rule run_kallisto_multiqc:
    input:
        kallisto_results = expand(RESULTS_DIR+"/kallisto/{sample}/abundance.h5", sample=SAMPLES),
        log_files = LOG_DIR
    output:
        name=LOG_DIR+"/multiqc_kallisto_align_report.html"
    params:
        proj_dir=PROJECT_DIR
    log:
        LOG_DIR+"/multiqc.html"
    version: 1.2
    #singularity:
    #    "docker://ewels/multiqc"
    shell:
        "multiqc --force {params.proj_dir} -n {output} {input.kallisto_results} {input.log_files}"

rule kallisto_with_qc:
    input: LOG_DIR+"/multiqc_kallisto_align_report.html"
    version: 1.1


rule salmon_quant:
    """Psuedoalign sequences using Salmon. MUCH faster than STAR and 
    I"m not convinced that STAR is any better at the alignment."""
    input:
        fq1=RESULTS_DIR+"/trimmed/{sample}.R1.fq.gz",
    output:
        RESULTS_DIR+"/salmon/{sample}/quant.sf"
    params:
        index=SALMON_INDEX,
        threads=THREADS,
        out_dir=RESULTS_DIR+"/salmon/{sample}/"
    log:
        LOG_DIR+"/salmon/salmon_{sample}.log"
    version: 1.0
    shell:
        """
        salmon quant \
            -l A \
            -p {params.threads} \
            -i {params.index} \
            --seqBias \
            --gcBias \
            --validateMappings \
            --fldMean 76 \
            --fldSD 5 \
            -r {input.fq1} \
            -o {params.out_dir} \
        """

rule salmon_quant_all:
    """Target rule to force alignement of all the samples. If aligning 
    with Salmon, use this as the target run since Salmon typically does 
    not make the bam files needed below."""
    input: expand(RESULTS_DIR+"/salmon/{sample}/quant.sf", sample=SAMPLES)

rule run_salmon_multiqc:
    input:
        alignment_results = expand(RESULTS_DIR+"/salmon/{sample}/quant.sf", sample=SAMPLES)
    output:
        LOG_DIR+"/multiqc_salmon_align_report.html"
    params:
        "-m fastqc",
        "-m bbmap",
        "-m salmon",
        "-ip"
    wrapper:
        "0.38.0/bio/multiqc"

rule salmon_with_qc:
    input: LOG_DIR+"/multiqc_salmon_align_report.html"
    version: 1.1

rule compress_salmon_results:
    input:
        quant=RESULTS_DIR+"/salmon/{sample}/quant.sf",
        summarized_qc=LOG_DIR+"/multiqc_salmon_align_report.html"
    output: RESULTS_DIR+"/salmon/{sample}/quant.sf.gz"
    params:
        threads=THREADS
    version: 1.0
    shell:
        """
        pigz -v -p {params.threads} {input.quant}
        """

rule can_fish:
    input: expand(RESULTS_DIR+"/salmon/{sample}/quant.sf.gz", sample=SAMPLES)