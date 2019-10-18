"""
Author: Miles Smith
Affiliation: OMRF
Aim: Snakemake workflow to process HyperPrep RNA-Seq data
Date: 2019/07/22
"""
version: 2.2

# libraries for Python functions
from pathlib import Path
from itertools import chain, combinations
from os.path import join
from os import getcwd, environ
import glob
import re

# Read in the values from the configuration file and build up the locations
# of files required in the pipeline
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
SUBREAD_INDEX = SEQUENCES_DIR + config["SUBREAD_INDEX"]
RESOURCE_DIR = REF_DIR + config["RESOURCE_DIR"]
GENOME_BUILD = config["GENOME_BUILD"]
POLY_A = RESOURCE_DIR + config["POLY_A"]
TRUSEQ_RNA = RESOURCE_DIR + config["TRUSEQ_RNA"]
TRUSEQ = RESOURCE_DIR + config["TRUSEQ"]
RRNAREF = RESOURCE_DIR + config["RRNAREF"]

USER = environ.get("USER")

THREADS = 4

# Find the fastq files to be processed
SAMPLES = glob.glob(f"{RAW_DATA_DIR}**/*.fastq.gz",
                    recursive=True)
SAMPLES = [sample.replace(f"{RAW_DATA_DIR}/","").replace(".fastq.gz","") 
           for sample
           in SAMPLES]
SAMPLES = [("_").join(sample.split("_")[:-2])
           for sample
           in SAMPLES]
SAMPLES = [_.split("/")[-1]
           for _
           in SAMPLES]

# Testing code, used when Snakemake seems unable to find the files.
# print(RAW_DATA_DIR)
# print(SAMPLES[0])

rule run_all:
    input:
        kallisto: LOG_DIR+"/multiqc_kallisto_align_report.html",
        salmon: expand(RESULTS_DIR+"/salmon/{sample}/quant.sf.gz", sample=SAMPLES)
        star_with_stringtie: LOG_DIR+"/multiqc_star_align_report.html"
        star_with_featureCounts: LOG_DIR+"/multiqc_star_align_report.html"

rule initial_qc:
    """Use Fastqc to examine the quality of the fastqs from the CGC."""
    input:
        R1=RAW_DATA_DIR+"/{sample}_R1_001.fastq.gz",
        R2=RAW_DATA_DIR+"/{sample}_R2_001.fastq.gz"
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
    version: 1.3
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
        R2=RAW_DATA_DIR+"/{sample}_R2_001.fastq.gz",
    params:
        out_dir="trimmed",
        phred_cutoff=5,
        polyA_ref=POLY_A,
        truseq_rna_adapter_ref=TRUSEQ_RNA,
        truseq_adapter_ref=TRUSEQ,
        rRNA_ref=RRNAREF
    output:
        filteredR1=RESULTS_DIR+"/trimmed/{sample}.R1.fq.gz",
        filteredR2=RESULTS_DIR+"/trimmed/{sample}.R2.fq.gz",
        wasteR1=RESULTS_DIR+"/trimmed/removed_{sample}.R1.fq.gz",
        wasteR2=RESULTS_DIR+"/trimmed/removed_{sample}.R2.fq.gz",
        # to collect metrics on how many ribosomal reads were eliminated
        contam=LOG_DIR+"/trimmed/contam_{sample}.csv"
    #singularity:
    #    "docker://milescsmith/bbmap"
    version: 2.1
    shell:
        """
        bbduk.sh \
                in={input.R1} \
                in2={input.R2} \
                outu={output.filteredR1} \
                out2={output.filteredR2} \
                outm={output.wasteR1} \
                outm2={output.wasteR2} \
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
        R2=expand(RESULTS_DIR+"/trimmed/{sample}.R2.fq.gz",
                  sample = SAMPLES),
        wasteR1=expand(RESULTS_DIR+"/trimmed/removed_{sample}.R1.fq.gz",
                       sample = SAMPLES),
        wasteR2=expand(RESULTS_DIR+"/trimmed/removed_{sample}.R2.fq.gz",
                       sample = SAMPLES),
        contam=expand(LOG_DIR+"/contam_{sample}.csv",
                      sample = SAMPLES)

rule kallisto_quant:
    """Psuedoalign sequences using Kallisto. MUCH faster than STAR and 
    I"m not convinced that STAR is any better at the alignment."""
    input:
        fq1=RESULTS_DIR+"/trimmed/{sample}.R1.fq.gz",
        fq2=RESULTS_DIR+"/trimmed/{sample}.R2.fq.gz",
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
        kallisto quant -t {params.threads} \
            -o {params.out_dir} \
            -i {params.index} \
            --rf-stranded \
            --bootstrap-samples=12 \
            --gtf {input.GTF} \
            --bias {input.fq1} {input.fq2}
        """
rule kallisto_quant_all:
    """Target rule to force alignement of all the samples. If aligning 
    with Kallisto, use this as the target run since Kallisto typically does 
    not make the bam files needed below."""
    input: expand(RESULTS_DIR+"kallisto/{sample}/abundance.h5", sample=SAMPLES)

rule run_kallisto_multiqc:
    input:
        kallisto_results = expand(RESULTS_DIR+"/kallisto/{sample}/abundance.h5", sample=SAMPLES),
        log_files = LOG_DIR
    output:
        LOG_DIR+"/multiqc_kallisto_align_report.html"
    params:
        "-m fastqc",
        "-m bbmap",
        "-m kallisto",
        "-ip"
    wrapper:
        "0.38.0/bio/multiqc"

rule kallisto_with_qc:
    input: LOG_DIR+"/multiqc_kallisto_align_report.html"
    version: 1.2

rule salmon_quant:
    """Psuedoalign sequences using Salmon. MUCH faster than STAR and 
    I"m not convinced that STAR is any better at the alignment.
    Also, Salmon results are easier to translate into gene-level counts than 
    Kallisto."""
    input:
        fq1=RESULTS_DIR+"/trimmed/{sample}.R1.fq.gz",
        fq2=RESULTS_DIR+"/trimmed/{sample}.R2.fq.gz",
    output:
        RESULTS_DIR+"/salmon/{sample}/quant.sf"
    params:
        index=SALMON_INDEX,
        threads=THREADS,
        out_dir=RESULTS_DIR+"/salmon/{sample}/"
    log:
        LOG_DIR+"/salmon/salmon_{sample}.log"
    version: 2.0
    shell:
        """
        salmon quant \
            --libtype MSR \
            -p {params.threads} \
            -i {params.index} \
            --seqBias \
            --gcBias \
            --validateMappings \
            -1 {input.fq1} \
            -2 {input.fq2} \
            -o {params.out_dir} \
        """

rule salmon_quant_all:
    """Target rule to force alignement of all the samples. If aligning 
    with Salmon, use this as the target run since Salmon typically does 
    not make the bam files needed below."""
    input: expand(RESULTS_DIR+"/salmon/{sample}/quant.sf", sample=SAMPLES)
    version: 1.0

rule run_salmon_multiqc:
    input:
        bcl2fastq_results = PROJECT_DIR+"/data/raw_data/bcls",
        alignment_results = expand(RESULTS_DIR+"/salmon/{sample}/quant.sf",
                                   sample=SAMPLES),
        fastqc_results = expand(RESULTS_DIR+"/qc/{sample}/{sample}_fastqc.html", 
                                sample=SAMPLES),
        contam=expand(LOG_DIR+"/trimmed/contam_{sample}.csv",
                      sample = SAMPLES)
    output:
        LOG_DIR+"/multiqc_salmon_align_report.html"
    params:
        "-m bcl2fastq",
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

rule star_align:
    input:
        fq1=RESULTS_DIR+"/trimmed/{sample}.R1.fq.gz",
        fq2=RESULTS_DIR+"/trimmed/{sample}.R2.fq.gz",
    output:
        RESULTS_DIR+"/star/{sample}/Aligned.sortedByCoord.out.bam"
    log:
        LOG_DIR+"/star/{sample}.log"
    params:
        threads=THREADS,
        index=STAR_INDEX,
        outdir=RESULTS_DIR+"/star/{sample}/"
    threads: THREADS
    shell:
    # The following options are based on the given ENCODE options
    # in the STAR manual
        """
        STAR \
            --outFilterType BySJout \
            --outFilterMultimapNmax 20 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --genomeDir {params.index} \
            --quantMode GeneCounts \
            --outFileNamePrefix {params.outdir}/ \
            --outSAMtype BAM SortedByCoordinate \
            --readFilesIn {input.fq1} {input.fq2}\
            --readFilesCommand gunzip -c \
            --runThreadN {threads}
        """

rule rename_star_results:
    input: RESULTS_DIR+"/star/{sample}/Aligned.sortedByCoord.out.bam"
    output: RESULTS_DIR+"/star/{sample}.bam"
    shell:
        "mv {input} {output}"

rule star_align_all:
    input:
        expand(RESULTS_DIR+"/star/{sample}.bam", sample=SAMPLES)

rule stringtie_quant:
    input:
        # merged_gtf="stringtie/merged.gtf",
        genome_gtf=GTF,
        sample_bam=RESULTS_DIR+"/star/{sample}.bam"
    output:
        gtf=RESULTS_DIR+"/stringtie/{sample}/{sample}.gtf",
        ctabs=expand(
            RESULTS_DIR+"/stringtie/{{sample}}/{name}.ctab",
            name=["i2t", "e2t", "i_data", "e_data", "t_data"]
        )
    threads: THREADS
    shell:
        """
        stringtie \
            -e \
            -B \
            -p {threads} \
            --fr \
            -G {input.genome_gtf} \
            -o {output.gtf} \
            {input.sample_bam}
        """

rule stringtie_quant_all:
    input: expand(RESULTS_DIR+"/stringtie/{sample}/{sample}.gtf", sample=SAMPLES)

rule run_star_stringtie_multiqc:
    input:
        fastqc_results=expand(RESULTS_DIR+"/qc/{sample}/{sample}_fastqc.html",
                              sample=SAMPLES),
        alignment_results=expand(RESULTS_DIR+"/star/{sample}.bam",
                                 sample=SAMPLES),
        contamination=expand(LOG_DIR+"/trimmed/contam_{sample}.csv",
                             sample = SAMPLES),
	    quant_results=expand(RESULTS_DIR+"/stringtie/{sample}/{sample}.gtf",
                             sample=SAMPLES)
    output:
        LOG_DIR+"/multiqc_star_align_report.html"
    params:
        "-m bcl2fastq",
        "-m fastqc",
        "-m bbmap",
        "-m star",
        "-m"
        "-ip"
    wrapper:
        "0.38.0/bio/multiqc"

rule star_with_qc:
    input: LOG_DIR+"/multiqc_star_align_report.html"
    version: 1.1

rule featureCounts:
    input: expand(RESULTS_DIR+"/star/{sample}.bam", sample=SAMPLES)
    output: 
        counts = RESULTS_DIR+"/featureCounts/counts.txt",
        summary = RESULTS_DIR+"/featureCounts/counts.txt.summary"
    threads: THREADS
    params:
        annotation = GTF
    shell:
        """
        featureCounts \
            -a {params.annotation} \
            -F GTF \
            -g gene_name \
            -p \
            -s 2 \
            -T {threads} \
            -o {output.counts} \
            -B \
            -C \
            {input}
        """

rule compress_featureCounts:
    input: RESULTS_DIR+"/featureCounts/counts.txt""
    output: RESULTS_DIR+"/featureCounts/counts.txt.gz"
    threads:
    shell: "pigz -p {threads} {input}"


rule run_star_featureCounts_multiqc:
    input:
        fastqc_results=expand(RESULTS_DIR+"/qc/{sample}/{sample}_fastqc.html",
                              sample=SAMPLES),
        alignment_results=expand(RESULTS_DIR+"/star/{sample}.bam",
                                 sample=SAMPLES),
        contamination=expand(LOG_DIR+"/trimmed/contam_{sample}.csv",
                             sample = SAMPLES),
	    quant_results=RESULTS_DIR+"/featureCounts/counts.txt.summary"
    output:
        LOG_DIR+"/multiqc_star_align_report.html"
    params:
        "-m bcl2fastq",
        "-m fastqc",
        "-m bbmap",
        "-m star",
        "-m featureCounts",
        "-ip"
    wrapper:
        "0.38.0/bio/multiqc"

rule featureCounts_all:
    input: 
        counts = RESULTS_DIR+"/featureCounts/counts.txt.gz",
        summary = LOG_DIR+"/multiqc_star_align_report.html"