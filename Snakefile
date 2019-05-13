"""
Author: Miles Smith
Affiliation: OMRF
Aim: Snakemake workflow to process HyperPrep RNA-Seq data
Date: 2019/02/28
"""
version: 2.1

from pathlib import Path
from itertools import chain, combinations
from os.path import join
from os import getcwd
import glob
import re
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
GS = GSRemoteProvider()

#configfile: "config.yaml"
# this is entirely because it does not seem to be possible to concatenate values
# in either YAML or JSON
BASE_DIR = config["BASE_DIR"]
SOURCE_DIR = BASE_DIR + config["SOURCE_DIR"]
PROJECT_DIR = SOURCE_DIR + config["PROJECT_DIR"]
RAW_DATA_DIR = PROJECT_DIR + config["RAW_DATA_DIR"]
OUT_DIR = PROJECT_DIR + config['OUT_DIR']
REF_DIR = BASE_DIR + config["REF_DIR"]
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

THREADS = 4

# The list of samples to be processed
SAMPLES = glob.glob(f'{RAW_DATA_DIR}**/*.fastq.gz', recursive=False)
#SAMPLES = GS.glob_wildcards(RAW_DATA_DIR + "{samplename}.fastq.gz")
SAMPLES = [sample.replace(f'{RAW_DATA_DIR}/','').replace('.fastq.gz','') for sample in SAMPLES]
SAMPLES = [('_').join(sample.split('_')[:-2]) for sample in SAMPLES]
SAMPLES = [_.split('/')[-1] for _ in SAMPLES]
print(f"raw_data_dir = {RAW_DATA_DIR}")
rule initial_qc:
    """Use Fastqc to examine the quality of the fastqs from the CGC."""
    input:
        R1=RAW_DATA_DIR+'/{sample}_R1_001.fastq.gz',
        R2=RAW_DATA_DIR+'/{sample}_R2_001.fastq.gz'
    params:
        threads=f"--threads {THREADS}",
        outdir=OUT_DIR+"/qc/initial/{sample}"
    output:
        html1=OUT_DIR+'/qc/initial/{sample}/{sample}_R1_001_fastqc.html',
        zip1=OUT_DIR+'/qc/initial/{sample}/{sample}_R1_001_fastqc.zip',
        html2=OUT_DIR+'/qc/initial/{sample}/{sample}_R2_001_fastqc.html',
        zip2=OUT_DIR+'/qc/initial/{sample}/{sample}_R2_001_fastqc.zip'
    log:
        "logs/fastqc/fastqc_{sample}.log"
    shell:
        """
        fastqc {params.threads} \
            --quiet \
            --outdir {params.outdir} \
            {input.R1} {input.R2}
        """

rule initial_qc_all:
    """Target rule to run just the inital Fastqc"""
    input: 
        R1=expand("qc/initial/{sample}/{sample}_R1_001_fastqc.html", sample=SAMPLES),
        R2=expand("qc/initial/{sample}/{sample}_R2_001_fastqc.html", sample=SAMPLES)
    version: 2.0

rule perfom_trimming:
    """Use BBmap to trim known adaptors, low quality reads, and polyadenylated sequences and filter out ribosomal reads"""
    input:
        R1=RAW_DATA_DIR+'/{sample}_R1_001.fastq.gz',
        R2=RAW_DATA_DIR+'/{sample}_R2_001.fastq.gz',
        wait=OUT_DIR+'/qc/initial/{sample}/{sample}_R1_001_fastqc.zip'
    params:
        out_dir='trimmed',
        phred_cutoff=5,
        polyA_ref=POLY_A,
        truseq_rna_adapter_ref=TRUSEQ_RNA,
        truseq_adapter_ref=TRUSEQ,
        rRNA_ref=RRNAREF
    output:
        filteredR1='trimmed/{sample}.R1.fq.gz',
        filteredR2='trimmed/{sample}.R2.fq.gz',
        wasteR1='trimmed/removed_{sample}.R1.fq.gz',
        wasteR2='trimmed/removed_{sample}.R2.fq.gz',
        contam='logs/contam_{sample}.csv' # to collect metrics on how many ribosomal reads were eliminated
    #singularity:
    #    "docker://milescsmith/bbmap"
    version: 2.0
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
        R1=expand('trimmed/{sample}.R1.fq.gz', sample = SAMPLES),
        R2=expand('trimmed/{sample}.R2.fq.gz', sample = SAMPLES),
        wasteR1=expand('trimmed/removed_{sample}.R1.fq.gz', sample = SAMPLES),
        wasteR2=expand('trimmed/removed_{sample}.R2.fq.gz', sample = SAMPLES),
        contam=expand('logs/contam_{sample}.csv', sample = SAMPLES)

rule kallisto:
    """Psuedoalign sequences using Kallisto. MUCH faster than STAR and I'm not convinced that STAR is any better at the alignment."""
    input:
        fq1='trimmed/{sample}.R1.fq.gz',
        fq2='trimmed/{sample}.R2.fq.gz',
        GTF=GTF
    output: 
        'kallisto/{sample}/abundance.h5',
        'kallisto/{sample}/abundance.tsv',
        'kallisto/{sample}/run_info.json'
    params:
        index=KALLISTO_INDEX,
        threads=THREADS,
        out_dir='kallisto/{sample}/'
    log:
        "logs/kallisto/kallisto_{sample}.log"
    #singularity:
    #    "docker://milescsmith/kallisto"
    version: 1.0
    shell:
        "kallisto quant -t {params.threads} -o {params.out_dir} -i {params.index} --rf-stranded --genomebam --gtf {input.GTF} --bias {input.fq1} {input.fq2}"

rule kallisto_quant_all:
    """Target rule to force alignement of all the samples. If aligning with Kallisto, use this as the target run since Kallisto typically does not make the bam files needed below."""
    input: expand("kallisto/{sample}/abundance.h5", sample=SAMPLES)

rule run_kallisto_multiqc:
    input:
        expand("kallisto/{sample}/abundance.h5", sample=SAMPLES)
    output:
        name="multiqc_kallisto_align_report.html"
    params:
        proj_dir=PROJECT_DIR
    log:
        "logs/multiqc.html"
    version: 1.1
    #singularity:
    #    "docker://ewels/multiqc"
    shell:
        "multiqc --force {params.proj_dir} -n {output}"

rule kallisto_with_qc:
    input: "multiqc_kallisto_align_report.html"
    version: 1.1


