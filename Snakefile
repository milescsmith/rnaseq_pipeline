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
from os import getcwd, environ
import glob
import re

#configfile: "config.yaml"
# this is entirely because it does not seem to be possible to concatenate values
# in either YAML or JSON

USER = environ.get("USER")
RAW_DATA_DIR = "/s/guth-aci/20190605-cgc-run-archives/raw_data"
OUT_DIR = "/s/guth-aci/20190605-cgc-run-archives/output"

THREADS = 4

# The list of samples to be processed
SAMPLES = glob.glob(f'{RAW_DATA_DIR}**/*.fastq.gz', recursive=False)
#SAMPLES = GS.glob_wildcards(RAW_DATA_DIR + "{samplename}.fastq.gz")
SAMPLES = [sample.replace(f'{RAW_DATA_DIR}/','').replace('.fastq.gz','') for sample in SAMPLES]
SAMPLES = [('_').join(sample.split('_')[:-2]) for sample in SAMPLES]
SAMPLES = [_.split('/')[-1] for _ in SAMPLES]

rule initial_qc:
    """Use Fastqc to examine the quality of the fastqs from the CGC."""
    input:
        R1=RAW_DATA_DIR+'/{sample}_R1_001.fastq.gz',
        R2=RAW_DATA_DIR+'/{sample}_R2_001.fastq.gz'
    params:
        threads=f"--threads {THREADS}",
        outdir=OUT_DIR+"/qc/initial/{sample}"
    output:
        html=OUT_DIR+'/qc/initial/{sample}/{sample}_fastqc.html',
        zip=OUT_DIR+'/qc/initial/{sample}/{sample}_fastqc.zip'
    # singularity:
    #     "docker://milescsmith/fastqc"
    log:
        "logs/fastqc/fastqc_{sample}.log"
    wrapper:
        "0.35.0/bio/fastqc"
    # shell:
    #     """
    #     fastqc {params.threads} \
    #         --quiet \
    #         --outdir {params.outdir} \
    #         {input.R1} {input.R2}
    #     """

rule initial_qc_all:
    """Target rule to run just the inital Fastqc"""
    input:
        expand(OUT_DIR+"/qc/initial/{sample}/{sample}_fastqc.html", sample=SAMPLES)
    version: 2.0

rule run_initial_multiqc:
    input:
        expand(OUT_DIR+"/qc/initial/{sample}/{sample}_fastqc.html", sample=SAMPLES)
    output:
        name=OUT_DIR+"/multiqc_initial_report.html"
    log:
        OUT_DIR+"/logs/multiqc.html"
    version: 1.0
    singularity:
        "docker://ewels/multiqc"
    shell:
        "multiqc --force {params.proj_dir} -n {output}"

rule just_initial_qc:
    input: OUT_DIR+"/multiqc_initial_report.html"
    version: 1.0