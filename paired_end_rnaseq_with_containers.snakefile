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

configfile: "config.yaml"
# this is entirely because it does not seem to be possible to concatenate values
# in either YAML or JSON
BASE_DIR = config["BASE_DIR"]
SOURCE_DIR = BASE_DIR + config["SOURCE_DIR"]
# PROJECT_DIR = getcwd()
RAW_DATA_DIR = PROJECT_DIR + config["RAW_DATA_DIR"]
OUT_DIR = PROJECT_DIR + config["OUT_DIR"]
REF_DIR = BASE_DIR + config["REF_DIR"]
SEQUENCES_DIR = REF_DIR + config["SEQUENCES_DIR"]
GTF = SEQUENCES_DIR + config["GTF"]
FASTA = SEQUENCES_DIR + config["FASTA"]
STAR_INDEX = SPECIES + config["STAR_INDEX"]
KALLISTO_INDEX = SPECIES + config["KALLISTO_INDEX"]
RESOURCE_DIR = REF_DIR + config["RESOURCE_DIR"]
GENOME_BUILD = config["GENOME_BUILD"]
POLY_A = RESOURCE_DIR + config["POLY_A"]
TRUSEQ_RNA = RESOURCE_DIR + config["TRUSEQ_RNA"]
TRUSEQ = RESOURCE_DIR + config["TRUSEQ"]
RRNAREF = RESOURCE_DIR + config["RRNAREF"]

THREADS = 8

# The list of samples to be processed
SAMPLES = glob.glob(f'{GS.remote(RAW_DATA_DIR)}**/*.fastq.gz', recursive=False)
SAMPLES = [sample.replace(f'{GS.remote(RAW_DATA_DIR)}/','').replace('.fastq.gz','') for sample in SAMPLES]
SAMPLES = [('_').join(sample.split('_')[:-2])  for sample in SAMPLES]
SAMPLES = [_.split('/')[-1] for _ in SAMPLES]

rule initial_qc:
    """Use Fastqc to examine the quality of the fastqs from the CGC."""
    input:
        R1=GS.remote(RAW_DATA_DIR+'/{sample}_R1_001.fastq.gz'),
        R2=GS.remote(RAW_DATA_DIR+'/{sample}_R2_001.fastq.gz')
    params:
        f'--threads {THREADS}'
    output:
        html='qc/initial/{sample}_fastqc.html',
        zip='qc/initial/{sample}_fastqc.zip'
    log:
        "logs/fastqc/fastqc_{sample}.log"
    wrapper: "0.31.0/bio/fastqc"

rule initial_qc_all:
    """Target rule to run just the inital Fastqc"""
    input: expand("qc/initial/{sample}_fastqc.html", sample=SAMPLES)
    version: 1.0

rule perfom_trimming:
    """Use BBmap to trim known adaptors, low quality reads, and polyadenylated sequences and filter out ribosomal reads"""
    input:
        R1=GS.remote(RAW_DATA_DIR+'/{sample}_R1_001.fastq.gz'),
        R2=GS.remote(RAW_DATA_DIR+'/{sample}_R2_001.fastq.gz'),
    	wait='qc/initial/{sample}_fastqc.zip'
    params:
        out_dir='trimmed',
        phred_cutoff=5,
        polyA_ref=GS.remote(POLY_A),
        truseq_rna_adapter_ref=GS.remote(TRUSEQ_RNA),
        truseq_adapter_ref=GS.remote(TRUSEQ),
        rRNA_ref=GS.remote(RRNAREF)
    output:
        filteredR1=GS.remote('trimmed/{sample}.R1.fq.gz')
        filteredR2=GS.remote('trimmed/{sample}.R2.fq.gz'),
        wasteR1=GS.remote('trimmed/removed_{sample}.R1.fq.gz'),
        wasteR2=GS.remote('trimmed/removed_{sample}.R2.fq.gz'),
        contam=GS.remote('trimmed/contam_{sample}.csv' # to collect metrics on how many ribosomal reads were eliminated
    singularity:
        "docker://milescsmith/bbmap"
    version: 1.0
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
        R1=GS.remote(expand('trimmed/{sample}.R1.fq.gz', sample = SAMPLES)),
        R2=GS.remote(expand('trimmed/{sample}.R2.fq.gz', sample = SAMPLES)),
        wasteR1=GS.remote(expand('trimmed/removed_{sample}.R1.fq.gz', sample = SAMPLES)),
        wasteR2=GS.remote(expand('trimmed/removed_{sample}.R2.fq.gz', sample = SAMPLES)),
        contam=GS.remote(expand('trimmed/contam_{sample}.txt', sample = SAMPLES))

# Examine the quality of the trimmed fastqs
# Typically not run.
rule trimmed_qc:
    input:
        R1='trimmed/{sample}.R1.fq.gz',
        R2='trimmed/{sample}.R2.fq.gz',
        wasteR1='trimmed/removed_{sample}.R1.fq.gz',
        wasteR2='trimmed/removed_{sample}.R2.fq.gz',
        contam='trimmed/contam_{sample}.csv'
    params:
        f'--threads {THREADS}'
    shadow:
        "shallow"
    output:
        html='qc/trimmed/{sample}_fastqc.html',
        zip='qc/trimmed/{sample}_fastqc.zip'
    wrapper: "0.31.0/bio/fastqc"
    
rule star_align:
    """STAR alignment of paired-end data using the STAR wrapper provided by http://snakemake-wrappers.readthedocs.io/"""
    input:
        fq1='trimmed/{sample}.R1.fq.gz',
        fq2='trimmed/{sample}.R2.fq.gz',
        annotation=GTF
        # wait='qc/trimmed/{sample}_fastqc.html'
    output: "star/{sample}/Aligned.sortedByCoord.out.bam"
        # see STAR manual for additional output files
        # "star/{sample}/Aligned.out.bam"
        
        # "star/{sample}/Log.final.out"
        # "star/{sample}/ReadsPerGene.out.tab"
    log:
        "logs/star/{sample}/Log.final.out"
    params:
        # path to STAR reference genome index
        index=STAR_INDEX,
        extra="--runThreadN {params.threads} \
        --genomeDir {params.index} \
        --readFilesIn {input.fq1} {input.fq2} \
        --outFilterType BySJout \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --outFilterMultiNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outSAMattrIHstart 0 \
        --outSAMattributes NH HI AS nM XS NM MD \
        --outReadsUnmapped Fastx \
        --outTmpDir {params.tempdir} \
        --outFileNamePrefix {sample}"
    threads: 8
    wrapper:
        "0.31.0/bio/star/align"

rule star_align_all:
    input:
        expand("star/{sample}/Aligned.sortedByCoord.out.bam", sample = SAMPLES)
    version: 1.0

rule kallisto:
    """Psuedoalign sequences using Kallisto. MUCH faster than STAR and I'm not convinced that STAR is any better at the alignment."""
    input:
        fq1='trimmed/{sample}.R1.fq.gz',
        fq2='trimmed/{sample}.R2.fq.gz',
        GTF=GTF
        # wait='qc/trimmed/{sample}_fastqc.html'
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
    singularity:
        "docker://milescsmith/kallisto"
    version: 1.0
    shell:
        "kallisto quant -t {params.threads} -o {params.out_dir} -i {params.index} --rf-stranded --genomebam --gtf {input.GTF} --bias {input.fq1} {input.fq2}"

rule kallisto_quant_all:
    """Target rule to force alignement of all the samples. If aligning with Kallisto, use this as the target run since Kallisto typically does not make the bam files needed below."""
    input: expand("kallisto/{sample}/abundance.h5", sample=SAMPLES)

rule sort:
    """Sort STAR-aligned sequences to allow Stringtie quantification."""
    input: 
        "star/{sample}/Aligned.sortedByCoord.out.bam"
    output: 
        'sorted/{sample}.sorted.bam'
    params:
        ''
    threads: THREADS
    wrapper:
        "0.31.0/bio/samtools/sort"

rule stringtie_quant:
    """Quantify reads at the gene level using Stringtie. May need to be adjusted should we later want to start looking at isoform-specific transcription."""
    input:
        genome_gtf=GTF,
        sample_bam="sorted/{sample}.sorted.bam"
    output:
        gtf="stringtie/{sample}/{sample}.gtf",
        ctabs=expand(
            "stringtie/{{sample}}/{name}.ctab",
            name=['i2t', 'e2t', 'i_data', 'e_data', 't_data'])
    threads: THREADS
    shell:
        "stringtie -e -B --rf -p {threads} -G {input.genome_gtf} -o {output.gtf} {input.sample_bam}"

rule stringtie_quant_all:
    input: expand("stringtie/{sample}/{sample}.gtf", sample=SAMPLES)
    version: 1.0

rule generate_count_matrices:
    """Run the python script that compiles all the count data into an expression matrix."""
    input: expand("stringtie/{sample}/{sample}.gtf", sample=SAMPLES)
    output: 
        gcm="stringtie/gene_count_matrix.csv",
        tcm="stringtie/transcript_count_matrix.csv"
    params: 
        inputdir = "stringtie/"
    version: 1.0
    shell:
        "prepDE.py -i {params.inputdir} -g {output.gcm} -t {output.tcm}"

rule quant_all:
    """Target rule to just quantify samples with STAR."""
    input: "stringtie/transcript_count_matrix.csv"
    version: 1.0

rule perform_qualimap_qc:
    """The next few rules are for extra quanity metrics.  I'm not sure most add anything other than run time."""
    input:  "sorted/{sample}.sorted.bam"
    output: 'mapped/post_mapping_qualimap/{sample}/qualimapReport.html'
    resources: mem_mb=96000
    params:
        outdir='mapped/post_mapping_qualimap/{sample}',
        gtf=GTF
    version: 1.0
    shell:
        "qualimap rnaseq -bam {input} -gtf {params.gtf} --outdir {params.outdir} --java-mem-size=8G"

rule alignmentsummarymetrics:
    input: 
        R=FASTA,
        I="sorted/{sample}.sorted.bam"
    output: "qc/alignmentsummarymetrics/{sample}/output.txt"
    version: 1.0
    shell:
        "picard CollectAlignmentSummaryMetrics R={input.R} I={input.I} O={output}"

rule featurecounts:
    input:
        gtf=GTF,
        bams=expand("sorted/{sample}.sorted.bam", sample=SAMPLES)
    output:
        "qc/featurecounts/counts.txt"
    params:
        threads=THREADS
    version: 1.0
    shell:
        "featureCounts -T {params.threads} -t exon -g gene_id -a {input.gtf} -o {output} {input.bams}"

rule samtools_index:
    input: "sorted/{sample}.sorted.bam"
    output: "sorted/{sample}.sorted.bam.bai"
    params:
        "" # optional params string
    wrapper:
        "0.31.0/bio/samtools/index"

rule deeptools_coverage:
    input:
        gtf=GTF,
        bams=expand("sorted/{sample}.sorted.bam", sample=SAMPLES),
        indices=expand("sorted/{sample}.sorted.bam.bai", sample=SAMPLES)
    output:
        pCrawcounts="qc/deeptools/plotCoverage.rawcounts.txt",
        pCplotfile="qc/deeptools/plotCoverage.plotfile.plotly.html"
    params:
        threads=THREADS
    version: 1.0
    shell:
        "plotCoverage --numberOfProcessors {params.threads} --bamfiles {input.bams} --plotFile {output.pCplotfile} --plotFileFormat plotly --outRawCounts {output.pCrawcounts}"

rule qorts:
    input:
        unalignedR1=RAW_DATA_DIR+'/{sample}_R1_001.fastq.gz',
        unalignedR2=RAW_DATA_DIR+'/{sample}_R2_001.fastq.gz',
        bam="sorted/{sample}.sorted.bam",
        gtf=GTF,
        genome=FASTA
    output: "qc/qorts/{sample}/"
    log: "logs/qorts/qorts_{sample}.log"
    version: 1.0
    shell:
        "qorts QC --generatePlots --stranded --rawfastq {input.unalignedR1},{input.unalignedR2} {input.bam} {input.gtf} {output}"

# It is easier to just run Multiqc manually
rule run_star_multiqc:
    input:
        "stringtie/gene_count_matrix.csv",
        "stringtie/transcript_count_matrix.csv",
        expand("mapped/post_mapping_qualimap/{sample}/qualimapReport.html", sample=SAMPLES),
        expand("qc/alignmentsummarymetrics/{sample}/output.txt", sample=SAMPLES),
        expand("qc/featurecounts/counts.txt", sample=SAMPLES),
        expand("sorted/{sample}.sorted.bam.bai", sample=SAMPLES),
        expand("qc/deeptools/plotCoverage.rawcounts.txt", sample=SAMPLES),
        expand("qc/deeptools/plotCoverage.plotfile.plotly.html", sample=SAMPLES)
    output:
        name="multiqc_star_align_report.html"
    log:
        "logs/multiqc.html"
    version: 1.1
    shell:
        "multiqc --force . -n {output}"

rule star_with_qc:
    input: "multiqc_star_align_report.html"
    version: 1.1

rule run_kallisto_multiqc:
    input:
        expand("kallisto/{sample}/abundance.h5", sample=SAMPLES),
    output:
        name="multiqc_kallisto_align_report.html"
    log:
        "logs/multiqc.html"
    version: 1.1
    shell:
        "multiqc --force . -n {output}"

rule kallisto_with_qc:
    input: "multiqc_kallisto_align_report.html"
    version: 1.1
