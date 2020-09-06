"""
Author: Miles Smith
Affiliation: OMRF
Aim: Snakemake workflow to process HyperPrep RNA-Seq data
Date: 2019/11/14
"""
version: 3.0

# libraries for Python functions
from os.path import join
from os import getcwd, environ, path
import glob
import re

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

# Read in the values from the configuration file and build up the locations
# of files required in the pipeline
PROJECT = config["PROJECT"]
RAW_DATA = PROJECT + config["RAW_DATA"]
RESULTS = PROJECT + config["RESULTS"]
LOGS = PROJECT + config["LOGS"]

REFERENCES = config["REFERENCES"]
GENOMIC_FILES = REFERENCES + config["GENOMIC_FILES"]
SPECIES = GENOMIC_FILES + config["SPECIES"]

SEQUENCES = SPECIES + config["SEQUENCES"]
GTF = SEQUENCES + config["GTF"]
FASTA = SEQUENCES + config["FASTA"]

INDICES = SPECIES + config["INDICES"]
STAR_INDEX = INDICES + config["STAR_INDEX"]
SALMON_INDEX = INDICES + config["SALMON_INDEX"]
SUBREAD_INDEX = INDICES + config["SUBREAD_INDEX"]
RSEM_INDEX = INDICES + config["RSEM_INDEX"]

MISCELLANEOUS = REFERENCES + config["MISCELLANEOUS"]
POLY_A = MISCELLANEOUS + config["POLY_A"]
TRUSEQ_RNA = MISCELLANEOUS + config["TRUSEQ_RNA"]
TRUSEQ = MISCELLANEOUS + config["TRUSEQ"]
RRNAREF = MISCELLANEOUS + config["RRNAREF"]

TRANSCRIPTOME = config["TRANSCRIPTOME"]
GENOME = config["GENOME"]

TRANSCRIPTOME_VERSION = TRANSCRIPTOME.split("/")[-1]
TRANSCRIPTOME_NAME = ".".join(TRANSCRIPTOME_VERSION.split(".")[:2])
GENOME_VERSION = GENOME.split("/")[-1]

USER = environ.get("USER")

THREADS = 8

# Find the fastq files to be processed
SAMPLES = glob.glob(f"{RAW_DATA}**/*.fastq.gz",
                    recursive=True)
SAMPLES = [sample.replace(f"{RAW_DATA}/","").replace(".fastq.gz","") 
           for sample
           in SAMPLES]
SAMPLES = [("_").join(sample.split("_")[:-2])
           for sample
           in SAMPLES]
SAMPLES = [_.split("/")[-1]
           for _
           in SAMPLES]

# Testing code, used when Snakemake seems unable to find the files.

rule run_all:
    input:
        # kallisto=LOGS+"/multiqc_kallisto_align_report.html",
        salmon=expand(RESULTS+"/salmon/{sample}/quant.sf.gz", sample=SAMPLES),
        star_with_stringtie=LOGS+"/multiqc_star_stringtie_align_report.html",
        star_with_featureCounts=LOGS+"/multiqc_star_featureCounts_align_report.html"

rule initial_qc:
    """Use Fastqc to examine fastq quality."""
    input:
        R1=RAW_DATA+"/{sample}_R1_001.fastq.gz",
        R2=RAW_DATA+"/{sample}_R2_001.fastq.gz"
    params:
        outdir=RESULTS+"/qc/{sample}"
    output:
        html1=RESULTS+"/qc/{sample}/{sample}_R1_001_fastqc.html",
        html2=RESULTS+"/qc/{sample}/{sample}_R2_001_fastqc.html",
        zip1=RESULTS+"/qc/{sample}/{sample}_R1_001_fastqc.zip",
        zip2=RESULTS+"/qc/{sample}/{sample}_R2_001_fastqc.zip"
    singularity:
        "docker://registry.gitlab.com/milothepsychic/rnaseq_pipeline/fastqc:0.11.9"
    log:
        LOGS+"/fastqc/fastqc_{sample}.log"
    threads:
        THREADS
    version: 2.0
    shell:
        """
        mkdir -p {params.outdir} && \
        fastqc \
            --noextract \
            --outdir {params.outdir} \
            --format fastq \
            --threads {threads} \
            {input.R1} {input.R2}
        """

rule initial_qc_all:
    """Target rule to run just the inital Fastqc"""
    input:
        expand(RESULTS+"/qc/{sample}/{sample}_fastqc.html", 
               sample=SAMPLES)
    version: 2.0

rule perform_trimming:
    """Use BBmap to trim known adaptors, low quality reads, 
    and polyadenylated sequences and filter out ribosomal reads"""
    input:
        R1=RAW_DATA+"/{sample}_R1_001.fastq.gz",
        R2=RAW_DATA+"/{sample}_R2_001.fastq.gz",
        make_qc_run=rules.initial_qc.output
    params:
        out_dir="trimmed",
        phred_cutoff=5,
        polyA_ref=POLY_A,
        truseq_rna_adapter_ref=TRUSEQ_RNA,
        truseq_adapter_ref=TRUSEQ,
        rRNA_ref=RRNAREF
    output:
        filteredR1=RESULTS+"/trimmed/{sample}.R1.fq.gz",
        filteredR2=RESULTS+"/trimmed/{sample}.R2.fq.gz",
        wasteR1=RESULTS+"/trimmed/removed_{sample}.R1.fq.gz",
        wasteR2=RESULTS+"/trimmed/removed_{sample}.R2.fq.gz",
        # to collect metrics on how many ribosomal reads were eliminated
        contam=LOGS+"/trimmed/contam_{sample}.csv"
    singularity:
        "docker://registry.gitlab.com/milothepsychic/rnaseq_pipeline/bbmap:38.86"
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
        R1=expand(RESULTS+"/trimmed/{sample}.R1.fq.gz",
                  sample = SAMPLES),
        R2=expand(RESULTS+"/trimmed/{sample}.R2.fq.gz",
                  sample = SAMPLES),
        wasteR1=expand(RESULTS+"/trimmed/removed_{sample}.R1.fq.gz",
                       sample = SAMPLES),
        wasteR2=expand(RESULTS+"/trimmed/removed_{sample}.R2.fq.gz",
                       sample = SAMPLES),
        contam=expand(LOGS+"/contam_{sample}.csv",
                      sample = SAMPLES)

rule retrieve_source_sequences:
    """If the files to build a reference are not available, retrieve them"""
    input:
        transcriptome=FTP.remote(TRANSCRIPTOME),
        genome=FTP.remote(GENOME)
    output:
        transcriptome=SEQUENCES+"/"+TRANSCRIPTOME_VERSION,
        genome=SEQUENCES+"/"+GENOME_VERSION
    params:
    log:
    singularity:
        "docker://registry.gitlab.com/milothepsychic/rnaseq_pipeline/snakemake"
    shell:
        """
        mv {input.transcriptome} {output.transcriptome}
        mv {input.genome} {output.genome}
        """

rule build_genome_decoys:
    input:
        genome=SEQUENCES+"/"+GENOME_VERSION,
        transcriptome=SEQUENCES+"/"+TRANSCRIPTOME_VERSION
    output:
        decoys=SEQUENCES+"/"+TRANSCRIPTOME_NAME+".decoys.txt",
        gentrome=SEQUENCES+"/"+TRANSCRIPTOME_NAME+".gentrome.fa.gz"
    singularity:
        "docker://registry.gitlab.com/milothepsychic/rnaseq_pipeline/snakemake"
    shell:
        """
        grep '^>' <(zcat {input.genome}) | cut -d ' ' -f 1 > {output.decoys}
        sed -i -e 's/>//g' {output.decoys}
        cat {input.transcriptome} {input.genome} > {output.gentrome}
        """

rule build_salmon_index:
    """If the index does not exist, build it."""
    input:
        decoys=SEQUENCES+"/"+TRANSCRIPTOME_NAME+".decoys.txt",
        gentrome=SEQUENCES+"/"+TRANSCRIPTOME_NAME+".gentrome.fa.gz"
    output:
        salmon_index=directory(SALMON_INDEX)
    threads:
        THREADS
    log: LOGS+"/salmon_index"
    singularity:
        "docker://registry.gitlab.com/guthridge_informatics/salmon:1.3.0"
    version: 1.0
    shell:
        """
        salmon \
            index \
            --transcripts {input.gentrome} \
            --decoys {input.decoys} \
            --threads {threads} \
            --index {output.salmon_index} \
            --gencode
        """

rule salmon_quant:
    """Psuedoalign sequences using Salmon. MUCH faster than STAR and 
    I"m not convinced that STAR is any better at the alignment.
    Also, Salmon results are easier to translate into gene-level counts than 
    Kallisto."""
    input:
        fq1=RESULTS+"/trimmed/{sample}.R1.fq.gz",
        fq2=RESULTS+"/trimmed/{sample}.R2.fq.gz",
    output:
        quant=RESULTS+"/salmon/{sample}/quant.sf",
    params:
        index=directory(SALMON_INDEX),
        threads=THREADS,
        outdir=directory(RESULTS+"/salmon/{sample}/")
    log:
        LOGS+"/salmon/salmon_{sample}.log"
    singularity:
        "docker://combinelab/salmon:1.3.0"
    version: 4.0
    shell:
        """
        salmon quant \
            --libType A \
            --threads {params.threads} \
            --index {params.index} \
            --seqBias \
            --gcBias \
            --validateMappings \
            --mates1 {input.fq1} \
            --mates2 {input.fq2} \
            --output {params.outdir}
        """

rule salmon_quant_all:
    """Target rule to force alignement of all the samples. If aligning 
    with Salmon, use this as the target run since Salmon typically does 
    not make the bam files needed below."""
    input: expand(RESULTS+"/salmon/{sample}/quant.sf", sample=SAMPLES)
    version: 1.0

rule run_salmon_multiqc:
    input:
        alignment_results=[f"{RESULTS}/salmon/{sample}/quant.sf"
                           for sample
                           in SAMPLES
                           ],
        fastqc_results_read1=[f"{RESULTS}/qc/{sample}/{sample}_R1_001_fastqc.html"
                              for sample
                              in SAMPLES
                              ],
        fastqc_results_read2=[f"{RESULTS}/qc/{sample}/{sample}_R2_001_fastqc.html"
                              for sample
                              in SAMPLES
                              ],
        contam=[f"{LOGS}/trimmed/contam_{sample}.csv"
                for sample
                in SAMPLES
                ]
    output:
        LOGS+"/multiqc_salmon_align_report.html"
    params:
        outdir=LOGS,
        reportname="/multiqc_salmon_align_report.html"
    singularity:
        "docker://ewels/multiqc:1.9"
    shell:
        """
        multiqc \
            -m fastqc \
            -m bbmap \
            -m salmon \
            -ip \
            --outdir {params.outdir} \
            --filename {params.reportname} \
            {input.alignment_results} \
            {input.fastqc_results_read1} \
            {input.fastqc_results_read2} \
            {input.contam}
        """

rule salmon_with_qc:
    input: LOGS+"/multiqc_salmon_align_report.html"
    version: 1.1

rule compress_salmon_results:
    input:
        quant=RESULTS+"/salmon/{sample}/quant.sf",
        #summarized_qc=LOGS+"/multiqc_salmon_align_report.html"
    output: RESULTS+"/salmon/{sample}/quant.sf.gz"
    params:
        threads=THREADS
    version: 1.0
    shell:
        """
        pigz \
            -v \
            -p {params.threads} \
            {input.quant}
        """

rule can_fish:
    input:
        alignments=expand(RESULTS+"/salmon/{sample}/quant.sf.gz", sample=SAMPLES),
        qc=rules.run_salmon_multiqc.output

rule star_align:
    input:
        fq1=RESULTS+"/trimmed/{sample}.R1.fq.gz",
        fq2=RESULTS+"/trimmed/{sample}.R2.fq.gz",
    output:
        RESULTS+"/star/{sample}/Aligned.sortedByCoord.out.bam"
    log:
        LOGS+"/star/{sample}.log"
    params:
        index=STAR_INDEX,
        outdir=RESULTS+"/star/{sample}/"
    threads: THREADS
    singularity:
        "docker://registry.gitlab.com/milothepsychic/rnaseq_pipeline/star:2.7.5c"
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
            --outFileNamePrefix {params.outdir} \
            --outSAMtype BAM SortedByCoordinate \
            --readFilesIn {input.fq1} {input.fq2}\
            --readFilesCommand gunzip -c \
            --runThreadN {threads}
        """

rule rename_star_results:
    input: RESULTS+"/star/{sample}/Aligned.sortedByCoord.out.bam"
    output: RESULTS+"/star/{sample}.bam"
    shell:
        "mv {input} {output}"

rule star_align_all:
    input:
        expand(RESULTS+"/star/{sample}.bam", sample=SAMPLES)

rule stringtie_quant:
    input:
        # merged_gtf="stringtie/merged.gtf",
        genome_gtf=GTF,
        sample_bam=RESULTS+"/star/{sample}.bam"
    output:
        gtf=RESULTS+"/stringtie/{sample}/{sample}.gtf",
        ctabs=expand(
            RESULTS+"/stringtie/{{sample}}/{name}.ctab",
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
    input: expand(RESULTS+"/stringtie/{sample}/{sample}.gtf", sample=SAMPLES)

rule run_star_stringtie_multiqc:
    input:
        alignment_results=set([path.dirname(i)
                               for i
                               in expand(RESULTS+"/star/{sample}.bam",
                               sample=SAMPLES)
        ]),
        fastqc_results=set([path.dirname(j)
                            for j
                            in expand(RESULTS+"/qc/{sample}/{sample}_fastqc.html",sample=SAMPLES)
        ]),
        contamination=set([path.dirname(k)
                           for k
                           in expand(LOGS+"/trimmed/contam_{sample}.csv",
                           sample = SAMPLES)
        ]),
        quant_results=set([path.dirname(m)
                           for m
                           in expand(RESULTS+"/stringtie/{sample}/{sample}.gtf",sample=SAMPLES)
        ])
    output:
        LOGS+"/multiqc_star_stringtie_align_report.html"
    params:
        outdir=LOGS,
        reportname="/multiqc_star_stringtie_align_report.html"
    singularity:
        "docker://ewels/multiqc:1.9"
    shell:
        """
        multiqc \
            -m fastqc \
            -m bbmap \
            -m star \
            -ip \
            --outdir {params.outdir} \
            --filename {params.reportname} \
            {input.alignment_results} \
            {input.fastqc_results} \
            {input.contamination} \
            {input.quant_results}
        """

rule star_with_stringtie_qc:
    input: LOGS+"/multiqc_star_stringtie_align_report.html"
    version: 1.1

rule featureCounts:
    input: expand(RESULTS+"/star/{sample}.bam", sample=SAMPLES)
    output: 
        counts=RESULTS+"/featureCounts/counts.txt",
        summary=RESULTS+"/featureCounts/counts.txt.summary"
    threads: THREADS
    params:
        annotation = GTF
    singularity:
        "docker://registry.gitlab.com/milothepsychic/rnaseq_pipeline/subread:2.0.0"
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
    input: RESULTS+"/featureCounts/counts.txt"
    output: RESULTS+"/featureCounts/counts.txt.gz"
    threads: THREADS
    shell: "pigz -p {threads} {input}"

rule run_star_with_featureCounts_qc:
    input:
        alignment_results=set([path.dirname(i)
                               for i
                               in expand(RESULTS+"/star/{sample}.bam",
                               sample=SAMPLES)
        ]),
        fastqc_results=set([path.dirname(j)
                            for j
                            in expand(RESULTS+"/qc/{sample}/{sample}_fastqc.html",sample=SAMPLES)
        ]),
        contamination=set([path.dirname(k)
                           for k
                           in expand(LOGS+"/trimmed/contam_{sample}.csv",
                           sample = SAMPLES)
        ]),
	    quant_results=RESULTS+"/featureCounts/counts.txt.summary"
    output:
        LOGS+"/multiqc_star_featureCounts_align_report.html"
    params:
        outdir=LOGS,
        reportname="/multiqc_star_featureCounts_align_report.html"
    singularity:
        "docker://ewels/multiqc:1.9"
    shell:
        """
        multiqc \
            -m fastqc \
            -m bbmap \
            -m star \
            -m featureCounts \
            -ip \
            --outdir {params.outdir} \
            --filename {params.reportname} \
            {input.alignment_results} \
            {input.fastqc_results} \
            {input.contamination} \
            {input.quant_results}
        """

rule featureCounts_all:
    input: 
        counts=RESULTS+"/featureCounts/counts.txt.gz",
        summary=LOGS+"/multiqc_star_featureCounts_align_report.html"

rule rsem_analysis:
    version: 1.0
    input:
        R1=RESULTS+"/trimmed/{sample}.R1.fq.gz",
        R2=RESULTS+"/trimmed/{sample}.R2.fq.gz",
    output:
        genes=RESULTS+"/rsem/{sample}/{sample}.genes.results",
        isoforms=RESULTS+"/rsem/{sample}/{sample}.isoforms.results",
        alignment=RESULTS+"/rsem/{sample}/{sample}.transcript.bam",
        stats=RESULTS+"/rsem/{sample}/{sample}.stat/{sample}.cnt",
    threads:
        THREADS
    singularity:
        "docker://dceoy/rsem"
    params:
        index=RSEM_INDEX,
        prefix=RESULTS+"/rsem/{sample}/{sample}",
        mem_per_thread=8
    resources:
        mem_mb=32768
    log:
        LOGS+"/rsem/{sample}.log"
    shell:
        """
        rsem-calculate-expression \
            --ci-memory {resources.mem_mb} \
            --calc-ci \
            --calc-pme \
            --strandedness reverse \
            --star \
            --star-gzipped-read-file \
            --num-threads {threads} \
            --paired-end \
            {input.R1} \
            {input.R2} \
            {params.index} \
            {params.prefix} \
            2> {log}
        """

rule run_rsem:
    input: 
        expand(RESULTS+"/rsem/{sample}/{sample}.genes.results", sample=SAMPLES)

rule gather_rsem:
    version: 1.0
    input:
        matrices=expand(RESULTS+"/rsem/{sample}/{sample}.genes.results", sample=SAMPLES)
    output:
        RESULTS+"/rsem/counts.matrix"
    singularity:
        "docker://dceoy/rsem"
    shell:
        """
        rsem-generate-data-matrix {input.matrices} > {output}
        """

rule run_rsem_multiqc:
    input:
        alignment_results=set([path.dirname(i)
                               for i
                               in expand(RESULTS+"/star/{sample}.bam",
                               sample=SAMPLES)
        ]),
        fastqc_results=set([path.dirname(j)
                            for j
                            in expand(RESULTS+"/qc/{sample}/{sample}_fastqc.html",sample=SAMPLES)
        ]),
        contamination=set([path.dirname(k)
                           for k
                           in expand(LOGS+"/trimmed/contam_{sample}.csv",
                           sample = SAMPLES)
        ]),
	    rsem_stats=set([path.dirname(m)
                        for m
                        in expand(RESULTS+"/rsem/{sample}/{sample}.stat/{sample}.cnt",
                                                     sample=SAMPLES)
        ])
    output:
        LOGS+"/multiqc_rsem_report.html"
    params:
        outdir=LOGS,
        reportname="/multiqc_rsem_report.html"
    singularity:
        "docker://ewels/multiqc:1.9"
    shell:
        """
        multiqc \
            -m fastqc \
            -m bbmap \
            -m star \
            -m rsem \
            -ip \
            --outdir {params.outdir} \
            --filename {params.reportname} \
            {input.alignment_results} \
            {input.fastqc_results} \
            {input.contamination} \
            {input.rsem_stats}
        """

rule rsem_qc:
    version: 1.1
    input:
        QC=LOGS+"/multiqc_rsem_report.html",
        COUNTS= RESULTS+"/rsem/counts.matrix"
