#!/usr/bin/env nextflow
/*
 * Copyright (c) 2020, Oklahoma Medical Research Foundation (OMRF).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * This Source Code Form is "Incompatible With Secondary Licenses", as
 * defined by the Mozilla Public License, v. 2.0.
 *
 */

// Nextflow pipeline for processing PacBio IsoSeq runs
// Author: Miles Smith <miles-smith@omrf.org>
// Date: 2020/09/04
// Version: 0.1.0

// File locations and program parameters

def helpMessage() {
    log.info nfcoreHeader()
    log.info """
    Usage:

      The typical command for running the pipeline is as follows:
 
      nextflow run milescsmith/nf-rnaseq \\
        --project /path/to/project \\
        --profile slurm

    Mandatory arguments:
      --project             Directory to use as a base where raw files are 
                            located and results and logs will be written
      --profile             Configuration profile to use for processing.
                            Available: slurm, gcp, standard
      
    Optional parameters:

    Reference locations:
      --species             By default, uses "chimera" for a mixed human/viral
                            index.
      --genome
      --truseq_adapter
      --truseq_rna_adapter
      --rRNAs
      --polyA
      --salmon_index
    
    Results locations:      If not specified, these will be created relative to
                            the `project` argument
                            Note: if undefined, the pipeline expects the raw BAM
                            file to be located in `/path/to/project/01_raw`
      --logs                Default: project/logs
      --raw_data            Default: project/data/raw_data
      --bcls                Default: project/data/raw_data/bcls
      --raw_fastqs          Default: project/data/raw_data/fastqs
      --results             Default: project/data/raw_data/data/results
      --qc                  Default: project/data/raw_data/data/results/qc
      --trimmed             Default: project/data/raw_data/data/results/trimmed
      --aligned             Default: project/data/raw_data/data/results/aligned

    Other:
      --help                Show this message
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

params.input = "${params.fastq}/*_S*_L00*_R{1,2}_00*.fastq.gz"

Channel
    .fromFilePairs( params.reads, checkIfExists: true, flat: true )
    .into{ raw_fastq_fastqc_reads_ch; raw_fastq_to_trim_ch }

Channel
    .fromPath( params.polyA, checkIfExists: true)
    .set{ polyA_ref }

Channel
    .fromPath( params.truseq_adapter, checkIfExists: true)
    .set{ truseq_adapter_ref }

Channel
    .fromPath( params.truseq_rna_adapter, checkIfExists: true)
    .set{ truseq_rna_adapter_ref }

Channel
    .fromPath( params.rRNAs, checkIfExists: true)
    .set{ rRNA_ref }

process initial_qc {
    //Use Fastqc to examine fastq quality
    
    tag "FastQC"
    container "registry.gitlab.com/milothepsychic/rnaseq_pipeline/fastqc:0.11.9"
    queue "highmem"
    
    publishDir "${params.qc}/${sample_id}", mode: "copy", pattern: "*.html", overwrite: true
    publishDir "${params.qc}/${sample_id}", mode: "copy", pattern: "*.zip", overwrite: true
    
    input:
        tuple val(sample_id), file(read1), file(read2) from raw_fastq_fastqc_reads_ch

    output:
         file "*.html" into fastqc_results_ch
    
    script:
        """
        fastqc \
            --noextract \
            --format fastq \
            --threads ${task.cpus} \
            ${read1} ${read2}
        """
}

process perfom_trimming {
    /* Use BBmap to trim known adaptors, low quality reads, and
       polyadenylated sequences and filter out ribosomal reads */
    
    tag "trim"
    container "registry.gitlab.com/milothepsychic/rnaseq_pipeline/bbmap:38.86"
    queue "highmem"
    
    input:
        tuple val(sample_id), file(read1), file(read2) from raw_fastq_to_trim_ch
        file polyA from polyA_ref
        file truseq_adapter from truseq_adapter_ref
        file truseq_rna_adapter from truseq_rna_adapter_ref
        file rRNAs from rRNA_ref

    output:
        file "*.trimmed.R1.fq.gz" into trimmed_read1_ch
        file "*.trimmed.R2.fq.gz" into trimmed_read2_ch
        val sample_id into trimmed_sample_name_ch
        file "*.csv" into contamination_ch
    
    script:
    """
    bbduk.sh \
        in=${read1} \
        in2=${read2} \
        outu=${sample_id}.trimmed.R1.fq.gz \
        out2=${sample_id}.trimmed.R2.fq.gz \
        outm=removed_${sample_id}.R1.fq.gz \
        outm2=removed_${sample_id}.R2.fq.gz \
        ref=${polyA},${truseq_adapter},${truseq_rna_adapter},${rRNAs} \
        stats=contam_${sample_id}.csv \
        statscolumns=3 \
        k=13 \
        ktrim=r \
        useshortkmers=t \
        mink=5 \
        qtrim=r \
        trimq=10 \
        minlength=20
    """
}

process salmon_quant {
    tag "salmon quant"
    container "docker://combinelab/salmon:1.3.0"

    input:
        val sample_id from trimmed_sample_name_ch
        file trimmed_read1 from trimmed_read1_ch
        file trimmed_read2 from trimmed_read2_ch

    output:
        file "${sample_id}/quant.sf" into pseudoquant_ch, pseudoquant_to_compress_ch
        file "${sample_id}/logs/salmon_quant.log" into pseudoquant_log_ch
        val sample_name into pseudoquant_name

    script:
    """
    salmon quant \
        --libType A \
        --threads ${task.cpus} \
        --index ${params.salmon_index} \
        --seqBias \
        --gcBias \
        --validateMappings \
        --mates1 ${trimmed_read1} \
        --mates2 ${trimmed_read2} \
        --output ${sample_id} \
    """
}

process multiqc {
    /* collect all qc metric into one report */
    
    tag "multiqc"
    container "docker://ewels/multiqc:1.9"
    
    publishDir "{params.qc}", mode: "copy", pattern: "multiqc_report.html", overwrite: true
    
    input:
        file alignment_results from pseudoquant_ch.collect()
        file fastqc_results from fastqc_results_ch.collect()
        file contam from contamination_ch.collect()
        file pseudoquant_log from pseudoquant_log_ch.collect()
    
    output:
        file "/multiqc_report.html" into multiqc_ch
    
    script:
    """
    multiqc \
        -m fastqc \
        -m bbmap \
        -m salmon \
        -ip \
        ${alignment_results} \
        ${fastqc_results} \
        ${contam} \
    """
}

process compress_salmon_results {
    /* No reason not to compress these results since tximport
       can read gzipped files */
    
    tag "compress results"
    container "registry.gitlab.com/milothepsychic/rnaseq_pipeline/pigz:2.4"

    publishDir "${params.aligned}", mode: "copy", pattern: "*.quant.sf.gz", overwrite: false

    input:
        file quant from pseudoquant_to_compress_ch
        val sample_id from pseudoquant_name
        
    output:
        file("${sample_id}/quant.sf.gz")
    
    script:
    """
    pigz -v -p ${task.cpus} ${sample_id}/${quant}
    """
}

def nfcoreHeader() {
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  milescsmith/rnaseq v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}