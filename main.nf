/*
 * Copyright (c) 2019, Oklahoma Medical Research Foundation (OMRF).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * This Source Code Form is "Incompatible With Secondary Licenses", as
 * defined by the Mozilla Public License, v. 2.0.
 *
 */


// RNAseq pipeline for processing Novaseq runs using Google Life Sciences
// pipelines.

// Author: Miles Smith <miles-smith@omrf.org>

// Originally derived from the nextflow-io/rnaseq-nf proof of concept by 
// Paolo Di Tommaso <paolo.ditommaso@gmail.com>, Emilio Palumbo
//  <emiliopalumbo@gmail.com>, and Evan Floden <evanfloden@gmail.com>


params.bucket                   = "gs://memory-beta"
params.reads_directory          = "test-data"
params.reads_pattern            = "*_S*_L00*_R{1,2}_00*.fastq.gz"
params.reads                    = "${params.bucket}/${params.reads_directory}/${params.reads_pattern}"
params.outdir                   = "${params.bucket}/results"
params.transcriptome            = "${params.bucket}/references/genomic/homo_sapiens/sequences/GRCh38.primary_assembly.genome.fa"
params.miscellaneous_sequences  = "${params.bucket}/references/miscellaneous"
params.polyA_ref                = "${params.miscellaneous_sequences}/polyA.fa.gz"
params.truseq_adapter_ref       = "${params.miscellaneous_sequences}/truseq.fa.gz"
params.truseq_rna_adapter_ref   = "${params.miscellaneous_sequences}/truseq_rna.fa.gz"
params.rRNA_ref                 = "${params.miscellaneous_sequences}/human_ribosomal.fa"
params.salmon_index             = "${params.bucket}/references/genomic/homo_sapiens/indices/salmon/ensembl97"

log.info """\
 G O O G L E - R N A S E Q - N F   P I P E L I N E
 ===================================
 reads        : ${params.reads}
 outdir       : ${params.outdir}
 """

polyA_reference               = file( params.polyA_ref )
truseq_adapters_reference     = file( params.truseq_adapter_ref )
truseq_rna_adapters_reference = file( params.truseq_rna_adapter_ref )
rRNA_reference                = file( params.rRNA_ref )

Channel
    .fromFilePairs( params.reads, checkExists:true )
    .into { qc_read_pairs; read_pairs}

salmon_index_ch = Channel.fromPath( params.salmon_index )

process fastqc {
    tag "FASTQC on $sample"
    publishDir "${params.outdir}/qc", mode: "copy", overwrite: true

    input:
    set sample, file(reads) from qc_read_pairs

    output:
    file "*_fastqc.{html,zip}" into fastqc_channel


    script:
    """
    fastqc -t ${task.cpus} -f fastq -q ${reads}
    """
}

process bbduk_trim {
    tag "Trimming $sample"
    publishDir "${params.outdir}/trimmed", mode: 'copy', pattern: '*.fq.gz', overwrite: true
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.csv', overwrite: true
    machineType 'n1-standard-8'
    cpus 8

    input:
    set sample, file(reads) from read_pairs
    file polyA_reference
    file truseq_adapters_reference
    file truseq_rna_adapters_reference
    file rRNA_reference

    output:
    file "*.trimmed.fq.gz" into trimmed_reads_salmon_channel
    file "*.trimmed.fq.gz" into trimmed_reads_star_channel
    file "*.trimmed.fq.gz" into trimmed_reads_qc_channel
    file "*.waste.fq.gz" into waste_channel
    file "*.contamination.csv" into contamination_channel
    val sample into sample_name_for_salmon_channel
    val sample into sample_name_for_star_channel
    

    script:
    """
    bbduk.sh \
        in=${reads[0]} \
        in2=${reads[1]} \
        outu=${sample}.R1.trimmed.fq.gz \
        out2=${sample}.R2.trimmed.fq.gz \
        outm=${sample}.R1.waste.fq.gz \
        outm2=${sample}.R2.waste.fq.gz \
        ref=${polyA_reference},${truseq_adapters_reference},${truseq_rna_adapters_reference},${rRNA_reference} \
        stats=${sample}.contamination.csv \
        statscolumns=3 \
        k=13 \
        ktrim=r \
        useshortkmers=t \
        mink=5 \
        qtrim=r \
        trimq=10 \
        minlength=20 \
        threads=${task.cpus}
    """
}

process salmon_align {
    tag "Aligning $sample"
    publishDir "${params.outdir}/aligned/salmon", mode: 'copy', overwrite: true
    machineType 'n1-highmem-8'
    cpus 8

    input:
    file salmon_index from salmon_index_ch.collect()
    val sample from sample_name_for_salmon_channel
    file reads name "*.R?.trimmed.fq.gz" from trimmed_reads_salmon_channel
    
    output:
    file "${sample}" into aligned_channel

    script:
    """
    salmon quant \
        --libType MSR \
        --threads {task.cpus} \
        --index {salmon_index} \
        --seqBias \
        --gcBias \
        --validateMappings \
        --mates1 ${reads[0]} \
        --mates2 ${reads[1]} \
        --output ${sample} \
    """
}

process star {
    tag "Aligning $sample"
    publishDir "${params.outdir}/aligned/salmon", mode: 'copy', overwrite: true
    machineType 'n1-highmem-8'
    cpus 8

    input:
    file star_index from star_index_ch.collect()
    val sample from sample_name_for_star_channel
    file reads name "*.R?.trimmed.fq.gz" from trimmed_reads_star_channel
    
    output:
    file "${sample}" into aligned_channel

    script:
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
            --genomeDir ${star_index} \
            --quantMode GeneCounts \
            --outFileNamePrefix ${sample} \
            --outSAMtype BAM SortedByCoordinate \
            --readFilesIn ${reads[0]} ${reads[1]} \
            --readFilesCommand gunzip -c \
            --runThreadN ${task.cpus}
        """
}

process multiqc {
    publishDir params.outdir, mode:'copy'

    input:
    file "${sample}" from aligned_channel.collect()
    file "*_fastqc.zip" from fastqc_channel.collect()
    file "*.trimmed.fq.gz" from trimmed_reads_qc_channel.collect()
    file "* .waste.fq.gz" from waste_channel.collect()
    file "${sample}.contamination.csv" from contamination_channel.collect()

    output:
    file('multiqc_report.html') optional true

    script:
    """multiqc -m bcl2fastq -m fastqc -m bbmap -m salmon -ip -v .
    """

}

// process compress_results{

//     script:
//     """
//     pigz -v -p {params.threads} {input.quant}
//     """
// }

workflow.onComplete {
	log.info ( workflow.success ? "\nDone!\n" : "Oops .. something went wrong" )
}
