#! /bin/bash -l

#SBATCH -J demux
#SBATCH -o 3_5_rerun_mixed_mkfastq.log
#SBATCH --mail-user=kevin-thomas@omrf.org
#SBATCH --mail-type=ALL
#SBATCH --mem=128
#SBATCH --partition=serial
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16

unset _JAVA_OPTIONS
nextflow rnaseq.nf \
    --project /s/guth-aci/narch_advanta_sequencing/ \
    --raw_fastqs /s/guth-aci/narch_advanta_sequencing/data/raw_data/fastqs_advanta/advanta \
    --with-dag flowchart.pdf \
    --with-report narch_advanta.html \
    -resume