#! /bin/bash -l

#SBATCH -J pipeline_runner
#SBATCH -o runner.log
#SBATCH --mail-user=miles-smith@omrf.org
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=16
#SBATCH --partition=serial
#SBATCH --cpus-per-task=4

# module load cellranger
conda activate snakemake
snakemake \
        --configfile /s/guth-aci/narch_advanta_sequencing/scripts/rnaseq_pipeline/config.yaml \
        --cluster-config /s/guth-aci/narch_advanta_sequencing/scripts/rnaseq_pipeline/cluster.json \
        --cluster \
        "sbatch -p {cluster.partition} \
        --cpus-per-task={cluster.cpus-per-task} \
        --mem={cluster.mem} \
        --nodes={cluster.nodes} \
        --export={cluster.path},{cluster.java_opts}"  \
        --jobs 30 \
        --snakefile /s/guth-aci/narch_advanta_sequencing/scripts/rnaseq_pipeline/Snakefile \
        --use-conda \
        --latency-wait 60 \
        --rerun-incomplete \
        --use-singularity \
        --singularity-args "--bind /s/guth-aci/:/s/guth-aci/ --bind /Volumes/guth_aci_informatics/references:/Volumes/guth_aci_informatics/references" \
        can_fish
