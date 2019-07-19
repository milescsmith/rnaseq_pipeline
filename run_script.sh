#! /usr/bin/env bash
snakemake \
--configfile /s/guth-aci/rnaseq_pipeline/config.yaml \
--cluster-config /s/guth-aci/rnaseq_pipeline/cluster.json \
--cluster \
"sbatch -p {cluster.partition} \
--cpus-per-task={cluster.cpus-per-task} \
--mem={cluster.mem} \
--nodes={cluster.nodes} \
--export={cluster.path},{cluster.java_opts}"  \
--jobs 20 \
--snakefile /s/guth-aci/rnaseq_pipeline/Snakefile \
--use-conda \
--use-singularity \
--singularity-args "-H $PWD" \
--latency-wait 30 \
--rerun-incomplete just_initial_qc
