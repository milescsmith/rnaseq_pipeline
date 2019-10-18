#! /usr/bin/env bash
snakemake \
--configfile ./star_pipeline/config.yaml \
--cluster-config ./star_pipeline/cluster.json \
--cluster \
"sbatch -p {cluster.partition} \
--cpus-per-task={cluster.cpus-per-task} \
--mem={cluster.mem} \
--nodes={cluster.nodes} \
--export={cluster.path},{cluster.java_opts}"  \
--jobs 30 \
--snakefile ./star_pipeline/Snakefile \
--use-conda \
--use-singularity \
--latency-wait 60 \
--rerun-incomplete \
featureCounts_all
