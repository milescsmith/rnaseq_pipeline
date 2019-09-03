#! /usr/bin/env bash
snakemake \
--configfile ./pipeline/config.yaml \
--cluster-config ./pipeline/cluster.json \
--cluster \
"sbatch -p {cluster.partition} \
--cpus-per-task={cluster.cpus-per-task} \
--mem={cluster.mem} \
--nodes={cluster.nodes} \
--export={cluster.path},{cluster.java_opts}"  \
--jobs 30 \
--snakefile ./pipeline/Snakefile \
--use-conda \
--use-singularity \
--latency-wait 60 \
--rerun-incomplete \
can_fish
