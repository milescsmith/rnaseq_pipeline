#! /usr/bin/env bash
snakemake \
--configfile ./config.yaml \
--cluster-config ./cluster.json \
--cluster \
"sbatch -p {cluster.partition} \
--cpus-per-task={cluster.cpus-per-task} \
--mem={cluster.mem} \
--nodes={cluster.nodes} \
--export={cluster.path},{cluster.java_opts}"  \
--jobs 60 \
--snakefile ./Snakefile \
--use-conda \
--latency-wait 60 \
--rerun-incomplete \
kallisto_with_qc