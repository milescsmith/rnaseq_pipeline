#! /usr/bin/env bash
snakemake \
--configfile $PWD/config.yaml \
--cluster-config $HOME/workspace/rnaseq_pipeline/cluster.json \
--cluster \
"sbatch -p {cluster.partition} \
--cpus-per-task={cluster.cpus-per-task} \
--mem={cluster.mem} \
--nodes={cluster.nodes} \
--export={cluster.path},{cluster.java_opts}"  \
--jobs 32 \
--snakefile $HOME/workspace/rnaseq_pipeline/Snakefile \
--use-conda \
--latency-wait 30 \
--rerun-incomplete kallisto_with_qc
