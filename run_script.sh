#! /usr/bin/env bash
snakemake \
--configfile $HOME/workspace/rnaseq_pipeline/config.yaml \
--cluster-config $HOME/workspace/rnaseq_pipeline/cluster.json \
--cluster \
"sbatch -p {cluster.partition} \
--cpus-per-task={cluster.cpus-per-task} \
--mem={cluster.mem} \
--nodes={cluster.nodes} \
--export={cluster.path},{cluster.java_opts}"  \
--jobs 12 \
--snakefile $HOME/workspace/rnaseq_pipeline/Snakefile \
--use-conda \
--rerun-incomplete kallisto_with_qc
