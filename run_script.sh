#! /usr/bin/env bash
snakemake --cluster-config cluster.json \
--cluster \
"sbatch -p {cluster.partition} \
--cpus-per-task={cluster.cpus-per-task} \
--mem={cluster.mem} \
--nodes={cluster.nodes} \
--export={cluster.path},{cluster.java_opts}"  \
--jobs 4 \
--snakefile paired_end_rnaseq.snakefile \
--use-conda \
--shadow-prefix ./ \
--rerun-incomplete kallisto_with_qc
