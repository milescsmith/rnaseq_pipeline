# rnaseq_pipeline

Snakemake scripts and associated resources for processing RNAseq data via Google Cloud services.

This script in particular is designed to handle data generated from KAPA Hyperprep RNA-seq libraries sequenced on Illumina Nextseq/Novaseq sequencers.  Very much still a work in progress.

Cluster job settings are handled by `cluster.json` while the specific locations of folders and files are set in `config.yaml`.