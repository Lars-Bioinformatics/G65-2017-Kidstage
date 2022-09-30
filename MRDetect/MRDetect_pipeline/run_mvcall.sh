#!/bin/bash
snakemake -s mvcall.smk --configfile config.yaml --cluster "sbatch --cpus-per-task {resources.cpus} --mem {resources.mem_mb}" --jobs all --rerun-incomplete
