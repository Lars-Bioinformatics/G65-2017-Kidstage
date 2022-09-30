# conda install mamba
# mamba install -c bioconda snakemake

snakemake -s /work/G65-2017-Kidstage/scripts/ngCGH.smk --cores 63 --use-conda --conda-frontend mamba