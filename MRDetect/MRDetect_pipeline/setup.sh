#!/bin/bash
# Snakemake
pip install snakemake
# Java
conda install -c conda-forge openjdk==8.0.152 -y
# R
conda create --name r -y
conda install --name r -c r r -
conda activate r
R
install.packages(c("stringr","rebus","yaml"))
q("no")
# Samtools
conda create --name samtools -y
conda install --name samtools -c bioconda samtools sambamba openssl=1.0.2p bam-readcount -y
# BCFtools
conda create --name bcftools -y
conda install --name bcftools -c bioconda bcftools htslib bedtools -y
conda install --name bcftools -c r r -y
conda activate bcftools
R
install.packages(c("stringr","data.table"), repos = c("https://mirrors.dotsrc.org/cran/", "http://cran.us.r-project.org"))
q("no")
# bam-readcount
conda create --name readcount -y
conda install --name readcount -c bioconda  bam-readcount -y
curl -L https://github.com/genome/bam-readcount/archive/master.zip -o scripts/master.zip
unzip scripts/master.zip -d scripts/
cd scripts/bam-readcount-master/
cmake .
make
cd ../..
# BEDOPS
curl -L https://github.com/bedops/bedops/archive/refs/heads/master.zip -o scripts/master.zip
unzip scripts/master.zip -d scripts/
cd scripts/bedops-master/
make
cd ../..
cp scripts/bedops-master/applications/bed/conversion/bin/convert2bed-typical scripts/bedops-master/applications/bed/conversion/bin/convert2bed
export PATH=$PATH:$PWD/scripts/bedops-master/applications/bed/conversion/bin
cp scripts/bedops-master/applications/bed/sort-bed/bin/sort-bed-typical scripts/bedops-master/applications/bed/sort-bed/bin/sort-bed
export PATH=$PATH:$PWD/scripts/bedops-master/applications/bed/sort-bed/bin
cp scripts/bedops-master/applications/bed/bedops/bin/bedops-typical scripts/bedops-master/applications/bed/bedops/bin/bedops
export PATH=$PATH:$PWD/scripts/bedops-master/applications/bed/bedops/bin
# vcf-annotation-tools
pip install vcf-annotation-tools
# vt
conda create --name vt -y
conda install --name vt -c bioconda vt -y
# leiningen
mkdir -p ~/bin/
curl -L https://raw.githubusercontent.com/technomancy/leiningen/stable/bin/lein -o ~/bin/lein
chmod +x ~/bin/lein
~/bin/lein
# bcbio-variation-recall
curl -L https://github.com/bcbio/bcbio.variation.recall/archive/master.zip -o scripts/master.zip
unzip scripts/master.zip -d scripts/
cd scripts/bcbio.variation.recall-master/
make
cd ../..
# bcbio-variation-recall dependencies
conda create --name bcbio -y
conda install --name bcbio -c bioconda -c conda-forge tabix bcftools htslib openjdk==8.0.152 -y
# GATK4
conda create --name gatk -y
conda install --name gatk -c bioconda gatk4 -y
# Strelka
conda create --name strelka -y
conda install --name strelka -c bioconda strelka -y
# VarScan
conda create --name varscan -y
conda install --name varscan -c bioconda varscan -y
# SomaticSniper
conda create --name somatic-sniper -y
conda install --name somatic-sniper -c bioconda somatic-sniper samtools=0.1.16 openssl=1.0.2p bam-readcount -y
# VarDict
conda create --name vardict -y
conda install --name vardict -c bioconda -c anaconda vardict vardict-java samtools openssl=1.0.2p perl -y
# MuSE
conda create --name muse -y
conda install --name muse -c bioconda muse -y
# LoFreq
conda create --name lofreq -y
conda install --name lofreq -c bioconda lofreq -y
# Scalpel
conda create --name scalpel -y
conda install --name scalpel -c bioconda scalpel -y
# SomaticSeq
conda create --name somaticseq -y
conda install --name somaticseq -c bioconda python=3.6.3 bedtools -y
conda activate somaticseq
pip install pysam numpy scipy pandas xgboost regex
curl -L https://github.com/bioinform/somaticseq/archive/master.zip -o scripts/master.zip
unzip scripts/master.zip -d scripts/
cd scripts/somaticseq-master/
./setup.py install
cd ../..
conda deactivate
# VCFlib and VCFtools
conda create --name vcf -y
conda install --name vcf -c conda-forge -c bioconda -c defaults vcflib vcftools -y
# VCF-kit
conda create --name vcf-kit python -y
conda activate vcf-kit
pip install numpy==1.19.5
pip install VCF-kit
conda deactivate
# MRDetect deps
conda create --name mrdetect python=2.7 -y 
conda run --name mrdetect pip install Cython==0.29.15 numpy==1.16.4 pandas==0.24.2 pysam==0.15.4 tqdm==4.43.0 scipy scikit-learn
# Varscan2cna
conda create --name varscan2cna -c conda-forge -c bioconda varscan -y 
conda activate varscan2cna
pip install numpy
conda deactivate
