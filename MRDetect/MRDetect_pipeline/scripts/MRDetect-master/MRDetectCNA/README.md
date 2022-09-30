## MRDetectCNA

### Dependencies 
1. GATK v3.4
2. Pandas
3. Numpy


### 1. Sample processing
This step takes as input a sample bam, CNV bed and reference fasta file to generate binned read counts.

E.g. CNA bed file
```
#chr  start end type  log2
18  31047303  34298225  NEU 0.073
18  35665723  40053938  NEU 0.138
14  19136137  22397429  DEL -0.740
14  22397430  22973824  DUP 0.949
```

E.g. usage
```
sh MRDetectCNA_sample_runner.sh
Usage: MRDetectCNA_sample_runner.sh BAM REF WIN_SIZE CNV_BED OUT_DIR SAMPLE_NAME INTERVAL_NAME

sh MRDetectCNA_sample_runner.sh sample1.bam ref.fa 500 cnv.bed <path to output directory> sample1 interval1 
```

### 2. Case vs Control processing 
This step takes pairwise case-control samples to generate detection metrics and Zscores. The control samples are typically a panel of normal (PON) or PBMC.

E.g. case-control pairs file
```
#case control case_average_coverage control_average_coverage  case_type
case1 PON case1_avg_cvg PON_avg_cvg pre
case2 PON case2_avg_cvg PON_avg_cvg post
control1  PON control1_avg_cvg PON_avg_cvg control
control2  PON control2_avg_cvg PON_avg_cvg control
```

E.g. usage
```
sh MRDetectCNA_pair_runner.sh
Usage: MRDetectCNA_pair_runner.sh INP_DIR OUT_DIR TN_PAIRS MEDIAN_THRESH PATIENT_NAME INTERVAL_NAME WIN_SIZE

sh MRDetectCNA_pair_runner.sh <path to input directory> <path to output directory> case_control.pairs.txt 1.5 patient1 interval1 500
```
