## MRDetectSNV

This implementation of MRDetectSNV has three modular steps: 
1) extracting individual reads from an input plasma, selecting only from sites present in an input set of tumor mutations
2) giving a quality score to to each individual read
3) filtering out reads (based on per-read quality score or per-site genomic blacklist), and then calculating the detection rate.

This should be run on individual pairings of plasma samples (BAM) and mutation lists (VCF). If the user has a cohort of multiple VCFs and BAMs (which could include both case and control), this protocol should be run on all BAM-VCF pairings, and then detection metrics can be evaluated across samples.


### Python Dependencies 
1. Cython==0.29.15
2. numpy==1.16.4
3. pandas==0.24.2
4. pysam==0.15.4
5. tqdm==4.43.0
6. sklearn
7. pickle

### MRDetectSNV_runner.sh

```
Usage: MRDetectCNA_sample_runner.sh [BAM] [VCF]

bash MRDetectSNV_runner.sh sample.bam sample.vcf 
```

This script runs the 3 modular steps in sequence.


### 1. pull_reads.py
pulls reads from the plasma bam that overlap any of the sites of interest in the VCF. 
- The VCF must include only SNVs (single nucleotide variations)
- Make sure that the BAM and VCF are both aligned to the same reference genome
- Make sure the BAM and VCF have consistent chromosome naming practices (i.e. "chr12" or "12")

### 2. quality_score.py
Assesses the quality of each read, and creates a new file with same output as pull_reads.py, with an additional quality score column appended. The implementation presented here scores each read using a support vector machine (SVM) that integrates base quality, mean read base quality, mapping quality, position in read, and orientation (R1 or R2) of the variant of interest. The pretrained SVM is provided as trained_SVM.pkl. 

### 3. filterAndDetect.py
Filters reads, then reports final detection metrics. The filter includes a read-specific filter, which uses the scores appended in step2, as well as a locus-specific blacklist filter. A sample blacklist is provided here as blacklist.txt.gz.  The final detection metrics reported include the detection rate. This value indicates the frequency at which tumor specific mutations are seen in the plasma DNA. 


