#$ -V
#$ -cwd

######
#Help Message
######

desc="This script runs a sample implementation of the MRDetectSNV algorithm. See readme for further details" 
usage="Usage: MRDetectSNV_runner.sh [Plasma BAM] [Tumor VCF] "
case ${1} in 
	-h) echo ~~~; echo $desc; echo $usage; echo ~~~; exit ;;
	--help) echo ~~~; echo $desc; echo $usage; echo ~~~; exit ;;

esac


######
#Set up variable file names
######

PLASMA_BAM=${1}
TUMOR_VCF=${2}

PLASMA_BAM_NAME=`basename $PLASMA_BAM`
TUMOR_VCF_NAME=`basename $TUMOR_VCF`
PLASMAvTUMOR_READS_FILE=${PLASMA_BAM_NAME}_VS_${TUMOR_VCF_NAME}.tsv
PLASMAvTUMOR_READS_FILE_SVMSCORED=${PLASMA_BAM_NAME}_VS_${TUMOR_VCF_NAME}.svm.tsv
PRETRAINED_SVM=./trained_SVM.pkl
OUTPUT=${PLASMA_BAM_NAME}_VS_${TUMOR_VCF_NAME}_RESULT.csv

echo ~~~~~~~~~~~~~~~
echo Running MRDetectSNV: Analysing  ${PLASMA_BAM_NAME} plasma using ${TUMOR_VCF_NAME} tumor signature

######
#STEP 1: Pull all reads from the Plasma BAM that overlap the sites in the Tumor VCF. Last argument is output file that will contain per-read data.
######
echo ~~~
echo Pulling candidate reads from plasma
python pull_reads.py --bam $PLASMA_BAM --vcf $TUMOR_VCF --out $PLASMAvTUMOR_READS_FILE

######
#STEP 2: Use SVM to score each read
######

echo ~~~
echo Quality scoring all reads
python quality_score.py --pickle-name ./trained_SVM.pkl --detections $PLASMAvTUMOR_READS_FILE --output_file $PLASMAvTUMOR_READS_FILE_SVMSCORED
######
#STEP 3: Filter Plasma Reads (Read-specific SVM score, Locus-specific Blacklist, R1R2). Then check for reads with mutations that match Tumor VCF. Outputs detection rate (# of tumor sites detected / total number of plasma reads at tumor sites) 
######
echo ~~~
echo Identifying tumor reads and calculating detection rate:
python filterAndDetect.py $TUMOR_VCF $PLASMAvTUMOR_READS_FILE_SVMSCORED $OUTPUT
echo Done. 
echo ~~~
echo Results written to ${OUTPUT}
echo ~~~~~~~~~~~~~~~

