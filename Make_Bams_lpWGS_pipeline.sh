#!/bin/sh
#$ -cwd
#$ -V
#$ -e sh_output
#$ -o sh_output
#$ -l h_rt=10:0:0 # Request 10 hour runtime
#$ -t 1-81
#$ -pe smp 6   
#$ -l h_vmem=10G # Request 10GB RAM / core

longLine="--------------------"

##############################################################################
# Go from fastq files to BAM files.
# QC fatsqs using fatsqc, Map files using BWA, mark duplicates using picard
# calculate coverage statistica dn check BAM using bamQC
##############################################################################

#files list should have 6 columns, read1name, read2name, patient, samplename, tissuetype (eg. HGD, LGD, C), and lesion # (to identify samples from same lesion)
filelist= LGD_HGD_81_samples.txt

##############################################################################
#set up directories and filenames
msg="Set up directories and filenames"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
fastqdir=/data/Raw
results_dir=/data/Results

mkdir -p ${results_dir}/processingbams
mkdir -p ${results_dir}/finalbams
mkdir -p ${results_dir}/coverage/plots
mkdir -p ${results_dir}/fastQC
mkdir -p ${results_dir}/bamQC
mkdir -p ${results_dir}/info
mkdir -p ${results_dir}/first_trimming
mkdir -p ${results_dir}/second_trimming
mkdir -p ${results_dir}/QDNAseq

read1=$(awk -v var="$SGE_TASK_ID" 'NR ==var { OFS="\t";print $1}' $filelist)
read2=$(awk -v var="$SGE_TASK_ID" 'NR ==var { OFS="\t";print $2}' $filelist)
patient=$(awk -v var="$SGE_TASK_ID" 'NR ==var { OFS="\t";print $3}' $filelist)
samplename=$(awk -v var="$SGE_TASK_ID" 'NR ==var { OFS="\t";print $4}' $filelist)
tissuetype=$(awk -v var="$SGE_TASK_ID" 'NR ==var { OFS="\t";print $5}' $filelist)

readgroup="@RG\tID:${patient}.${samplename}\tLB:${patient}.${samplename}\tSM:${patient}.${samplename}\tPL:ILLUMINA"

msg="Sample: ${patient}.${samplename}"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"

##############################################################################
msg="trim adapter ends with skewer"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
module load use.dev skewer
skewer -x AGATCGGAAGAGC -m any -l 35 -r 0.1 -d 0.03 -Q 10 -n ${fastqdir}/${read1} ${fastqdir}/${read2} -o ${results_dir}/first_trimming/${patient}.${samplename}
skewer -x ACGCTCTTCCGATCT -m head -l 35 -r 0.1 -d 0.03 -Q 10 -n ${results_dir}/first_trimming/${patient}.${samplename}-trimmed-pair1.fastq ${results_dir}/first_trimming/${patient}.${samplename}-trimmed-pair2.fastq -o ${results_dir}/second_trimming/${patient}.${samplename}

##############################################################################
msg="run fastqc"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
#perform fastqc quality control
module load fastqc/0.11.5
fastqc ${results_dir}/second_trimming/${patient}.${samplename}-trimmed-pair1.fastq ${results_dir}/second_trimming/${patient}.${samplename}-trimmed-pair2.fastq --outdir=${results_dir}/fastQC/

##############################################################################

#add read group headers and align using bwa mem, pipe to samtools to convert to bam
msg="map with BWA and pipe to samtools to convert to bam"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
module load bwa/0.7.8
bwa mem -M -t 6 -R ${readgroup} /data/refs/hg19/ucsc.hg19.fasta ${results_dir}/second_trimming/${patient}.${samplename}-trimmed-pair1.fastq ${results_dir}/second_trimming/${patient}.${samplename}-trimmed-pair2.fastq  > ${results_dir}/processingbams/${patient}.${samplename}.sam	


module load samtools/1.3.1
samtools view -S -b -q 37 ${results_dir}/processingbams/${patient}.${samplename}.sam > ${results_dir}/processingbams/${patient}.${samplename}_unsort.bam

##############################################################################

#sort bam file an index using picard
msg="sort bam with picard"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
module load java/oracle/1.8.0_11 
java -Xmx8G -jar /data/home/picard.jar SortSam INPUT=${results_dir}/processingbams/${patient}.${samplename}_unsort.bam \
OUTPUT=${results_dir}/processingbams/${patient}.${samplename}.bam \
SORT_ORDER=coordinate CREATE_INDEX=true

##############################################################################

#remove unsorted bam file
msg="remove unsorted bam"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
rm ${results_dir}/processingbams/${patient}.${samplename}_unsort.bam

##############################################################################

#mark duplicates
msg="mark duplicates with picard"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"

java -Xmx8G -jar /data/home/picard.jar MarkDuplicates \
INPUT=${results_dir}/processingbams/${patient}.${samplename}.bam \
OUTPUT=${results_dir}/processingbams/${patient}.${samplename}.dedup.bam \
METRICS_FILE=${results_dir}/processingbams/${patient}.${samplename}.metrics.file CREATE_INDEX=true

##############################################################################

#remove duplicates using samtools
msg="remove duplicates"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
samtools rmdup ${results_dir}/processingbams/${patient}.${samplename}.dedup.bam ${results_dir}/processingbams/${patient}.${samplename}.dedup2.bam

##############################################################################

#index using picard
msg="index with picard"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"

java -Xmx8G -jar /data/home/picard.jar BuildBamIndex INPUT=${results_dir}/processingbams/${patient}.${samplename}.dedup2.bam \
OUTPUT=${results_dir}/processingbams/${patient}.${samplename}.dedup2.bai

##############################################################################

#copy to final bams directory
cp ${results_dir}/processingbams/${patient}.${samplename}.dedup2.bam ${results_dir}/finalbams/${patient}.${samplename}.bam
cp ${results_dir}/processingbams/${patient}.${samplename}.dedup2.bai ${results_dir}/finalbams/${patient}.${samplename}.bam.bai

##############################################################################

#check with bamqc
msg="run bamqc"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
/data/home/bamqc -o ${results_dir}/bamQC/ ${results_dir}/finalbams/${patient}.${samplename}.bam

##############################################################################

#calculate coverage statistics
msg="run GATK coverage"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
module load java/oracle/1.8.0_11
java -jar -Xmx4G /data/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -T DepthOfCoverage \
-R /data/refs/hg19/ucsc.hg19.fasta \
-I ${results_dir}/finalbams/${patient}.${samplename}.bam \
--start 1 \
--stop 1000 \
--nBins 500 \
--omitDepthOutputAtEachBase \
-o ${results_dir}/coverage/${patient}.${samplename}.coverage

msg="finished"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"