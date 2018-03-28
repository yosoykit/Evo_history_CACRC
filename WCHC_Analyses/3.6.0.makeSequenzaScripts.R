sampleList <- read.csv(file="~/Projects/IBDproject/masterSampleList.csv", header=TRUE, stringsAsFactors=FALSE)
setNames <- unique(sampleList[["setID"]])

bamFiles <- "~/Projects/IBDproject/2.bamFiles/"

seqFiles <- "~/Projects/IBDproject/3.6.sequenza/"

scriptDir <- "~/Projects/IBDproject/A.runScripts/3.6.sequenza/"

GCcontentFile <- "~/Projects/IBDproject/hg19.gc50Base.txt.gz"
refFile <- "~/Projects/IBDproject/ucsc.hg19.fasta"

system(paste("mkdir ", scriptDir, sep=""))

#make sequenza prep scripts (for each sample)
for(currSam in 1:nrow(sampleList)){
  print(paste("#### making sequenza script for sample ", sampleList[currSam, "sampleID"], " ####",sep=""))
  
  sampleID <- sampleList[currSam, "sampleID"]
  
  #current sample
  setID <- sampleList[currSam, "setID"]
  
  #current normal
  subSample <- sampleList[sampleList[["setID"]]==setID, ]
  
  #normal name
  normalName <- subSample[1,"normalID"]
  
  if(sampleID == normalName){
    next()
  }
  
  #bam file names
  tumourBam <- paste(bamFiles, setID, "/", sampleID, "/", sampleID, ".mkdub.bam", sep="")
  normBam <- paste(bamFiles, setID, "/", normalName, "/", normalName, ".mkdub.bam", sep="")

  #sequenza files
  seqFile <- paste(seqFiles, setID, "/", sampleID, ".seqz.gz", sep="")
  binnedSeqFile <- paste(seqFiles, setID, "/", sampleID, ".seqz.binned.gz", sep="")
 
  #prep .sh file name
  SHfile <- paste(scriptDir, "sequenzaPrep.", setID, ".", sampleID, ".sh", sep="")
  
  prepSHstrings <- paste("#!/bin/sh
#$ -cwd
#$ -V
#$ -pe smp 1        # Request 1 CPU core1
#$ -l h_vmem=20G    # Request 20GB RAM / core, i.e. 24GB total
#$ -l h_rt=48:0:0   # Request 48 hour runtime

#make analysis dir
mkdir ", seqFiles, setID,"

#get .seqz file
~/bin/sequenza-utils.py bam2seqz -gc ", GCcontentFile," --fasta ", refFile," -n ", normBam," -t ", tumourBam," | gzip > ", seqFile,"

#bin .seqz file to shorten analysis time
~/bin/sequenza-utils.py seqz-binning -w 250 -s ", seqFile," | gzip > ", binnedSeqFile,"


#run pre-processed files using sequenza in R locally
", sep="")

  #write prep .sh file
  lapply(prepSHstrings, write, SHfile, append=FALSE)
     
}

