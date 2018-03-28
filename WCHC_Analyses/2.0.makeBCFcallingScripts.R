#get sample list information
sampleList <- read.csv(file="~/Projects/IBDproject/masterSampleList.csv", header=TRUE, stringsAsFactors=FALSE)
setNames <- unique(sampleList[["setID"]])

bamDir <- "~/Projects/IBDproject/2.processedBams/"

vcfOut <- "~/Projects/IBDproject/3.0.bcfCalls/"

scriptsOut <- "~/Projects/IBDproject/A.runScripts/"

system(command = paste("mkdir ", scriptsOut, "3.0.bcfCalls/", sep=""))

for(currSet in 1:length(setNames)){
  subSample <- sampleList[sampleList[["setID"]]==setNames[currSet], ]
  normalSample <- subSample[1,"normalID"]
  biopsyNames <- subSample[["sampleID"]]

  #seq file list
  vcfFiles <- c()
  vcfCounter <- 1
  
  #make mutect scripts for each biopsy
  for(currBio in 1:length(biopsyNames)){
    print(paste("####### making bcftools scripts for sample", subSample[currBio, "sampleID"], "#########"))
    bamFileIn <- paste(bamDir, setNames[currSet], "/", biopsyNames[currBio], "/", biopsyNames[currBio], ".mkdub.bam", sep="")
    
    #save vcfs for by chromosome concatenation
    vcfsByChromsome <- c()
    chrCounter <- 1
    
    for(currChr in c(1:22, "X", "Y")){
      outputvcf <- paste(vcfOut, setNames[currSet], "/", setNames[currSet], "_", biopsyNames[currBio], "_chr", currChr, "_bcftools.vcf", sep="")
      
      vcfsByChromsome[chrCounter] <- paste(outputvcf, ".gz", sep="")
      chrCounter <- chrCounter + 1
      
      #set regions file, if applicable
      if(sampleList[currSet, "regions"] == "WGS"){
        samtoolsCommand <- paste("samtools mpileup -d 250 -r chr", currChr," -u -f ", subSample[currBio, "reference"], " ", bamFileIn," | ", sep="")
      }else{
        samtoolsCommand <- paste("samtools mpileup -d 250 -u -r ", currChr," -l ", subSample[currBio, "regions"]," -f ", subSample[currBio, "reference"], " ", bamFileIn," | ", sep="")
      }
      
      outName <- paste(scriptsOut, "3.bcfCalls/run_bcftools_", setNames[currSet], "_", biopsyNames[currBio], "_chr", currChr, ".sh", sep="") 
      totalStrings <- paste("#!/bin/sh
#$ -cwd
#$ -V
#$ -pe smp 1
#$ -l h_vmem=6G
#$ -l h_rt=48:0:0

source activate ccevo
#module load bcftools/1.4
#module load samtools

mkdir ", vcfOut, setNames[currSet], "/

#run bcftools 
",samtoolsCommand," \\
bcftools call --skip-variants indels -cv -P 1.1e-1 -O v > ", outputvcf, "

bgzip ", outputvcf," 
tabix -p vcf ", outputvcf,".gz
", sep="")
      
      #write seq script file
      lapply(totalStrings, write, outName, append=FALSE)
    }
    
    
    #make vcf merging script
    outNameMerge <- paste(scriptsOut, "3.bcfCalls/merge_vcftools_byChr_", biopsyNames[currBio], ".sh", sep="") 
    mergedVCF <- paste(vcfOut, setNames[currSet], "/", biopsyNames[currBio], ".bcftools.calls.vcf", sep="")
    
    vcfFiles[vcfCounter] <- paste(mergedVCF, ".gz", sep="")
    vcfCounter <- vcfCounter + 1
    
    totalStringsMerge <- paste("#!/bin/sh
#$ -cwd
#$ -V
#$ -pe smp 1
#$ -l h_vmem=5G
#$ -l h_rt=1:0:0      # Request 100 hour runtime

source activate ccevo

vcf-concat ", paste(vcfsByChromsome, collapse = " ")," | bgzip > ", mergedVCF, ".gz

tabix -p vcf ", mergedVCF,".gz", sep="")

    #write seq script file
    lapply(totalStringsMerge, write, outNameMerge, append=FALSE)
    
  }
  
    
    
    
  #make vcf merging script
  outNameMerge <- paste(scriptsOut, "3.bcfCalls/merge_vcftools_", setNames[currSet], ".sh", sep="") 
  mergedVCF <- paste(vcfOut, setNames[currSet], "/", setNames[currSet], ".bcftools.calls.vcf", sep="")
  totalStringsMerge <- paste("#!/bin/sh
#$ -cwd
#$ -V
#$ -pe smp 1
#$ -l h_vmem=5G
#$ -l h_rt=1:0:0      # Request 100 hour runtime

source activate ccevo

vcf-merge ", paste(vcfFiles, collapse = " ")," | bgzip > ", mergedVCF, ".gz

tabix -p vcf ", mergedVCF,".gz", sep="")
  
  #write seq script file
  lapply(totalStringsMerge, write, outNameMerge, append=FALSE)
    
}






