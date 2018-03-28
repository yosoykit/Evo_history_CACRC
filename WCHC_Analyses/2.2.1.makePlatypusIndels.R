#get sample list information
sampleList <- read.csv(file="~/Projects/IBDproject/masterSampleList.csv", header=TRUE, stringsAsFactors=FALSE)
setNames <- unique(sampleList[["setID"]])

vcfDir <- "~/Projects/IBDproject/4.platypusIndels/"

bcfDir <-  "~/Projects/IBDproject/3.bcfIndels/"

scriptsOut <- "~/Projects/IBDproject/"

system(command = paste("mkdir ", scriptsOut, "runPlatypusIndel/", sep=""))

for(currSet in 1:length(setNames)){
  currSetId <- setNames[currSet]
  print(paste("#### making scripts for sample ", currSetId, " ####",sep=""))
  
  subSample <- sampleList[sampleList[["setID"]]==setNames[currSet], ]
  vcfFiles <- c()
  vcfCounter <- 1
  
  #for each chromosome make a platypus calling script
  for(currChr in c(1:22, "X", "Y")){
    platLogFile <- paste(vcfDir, setNames[currSet], "/", setNames[currSet], "_platypus.chr", currChr,".indel.log", sep="")
    platVCFFile <- paste(vcfDir, setNames[currSet], "/", setNames[currSet], "_platypus.chr", currChr,".indel.vcf", sep="")
    vcfFiles[vcfCounter] <- platVCFFile
    vcfCounter <- vcfCounter + 1
    bcftoolsFile <- paste(bcfDir, setNames[currSet], "/", setNames[currSet], ".bcftools.indels.vcf.gz", sep="")
    
    if(subSample[1, "regions"] == "WGS"){
      regionsString <- paste("--regions=chr", currChr, sep="")
    }else{
      regionsString <- paste("--regions=", subSample[1, "regions"], "_chr", currChr,".txt, sep=")
    }
    
    outName <- paste(scriptsOut, "run_Platypus_indels_", setNames[currSet], "_chr", currChr, ".sh", sep="") 
    totalStrings <- paste("#!/bin/sh
#$ -cwd
#$ -V
#$ -pe smp 6
#$ -l h_vmem=10G
#$ -l h_rt=100:0:0
   
source activate ccevo
 
mkdir ", vcfDir, setNames[currSet],"/

platypus callVariants \\
--bamFiles=", vcfDir, setNames[currSet],".bamList.txt \\
--output=", platVCFFile," \\
--refFile=", subSample[1, "reference"]," \\
--nCPU=6 \\
--minReads=25 \\
--maxVariants=100 \\
--mergeClusteredVariants=1 \\
--minMapQual=1 \\
--minPosterior=0 \\
--minVarFreq=0.1 \\
--maxSize=100 \\
--genSNPs=0 \\
--logFileName=", platLogFile," \\
--source=", bcftoolsFile," \\
",regionsString,"
", sep="")
    lapply(totalStrings, write, outName, append=FALSE)
  }
    
    
  #make vcf merging script
  mergedVCF <- paste(vcfDir, setNames[currSet], "/", setNames[currSet], ".merged_indels.vcf", sep="")
  
  outNameMerge <- paste(scriptsOut, "merge_indelVCF_", setNames[currSet], ".sh", sep="") 
  totalStringsMerge <- paste("#!/bin/sh
#$ -cwd
#$ -V
#$ -pe smp 1
#$ -l h_rt=100:0:0      # Request 100 hour runtime
#$ -l h_vmem=2G

source activate ccevo

vcf-concat ", paste(vcfFiles, collapse = " "), " > ", mergedVCF, sep="")
  
  #write seq script file
  lapply(totalStringsMerge, write, outNameMerge, append=FALSE)
    
}

