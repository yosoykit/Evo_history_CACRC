#get sample list information
sampleList <- read.csv(file="~/Projects/IBDproject/masterSampleList.csv", header=TRUE, stringsAsFactors=FALSE)
setNames <- unique(sampleList[["setID"]])

platDir <- "~/Projects/IBDproject/annoPlatypusCalls/exomes/"

vcfDir <- "~/Projects/IBDproject/3.1.platypusCalls/exomes/"

bcfDir <-  "~/Projects/IBDproject/3.0.bcfCalls/"

scriptsOut <- "~/Projects/IBDproject/A.runScripts/"

system(command = paste("mkdir ", scriptsOut, "3.1.runPlatypus/", sep=""))

for(currSet in 1:length(setNames)){
  
  subSample <- sampleList[sampleList[["setID"]]==setNames[currSet], ]
  
  #output bam file list
  bamFileVect <- paste(subSample[1, "FileHolder2"], "2.processedBams/exomes/", setNames[currSet], "/", subSample[["sampleID"]], "/", subSample[["sampleID"]], ".mkdub.bam", sep="")
  bamListFileName <- paste(scriptsOut, "3.1.runPlatypus/", setNames[currSet], ".bamList.txt", sep="")
  write.table(matrix(bamFileVect, nrow = length(bamFileVect), ncol=1), file = bamListFileName, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  vcfFiles <- c()
  vcfCounter <- 1
  
  #for each chromosome make a platypus calling script
  for(currChr in c(1:22, "X", "Y")){
    platLogFile <- paste(vcfDir, setNames[currSet], "/", setNames[currSet], "_platypus.chr", currChr,".log", sep="")
    platVCFFile <- paste(vcfDir, setNames[currSet], "/", setNames[currSet], "_platypus.chr", currChr,".vcf", sep="")
    vcfFiles[vcfCounter] <- platVCFFile
    vcfCounter <- vcfCounter + 1
    bcftoolsFile <- paste(bcfDir, setNames[currSet], "/", setNames[currSet], ".bcftools.calls.vcf.gz", sep="")
    
    if(subSample[1, "regions"] == "WGS"){
      regionsString <- paste(" \\
--regions=chr", currChr, sep="")
    }else{
      regionsString <- paste(" \\
--regions=", subSample[1, "regions"], ".chr", currChr,".txt", sep="")
    }

    outName <- paste(scriptsOut, "3.1.runPlatypus/", "run_Platypus_", setNames[currSet], "_chr", currChr, ".sh", sep="") 
    totalStrings <- paste("#!/bin/sh
#!/bin/sh
#$ -cwd
#$ -V
#$ -pe smp 1
#$ -l h_rt=100:0:0
#$ -l h_vmem=4G

module load use.sl6
module load python
module load htslib
module load tabix

mkdir ", vcfDir, setNames[currSet],"/

python2.7 /data/BCI-EvoCa2/jacob/Central/software/Platypus_0.8.1/Platypus.py callVariants \\
--bamFiles=", vcfDir, setNames[currSet],".bamList.txt \\
--output=", platVCFFile," \\
--refFile=", subSample[1, "reference"]," \\
--nCPU=1 \\
--minReads=2 \\
--maxVariants=100 \\
--mergeClusteredVariants=1 \\
--minMapQual=1 \\
--genIndels=0 \\
--minPosterior=0 \\
--logFileName=", platLogFile," \\
--source=", bcftoolsFile, regionsString,"

module unload python
module unload htslib
module unload tabix
module unload use.sl6", sep="")
    lapply(totalStrings, write, outName, append=FALSE)
  }
    
    
  #make vcf merging script
  mergedVCF <- paste(vcfDir, setNames[currSet], "/", setNames[currSet], ".merged.vcf", sep="")
  
  outNameMerge <- paste(scriptsOut, "merge_finalVCF_", setNames[currSet], ".sh", sep="") 
  totalStringsMerge <- paste("#!/bin/sh
#$ -cwd
#$ -V
#$ -pe smp 1
#$ -l h_rt=10:0:0      # Request 100 hour runtime
#$ -l h_vmem=2G

module load vcftools

vcf-concat ", paste(vcfFiles, collapse = " "), " > ", mergedVCF, sep="")
  
  #write seq script file
  lapply(totalStringsMerge, write, outNameMerge, append=FALSE)
    
}

