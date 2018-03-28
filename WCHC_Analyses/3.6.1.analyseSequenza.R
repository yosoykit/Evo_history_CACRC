#get sample list information
sampleList <- read.csv(file="~/Projects/IBDproject/masterSampleList.csv", header=TRUE, stringsAsFactors=FALSE)
setNames <- unique(sampleList[["setID"]])

#seqFiles <- "/Users/cross01/Projects/MseqVol2/3.6.sequenza/"
seqFiles <- argument[2]


for(currSam in 1:nrow(sampleList)){
  
  #current set
  currSet <- sampleList[currSam, "setID"]
  currID <- sampleList[currSam, "sampleID"]
  
  subSample <- sampleList[sampleList[[1]]==currSet, ]
  normalName <- subSample[1, "normalID"]
  
  if(currID == normalName){
    next()
  }
  
  print(paste("#### performing sequenza analysis for ", sampleList[currSam, "sampleID"], " ####",sep=""))


  #read data from tumour seqz file name
  seqFileName <- paste(seqFiles, currSet, "/", currID, ".seqz.gz", sep="")
  
  #analyze data using sequenza
  chromosomes <- paste("chr", c(1:22), sep="")
  seqzExt <- sequenza.extract(file=seqFileName, chromosome.list = chromosomes, gamma = 40, kmin = 50, min.reads.baf = 15)
  
  #infer cellularity and ploidy
  paraSpace <- sequenza.fit(seqzExt, cellularity = seq(0.5, 1, 0.1), ploidy = seq(2, 5, 0.1))
  
  #sequenza analysis
  sequenza.results(sequenza.extract = seqzExt, cp.table = paraSpace, sample.id = sampleList[currSam, "sampleID"], out.dir=paste(seqFiles, currSet, "/", sep="") )
  
}



