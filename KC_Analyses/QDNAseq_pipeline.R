source("https://bioconductor.org/biocLite.R")
biocLite("QDNAseq.hg19")
library(QDNAseq)
#directory of bams

setwd("/data/Results/finalbams")

#directory to store plots

plotdir<-"/data/Results/QDNAseq/plots500kbp"

## output of list with all copy number calling information
outputdir<-"/data/Results/QDNAseq/CNAcalling"

run_size = 81

#download bin annotations
bins500 <- getBinAnnotations(binSize=500)


#process each bam individually, first read in bam filenames (n=81)
samples<-read.table('../../LGD_HGD_81_samples_final_withlesion.txt', colClasses = c("character"))[,3:4]
filenames<-paste(samples[,1],samples[,2],sep=".")[1:run_size]
callnames=c("copynumber", "probgain"   ,"probdloss",  "probamp" ,   "segmented" , "probnorm" ,  "calls"   ,   "probloss") 

## Make list for all call values for each patient
call_stats <- list(matrix(rep(0,8*4401),nrow= 8))
names(call_stats)[1]=filenames[1]
rownames(call_stats[[filenames[1]]])= callnames 

for (i in 2:run_size){
  call_stats[[i]]=matrix(rep(0,8*4401),nrow= 8)
  names(call_stats)[i]=filenames[i]
  rownames(call_stats[[filenames[i]]])= callnames 
}

for (i in filenames[1:run_size]){

  #load sequencing data
  readCounts <- binReadCounts(bins500,bamfiles=paste(i,".bam",sep=""), pairedEnds = TRUE)
  
  #plot raw readcounts
  pdf(paste(plotdir,"/","raw_profile",i,".pdf",sep=""),6,4)
  plot(readCounts, logTransform=FALSE, ylim=c(-10, 2000),main=paste("Raw Read Counts ",i,sep=""))
  highlightFilters(readCounts, logTransform=FALSE,residual=TRUE, blacklist=TRUE)
  dev.off()
  
  #apply filters and plot median read counts per bin as a function of GC content and mappability
  readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)

  ## Uncomment next two lines when using 500 kbp bins, spurious bins found also in normal blood
  bad_bins = c(1750, 1751, 1752, 1753 ,1754, 1755)
  readCountsFiltered@featureData$use[readCountsFiltered@featureData$use==T][bad_bins]=FALSE
  pdf(paste(plotdir,"/","isobar_",i,".pdf",sep=""),5,5)
  isobarPlot(readCountsFiltered)
  dev.off()
  
  #Estimate the correction for GC content and mappability, and make a plot for the relationship between the 
  #observed standard deviation in the data and its read depth
  readCountsFiltered <- estimateCorrection(readCountsFiltered)
  
  pdf(paste(plotdir,"/","noise_",i,".pdf",sep=""),5,5)
  noisePlot(readCountsFiltered)
  dev.off()
  
  #apply the correction for GC content and mappability which we then normalize, smooth outliers, calculate segmentation 
  #and plot the copy number profile
  copyNumbers <- correctBins(readCountsFiltered)
  copyNumbersNormalized <- normalizeBins(copyNumbers)
  copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
  
  pdf(paste(plotdir,"/","copynumber_profile_",i,".pdf",sep=""),7,4)
  plot(copyNumbersSmooth, ylim=c(-2,2))
  dev.off()
  
  # Segment copy number profile using DNAcopy 
  copyNumbersSegmented <- segmentBins(copyNumbersSmooth, segmentStatistic="seg.median")
  copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
  
  
  pdf(paste(plotdir,"/","segments",i,".pdf",sep=""),7,4)
  plot(copyNumbersSegmented, ylim=c(-2,2))
  dev.off()

  
  # Call copy number aberrations using CGHcall
  copyNumbersCalled <- callBins(copyNumbersSegmented)
  
  pdf(paste(plotdir,"/","copy_number_calls_",i,".pdf",sep=""),7,4)
  plot(copyNumbersCalled, ylim=c(-2,2))
  dev.off()

  ind = which(!is.na(copyNumbersCalled@assayData$copynumber))

  call_stats[[i]]['copynumber',]= copyNumbersCalled@assayData$copynumber[ind]
  call_stats[[i]]['probgain',]= copyNumbersCalled@assayData$probgain[ind]
  call_stats[[i]]['probdloss',]= copyNumbersCalled@assayData$probdloss[ind]
  call_stats[[i]]['probamp',]= copyNumbersCalled@assayData$probamp[ind]
  call_stats[[i]]['segmented',]= copyNumbersCalled@assayData$segmented[ind]
  call_stats[[i]]['probnorm',]= copyNumbersCalled@assayData$probnorm[ind]
  call_stats[[i]]['calls',]= copyNumbersCalled@assayData$calls[ind]
  call_stats[[i]]['probloss',]= copyNumbersCalled@assayData$probloss[ind]
}

call_stats_LGDHGD_500kbp = call_stats 

save(call_stats_LGDHGD_500kbp,file = paste(outputdir,"/","CNAcall_list_LGDHGD_500kbp.RData",sep=""))

