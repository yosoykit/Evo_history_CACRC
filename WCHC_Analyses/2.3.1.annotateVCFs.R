###### subroutine #####

#function to annotate with annovar
addAnnotations <- function(annoSub, exomeSub, mutTypeSub){
  #add mutation annotations for filtered variants
  annoSub <- data.frame(append(x = annoSub, values = NA, after = 2), stringsAsFactors = FALSE)
  names(annoSub)[3] <- "mutation"
  annoSub <- data.frame(append(x = annoSub, values = NA, after = 3), stringsAsFactors = FALSE)
  names(annoSub)[4] <- "mutationtype"
  
  tracker <- round(seq(1, nrow(annoSub), length.out = 10), digits = 0)
  tracker <- tracker[-1]
  
  for(currRow in 1:nrow(annoSub)){
    if(currRow %in% tracker){
      print(paste("##### annotating ", (which(currRow==tracker))*10, "% complete ####", sep=""))
    }
    if(annoSub[currRow, "start"] %in% exomeSub[[5]]){
      if(nrow(exomeSub[exomeSub[[4]]==annoSub[currRow, "chr"] & exomeSub[[5]]==annoSub[currRow, "start"], ]) != 0){
        annoSub[currRow, "mutation"] <- exomeSub[exomeSub[[4]]==annoSub[currRow, "chr"] & exomeSub[[5]]==annoSub[currRow, "start"], 3][1]
        annoSub[currRow, "mutationtype"] <- exomeSub[exomeSub[[4]]==annoSub[currRow, "chr"] & exomeSub[[5]]==annoSub[currRow, "start"], 2][1]
      }
    }else{
      tempGeneId <- strsplit(annoSub[currRow, "gene"], split = "[(]")
      annoSub[currRow, "gene"] <- tempGeneId[[1]][1]
      annoSub[currRow, "mutation"] <- tempGeneId[[1]][2]
      annoSub[currRow, "mutationtype"] <- mutTypeSub
    }
  }
  return(annoSub)
}


#function to add protein annotations
annoBiomart <- function(annoSub, bioMartSub){
  annoSub <- data.frame(append(x = annoSub, values = NA, after = 4), stringsAsFactors = FALSE)
  names(annoSub)[5] <- "proteinID"
  annoSub <- data.frame(append(x = annoSub, values = NA, after = 5), stringsAsFactors = FALSE)
  names(annoSub)[6] <- "band"
  
  tracker <- round(seq(1, nrow(annoSub), length.out = 10), digits = 0)
  tracker <- tracker[-1]
  
  for(currRow in 1:nrow(annoSub)){
    if(currRow %in% tracker){
      print(paste("##### annotating ", (which(currRow==tracker))*10, "% complete ####", sep=""))
    }
    
    if(annoSub[currRow, "gene"] %in% bioMartSub[["geneID"]]){
      bioTabSub <- bioMartSub[annoSub[currRow, "gene"] == bioMartSub[["geneID"]], ]
      protID <- bioTabSub[["uniProtSwissID"]]
      protID <- protID[protID!=""]
      if(length(protID)==0){
        annoSub[currRow, "proteinID"] <- "noID"
      }else if(length(unique(protID))>1){
        annoSub[currRow, "proteinID"] <- paste(unique(protID), collapse=":")
      }else{
        annoSub[currRow, "proteinID"] <- protID[1]
      }
      annoSub[currRow, "band"] <- bioTabSub[1,"band"]
    }else{
      annoSub[currRow, "proteinID"] <- "noID"
      annoSub[currRow, "band"] <- "noID"
    }
  }
  return(annoSub)
}



###### begin ########
#get sample list information
sampleList <- read.csv(file="~/Projects/IBDproject/masterSampleList.csv", header=TRUE, stringsAsFactors=FALSE)
setNames <- unique(sampleList[["setID"]])

platDir <- "~/Projects/IBDproject/annoPlatypusCalls/exomes/"

#mutType <- "indel"
mutType <- "snv"

#summary table
summTab <- data.frame(matrix(data = NA, nrow = length(setNames), ncol = 4))
names(summTab) <- c("set", "noDrivers", "noSom", "noGerm")
summTab["set"] <- setNames

fileName <- paste(".", mutType, ".somatic.annoVar.variant_function.txt", sep="")
fileNameExome <- paste(".", mutType, ".somatic.annoVar.exonic_variant_function.txt", sep="")

fileNameGermline <- paste(".germline.annoVar.variant_function.txt", sep="")
fileNameGermlineExome <- paste(".germline.annoVar.exonic_variant_function.txt", sep="")


#get biomart table
bioMart <- read.table(file="~/Projects/ReferenceGenome/biomartExome.hg19.IDs.bed", header = TRUE, sep="\t", stringsAsFactors = FALSE)
bioMart <- bioMart[bioMart[["chrom"]] %in% c(1:22, "X", "Y"),]

#perform pipeline for each sample
for(j in 1:length(setNames)){
  currSetId <- setNames[j]
	print(paste("#### analysing variants for sample ", currSetId, " ####",sep=""))
	
	#setup input/output names
  subSample <- sampleList[sampleList[[1]]==currSetId, ]
  
  normName <- subSample[1, "normalID"]
  
  #get biopsy names from vcf file
  biopList <-  read.table(file=paste(platDir, currSetId, "/", currSetId, ".vcfbioID.txt", sep=""), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  biopList <- biopList[[1]]
  
  # somatic call file
  annoSom <- read.table(file=paste(platDir, currSetId, "/", currSetId, fileName, sep=""), sep="\t", header = FALSE, stringsAsFactors = FALSE)
  names(annoSom) <- c("type", "gene", "chr", "start", "end", "ref", "alt", paste(biopList, ".NR", sep=""), paste(biopList, ".NV", sep=""))
  
  exomeSom <- read.table(file=paste(platDir, currSetId, "/", currSetId, fileNameExome, sep=""), sep="\t", header = FALSE, stringsAsFactors = FALSE)
  
 
  # germline variants file
  totalGermline <- read.table(file=paste(platDir, currSetId, "/", currSetId, fileNameGermline, sep=""), sep="\t", header = FALSE, stringsAsFactors = FALSE)
  names(totalGermline) <- c("type", "gene", "chr", "start", "end", "ref", "alt", paste(biopList, ".NR", sep=""), paste(biopList, ".NV", sep=""))
  
  totalGermlineExome <- read.table(file=paste(platDir, currSetId, "/", currSetId, fileNameGermlineExome, sep=""), sep="\t", header = FALSE, stringsAsFactors = FALSE)
  
  
  ######## 1. rearranges vcf tables to include mutations ######## 
  
  for(currCol in 1:length(biopList)){
    annoSom[ncol(annoSom)+1] <- as.numeric(annoSom[[paste(biopList[currCol], ".NV", sep="")]]) / as.numeric(annoSom[[paste(biopList[currCol], ".NR", sep="")]])
    names(annoSom)[ncol(annoSom)] <- paste(biopList[currCol], ".VAF", sep="")
    
    totalGermline[ncol(totalGermline)+1] <- as.numeric(totalGermline[[paste(biopList[currCol], ".NV", sep="")]]) / as.numeric(totalGermline[[paste(biopList[currCol], ".NR", sep="")]])
    names(totalGermline)[ncol(totalGermline)] <- paste(biopList[currCol], ".VAF", sep="")
  }
  
  #add mutation annotations for exonic somatic calls
  print("#### annotating exome somatic calls ####")
  annoSom <- addAnnotations(annoSom, exomeSom, mutType)
  
  #print("#### annotating germline somatic calls ####")
  #totalGermline <- addAnnotations(totalGermline, totalGermlineExome, mutType)
  
  
  ######### 2. gets protein and transcript IDs from biomart table ############
  
  print("#### getting biomart IDs for somatic calls ####")
  annoSom <- annoBiomart(annoSom, bioMart)
  
  #print("#### getting biomart IDs for germline calls ####")
  #totalGermline <- annoBiomart(totalGermline, bioMart)
  
  
  ######## 3.plot VAF distributions ##########
  
  pdf(file = paste(platDir, currSetId, "/", currSetId, ".", mutType, ".VAFplots.pdf", sep=""), width = 10, height = (nrow(subSample)*5))
    par(mfrow=c(nrow(subSample), 2), mar=c(5,3,3,3), cex=1.2)
    
    for(currSam in 1:length(biopList)){
      #plot somatic muts
      tempVAR <- annoSom[[paste(biopList[currSam], ".VAF", sep="")]]
      tempVAR <- tempVAR[tempVAR > 0]
      hist(tempVAR, breaks = seq(0,1,0.01), xlab = "", main = paste(biopList[currSam], "somatic mut distribution"), col="darkgreen", border = "white")
      meanVAR <- round(mean(tempVAR), digits = 2)
      noVAR <- length(tempVAR)
      title(sub = paste("mean VAF =", meanVAR, "\n", "no muts =", noVAR))
      
      #plot GL muts
      tempVAR <- totalGermline[[paste(biopList[currSam], ".VAF", sep="")]]
      tempVAR <- tempVAR[!is.na(tempVAR)]
      tempVAR <- tempVAR[tempVAR > 0]
      hist(tempVAR, breaks = seq(0,1,0.01), xlab = "", main = paste(biopList[currSam], "GL distribution"), col="darkblue", border = "white")
      meanVAF <- round(mean(tempVAR), digits = 2)
      noVAR <- length(tempVAR)
      title(sub = paste("mean VAF =", meanVAF, "\n", "no muts =", noVAR))
    }
  dev.off()
  
  
  #output somatic variant table
  somOut <- paste(platDir, currSetId, "/", currSetId, ".", mutType, ".somatic.annoVar.comp.txt", sep="")
  write.table(x=annoSom, file = somOut, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  #output somatic variant table
  #germOut <- paste(platDir, currSetId, "/", currSetId, ".", mutType, ".germline.annoVar.comp.txt", sep="")
  #write.table(x=totalGermline, file = germOut, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  ###### count variants and add to table #######
  summTab[j, "noDrivers"] <- NA
  summTab[j, "noSom"] <- nrow(annoSom)
  summTab[j, "noGerm"] <- nrow(totalGermline)
}

#write summary table
summaryOut <- paste(platDir, "Platypus.", mutType, ".summary.txt", sep="")
write.table(summTab, file = summaryOut, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

