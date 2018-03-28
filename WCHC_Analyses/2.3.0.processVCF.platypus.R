#get sample list information
sampleList <- read.csv(file="~/Projects/IBDproject/masterSampleList.csv", header=TRUE, stringsAsFactors=FALSE)
setNames <- unique(sampleList[["setID"]])

platDir <- "~/Projects/IBDproject/annoPlatypusCalls/exomes/"

vcfName <- ".merged_indels.vcf"

mutType <- "snv"
#mutType <- "indel"

outAnno <- "~/Projects/IBDproject/6.indelPlatypusCalls/annotatedIndels/"

#concatenate names in table and delete unwanted columns
sampleList["sampleInfo"] <- paste(platDir, sampleList[[1]],"/", sampleList[[1]], vcfName, sep="")


#perform pipeline for each sample
for(j in 1:length(setNames)){
  currSetId <- setNames[j]
	print(paste("#### filtering sample ", currSetId, " ####",sep=""))
	
	#setup input/output names
  subSample <- sampleList[sampleList[[1]]==currSetId, ]
  
	#### these file are temporary and later deleted ####
	FILTName <- paste(platDir, currSetId,"/", currSetId, ".FILT.vcf", sep="")
	VARName <- paste(platDir, currSetId,"/", currSetId, ".VAR.vcf", sep="")
	SOMAName <- paste(platDir, currSetId,"/", currSetId, ".SOMA.vcf", sep="")
	
	#somatic files
	confTotalName <- paste(platDir, currSetId,"/", currSetId, ".", mutType,".somatic.vcf", sep="")
	confTotalOutput <- paste(platDir, currSetId,"/", currSetId, ".", mutType, ".somatic.txt", sep="")
	
	#SNP (germline variant) files
	germTotalName <- paste(platDir, currSetId,"/", currSetId, ".germline.vcf", sep="")
	germTotalOutput <- paste(platDir, currSetId,"/", currSetId, ".germline.txt", sep="")
	
	#biopsy order lists and therefore indexes to remove
	biopList <- system(command = paste("grep \"#CHROM\" ", platDir, currSetId, "/", currSetId, vcfName, sep=""), intern = TRUE, wait = TRUE)
	biopList <- strsplit(biopList, split = "\t")
	biopList <- biopList[[1]]
	biopList <- biopList[10:length(biopList)]
	write.table(biopList, file=paste(platDir, currSetId, "/", currSetId, ".bioIDorder.txt", sep=""), sep = "\t", row.names = FALSE, quote = FALSE)
	normalIndex <- (which(biopList==subSample[1, "normalID"]))-1
	
	system(command = paste("mkdir ", outAnno, currSetId, sep=""))
	
	#subset exclusion list
	colIndexes <- c(0:(nrow(subSample)-1))
	remList <- subSample[subSample[["retain"]]==2, "sampleID"]
	if(length(remList)==0){
	  keepList <- colIndexes
	  write.table(biopList, file=paste(outAnno, currSetId, "/", currSetId, ".vcfbioID.txt", sep=""), sep = "\t", row.names = FALSE, quote = FALSE)
	}else{
	  keepList <- colIndexes[-(which(biopList %in% remList))]
	  biopListTemp <- biopList[-(which(biopList %in% remList))]
	  write.table(biopListTemp, file=paste(outAnno, currSetId, "/", currSetId, ".vcfbioID.txt", sep=""), sep = "\t", row.names = FALSE, quote = FALSE)
	}
	
	#remove normal column from keep list
	keepListNoNorm <- keepList[-which(keepList == normalIndex)]

	#prepare indexes for total variant output
	counterTemp <- 1
	totalIndexStrings <- as.list(NA)
	for(k in keepListNoNorm){
		totalIndexStrings[[counterTemp]] <- paste(" ((GEN[", k,"].NR > 9) & (GEN[", k,"].NV > 0)) | ", sep="")
		counterTemp <- counterTemp +1
	}
	totalIndexStrings[[length(totalIndexStrings)]] <- substr(totalIndexStrings[[length(totalIndexStrings)]], 1,  (nchar(totalIndexStrings[[length(totalIndexStrings)]]) - 2 ))
	
	#prepare indexes for extractFields command
	counterTemp <- 1
	extractStringsNR <- as.list(NA)
	for(k in keepList){
		extractStringsNR[[counterTemp]] <- paste(" \"GEN[", k,"].NR\" ", sep="")
		counterTemp <- counterTemp +1

	}
	counterTemp <- 1
	extractStringsNV <- as.list(NA)
	for(k in keepList){
		extractStringsNV[[counterTemp]] <- paste(" \"GEN[", k,"].NV\" ", sep="")
		counterTemp <- counterTemp +1
	}
	

	# 1 .filter by FILTER field
	filtVarCommand <- paste("cat ", subSample[1, "sampleInfo"], " | java -jar ~/bin/SnpSift.jar filter \"( ( (FILTER ='PASS') | (FILTER ='alleleBias') | (FILTER ='HapScore') | (FILTER ='SC') | (FILTER ='badReads') | (FILTER ='SC;alleleBias') | (FILTER ='HapScore;alleleBias') | (FILTER ='HapScore;SC') ) )\" > ", FILTName, sep="")
	system(command=filtVarCommand)
	
	# 2. annotate variant types
	annoCommand <- paste("java -jar bin/SnpSift.jar varType ", FILTName," > ", VARName, sep="")
	system(command=annoCommand)
  
  
	#### somatic files ####
	
	# 3. filter for somatic mutations (not in normal)
	if(mutType == "snv"){
	  somaticCommand <- paste("cat ", VARName, " | java -jar ~/bin/SnpSift.jar filter \"( ( (exists SNP) & (GEN[", normalIndex,"].NV < 2) & (GEN[", normalIndex,"].NR < 100) & (GEN[", normalIndex,"].NR > 9) ) | ( (exists SNP) & (GEN[", normalIndex,"].NV < 3) & (GEN[", normalIndex,"].NR > 99) ) )\" > ", SOMAName, sep="")
	  system(command=somaticCommand)
	}else if(mutType == "indel"){
	  somaticCommand <- paste("cat ", VARName, " | java -jar ~/bin/SnpSift.jar filter \"( ( (exists DEL) & (GEN[", normalIndex,"].NV = 0)) | ( (exists INS) & (GEN[", normalIndex,"].NV = 0))  )\" > ", SOMAName, sep="")
	  system(command=somaticCommand)
	}
	
	# 4. filter by read depth (>9X for ANY samples), this is the final somatic variants vcf
	depthCommand5 <- paste("cat ", SOMAName, " | java -jar ~/bin/SnpSift.jar filter \"( ( ", paste(totalIndexStrings, collapse=" ") ,"))\" > ", confTotalName, sep="")
	system(command=depthCommand5)
	

	#### germline files ####

	# 4b. produce non-fitered germline file, this is the final germline vcf
	if(mutType == "snv"){
	  depthCommand3 <- paste("cat ", VARName, " | java -jar ~/bin/SnpSift.jar filter \"( (exists SNP) & (GEN[", normalIndex,"].NV > 0) & (GEN[", normalIndex,"].NR > 9) )\" > ", germTotalName, sep="")
	  system(command=depthCommand3)
	}else if(mutType == "indel"){
	  depthCommand3 <- paste("cat ", VARName, " | java -jar ~/bin/SnpSift.jar filter \"( ((exists DEL) & (GEN[", normalIndex,"].NV > 20) & (GEN[", normalIndex,"].NR > 49)) | ((exists INS) & (GEN[", normalIndex,"].NV > 20) & (GEN[", normalIndex,"].NR > 49)))\" > ", germTotalName, sep="")
	  system(command=depthCommand3)
	}
	
	
	#### vcf to text file conversion ####
	
	extractCommand2 <- paste("java -jar ~/bin/SnpSift.jar extractFields ", confTotalName, " \"CHROM\" \"POS\" \"REF\" \"ALT\" ", paste(extractStringsNR, collapse=" "), " ", paste(extractStringsNV, collapse=" "), " > ", confTotalOutput, sep="")
	system(command=extractCommand2)
		
	extractCommand5 <- paste("java -jar ~/bin/SnpSift.jar extractFields ", germTotalName, " \"CHROM\" \"POS\" \"REF\" \"ALT\" ", paste(extractStringsNR, collapse=" "), " ", paste(extractStringsNV, collapse=" "), " > ", germTotalOutput, sep="")
	system(command=extractCommand5)

	#tidy up
	system(command=paste("rm ", FILTName, sep=""))
	system(command=paste("rm ", VARName, sep=""))
	system(command=paste("rm ", SOMAName, sep=""))
}
	
  #annotate somatic file with annoVar
  makeAnnovar(subSample, paste(".", mutType, ".somatic.txt", sep=""), paste(".", mutType, ".somatic", sep=""), platDir)
  makeAnnovar(subSample, ".germline.txt", ".germline", platDir)
	
  system(command=paste("rm ", confTotalOutput, sep=""))
  system(command=paste("rm ", germTotalOutput, sep=""))
  
  #move annotated files to new directory
  annoVcf <- paste(platDir, currSetId,"/", currSetId, ".", mutType, ".somatic.annoVar.variant_function.txt", sep="")
  annoVcfNew <- paste(outAnno, currSetId,"/", currSetId, ".", mutType, ".somatic", ".annoVar.variant_function.txt", sep="")
  
  annoVcfEx <- paste(platDir, currSetId,"/", currSetId, ".", mutType, ".somatic.annoVar.exonic_variant_function.txt", sep="")
  annoVcfExNew <- paste(outAnno, currSetId,"/", currSetId, ".", mutType, ".somatic.annoVar.exonic_variant_function.txt", sep="")
   
  annoSNPTotal <- paste(platDir, currSetId,"/", currSetId, ".germline.annoVar.variant_function.txt", sep="")
  annoSNPTotalNew <- paste(outAnno, currSetId,"/", currSetId, ".germline.annoVar.variant_function.txt", sep="")
   
  annoSNPTotalEx <- paste(platDir, currSetId,"/", currSetId, ".germline.annoVar.exonic_variant_function.txt", sep="")
  annoSNPTotalExNew <- paste(outAnno, currSetId,"/", currSetId, ".germline.annoVar.exonic_variant_function.txt", sep="")
   
  #rename files with .txt prepend
  system(command=paste("mv", annoVcf, annoVcfNew))
  system(command=paste("mv", annoVcfEx, annoVcfExNew))
  system(command=paste("mv", annoSNPTotal, annoSNPTotalNew))
  system(command=paste("mv", annoSNPTotalEx, annoSNPTotalExNew))
}



