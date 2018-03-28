#### subroutines ####

#parse log file
sedLogFile <- function(logFile, parseFile){
	system(command=paste("gsed -n '/Apomorphy lists:/,/Character\ diagnostics:/p' ", logFile, " > ", parseFile, sep=""))
}

#adjust/correct parsed file format
correctFile <- function(inFile, outFile){
	dataIn <- read.table(file=inFile, stringsAsFactors=FALSE, fill=TRUE)
	dataIn <- dataIn[-c(1:3),]
	dataIn <- dataIn[-nrow(dataIn), ]
	
	for(currRow in 1:nrow(dataIn)){
		if(!grepl("node", dataIn[currRow,1])){
			dataIn[currRow,c(4:9)] <- dataIn[currRow,c(1:6)]
			dataIn[currRow,c(1:3)] <- currentNode
		}else{currentNode <- as.vector(dataIn[currRow,1:3])}
	}
	write.table(dataIn, file=outFile, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
}



#### main program ####

sampleList <- read.csv(file="~/Projects/IBDproject/masterSampleList.csv", header=TRUE, stringsAsFactors=FALSE, fill=TRUE)
setNames <- unique(sampleList[["setID"]])

holdingDir <- "~/Projects/IBDproject/5.phylogenetics/"


for(i in 1:length(setNames)){
	#subset current sample names from list
	subSample <- subset(sampleList, sampleList[1]== setNames[i])
	
	#assign file names
	logFileName <- paste(holdingDir, subSample[1,1],"/", subSample[1,1],".log", sep="")
	logFileNameCorrected <- paste(holdingDir, subSample[1,1],"/", subSample[1,1],".corr.log", sep="")
	parseFileName <- paste(holdingDir, subSample[1,1],"/", subSample[1,1],".HomoTab.txt", sep="")
	correctFileName <- paste(holdingDir, subSample[1,1],"/", subSample[1,1],".HomoTabCorr.txt", sep="")
	outFinal <- paste(holdingDir, subSample[1,1],"/", subSample[1,1],".finalHomo.txt", sep="")
	outFinalAnno <- paste(holdingDir, subSample[1,1],"/", subSample[1,1],".homoAnno.txt", sep="")

	#rename set names in log file to remove spaces in names
	sedStrings <- c(0)
	for(j in 1:nrow(subSample)){
		strSample <- strsplit(subSample[j,2], "_")
		if(length(strSample[[1]])==2){
      strSample <- paste(strSample[[1]][1], strSample[[1]][2])
		}
		if(length(strSample[[1]])==3){
		  strSample <- paste(strSample[[1]][1], strSample[[1]][2], strSample[[1]][3])
		}
		if(length(strSample[[1]])==4){
		  strSample <- paste(strSample[[1]][1], strSample[[1]][2], strSample[[1]][3], strSample[[1]][4])
		}
		if(length(strSample[[1]])==5){
		  strSample <- paste(strSample[[1]][1], strSample[[1]][2], strSample[[1]][3], strSample[[1]][4], strSample[[1]][5])
		}
    sedStrings[j] <- paste("-e 's/ ", strSample," / ",subSample[j,2]," /g'", sep="")
	}
	sedCommand <- paste("gsed -E ", paste(sedStrings, collapse=" "), " ", logFileName," > ", logFileNameCorrected, sep="")
	system(command=sedCommand)

	#parse out character table from log file	
	sedLogFile(logFileNameCorrected, parseFileName)
		
	#now correct file for spacing and headers
	correctFile(parseFileName, correctFileName)
	
	system(command=paste("rm", logFileNameCorrected))
	system(command=paste("rm", parseFileName))


	#now analyse corrected homoplasy table
	corrHomoplasy <- read.table(file=correctFileName, sep="\t", stringsAsFactors=FALSE, header=FALSE)

	#remove unwanted columns
	corrHomoplasy <- corrHomoplasy[-c(2,5,8)]
	corrHomoplasy <- subset(corrHomoplasy, corrHomoplasy$V6!=1)
	if(nrow(corrHomoplasy)==0){
	  print(paste("#### no homoplasy in sample", subSample[1,1], "####"))
	  next()
	}
	
	#output file
	write.table(corrHomoplasy, file=outFinal, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
	
	#now get unique variants by removing duplicates, annotate and export file
	corrHomoplasy <- corrHomoplasy[order(corrHomoplasy[3]),]
	removeVect <- c()
	counter <- 1
	for(currRow in 2:nrow(corrHomoplasy)){
		if(corrHomoplasy[currRow, 3] == corrHomoplasy[(currRow-1), 3]){
			removeVect[counter] <- currRow
			counter <- counter + 1
		}
	}
	corrHomoplasy <- corrHomoplasy[-removeVect,]
	corrHomoplasy <- corrHomoplasy[-c(1:2)]
	
	#get variant list
	varInName <- paste(holdingDir, subSample[1,1],"/", subSample[1,1], ".snv.annoVar.retained.txt", sep="")
	varList <- read.table(file=varInName, sep="\t", stringsAsFactors=FALSE, header=TRUE)
	varList[1] <- c(1:nrow(varList))
	getGeneCol <- which(names(varList)=="gene")
  
  
	#find gene names, positions etc from var list and make homoplasy variant list
	for(currRow in 1:nrow(corrHomoplasy)){
	  tempGene <- varList[as.numeric(corrHomoplasy[currRow, 1]), getGeneCol]
	  tempGene <- strsplit(tempGene, ":")
    corrHomoplasy[currRow, 2] <- tempGene[[1]][1]
	 	corrHomoplasy[currRow, 3] <- varList[as.numeric(corrHomoplasy[currRow, 1]), "chr"]
	  corrHomoplasy[currRow, 4] <- varList[as.numeric(corrHomoplasy[currRow, 1]), "start"]
	}
	
	system(command=paste("rm", correctFileName))
	system(command=paste("rm", outFinal))
  
	write.table(corrHomoplasy, file=outFinalAnno, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
}

