#get sample list information
sampleList <- read.csv(file="~/Projects/IBDproject/masterSampleList.csv", header=TRUE, stringsAsFactors=FALSE)
setNames <- unique(sampleList[["setID"]])

vcfDir <- "~/Projects/IBDproject/4.platypusIndels/"

platDir <- "~/Projects/IBDproject/4.platypusCalls/annotatedSNVs/"

phyloDir <- "~/Projects/IBDproject/4.phylogenetics/subsetBiopsies/"

#vafFilter <- 0.1

#perform pipeline for each sample
for(j in 1:length(setNames)){
  currSetId <- setNames[j]
	print(paste("#### analysing sample ", currSetId, " ####",sep=""))
	
	system(command = paste("mkdir ", phyloDir, currSetId, "/", sep=""))
	
	#setup input/output names
  subSample <- sampleList[sampleList[[1]]==currSetId, ]
  normID <- subSample[1, "normalID"]
  
  biopList <-  read.table(file=paste(platDir, currSetId, "/", currSetId, ".vcfbioID.txt", sep=""), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  biopList <- biopList[[1]]
  bioTotal <- biopList 
  normName <- subSample[1, "normalID"]
  biopList <- biopList[-(which(biopList == normName))]
  noSamples <- length(biopList)
  
  annoTable <- read.table(file=paste(platDir, currSetId, "/", currSetId, namePrep, sep=""), sep="\t", header = TRUE, stringsAsFactors = FALSE)
  
  #mark each variant by phylogenetic locations and VAF
  annoTable[ncol(annoTable)+1] <- 0
  names(annoTable)[ncol(annoTable)] <- "phyloLoc"
  
  annoTable[ncol(annoTable)+1] <- 0
  names(annoTable)[ncol(annoTable)] <- "VAForder"
  
  for(currRow in 1:nrow(annoTable)){
    #get phylogenetic location
    assessRow <- annoTable[currRow, paste(biopList, ".VAF", sep="")] !=0
    assessTable <- table(assessRow)
    
    #add VAF information (sum)
    tempVal <- annoTable[currRow, paste(biopList, ".VAF", sep="")]
    tempVal <- tempVal[!is.na(tempVal)]
    tempVal <- tempVal[tempVal!=0]
    annoTable[currRow, "VAForder"] <- sum(tempVal)
    
    #assess if variant present at low VAF in normal
    normAssess <- annoTable[currRow, paste(normName, ".VAF", sep="")]
    if(!is.na(normAssess)){
      if(normAssess > 0){
        annoTable[currRow, "phyloLoc"] <- 0
        next()
      }
    }
    
    #if trunkal mark as 1 etc
    if(as.integer(assessTable["TRUE"])==noSamples){
      annoTable[currRow, "phyloLoc"] <- 1
    }
    if(as.integer(assessTable["TRUE"])==1){
      annoTable[currRow, "phyloLoc"] <- 3
    }
    if(as.integer(assessTable["TRUE"])!=1 & as.integer(assessTable["TRUE"])!=noSamples){
      annoTable[currRow, "phyloLoc"] <- 2
    }
    
  }
  
  sortTable <- annoTable[order( annoTable[["phyloLoc"]], -annoTable[["VAForder"]] ), ]
  
  #output sorted table
  write.table(sortTable, file=paste(platDir, currSetId, "/", currSetId, ".snv.annoVar.sorted.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  
  ########## prepare data for phylogenetic tree construction
  
  #convert VAFs to binary table
  confData <- sortTable[, paste(bioTotal, ".VAF", sep="")]
  confData[confData < vafFilter] <- 0
  confData[confData >= vafFilter] <- 1
  confData[is.na(confData)] <- 0
  
  keepRows <- c()
  counter <- 1
  for(currRow in 1:nrow(confData)){
    sumVAF <- sum(confData[currRow, ])
    normVAF <- confData[currRow, paste(normID, ".VAF", sep="")]
    if(sumVAF != 0 & normVAF != 1){
      keepRows[counter] <- currRow
      counter <- counter + 1
    }
  }
  
  #remove filtered variants and save filtered tables
  confData <- confData[keepRows, ]
  tempTable <- sortTable[-keepRows, ]
  write.table(tempTable, file=paste(phyloDir, currSetId, "/", currSetId, ".snv.annoVar.filtered.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  keepTable <- sortTable[keepRows, ]
  write.table(keepTable, file=paste(phyloDir, currSetId, "/", currSetId, ".snv.annoVar.retained.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  
  #store strings with names in list
  sampleData <- as.list(0)
  for(z in 1:ncol(confData)){	
    sampleData[[z]] <- confData[[z]]
    names(sampleData)[[z]] <- bioTotal[z]
  }
  noSamples <- noSamples + 1
  counter<-1
  asStrings <- as.list(0)
  for(k in seq(1,(2*(noSamples)),2)){
    asStrings[[k]] <- as.character(bioTotal[counter])
    asStrings[[k+1]] <- paste(sampleData[[counter]], collapse='')
    counter <- counter+1
  }
  
  #get assign all file locations and parsomony parameters
  
  noTaxa <- (noSamples)
  taxaNames <- paste(bioTotal, collapse=" ")
  noCharacters <- nchar(asStrings[[2]][[1]])
  matrixStrings <- paste(asStrings, collapse="\n")
  ascestState <- 0
  for(w in 1:(noCharacters-1)){ascestState <- paste(ascestState,0,sep="")}
  
  upperHIvalue <- round(noCharacters + ((noCharacters / 100) * 55), 0)
  
  #dir names for nexux control file
  outputLoc <- paste(subSample[1, "FileHolder2"], "7.phylogenetics/", currSetId, "/", currSetId, sep="")
  
  treeFile <- paste(outputLoc, ".tre", sep="")
  treeFileAll <- paste(outputLoc, ".allTreeSearch.tre", sep="")
  logFile <- paste(outputLoc, ".log", sep="")
  logFileAll <- paste(outputLoc, ".allTrees.log", sep="")
  logFileBoot <- paste(outputLoc, ".boot.log", sep="")
  histogramFile <- paste(outputLoc, ".hist.txt", sep="")
  
  nexusTreeFileLoc <- paste(outputLoc, ".nex", sep="")
  nexusBootFileLoc <- paste(outputLoc, ".boot.nex", sep="")
  nexusAllTreesFile <- paste(outputLoc, ".allTrees.nex", sep="")
  
  outputLocNex <- paste(phyloDir, currSetId, "/", currSetId, sep="")
  
  #put all strings into nexus format for tree production
  totalStringsTree <- paste("#NEXUS
                            begin paup;
                            set criterion=parsimony autoclose=yes warntree=no warnreset=no monitor=no warntsave=no maxtrees=1000 increase=no;
                            log start file= ", logFile,"  replace;
                            
                            begin taxa;
                            dimensions ntax=", noTaxa,";
                            taxlabels
                            ", taxaNames,";
                            end;
                            
                            begin characters;
                            dimensions nchar= ", noCharacters," ;
                            format symbols = \"01\";
                            matrix\n", matrixStrings,";
                            
                            end;
                            
                            begin assumptions;
                            options deftype = unord;
                            
                            usertype statetransitionsmatrix (stepmatrix)= 2
                            0 1 
                            . i 
                            1 . 
                            ;
                            
                            end;
                            
                            hsearch addseq=simple nreps=20 hold=1 rearrlimit=10000000 nbest=1000;
                            
                            outgroup ", normName," /only;
                            
                            roottrees outroot=monophyl rootmethod=outgroup;
                            
                            describetrees /plot=both root=outgroup brlens=yes chglist=yes apolist=yes diag=yes homoplasy=yes xout=both;
                            
                            savetrees from=1 to=1000 file= ", treeFile,"  savebootp=Both format=altnex brlens=yes root=yes append=yes;
                            
                            log stop;
                            
                            end;
                            
                            quit;", sep="")
  
  
  #write nexus file
  lapply(totalStringsTree, write, paste(outputLocNex,".nex",sep=""), append=FALSE)
  
  
  #put all strings into nexus format for all tree production
  totalStringsTreeAll <- paste("#NEXUS
                               begin paup;
                               set criterion=parsimony autoclose=yes warntree=no warnreset=no monitor=no warntsave=no maxtrees=1000 increase=no;
                               log start file= ", logFileAll,"  replace;
                               
                               begin taxa;
                               dimensions ntax=", noTaxa,";
                               taxlabels
                               ", taxaNames,";
                               end;
                               
                               begin characters;
                               dimensions nchar= ", noCharacters," ;
                               format symbols = \"01\";
                               matrix\n", matrixStrings,";
                               
                               end;
                               
                               begin assumptions;
                               options deftype = unord;
                               
                               usertype statetransitionsmatrix (stepmatrix)= 2
                               0 1 
                               . i 
                               1 . 
                               ;
                               
                               end;
                               
                               alltrees fdfile=", histogramFile," keep=", upperHIvalue,";
                               
                               outgroup ", normName," /only;
                               
                               roottrees outroot=monophyl rootmethod=outgroup;
                               
                               describetrees /plot=both root=outgroup brlens=yes chglist=yes apolist=yes diag=yes homoplasy=yes xout=both;
                               
                               savetrees from=1 to=1000 file= ", treeFileAll,"  savebootp=Both format=altnex brlens=yes root=yes append=yes;
                               
                               log stop;
                               
                               end;
                               
                               quit;", sep="")
  
  
  #write nexus file
  lapply(totalStringsTreeAll, write, paste(outputLocNex,".allTrees.nex",sep=""), append=FALSE)
  
  
  totalStringsBoot <- paste("#NEXUS
begin paup;
                            set criterion=parsimony autoclose=yes warntree=no warnreset=no monitor=no warntsave=no maxtrees=1000 increase=no;
                            log start file= ", logFileBoot,"  replace;
                            
                            begin taxa;
                            dimensions ntax=", noTaxa,";
                            taxlabels
                            ", taxaNames,";
                            end;
                            
                            begin characters;
                            dimensions nchar= ", noCharacters," ;
                            format symbols = \"01\";
                            matrix\n",matrixStrings,";
                            
                            end;
                            
                            begin assumptions;
                            options deftype = unord;
                            
                            usertype statetransitionsmatrix (stepmatrix)= 2
                            0 1 
                            . i 
                            1 . 
                            ;
                            
                            end;
                            
                            hsearch addseq=simple nreps=20 hold=1 rearrlimit=10000000;
                            
                            outgroup ", normName," /only;
                            
                            roottrees outroot=monophyl rootmethod=outgroup;
                            
                            bootstrap nreps=10000;
                            
                            log stop;
                            
                            end;
                            
                            quit;", sep="")
	
	
	#write nexus file
	lapply(totalStringsBoot, write, paste(outputLocNex,".boot.nex",sep=""), append=FALSE)

	#create shell scripts to run nexus files
	shellStrings <- paste("#!/bin/sh
	                      #$ -cwd
	                      #$ -V
	                      #$ -pe smp 1            # Request 2 CPU cores
	                      #$ -l h_rt=48:0:0      # Request 48 hour runtime
	                      #$ -l h_vmem=2G    # Request 12GB RAM / core, i.e. 24GB total
	                      
	                      ./bin/paup4b10-x86-linux -n ", nexusTreeFileLoc, "
	                      ./bin/paup4b10-x86-linux -n ",nexusBootFileLoc, "
	                      ./bin/paup4b10-x86-linux -n ", nexusAllTreesFile
, sep="")

	#write nexus file
	lapply(shellStrings, write, paste(outputLocNex,".runPAUP.sh",sep=""), append=FALSE)
}

