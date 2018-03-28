############## libraries ################
library(apTreeshape)
library(ape)
library(phangorn)

#get sample list information
sampleList <- read.csv(file="~/Projects/IBDproject/masterSampleList.csv", header=TRUE, stringsAsFactors=FALSE)
setNames <- unique(sampleList[["setID"]])
sampleNames <- unique(sampleList[[1]])
noSets <- length(sampleNames)

holdingDir <- "5.phylogenetics/CNAfilt/"
namePrepended <- ".tre"
namePrepended2 <- ".allTreeSearch.tre"

######## get the most parsimonious and shape stats ########
for(i in 1:noSets){
  subList <- subset(sampleList, sampleList[1]== sampleNames[i])
  
  treFile <- paste(subList[1, "FileHolder1"], holdingDir, subList[1,1],"/", subList[1,1], namePrepended, sep="")
  treFileTotal <- paste(subList[1, "FileHolder1"], holdingDir, subList[1,1],"/", subList[1,1], namePrepended2, sep="")
  
  normalIndex <- subList[1,"normalIndex"]+1
  normalName <- subList[normalIndex,"sampleID"]
  histFile <- paste(subList[1, "FileHolder1"], holdingDir, subList[1,1],"/", subList[1,1], ".hist.txt", sep="")
  
  treeList <- as.list(0)
  
  #get specific tree set from file
  if(file.exists(treFileTotal)){
    treeList <- read.nexus(file=treFileTotal)
    histTab <- read.table(file=histFile, header=FALSE, sep="\t")
  }else if(file.exists(treFile)){
    treeList <- read.nexus(file=treFile)
    histTab <- NA
  }else{
    next
  }
  
  #no of trees
  if(length(is.binary.tree(treeList)) > 1){
    notrees <- length(treeList)
  }else{
    print(paste("#### only one tree for set", sampleNames[i],"####"))
    next()
  }
  
  #setup results table
  resultsTable <- as.data.frame(matrix(NA, nrow=notrees, ncol=3))
  names(resultsTable) <- c("tree length", "colless test (yule)", "colless test (PDA)")
  row.names(resultsTable) <- paste("tree", c(1:notrees))
  
  #for each tree get stats
  for(j in 1:notrees){
    print(paste("calculating stats for tree", j, "of set", subList[1,1])) 
    
    #convert to single tree object
    phyloTree <- treeList[[j]]
    
    #resolve polytomies
    phyloTree <- multi2di(phyloTree, random = TRUE)
    
    #get tree length
    resultsTable[j,1] <- sum(phyloTree$edge.length)
    
    #drop normal sample
    treeTemp <- drop.tip(phyloTree, normalName)
    
    #convert to treeshape object
    phyloShape <- as.treeshape(treeTemp)
    
    #perform stats (silently)
    dummyFile <- file()
    sink(file=dummyFile)
    if(length(phyloShape$names) < 5){
      yuleTemp <- colless.test(phyloShape, model="yule", n.mc=1000, alternative="greater")
    }else{
      yuleTemp <- likelihood.test(phyloShape, model="yule", alternative = "greater")
    }
    pdaTemp <- colless.test(phyloShape, model="pda", n.mc=1000, alternative="greater")
    sink()
    close(dummyFile)
    
    #populate table
    resultsTable[j,2] <- yuleTemp$p.value
    resultsTable[j,3] <- pdaTemp$p.value
  }
  
  #reorder table
  resultsTable <- resultsTable[order(resultsTable[1]), ]
  #parseTreeList[i, treeTypes[eventType]] <- row.names(resultsTable)[1]
  
  write.table(resultsTable, file=paste(subList[1, "FileHolder1"], holdingDir, subList[1,1], "/", subList[1,1], ".treeStats.txt", sep=""), sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
  
  #output most parsimonous tree
  parseFile <- paste(subList[1, "FileHolder1"], holdingDir, subList[1,1],"/", subList[1,1], ".parseTree.tre", sep="")
  mostParseTree <- as.numeric(strsplit(rownames(resultsTable)[1], " ")[[1]][2])
  parseTree <- treeList[[mostParseTree]]
  write.tree(parseTree, file=parseFile)
  
  #plot graphs of distribution
  pdf(file=paste(subList[1, "FileHolder1"], holdingDir, subList[1,1],"/", subList[1,1], ".tree_p.value.pdf", sep=""), onefile=TRUE, width=5, height=5)
    plot(resultsTable[[2]], resultsTable[[1]], main="tree length vs p-value from yule", xlab="p-value from yule", ylab="tree length")
  dev.off()
}



