require(gtools)
require(nnls)
require(SomaticSignatures)
require(RColorBrewer)
require(vioplot)

# Function to classify mutations in terms of timing based on their presence/absence in individual samples
SampleType <- function(early.samples, late.samples, samples){
  tmp <- strsplit(samples, ':')[[1]]
  type <- 'NA'
  if(length(intersect(tmp, early.samples)) > 0){
    type <- 'pre'
  } else if ((length(intersect(tmp, late.samples)) < length(late.samples)) & (length(intersect(samples, late.samples)) > 0)){
    type <- 'early'  
  } else if (length(intersect(tmp, late.samples)) == length(late.samples)){
    type <- 'late'
  }
  return(type)
}

#########################################################################################################################
### User specified parameters; Update as applicable ###
setwd('/PROJECT_PATH/')
# Replace with paths to the relevant files in local #
hg19.path <- "/USER_PATH/ucsc.hg19.fasta"
cosmic.signatures.path <- "/USER_PATH/signatures_probabilities.txt"
### User specified parameters; Update as applicable ###
#########################################################################################################################

msl <- read.table(file = 'data/master.sample.list.csv', sep = ',', header = T)
sets <- sort(unique(msl$sampleInfo[which((msl$retain == 1) & (msl$noSamples > 2))]))

# calling criteria
variant.read.cutoff <- 2
vaf.cutoff <- 0.05
coverage.cutoff <- 10
fdr.threshold <- 0.05

# mutation filtering criteria
bases <- c('A', 'C', 'G', 'T')

# signature information
contexts <- paste(rep(c('CA', 'CG', 'CT', 'TA', 'TC', 'TG'), each = 16), 
                  rep(apply(permutations(n = 4, v=c('A', 'C', 'G', 'T'), r = 2, repeats.allowed = T), 1, function(x) paste(x, collapse = '.')), 6))
context.names <- sapply(as.character(contexts), function(x) paste(gsub('\\.', substr(x, 1, 1), strsplit(x, ' ')[[1]][2]), gsub('\\.', substr(x, 2, 2), strsplit(x, ' ')[[1]][2]), sep = '>'))
prc <- read.table(cosmic.signatures.path, sep = '\t', header = T, stringsAsFactors = F)
unordered.contexts <- paste(gsub('>', '', prc$Substitution.Type), paste(substr(prc$Trinucleotide, 1, 1), substr(prc$Trinucleotide, 3, 3), sep = '.'))
prc <- prc[order(unordered.contexts), ]
rownames(prc) <- contexts
prc.proven <- prc[, paste('Signature.', 1:30, sep = '')]
tnf <- read.table("data/trinucleotide.frequencies.tsv", header = T, stringsAsFactors = F)
tnf$genome.to.exome <- (tnf$exome / tnf$genome)
frequencies <- tnf$genome.to.exome[match(as.character(sapply(contexts, function(x) paste(substr(x,4,4), substr(x,1,1), substr(x,6,6), sep = ''))),
                                         tnf$type)]
prc.proven.exome <- prc.proven
for(i in ncol(prc.proven.exome)){
  prc.proven.exome[, i] <- prc.proven.exome[, i] * frequencies
  prc.proven.exome[, i] <- prc.proven.exome[, i] / sum(prc.proven.exome[, i])
}

# active signatures in CRC
active.signature.numbers <-  c(1, 2, 5, 6, 10, 13, 17)
active.signatures <- paste('Signature.', active.signature.numbers, sep = '')

# signature annotations
signature.annotations <- read.table(file = "data/signature.annotations.csv", 
                                    sep = ',', col.names = c('process', 'association', 'type', 'clock.like'))
signature.annotations$type <- sapply(1:nrow(signature.annotations), function(x) if(signature.annotations$association[x] == 'Unknown aetiology') 
{paste(signature.annotations$process[x], signature.annotations$association[x])} else {paste(signature.annotations$association[x])})

# signature colour mapping
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector[c(9, 10, 15)] <- col_vector[26 + c(3, 4, 7)]
color.map <- data.frame(signature = colnames(prc.proven), 
                        col = col_vector[match(signature.annotations$type, unique(signature.annotations$type))])

sample.types <- c('pre', 'early', 'late')

#########################################################################################################################
## A) Read in mutation data and assign mutations to individual samples and classify mutations in terms of timing

# read in CA-CRC data
for(set in sets){
  print(set)
  
  msl.sub <- msl[which((msl$sampleInfo == set) & (msl$retain == 1)), ]
  early.samples <- msl.sub$sampleID[which(msl.sub$type %in% c('proximal', 'polyp'))]
  late.samples <- msl.sub$sampleID[which(msl.sub$type == 'cancer')]
  somatic.samples <- union(early.samples, late.samples)
  
  set.name <- msl$setID[match(set, msl$sampleInfo)]
  df <- read.table(file = paste('data/', set.name, '.snv.total.annotated.txt', sep = ''),
                   sep = '\t',
                   header = T)
  
  # assign mutations to samples
  variant.reads <- as.data.frame(df[, match(paste(somatic.samples, '.NV', sep = ''), colnames(df))])
  vafs <- as.data.frame(df[, match(paste(somatic.samples, '.VAF', sep = ''), colnames(df))])
  coverage <- as.data.frame(df[, match(paste(somatic.samples, '.NR', sep = ''), colnames(df))])
  df$samples <- sapply(1:nrow(df), 
                       function(x) paste(somatic.samples[which((variant.reads[x,] > variant.read.cutoff) & (vafs[x,] > vaf.cutoff) & (coverage[x,] > coverage.cutoff))], collapse = ':'))
  
  # add mutation contexts
  vr <- VRanges(seqnames = Rle(paste(df$chr, sep = '')), 
                ranges = IRanges(start = df$start, end = df$end),
                ref = as.character(df$ref),
                alt = as.character(df$alt))
  file = hg19.path
  fa <- open(FaFile(file, sprintf("%s.fai", file)))
  vr.tmp <- mutationContext(vr = vr, ref = fa, k = 3)
  df$context <- paste(vr.tmp$alteration, vr.tmp$context)
  df$sample.type <- sapply(1:nrow(df), function(x) SampleType(early.samples = early.samples,
                                                              late.samples = late.samples,
                                                              samples = df$samples[x]))
  
  assign(paste('df.', set, sep = ''), df)
}

# Count mutations of each context for each mutation timing type (i.e. 'pre', 'early' and 'late')
for(sample.type in sample.types){
  sample.contexts <- data.frame(row.names = contexts)
  for(set in sets){
    print(set)
    df <- get(paste('df.', set, sep = ''))
    
    sub <- df[which(df$sample.type == sample.type),]
    sample.contexts <- cbind(sample.contexts,
                             as.numeric(table(factor(sub$context, levels = contexts))))
    colnames(sample.contexts)[ncol(sample.contexts)] <- set
  }
  assign(paste('sample.contexts.', sample.type, sep = ''), sample.contexts)
}

#########################################################################################################################
## B) Assign mutational signatures

# Assign signatures
sample.types <- c('pre', 'early', 'late')
for(sample.type in sample.types){
  sample.contexts <- get(paste('sample.contexts.', sample.type, sep = ''))
  
  signatures <- prc.proven.exome[, match(active.signatures, colnames(prc.proven.exome))]
  sample.signatures <- as.data.frame(lapply(1:ncol(sample.contexts), function(y) nnls(as.matrix(signatures), sample.contexts[,y])$x))
  fitted.signatures <- as.data.frame(lapply(1:ncol(sample.contexts), function(y) nnls(as.matrix(signatures), sample.contexts[,y])$fitted))
  
  rownames(sample.signatures) <- active.signatures
  colnames(sample.signatures) <- colnames(sample.contexts)
  sample.normalised.signatures <- apply(sample.signatures, 2, function(x) x / sum(x))
  assign(paste('sample.signatures.', sample.type, sep = ''), sample.signatures)
  assign(paste('fitted.signatures.', sample.type, sep = ''), fitted.signatures)
  assign(paste('sample.normalised.signatures.', sample.type, sep = ''), sample.normalised.signatures)
}

#########################################################################################################################
## C) Compare mutational signautres for mutations in different timing classes

# tests
tests <- as.data.frame(do.call('rbind', lapply(1:length(active.signatures), function(x) t(combn(sample.types, 2)))))
tests$signature <- rep(active.signatures, each = 3)
tests$p.value = 1
for(i in 1:nrow(tests)){
  df.1 <- get(paste('sample.normalised.signatures.', tests[i,1], sep = ''))
  df.2 <- get(paste('sample.normalised.signatures.', tests[i,2], sep = ''))
  tests$p.value[i] <- wilcox.test(df.1[match(tests$signature[i], rownames(df.1)), ],
                                  df.2[match(tests$signature[i], rownames(df.2)), ],
                                  paired = T)$p.value
}
tests <- tests[order(tests$p.value),]
tests$stat <- sapply(1:nrow(tests), function(x) (tests$p.value[x] * nrow(tests)) / length(which(tests$p.value <= tests$p.value[x])))
tests$q.value <- sapply(1:nrow(tests), function(x) min(tests$stat[x:nrow(tests)]))
time.evolution <- tests

# plots
pdf(file = 'time_evolution_plots/absolute.signature.evolution.pdf',
    width = 11.69, height = 8.27)
par(mfrow = c(2, 4), oma = c(4,4,4,4))
for(sig in active.signatures){
  vioplot(as.numeric(sample.signatures.pre[match(sig, rownames(sample.signatures.pre)), ]), 
          as.numeric(sample.signatures.early[match(sig, rownames(sample.signatures.early)), ]),
          as.numeric(sample.signatures.late[match(sig, rownames(sample.signatures.late)), ]),
          col = as.character(color.map$col[match(sig, color.map$signature)]),
          names = c('Pre', 'Early', 'Late'))
  title(ylab = 'Signature Exposure', main = gsub('\\.', ' ', sig))
}
dev.off()


pdf(file = 'time_evolution_plots/relative.signature.evolution.pdf',
    width = 11.69, height = 8.27)
par(mfrow = c(2, 4), oma = c(4,4,4,4))
for(sig in active.signatures){
  vioplot(sample.normalised.signatures.pre[match(sig, rownames(sample.normalised.signatures.pre)), ], 
          sample.normalised.signatures.early[match(sig, rownames(sample.normalised.signatures.early)), ],
          sample.normalised.signatures.late[match(sig, rownames(sample.normalised.signatures.late)), ],
          col = as.character(color.map$col[match(sig, color.map$signature)]),
          names = c('Pre', 'Early', 'Late'))
  title(ylab = 'Relative Signature Exposure', main = gsub('\\.', ' ', sig))
}
dev.off()



