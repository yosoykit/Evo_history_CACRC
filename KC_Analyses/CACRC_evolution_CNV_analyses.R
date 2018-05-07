################################################
## Baker et al. 2018 												
## This script performs the analyses to create:						
## Figures 3B, 3C, 4 												
## Supplementary Tables 9, 10 										
## All corresponding Results in the Main text 
################################################

## Load LGD/HGD call data for 500 kbp binning, 6 bins blacklisted - total of 4401 bins. All specific folder names have been changed to generic placeholders  
## 81 samples total (LGD = 38, HGD = 23, CACRC = 7, U = 12, pseudopolyp = 1) - to be used lpWGS in CA-CRC paper Baker et al. 2018

## For final plots, remove all MSI+ patients' samples: STM005, Oxford_IBD5 (WES), IF3, IF5, IF7, IF8 (SNP), and 3611 (LPWGS) and sample 3305.H2 (bulk of A and B of that lesion)
library(vioplot)
library(gplots)
library(beeswarm)
library(forestplot)

run_size=81
total_binslowpass = 4401

setwd("/data/Results/finalbams")
outputdir<-"/data/Results/QDNAseq/CNAcalling"
plotdir<-"/data/Results/QDNAseq/CACRC_paper_Figures"
samples<-read.table('../../LGD_HGD_81_samples_final_withlesion.txt', colClasses = c("character"))[,3:6]
filenames<-paste(samples[,1],samples[,2],sep=".")[1:run_size]
MSIind = which(substr(filenames,1,4)=="3611")
bulk3305_h2 =which(filenames=="3305.H2") 
remove_LPWGS = c(MSIind,bulk3305_h2)
filenames = filenames[-remove_LPWGS]
samples = samples[-remove_LPWGS,]
patients = unique(substr(filenames,1,4))
## Load bin locations (chrom, start, end bp) for low pass data
bins_4401 = read.table("../../bin_locations.txt",header=T)[,2:4]
median.bp = rep(0,length(bins_4401[,1]))
for (i in 1:length(bins_4401[,1])){
	median.bp[i]=median(as.numeric(bins_4401[i,2:3]))
}
chroms_4401 = bins_4401[,1]


## load low pass data for 81 samples (n=77, removed MSI + patient 3611 for a total of 18 MSI- patients included in paper results)
load(file = paste(outputdir,"/","bin_locations500kbp.Rdata",sep=""),verbose=T)
load(file = paste(outputdir,"/","CNAcall_list_LGDHGD_500kbp.RData",sep=""),verbose=T)
tissue_type = read.table('../../LGD_HGD_81_samples_final_withlesion.txt')[-remove_LPWGS,5]

HGD_ind = which(tissue_type=="H")
C_ind = which(tissue_type=="C")
LGD_ind = which(tissue_type=='L')
### Compare LGD/HGD with U2 samples, pathologists gave differing diagnoses of LGD or HGD, and confirmed U
U_ind = c(which(tissue_type=="U2"), which(tissue_type=="U"))

CNA_freq = rep(0,length(filenames))

call_stats= call_stats_LGDHGD_500kbp
count=1
for (i in filenames){
	CNA_freq[count] = length(which(call_stats[[i]]['calls',]!=0))/length(call_stats[[i]]['calls',])
	count=count+1
}


### Analysess with Sequenza tumors (n=70 samples from 12 tumors) - use MSI- samples ONLY: (n=47 samples from 11 tumors/patients)

Exome_tumors = read.csv('../QDNAseq/WILLmasterSampleList.csv',header=T)

## Make list for all call values for each patient
MSIneg_ind=c(which(substr(Exome_tumors$Dx,1,1)=="B"), which(substr(Exome_tumors$Dx,1,3)=="MSI"))

exome_names = Exome_tumors$sampleID[-MSIneg_ind]
exome_dx = Exome_tumors$Dx[-MSIneg_ind ]
exome_set = Exome_tumors$setID[-MSIneg_ind]
exome_patients = unique(exome_set)
Exome_data_call_stats = list()
## these are the columns for the Exome segmented data"
#callnames=c("chromosome","start.pos",	"end.pos"	,"Bf",	"N.BAF",	"sd.BAF","depth.ratio",	"N.ratio",	"sd.ratio",	"CNt",	"A",	"B",	"LPP") 

for (i in 1:length(exome_names)){
	txt_dir = exome_set[i] 
	sample = exome_names[i]
	Exome_data_call_stats[[i]] = read.table(paste('/LGD_HGD_CACRC/sequenzaFiles/',txt_dir,'/',sample,'_segments.txt',sep=""),head=TRUE)	
}
names(Exome_data_call_stats) = exome_names

N_indE = which(exome_dx=="N")
HGD_indE = which(exome_dx=="H")
C_indE = which(exome_dx=="C")
LGD_indE = which(substr(exome_dx,1,1)=="L")
CNA_freqExome = rep(0,length(exome_names))

count=1
for (i in 1:length(exome_names)){
	segs=Exome_data_call_stats[[i]]$CNt	
	medianploidy = median(segs)
	total_size=sum(as.numeric(Exome_data_call_stats[[i]]$end.pos-Exome_data_call_stats[[i]]$start.pos))
	CNA_size=0 
	for ( j in 1: length(segs)){
		if (abs(segs[j]-medianploidy)>=1){
			CNA_size= CNA_size + (Exome_data_call_stats[[i]]$end.pos[j]-Exome_data_call_stats[[i]]$start.pos[j])
		}
	}
	CNA_freqExome[count] = CNA_size/total_size
	count=count+1
}


## SNP array tumors: 13 all cancer, ASCAT (n=13 from 13 tumors) ; Remove MSI + IF5, IF7, IF8, n=9 from 9 tumors/patients

load('/data/CGH/IBD/ascat/ascat_segments.RData',verbose=T)
snp_names = unique(ryan_ascat_segments$SampleID)
MSIneg_ind = c(which(snp_names=="IF5-T"),which(snp_names=="IF7-T"),which(snp_names=="IF8-T"), which(snp_names=="IF3-T"))
snp_names= snp_names[-MSIneg_ind]
CNA_freqSNP = rep(0,length(snp_names))
count=1
for (i in 1:length(snp_names)){
	pat_id=which(ryan_ascat_segments$SampleID==snp_names[i])
	segs = ryan_ascat_segments$nA[pat_id]+ryan_ascat_segments$nB[pat_id]
	medianploidy = median(segs)
	total_size=sum(as.numeric(ryan_ascat_segments$End[pat_id]-ryan_ascat_segments$Start[pat_id]))
	CNA_size=0 
	for ( j in 1: length(segs)){
		if (abs(segs[j]-medianploidy)>=1){
			CNA_size= CNA_size + (ryan_ascat_segments$End[pat_id[j]]-ryan_ascat_segments$Start[pat_id[j]])
		}
	}
	CNA_freqSNP[count] = CNA_size/total_size
	count=count+1
}

## ASCAT segments for 25 adenomatous polyps 
load('/data/ascat/AdInCa_unmatched_ascat_segments.RData',verbose=T)
polyp_names = c("1-D1","1-D2","2-D1","2-D2","3-D1","3-D2","4-D1","4-D2","5-D1","5-D2","6-D1","6-D2","7-D1","7-D2","8-D1","8-D2","9-D1","9-D2","10-D1","10-D2",
					"11-D1","11-D2","12-D1","12-D2","13-D1","13-D2","14-D1","14-D2","15-D1","15-D2","16-D1","16-D2","17-D1","17-D2","18-D1","18-D2","19-D1","19-D2","20-D1","20-D2",	
								"21-D1","21-D2","22-D1","22-D2","23-D1","23-D2","24-D1","24-D2","25-D1","25-D2")
CNA_freqpolyps = rep(0,50)
count=1
for (i in 1:length(polyp_names)){
	pat_id=which(adinca_ascat_segments$SampleID==polyp_names[i])
	segs = adinca_ascat_segments$nA[pat_id]+adinca_ascat_segments$nB[pat_id]
	medianploidy = median(segs)
	total_size=sum(as.numeric(adinca_ascat_segments$End[pat_id]-adinca_ascat_segments$Start[pat_id]))
	CNA_size=0 
	for ( j in 1: length(segs)){
		if (abs(segs[j]-medianploidy)>=1){
			CNA_size= CNA_size + (adinca_ascat_segments$End[pat_id[j]]-adinca_ascat_segments$Start[pat_id[j]])
		}
	}
	CNA_freqpolyps[count] = CNA_size/total_size
	count=count+1
}

CNA_freqTCGA= rep(0,length(TCGA_names))

count=1
for (i in TCGA_names){
	TCGA_data_call_stats[[i]]$CNt = TCGA_data_call_stats[[i]]$nMaj + TCGA_data_call_stats[[i]]$nMin
	segs=TCGA_data_call_stats[[i]]$CNt	
	medianploidy = median(segs)
	total_size=sum(as.numeric(TCGA_data_call_stats[[i]]$End-TCGA_data_call_stats[[i]]$Start))
	CNA_size=0 
	for ( j in 1: length(segs)){
		if (abs(segs[j]-medianploidy)>=1){
			CNA_size= CNA_size + (TCGA_data_call_stats[[i]]$End[j]-TCGA_data_call_stats[[i]]$Start[j])
		}
	}
	CNA_freqTCGA[count] = CNA_size/total_size
	count=count+1
}



## To do analysis by lesion, need lesion_id

## ID FOR LPWGS
lesion_id = paste(samples[,3],samples[,4],sep=".")
lesions=unique(lesion_id)
lesion_dx= NULL
CNA_freq_by_lesion= rep(0,length(patients))
for (i in 1:length(lesions)){
  lesion_ind = which(lesion_id==lesions[i])
  CNA_freq_by_lesion[i]=mean(CNA_freq[lesion_ind])
  lesion_dx = c(lesion_dx,as.character(tissue_type[lesion_ind[1]]))
}


## ID FOR EXOME
lesion_idExome = paste(exome_set,exome_dx,sep=".")
lesionsExome=unique(lesion_idExome)

lesion_dxExome= lesion_nameExome=rep("",length(lesionsExome))
CNA_freq_by_lesionE= rep(0,length(lesionsExome))
for (i in 1:length(lesionsExome)){
  lesion_ind = which(lesion_idExome==lesionsExome[i])
  CNA_freq_by_lesionE[i]=mean(CNA_freqExome[lesion_ind])
  lesion_dxExome[i] = as.character(exome_dx[lesion_ind[1]])
  lesion_nameExome[i]= as.character(exome_names[lesion_ind[1]])
}

### analysis by mean (or max) % of genome CNA across lesions
HGD_ind_lesion = which(lesion_dx=="H")
C_ind_lesion = which(lesion_dx=="C")
LGD_ind_lesion = which(lesion_dx=='L')

### Compare LGD/HGD with U2 samples, pathologists gave differing diagnoses of LGD or HGD, and confirmed U
U_ind_lesion = c(which(lesion_dx=="U2"), which(lesion_dx=="U"))

N_indE_lesion = which(lesion_dxExome=="N")
HGD_indE_lesion = which(lesion_dxExome=="H")
C_indE_lesion = which(lesion_dxExome=="C")
LGD_indE_lesion = which(substr(lesion_dxExome,1,1)=="L")

CNA_freqAll = c(CNA_freq_by_lesionE[N_indE_lesion], CNA_freq_by_lesion[LGD_ind_lesion],CNA_freq_by_lesionE[LGD_indE_lesion], CNA_freq_by_lesion[U_ind_lesion], CNA_freq_by_lesion[HGD_ind_lesion], CNA_freq_by_lesionE[HGD_indE_lesion], CNA_freq_by_lesion[C_ind_lesion], CNA_freq_by_lesionE[C_indE_lesion], CNA_freqSNP)
grades = factor(c(rep("N",9),rep("L",28),rep("U",7),rep("H",13),rep("C",25)),levels=c("N","L","U","H","C"))
DF <- data.frame(CNA= CNA_freqAll, grade = grades)

## Kruskal-Wallis test p-value
kwtest=kruskal.test(CNA~grade, data=DF)
KWp = kwtest$p.value
fit=lm(CNA~grade, data=DF)
tabletext<-c("LGD","Mixed","HGD","CA-CRC")


##################################### PLOT 1 #######################################################################################################
## N= 82 Lesions included: N (n=9) all from Exome, L (n=28) 24 from LPWGS and 4 from Exome, U (n=7) all from LPWGS, H (n=13) 11 from LPWGS, 2 from Exome, C (n=25) 11 from Exome, 9 from SNP, 5 from LPWGS
##
percent= 100
#dev.new()
setwd("/data/CA-CRC\ exome\ manuscript/Figure\ drafts/LPWGS_Figures")
pdf("CNAburden_vioplot_v2.pdf",width=9,height=6)
vioplot(CNA_freq_by_lesionE[N_indE_lesion]*percent,c(CNA_freq_by_lesionE[LGD_indE_lesion], CNA_freq_by_lesion[LGD_ind_lesion])*percent, CNA_freq_by_lesion[U_ind_lesion]*percent, c(CNA_freq_by_lesionE[HGD_indE_lesion],CNA_freq_by_lesion[HGD_ind_lesion])*percent,c(CNA_freq_by_lesion[C_ind_lesion],CNA_freq_by_lesionE[C_indE_lesion],CNA_freqSNP)*percent,names=c("Normal (n=9)","LGD (n=28)", "Mixed (n=7)","HGD (n=13)",  'CA-CRC (n=25)'), col=c("lightblue", "darkblue"),border=c("lightblue","darkblue"))
title(ylab="% of genome with CNA", cex.lab=.95)
beeswarm(percent*c(CNA_freq_by_lesionE[N_indE_lesion],c(CNA_freq_by_lesionE[LGD_indE_lesion], CNA_freq_by_lesion[LGD_ind_lesion]), CNA_freq_by_lesion[U_ind_lesion], c(CNA_freq_by_lesionE[HGD_indE_lesion],CNA_freq_by_lesion[HGD_ind_lesion]),c(CNA_freq_by_lesion[C_ind_lesion],CNA_freq_by_lesionE[C_indE_lesion],CNA_freqSNP))~grades, add=TRUE, pch=1,col="dimgrey")
text(1,80, sprintf("p =%11.2e", KWp))
dev.off()
##
####################################################################################################################################################


##################################### PLOT 2 #######################################################################################################
## N= 83 Lesions included: N (n=9) all from Exome, L (n=28) 24 from LPWGS and 4 from Exome, U (n=7) all from LPWGS, H (n=13) 11 from LPWGS, 2 from Exome, C (n=25) 11 from Exome, 9 from SNP, 5 from LPWGS
##
#dev.new()
pdf("CNAregression_forestplot.pdf",width=9,height=6)
forestplot(tabletext,fit$coefficients[2:5]*100,confint(fit, level=0.95)[2:5,1]*100,confint(fit, level=0.95)[2:5,2]*100,xlab= "% genome with CNA relative to Normal",clip =c(-.5, 60),col=fpColors(box= "darkred"),xticks = c(0, 10, 20, 30, 40 ,50, 60))
dev.off()
##
####################################################################################################################################################



## Analysis by lesion:
CNA_freqpolyps_by_lesion = rep(0,25)
count=1
for (i in 1:25){
	CNA_freqpolyps_by_lesion[i] = mean(c(CNA_freqpolyps[count],CNA_freqpolyps[count+1]))
	count = count+2 
	print(count)
}
Sgrades = factor(c(rep("A",25),rep("C",127)),levels=c("A","C"))

Bothgrades = factor(c(rep("N",9),rep("L",28),rep("U",7),rep("H",13),rep("C",25),rep("A",25),rep("SCRC",127)),levels=c("N","L","U","H","C","A","SCRC"))

##################################### PLOT for Sporadic CRC #######################################################################################################
## N= 83 Lesions included: N (n=9) all from Exome, L (n=28) 24 from LPWGS and 4 from Exome, U (n=7) all from LPWGS, H (n=13) 11 from LPWGS, 2 from Exome, C (n=25) 11 from Exome, 9 from SNP, 5 from LPWGS
##
#dev.new()
pdf("Fig3c_S.pdf",width=5,height=6)
vioplot(CNA_freqpolyps_by_lesion*percent,CNA_freqTCGA*percent,names=c("Adenomas (n=25)","S-CRC (n=127)"), col=c("lightblue", "darkblue"),border=c("lightblue","darkblue"))
title(ylab="% of genome with CNA", cex.lab=.95)
beeswarm(percent*c(CNA_freqpolyps_by_lesion,CNA_freqTCGA)~Sgrades, add=TRUE, pch=1,col="dimgrey")
text(1,80, sprintf("p =%11.2e", wilcox.test(CNA_freqpolyps_by_lesion,CNA_freqTCGA)$p.value))
dev.off()


percent= 100
#dev.new()
setwd("/data/CA-CRC\ exome\ manuscript/Figure\ drafts/LPWGS_Figures")
pdf("CNAburden_vioplot_v2_All.pdf",width=13,height=6)
vioplot(CNA_freq_by_lesionE[N_indE_lesion]*percent,c(CNA_freq_by_lesionE[LGD_indE_lesion], CNA_freq_by_lesion[LGD_ind_lesion])*percent, CNA_freq_by_lesion[U_ind_lesion]*percent, c(CNA_freq_by_lesionE[HGD_indE_lesion],CNA_freq_by_lesion[HGD_ind_lesion])*percent,c(CNA_freq_by_lesion[C_ind_lesion],CNA_freq_by_lesionE[C_indE_lesion],CNA_freqSNP)*percent,CNA_freqpolyps_by_lesion*percent,CNA_freqTCGA*percent,names=c("Normal (n=9)","LGD (n=28)", "Mixed (n=7)","HGD (n=13)",  'CA-CRC (n=25)',"Adenomas (n=25)","S-CRC (n=127)"), col=c("lightblue", "darkblue"),border=c("lightblue","darkblue"))
title(ylab="% of genome with CNA", cex.lab=.95)
beeswarm(percent*c(CNA_freq_by_lesionE[N_indE_lesion],c(CNA_freq_by_lesionE[LGD_indE_lesion], CNA_freq_by_lesion[LGD_ind_lesion]), CNA_freq_by_lesion[U_ind_lesion], c(CNA_freq_by_lesionE[HGD_indE_lesion],CNA_freq_by_lesion[HGD_ind_lesion]),c(CNA_freq_by_lesion[C_ind_lesion],CNA_freq_by_lesionE[C_indE_lesion]),CNA_freqSNP,CNA_freqpolyps_by_lesion,CNA_freqTCGA)~Bothgrades, add=TRUE, pch=1,col="dimgrey")
text(1,80, sprintf("p =%11.2e", KWp))
text(6,80, sprintf("p =%11.2e", wilcox.test(CNA_freqpolyps_by_lesion,CNA_freqTCGA)$p.value))
dev.off()
##
##
####################################################################################################################################################
vioplot_data = NULL
vioplot_data$Normal = CNA_freq_by_lesionE[N_indE_lesion]*percent
vioplot_data$LGD = c(CNA_freq_by_lesionE[LGD_indE_lesion], CNA_freq_by_lesion[LGD_ind_lesion])*percent
vioplot_data$U= CNA_freq_by_lesion[U_ind_lesion]*percent
vioplot_data$HGD = c(CNA_freq_by_lesionE[HGD_indE_lesion],CNA_freq_by_lesion[HGD_ind_lesion])*percent
vioplot_data$CACRC = c(CNA_freq_by_lesion[C_ind_lesion],CNA_freq_by_lesionE[C_indE_lesion],CNA_freqSNP)*percent
vioplot_data$polyp= CNA_freqpolyps_by_lesion*percent
vioplot_data$TCGA = CNA_freqTCGA*percent

write.xlsx(vioplot_data$Normal, file="violin_plot_v2.xlsx", sheetName="Normal", row.names=FALSE,col.names = FALSE)
write.xlsx(vioplot_data$LGD, file="violin_plot_v2.xlsx", sheetName="LGD", row.names=FALSE,append=TRUE,col.names=FALSE)
write.xlsx(vioplot_data$U, file="violin_plot_v2.xlsx", sheetName="U", row.names=FALSE,append=TRUE,col.names=FALSE)
write.xlsx(vioplot_data$HGD, file="violin_plot_v2.xlsx", sheetName="HGD", row.names=FALSE,append=TRUE,col.names=FALSE)
write.xlsx(vioplot_data$CACRC, file="violin_plot_v2.xlsx", sheetName="CACRC", row.names=FALSE,append=TRUE,col.names=FALSE)
write.xlsx(vioplot_data$polyp, file="violin_plot_v2.xlsx", sheetName="Polyp", row.names=FALSE,append=TRUE,col.names=FALSE)
write.xlsx(vioplot_data$TCGA, file="violin_plot_v2.xlsx", sheetName="SCRC", row.names=FALSE,append=TRUE,col.names=FALSE)






call_stats_N = list()
for (i in filenames[N_ind]){
	call_stats_N[[i]]=call_stats[[i]]
}
call_stats_LGD = list()
for (i in filenames[LGD_ind]){
	call_stats_LGD[[i]]=call_stats[[i]]
}
call_stats_HGD = list()
for (i in filenames[HGD_ind]){
	call_stats_HGD[[i]]=call_stats[[i]]
}
call_stats_U = list()
for (i in filenames[U_ind]){
	call_stats_U[[i]]=call_stats[[i]]
}
call_stats_C = list()
for (i in filenames[C_ind]){
	call_stats_C[[i]]=call_stats[[i]]
}

chr_start_spot=c(0,22)
for (i in 1:22){
	chr_start_spot[i]=which(chroms_4401==i)[1]

}




#### Version with Exome and SNP samples (cnLOH events not included as change)
### N (n=9) all from Exome
run_sizeN = length(N_indE)
all_gain_counts = NULL
for (i in 1:total_binslowpass){
  	## Exome N data
	chrom = bins_4401[i,1]
  	start = bins_4401[i,2]
  	end = bins_4401[i,3]
  	for (j in N_indE){
  		exome_chroms = substr(Exome_data_call_stats[[j]]$chromosome,4,length(Exome_data_call_stats[[j]]$chromosome))
  		medianploidy = median(Exome_data_call_stats[[j]]$CNt)
  		if (sum(exome_chroms==chrom)>0){
			chromseg_ind = which(exome_chroms==chrom)
			for (k in chromseg_ind){
				if (Exome_data_call_stats[[j]]$start.pos[k]<start && Exome_data_call_stats[[j]]$end.pos[k]> end ){
					#if (Exome_data_call_stats[[j]]$CNt[k]>2){
					if ((Exome_data_call_stats[[j]]$CNt[k]-medianploidy)>=1){
						all_gain_counts=c(all_gain_counts,i)
						break
					}
				}
				else if(Exome_data_call_stats[[j]]$start.pos[k]>start && Exome_data_call_stats[[j]]$start.pos[k]< end ){
					#if (Exome_data_call_stats[[j]]$CNt[k]>2){
					if ((Exome_data_call_stats[[j]]$CNt[k]-medianploidy)>=1){
						all_gain_counts=c(all_gain_counts,i)
						break
					}
				}
				else if(Exome_data_call_stats[[j]]$end.pos[k]>start && Exome_data_call_stats[[j]]$end.pos[k]< end ){
					#if (Exome_data_call_stats[[j]]$CNt[k]>2){
					if ((Exome_data_call_stats[[j]]$CNt[k]-medianploidy)>=1){
						all_gain_counts=c(all_gain_counts,i)
						break
					}
				}
			}
  		}
  	}
}

gainfreq_each_bin = rep(0,total_binslowpass)
for (i in 1:total_binslowpass){
  gainfreq_each_bin[i] = sum(all_gain_counts==i)
}

gainfreqN = gainfreq_each_bin/run_sizeN*100


all_loss_counts = NULL

for (i in 1:total_binslowpass){
 	chrom = bins_4401[i,1]
  	start = bins_4401[i,2]
  	end = bins_4401[i,3]
  	for (j in N_indE){
  		exome_chroms = substr(Exome_data_call_stats[[j]]$chromosome,4,length(Exome_data_call_stats[[j]]$chromosome))
  		medianploidy = median(Exome_data_call_stats[[j]]$CNt)
  		if (sum(exome_chroms==chrom)>0){
			chromseg_ind = which(exome_chroms==chrom)
			for (k in chromseg_ind){
				if (Exome_data_call_stats[[j]]$start.pos[k]<start && Exome_data_call_stats[[j]]$end.pos[k]> end ){
					#if (Exome_data_call_stats[[j]]$CNt[k]<2){
					if ((medianploidy-Exome_data_call_stats[[j]]$CNt[k])>=1){
						all_loss_counts=c(all_loss_counts,i)
						break
					}
				}
				else if(Exome_data_call_stats[[j]]$start.pos[k]>start && Exome_data_call_stats[[j]]$start.pos[k]< end ){
					#if (Exome_data_call_stats[[j]]$CNt[k]<2){
					if ((medianploidy-Exome_data_call_stats[[j]]$CNt[k])>=1){
						all_loss_counts=c(all_loss_counts,i)
						break
					}
				}
				else if(Exome_data_call_stats[[j]]$end.pos[k]>start && Exome_data_call_stats[[j]]$end.pos[k]< end ){
					#if (Exome_data_call_stats[[j]]$CNt[k]<2){
					if ((medianploidy-Exome_data_call_stats[[j]]$CNt[k])>=1){
						all_loss_counts=c(all_loss_counts,i)
						break
					}
				}
			}
  		}
  	}
}
lossfreq_each_bin = rep(0,total_binslowpass)
for (i in 1:total_binslowpass){
  lossfreq_each_bin[i] = sum(all_loss_counts==i)
}

lossfreqN = lossfreq_each_bin/run_sizeN*100


### LGD: (n=28) 24 from LPWGS and 4 from Exome,

call_stats=call_stats_LGD
run_sizeL = length(LGD_ind)+length(LGD_indE)

all_gain_counts = NULL
for (i in 1:total_binslowpass){
	## lowpass data 
  	all_gain_counts = c(all_gain_counts,rep(i,sum(lapply(call_stats,'[','calls',i)>=1)))
  	## Exome LGD data
  	chrom = bins_4401[i,1]
  	start = bins_4401[i,2]
  	end = bins_4401[i,3]
  	for (j in LGD_indE){
  		exome_chroms = substr(Exome_data_call_stats[[j]]$chromosome,4,length(Exome_data_call_stats[[j]]$chromosome))
  		medianploidy=median(Exome_data_call_stats[[j]]$CNt)
  		if (sum(exome_chroms==chrom)>0){
			chromseg_ind = which(exome_chroms==chrom)
			for (k in chromseg_ind){
				if (Exome_data_call_stats[[j]]$start.pos[k]<start && Exome_data_call_stats[[j]]$end.pos[k]> end ){
					#if (Exome_data_call_stats[[j]]$CNt[k]>2){
					if ((Exome_data_call_stats[[j]]$CNt[k]-medianploidy)>=1){
						all_gain_counts=c(all_gain_counts,i)
						break
					}
				}
				else if(Exome_data_call_stats[[j]]$start.pos[k]>start && Exome_data_call_stats[[j]]$start.pos[k]< end ){
					#if (Exome_data_call_stats[[j]]$CNt[k]>2){
					if ((Exome_data_call_stats[[j]]$CNt[k]-medianploidy)>=1){
						all_gain_counts=c(all_gain_counts,i)
						break
					}
				}
				else if(Exome_data_call_stats[[j]]$end.pos[k]>start && Exome_data_call_stats[[j]]$end.pos[k]< end ){
					#if (Exome_data_call_stats[[j]]$CNt[k]>2){
					if ((Exome_data_call_stats[[j]]$CNt[k]-medianploidy)>=1){
						all_gain_counts=c(all_gain_counts,i)
						break
					}
				}
			}
  		}
  	}
}

gainfreq_each_bin = rep(0,total_binslowpass)
for (i in 1:total_binslowpass){
  gainfreq_each_bin[i] = sum(all_gain_counts==i)
}

gainfreqL = gainfreq_each_bin/run_sizeL*100


all_loss_counts = NULL

for (i in 1:total_binslowpass){
  	all_loss_counts = c(all_loss_counts,rep(i,sum(lapply(call_stats,'[','calls',i)<=-1)))
  ## Exome LGD data
  	chrom = bins_4401[i,1]
  	start = bins_4401[i,2]
  	end = bins_4401[i,3]
  	for (j in LGD_indE){
  		exome_chroms = substr(Exome_data_call_stats[[j]]$chromosome,4,length(Exome_data_call_stats[[j]]$chromosome))
  		medianploidy= median(Exome_data_call_stats[[j]]$CNt)
  		if (sum(exome_chroms==chrom)>0){
			chromseg_ind = which(exome_chroms==chrom)
			for (k in chromseg_ind){
				if (Exome_data_call_stats[[j]]$start.pos[k]<start && Exome_data_call_stats[[j]]$end.pos[k]> end ){
					#if (Exome_data_call_stats[[j]]$CNt[k]<2){
					if ((medianploidy-Exome_data_call_stats[[j]]$CNt[k])>=1){
						all_loss_counts=c(all_loss_counts,i)
						break
					}
				}
				else if(Exome_data_call_stats[[j]]$start.pos[k]>start && Exome_data_call_stats[[j]]$start.pos[k]< end ){
					#if (Exome_data_call_stats[[j]]$CNt[k]<2){
					if ((medianploidy-Exome_data_call_stats[[j]]$CNt[k])>=1){
						all_loss_counts=c(all_loss_counts,i)
						break
					}
				}
				else if(Exome_data_call_stats[[j]]$end.pos[k]>start && Exome_data_call_stats[[j]]$end.pos[k]< end ){
					#if (Exome_data_call_stats[[j]]$CNt[k]<2){
					if ((medianploidy-Exome_data_call_stats[[j]]$CNt[k])>=1){
						all_loss_counts=c(all_loss_counts,i)
						break
					}
				}
			}
  		}
  	}
}
lossfreq_each_bin = rep(0,total_binslowpass)
for (i in 1:total_binslowpass){
  lossfreq_each_bin[i] = sum(all_loss_counts==i)
}

lossfreqL = lossfreq_each_bin/run_sizeL*100




### U (n=7) all from LPWGS

call_stats=call_stats_U
run_sizeU = length(U_ind)

all_gain_counts = NULL
for (i in 1:total_binslowpass){
    all_gain_counts = c(all_gain_counts,rep(i,sum(lapply(call_stats,'[','calls',i)>=1))) 
}
gainfreq_each_bin = rep(0,total_binslowpass)
for (i in 1:total_binslowpass){
    gainfreq_each_bin[i] = sum(all_gain_counts==i)
}

gainfreqU = gainfreq_each_bin/run_sizeU*100


all_loss_counts = NULL

for (i in 1:total_binslowpass){
    all_loss_counts = c(all_loss_counts,rep(i,sum(lapply(call_stats,'[','calls',i)<=-1)))
}
lossfreq_each_bin = rep(0,total_binslowpass)
for (i in 1:total_binslowpass){
    lossfreq_each_bin[i] = sum(all_loss_counts==i)
}

lossfreqU = lossfreq_each_bin/run_sizeU*100






### H (n=13) 11 from LPWGS, 2 from Exome
call_stats=call_stats_HGD
run_sizeH = length(HGD_ind)+length(HGD_indE)

all_gain_counts = NULL
for (i in 1:total_binslowpass){
	## lowpass data 
  	all_gain_counts = c(all_gain_counts,rep(i,sum(lapply(call_stats,'[','calls',i)>=1)))
  	## Exome HGD data
  	chrom = bins_4401[i,1]
  	start = bins_4401[i,2]
  	end = bins_4401[i,3]
  	for (j in HGD_indE){
  		exome_chroms = substr(Exome_data_call_stats[[j]]$chromosome,4,length(Exome_data_call_stats[[j]]$chromosome))
  		medianploidy=median(Exome_data_call_stats[[j]]$CNt)
  		if (sum(exome_chroms==chrom)>0){
			chromseg_ind = which(exome_chroms==chrom)
			for (k in chromseg_ind){
				if (Exome_data_call_stats[[j]]$start.pos[k]<start && Exome_data_call_stats[[j]]$end.pos[k]> end ){
					#if (Exome_data_call_stats[[j]]$CNt[k]>2){
					if ((Exome_data_call_stats[[j]]$CNt[k]-medianploidy)>=1){
						all_gain_counts=c(all_gain_counts,i)
						break
					}
				}
				else if(Exome_data_call_stats[[j]]$start.pos[k]>start && Exome_data_call_stats[[j]]$start.pos[k]< end ){
					#if (Exome_data_call_stats[[j]]$CNt[k]>2){
					if ((Exome_data_call_stats[[j]]$CNt[k]-medianploidy)>=1){
						all_gain_counts=c(all_gain_counts,i)
						break
					}
				}
				else if(Exome_data_call_stats[[j]]$end.pos[k]>start && Exome_data_call_stats[[j]]$end.pos[k]< end ){
					#if (Exome_data_call_stats[[j]]$CNt[k]>2){
					if ((Exome_data_call_stats[[j]]$CNt[k]-medianploidy)>=1){
						all_gain_counts=c(all_gain_counts,i)
						break
					}
				}
			}
  		}
  	}
}

gainfreq_each_bin = rep(0,total_binslowpass)
for (i in 1:total_binslowpass){
  gainfreq_each_bin[i] = sum(all_gain_counts==i)
}

gainfreqH = gainfreq_each_bin/run_sizeH*100


all_loss_counts = NULL

for (i in 1:total_binslowpass){
  all_loss_counts = c(all_loss_counts,rep(i,sum(lapply(call_stats,'[','calls',i)<=-1)))
  ## Exome HGD data
  	chrom = bins_4401[i,1]
  	start = bins_4401[i,2]
  	end = bins_4401[i,3]
  	for (j in HGD_indE){
  		exome_chroms = substr(Exome_data_call_stats[[j]]$chromosome,4,length(Exome_data_call_stats[[j]]$chromosome))
  		medianploidy =median(Exome_data_call_stats[[j]]$CNt)
  		if (sum(exome_chroms==chrom)>0){
			chromseg_ind = which(exome_chroms==chrom)
			for (k in chromseg_ind){
				if (Exome_data_call_stats[[j]]$start.pos[k]<start && Exome_data_call_stats[[j]]$end.pos[k]> end ){
					#if (Exome_data_call_stats[[j]]$CNt[k]<2){
					if ((medianploidy-Exome_data_call_stats[[j]]$CNt[k])>=1){
						all_loss_counts=c(all_loss_counts,i)
						break
					}
				}
				else if(Exome_data_call_stats[[j]]$start.pos[k]>start && Exome_data_call_stats[[j]]$start.pos[k]< end ){
					#if (Exome_data_call_stats[[j]]$CNt[k]<2){
					if ((medianploidy-Exome_data_call_stats[[j]]$CNt[k])>=1){
						all_loss_counts=c(all_loss_counts,i)
						break
					}
				}
				else if(Exome_data_call_stats[[j]]$end.pos[k]>start && Exome_data_call_stats[[j]]$end.pos[k]< end ){
					#if (Exome_data_call_stats[[j]]$CNt[k]<2){
					if ((medianploidy-Exome_data_call_stats[[j]]$CNt[k])>=1){
						all_loss_counts=c(all_loss_counts,i)
						break
					}
				}
			}
  		}
  	}
}
lossfreq_each_bin = rep(0,total_binslowpass)
for (i in 1:total_binslowpass){
  lossfreq_each_bin[i] = sum(all_loss_counts==i)
}

lossfreqH = lossfreq_each_bin/run_sizeH*100



### C(n=25) 11 from Exome, 9 from SNP, 5 from LPWGS
call_stats=call_stats_C
run_sizeC = length(C_ind)+length(C_indE) +length(snp_names)

all_gain_counts = NULL
for (i in 1:total_binslowpass){
	## lowpass data 
  	all_gain_counts = c(all_gain_counts,rep(i,sum(lapply(call_stats,'[','calls',i)>=1)))
  	## Exome and SNP CA-CRC data
  	chrom = bins_4401[i,1]
  	start = bins_4401[i,2]
  	end = bins_4401[i,3]
  	for (j in C_indE){
  		exome_chroms = substr(Exome_data_call_stats[[j]]$chromosome,4,length(Exome_data_call_stats[[j]]$chromosome))
  		medianploidy = median(Exome_data_call_stats[[j]]$CNt)
  		if (sum(exome_chroms==chrom)>0){
			chromseg_ind = which(exome_chroms==chrom)
			for (k in chromseg_ind){
				if (Exome_data_call_stats[[j]]$start.pos[k]<start && Exome_data_call_stats[[j]]$end.pos[k]> end ){
					#if (Exome_data_call_stats[[j]]$CNt[k]>2){
					if ((Exome_data_call_stats[[j]]$CNt[k]-medianploidy)>=1){
						all_gain_counts=c(all_gain_counts,i)
						break
					}
				}
				else if(Exome_data_call_stats[[j]]$start.pos[k]>start && Exome_data_call_stats[[j]]$start.pos[k]< end ){
					#if (Exome_data_call_stats[[j]]$CNt[k]>2){
					if ((Exome_data_call_stats[[j]]$CNt[k]-medianploidy)>=1){
						all_gain_counts=c(all_gain_counts,i)
						break
					}
				}
				else if(Exome_data_call_stats[[j]]$end.pos[k]>start && Exome_data_call_stats[[j]]$end.pos[k]< end ){
					#if (Exome_data_call_stats[[j]]$CNt[k]>2){
					if ((Exome_data_call_stats[[j]]$CNt[k]-medianploidy)>=1){
						all_gain_counts=c(all_gain_counts,i)
						break
					}
				}
			}
  		}
  	}
  	for (j in 1:length(snp_names)){
  		pat_id=which(ryan_ascat_segments$SampleID==snp_names[j])
  		segs = ryan_ascat_segments$nA[pat_id]+ryan_ascat_segments$nB[pat_id]
  		medianploidy=median(segs)
  		SNP_chroms = ryan_ascat_segments$Chr[pat_id]
  		if (sum(SNP_chroms==chrom)>0){
			chromseg_ind = which(SNP_chroms==chrom)
			for (k in chromseg_ind){
				if (ryan_ascat_segments$Start[pat_id][k]<start && ryan_ascat_segments$End[pat_id][k]> end ){
					#if (segs[k]>2){
					if ((segs[k]-medianploidy)>=1){
						all_gain_counts=c(all_gain_counts,i)
						break
					}
				}
				else if(ryan_ascat_segments$Start[pat_id][k]>start && ryan_ascat_segments$Start[pat_id][k]< end ){
					#if (segs[k]>2){
					if ((segs[k]-medianploidy)>=1){
						all_gain_counts=c(all_gain_counts,i)
						break
					}
				}
				else if(ryan_ascat_segments$End[pat_id][k]>start && ryan_ascat_segments$End[pat_id][k]<end ){
					#if (segs[k]>2){
					if ((segs[k]-medianploidy)>=1){
						all_gain_counts=c(all_gain_counts,i)
						break
					}
				}	
			}
  		}
  	}
}

gainfreq_each_bin = rep(0,total_binslowpass)
for (i in 1:total_binslowpass){
  gainfreq_each_bin[i] = sum(all_gain_counts==i)
}

gainfreqC = gainfreq_each_bin/run_sizeC*100


all_loss_counts = NULL

for (i in 1:total_binslowpass){
  all_loss_counts = c(all_loss_counts,rep(i,sum(lapply(call_stats,'[','calls',i)<=-1)))
  ## Exome + SNP CA-CRC data
  	chrom = bins_4401[i,1]
  	start = bins_4401[i,2]
  	end = bins_4401[i,3]
  	for (j in C_indE){
  		exome_chroms = substr(Exome_data_call_stats[[j]]$chromosome,4,length(Exome_data_call_stats[[j]]$chromosome))
  		medianploidy = median(Exome_data_call_stats[[j]]$CNt)
  		if (sum(exome_chroms==chrom)>0){
			chromseg_ind = which(exome_chroms==chrom)
			for (k in chromseg_ind){
				if (Exome_data_call_stats[[j]]$start.pos[k]<start && Exome_data_call_stats[[j]]$end.pos[k]> end ){
					#if (Exome_data_call_stats[[j]]$CNt[k]<2){
					if ((medianploidy-Exome_data_call_stats[[j]]$CNt[k])>=1){
						all_loss_counts=c(all_loss_counts,i)
						break
					}
					#else if(Exome_data_call_stats[[j]]$A[k]==0 || Exome_data_call_stats[[j]]$B[k]==0){
					#	if (Exome_data_call_stats[[j]]$A[k]>2 || Exome_data_call_stats[[j]]$B[k]>2){
					#		all_loss_counts=c(all_loss_counts,i)
					#		break
					#	}
					#}
				}
				else if(Exome_data_call_stats[[j]]$start.pos[k]>start && Exome_data_call_stats[[j]]$start.pos[k]< end ){
					#if (Exome_data_call_stats[[j]]$CNt[k]<2){
					if ((medianploidy-Exome_data_call_stats[[j]]$CNt[k])>=1){
						all_loss_counts=c(all_loss_counts,i)
						break
					}
					#else if(Exome_data_call_stats[[j]]$A[k]==0 || Exome_data_call_stats[[j]]$B[k]==0){
					#	if (Exome_data_call_stats[[j]]$A[k]>2 || Exome_data_call_stats[[j]]$B[k]>2){
					#		all_loss_counts=c(all_loss_counts,i)
					#		break
					#	}
					#}
				}
				else if(Exome_data_call_stats[[j]]$end.pos[k]>start && Exome_data_call_stats[[j]]$end.pos[k]< end ){
					#if (Exome_data_call_stats[[j]]$CNt[k]<2){
					if ((medianploidy-Exome_data_call_stats[[j]]$CNt[k])>=1){
						all_loss_counts=c(all_loss_counts,i)
						break
					}
					#else if(Exome_data_call_stats[[j]]$A[k]==0 || Exome_data_call_stats[[j]]$B[k]==0){
					#	if (Exome_data_call_stats[[j]]$A[k]>2 || Exome_data_call_stats[[j]]$B[k]>2){
					#		all_loss_counts=c(all_loss_counts,i)
					#		break
					#	}
					#}
				}
			}
  		}
  	}
  	for (j in 1:length(snp_names)){
  		pat_id=which(ryan_ascat_segments$SampleID==snp_names[j])
  		segs = ryan_ascat_segments$nA[pat_id]+ryan_ascat_segments$nB[pat_id]
  		medianploidy= median(segs)
  		SNP_chroms = ryan_ascat_segments$Chr[pat_id]
  		if (sum(SNP_chroms==chrom)>0){
			chromseg_ind = which(SNP_chroms==chrom)
			for (k in chromseg_ind){
				if (ryan_ascat_segments$Start[pat_id][k]<start && ryan_ascat_segments$End[pat_id][k]> end ){
					#if (segs[k]<2){
					if ((medianploidy-segs[k])>=1){
						all_loss_counts=c(all_loss_counts,i)
						break
					}
					#else if(ryan_ascat_segments$nA[pat_id][k]==0 || ryan_ascat_segments$nB[pat_id][k]==0){
					#	if (ryan_ascat_segments$nA[pat_id][k]>2 || ryan_ascat_segments$nB[pat_id][k]>2){
					#		all_loss_counts=c(all_loss_counts,i)
					#		break
					#	}
					#}
				}
				else if(ryan_ascat_segments$Start[pat_id][k]>start && ryan_ascat_segments$Start[pat_id][k]< end ){
					#if (segs[k]<2)
					if ((medianploidy-segs[k])>=1){
						all_loss_counts=c(all_loss_counts,i)
						break
					}
					#else if(ryan_ascat_segments$nA[pat_id][k]==0 || ryan_ascat_segments$nB[pat_id][k]==0){
					#	if (ryan_ascat_segments$nA[pat_id][k]>2 || ryan_ascat_segments$nB[pat_id][k]>2){
					#		all_loss_counts=c(all_loss_counts,i)
					#		break
					#	}
					#}
				}
				else if(ryan_ascat_segments$End[pat_id][k]>start && ryan_ascat_segments$End[pat_id][k]<end ){
					#if (segs[k]<2){
					if ((medianploidy-segs[k])>=1){
						all_loss_counts=c(all_loss_counts,i)
						break
					}
					#else if(ryan_ascat_segments$nA[pat_id][k]==0 || ryan_ascat_segments$nB[pat_id][k]==0){
					#	if (ryan_ascat_segments$nA[pat_id][k]>2 || ryan_ascat_segments$nB[pat_id][k]>2){
					#		all_loss_counts=c(all_loss_counts,i)
					#		break
					#	}
					#}
				}
			}
  		}
  	}
}
lossfreq_each_bin = rep(0,total_binslowpass)
for (i in 1:total_binslowpass){
  lossfreq_each_bin[i] = sum(all_loss_counts==i)
}

lossfreqC = lossfreq_each_bin/run_sizeC*100

### Sporadic Adenomas (n= 25, 2 samples from each for total of 50 samples )

run_sizepolyps = length(polyp_names)

all_gain_counts = NULL
for (i in 1:total_binslowpass){
	chrom = bins_4401[i,1]
  	start = bins_4401[i,2]
  	end = bins_4401[i,3]
  	for (j in 1:length(polyp_names)){
  		pat_id=which(adinca_ascat_segments$SampleID==polyp_names[j])
  		segs = adinca_ascat_segments$nA[pat_id]+adinca_ascat_segments$nB[pat_id]
  		medianploidy=median(segs)
  		SNP_chroms = adinca_ascat_segments$Chr[pat_id]
  		if (sum(SNP_chroms==chrom)>0){
			chromseg_ind = which(SNP_chroms==chrom)
			for (k in chromseg_ind){
				if (adinca_ascat_segments$Start[pat_id][k]<start && adinca_ascat_segments$End[pat_id][k]> end ){
					#if (segs[k]>2){
					if ((segs[k]-medianploidy)>=1){
						all_gain_counts=c(all_gain_counts,i)
						break
					}
				}
				else if(adinca_ascat_segments$Start[pat_id][k]>start && adinca_ascat_segments$Start[pat_id][k]< end ){
					#if (segs[k]>2){
					if ((segs[k]-medianploidy)>=1){
						all_gain_counts=c(all_gain_counts,i)
						break
					}
				}
				else if(adinca_ascat_segments$End[pat_id][k]>start && adinca_ascat_segments$End[pat_id][k]<end ){
					#if (segs[k]>2){
					if ((segs[k]-medianploidy)>=1){
						all_gain_counts=c(all_gain_counts,i)
						break
					}
				}	
			}
  		}
  	}
}

gainfreq_each_bin = rep(0,total_binslowpass)
for (i in 1:total_binslowpass){
  gainfreq_each_bin[i] = sum(all_gain_counts==i)
}

gainfreqpolyp = gainfreq_each_bin/run_sizepolyps*100


all_loss_counts = NULL

for (i in 1:total_binslowpass){
  	chrom = bins_4401[i,1]
  	start = bins_4401[i,2]
  	end = bins_4401[i,3]
  	for (j in 1:length(polyp_names)){
  		pat_id=which(adinca_ascat_segments$SampleID==polyp_names[j])
  		segs = adinca_ascat_segments$nA[pat_id]+adinca_ascat_segments$nB[pat_id]
  		medianploidy= median(segs)
  		SNP_chroms = adinca_ascat_segments$Chr[pat_id]
  		if (sum(SNP_chroms==chrom)>0){
			chromseg_ind = which(SNP_chroms==chrom)
			for (k in chromseg_ind){
				if (adinca_ascat_segments$Start[pat_id][k]<start && adinca_ascat_segments$End[pat_id][k]> end ){
					#if (segs[k]<2){
					if ((medianploidy-segs[k])>=1){
						all_loss_counts=c(all_loss_counts,i)
						break
					}
				}
				else if(adinca_ascat_segments$Start[pat_id][k]>start && adinca_ascat_segments$Start[pat_id][k]< end ){
					#if (segs[k]<2)
					if ((medianploidy-segs[k])>=1){
						all_loss_counts=c(all_loss_counts,i)
						break
					}
				}
				else if(adinca_ascat_segments$End[pat_id][k]>start && adinca_ascat_segments$End[pat_id][k]<end ){
					#if (segs[k]<2){
					if ((medianploidy-segs[k])>=1){
						all_loss_counts=c(all_loss_counts,i)
						break
					}
				}
			}
  		}
  	}
}
lossfreq_each_bin = rep(0,total_binslowpass)
for (i in 1:total_binslowpass){
  lossfreq_each_bin[i] = sum(all_loss_counts==i)
}

lossfreqpolyp= lossfreq_each_bin/run_sizepolyps*100


### TCGA Sporadic CRC (n = 127)
TCGA_names = as.character(read.table('/data/LGD_HGD_CACRC/CRCids.txt',header=T)[,1])
missing_IDs = c(14,17,18,28,64,84,90,91,93,111,117,126)
TCGA_names = TCGA_names[-missing_IDs]
run_sizeTCGA = length(TCGA_names)


TCGA_data_call_stats = NULL

for (i in 1:length(TCGA_names)){
	sample = TCGA_names[i]
	TCGA_data_call_stats[[i]] = read.table(paste('/data/LGD_HGD_CACRC/TCGA-CN/',sample,'.txt',sep=""),head=TRUE)	
}
names(TCGA_data_call_stats) = TCGA_names

all_gain_counts = NULL
for (i in 1:total_binslowpass){
	chrom = bins_4401[i,1]
  	start = bins_4401[i,2]
  	end = bins_4401[i,3]
  	for (j in TCGA_names){
  		exome_chroms = TCGA_data_call_stats[[j]]$Chr
  		TCGA_data_call_stats[[j]]$CNt = TCGA_data_call_stats[[j]]$nMaj + TCGA_data_call_stats[[j]]$nMin
  		medianploidy = median(TCGA_data_call_stats[[j]]$CNt)
  		if (sum(exome_chroms==chrom)>0){
			chromseg_ind = which(exome_chroms==chrom)
			for (k in chromseg_ind){
				if (TCGA_data_call_stats[[j]]$Start[k]<start && TCGA_data_call_stats[[j]]$End[k]> end ){
					if ((TCGA_data_call_stats[[j]]$CNt[k]-medianploidy)>=1){
						all_gain_counts=c(all_gain_counts,i)
						break
					}
				}
				else if(TCGA_data_call_stats[[j]]$Start[k]>start && TCGA_data_call_stats[[j]]$Start[k]< end ){
					#if (Exome_data_call_stats[[j]]$CNt[k]>2){
					if ((TCGA_data_call_stats[[j]]$CNt[k]-medianploidy)>=1){
						all_gain_counts=c(all_gain_counts,i)
						break
					}
				}
				else if(TCGA_data_call_stats[[j]]$End[k]>start && TCGA_data_call_stats[[j]]$End[k]< end ){
					#if (Exome_data_call_stats[[j]]$CNt[k]>2){
					if ((TCGA_data_call_stats[[j]]$CNt[k]-medianploidy)>=1){
						all_gain_counts=c(all_gain_counts,i)
						break
					}
				}
			}
  		}
  	}
}

gainfreq_each_bin = rep(0,total_binslowpass)
for (i in 1:total_binslowpass){
  gainfreq_each_bin[i] = sum(all_gain_counts==i)
}

gainfreqTCGA = gainfreq_each_bin/run_sizeTCGA*100


all_loss_counts = NULL

for (i in 1:total_binslowpass){
  	chrom = bins_4401[i,1]
  	start = bins_4401[i,2]
  	end = bins_4401[i,3]
  	for (j in TCGA_names){
  		exome_chroms = TCGA_data_call_stats[[j]]$Chr
  		TCGA_data_call_stats[[j]]$CNt = TCGA_data_call_stats[[j]]$nMaj + TCGA_data_call_stats[[j]]$nMin
  		medianploidy = median(TCGA_data_call_stats[[j]]$CNt)
  		if (sum(exome_chroms==chrom)>0){
			chromseg_ind = which(exome_chroms==chrom)
			for (k in chromseg_ind){
				if (TCGA_data_call_stats[[j]]$Start[k]<start && TCGA_data_call_stats[[j]]$End[k]> end ){
					if ((medianploidy-TCGA_data_call_stats[[j]]$CNt[k])>=1){
						all_loss_counts=c(all_loss_counts,i)
						break
					}
				}
				else if(TCGA_data_call_stats[[j]]$Start[k]>start && TCGA_data_call_stats[[j]]$Start[k]< end ){
					if ((medianploidy-TCGA_data_call_stats[[j]]$CNt[k])>=1){
						all_loss_counts=c(all_loss_counts,i)
						break
					}
				}
				else if(TCGA_data_call_stats[[j]]$End[k]>start && TCGA_data_call_stats[[j]]$End[k]< end ){
					if ((medianploidy-TCGA_data_call_stats[[j]]$CNt[k])>=1){
						all_loss_counts=c(all_loss_counts,i)
						break
					}
				}
			}
		}	
  	}
}
lossfreq_each_bin = rep(0,total_binslowpass)
for (i in 1:total_binslowpass){
  lossfreq_each_bin[i] = sum(all_loss_counts==i)
}

lossfreqTCGA= lossfreq_each_bin/run_sizeTCGA*100

dev.new()
plot(c(0,total_binslowpass),c(-100,100),type = "n", ylab = "Percentage of TCGA samples", xlab="Chromosome", axes=F)
legend('topleft',c("Gains","Losses"),col=c("red","blue"), lty=1,lwd=2)
axis(2,at=seq(-100,100,20),labels=c(100,80, 60, 40,20,0,20,40,60,80,100))
axis(1,at=chr_start_spot, labels=seq(1,22))
lines( gainfreqTCGA,col="red",lwd=2)
lines(-lossfreqTCGA,col="blue",lwd=2)
abline(0,0)
box()

##################################### PLOT 3 #######################################################################################################

dev.new()
par(mfrow=c(5,1))
plot(c(0,total_binslowpass),c(-100,100),type = "n", ylab = "Percentage of normal samples", xlab="Chromosome", axes=F)
legend('topleft',c("Gains","Losses"),col=c("red","blue"), lty=1,lwd=2)
axis(2,at=seq(-100,100,20),labels=c(100,80, 60, 40,20,0,20,40,60,80,100))
axis(1,at=chr_start_spot, labels=seq(1,22))
lines( gainfreqN,col="red",lwd=2)
lines(-lossfreqN,col="blue",lwd=2)
abline(0,0)
box()

plot(c(0,total_binslowpass),c(-100,100),type = "n", ylab = "Percentage of LGD samples", xlab="Chromosome", axes=F)
legend('topleft',c("Gains","Losses"),col=c("red","blue"), lty=1,lwd=2)
axis(2,at=seq(-100,100,20),labels=c(100,80, 60, 40,20,0,20,40,60,80,100))
axis(1,at=chr_start_spot, labels=seq(1,22))
lines( gainfreqL,col="red",lwd=2)
lines(-lossfreqL,col="blue",lwd=2)
abline(0,0)
box()



plot(c(0,total_binslowpass),c(-100,100),type = "n", ylab = "Percentage of U samples", xlab="Chromosome", axes=F)
legend('topleft',c("Gains","Losses"),col=c("red","blue"), lty=1,lwd=2)
axis(2,at=seq(-100,100,20),labels=c(100,80, 60, 40,20,0,20,40,60,80,100))
axis(1,at=chr_start_spot, labels=seq(1,22))
lines( gainfreqU,col="red",lwd=2)
lines(-lossfreqU,col="blue",lwd=2)
abline(0,0)
box()


plot(c(0,total_binslowpass),c(-100,100),type = "n", ylab = "Percentage of HGD samples", xlab="Chromosome", axes=F)
legend('topleft',c("Gains","Losses"),col=c("red","blue"), lty=1,lwd=2)
axis(2,at=seq(-100,100,20),labels=c(100,80, 60, 40,20,0,20,40,60,80,100))
axis(1,at=chr_start_spot, labels=seq(1,22))
lines( gainfreqH,col="red",lwd=2)
lines(-lossfreqH,col="blue",lwd=2)
abline(0,0)
box()

plot(c(0,total_binslowpass),c(-100,100),type = "n", ylab = "Percentage of CA-CRC samples", xlab="Chromosome", axes=F)
legend('topleft',c("Gains","Losses"),col=c("red","blue"), lty=1,lwd=2)
axis(2,at=seq(-100,100,20),labels=c(100,80, 60, 40,20,0,20,40,60,80,100))
axis(1,at=chr_start_spot, labels=seq(1,22))
lines( gainfreqC,col="red",lwd=2)
lines(-lossfreqC,col="blue",lwd=2)
abline(0,0)
box()

par(mfrow=c(5,1))
plot(c(0,total_binslowpass),c(-100,100),type = "n", ylab = "Percentage of polyp samples", xlab="Chromosome", axes=F)
legend('topleft',c("Gains","Losses"),col=c("red","blue"), lty=1,lwd=2)
axis(2,at=seq(-100,100,20),labels=c(100,80, 60, 40,20,0,20,40,60,80,100))
axis(1,at=chr_start_spot, labels=seq(1,22))
lines( gainfreqpolyp,col="red",lwd=2)
lines(-lossfreqpolyp,col="blue",lwd=2)
abline(0,0)
box()

plot(c(0,total_binslowpass),c(-100,100),type = "n", ylab = "Percentage of TCGA samples", xlab="Chromosome", axes=F)
legend('topleft',c("Gains","Losses"),col=c("red","blue"), lty=1,lwd=2)
axis(2,at=seq(-100,100,20),labels=c(100,80, 60, 40,20,0,20,40,60,80,100))
axis(1,at=chr_start_spot, labels=seq(1,22))
lines( gainfreqTCGA,col="red",lwd=2)
lines(-lossfreqTCGA,col="blue",lwd=2)
abline(0,0)
box()

CACRC_data = data.frame(chrom=bins_4401[,1],start=bins_4401[,2], end=bins_4401[,3],gainfreq=gainfreqC,lossfreq=lossfreqC)
write.csv(CACRC_data,file="CA_CRC_CNAfreq.csv")

CNAfreq_data = data.frame(chrom=bins_4401[,1],start=bins_4401[,2], end=bins_4401[,3],gainfreqN=gainfreqN,lossfreqN=lossfreqN,gainfreqL=gainfreqL,lossfreqL=lossfreqL,gainfreqU=gainfreqU,lossfreqU=lossfreqU,gainfreqH=gainfreqH,lossfreqH=lossfreqH,gainfreqC=gainfreqC,lossfreqC=lossfreqC)
write.csv(CNAfreq_data,file="CNAfreq.csv")

### Add TCGA focal amplifications and deletions locations (45 total) to cancer plot above:
cytobands = read.table('../QDNAseq/CACRC_paper_Figures/cytoBand.txt')
## 17 focal gains that were found to be significant p<0.05
focg<-c( "17q12" ,  "8q24.21" , "11p15.5" , "20q13.12" ,"13q12.13" ,"8p12","20q11.21" ,"13q12.2" , "12p13.32" ,"13q22.1" , "8p11.21"  ,"10q22.3" , "17q24.1" , "16p11.2" , "20p11.23", "6p21.1"  , "5p12"   ) 
focgchrom=c(2,1,2,2,2,1,2,2,2,2,1,2,2,2,2,1,1)
## 28 focal losses that were found to be significant p<0.05
focl <- c("3p14.2",   "16p13.3",  "20p12.1" ,"6q26"    , "1p36.11" , "16q23.1" , "3q26.31", "4q22.1"  , "10q11.23", "4q35.1" ,  "6p25.3",   "18q21.2" ,
		"5q11.2" ,  "8p23.2"  , "5q22.2"   ,"1p33"    , "7q31.1" ,  "10q23.2" , "21q21.1" , "15q22.33", "4p16.1"  , "17p12"   , "10q25.2",  "11q22.3" ,
			"1p13.1" ,  "15q21.1" , "15q12"   ,"2q37.3" )
foclchrom = c(1,2,2,1,1,2,1,1,2,1,1,2,1,1,1,1,1,2,2,2,1,2,2,2,1,2,2,1)

gain_regions_ind = list()
for (i in 1:length(focg)){
	chrom = as.character(substr(focg[i],1,focgchrom[i]))
	chrom_ind=which(as.character(substr(cytobands[,1],4,(4+length(chrom))))==chrom)
	temp = which(as.character(substr(cytobands[chrom_ind,4],1,length(cytobands[chrom_ind,4])))==substr(focg[i],(focgchrom[i]+1),nchar(focg[i])))
	gain_regions_ind[[i]]=cytobands[chrom_ind[temp],]
}
names(gain_regions_ind)=focg

loss_regions_ind = list()
for (i in 1:length(focl)){
	chrom = as.character(substr(focl[i],1,foclchrom[i]))
	chrom_ind=which(as.character(substr(cytobands[,1],4,(4+length(chrom))))==chrom)
	temp = which(as.character(substr(cytobands[chrom_ind,4],1,length(cytobands[chrom_ind,4])))==substr(focl[i],(foclchrom[i]+1),nchar(focl[i])))
	loss_regions_ind[[i]]=cytobands[chrom_ind[temp],]
}
names(loss_regions_ind)=focl

gain_bins = rep(0,length(focg))
loss_bins = rep(0,length(focl))


### Plot gains and losses from TCGA 
for (j in 1:length(focg)){
	chrom_focg = as.character(substr(focg[j],1,focgchrom[j]))
	start_focg=gain_regions_ind[[j]][2]
	for (i in 1:(total_binslowpass-1)){
		## Check each of the LPWGS bins if there is a focal gain or a focal loss in that bin
	  	chrom = bins_4401[i,1]
	  	start = bins_4401[i,2]
	  	next_start= bins_4401[(i+1),2]
	  	end = bins_4401[i,3]
		if (chrom_focg==chrom){
			if (start<=start_focg && start_focg <next_start ){
				gain_bins[j]=i
				lines(c(gain_bins[j],gain_bins[j]),c(0,100),col="red",lty=2)
				points(gain_bins[j],100,col="green",pch=4)
				break
			}
		}
	}
}
for (j in 1:length(focl)){
	chrom_focl = as.character(substr(focl[j],1,foclchrom[j]))
	start_focl=loss_regions_ind[[j]][2]
	for (i in 1:(total_binslowpass-1)){
		## Check each of the LPWGS bins if there is a focal gain or a focal loss in that bin
	  	chrom = bins_4401[i,1]
	  	start = bins_4401[i,2]
	  	next_start= bins_4401[(i+1),2]
	  	end = bins_4401[i,3]
		if (chrom_focl==chrom){
			if (start<=start_focl && start_focl <next_start ){
				loss_bins[j]=i
				lines(c(loss_bins[j],loss_bins[j]),c(-100,0),col="blue",lty=2)
				points(loss_bins[j],-100,col="green",pch=4)
				break
			}
		}
	}
}

### Fisher's tests for arm level changes between TCGA sporadic and CA-CRC


### PHYLOGENETIC ANALYSIS
library(phangorn)

call_stats = call_stats_LGDHGD_500kbp
segment_data = matrix(rep(0,4401*5),ncol=4401)
## Patient 1234 from LPWGS
pat_ind=1:5
for (i in 1:4401){
	for (j in pat_ind){
		segment_data[j,i]=call_stats[[j]]['segmented',i]
	}
}

event_data = matrix(rep(0,4401*5),ncol=4401)
for (i in 1:4401){
	for (j in pat_ind){
		event_data[j,i]=call_stats[[j]]['calls',i]
	}
}

event_data=t(event_data)
non <- rep(0,nrow(event_data))
event_data <- cbind(event_data,non)
colnames(event_data) <- c(filenames[1:5], 'Unaltered Genome')


segment_data=t(segment_data)
non <- rep(0,nrow(segment_data))
segment_data <- cbind(segment_data,non)
colnames(segment_data) <- c(filenames[1:5], 'Unaltered Genome')

temp = event_data
tempphy <- as.phyDat(as.data.frame(temp),type="USER",levels=c(-2,-1,0,1))
CNAtreeMP <- pratchet(tempphy)
CNAtreeMP <- root(CNAtreeMP,"Unaltered Genome",resolve.root=TRUE)
CNAtreeMP <- acctran(CNAtreeMP,tempphy);set.seed(120)
CNABStrees <- bootstrap.phyDat(tempphy, pratchet, bs = 100)

##################################### PLOT 4: EVOLUTIONARY TREES FOR PATIENT 1234 (LPWGS) and PATIENT UC06 aka Oxford_IBD1 (WES) #######################################################################################################
dev.new()

plot(CNAtreeMP,edge.width=2,font=2,label.offset=0,tip.color='blue')
add.scale.bar(1500,1.5,lwd=2)




## Patient UC06 from WES
pat_ind=which(substr(names(Exome_data_call_stats),1,11)=="Oxford_IBD1")
#segmentsUC06=length(Exome_data_call_stats[[pat_ind[1]]]$CNt) 
event_dataUC06 = matrix(rep(0,4401*length(pat_ind)),ncol=4401)
for (i in 1:4401){
	chrom = bins_4401[i,1]
  	start = bins_4401[i,2]
  	end = bins_4401[i,3]
	count=1
	for (j in pat_ind){
		exome_chroms = substr(Exome_data_call_stats[[j]]$chromosome,4,length(Exome_data_call_stats[[j]]$chromosome))
  		medianploidy = median(Exome_data_call_stats[[j]]$CNt)
  		if (sum(exome_chroms==chrom)>0){
			chromseg_ind = which(exome_chroms==chrom)
			for (k in chromseg_ind){
				if (Exome_data_call_stats[[j]]$start.pos[k]<start && Exome_data_call_stats[[j]]$end.pos[k]> end ){
					if ((medianploidy-Exome_data_call_stats[[j]]$CNt[k])>=1){
						# loss
						event_dataUC06[count,i]=-1
						break
					}
					else if((Exome_data_call_stats[[j]]$CNt[k]-medianploidy)>=1){
						# gain
						event_dataUC06[count,i]=1
						break
					}
				}
				else if(Exome_data_call_stats[[j]]$start.pos[k]>start && Exome_data_call_stats[[j]]$start.pos[k]< end ){
					if ((medianploidy-Exome_data_call_stats[[j]]$CNt[k])>=1){
						# loss
						event_dataUC06[count,i]=-1
						break
					}
					else if((Exome_data_call_stats[[j]]$CNt[k]-medianploidy)>=1){
						# gain
						event_dataUC06[count,i]=1
						break
					}
				}
				else if(Exome_data_call_stats[[j]]$end.pos[k]>start && Exome_data_call_stats[[j]]$end.pos[k]< end ){
					if ((medianploidy-Exome_data_call_stats[[j]]$CNt[k])>=1){
						# loss
						event_dataUC06[count,i]=-1
						break
					}
					else if((Exome_data_call_stats[[j]]$CNt[k]-medianploidy)>=1){
						# gain
						event_dataUC06[count,i]=1
						break
					}
				}
			}
  		}
  		count=count+1
	}
}

event_dataUC06=t(event_dataUC06)
non <- rep(0,nrow(event_dataUC06))
event_dataUC06 <- cbind(event_dataUC06,non)
colnames(event_dataUC06) <- c(names(Exome_data_call_stats)[pat_ind], 'Unaltered Genome')


temp = event_dataUC06
tempphy <- as.phyDat(as.data.frame(temp),type="USER",levels=c(-1,0,1))
CNAtreeMP <- pratchet(tempphy)
CNAtreeMP <- root(CNAtreeMP,"Unaltered Genome",resolve.root=TRUE)
CNAtreeMP <- acctran(CNAtreeMP,tempphy);set.seed(120)
CNABStrees <- bootstrap.phyDat(tempphy, pratchet, bs = 100)
plotBS(CNAtreeMP,CNABStrees,"phylogram")


dev.new()

plot(CNAtreeMP,edge.width=2,font=2,adj=.5,tip.color='blue')
add.scale.bar(1500,1.5,lwd=2)


## Patient STM008 from WES
pat_ind=which(substr(names(Exome_data_call_stats),1,6)=="STM008")
#segmentsUC06=length(Exome_data_call_stats[[pat_ind[1]]]$CNt) 
event_dataSTM008 = matrix(rep(0,4401*length(pat_ind)),ncol=4401)
for (i in 1:4401){
	chrom = bins_4401[i,1]
  	start = bins_4401[i,2]
  	end = bins_4401[i,3]
	count=1
	for (j in pat_ind){
		exome_chroms = substr(Exome_data_call_stats[[j]]$chromosome,4,length(Exome_data_call_stats[[j]]$chromosome))
  		medianploidy = median(Exome_data_call_stats[[j]]$CNt)
  		if (sum(exome_chroms==chrom)>0){
			chromseg_ind = which(exome_chroms==chrom)
			for (k in chromseg_ind){
				if (Exome_data_call_stats[[j]]$start.pos[k]<start && Exome_data_call_stats[[j]]$end.pos[k]> end ){
					if ((medianploidy-Exome_data_call_stats[[j]]$CNt[k])>=1){
						# loss
						event_dataSTM008[count,i]=-1
						break
					}
					else if((Exome_data_call_stats[[j]]$CNt[k]-medianploidy)>=1){
						# gain
						event_dataSTM008[count,i]=1
						break
					}
				}
				else if(Exome_data_call_stats[[j]]$start.pos[k]>start && Exome_data_call_stats[[j]]$start.pos[k]< end ){
					if ((medianploidy-Exome_data_call_stats[[j]]$CNt[k])>=1){
						# loss
						event_dataSTM008[count,i]=-1
						break
					}
					else if((Exome_data_call_stats[[j]]$CNt[k]-medianploidy)>=1){
						# gain
						event_dataSTM008[count,i]=1
						break
					}
				}
				else if(Exome_data_call_stats[[j]]$end.pos[k]>start && Exome_data_call_stats[[j]]$end.pos[k]< end ){
					if ((medianploidy-Exome_data_call_stats[[j]]$CNt[k])>=1){
						# loss
						event_dataSTM008[count,i]=-1
						break
					}
					else if((Exome_data_call_stats[[j]]$CNt[k]-medianploidy)>=1){
						# gain
						event_dataSTM008[count,i]=1
						break
					}
				}
			}
  		}
  		count=count+1
	}
}

event_dataSTM008=t(event_dataSTM008)
non <- rep(0,nrow(event_dataSTM008))
event_dataSTM008 <- cbind(event_dataSTM008,non)
colnames(event_dataSTM008) <- c(names(Exome_data_call_stats)[pat_ind], 'Unaltered Genome')


temp = event_dataSTM008
tempphy <- as.phyDat(as.data.frame(temp),type="USER",levels=c(-1,0,1))
CNAtreeMP <- pratchet(tempphy)
CNAtreeMP <- root(CNAtreeMP,"Unaltered Genome",resolve.root=TRUE)
CNAtreeMP <- acctran(CNAtreeMP,tempphy);set.seed(120)
CNABStrees <- bootstrap.phyDat(tempphy, pratchet, bs = 100)

dev.new()

plot(CNAtreeMP,edge.width=2,font=2,adj=.5,tip.color='blue')
add.scale.bar(1500,1.5,lwd=2)





call_stats = call_stats_LGDHGD_500kbp
## Patient 3305 from LPWGS
pat_ind=which(substr(filenames,1,4)=="3305")

event_data3305 = matrix(rep(0,4401*length(pat_ind)),ncol=4401)
for (i in 1:4401){
	count=1
	for (j in pat_ind){
		event_data3305[count,i]=call_stats[[filenames[j]]]['calls',i]
		count=count+1
	}
}

event_data3305=t(event_data3305)
non <- rep(0,nrow(event_data3305))
event_data3305 <- cbind(event_data3305,non)
colnames(event_data3305) <- c(filenames[pat_ind], 'Unaltered Genome')


temp = event_data3305

tempphy <- as.phyDat(as.data.frame(temp),type="USER",levels=c(-2,-1,0,1,2))
CNAtreeMP <- pratchet(tempphy)
CNAtreeMP <- root(CNAtreeMP,"Unaltered Genome",resolve.root=TRUE)
CNAtreeMP <- acctran(CNAtreeMP,tempphy);set.seed(1210)
CNABStrees <- bootstrap.phyDat(tempphy, pratchet, bs = 10000)

##################################### PLOT 4: EVOLUTIONARY TREES FOR PATIENT 1234 (LPWGS) and PATIENT UC06 aka Oxford_IBD1 (WES) #######################################################################################################
dev.new()

plot(CNAtreeMP,edge.width=2,font=2,adj=0.5,tip.color='blue')
add.scale.bar(1500,1.5,lwd=2)


##################################### FISHER'S TESTS FOR ARM LEVEL FREQUENCIES #######################################################################################################


chroms = bins_4401[,1]

centromeres = read.csv("../QDNAseq/centromerePositions.csv", header=T)
centro_starts = centromeres[,3]
centro_ends= centromeres[,4]
last_bin_til_centro = rep(0,22)
for (i in 1:22){
	chr_ind = which(chroms==i)
	last_bin_til_centro[i]=which.min(abs(bins_4401[chr_ind,3]-centro_starts[i]))
}
chroms = bins_4401[,1]

## columns: 22 total signif changes found in TCGA

colnames(CACRC_CNAs)= c( "Loss on 1p","Gain on 1q",  "Loss on 4q","Loss on 5q","Gain on 7p", "Gain on 7q", "Gain on 8p", "Loss on 8p", "Gain on 8q", 
							"Gain on 12q",  "Gain on 13q","Loss on 14q","Loss on 15q", "Loss on 17p",  "Loss on 17q", "Loss on 18p", "Loss on 18q",
								 "Gain on 19q", "Gain on 20p", "Loss on 20p", "Gain on 20q", "Loss on 22q")
rownames(CACRC_CNAs)=c("CACRC_CNA_freq","TCGA_CNA_freq")
CACRC_CNAs[2,]=c(0.19,0.17,0.24,0.17,0.47,0.41,0.28,0.5,0.46,0.18,0.56,0.3,0.32,0.56,0.15,0.61,0.66,0.15,0.58,0.32,0.72,0.26)
TCGA_chroms = c(1,1,4,5,7,7,8,8,8,12,13,14,15,17,17,18,18,19,20,20,20,22)
TCGA_chromArms = c('p','q','q','q','p','q','p','p','q','q','q','q','q','p','q','p','q','q','p','p','q','q')


for (j in 1:length(TCGA_chroms)){
  	chr_ind = which(chroms==TCGA_chroms[j])
  	### p arm changes
  	if (TCGA_chromArms[j]=="p"){
  		if (substr(colnames(CACRC_CNAs)[j],1,4)=="Gain"){
	  		CACRC_CNAs[1,j] = mean(gainfreqC[chr_ind[1:last_bin_til_centro[TCGA_chroms[j]]]])
	  	}
	  	else if (substr(colnames(CACRC_CNAs)[j],1,4)=="Loss"){
	  		CACRC_CNAs[1,j] = mean(lossfreqC[chr_ind[1:last_bin_til_centro[TCGA_chroms[j]]]])
	  	}
  	}
  	else if (TCGA_chromArms[j]=="q"){
  	### q arm changes
  		if (substr(colnames(CACRC_CNAs)[j],1,4)=="Gain"){
	  		CACRC_CNAs[1,j] = mean(gainfreqC[chr_ind[(last_bin_til_centro[TCGA_chroms[j]]+1):length(chr_ind)]])
	  	}
	  	else if (substr(colnames(CACRC_CNAs)[j],1,4)=="Loss"){
	  		if (last_bin_til_centro[TCGA_chroms[j]]==1){
	  			CACRC_CNAs[1,j] = mean(lossfreqC[chr_ind[last_bin_til_centro[TCGA_chroms[j]]:length(chr_ind)]])
	  		}
	  		else{
	  			CACRC_CNAs[1,j] = mean(lossfreqC[chr_ind[(last_bin_til_centro[TCGA_chroms[j]]+1):length(chr_ind)]])
	  		}
	  	}
	}
}


CACRC_CNAs[1,]=CACRC_CNAs[1,]/100

total_TCGA = 257
total_CACRC = length(c(C_ind,C_indE,snp_names))
fisherexactp = fisherexactOR= rep(0,length(CACRC_CNAs[1,]))

CNA_details = matrix(rep(0,(length(CACRC_CNAs[1,])*4)),ncol=4)
for (i in 1:length(CACRC_CNAs[1,])){
	temp = matrix(c(CACRC_CNAs[1,i]*total_CACRC,(total_CACRC-CACRC_CNAs[1,i]*total_CACRC),CACRC_CNAs[2,i]*total_TCGA,(total_TCGA-CACRC_CNAs[2,i]*total_TCGA)),byrow=T,nrow = 2)
	print(temp)
	t = fisher.test(round(temp))
	fisherexactp[i]=t$p.value
	fisherexactOR[i] = t$estimate
	CNA_details[i,] = c(CACRC_CNAs[1,i],CACRC_CNAs[2,i],t$estimate,t$p.value )
}
row.names(CNA_details)= colnames(CACRC_CNAs)

#row.names(CNA_details)= c('gain20q','gain13q','gain20p','gain8q','gain5p','gain7p','loss17p','loss5q','loss4p','loss18q','loss4q','loss14q')
colnames(CNA_details)=c('Fraction in CA-CRC','Fraction in S-CRC', 'OR', 'Fisher p')
#write.table(CNA_details, file='Fishertest_CNASinCACRC.xlsx')


## Benjamini-Hochberg (FDR) adjusted p value for multiple comparisons
fisheradjustedp = p.adjust(fisherexactp, method="BH")
ind_signif=which(fisheradjustedp<.05)
print(ind_signif)
CNA_details_a = cbind(CNA_details,fisheradjustedp)
colnames(CNA_details_a)=c('Fraction in CA-CRC','Fraction in S-CRC', 'OR', 'Fisher p', 'Adjusted p')
write.csv(CNA_details_a, file='Fishertest_CNASinCACRC_adj.csv')


### Version with all arm frequencies

## columns: gain on p, loss on p, gain on q, loss on q = 4*22 total -10 no p arms of 21 22 13 14 15
q_only_set=c(13, 14 ,15,21,22)
total_changes = 22*4-10
CACRC_CNAs = matrix(rep(0,total_changes), nrow=2,ncol=total_changes)
colnames(CACRC_CNAs)= c("Gain on 1p", "Loss on 1p", "Gain on 1q", "Loss on 1q","Gain on 2p", "Loss on 2p", "Gain on 2q", "Loss on 2q","Gain on 3p", "Loss on 3p", "Gain on 3q", "Loss on 3q","Gain on 4p", "Loss on 4p", "Gain on 4q", "Loss on 4q",
							"Gain on 5p", "Loss on 5p", "Gain on 5q", "Loss on 5q","Gain on 6p", "Loss on 6p", "Gain on 6q", "Loss on 6q","Gain on 7p", "Loss on 7p", "Gain on 7q", "Loss on 7q","Gain on 8p", "Loss on 8p", "Gain on 8q", "Loss on 8q",
								"Gain on 9p", "Loss on 9p", "Gain on 9q", "Loss on 9q","Gain on 10p", "Loss on 10p", "Gain on 10q", "Loss on 10q",
								"Gain on 11p", "Loss on 11p", "Gain on 11q", "Loss on 11q","Gain on 12p", "Loss on 12p", "Gain on 12q", "Loss on 12q", "Gain on 13q", "Loss on 13q","Gain on 14q", "Loss on 14q",
							 "Gain on 15q", "Loss on 15q","Gain on 16p", "Loss on 16p", "Gain on 16q", "Loss on 16q","Gain on 17p", "Loss on 17p", "Gain on 17q", "Loss on 17q","Gain on 18p", "Loss on 18p", "Gain on 18q", "Loss on 18q",
								"Gain on 19p", "Loss on 19p", "Gain on 19q", "Loss on 19q","Gain on 20p", "Loss on 20p", "Gain on 20q", "Loss on 20q", "Gain on 21q", "Loss on 21q", "Gain on 22q", "Loss on 22q")
rownames(CACRC_CNAs)=c("CACRC_CNA_freq","TCGA_CNA_freq")
TCGA_chroms =1:22
TCGA_numArms=c(rep(2,12),1,1,1,rep(2,5),1,1)
TCGA_chromArms=c(rep(c('p','q'),12),'q','q','q',rep(c('p','q'),5),'q','q')
TCGA_typeArms=rep(c('Gain','Loss'),total_changes/2)
TCGA_arms=read.delim('amp_del_inCRC_TCGA.txt',header=T)
TCGAampfreqs = TCGA_arms[[2]]
TCGAdelfreqs = TCGA_arms[[4]]
CACRC_CNAs[2,]=c(rbind(TCGAampfreqs,TCGAdelfreqs))


count=1
count_arms=1
for (j in 1:length(TCGA_chroms)){
  	chr_ind = which(chroms==TCGA_chroms[j])
  	num_arms = TCGA_numArms[j]
  	### p arm changes
  	for (i in 1:num_arms){
	  	if (TCGA_chromArms[count_arms]=="p"){
		  	CACRC_CNAs[1,count] = mean(gainfreqC[chr_ind[1:last_bin_til_centro[TCGA_chroms[j]]]])
		  	CACRC_CNAs[1,(count+1)] = mean(lossfreqC[chr_ind[1:last_bin_til_centro[TCGA_chroms[j]]]])
		  	count = count + 2
		  	count_arms = count_arms + 1
	  	}
	  	else if (TCGA_chromArms[count_arms]=="q"){
	  	### q arm changes
		  	CACRC_CNAs[1,count] = mean(gainfreqC[chr_ind[(last_bin_til_centro[TCGA_chroms[j]]+1):length(chr_ind)]])
	  		if (last_bin_til_centro[TCGA_chroms[j]]==1){
	  			CACRC_CNAs[1,(count+1)] = mean(lossfreqC[chr_ind[last_bin_til_centro[TCGA_chroms[j]]:length(chr_ind)]])
	  		}
	  		else{
	  			CACRC_CNAs[1,(count+1)] = mean(lossfreqC[chr_ind[(last_bin_til_centro[TCGA_chroms[j]]+1):length(chr_ind)]])
	  		}
	  		count = count + 2
		  	count_arms = count_arms + 1
		  	print(count)
		}
	}
}


CACRC_CNAs[1,]=CACRC_CNAs[1,]/100

total_TCGA = 257
total_CACRC = length(c(C_ind,C_indE,snp_names))
fisherexactp = fisherexactOR= rep(0,length(CACRC_CNAs[1,]))

CNA_details = matrix(rep(0,(length(CACRC_CNAs[1,])*4)),ncol=4)
for (i in 1:length(CACRC_CNAs[1,])){
	temp = matrix(c(CACRC_CNAs[1,i]*total_CACRC,(total_CACRC-CACRC_CNAs[1,i]*total_CACRC),CACRC_CNAs[2,i]*total_TCGA,(total_TCGA-CACRC_CNAs[2,i]*total_TCGA)),byrow=T,nrow = 2)
	print(temp)
	t = fisher.test(round(temp))
	fisherexactp[i]=t$p.value
	fisherexactOR[i] = t$estimate
	CNA_details[i,] = c(CACRC_CNAs[1,i],CACRC_CNAs[2,i],t$estimate,t$p.value )
}
row.names(CNA_details)= colnames(CACRC_CNAs)
colnames(CNA_details)=c('Fraction in CA-CRC','Fraction in S-CRC', 'OR', 'Fisher p')
#write.table(CNA_details, file='Fishertest_CNASinCACRC.xlsx')

## Benjamini-Hochberg (FDR) adjusted p value for multiple comparisons
fisheradjustedp = p.adjust(fisherexactp, method="BH")
ind_signif=which(fisheradjustedp<.05)
print(ind_signif)
CNA_details_a = cbind(CNA_details,fisheradjustedp)
colnames(CNA_details_a)=c('Fraction in CA-CRC','Fraction in S-CRC', 'OR', 'Fisher p', 'Adjusted p')
write.csv(CNA_details_a, file='Fishertest_CNASinCACRC_adj.csv')



### Fisher's tests and Benjamini-Hochberg (FDR) adjusted p values for multiple comparisons for # mutants in TCGA vs. meta-analysis and our CA-CRCs avs TCGA
require(gdata)
EvoCa_v_TCGA=read.xls("/Users/curtiu01/Dropbox/QMUL_2016/LGD_HGD_CACRC/Raw_data_for_BH-FDR.xlsx",header=T)
Meta_v_TCGA=read.xls("/Users/curtiu01/Dropbox/QMUL_2016/LGD_HGD_CACRC/Meta_analysis_raw_data_for_BH-FDR.xlsx",header=T)

total_genes = length(EvoCa_v_TCGA[,1])
fisherexactp = fisherexactOR= rep(0,total_genes)

Mutation_details = matrix(rep(0,(total_genes*4)),ncol=4)
for (i in 1:total_genes){
	temp = matrix(c(EvoCa_v_TCGA[i,2],EvoCa_v_TCGA[i,3],EvoCa_v_TCGA[i,5],EvoCa_v_TCGA[i,6]),byrow=T,nrow = 2)
	print(temp)
	t = fisher.test(round(temp))
	fisherexactp[i]=t$p.value
	fisherexactOR[i] = t$estimate
	Mutation_details[i,] = c(EvoCa_v_TCGA[i,2]/EvoCa_v_TCGA[i,4],EvoCa_v_TCGA[i,5]/EvoCa_v_TCGA[i,7],t$estimate,t$p.value )
}
row.names(Mutation_details)= EvoCa_v_TCGA[,1]
colnames(Mutation_details)=c('Fraction in CA-CRC','Fraction in S-CRC', 'OR', 'Fisher p')
#write.table(CNA_details, file='Fishertest_CNASinCACRC.xlsx')

## Benjamini-Hochberg (FDR) adjusted p value for multiple comparisons
fisheradjustedp = p.adjust(fisherexactp, method="BH")
ind_signif=which(fisheradjustedp<.05)
print(ind_signif)
Mutation_details_a = cbind(Mutation_details,fisheradjustedp)
colnames(Mutation_details_a)=c('Fraction in CA-CRC','Fraction in S-CRC', 'OR', 'Fisher p', 'Adjusted p')
write.csv(Mutation_details_a, file='Fishertest_MutationsinCACRCvTCGA_adj.csv')

## Meta vs TCGA
fisherexactp = fisherexactOR= rep(0,total_genes)

Mutation_details2 = matrix(rep(0,(total_genes*4)),ncol=4)
for (i in 1:total_genes){
	temp = matrix(c(Meta_v_TCGA[i,2],Meta_v_TCGA[i,3],Meta_v_TCGA[i,5],Meta_v_TCGA[i,6]),byrow=T,nrow = 2)
	print(temp)
	t = fisher.test(round(temp))
	fisherexactp[i]=t$p.value
	fisherexactOR[i] = t$estimate
	Mutation_details2[i,] = c(Meta_v_TCGA[i,2]/Meta_v_TCGA[i,4],Meta_v_TCGA[i,5]/Meta_v_TCGA[i,7],t$estimate,t$p.value )
}
row.names(Mutation_details2)= Meta_v_TCGA[,1]
colnames(Mutation_details2)=c('Fraction in Meta analyses','Fraction in S-CRC', 'OR', 'Fisher p')
#write.table(CNA_details, file='Fishertest_CNASinCACRC.xlsx')

## Benjamini-Hochberg (FDR) adjusted p value for multiple comparisons
fisheradjustedp = p.adjust(fisherexactp, method="BH")
ind_signif=which(fisheradjustedp<.05)
print(ind_signif)
Mutation_details2_a = cbind(Mutation_details2,fisheradjustedp)
colnames(Mutation_details2_a)=c('Fraction in Meta analyses','Fraction in S-CRC', 'OR', 'Fisher p', 'Adjusted p')
write.csv(Mutation_details2_a, file='Fishertest_MutationsinMetavTCGA_adj.csv')

## REVISIONS:

## FIGURE 3 with 9 control samples from UC case control study 'Normal normal' = colitis tissue from patients without cancer and who did not progress

normalUCnames = c("29938.S3","51144.D1","16239.D1","49806.D3","16399.R2",'02275.S2',"09663.SP2","29348.R4","34071.S1")
UC_ind = match(normalUCnames,names(all_260_calls))
Control_call_stats = NULL
for (i in normalUCnames){
	Control_call_stats[[i]] = all_260_calls[[i]]
}


all_gain_counts = NULL

for (i in 1:4401){
	all_gain_counts = c(all_gain_counts,rep(i,sum(lapply(Control_call_stats,'[','calls',i)>0)))
}
gainfreq_each_bin = rep(0,4401)
for (i in 1:4401){
	gainfreq_each_bin[i] = sum(all_gain_counts==i)
}

gainfreq = gainfreq_each_bin/length(normalUCnames)*100
#dev.new()
#barplot(gainfreq, ylab= "% of patients with gain")

all_loss_counts = NULL

for (i in 1:4401){
	all_loss_counts = c(all_loss_counts,rep(i,sum(lapply(Control_call_stats,'[','calls',i)<0)))
}
lossfreq_each_bin = rep(0,4401)
for (i in 1:4401){
	lossfreq_each_bin[i] = sum(all_loss_counts==i)
}

lossfreq = lossfreq_each_bin/length(normalUCnames)*100

### PLOT FOR NORMAL NORMAL CONTROLS FROM UC to compare with normal Exomes from above plot in manuscript
plot(c(0,total_binslowpass),c(-100,100),type = "n", ylab = "Percentage of normal samples", xlab="Chromosome", axes=F)
legend('topleft',c("Gains","Losses"),col=c("red","blue"), lty=1,lwd=2)
axis(2,at=seq(-100,100,20),labels=c(100,80, 60, 40,20,0,20,40,60,80,100))
axis(1,at=chr_start_spot, labels=seq(1,22))
lines( gainfreqN,col="red",lwd=2)
lines(-lossfreqN,col="blue",lwd=2)
abline(0,0)
box()










### DO analysis for APC (5q22.2) and ERBB2 (17q12 in general) in cytobands... and more TSG and Oncogenes
APCindcyto = which(cytobands[,1]=="chr5"& cytobands[,4]=="q22.2")
ERBB2indcyto  =which(cytobands[,1]=="chr17"& cytobands[,4]=="q12")
p53indcyto  =which(cytobands[,1]=="chr17"& cytobands[,4]=="p13.1")
KRASindcyto  =which(cytobands[,1]=="chr12"& cytobands[,4]=="p12.1")
PIK3CAindcyto  =which(cytobands[,1]=="chr3"& cytobands[,4]=="q26.32")
FBXW7indcyto  =which(cytobands[,1]=="chr4"& cytobands[,4]=="q31.3")
ARID1Aindcyto  =which(cytobands[,1]=="chr1"& cytobands[,4]=="p36.11")
Bcatindcyto  =which(cytobands[,1]=="chr3"& cytobands[,4]=="p22.1")

start5q22.2 = cytobands[APCindcyto,2]
end5q22.2 = cytobands[APCindcyto,3]
start17q12 = cytobands[ERBB2indcyto,2]
end17q12= cytobands[ERBB2indcyto,3]
startp53 = cytobands[p53indcyto,2]
endp53 = cytobands[p53indcyto,3]
startKRAS = cytobands[KRASindcyto,2]
endKRAS= cytobands[KRASindcyto,3]
startPIK3CA = cytobands[PIK3CAindcyto,2]
endPIK3CA= cytobands[PIK3CAindcyto,3]
startFBXW7 = cytobands[FBXW7indcyto,2]
endFBXW7= cytobands[FBXW7indcyto,3]
startARID1A = cytobands[ARID1Aindcyto,2]
endARID1A = cytobands[ARID1Aindcyto,3]
startBcat = cytobands[Bcatindcyto,2]
endBcat= cytobands[Bcatindcyto,3]
#chr5 111500000 113100000   q22.2  gpos50
#chr17  28800000  35400000     q12  gpos50
## Make list for all call values for each patient
MSIneg_ind=c(which(substr(Exome_tumors$Dx,1,1)=="B"), which(substr(Exome_tumors$Dx,1,3)=="MSI"))

exome_names = Exome_tumors$sampleID[-MSIneg_ind]
exome_dx = Exome_tumors$Dx[-MSIneg_ind ]
exome_set = Exome_tumors$setID[-MSIneg_ind]
exome_patients = unique(exome_set)
Exome_data_call_stats = list()
## these are the columns for the Exome segmented data"
#callnames=c("chromosome","start.pos",	"end.pos"	,"Bf",	"N.BAF",	"sd.BAF","depth.ratio",	"N.ratio",	"sd.ratio",	"CNt",	"A",	"B",	"LPP") 

for (i in 1:length(exome_names)){
	txt_dir = exome_set[i] 
	sample = exome_names[i]
	Exome_data_call_stats[[i]] = read.table(paste('/data/LGD_HGD_CACRC/sequenzaFiles/',txt_dir,'/',sample,'_segments.txt',sep=""),head=TRUE)	
}
names(Exome_data_call_stats) = exome_names

C_indE = which(exome_dx=="C")




exome_C_names = exome_names[C_indE]
cancerLesions_specific_cyto_counts = NULL

count=1
for (i in exome_C_names){
	temp = NULL
	segsA=Exome_data_call_stats[[i]]$A	
	segsB=Exome_data_call_stats[[i]]$B
	medianploidy = median(segsA+segsB)
	#total_size=sum(as.numeric(Exome_data_call_stats[[i]]$end.pos-Exome_data_call_stats[[i]]$start.pos))
	CNA_size=0 
	for ( j in 1: length(segsA)){
		if (Exome_data_call_stats[[i]]$chromosome[j]=="chr1"){
			if (Exome_data_call_stats[[i]]$start.pos[j] > startARID1A && Exome_data_call_stats[[i]]$end.pos[j] <endARID1A){
				temp[["ARID1AA"]] = Exome_data_call_stats[[i]]$A[j]
				temp[["ARID1BB"]] = Exome_data_call_stats[[i]]$B[j]
				print('1 seg fits inside!')
				print(i)
			} 
			else if (Exome_data_call_stats[[i]]$start.pos[j] < startARID1A && Exome_data_call_stats[[i]]$end.pos[j] >endARID1A){
				temp[["ARID1AA"]] = Exome_data_call_stats[[i]]$A[j]
				temp[["ARID1AB"]]= Exome_data_call_stats[[i]]$B[j]
				print('1 segment covers cyto!')
				print(i)
			} 
			else if (Exome_data_call_stats[[i]]$start.pos[j] > startARID1A && Exome_data_call_stats[[i]]$start.pos[j] <endARID1A){
				temp[["ARID1AA"]] = Exome_data_call_stats[[i]]$A[j]
				temp[["ARID1AB"]]= Exome_data_call_stats[[i]]$B[j]
				print('1 segment starts in cyto!')
				print(i)
			} 
		}
		if (Exome_data_call_stats[[i]]$chromosome[j]=="chr3"){
			if (Exome_data_call_stats[[i]]$start.pos[j] > startPIK3CA && Exome_data_call_stats[[i]]$end.pos[j] <endPIK3CA){
				temp[["PIK3CAA"]] = Exome_data_call_stats[[i]]$A[j]
				temp[["PIK3CAB"]] = Exome_data_call_stats[[i]]$B[j]
				print('3q seg fits inside!')
				print(i)
			} 
			else if (Exome_data_call_stats[[i]]$start.pos[j] < startPIK3CA && Exome_data_call_stats[[i]]$end.pos[j] >endPIK3CA){
				temp[["PIK3CAA"]] = Exome_data_call_stats[[i]]$A[j]
				temp[["PIK3CAB"]] = Exome_data_call_stats[[i]]$B[j]
				print('3q segment covers cyto!')
				print(i)
			} 
			else if (Exome_data_call_stats[[i]]$start.pos[j] > startPIK3CA && Exome_data_call_stats[[i]]$start.pos[j] <endPIK3CA){
				temp[["PIK3CAA"]] = Exome_data_call_stats[[i]]$A[j]
				temp[["PIK3CAB"]] = Exome_data_call_stats[[i]]$B[j]
				print('3q segment starts in cyto!')
				print(i)
			} 
			if (Exome_data_call_stats[[i]]$start.pos[j] > startBcat && Exome_data_call_stats[[i]]$end.pos[j] <endBcat){
				temp[["BcatA"]] = Exome_data_call_stats[[i]]$A[j]
				temp[["BcatB"]] = Exome_data_call_stats[[i]]$B[j]
				print('3p seg fits inside!')
				print(i)
			} 
			else if (Exome_data_call_stats[[i]]$start.pos[j] < startBcat && Exome_data_call_stats[[i]]$end.pos[j] >endBcat){
				temp[["BcatA"]] = Exome_data_call_stats[[i]]$A[j]
				temp[["BcatB"]] = Exome_data_call_stats[[i]]$B[j]
				print('3p segment covers cyto!')
				print(i)
			} 
			else if (Exome_data_call_stats[[i]]$start.pos[j] > startBcat && Exome_data_call_stats[[i]]$start.pos[j] <endBcat){
				temp[["BcatA"]] = Exome_data_call_stats[[i]]$A[j]
				temp[["BcatB"]] = Exome_data_call_stats[[i]]$B[j]
				print('3p segment starts in cyto!')
				print(i)
			}
		}
		if (Exome_data_call_stats[[i]]$chromosome[j]=="chr4"){
			if (Exome_data_call_stats[[i]]$start.pos[j] > startFBXW7 && Exome_data_call_stats[[i]]$end.pos[j] <endFBXW7){
				temp[["FBXW7A"]] = Exome_data_call_stats[[i]]$A[j]
				temp[["FBXW7B"]] = Exome_data_call_stats[[i]]$B[j]
				print('4 seg fits inside!')
				print(i)
			} 
			else if (Exome_data_call_stats[[i]]$start.pos[j] < startFBXW7 && Exome_data_call_stats[[i]]$end.pos[j] >endFBXW7){
				temp[["FBXW7A"]] = Exome_data_call_stats[[i]]$A[j]
				temp[["FBXW7B"]]= Exome_data_call_stats[[i]]$B[j]
				print('4 segment covers cyto!')
				print(i)
			} 
			else if (Exome_data_call_stats[[i]]$start.pos[j] > startFBXW7 && Exome_data_call_stats[[i]]$start.pos[j] <endFBXW7){
				temp[["FBXW7A"]] = Exome_data_call_stats[[i]]$A[j]
				temp[["FBXW7B"]]= Exome_data_call_stats[[i]]$B[j]
				print('4 segment starts in cyto!')
				print(i)
			} 
		}
		if (Exome_data_call_stats[[i]]$chromosome[j]=="chr5"){
			if (Exome_data_call_stats[[i]]$start.pos[j] > start5q22.2 && Exome_data_call_stats[[i]]$end.pos[j] <end5q22.2){
				temp[["5q22.2A"]] = Exome_data_call_stats[[i]]$A[j]
				temp[["5q22.2B"]] = Exome_data_call_stats[[i]]$B[j]
				print('5 seg fits inside!')
				print(i)
			} 
			else if (Exome_data_call_stats[[i]]$start.pos[j] < start5q22.2 && Exome_data_call_stats[[i]]$end.pos[j] >end5q22.2){
				temp[["5q22.2A"]] = Exome_data_call_stats[[i]]$A[j]
				temp[["5q22.2B"]]= Exome_data_call_stats[[i]]$B[j]
				print('5 segment covers cyto!')
				print(i)
			} 
			else if (Exome_data_call_stats[[i]]$start.pos[j] > start5q22.2 && Exome_data_call_stats[[i]]$start.pos[j] <end5q22.2){
				temp[["5q22.2A"]] = Exome_data_call_stats[[i]]$A[j]
				temp[["5q22.2B"]]= Exome_data_call_stats[[i]]$B[j]
				print('5 segment starts in cyto!')
				print(i)
			} 
		}
		if (Exome_data_call_stats[[i]]$chromosome[j]=="chr12"){
			if (Exome_data_call_stats[[i]]$start.pos[j] > startKRAS&& Exome_data_call_stats[[i]]$end.pos[j] <endKRAS){
				temp[["KRASA"]] = Exome_data_call_stats[[i]]$A[j]
				temp[["KRASB"]] = Exome_data_call_stats[[i]]$B[j]
				print('12 seg fits inside!')
				print(i)
			} 
			else if (Exome_data_call_stats[[i]]$start.pos[j] < startKRAS && Exome_data_call_stats[[i]]$end.pos[j] >endKRAS){
				temp[["KRASA"]] = Exome_data_call_stats[[i]]$A[j]
				temp[["KRASB"]]= Exome_data_call_stats[[i]]$B[j]
				print('12 segment covers cyto!')
				print(i)
			} 
			else if (Exome_data_call_stats[[i]]$start.pos[j] > startKRAS && Exome_data_call_stats[[i]]$start.pos[j] <endKRAS){
				temp[["KRASA"]] = Exome_data_call_stats[[i]]$A[j]
				temp[["KRASB"]]= Exome_data_call_stats[[i]]$B[j]
				print('12 segment starts in cyto!')
				print(i)
			} 
		}
		if (Exome_data_call_stats[[i]]$chromosome[j]=="chr17"){
			if (Exome_data_call_stats[[i]]$start.pos[j] > start17q12 && Exome_data_call_stats[[i]]$end.pos[j] <end17q12){
				temp[["17q12A"]] = Exome_data_call_stats[[i]]$A[j]
				temp[["17q12B"]] = Exome_data_call_stats[[i]]$B[j]
				print('17q seg fits inside!')
				print(i)
			} 
			else if (Exome_data_call_stats[[i]]$start.pos[j] < start17q12 && Exome_data_call_stats[[i]]$end.pos[j] >end17q12){
				temp[["17q12A"]] = Exome_data_call_stats[[i]]$A[j]
				temp[["17q12B"]] = Exome_data_call_stats[[i]]$B[j]
				print('17q segment covers cyto!')
				print(i)
			} 
			else if (Exome_data_call_stats[[i]]$start.pos[j] > start17q12 && Exome_data_call_stats[[i]]$start.pos[j] <end17q12){
				temp[["17q12A"]] = Exome_data_call_stats[[i]]$A[j]
				temp[["17q12B"]] = Exome_data_call_stats[[i]]$B[j]
				print('17q segment starts in cyto!')
				print(i)
			} 
			if (Exome_data_call_stats[[i]]$start.pos[j] > startp53 && Exome_data_call_stats[[i]]$end.pos[j] <endp53){
				temp[["p53A"]] = Exome_data_call_stats[[i]]$A[j]
				temp[["p53B"]] = Exome_data_call_stats[[i]]$B[j]
				print('17p seg fits inside!')
				print(i)
			} 
			else if (Exome_data_call_stats[[i]]$start.pos[j] < startp53 && Exome_data_call_stats[[i]]$end.pos[j] >endp53){
				temp[["p53A"]] = Exome_data_call_stats[[i]]$A[j]
				temp[["p53B"]] = Exome_data_call_stats[[i]]$B[j]
				print('17p segment covers cyto!')
				print(i)
			} 
			else if (Exome_data_call_stats[[i]]$start.pos[j] > startp53 && Exome_data_call_stats[[i]]$start.pos[j] <endp53){
				temp[["p53A"]] = Exome_data_call_stats[[i]]$A[j]
				temp[["p53B"]] = Exome_data_call_stats[[i]]$B[j]
				print('17p segment starts in cyto!')
				print(i)
			}
		}
	}
	cancerLesions_specific_cyto_counts[[i]]=temp
	#CNA_freqExome[count] = CNA_size/total_size
	count=count+1
}

for (i in snp_names){
	temp = NULL
	pat_id=which(ryan_ascat_segments$SampleID==i)
	segs = ryan_ascat_segments$nA[pat_id]+ryan_ascat_segments$nB[pat_id]
	medianploidy = median(segs)
	for ( j in 1: length(segs)){
		if (ryan_ascat_segments$Chr[pat_id[j]]=="1"){
			if (ryan_ascat_segments$Start[pat_id[j]] > startARID1A && ryan_ascat_segments$End[pat_id[j]]  <endARID1A){
				temp[["ARID1AA"]] = ryan_ascat_segments$nA[pat_id[j]]
				temp[["ARID1AB"]] = ryan_ascat_segments$nB[pat_id[j]]
				print('1 seg fits inside!')
				print(i)
			} 
			else if (ryan_ascat_segments$Start[pat_id[j]] < startARID1A && ryan_ascat_segments$End[pat_id[j]]  >endARID1A){
				temp[["ARID1AA"]] = ryan_ascat_segments$nA[pat_id[j]]
				temp[["ARID1AB"]]= ryan_ascat_segments$nB[pat_id[j]]
				print('1 segment covers cyto!')
				print(i)
			} 
			else if (ryan_ascat_segments$Start[pat_id[j]] < startARID1A && ryan_ascat_segments$Start[pat_id[j]]  <endARID1A){
				temp[["ARID1AA"]] = ryan_ascat_segments$nA[pat_id[j]]
				temp[["ARID1AB"]]= ryan_ascat_segments$nB[pat_id[j]]
				print('1 segment starts in cyto!')
				print(i)
			} 
		}
		if (ryan_ascat_segments$Chr[pat_id[j]]=="3"){
			if (ryan_ascat_segments$Start[pat_id[j]] > startPIK3CA&& ryan_ascat_segments$End[pat_id[j]]  <endPIK3CA){
				temp[["PIK3CAA"]] = ryan_ascat_segments$nA[pat_id[j]]
				temp[["PIK3CAB"]] = ryan_ascat_segments$nB[pat_id[j]]
				print('3q seg fits inside!')
				print(i)
			} 
			else if (ryan_ascat_segments$Start[pat_id[j]] < startPIK3CA && ryan_ascat_segments$End[pat_id[j]]  >endPIK3CA){
				temp[["PIK3CAA"]] = ryan_ascat_segments$nA[pat_id[j]]
				temp[["PIK3CAB"]] = ryan_ascat_segments$nA[pat_id[j]]
				print('3q segment covers cyto!')
				print(i)
			} 
			else if (ryan_ascat_segments$Start[pat_id[j]] < startPIK3CA && ryan_ascat_segments$Start[pat_id[j]]  <endPIK3CA){
				temp[["PIK3CAA"]] = ryan_ascat_segments$nA[pat_id[j]]
				temp[["PIK3CAB"]] = ryan_ascat_segments$nA[pat_id[j]]
				print('3q segment starts in cyto!')
				print(i)
			} 
			if (ryan_ascat_segments$Start[pat_id[j]] > startBcat && ryan_ascat_segments$End[pat_id[j]]  <endBcat){
				temp[["BcatA"]] = ryan_ascat_segments$nA[pat_id[j]]
				temp[["BcatB"]] = ryan_ascat_segments$nB[pat_id[j]]
				print('3p seg fits inside!')
				print(i)
			} 
			else if (ryan_ascat_segments$Start[pat_id[j]] < startBcat&& ryan_ascat_segments$End[pat_id[j]]  >endBcat){
				temp[["BcatA"]] = ryan_ascat_segments$nA[pat_id[j]]
				temp[["BcatB"]] = ryan_ascat_segments$nA[pat_id[j]]
				print('3p segment covers cyto!')
				print(i)
			} 
			else if (ryan_ascat_segments$Start[pat_id[j]] < startBcat && ryan_ascat_segments$Start[pat_id[j]]  <endBcat){
				temp[["BcatA"]] = ryan_ascat_segments$nA[pat_id[j]]
				temp[["BcatB"]] = ryan_ascat_segments$nA[pat_id[j]]
				print('3p segment starts in cyto!')
				print(i)
			} 
		}
		if (ryan_ascat_segments$Chr[pat_id[j]]=="4"){
			if (ryan_ascat_segments$Start[pat_id[j]] > startFBXW7 && ryan_ascat_segments$End[pat_id[j]]  <endFBXW7){
				temp[["FBXW7A"]] = ryan_ascat_segments$nA[pat_id[j]]
				temp[["FBXW7B"]] = ryan_ascat_segments$nB[pat_id[j]]
				print('4 seg fits inside!')
				print(i)
			} 
			else if (ryan_ascat_segments$Start[pat_id[j]] < startFBXW7 && ryan_ascat_segments$End[pat_id[j]]  >endFBXW7){
				temp[["FBXW7A"]] = ryan_ascat_segments$nA[pat_id[j]]
				temp[["FBXW7B"]]= ryan_ascat_segments$nB[pat_id[j]]
				print('4 segment covers cyto!')
				print(i)
			} 
			else if (ryan_ascat_segments$Start[pat_id[j]] < startFBXW7 && ryan_ascat_segments$Start[pat_id[j]]  <endFBXW7){
				temp[["FBXW7A"]] = ryan_ascat_segments$nA[pat_id[j]]
				temp[["FBXW7B"]]= ryan_ascat_segments$nB[pat_id[j]]
				print('4 segment starts in cyto!')
				print(i)
			} 
		}
		if (ryan_ascat_segments$Chr[pat_id[j]]=="5"){
			if (ryan_ascat_segments$Start[pat_id[j]] > start5q22.2 && ryan_ascat_segments$End[pat_id[j]]  <end5q22.2){
				temp[["5q22.2A"]] = ryan_ascat_segments$nA[pat_id[j]]
				temp[["5q22.2B"]] = ryan_ascat_segments$nB[pat_id[j]]
				print('5 seg fits inside!')
				print(i)
			} 
			else if (ryan_ascat_segments$Start[pat_id[j]] < start5q22.2 && ryan_ascat_segments$End[pat_id[j]]  >end5q22.2){
				temp[["5q22.2A"]] = ryan_ascat_segments$nA[pat_id[j]]
				temp[["5q22.2B"]]= ryan_ascat_segments$nB[pat_id[j]]
				print('5 segment covers cyto!')
				print(i)
			} 
			else if (ryan_ascat_segments$Start[pat_id[j]] < start5q22.2 && ryan_ascat_segments$Start[pat_id[j]]  <end5q22.2){
				temp[["5q22.2A"]] = ryan_ascat_segments$nA[pat_id[j]]
				temp[["5q22.2B"]]= ryan_ascat_segments$nB[pat_id[j]]
				print('5 segment starts in cyto!')
				print(i)
			} 
		}
		if (ryan_ascat_segments$Chr[pat_id[j]]=="12"){
			if (ryan_ascat_segments$Start[pat_id[j]] > startKRAS && ryan_ascat_segments$End[pat_id[j]]  <endKRAS){
				temp[["KRASA"]] = ryan_ascat_segments$nA[pat_id[j]]
				temp[["KRASB"]] = ryan_ascat_segments$nB[pat_id[j]]
				print('12 seg fits inside!')
				print(i)
			} 
			else if (ryan_ascat_segments$Start[pat_id[j]] < startKRAS && ryan_ascat_segments$End[pat_id[j]]  >endKRAS){
				temp[["KRASA"]] = ryan_ascat_segments$nA[pat_id[j]]
				temp[["KRASB"]]= ryan_ascat_segments$nB[pat_id[j]]
				print('12 segment covers cyto!')
				print(i)
			} 
			else if (ryan_ascat_segments$Start[pat_id[j]] < startKRAS && ryan_ascat_segments$Start[pat_id[j]]  <endKRAS){
				temp[["KRASA"]] = ryan_ascat_segments$nA[pat_id[j]]
				temp[["KRASB"]]= ryan_ascat_segments$nB[pat_id[j]]
				print('12 segment starts in cyto!')
				print(i)
			} 
		}
		if (ryan_ascat_segments$Chr[pat_id[j]]=="17"){
			if (ryan_ascat_segments$Start[pat_id[j]] > start17q12 && ryan_ascat_segments$End[pat_id[j]]  <end17q12){
				temp[["17q12A"]] = ryan_ascat_segments$nA[pat_id[j]]
				temp[["17q12B"]] = ryan_ascat_segments$nB[pat_id[j]]
				print('17q seg fits inside!')
				print(i)
			} 
			else if (ryan_ascat_segments$Start[pat_id[j]] < start17q12 && ryan_ascat_segments$End[pat_id[j]]  >end17q12){
				temp[["17q12A"]] = ryan_ascat_segments$nA[pat_id[j]]
				temp[["17q12B"]] = ryan_ascat_segments$nA[pat_id[j]]
				print('17q segment covers cyto!')
				print(i)
			} 
			else if (ryan_ascat_segments$Start[pat_id[j]] < start17q12 && ryan_ascat_segments$Start[pat_id[j]]  <end17q12){
				temp[["17q12A"]] = ryan_ascat_segments$nA[pat_id[j]]
				temp[["17q12B"]] = ryan_ascat_segments$nA[pat_id[j]]
				print('17q segment starts in cyto!')
				print(i)
			} 
			if (ryan_ascat_segments$Start[pat_id[j]] > startp53 && ryan_ascat_segments$End[pat_id[j]]  <endp53){
				temp[["p53A"]] = ryan_ascat_segments$nA[pat_id[j]]
				temp[["p53B"]] = ryan_ascat_segments$nB[pat_id[j]]
				print('17p seg fits inside!')
				print(i)
			} 
			else if (ryan_ascat_segments$Start[pat_id[j]] < startp53&& ryan_ascat_segments$End[pat_id[j]]  >endp53){
				temp[["p53A"]] = ryan_ascat_segments$nA[pat_id[j]]
				temp[["p53B"]] = ryan_ascat_segments$nA[pat_id[j]]
				print('17p segment covers cyto!')
				print(i)
			} 
			else if (ryan_ascat_segments$Start[pat_id[j]] < startp53 && ryan_ascat_segments$Start[pat_id[j]]  <endp53){
				temp[["p53A"]] = ryan_ascat_segments$nA[pat_id[j]]
				temp[["p53B"]] = ryan_ascat_segments$nA[pat_id[j]]
				print('17p segment starts in cyto!')
				print(i)
			} 
		}
	}
	cancerLesions_specific_cyto_counts[[i]]=temp
	count=count+1
}

write.csv(t(cancerLesions_specific_cyto_counts),file="/data/CA-CRC_Gut/Specific_gene_counts.csv")

