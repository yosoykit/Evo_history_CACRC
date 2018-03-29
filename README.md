# The evolutionary history of human colitis-associated colorectal cancer

These are the scripts written and used to implement all bioinformatic and statistical analyses by the authors of manuscript:

"The evolutionary history of human colitis-associated colorectal cancer"

Authors: Ann-Marie Baker^\*, William Cross\*, Kit Curtius\*, Ibrahim Al-Bakir\*, Chang-Ho Ryan Choi\*, Hayley Davis, Daniel Temko, Sujata Biswas, Pierre Martinez, Marc Williams, James O Lindsay, Roger Feakins, Roser Vega6, Stephen J Hayes, Ian PM Tomlinson, Stuart AC McDonald, Morgan Moorghen, Andrew Silver, James E East, Nicholas A Wright, Lai Mun Wang, Manuel Rodriguez-Justo, Marnix Jansen, Ailsa L Hart^, Simon J Leedham^ and Trevor A Graham^

\* joint first authors

^For correspondence:
- Ann-Marie Baker

Email: a.m.c.baker@qmul.ac.uk

- Simon Leedham

Email: simonl@well.ox.ac.uk

- Trevor Graham

Email: t.graham@qmul.ac.uk

## Folder contents:

### KC_Analyses

### Make_Bams_lpWGS_pipeline.sh
- Bash script that contains bioinformatics pipeline to create bam files from 81 low-pass whole genome sequenced tissue samples
- takes fastq file format and .txt with sample info as input (see workflow)


### QDNAseq_pipeline.R
- R script to take processed bam files
- Saves segmentation and CNA call data in a .Rdata file for downstream analyses (see workflow below)
  
### CACRC_evolution_CNV_analyses.R
-  R script performs the analyses to create:
						
	Figures 3B, 3C, 4 		
											
	Supplementary Tables 9, 10 		
									
  	All corresponding Results in the Main text regarding the above	

<img src="LPWGS_workflow.png" height="60%" width="60%">
  
### WCHC_Analyses
   
   
### 2.0.makeBCFcallingScripts.R

- Makes a shell script for each sample in a list. Shell script is designed to run on a cluster.
- Runs BCFtools to call a latent list of potential SNVs
    
    
### 2.1.makeBCFindelScripts.R

- Same as above but for indels


### 2.2.0.makePlatypusScripts.R

- Makes a shell script to run Platypus jointly for a sample set (defined in the sample list)
- Each chromosome is run separately
- The source command is used to assess the BCFtools variant proposals 
- The mergeFinalVCF script concatenates the resulting vcfs to one


### 2.2.1.makePlatypusIndels.R

- Same as above but for indels


### 2.3.0.processVCF.platypus.R

- Script runs locally
- Filters the Platypus derived variants by coverage (min 10X for each variant)
- Annotates variants using AnnoVar
- Outputs a .txt version of the vcf for germline and somatic variants separately


### 3.6.0.makeSequenzaScripts.R

- Script produces .seqz sequenza files from bams


### 3.6.1.analyseSequenza.R

- Script runs locally to analyse .seqz files as per manual


### 5.0.phylogeneticsPrep.R

- Takes a .txt vcf files and outputs a .nexux file for phylogenetic analyses of variants
- Nexus file can be run in PAUP, Phylip or other compatible software


### 5.1.HomoplasyParse.R

- Runs locally
- Reads in a multiple phylogeny, .tre file and outputs the most parsimonious tree along with statistics


### 5.3.phylogeneticsStats.R

- Runs locally
- Produces a table of tree shape statistics and homoplasy indexes
    
    
