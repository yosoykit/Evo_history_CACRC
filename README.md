# The evolutionary history of human colitis-associated colorectal cancer

These are the scripts written and used to implement all bioinformatic and statistical analyses by the authors of manuscript:

"The evolutionary history of human colitis-associated colorectal cancer"

Authors: Ann-Marie Baker*^, William Cross*, Kit Curtius1*, Ibrahim Al-Bakir*, Chang-Ho Ryan Choi*, Hayley Davis, Daniel Temko, Sujata Biswas, Pierre Martinez, Marc Williams, James O Lindsay, Roger Feakins, Roser Vega6, Stephen J Hayes, Ian PM Tomlinson, Stuart AC McDonald, Morgan Moorghen, Andrew Silver, James E East, Nicholas A Wright, Lai Mun Wang, Manuel Rodriguez-Justo, Marnix Jansen, Ailsa L Hart^, Simon J Leedham^ and Trevor A Graham^

* joint first authors

^For correspondence:- Ann-Marie BakerEmail: a.m.c.baker@qmul.ac.uk- Simon LeedhamEmail: simonl@well.ox.ac.uk- Trevor GrahamEmail: t.graham@qmul.ac.uk

Make_Bams_lpWGS_pipeline.sh
- Bash script that contains bioninformatics pipeline to create bam files from 81 low-pass whole genome sequenced tissue samples
- takes fastq file format and .txt with sample info as input (see workflow)


QDNAseq_pipeline.R
- R sript to take processed bam files
- Saves segmentation and CNA calll data in a .Rdata file for downstream analyses (see workflow)
  
 ![plot](LPWGS_workflow.png)
  
CACRC_evolution_CNV_analyses.R
-  R script performs the analyses to create:						
	Figures 3B, 3C, 4 												
	Supplementary Tables 9, 10 										
  	All corresponding Results in the Main text regarding the above	


    
