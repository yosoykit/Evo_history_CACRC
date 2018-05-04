# Read Me

There are two scripts in the 'scripts' directory. 

1. quasi_bulk_signatures

This script is in three parts: 

A) Code is provided to read in mutation data for CA-CRC and (TCGA) sporadic CRC tumours,
and classify mutations in each tumour among 96 mutation channels. For the CA-CRC data, 
since there are multiple sequenced regions for the same tumour, mutations are also 
assigned to individual tumour samples/regions in each sample.

B) Code is provided to assign mutational signature activities to each CA-CRC and sporadic 
CRC sample using non-negative least squares regression, implemented in the R package 
'nnls'. This section also contains code for visualising the mutational signature 
assignments in CA-CRC samples.

C) Code is provided to assess the differences in inferred mutational signature composition 
between CA-CRC and sporadic CRC samples.

2. time_evolution

This script is in three parts:

A) Code is provided to read in mutation data for CA-CRC tumours and classify mutations in 
each tumour among 96 mutation channels. Since there are multiple sequenced regions for 
the same tumour, mutations are also assigned to individual tumour samples/regions in each 
sample. Based on the regional mutation assignments mutations are also classified in terms
of timing as 'pre', 'early', or 'late'

B) Code is provided to assign mutational signatures to the mutations in each timing class
for each tumour

C) Code is provided to test and visualise differences in each mutational signature between 
timing classes across samples. 