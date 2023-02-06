# KindeyCancerThesisMaster

This project builds a Random Survival Forest (RSF) to predict the Papillary Renal Cell Carcinoma prognosis using a combination of image features extracted from H&E whole-slide tissue images and gene expression profiles. Patient data is acquired from The Cancer Genome Atlas.  

## Code organization and directory description

Not all patients from TCGA have image data and gene expression data simultaneously, so these two types of data are first prepared separately and then intersected using the cBioPortal and Python Data Processing libraries, resulting in 282 patients with all two types of data available.The execution order of programs in each folder is indicated by their filenames, beginning with "01". R and Python programming languages are used in this study.

The datasets used in the study are uploaded in a separate directory. Not all datasets (image and gene sequence data) could be uploaded due to file size, but computed image features using CellProfiler are included. 

Steps in download the Gene Sequence data
1) Proceed to https://www.cbioportal.org/study/summary?id=kirp_tcga
2) Select the intersection of Mutations and mRNA expression (RNA Seq V2 RSEM)
3) Download the clinical ang genomic data of Kidney Renal Papillary Cell Carcinoma (TCGA, Firehose Legacy)
4) Use preferred extraction tools for tar.gz file
5) Extract data_mrna_seq_v2_rsem.txt and convert to desired format

