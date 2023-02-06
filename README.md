# KindeyCancerThesisMaster
note: not all datasets could be uploaded.

This project builds a Random Survival Forest (RSF) to predict the Papillary Renal Cell Carcinoma prognosis using a combination of image features extracted from H&E whole-slide tissue images and gene expression profiles. Patient data is acquired from The Cancer Genome Atlas.  

## Code organization and directory description

Not all patients from TCGA have image data and gene expression data simultaneously, so these two types of data are first prepared separately and then intersected using the cBioPortal and Python Data Processing libraries, resulting in 282 patients with all two types of data available.The execution order of programs in each folder is indicated by their filenames, beginning with "01". R and Python programming languages are used in this study.
