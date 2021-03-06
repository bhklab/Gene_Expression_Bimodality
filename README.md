# Gene_Expression_Bimodality

## Introduction
This is a repo that contains the code used to do the computational analyses in this paper:

Ba-Alawi, W., Nair, S. K., Li, B., Mammoliti, A., Smirnov, P., Mer, A. S., Penn, L., & Haibe-Kains, B. (2020). **Bimodality of gene expression in cancer patient tumors as interpretable biomarkers for drug sensitivity**. BioRxiv, https://doi.org/10.1101/2020.09.08.288688

The script for the analyses is written in R.


----

## Dependencies:
These R packages need to be installed in order to run the analysis script:
- Biobase
- PharmacoGx
- ggplot2
- GSA
- piano
- fgsea
- dplyr
- RLOBICO


----
## Reproducibility of the Analysis:
- Once the project is downloaded to the user computer, the user needs to navigate to the main directory of the project "Gene_Expression_Bimodality-master".
- Inside the main directory, there is an R script file named "Bimodality_project.R". Running this script will regenerate the different plots used in the paper

<br>
*Important Note 1:* the user needs to set the working directory inside the script file before running it, i.e. setting the working directory to "Gene_Expression_Bimodality-master"
<br>
<br>
*Important Note 2:* RLOBICO package is used to create the logical models. This package requirs the installation of the CPLEX solver by IBM. Please follow instructions on the package's github repo for installation (https://github.com/bhklab/RLOBICO)
<br>
<br>
*Important Note 3:* Data needed to run these analyses can be downloaded from 
https://figshare.com/projects/Bimodality_of_genes_expression_in_cancer_patients_as_interpretable_biomarkers_for_drug_sensitivity/98693
Please download these PSets to folder
<br>


