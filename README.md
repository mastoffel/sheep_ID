### Genetic architecture and lifetime dynamics of inbreeding depression in a wild mammal
Stoffel, M.A, Johnston, S.E., Pilkington, J.G., Pemberton, J.M  
*bioRxiv* [https://doi.org/10.1101/2020.05.27.118877](https://doi.org/10.1101/2020.05.27.118877)  

#### Overview   
This repository contains the analysis code for our paper, in order 1-7.  

Script 1-3 process data, specifically:  
1_ROH_calling: Calls ROH  
2_calculate_fitness_parameters_v2: Transform tables from the Soay sheep database into a table with annual fitness measures  
3_combine_ROH_and_Fit: Calculates FROH and combines fitness and FROH data

Scripts 4-7b contain the main analyses in the paper:  
4_ROH_patterns: First part of the paper, including Figure 1  
5a_survival_models: INLA animal models to quantify inbreeding depression in survival  
5b_survival_models_figure: Creates Figure 2 based on INLA output  
6_alt_gwas_annual_survival: ROH GWAS, needs to run on a cluster  
7a_gwas_postprocessing: Takes GWAS results, checks for errors and makes a Manhattan plot    
7b_gwas_top_snps: Visualises genetic diversity and GWAS estimates in the genomic vicinity of GWAS peaks  
99_add_suppl_figures: Creates some additional figures for Supplementary Material  
99_make_pcs_for_gwas: Get PCs in the right format for use in GWAS.  

#### Data
Full data will be uploaded upon publication.  
All analysis scripts (4-7b) will run with data provided in example_data and example_output. 
To run the code, you either have to change all data/ and output/ filepaths
to example_data/ and example_output/ or you change the folder names to data/ and output/
after downloading the repository. The example data consists of 100 randomly chosen
individuals and all analysis will thus give different results from those reported 
in the manuscript.

#### Related repositories
The scripts for some of the data preprocessing are stored in the following repositories:
1) [SNP chip merging and LD decay](https://github.com/mastoffel/sheep)
2) [Genotype imputation and cross-validation](https://github.com/mastoffel/imputation_eddie)
3) [Imputation output postprocessing](https://github.com/mastoffel/imputation_mac)

