![alt text](https://github.com/mastoffel/uncarved_blog/tree/master/public/images/mum_lamb.jpg)


### Genetic architecture and lifetime dynamics of inbreeding depression in a wild mammal
Stoffel, M.A, Johnston, S.E., Pilkington, J.G., Pemberton, J.M  
*bioRxiv* [https://doi.org/10.1101/2020.05.27.118877](https://doi.org/10.1101/2020.05.27.118877)  

#### Overview   
This repository contains the analysis code for our paper, in order 1-7.  

*Script 1-3 process data, specifically:*  
**1_ROH_calling:** Calls ROH  
**2_calculate_fitness_parameters_v2:** Transform tables from the Soay sheep database into a table with annual fitness measures  
**3_combine_ROH_and_Fit:** Calculates FROH and combines fitness and FROH data

*Scripts 4-7b contain the main analyses in the paper:*  
**4_ROH_patterns:** First part of the paper, including Figure 1  
**5a_survival_models:** INLA animal models to quantify inbreeding depression in survival  
**5b_survival_models_figure:** Creates Figure 2 based on INLA output  
**6_alt_gwas_annual_survival:** ROH GWAS, needs to run on a cluster  
**7a_gwas_postprocessing:** Takes GWAS results, checks for errors and makes a Manhattan plot    
**7b_gwas_top_snps:** Visualises genetic diversity and GWAS estimates in the genomic vicinity of GWAS peaks  
**99_add_suppl_figures:** Creates some additional figures for Supplementary Material  
**99_make_pcs_for_gwas:** Get PCs in the right format for use in GWAS.  

#### Data
All analysis scripts (4-7b) can be run with data provided in example_data and example_output. For a smooth experience with running these scripts, download the complete repository and change the names of example_data and example_output to data and output, respectively. Examples are based on a random subset of 100 individuals, so results will be different from those reported in the paper.

Complete data will be uploaded upon publication.  

#### Versions and dependencies
All R package dependencies and their versions to run the code are stored in the `renv.lock` file. You can either install packages yourself while running the code, or you can use the `renv` package to setup everything for you. To do so, download or clone this repository somewhere onto your computer. Then install `renv` with:

```r
if (!requireNamespace("remotes"))
  install.packages("remotes")

remotes::install_github("rstudio/renv")
```

You can then simply run `renv::init()` from the directory. This will find the `renv.lock` file and create an renv folder. The folder will contain a private R library containing all the packages used in these scripts. This might take a while as there are a lot of packages to download and install. More info here: [https://rstudio.github.io/renv/](https://rstudio.github.io/renv/)


#### Related repositories
The scripts for some of the data preprocessing are stored in the following repositories:
1) [SNP chip merging and LD decay](https://github.com/mastoffel/sheep)
2) [Genotype imputation and cross-validation](https://github.com/mastoffel/imputation_eddie)
3) [Imputation output postprocessing](https://github.com/mastoffel/imputation_mac)

