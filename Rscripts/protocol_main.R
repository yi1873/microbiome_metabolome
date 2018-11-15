### 
### Procedure for three-pronged host-microbiome analysis pipeline
###
### This R script essentially reproduces results from the analysis of MetaHIT
### data in Pedersen, H. K. et al. 'Human gut microbes impact host serum
### metabolome and insulin sensitivity'. Nature 535, 376–381 (2016)., and is
### intended to allow adaptations for analogous work on other datasets.
###
### The sections in this script corresponds to procedural sections in the parent
### protocol: 'A computational framework for conducting a three-pronged
### association study to integrate high-throughput ‘-omics’ datasets for the
### identification of potential mechanistic links'.
### It was adapted from initial scripts by Helle Krogh Pedersen (HKP), Valborg
### Gudmundsdottir (VG) and Henrik Bjørn Nielsen (BN) by HKP, Sofia K. Forslund
### (SKF), VG and Anders Petersen (AP). 
### Last version 26/07/2018.
###

###
### Step 1 - Set working directory
###
### The code will look for input files and deposit output files in
### subdirectories relative to a main directory. On your computer, create a
### directory for the analysis, to remain generic (and without assumption on
### user operating system) we here call it top/. Then create subdirectories to
### obtain the following directory hierarchy:
###
###   top/ 
###     r-code/   // containing this script and any others
###     data/     // containing input files (i.e. all files (unzipped) from the 'example_input.zip'-file).
###     results/  // location of output files
###
### Copy all R-code from the git repository to your top/r-code/ subdirectory (i.e.
### all .R files). Then open the main R-script called protocol_main.R in e.g.
### RStudio. 
### Finally, modify the following command to set the working directory to
### your top/ directory: setwd(“~/top”)
###

setwd ("~/top") ### Change the path to your main working directory. 

###
### Step 2 - Ensure availability of software packages and satisfaction of
### dependencies
###

source ("r-code/step2_load.libraries.R")

###
### Step 3 - Import input files
### 

source ("r-code/step3_import.input.files.R")

### 
### Step 4 - Preprocessing/cleanup of loaded data for sparsity and domain
### limitation
###

source ("r-code/step4_preprocessing.data.for.sparsity.R")

###
### Step 5 - Identify clusters of polar metabolites
###

source ("r-code/step5_identify.WGCNA.clusters.for.metabolites.R")

###
### Step 6 - Link individual polar metabolites to phenotype of interest
###

source ("r-code/step6_associate.metabolites.with.phenotype.R")

###
### Step 7 - Identify clusters of molecular lipids and link individual molecular
### lipids to phenotype of interest
### 

source ("r-code/step7a_identify.WGCNA.clusters.for.lipids.R")
source ("r-code/step7b_associate.lipids.with.phenotype.R")

###
### Step 8 - Link metabolite clusters to phenotype of interest
###

source ("r-code/step8_associate.metabolite.clusters.with.phenotype.R")

###
### Step 9 - Link MGS metagenomic entities to phenotype of interest
###

source ("r-code/step9_associate.MGSs.with.phenotype.R")

###
### Step 10 - Link KEGG functions to phenotype of interest
###

source ("r-code/step10_associate.KEGG.modules.with.phenotype.R")

###
### Step 11 - Save phenotype associations
###

source ("r-code/step11_save.phenotype.associations.R")

###
### Step 12 - Select features with significant differences
###

source ("r-code/step12_select.significant.features.R")

###
### Step 13 - Correlate metabolite clusters to functional metagenomic potentials
### (KOs)
###

source ("r-code/step13_associate.metabolite.clusters.with.KOs.R")

###
### Step 14 - Associate metabolite clusters to functional metagenomic potentials
### (KEGG modules)
###

source ("r-code/step14_associate.metabolite.clusters.with.KEGG.modules.R")

###
### Step 15 - Plot metabolome-microbiome functional analysis results
###

source ("r-code/step15_plot.metabolome.microbiome.functional.analysis.R")

###
### Step 16 (optional) - Export metabolome-microbiome associations for network
### analysis
###

source ("r-code/step16_export.edge.node.files.R")

###
### Step 17 - Leave-one-MGS-out analysis
###

source ("r-code/step17_leave.one.MGS.out.analysis.R")

###
### Step 18 - Extracting top driver species from leave-one-MGS-out analysis for
### each microbiome functional module
###

source ("r-code/step18_extract.top.driver.species.R")

###
### Step 19 - Plotting leave-one-MGS-out results
###

source ("r-code/step19_plot.leave.one.MGS.out.results.R")


print (sessionInfo ())

# R version 3.3.3 (2017-03-06)
# Platform: x86_64-apple-darwin13.4.0 (64-bit)
# Running under: macOS Sierra 10.12.6
# 
# locale:
# [1] da_DK.UTF-8/da_DK.UTF-8/da_DK.UTF-8/C/da_DK.UTF-8/da_DK.UTF-8
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] plyr_1.8.4            cowplot_0.9.1         ggplot2_2.2.1         gplots_3.0.1          ppcor_1.1            
# [6] MASS_7.3-47           flashClust_1.01-2     WGCNA_1.61            fastcluster_1.1.24    dynamicTreeCut_1.63-1
# [11] data.table_1.10.4-3  xlsx_0.5.7            xlsxjars_0.6.1        rJava_0.9-9          
# 
# loaded via a namespace (and not attached):
# [1]  Biobase_2.34.0        bit64_0.9-7           splines_3.3.3         foreach_1.4.3         gtools_3.5.0         
# [6]  Formula_1.2-2         stats4_3.3.3          latticeExtra_0.6-28   blob_1.1.0            fit.models_0.5-14    
# [11] yaml_2.1.14           robustbase_0.92-8     impute_1.48.0         RSQLite_2.0           backports_1.1.1      
# [16] lattice_0.20-35       digest_0.6.12         RColorBrewer_1.1-2    checkmate_1.8.5       colorspace_1.3-2     
# [21] htmltools_0.3.6       preprocessCore_1.36.0 Matrix_1.2-12         pcaPP_1.9-72          pkgconfig_2.0.1      
# [26] GO.db_3.4.0           mvtnorm_1.0-6         scales_0.5.0          gdata_2.18.0          htmlTable_1.11.2     
# [31] tibble_1.3.4          IRanges_2.8.2         nnet_7.3-12           BiocGenerics_0.20.0   lazyeval_0.2.1       
# [36] survival_2.41-3       magrittr_1.5          memoise_1.1.0         doParallel_1.0.11     foreign_0.8-69       
# [41] tools_3.3.3           matrixStats_0.52.2    stringr_1.2.0         S4Vectors_0.12.2      munsell_0.4.3        
# [46] cluster_2.0.6         AnnotationDbi_1.36.2  caTools_1.17.1        rlang_0.1.4           grid_3.3.3           
# [51] iterators_1.0.8       rstudioapi_0.7        htmlwidgets_0.9       robust_0.4-18         labeling_0.3         
# [56] bitops_1.0-6          base64enc_0.1-3       gtable_0.2.0          codetools_0.2-15      DBI_0.7              
# [61] reshape2_1.4.2        rrcov_1.4-3           gridExtra_2.3         knitr_1.17            bit_1.1-12           
# [66] Hmisc_4.1-1           KernSmooth_2.23-15    stringi_1.1.6         parallel_3.3.3        Rcpp_0.12.13         
# [71] rpart_4.1-11          acepack_1.4.1         DEoptimR_1.0-8       

