###
### This is step2_load.libraries.R, part of the Nature Protocol: 
### 'A computational framework for conducting a three-pronged association study
### to integrate high-throughput ‘-omics’ datasets for the identification of
### potential mechanistic links'.
### 26/7-2018
###

###
### Step 2 - Ensure availability of software packages and satisfaction of
### dependencies
###
### The code in this section simply makes all script-wise dependencies and
### settings available in the present workspace. The analysis was tested using R
### 3.3.3. If packages are available under another version, it should run, but
### specifics of the implementation of each package may change results slightly.
### The provided r-code completes in ~1 hour on a MacBook Pro (2.9GHz quad- core
### 7th-generation Intel Core i7 processor, 16GB 2133MHz LPDDR3 memory) and can
### be paused at any point.
###
### Ensure the following packages (including their indirect dependencies other
### packages needed for their compilation and operation, some of which should be
### installed via Bioconductor) are installed and loaded correctly correctly (in
### every case available via the built-in R and Bioconductor package managers).
###

library (xlsx) ### Saving to spreadsheet
library (data.table) ### Fast read of large files into R
library (WGCNA) ### -	Clustering software. Previously reported work done using v1.34
library (flashClust) ### Clustering software
library (ppcor) ### Partial Spearman correlations, for confounder analysis. Previously reported work done using v1.0
library (gplots) ### Plotting
library (cowplot) ### Plotting; to arrange several plots on the same page
library (ggplot2) ### Plotting
library (plyr) ### Data transformations

### If you have problems installing some packages (e.g. WGCNA), it is likely 
### because you are missing one or more of the Bioconductor packages. 
### First install the dependencies with the code below, then reattempt the 
### installation of WGCNA
# source("https://bioconductor.org/biocLite.R")
# biocLite("impute")
# biocLite("GO.db")
# install.packages("WGCNA)

### The script was testing using the following version of packages:
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
# [1] B iobase_2.34.0        bit64_0.9-7           splines_3.3.3         foreach_1.4.3         gtools_3.5.0         
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
