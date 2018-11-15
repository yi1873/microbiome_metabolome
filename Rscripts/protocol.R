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

###
### Step 3 - Import input files
### 
### The following input data is assumed to exist within the top/data subdirectory:
###
### Data-files: 
### -	phenotypes.tab: 
###   File with clinical phenotypes (columns) per
###   individual (rows). Used to test for associations with, or for confounder
###   analysis. Individuals are labeled ‘idv’ followed by a number, e.g. ‘idv001’.
###
### -	metabolomic.tab / lipidomic.tab: 
###   Input data matrix for abundance of 325
###   polar metabolites or 876 molecular lipids per individual. Note that no
###   additional normalization is done in this script, so data is assumed to be
###   comparable in these regards. Such different data types are eventually merged
###   into a single set of metabolite cluster abundances. Individual
###   metabolite/lipids are named M or L (for specifying a polar metabolite or
###   molecular lipid, respectively), followed by a number and lastly the
###   annotation or ‘unknown’ in case of unannotated metabolites/lipids, e.g.
###   ‘M_20_Valine’.
###
### -	MGS_abundance.tab: 
###   File with the abundance (e.g. median gene abundance) of
###   MGSs (columns) per individual (rows). These are assumed to have been
###   rarefied to comparable depth or otherwise normalized. For historical
###   reasons, the MGSs are labeled ‘T2DCAG’ followed by a number, e.g.
###   ‘T2DCAG00001’.
###
### -	KO_abundance.tab: 
###   File with the abundance of each KO (columns) per
###   individual (rows). The data is assumed to be rarefied to comparable depth or
###   otherwise normalized. For this, the software tool rtk58 can be used.
###
### - gene_abundance_sub.tab: 
###   File with the abundance of each catalog gene
###   (subset-version) in each individual assumed to be rarefied to comparable
###   depth or otherwise normalized.
### 
### Annotation-files:
### -	cluster_mapping_file.tab: 
###   Input file with annotation for metabolite
###   clusters, as available from curation of data in the specific dataset. The
###   WGCNA clustering algorithm names the generated clusters with color codes.
###   This mapping file simply facilitate renaming to more meaningful cluster
###   descriptions. Here, the serum polar metabolite and serum molecular lipid
###   clusters are labelled M01–M35 and L01–L39, respectively, and collectively
###   termed metabolite clusters.
###
### -	MGS_taxonomy.tab: 
###   File with taxonomic annotation of the MGSs (rows) used
###   in the analysis. Each row contains the following information:
###   species_taxonomy, species_pct, genus_taxonomy, genus_pct, family_taxonomy,
###   family_pct, order_taxonomy, order_pct, phylum_taxonomy, phylum_pct, where
###   x_pct is the percentage of the MGS genes that can be annotated (by sequence
###   similarity) to the taxonomy of the MGS. If no taxonomy can be assigned to
###   the MGS, the value will be NA.
###
### -	KEGG_modules.tab: 
###   File containing definition of KEGG gene functional
###   modules (rows) specifying which KO gene groups constitute each KEGG module.
###   The first two columns contain KEGG module entry (number) and name, the third
###   column list all KOs separated by semicolon. Any other functional annotation
###   used analogously could be swapped in instead of the name. This file can be
###   obtained by downloading KEGG modules from
###   www.genome.jp/kegg-bin/get_htext?ko00002.keg (www.genome.jp/kegg/ --> KEGG
###   MODULE --> KEGG modules); download htext and then running the provided
###   script ‘parse_kegg.pl’ (after changing input filename).
###
### -	KO_to_MGS.tab: 
###   File listing for each KO (rows) the MGSs (space-separated)
###   it is a member of. This is based on the KO annotation of the gene catalog
###   (gene_to_KO.tab) and the information, what gene from the gene catalog is
###   within each MGS as given in the list MGS_to_gene.tab.
###
### -	gene_to_KO.tab: 
###   File containing the KO annotation (if any) per gene (rows)
###   in the gene catalog (constituting 7,328,469 genes). Genes are labeled
###   ‘RefCat620’ followed by a number (1…7328469), e.g. ‘RefCat620.1’. Used to
###   create KO_to_MGS.tab. Note, for the purpose of this protocol and to
###   considerably reduce size of input data, only the subset of catalogue genes
###   with KO annotation (n = 2,205,769) are provided in files with gene abundance
###   or annotation (gene_to_KO.tab, MGS_to_gene.tab and gene_abundance_sub.tab)
###   as only those genes are used in the driver-species analysis.
###
### -	MGS_to_gene.tab: 
###   List of genes binned into a given MGS (rows). Used to
###   create KO_to_MGS.tab. Rather than gene names the file contains the index
###   value of the gene. i.e. position of the gene in the gene catalogue
###   (subset-version, thus 1…2205769).
###
### All demonstration files are tab-delimited text files, but other formats
### would equally work after modifying the respective file-import commands in
### the R-script. 
###
### For all demonstration data, pseudonymised sample names were re-randomised to
### generate anonymised data.
###

options (stringsAsFactors = FALSE)

### - phenotypes.tab

phenotypes = read.table (file = "data/phenotypes.tab", row.names = 1, header = T, sep = "\t")
ctrl = rownames (subset (phenotypes, Diabetes == "nonDiabetic")) ### set of control samples for analyses restricted to non-diabetic individuals.

### - metabolomic.tab

metabolomic = read.table (file = "data/metabolomic.tab", row.names = 1, header = T, sep = "\t")

### - lipidomic.tab

lipidomic = read.table (file = "data/lipidomic.tab", row.names = 1, header = T, sep = "\t")

### - cluster_mapping_file.tab

cluster_mapping_file = read.table (file = "data/cluster_mapping_file.tab", header = T, sep = "\t", row.names = 1)
cluster_mapping_file$label =  sapply (rownames (cluster_mapping_file),  function (x) paste (cluster_mapping_file [x, "New_Name"], cluster_mapping_file [x, "Description"], sep = ": "))

### - MGS_abundance.tab

mgs_abundance = read.table (file = "data/MGS_abundance.tab", sep = "\t", row.names = 1, header = T)

### - MGS_taxonomy.tab

mgs_taxonomy  = read.table (file = "data/MGS_taxonomy.tab", sep = "\t", row.names = 1, header = T)

### - KEGG_modules.tab

tmp = read.table ("data/KEGG_modules.tab", sep = "\t")
koann = strsplit (tmp[,3], split = ";")
names (koann) = tmp[,1]

module_mapping = tmp[,2] ### description of Kegg modules
names (module_mapping) = tmp[,1] ; rm (tmp)

### remove KEGG references in square brackets for more clean names for plotting
module_mapping_clean = sapply(module_mapping, function(x) strsplit(x, " \\[")[[1]][1])

### - KO_abundance.tab

ko_abundance = read.table (file = "data/KO_abundance.tab", sep = "\t", row.names = 1, header = T)

### - KO_to_MGS.tab

tmp = read.table ("data/KO_to_MGS.tab", sep = "\t", strip.white = T)
KO2MGS = strsplit (tmp[,2], split = " ")
names (KO2MGS) = tmp[,1] ; rm (tmp)

### - gene_to_KO.tab

gene2KO = read.table ("data/gene_to_KO.tab", sep = "\t", row.names = 1, header = T)
KO2gene = tapply (1:nrow (gene2KO), gene2KO[,1], c) ; rm (gene2KO)

### - MGS_to_gene.tab

tmp = read.table ("data/MGS_to_gene.tab", sep = "\t", strip.white = T)
MGS2gene = strsplit (tmp[,2], split = " ")
names (MGS2gene) = tmp[,1] ; rm (tmp)

### - gene_abundance_sub.tab

gene_abundance_sub = data.frame( fread ("data/gene_abundance_sub.tab", sep = "\t", header = T), row.names = 1)

### 
### Step 4 - Preprocessing/cleanup of loaded data for sparsity and domain
### limitation
###
### The following commands restrict the input data to account for method
### limitations.
###

### MGS sparsity filter step - exclude MGSs that occur in <3 of the control individuals
test = apply (mgs_abundance [intersect (ctrl, rownames (mgs_abundance)),], 2, function(x) length (x [x != 0])) 
incl = names (test [test >= 3])
length (incl)
MGSs = incl
mgs_abundance = mgs_abundance [,incl]
mgs_taxonomy = mgs_taxonomy [incl,]
rm (test, incl)
### MGS sparsity filter step done

### Filtering out KEGG modules that are eukaryotic only
euk = unique (c ("M00352", "M00355", "M00354", "M00285", "M00295", "M00341", 
                 "M00177", "M00160", "M00359", "M00391", "M00182", "M00340", 
                 "M00359", "M00182", "M00160", "M00177", "M00391", "M00180",  ### "eukaryotes" in name
                 "M00351", "M00352", "M00355", "M00354", "M00353", "M00427")) ### spliceosome or nuclear functions
incl = setdiff (names (koann), euk)
koann = koann [incl] ; rm (incl)
### Eukaryote filtering done

### KO sparsity filter step - exclude KOs that occur in <3 of the control individuals
test = apply (ko_abundance [intersect (ctrl, rownames (ko_abundance)),], 2, function (x) length (x [x != 0]))  
incl = names (test [test >= 3])
length (incl)
ko_abundance = ko_abundance [, incl]
rm (test, incl)
### KO sparsity filter step done

### Determine control sample set and ensure there are no missing values for
### phenotypes tested, as algorithms used cannot handle this well.

ctrl.no.na = ctrl [! is.na (phenotypes [ctrl, "Homa.IR"])]

###
### Step 5 - Identify clusters of polar metabolites
###
### A central part of the protocol is dimensionality reduction of the high-
### dimensional input data. Here the space of metabolite measurements is reduced
### by identifying clusters of co-abundant metabolite peaks, analogously to the
### identification of MGSs as co-abundant gene groups. This is done using
### algorithms developed for gene expression analysis, specifically WGCNA.
###
### Importantly, effective WGCNA performance depends on parameter settings for
### cluster reconstruction which will vary to some extent between datasets.
### Establishing these for a particular dataset requires some inspection of
### resulting clustering behavior under parameter exploration. The full
### procedure for this lies beyond the scope of the present Protocols, but the
### reader is referred to the WGCNA main documentation.
###
### Thus, for optimal performance, the steps below should be performed for
### different values and the resulting curves inspected, hereafter parameters
### should be specified as appropriate. The examples below corresponds to
### optimal choices for the dataset and analysis used in Pedersen et al., 2016.
###
### In this step, WGCNA clustering is performed separately on lipidomic and
### metabolomic measurements to detect cluster of densely connected
### metabolites/lipids. The metabolite/lipid profiles constituting a given
### cluster are summarized by the 1st principal component of the
### metabolite/lipid abundance matrix (‘Module Eigen-metabolite/lipid’); i.e.
### basically a weighted average abundance profile. Prior to this, optimal
### parameters for WGCNA should be established for the dataset being analyzed.
###

### Settings for WGCNA generally
cor_method          = "spearman" ### for association with clinical parameters
corFun_tmp          = "bicor"
cluster_method      = "average"
corOptions_list     = list (use = 'pairwise.complete.obs') 
corOptions_str      = "use = 'pairwise.complete.obs'"
BH_pval_asso_cutoff = 0.05
NetworkType         = "signed" ### Signed-network (as the PC1 and median profile does not make sense as a summary measure of a cluster with anticorrelated metabolites.)

### Specify data and parameters 
ID_common_all_samples_meta = intersect (rownames (phenotypes), rownames (metabolomic)) ### only include individuals with information for both domains.
dat = as.data.frame (metabolomic [ID_common_all_samples_meta,])
dat_tmp = log2 (dat) ### work in logarithmic space

### Settings for WGCNA on polar metabolite measurements
### Once these are established the steps below can be run
RsquareCut_val      = 0.89 ### usually ranges 0.80-0.95 but requires inspecting curves
mergingThresh       = 0.20 ### Maximum dissimilarity of module eigengenes (i.e. 1-correlation) for merging modules.
minModuleSize       = 3 ### minimum number of metabolites constituting a cluster
SoftPower           = 13 ### beta-value, main parameter to optimize

### Calculate weighted adjacency matrix
A = adjacency (dat_tmp, power = SoftPower, type = NetworkType, corFnc = corFun_tmp, corOptions = corOptions_str)
colnames (A) = rownames (A) = colnames (dat_tmp)
### Define dissimilarity based on topological overlap
dissTOM = TOMdist (A, TOMType = NetworkType)
colnames (dissTOM) = rownames (dissTOM) = colnames (dat_tmp)
### Hierarchical clustering
metaTree = flashClust (as.dist (dissTOM), method = cluster_method)
### Define modules by cutting branches
moduleLabels1 = cutreeDynamic (dendro = metaTree, distM = dissTOM, method = "hybrid", deepSplit = 4, pamRespectsDendro = T, minClusterSize = minModuleSize)
moduleLabels1 = labels2colors (moduleLabels1)
### Automatically merge highly correlated modules
merge = mergeCloseModules (dat_tmp, moduleLabels1, corFnc = corFun_tmp, corOptions = corOptions_list, cutHeight = mergingThresh)
### Determine resulting merged module colors
moduleLabels2 = merge$colors
### Establish eigengenes of the newly merged modules, used for cluster overall abundances
MEs = merge$newMEs
### Choose final module assignments
moduleColorsMeta = moduleLabels2
names (moduleColorsMeta) = colnames (dat_tmp)
MEsMeta = orderMEs (MEs)
rownames (MEsMeta) = rownames (dat_tmp)

### Determine relevant descriptive statistics of established clusters
### kIN: within-module connectivity, determined by summing connectivity with all
###      other metabolites in the given cluster.
### kME: bicor-correlation between the metabolite profile and module eigenvector; 
### both measures of intramodular hub-metabolite status.
kIN <-      vector (length = ncol (dat_tmp)); names (kIN) = colnames (dat_tmp)
kME <-      vector (length = ncol (dat_tmp)); names (kME) = colnames (dat_tmp)
modules <-  vector (length = ncol (dat_tmp)); names (modules) = colnames (dat_tmp)

for (module in names (table (moduleColorsMeta))) {   

	all.metabolites = names (dat_tmp)
	inModule = (moduleColorsMeta == module)
	module.metabolites = names (moduleColorsMeta [inModule])
  modules [module.metabolites] = module 
	kIN [module.metabolites] = sapply (module.metabolites, function (x) sum (A [x, module.metabolites]) - 1)
  datKME = signedKME (dat_tmp, MEsMeta, corFnc = corFun_tmp, corOptions = corOptions_str)
	rownames (datKME) = colnames (dat_tmp)
	kME [module.metabolites] = datKME [module.metabolites, paste ("kME", module, sep = "")]   

}
output = data.frame ("module" = modules, "kME" = kME, "kIN" = kIN, "cluster_name" = sapply (modules, function (m) cluster_mapping_file [paste0 ("M_ME", m), "New_Name"]))

###
### Step 6 - Link individual polar metabolites to phenotype of interest
###
### The dimensionality reduction approach (detailed in Stage 3 (Step 8-12))
### hinges on identifying features (e.g. metabolites and lipids clusters, etc.)
### that are associated with a host phenotype of interest. In the example work
### we describe here, this was insulin resistance as assessed by the HOMA-IR
### measurement.
###
### In this step, for reference, the individual polar metabolites are linked
### with HOMA-IR, both directly and under adjustment for a potential confounder
### variable. In this case, such de-confounding was done for body mass index
### (BMI). In case of a binary phenotype variable, one can, for example,
### substitute the Spearman correlation test with a Mann-Whitney U (MWU) test.
###

tmpMat = array (NA, c (ncol (metabolomic), 2, 2))
dimnames (tmpMat) [[1]] = colnames (metabolomic)
dimnames (tmpMat) [[2]] = c ("HOMA.ir", "HOMA.ir_BMI.adj.partial") 
dimnames (tmpMat) [[3]] = c ("estimate", "p.value")

### Associating individual polar metabolites with HOMA-IR without de-confounding
### for BMI
tmpMat [, "HOMA.ir", c ("estimate", "p.value")] =
	t (apply (metabolomic [ctrl.no.na,], MARGIN = 2, FUN = function (x) 
	  unlist (cor.test (phenotypes [ctrl.no.na, "Homa.IR"], x,
		method = cor_method, use = "pairwise.complete.obs")[c ("estimate", "p.value")])))       

### Associating individual polar metabolites with HOMA-IR while de-confounding
### for BMI
tmpMat [, "HOMA.ir_BMI.adj.partial", c ("estimate", "p.value")] =
	t (apply (metabolomic [ctrl.no.na,], MARGIN = 2, FUN = function (x) 
	  unlist (pcor.test (x = x, y = phenotypes [ctrl.no.na, "Homa.IR"], z = phenotypes [ctrl.no.na, "BMI.kg.m2"],
		method = cor_method) [c ("estimate", "p.value")])))       

### Sort by cluster_name and then decreasing values of kIN
output2 = cbind (tmpMat[, "HOMA.ir", c ("estimate", "p.value")], "p.adjust" = p.adjust (tmpMat [, "HOMA.ir", "p.value"], method = "BH"), 
                 tmpMat[, "HOMA.ir_BMI.adj.partial", c ("estimate", "p.value")], "p.adjust" = p.adjust (tmpMat [, "HOMA.ir_BMI.adj.partial", "p.value"], method = "BH"))
colnames (output2) = paste (rep (c ("HOMA.ir", "HOMA.ir_BMI.adj.partial"), each = 3), colnames (output2), sep = "_")
output3 = cbind (output, output2)
output3 = output3 [with (output3, order (cluster_name, -kIN)), ]

### Write to file
write.table (output3, file = "results/individual_metabolites.txt", sep = "\t", col.names = NA, quote = F, row.names = T)
rm (output, output2, output3, tmpMat, dat, dat_tmp)

###
### Step 7a - Identify clusters of molecular lipid
###
### Analogous to Step 5. First identify optimal parameters for WGCNA, then
### establish clusters.
###

### Specify data and parameters
dat = as.data.frame (lipidomic [rownames (phenotypes),])
dat_tmp = log2 (dat)

### Settings for WGCNA on molecular lipids measurements
### Once these are established the steps below can be run
RsquareCut_val = 0.90
mergingThresh = 0.25
minModuleSize = 5
SoftPower = 14

### Calculate weighted adjacency matrix
A = adjacency (dat_tmp, power = SoftPower, type = NetworkType, corFnc = corFun_tmp, corOptions = corOptions_str)
### Define dissimilarity based on topological overlap
dissTOM = TOMdist (A, TOMType = NetworkType)
### Hierarchical clustering
lipidTree = flashClust (as.dist (dissTOM), method = cluster_method)
### Define modules by cutting branches
moduleLabels1 = cutreeDynamic (dendro = lipidTree, distM = dissTOM, method = "hybrid", deepSplit = 4, pamRespectsDendro = T, minClusterSize = minModuleSize)
moduleLabels1 = labels2colors (moduleLabels1)
### Automatically merge highly correlated modules
merge = mergeCloseModules (dat_tmp, moduleLabels1, corFnc = corFun_tmp, corOptions = corOptions_list, cutHeight = mergingThresh)
### Determine resulting merged module colors
moduleLabels2 = merge$colors
### Establish eigengenes of the newly merged modules, used for cluster overall
### abundances
MEs = merge$newMEs
### Choose module assignments
moduleColorsLipid = moduleLabels2
names (moduleColorsLipid) = colnames (dat_tmp)
MEsLipid = orderMEs (MEs)
rownames (MEsLipid) = rownames (dat_tmp)

### Determine relevant descriptive statistics of established clusters
kIN <-      vector (length = ncol (dat_tmp)); names (kIN) = colnames (dat_tmp)
kME <-      vector (length = ncol (dat_tmp)); names (kME) = colnames (dat_tmp)
modules <-  vector (length = ncol (dat_tmp)); names (modules) = colnames (dat_tmp)
 
for (module in names (table (moduleColorsLipid))) {   

	all.lipids = names (dat_tmp)
	inModule = (moduleColorsLipid == module)
	module.lipids = names (moduleColorsLipid [inModule])
  modules [module.lipids] = module 
	kIN [module.lipids] = sapply (module.lipids, function (x) sum (A [x, module.lipids]) -1)
	datKME = signedKME (dat_tmp, MEsLipid, corFnc = corFun_tmp, corOptions = corOptions_str)
	rownames (datKME) = colnames (dat_tmp)
	kME [module.lipids] = datKME [module.lipids, paste ("kME", module, sep = "")]   

}
output = data.frame ("module" = modules, "kME" = kME, "kIN" = kIN, "cluster_name" = sapply (modules, function (l) cluster_mapping_file [paste0 ("L_ME", l), "New_Name"]))

###
### Step 7b - Link individual molecular lipids to phenotype of interest
###
### Analogous to Step 6.
### The resulting metabolite and lipid clusters are thereafter merged into a
### combined dataset, collectively termed 'metabolite clusters', for downstream
### analyses.
###

tmpMat = array (NA, c (ncol (lipidomic), 2,  2))
dimnames (tmpMat) [[1]] = colnames (lipidomic)
dimnames (tmpMat) [[2]] = c ("HOMA.ir", "HOMA.ir_BMI.adj.partial") 
dimnames (tmpMat) [[3]] = c ("estimate", "p.value")

### Associating individual molecular lipids with HOMA-IR without de-confounding
### for BMI
tmpMat [, "HOMA.ir", c ("estimate", "p.value")] =
	t (apply (lipidomic [ctrl.no.na,], MARGIN = 2, FUN = function (x) 
	  unlist (cor.test (phenotypes [ctrl.no.na, "Homa.IR"], x,
		method = cor_method, use = "pairwise.complete.obs") [c ("estimate", "p.value")])))       

### Associating individual molecular lipids with HOMA-IR while de-confounding
### for BMI
tmpMat[, "HOMA.ir_BMI.adj.partial", c ("estimate", "p.value")] =
	t (apply (lipidomic [ctrl.no.na,], MARGIN = 2, FUN = function (x) 
	  unlist (pcor.test (x = x, y = phenotypes [ctrl.no.na, "Homa.IR"], z = phenotypes [ctrl.no.na, "BMI.kg.m2"],
		method = cor_method) [c ("estimate", "p.value")])))       

### Sort by cluster_name and then decreasing values of kIN
output2 = cbind (tmpMat[, "HOMA.ir", c ("estimate", "p.value")], "p.adjust" = p.adjust (tmpMat [, "HOMA.ir", "p.value"], method = "BH"), 
                 tmpMat[, "HOMA.ir_BMI.adj.partial", c ("estimate", "p.value")], "p.adjust" = p.adjust (tmpMat [, "HOMA.ir_BMI.adj.partial", "p.value"], method = "BH"))
colnames (output2) = paste (rep (c ("HOMA.ir", "HOMA.ir_BMI.adj.partial"), each = 3), colnames (output2), sep = "_")
output3 = cbind (output, output2)
output3 = output3 [with (output3, order (cluster_name, -kIN)), ]

### Write to file
write.table (output3, file = "results/individual_lipids.txt", sep = "\t", col.names = NA, quote = F, row.names = T)
rm (output, output2, output3, tmpMat, dat, dat_tmp)

### Create a joint data frame of metabolite/lipid cluster eigengene
### equivalents/effective abundances ('MEsMetLip')
MEsMetLip = cbind (MEsMeta [rownames (MEsMeta),], MEsLipid [rownames (MEsMeta),])
### Rename module names (columnames) from 'colors' to numbers
colnames (MEsMetLip) = cluster_mapping_file [c (paste0 ("M_", colnames (MEsMeta)), paste0 ("L_", colnames (MEsLipid))), "New_Name"]
### Exclude the two bin-clusters (i.e. "M_remaining" and "L_remaining")
MEsMetLip = subset (MEsMetLip, select = c (-M_remaining, -L_remaining))
### Save MEsMetLip to file, order by module number
write.table (MEsMetLip [ , order (colnames (MEsMetLip))], file = "results/MEs_metabolite_clusters.txt", sep = "\t", row.names = T, col.names = NA, quote = F)

###
### Step 8 - Link metabolite clusters to phenotype of interest
###
### Analogous to Step 6 (and 7b). 
### This is a core analysis step generating associations between the
### integrated/clustered –omics data and a clinically interesting phenotype. In
### Pedersen et al., 2016, this was insulin resistance (HOMA-IR measurement),
### but any phenotype is possible, as is checking against other –omics spaces or
### overall –omics measurements such as gut diversity or enterotype. This
### analysis can further be conducted controlling for confounders such as BMI in
### the analysis we previously reported, by performing tests with partial
### correlations or extended to binary phenotype variables by substituting tests
### of Spearman correlation with e.g. MWU tests.
###

cor_HOMA.IR <- list () ### Data structure for storing results of correlation tests under different setups

tmpMat = array (NA, c (ncol (MEsMetLip), 2, 2))
dimnames (tmpMat) [[1]] = names (MEsMetLip)
dimnames (tmpMat) [[2]] = c ("HOMA.ir", "HOMA.ir_BMI.adj.partial") 
dimnames (tmpMat) [[3]] = c ("estimate", "p.value")

### Associating metabolite clusters with HOMA-IR without de-confounding for BMI    
tmpMat [, "HOMA.ir", c ("estimate", "p.value")] =
	t (apply (MEsMetLip [ctrl.no.na,], MARGIN = 2, FUN = function (x) 
	  unlist (cor.test (phenotypes [ctrl.no.na, "Homa.IR"], x,
		method = cor_method, use = "pairwise.complete.obs") [c ("estimate", "p.value")])))       

### Associating metabolite clusters with HOMA-IR while de-confounding for BMI    
tmpMat [, "HOMA.ir_BMI.adj.partial", c ("estimate", "p.value")] =
	t (apply (MEsMetLip [ctrl.no.na,], MARGIN = 2, FUN = function (x) 
	  unlist (pcor.test (x = x, y = phenotypes [ctrl.no.na, "Homa.IR"], z = phenotypes [ctrl.no.na, "BMI.kg.m2"],
		method = cor_method) [c ("estimate", "p.value")])))       

cor_HOMA.IR [["metlip"]] <- tmpMat
rm (tmpMat)

###
### Step 9 - Link MGS metagenomic entities to phenotype of interest
###
### Analogous to Step 6 but for metagenomic taxonomic data.
###

tmpMat = array (NA, c (ncol (mgs_abundance), 2, 2))
dimnames (tmpMat) [[1]] = colnames (mgs_abundance)
dimnames (tmpMat) [[2]] = c ("HOMA.ir", "HOMA.ir_BMI.adj.partial")
dimnames (tmpMat) [[3]] = c ("estimate", "p.value")

### MGS sparsity filter step - exclude MGSs that occur in <3 of the control
### individuals that are actually used in the analysis. Repeat the step since we
### have excluded control individuals with missing information for HOMA-IR
ctrl.no.na.4MGS = intersect (ctrl.no.na, rownames (mgs_abundance)) # length=275
test = apply (mgs_abundance [intersect (ctrl.no.na.4MGS, rownames (mgs_abundance)),], 2, function (x) length (x [x != 0]))
mgs.subset = names (test[test>=3])
rm (test)
### Sparsity filter done

### Associating MGSs with HOMA-IR without de-confounding for BMI    
tmpMat [mgs.subset, "HOMA.ir", c ("estimate", "p.value")] =
	t (apply (mgs_abundance [ctrl.no.na.4MGS, mgs.subset], MARGIN = 2, FUN = function (x) 
	  unlist (cor.test (phenotypes [ctrl.no.na.4MGS, "Homa.IR"], x,
		method = cor_method, use = "pairwise.complete.obs") [c ("estimate", "p.value")])))       

### Associating MGSs with HOMA-IR while de-confounding for BMI    
tmpMat [mgs.subset, "HOMA.ir_BMI.adj.partial", c ("estimate", "p.value")] =
	t (apply (mgs_abundance [ctrl.no.na.4MGS, mgs.subset], MARGIN = 2, FUN = function (x) 
	  unlist (pcor.test (x = x, y = phenotypes [ctrl.no.na.4MGS, "Homa.IR"], z = phenotypes [ctrl.no.na.4MGS, "BMI.kg.m2"],
		method = cor_method) [c ("estimate", "p.value")])))       

cor_HOMA.IR [["MGSs"]] <- tmpMat
rm (tmpMat)

###
### Step 10 - Link KEGG functions to phenotype of interest
###
### Analogous in goal to Step 6 but for metagenomics functional data. 
### Here we use KEGG modules, but any other groupings of genes into functional
### modules could similarly be used (see examples in Table 1 in the accompanying
### Protocol). However, it is more complex in that each KEGG module is
### constituted by multiple KOs. Thus, to generate results on the level of
### modules, a test is made if correlations between the phenotype and the
### abundances of KOs in the module is significantly higher or lower (MWU test)
### for the module member KOs than for all other KOs, thus also considering
### module completeness beyond the single gene level. In case of a binary
### phenotype variable, the KOs can instead (of Spearman correlation
### coefficients) be ranked based on Wald statistics for testing differentially
### abundant KOs with a negative binomial test with the DESeq2 R package using
### non-rarefied gene counts.
###

tmpMat = array (NA, c (length (koann), 3, 2))
dimnames (tmpMat) [[1]] = names (koann)
dimnames (tmpMat) [[2]] = c ("HOMA.ir", "HOMA.ir_BMI.adj.partial", "HOMA.ir_RichnessGenes7M.adj.partial")
dimnames (tmpMat) [[3]] = c ("estimate", "p.value")

### KOs sparsity filter step - exclude KOs that occur in <3 of the control
### individuals that are actually used in the analysis. Repeat the step since we
### have excluded control individuals with missing information for HOMA-IR
ctrl.no.na.4KOs = intersect (ctrl.no.na, rownames (ko_abundance)) # length=275
test = apply (ko_abundance [intersect (ctrl.no.na.4KOs, rownames (ko_abundance)),], 2, function (x) length (x [x != 0]))
ko.subset = names (test [test >= 3])
rm (test)
### Sparsity filtering done

### Associating KOs with HOMA-IR without de-confounding for BMI    
### First, obtain the correlation coefficients for Spearman correlation between
### KOs and HOMA-IR
KO_cor = apply (ko_abundance [ctrl.no.na.4KOs, ko.subset], MARGIN = 2, FUN = function (x) 
	cor (phenotypes [ctrl.no.na.4KOs, "Homa.IR"], x,
	method = cor_method, use = "pairwise.complete.obs"))      
KO_cor_HOMA = KO_cor
### Then, test for difference in correlation coefficients between KOs in the
### KEGG module and all other KOs.
for (k in names (koann)) {
	incat =    na.omit (KO_cor [   names (KO_cor) %in% koann [[k]] ]) ### select all correlations between HOMA-IR and KOs in the KEGG module
	notincat = na.omit (KO_cor [! (names (KO_cor) %in% koann [[k]])]) ### select all correlations between HOMA-IR and KOs NOT in the KEGG module
	if (length (incat) > 0 & length (notincat) > 0) {
		x = wilcox.test (incat, notincat)
		tmpMat [k, "HOMA.ir", "p.value"] = x$p.value
		tmpMat [k, "HOMA.ir", "estimate"] = (median (incat, na.rm = T) - median (notincat, na.rm = T))
	}
}
rm (KO_cor)

### Associating KOs with HOMA-IR while de-confounding for BMI    
### First, obtain the correlation coefficients for partial Spearman correlation
### between KOs and HOMA-IR
KO_cor = apply (ko_abundance [ctrl.no.na.4KOs, ko.subset], MARGIN = 2, FUN = function (x) 
	pcor.test (x = x, y = phenotypes [ctrl.no.na.4KOs, "Homa.IR"], z = phenotypes [ctrl.no.na.4KOs, "BMI.kg.m2"],
	method = cor_method) [1, "estimate"])      
KO_cor_HOMAadjBMI = KO_cor
### Then, test for difference in correlation coefficients between KOs in the 
### KEGG module and all other KOs.
for (k in names (koann)) {
	incat =    na.omit (KO_cor [   names (KO_cor) %in% koann [[k]] ]) ### select all correlations between HOMA-IR.BMI.adj and KOs in the KEGG module
	notincat = na.omit (KO_cor [! (names (KO_cor) %in% koann [[k]])]) ### select all correlations between HOMA-IR.BMI.adj and KOs NOT in the KEGG module
	if (length (incat) > 0 & length (notincat) > 0) {
		x = wilcox.test (incat, notincat)
		tmpMat [k, "HOMA.ir_BMI.adj.partial", "p.value"]  = x$p.value
		tmpMat [k, "HOMA.ir_BMI.adj.partial", "estimate"] = (median (incat, na.rm = T) - median (notincat, na.rm = T))
	}
}
rm (KO_cor)

### Associating KOs with HOMA-IR while de-confounding for gut gene richness    
### First, obtain the correlation coefficients for partial Spearman correlation
### between KOs and HOMA.IR
KO_cor = apply (ko_abundance [ctrl.no.na.4KOs, ko.subset], MARGIN = 2, FUN = function (x) 
	pcor.test (x = x, y = phenotypes [ctrl.no.na.4KOs, "Homa.IR"], z = phenotypes [ctrl.no.na.4KOs, "richnessGenes7M"],
	method = cor_method) [1, "estimate"])      
KO_cor_HOMAadjRichnessGenes7M = KO_cor
### Then, test for difference in correlation coefficients between KOs in the 
### KEGG module and all other KOs.
for (k in names (koann)) {
	incat =    na.omit (KO_cor [   names (KO_cor) %in% koann [[k]] ]) ### select all correlations between HOMA-IR.RichnessGenes.adj and KOs in the KEGG module
	notincat = na.omit (KO_cor [! (names (KO_cor) %in% koann [[k]])]) ### select all correlations between HOMA-IR.RichnessGenes.adj and KOs NOT in the KEGG module
	if (length (incat) > 0 & length (notincat) > 0) {
		x = wilcox.test (incat, notincat)
		tmpMat [k, "HOMA.ir_RichnessGenes7M.adj.partial", "p.value"]  = x$p.value
		tmpMat [k, "HOMA.ir_RichnessGenes7M.adj.partial", "estimate"] = (median (incat, na.rm = T) - median (notincat, na.rm = T))
	}
}
rm (KO_cor)

cor_HOMA.IR [["keggmodules"]] <- tmpMat
rm (tmpMat, ko.subset, ctrl.no.na.4KOs)

###
### Step 11 - Save phenotype associations
###
### Save the (BMI corrected) HOMA-IR association of metabolite clusters, MGSs
### and KEGG modules calculated in Step 8-10.
###

tmp.excel.file = "results/HOMA.IR_associations.xlsx"
write.xlsx (paste ("Associations with HOMA-IR and HOMA-IR adjusted for BMI"), 
            sheetName = "info", file = tmp.excel.file, row.names = F, col.names = F)

for (tmp_name in names (cor_HOMA.IR)) {
  tmpMat = cor_HOMA.IR [[paste (tmp_name)]]
  out = cbind (tmpMat [, "HOMA.ir", c("estimate", "p.value")], "p.adjust" = p.adjust (tmpMat [, "HOMA.ir","p.value"], method = "BH"), 
               tmpMat [, "HOMA.ir_BMI.adj.partial", c ("estimate", "p.value")], "p.adjust" = p.adjust (tmpMat [, "HOMA.ir_BMI.adj.partial", "p.value"], method = "BH"))
  colnames (out) = paste (rep (c ("HOMA.IR", "HOMA.IR_BMI.adj.partial"), each = 3), colnames (out), sep = "_")
  if (tmp_name == "metlip") { rownames (out) = cluster_mapping_file [match( rownames (out), cluster_mapping_file$New_Name), "label"] }
  write.xlsx (x = out, file = tmp.excel.file, append = T, sheetName = paste (tmp_name, sep = "_"), row.names = T, col.names = T)            
}

###
### Step 12 - Select features with significant differences
###
### Here, combine and integrate those functional, taxonomic and metabolomics
### features which reliably correspond to the host phenotype of interest, then
### later determine their inter-correlations.
###

final.fdr.cutoffs = 0.1 ### FDR thresholds, change to your likings.

vartotest_union <- list ()
vartotest_intersect <- list ()
for (tmp_name in names (cor_HOMA.IR)) {
	tmpMat = cor_HOMA.IR [[paste (tmp_name)]]
	vartotest_union [[paste ("fdr", final.fdr.cutoffs, sep = "_")]] [[tmp_name]] =
	unique (c (names (which (p.adjust (tmpMat [ , "HOMA.ir", "p.value"], method = "BH") < final.fdr.cutoffs)), 
		   names (which (p.adjust (tmpMat [ , "HOMA.ir_BMI.adj.partial", "p.value"], method = "BH") < final.fdr.cutoffs))))
	vartotest_intersect [[paste("fdr", final.fdr.cutoffs, sep = "_")]] [[tmp_name]] =
	intersect (names (which (p.adjust (tmpMat[ , "HOMA.ir", "p.value"], method = "BH") < final.fdr.cutoffs)),
		   names (which (p.adjust (tmpMat[ , "HOMA.ir_BMI.adj.partial", "p.value"], method = "BH") < final.fdr.cutoffs)))
}

###
### Step 13 - Correlate metabolite clusters to functional metagenomic potentials
### (KOs)
###
### The set of metabolite clusters associated with the phenotype of interest
### should next be tested for association with the set of functional metagenome
### features likewise so associated (identified in Step 12). Here once again the
### complex nature of KEGG modules (consisting of multiple KOs) must be taken
### into account.
###
### In this step, the correlations between each metabolite/lipid cluster and
### each KO are determined.
###

ctrl.4KOs = intersect (ctrl, rownames (ko_abundance))

KO_MetLip_cor = matrix (NA, nrow = ncol (ko_abundance), ncol = ncol (MEsMetLip))
rownames (KO_MetLip_cor) = colnames (ko_abundance)
colnames (KO_MetLip_cor) = colnames (MEsMetLip)

for (m in colnames (MEsMetLip)) {
		KO_MetLip_cor [ , m] = apply (ko_abundance [ctrl.4KOs, ], MARGIN = 2, FUN = function (x) 
		  cor (x, MEsMetLip [ctrl.4KOs, m],
		       method = cor_method, use = "pairwise.complete.obs"))
}

###
### Step 14 - Associate metabolite clusters to functional metagenomic potentials
### (KEGG modules)
###
### Using the KO level data generated in Step 13, to calculate module-level
### associations between a KEGG module and a metabolite cluster.
###

### Select variables to include in the heatmap
### note that for the metabolite clusters (metlip) we take the intersect, to 
### reduce the number of clusters in the heatmap, but for the KEGG modules 
### (keggmodules) we use the union.
metlip2test = vartotest_intersect [[paste("fdr", final.fdr.cutoffs, sep = "_")]][["metlip"]] 
keggmodules2test = vartotest_union [[paste("fdr", final.fdr.cutoffs, sep = "_")]][["keggmodules"]] 

kegg_metlip_est = matrix (data = NA, nrow = length (keggmodules2test), ncol = length (metlip2test))
rownames (kegg_metlip_est) = keggmodules2test
colnames (kegg_metlip_est) = metlip2test
kegg_metlip_p = kegg_metlip_est

for (m in metlip2test) {
	for (k in keggmodules2test) {
		incat =    na.omit (KO_MetLip_cor [   rownames (KO_MetLip_cor) %in% koann [[k]],  m]) ### select all correlations between metlip group and KOs in the KEGG module
		notincat = na.omit (KO_MetLip_cor [! (rownames (KO_MetLip_cor) %in% koann [[k]]), m]) ### select all correlations between metlip group and KOs NOT in the KEGG module
		if (length (incat) > 0 & length (notincat) > 0) {
			x = wilcox.test (incat, notincat)
			kegg_metlip_p [k, m] = x$p.value
			kegg_metlip_est [k, m] = median (incat, na.rm = T) - median (notincat, na.rm = T)
		}
	}
}
rm (KO_MetLip_cor)

###
### Step 15 - Plot metabolome-microbiome functional analysis results
###
### Create visual representation of the generated results.
###

### load modified heatmap function from the gplots R-package, which allows for
### both multi-column row-sidebar and showing text within cells (here
### significance). Note this function (heatmap.3) is different from the
### 'heatmap3' R package and the heatmap.3 function in GMD package The latest
### version can be obtained from https://gist.github.com/amcdavid/5439787 We
### have tested the script with the version loaded below.
source("r-code/heatmap_3.R")

### BH adjust p-values for associations between KEGG modules and metabolite
### modules
kegg_metlip_p_adj = p.adjust (kegg_metlip_p, method = "BH")
dim (kegg_metlip_p_adj) = dim (kegg_metlip_p)
rownames (kegg_metlip_p_adj) = rownames (kegg_metlip_p)
colnames (kegg_metlip_p_adj) = colnames (kegg_metlip_p)

### Select only significant associations to show in heatmap (to make it visually
### comprehensible) I.e. exclude KEGG modules and metabolite modules with not at
### least one significant association. In this example, all rows/columns are
### kept.
tmp = kegg_metlip_p_adj < final.fdr.cutoffs
issig = rowSums (tmp, na.rm = T)
kegg = na.omit (names (issig [issig > 0]))
length (kegg)
issig = colSums (tmp, na.rm=T)
metlip = na.omit (names (issig [issig > 0]))
length (metlip)

### create matrix for heatmap
plotmat = kegg_metlip_est [kegg, metlip]

### create matrix with significance stars
plotmat_p = kegg_metlip_p_adj [kegg, metlip] ### p-values to make stars for heatmap
stars = matrix ("", ncol = ncol (plotmat_p), nrow = nrow (plotmat_p))
rownames (stars) = rownames (plotmat_p)
colnames (stars) = colnames (plotmat_p)
for (z in 1:ncol (stars)) {
  for (j in 1:nrow (stars)) {
    if (plotmat_p [j, z] < 0.1) {
      stars [j, z] = "+"
    }
    if (plotmat_p [j, z] < 0.01) {
      stars [j, z] = "*"
    }
    if (plotmat_p [j, z] < 0.001) {
      stars [j, z] = "**"
    }
  }
}
rm (plotmat_p)

### make sidebar with phenotype associations
pheno = c ("HOMA.ir", "HOMA.ir_BMI.adj.partial")
phenobar = matrix (NA, ncol = length (pheno), nrow = length (kegg))
colnames (phenobar) = pheno
rownames (phenobar) = kegg

for (j in pheno) {
  tmp = p.adjust (cor_HOMA.IR$keggmodules [ , j, "p.value"], method="BH") ### adjust column wise on phenotypes
  for (k in kegg) {
    if (tmp [k] >= 0.1) {
      phenobar [k, j] = "grey"
    } else {
      if (cor_HOMA.IR$keggmodules [k , j, "estimate"] > 0) {
        phenobar [k, j] = "darkblue"
      } else if (cor_HOMA.IR$keggmodules [k , j, "estimate"] < 0){
        phenobar[k, j] = "darkred"
      } else {
        phenobar[k, j] = "white"
      }
    }
  }
}

### Rename column names of phenobar
colnames (phenobar) = c("HOMA-IR", "HOMA-IR.BMI.adj")

### Note, normally one would cluster the rows and columns in the heatmap as
### shown below but for the paper (Pedersen et al, 2016) we needed to group the
### KEGG modules by biological similarity to make a higher-level annotation for
### the figure and thus made a manual arrangement.

### cluster rows of matrix 
d = dist (plotmat)
dend = as.dendrogram (hclust (d, method = "average"))
rm (d)

### cluster columns of matrix 
d2 = dist (t (plotmat))
dend2 = as.dendrogram (hclust (d2, method = "average"))
rm (d2)

### make heatmap with all the selected KEGG modules
pdf ("results/heatmap_KEGG_vs_metabolite_clusters.pdf", height = 10, width = 10)
heatmap.3 (plotmat,
          Rowv = dend, 
          Colv = dend2, 
          dendrogram = "column",
          col = rev (bluered (100)),
          symbreaks = T,
          key = T,
          symkey = T,
          keysize = 1,
          KeyValueName = "SCC.bg.adj",
          lhei = c(0.6, 4),
          cellnote = stars, 
          notecol = "black",
          notecex = 1.1,
          trace = "none",
          labRow = module_mapping_clean [rownames (plotmat)],
          labCol = cluster_mapping_file [match( colnames (plotmat), cluster_mapping_file$New_Name), "label"],
          RowSideColors = t (phenobar [rownames (plotmat), ]), 
          side.height.fraction = 0.35,    
          NumColSideColors = dim (phenobar) [2], 
          margins = c (23, 33),
          cexRow = 1,
          cexCol = 1
)
dev.off () 

###
### Step 16 (optional) - Export metabolome-microbiome associations for network
### analysis
###
### Further exploration of such high dimensional association data can be
### performed using software for network analysis, e.g. Cytoscape
### (www.cytoscape.org) or igraph (www.igraph.org). Here we provide code for
### exporting an edge-file with pairwise association scores and FDR-values
### between metabolite clusters and KEGG modules as well as a corresponding node
### attribute file; both in .txt-format. Several tutorials for importing,
### visualizing and analyzing networks in Cytoscape can be found here:
### https://github.com/cytoscape/cytoscape-tutorials/wiki.
###

### Note, for simplicity, the node and edge files are made for a network
### representing the heatmap plotted in Step 15 (i.e. KEGG modules in
### "keggmodules2test" and metabolite clusters in "metlip2test")

### Make edge file with the following four columns: 
### - KEGG_module, 
### - Metabolite_cluster, 
### - association estimate (SCC.bg.adj), 
### - FDR adjusted p-values
tmp.edge.file.p_adj = melt (kegg_metlip_p_adj [keggmodules2test, metlip2test])
tmp.edge.file.est   = melt (kegg_metlip_est [keggmodules2test, metlip2test])
colnames (tmp.edge.file.p_adj) = c ("KEGG_module", "Metabolite_cluster", "p.adjust")
colnames (tmp.edge.file.est)   = c ("KEGG_module", "Metabolite_cluster", "SCC.bg.adj")
edge.file = cbind (tmp.edge.file.est, "p.adjust" = tmp.edge.file.p_adj$p.adjust)
rm (tmp.edge.file.p_adj, tmp.edge.file.est)

### Save the edge-file as a tab-seperated txt.file
write.table (edge.file, file = "results/edge_file.txt", sep = "\t", row.names = F, col.names = T, quote = F)

### Make node annotation file with the following three columns: 
### - Node_name (corresponds to names in the edge-file), 
### - Description (where existing), 
### - Examples (example metabolites, only for metabolite clusters)
### Get metabolite annotation
node.file = cluster_mapping_file [cluster_mapping_file$New_Name %in% metlip2test, c ("New_Name", "Description", "Examples")]
colnames (node.file) = c ("Node_name", "Description", "Examples")
### Add KEGG module annotation
node.file = rbind (node.file, data.frame ("Node_name" = keggmodules2test, "Description" = module_mapping_clean [keggmodules2test], "Examples" = rep ("", time = length (keggmodules2test))))

### Save the node-file as a tab-seperated txt.file
write.table (node.file, file = "results/node_file.txt", sep = "\t", row.names = F, col.names = T, quote = F)

###
### Step 17 - Leave-one-MGS-out analysis
###
### This part of the approach allows testing which bacterial taxa (sensu
### metagenomics species (MGS) as defined from the metagenomics datasets
### themselves) are driving the functional effects seen. I.e. enabling
### assessment of the extent to which different taxa explain a functional
### potential association to a phenotype of interest, such as insulin resistance
### in the case of Pedersen et al., 2016. It is done by testing for each
### functional feature (here: KEGG modules) to what extent leaving out each MGS
### and the genes it contains causes a change in the association between those
### modules and the target phenotype. Thus, the first step is setting which
### modules and taxa will be tested, then performing this test for each
### combination, after which results are plotted and reported.
###

### Specify the MGSs to be tested in the leave-on-MGS-out analysis. 
### Here it is set to all 788 MGSs.
MGS.list = MGSs
### Specify the KEGG modules to test in the leave-one-MGS-out analysis, and get
### the KOs for those modules.
### Here it is set to all modules defined as significant associated with HOMA-IR
### and/or HOMA-IR.bmi.adjusted (i.e. 'keggmodules2test').
KOsets = koann [keggmodules2test] 
### Manually add the bcaa_biosynthesis 'mega-module' (relevant only for Pedersen
### et al., 2016)
KOsets [["bcaa_biosyn"]] = unique (c (koann [["M00019"]], koann [["M00570"]], koann [["M00535"]], koann [["M00432"]]))

### The following code can take long time.
### for testing purposes, especially on older systems, we recommend subsetting 
### to a few KEGG modules and maybe also MGSs.
# KOsets = list ("M00060"= KOsets[["M00060"]],  "bcaa_biosyn"=KOsets[["bcaa_biosyn"]])
# MGS.list = c ("T2DCAG00004", "T2DCAG00005", "T2DCAG00385", " T2DCAG00011")

### Functions: 

### Calculate abundance sum of ALL genes annotated to a KO
KOprofile_FUN <- function (KO) {
	KOgenes_i <- unique (unlist (KO2gene [KO])) ### get all genes that are annotated to the KO
	return (colSums (gene_abundance_sub [KOgenes_i,]))
}

### Calculate abundance sum of genes annotated to a KO while excluding the genes
### that are part of 'LeftOutMGS'
LOOKOprofile <- function (KO, LeftOutMGS) {
	if ( any (LeftOutMGS == unlist (KO2MGS [KO])) ){	KOgenes_i <- setdiff (unique (unlist (KO2gene [KO])), MGS2gene [[LeftOutMGS]])} 
	else { KOgenes_i <- unique (unlist (KO2gene [KO])) }
	return (colSums (gene_abundance_sub [KOgenes_i,]))
}

### The main function calculating the effect of leaving-one-MGS-out on the 
### association between a KEGG module and a phenotype of interest
LeaveOneMGSOut_FUN <- function (KOsets = KOsets, MGS.list = MGS.list, sample.list = sample.list, y, z = NULL) {
	### number of KOs kinds per MGS
	KO_types_per_MGS <- lapply (KOsets, function (KOset) {
		mgses <- sapply (KOset, function (KO) { intersect (unique (unlist (KO2MGS [KO])), MGS.list)})
		table (unlist (mgses))})
  
	if  (is.null (z) == T) { 
    
		print ("No 'z' is provided, calculating spearman correlation")
    
		### Code for spearman correlation, i.e. NOT adjusting for a confounding factor
		y = y [sample.list,]
    
		### correlating the mean KO signal (calculated with 'KOprofile_FUN') to 
		### HOMA-IR, and then taking the median SCC for the module
		print ("calculating spearman correlation, using all genes")
		SCCallMGS<-sapply (KOsets, function (KOset) { 
      			KOprofiles<-t (sapply (KOset, KOprofile_FUN))
			return (median (apply (KOprofiles [rowSums (KOprofiles) > 0, sample.list], 1, function (x) { 
				cor.test (x = x, y = y, method = "spearman", use = "pairwise.complete.obs", exact = FALSE)$estimate }))) 
		})
    
		### correlating the mean KO signal while excluding contribution from 'MGS' 
		### (calculated with 'LOOKOprofile') to HOMA-IR, and then taking the median 
		### SCC for the module.
		### Repeating for all MGS in 'MGS.list'
		print ("calculating spearman correlation, leaving-one-MGS-out")
		SCComitingMGS <- lapply (KOsets, function (KOset) {
			sapply (intersect (unique (unlist (KO2MGS[KOset])),MGS.list), function (MGS){ 
				KOprofiles <- (sapply (KOset, function (KO) { LOOKOprofile (KO, MGS) }))
				return (median (cor (x = KOprofiles [sample.list, colSums (KOprofiles) > 0], y = y, method = "spearman", use = "pairwise.complete.obs"))) }) 
		})
    
		### Summarizing outout
		DeltaSCCperMGS <- lapply (names (SCComitingMGS), function (N) {
			SCC = SCCallMGS [[N]]
			SCC.bgadj =  (median (na.omit (KO_cor_HOMA [names (KO_cor_HOMA) %in% KOsets [[N]]]), na.rm = T) 
				          - median (na.omit (KO_cor_HOMA[! (names (KO_cor_HOMA) %in% KOsets [[N]])]), na.rm = T)) 
			# summary(KO_cor_HOMA); adjust % SCC effect for background distribution (i.e. the fact that median(SCC for HOMA vs modules) is negative and not = 0)
			data.frame (SCC = SCC, SCC.bgadj = SCC.bgadj, SCC_omiting_MGS = SCComitingMGS [[N]], DeltaMGS_SCC = SCC - SCComitingMGS [[N]], 
				pctSCCeffect = 100 * (SCC - SCComitingMGS [[N]]) / SCC, pctSCCeffect.bgadj = 100 * (SCC - SCComitingMGS [[N]]) / SCC.bgadj,
				Distinct_KOs_in_MGS = as.vector (KO_types_per_MGS [[N]][names (SCComitingMGS [[N]])]), row.names = names (SCComitingMGS [[N]])
			) [rev (order (pctSCCeffect = 100 * (SCC - SCComitingMGS[[N]]) / SCC)),] 
		})

		names (DeltaSCCperMGS) <- names (KOsets)
		return (DeltaSCCperMGS)

	} 

	else {
    
		print ("'z' is provided, calculating partial spearman correlation, adjusting for 'z'")

		### Code for partial spearman correlation, i.e. adjusting for a confounding factor
		y = y [sample.list,]
		z = z [sample.list,]

		### correlating the mean KO signal (calculated with 'KOprofile_FUN') to 
		### HOMA-IR, and then taking the median SCC for the module
		print ("calculating partial spearman correlation, using all genes")
		partialSCCallMGS <- sapply (KOsets, function (KOset) {	
			KOprofiles <- t (sapply (KOset, KOprofile_FUN))
			return (median (apply (KOprofiles [rowSums (KOprofiles) > 0,], 1, function (x) { pcor.test (x = x [sample.list], y = y, z = z, method = "spearman")$estimate })))
		})

		### correlating the mean KO signal while excluding contribution from 'MGS' 
		### (calculated with 'LOOKOprofile') to HOMA-IR, and then taking the median 
		### SCC for the module.
		### Repeating for all MGS in 'MGS.list'
		print ("calculating partial spearman correlation, leaving-one-MGS-out")
		partialSCComitingMGS <- lapply (KOsets, function (KOset) {
			sapply (intersect (unique (unlist (KO2MGS [KOset])), MGS.list), function (MGS) { # looping over all MGGs in 'MGS.list'
				KOprofiles <- t (sapply (KOset, function (KO) {LOOKOprofile (KO, MGS) }))
				return (median (apply (KOprofiles [rowSums (KOprofiles) > 0,], 1, function (x) { pcor.test (x = x [sample.list], y = y, z = z, method="spearman")$estimate })))
			})
		})
    
		### Summarizing outout
		partialDeltaSCCperMGS<-lapply (names (partialSCComitingMGS), function (N){
			SCC = partialSCCallMGS[[N]]
			SCC.bgadj =  (median (na.omit (KO_cor_HOMAadjBMI[names (KO_cor_HOMAadjBMI) %in% KOsets[[N]]]), na.rm=T) 
				          - median (na.omit (KO_cor_HOMAadjBMI[! (names (KO_cor_HOMAadjBMI) %in% KOsets[[N]])]), na.rm=T)) 
			### summary(KO_cor_HOMA); adjust % SCC effect for background distribution (i.e. the fact that median(SCC for HOMA vs modules) is negative and not = 0) 
			data.frame (SCC =SCC, SCC.bgadj = SCC.bgadj, SCC_omiting_MGS = partialSCComitingMGS [[N]], DeltaMGS_SCC = partialSCCallMGS [[N]] - partialSCComitingMGS [[N]], 
				pctSCCeffect = 100 * (SCC - partialSCComitingMGS [[N]]) / SCC, pctSCCeffect.bgadj = 100 * (SCC - partialSCComitingMGS [[N]]) / SCC.bgadj,
				Distinct_KOs_in_MGS = as.vector (KO_types_per_MGS [[N]][names (partialSCComitingMGS [[N]])]), 
				row.names = names (partialSCComitingMGS [[N]])
			) [rev (order (pctSCCeffect = 100 * (SCC - partialSCComitingMGS [[N]]) / SCC)),] 
		})

		names (partialDeltaSCCperMGS) <- names (KOsets)
		return (partialDeltaSCCperMGS)

	}  

}
### End functions

### Next, run these functions to compute the delta SCC values for the test in
### question - change here if partial correlation to account for a confounder is
### desired.

### Performing the Leave-one-MGS-out for HOMA-IR
DeltaSCCperMGS <- LeaveOneMGSOut_FUN (KOsets = KOsets, MGS.list = MGS.list,
	sample.list = ctrl.no.na.4MGS, ### the subset of individuals with complete information for phenotype and metagenomic data.
	y = phenotypes [,"Homa.IR", drop = F]) 
save (DeltaSCCperMGS, file = "results/delta_SCC_per_MGS.RData")

### Performing the Leave-one-MGS-out for HOMA-IR.bmi.adjusted
# partialDeltaSCCperMGS = LeaveOneMGSOut_FUN (KOsets = KOsets, MGS.list = MGS.list,
#	sample.list = ctrl.no.na.4MGS, y = phenotypes [,"Homa.IR", drop = F], z = phenotypes [,"BMI.kg.m2", drop = F])
# save (partialDeltaSCCperMGS, file = "results/partial_delta_SCC_per_MGS.RData")

###
### Step 18 - Extracting top driver species from leave-one-MGS-out analysis for
### each microbiome functional module
###
### Here, having computed the contribution per taxon to each module, we extract
### top five drivers for interpretation.
###

topX = 5
output = as.data.frame (matrix (NA, nrow = length (DeltaSCCperMGS), ncol= (2 + 4 * topX)))
colnames (output) = c ("Module name","Number of genes in KEGG_module",
			"top1", "top1_Number of module genes in MGS", "top1_DeltaMGS_SCC",  "top1_pctSCCeffect_bg.adj",
			"top2", "top2_Number of module genes in MGS", "top2_DeltaMGS_SCC",  "top2_pctSCCeffect_bg.adj",
			"top3", "top3_Number of module genes in MGS", "top3_DeltaMGS_SCC",  "top3_pctSCCeffect_bg.adj",
			"top4", "top4_Number of module genes in MGS", "top4_DeltaMGS_SCC",  "top4_pctSCCeffect_bg.adj",
			"top5", "top5_Number of module genes in MGS", "top5_DeltaMGS_SCC",  "top5_pctSCCeffect_bg.adj")
rownames (output) = names (DeltaSCCperMGS)

output [,"Number of genes in KEGG_module"] = sapply (rownames (output), function (i) length (KOsets [[i]]))
output [,"Module name"] = module_mapping [rownames (output)]

for (i in names (DeltaSCCperMGS)) {

	if (length (rownames (DeltaSCCperMGS[[i]])) <= topX ) {

		tmp.MGSs.top = rownames (DeltaSCCperMGS [[i]]) [1:topX]
		tmp.MGSs.top [DeltaSCCperMGS [[i]][tmp.MGSs.top, "pctSCCeffect"] <= 0] <- NA ### set top.MGSs to NA if removing them are NOT INCREASING the effect
		tmp.MGSs = c (tmp.MGSs.top)

	} 

	else {

		tmp.MGSs.top = rownames (DeltaSCCperMGS [[i]]) [1:topX]
		tmp.MGSs.top [DeltaSCCperMGS [[i]][tmp.MGSs.top, "pctSCCeffect"] <= 0] <- NA ### set top.MGSs to NA if removing them are NOT INCREASING the effect
		tmp.MGSs = c (tmp.MGSs.top)

	}

	## Add taxonomically annotation (for those MGSs that have one)
	tmp.MGSs.names = sapply (tmp.MGSs, function (i) ifelse (i %in% rownames (mgs_taxonomy),
		paste (" ", i, " : ", mgs_taxonomy [i, "genus"], " sp ", sep = ""), # true
		paste (" ", i, " ", sep = "") )) # false
	tmp.KOsInMGS      = DeltaSCCperMGS [[i]][ tmp.MGSs, "Distinct_KOs_in_MGS"]
	tmp.delta         = DeltaSCCperMGS [[i]][ tmp.MGSs, "DeltaMGS_SCC"]
	tmp.pct_bg.adj    = DeltaSCCperMGS [[i]][ tmp.MGSs, "pctSCCeffect.bgadj"]
	output [i, 3:ncol (output)] = c (rbind (tmp.MGSs.names, tmp.KOsInMGS, tmp.delta, tmp.pct_bg.adj))

	rm (tmp.delta, tmp.pct_bg.adj, tmp.MGSs.names, tmp.KOsInMGS, tmp.MGSs, tmp.MGSs.top)

}

### write 'output' to tab-delimited text file
write.table (x = output, file = "results/top_driver_species.txt", row.names = T, col.names = NA, quote = F, sep = "\t")         

###
### Step 19 - Plotting leave-one-MGS-out results
###
### Having computed these results, we plot the top driver taxa for the
### phenotype-associated gene functions (Figure 4).
### These include sub-plots for:
### - Distribution/density plots of correlations for KOs in a KEGG module vs 
###   all other KOs
### - Distribution of correlations when leave-one-MGS-out. 
### - Distribution of correlations when leave-one-MGS-out. bg.adj. SCC 
###   (i.e. what was shown in Figure 3c+d in Pedersen et al., 2016)
### Median SCC for KO within a module (red) and all other remaining KOs (green)
### are indicated in the first two plots.
###

### Specify pdf-file for output.
pdf (file = "results/density_plot_SCC_HOMA.IR.pdf", width = 6 * 1, height = 3 * 3)

for (m in names (KOsets)) {
  
  ### Distribution of correlations for KOs in modules vs all other KOs
  incat =     na.omit (KO_cor_HOMA [names (KO_cor_HOMA) %in% KOsets [[m]] ]) ### select all correlations between HOMA-IR and KOs in the KEGG module
  incat2 =    data.frame ("KO_cor_HOMA" = incat, "cat" = "KOs in module")
  notincat =  na.omit (KO_cor_HOMA [! (names (KO_cor_HOMA) %in% KOsets [[m]])]) ### select all correlations between HOMA-IR and KOs NOT in the KEGG module   
  notincat2 = data.frame ("KO_cor_HOMA" = notincat, "cat" = "KOs not in module")
  tmp.df = rbind (incat2, notincat2)
  cdf <- ddply (tmp.df, "cat", summarise, rating.median = median (KO_cor_HOMA)) ### find median for each group
  g1 <- ggplot (tmp.df, aes (x = KO_cor_HOMA, fill = cat)) + 
    geom_density (alpha = 0.4) +
    geom_vline (data = cdf, aes (xintercept = rating.median, colour = cat), linetype = "dashed", size = 0.8) +
    ggtitle (paste0 ("KEGG module: ", m, 
                     "\n", strsplit (module_mapping [m], split = "\\[|,")[[1]][1], 
                     "\n", length (KOsets [[m]]), " KOs in module vs all remaining ", (length (KO_cor_HOMA) - length (KOsets[[m]])), " KOs", "\n")) +
    xlab ("SCC for KOs and HOMA-IR") +
    theme_classic () + theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ())
  
  ### Plot Distribution of correlations when leave-one-MGS-out.   
  if (nrow (DeltaSCCperMGS [[m]]) > 1) {
    g2 <- ggplot (DeltaSCCperMGS [[m]], aes (x = SCC_omiting_MGS)) + 
      geom_density (alpha = 1, fill = "grey", aes (y = ..scaled..)) + 
      geom_segment (aes (y = -0.1, yend = -0.02, x = SCC_omiting_MGS, xend = SCC_omiting_MGS)) +
      geom_vline (data = cdf, aes (xintercept = rating.median, colour = cat), linetype = "dashed", size = 0.8) +
      ggtitle (paste0 ("KEGG module: ", m, 
                       "\n", strsplit (module_mapping [m], split = "\\[|,")[[1]][1], 
                       "\nSCC for leave-1-MGS-out",
                       "\nNumber of MGSs = ", nrow (DeltaSCCperMGS [[m]]))) +
      xlab ("SCC for KOs and HOMA-IR") +
      theme_classic () + theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ())
    
  } else {
    
    ### Special circumstance when there is only one MGS (then one cannot make a 
    ### density plot, instead a black line is plotted at the value when that 
    ### one MGS is left out).
    
    g2 <- ggplot () +
      scale_x_continuous (limits = range (c (DeltaSCCperMGS [[m]][, "SCC_omiting_MGS"], cdf$rating.median))) +
      scale_y_continuous (name = "", limits = c (0, 1)) +
      geom_vline (data = cdf, aes (xintercept = rating.median, colour = cat), linetype = "dashed", size = 0.8) +
      geom_vline (data = DeltaSCCperMGS [[m]], aes (xintercept = SCC_omiting_MGS), linetype = "longdash", color = "black",size = 1) + 
      ggtitle (paste0 ("KEGG module: ", m, 
                       "\n", strsplit (module_mapping [m], split = "\\[|,") [[1]][1], 
                       "\nSCC for leave-1-MGS-out",
                       "\nNumber of MGSs = ", nrow (DeltaSCCperMGS [[m]]))) +
      xlab ("SCC for KOs and HOMA-IR") +
      theme_classic () + theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ()) 
    
  }
  
  ### Plot distribution of correlations when leave-one-MGS-out. Bg.adjust SCC 
  ### (i.e. what we show in Figure 3c+d)
  tmp.DeltaSCCperMGS = DeltaSCCperMGS [[m]]
  ### Calculate delta-SCC in relation to the background adjusted SCC. 
  ### i.e the 'plottedvalue.s.m' shown as the second last equation in the 
  ### methods section in Pedersen et al., 2016.
  tmp.DeltaSCCperMGS$DeltaMGS_SCC.bgadj = tmp.DeltaSCCperMGS$SCC.bgadj - tmp.DeltaSCCperMGS$DeltaMGS_SCC
  ### specify the range of the x-axis to include 0, i.e. [min:0] for negative 
  ### correlations and [0:max] for positive correlations
  x.range = c (min (0, tmp.DeltaSCCperMGS$DeltaMGS_SCC.bgadj), max (0, tmp.DeltaSCCperMGS$DeltaMGS_SCC.bgadj))
  
  if (nrow (tmp.DeltaSCCperMGS) > 1 ) {
    
    g3 <- ggplot (tmp.DeltaSCCperMGS, aes (x = DeltaMGS_SCC.bgadj)) + 
      geom_density (alpha = 1, fill = "grey", aes (y = ..scaled..)) + 
      geom_segment (aes (y = -0.1, yend = -0.02, x = DeltaMGS_SCC.bgadj, xend = DeltaMGS_SCC.bgadj)) +
      ggtitle (paste0 ("KEGG module: ", m, 
                       "\n", strsplit (module_mapping [m], split = "\\[|,")[[1]][1], 
                       "\nbg.adj.SCC for leave-1-MGS-out",
                       "\nNumber of MGSs = ", nrow (tmp.DeltaSCCperMGS))) +
      xlab ("bg.adj.SCC for KOs and HOMA-IR") + xlim (x.range) +
      theme_classic () + theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ())
    
  } else { 
    
    ### Special circumstance when there is only one MGS (then one cannot make a 
    ### density plot)
    
    g3 <- ggplot () +
      ggtitle (paste0 (m, 
                       "\n", strsplit (module_mapping [m], split="\\[|,")[[1]][1], 
                       "\nbg.adj.SCC for leave-1-MGS-out",
                       "\nNumber of MGSs = ", nrow (tmp.DeltaSCCperMGS), " - consequently not showing a density/rug plot")) 
    
  }
  
  g = plot_grid (g1, g2, g3, ncol = 1, nrow = 3, align = "v", labels = c ("a", "b", "c"), axis = "rl")
  # g = plot_grid (g1, g3, ncol = 1, nrow = 2, align = "v", labels = c ("a", "b", "c"), axis = "rl") ### use this if deleting g2
  plot (g)
  
}

dev.off ()  

### END

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
