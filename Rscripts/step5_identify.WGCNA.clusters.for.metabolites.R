###
### This is step5_identify.WGCNA.clusters.for.metabolites.R, part of the Nature Protocol: 
### 'A computational framework for conducting a three-pronged association study
### to integrate high-throughput ‘-omics’ datasets for the identification of
### potential mechanistic links'.
### 26/7-2018
###

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
