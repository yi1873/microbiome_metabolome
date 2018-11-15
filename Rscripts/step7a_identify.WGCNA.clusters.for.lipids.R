###
### This is step7a_identify.WGCNA.clusters.for.lipids.R, part of the Nature Protocol: 
### 'A computational framework for conducting a three-pronged association study
### to integrate high-throughput ‘-omics’ datasets for the identification of
### potential mechanistic links'.
### 26/7-2018
###

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
