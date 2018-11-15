###
### This is step14_associate.metabolite.clusters.with.KEGG.modules.R, part of the Nature Protocol: 
### 'A computational framework for conducting a three-pronged association study
### to integrate high-throughput ‘-omics’ datasets for the identification of
### potential mechanistic links'.
### 26/7-2018
###

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
