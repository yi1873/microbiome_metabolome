###
### This is step18_extract.top.driver.species.R, part of the Nature Protocol: 
### 'A computational framework for conducting a three-pronged association study
### to integrate high-throughput ‘-omics’ datasets for the identification of
### potential mechanistic links'.
### 26/7-2018
###

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
