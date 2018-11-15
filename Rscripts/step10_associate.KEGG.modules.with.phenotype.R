###
### This is step10_associate.KEGG.modules.with.phenotype.R, part of the Nature Protocol: 
### 'A computational framework for conducting a three-pronged association study
### to integrate high-throughput ‘-omics’ datasets for the identification of
### potential mechanistic links'.
### 26/7-2018
###

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
