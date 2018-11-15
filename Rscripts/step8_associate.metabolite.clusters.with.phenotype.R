###
### This is step8_associate.metabolite.clusters.with.phenotype.R, part of the Nature Protocol: 
### 'A computational framework for conducting a three-pronged association study
### to integrate high-throughput ‘-omics’ datasets for the identification of
### potential mechanistic links'.
### 26/7-2018
###

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
