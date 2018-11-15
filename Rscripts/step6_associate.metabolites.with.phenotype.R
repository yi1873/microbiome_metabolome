###
### This is step6_associate.metabolites.with.phenotype.R, part of the Nature Protocol: 
### 'A computational framework for conducting a three-pronged association study
### to integrate high-throughput ‘-omics’ datasets for the identification of
### potential mechanistic links'.
### 26/7-2018
###

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
