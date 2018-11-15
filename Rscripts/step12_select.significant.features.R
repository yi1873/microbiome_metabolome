###
### This is step12_select.significant.features.R, part of the Nature Protocol: 
### 'A computational framework for conducting a three-pronged association study
### to integrate high-throughput ‘-omics’ datasets for the identification of
### potential mechanistic links'.
### 26/7-2018
###

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
