###
### This is , part of the Nature Protocol: 
### 'A computational framework for conducting a three-pronged association study
### to integrate high-throughput ‘-omics’ datasets for the identification of
### potential mechanistic links'.
### 26/7-2018
###

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
