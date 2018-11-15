###
### This is step4_preprocessing.data.for.sparsity.R, part of the Nature Protocol: 
### 'A computational framework for conducting a three-pronged association study
### to integrate high-throughput ‘-omics’ datasets for the identification of
### potential mechanistic links'.
### 26/7-2018
###

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
