###
### This is step17_leave.one.MGS.out.analysis.R, part of the Nature Protocol: 
### 'A computational framework for conducting a three-pronged association study
### to integrate high-throughput ‘-omics’ datasets for the identification of
### potential mechanistic links'.
### 26/7-2018
###

###
### Step 17 - Leave-one-MGS-out analysis
###
### This part of the approach allows testing which bacterial taxa (sensu
### metagenomics species (MGS) as defined from the metagenomics datasets
### themselves) are driving the functional effects seen. I.e. enabling
### assessment of the extent to which different taxa explain a functional
### potential association to a phenotype of interest, such as insulin resistance
### in the case of Pedersen et al., 2016. It is done by testing for each
### functional feature (here: KEGG modules) to what extent leaving out each MGS
### and the genes it contains causes a change in the association between those
### modules and the target phenotype. Thus, the first step is setting which
### modules and taxa will be tested, then performing this test for each
### combination, after which results are plotted and reported.
###

### Specify the MGSs to be tested in the leave-on-MGS-out analysis. 
### Here it is set to all 788 MGSs.
MGS.list = MGSs
### Specify the KEGG modules to test in the leave-one-MGS-out analysis, and get
### the KOs for those modules.
### Here it is set to all modules defined as significant associated with HOMA-IR
### and/or HOMA-IR.bmi.adjusted (i.e. 'keggmodules2test').
KOsets = koann [keggmodules2test] 
### Manually add the bcaa_biosynthesis 'mega-module' (relevant only for Pedersen
### et al., 2016)
KOsets [["bcaa_biosyn"]] = unique (c (koann [["M00019"]], koann [["M00570"]], koann [["M00535"]], koann [["M00432"]]))

### The following code can take long time.
### for testing purposes, especially on older systems, we recommend subsetting 
### to a few KEGG modules and maybe also MGSs.
# KOsets = list ("M00060"= KOsets[["M00060"]],  "bcaa_biosyn"=KOsets[["bcaa_biosyn"]])
# MGS.list = c ("T2DCAG00004", "T2DCAG00005", "T2DCAG00385", " T2DCAG00011")

### Functions: 

### Calculate abundance sum of ALL genes annotated to a KO
KOprofile_FUN <- function (KO) {
  KOgenes_i <- unique (unlist (KO2gene [KO])) ### get all genes that are annotated to the KO
  return (colSums (gene_abundance_sub [KOgenes_i,]))
}

### Calculate abundance sum of genes annotated to a KO while excluding the genes
### that are part of 'LeftOutMGS'
LOOKOprofile <- function (KO, LeftOutMGS) {
  if ( any (LeftOutMGS == unlist (KO2MGS [KO])) ){	KOgenes_i <- setdiff (unique (unlist (KO2gene [KO])), MGS2gene [[LeftOutMGS]])} 
  else { KOgenes_i <- unique (unlist (KO2gene [KO])) }
  return (colSums (gene_abundance_sub [KOgenes_i,]))
}

### The main function calculating the effect of leaving-one-MGS-out on the 
### association between a KEGG module and a phenotype of interest
LeaveOneMGSOut_FUN <- function (KOsets = KOsets, MGS.list = MGS.list, sample.list = sample.list, y, z = NULL) {
  ### number of KOs kinds per MGS
  KO_types_per_MGS <- lapply (KOsets, function (KOset) {
    mgses <- sapply (KOset, function (KO) { intersect (unique (unlist (KO2MGS [KO])), MGS.list)})
    table (unlist (mgses))})
  
  if  (is.null (z) == T) { 
    
    print ("No 'z' is provided, calculating spearman correlation")
    
    ### Code for spearman correlation, i.e. NOT adjusting for a confounding factor
    y = y [sample.list,]
    
    ### correlating the mean KO signal (calculated with 'KOprofile_FUN') to 
    ### HOMA-IR, and then taking the median SCC for the module
    print ("calculating spearman correlation, using all genes")
    SCCallMGS<-sapply (KOsets, function (KOset) { 
      KOprofiles<-t (sapply (KOset, KOprofile_FUN))
      return (median (apply (KOprofiles [rowSums (KOprofiles) > 0, sample.list], 1, function (x) { 
        cor.test (x = x, y = y, method = "spearman", use = "pairwise.complete.obs", exact = FALSE)$estimate }))) 
    })
    
    ### correlating the mean KO signal while excluding contribution from 'MGS' 
    ### (calculated with 'LOOKOprofile') to HOMA-IR, and then taking the median 
    ### SCC for the module.
    ### Repeating for all MGS in 'MGS.list'
    print ("calculating spearman correlation, leaving-one-MGS-out")
    SCComitingMGS <- lapply (KOsets, function (KOset) {
      sapply (intersect (unique (unlist (KO2MGS[KOset])),MGS.list), function (MGS){ 
        KOprofiles <- (sapply (KOset, function (KO) { LOOKOprofile (KO, MGS) }))
        return (median (cor (x = KOprofiles [sample.list, colSums (KOprofiles) > 0], y = y, method = "spearman", use = "pairwise.complete.obs"))) }) 
    })
    
    ### Summarizing outout
    DeltaSCCperMGS <- lapply (names (SCComitingMGS), function (N) {
      SCC = SCCallMGS [[N]]
      SCC.bgadj =  (median (na.omit (KO_cor_HOMA [names (KO_cor_HOMA) %in% KOsets [[N]]]), na.rm = T) 
                    - median (na.omit (KO_cor_HOMA[! (names (KO_cor_HOMA) %in% KOsets [[N]])]), na.rm = T)) 
      # summary(KO_cor_HOMA); adjust % SCC effect for background distribution (i.e. the fact that median(SCC for HOMA vs modules) is negative and not = 0)
      data.frame (SCC = SCC, SCC.bgadj = SCC.bgadj, SCC_omiting_MGS = SCComitingMGS [[N]], DeltaMGS_SCC = SCC - SCComitingMGS [[N]], 
                  pctSCCeffect = 100 * (SCC - SCComitingMGS [[N]]) / SCC, pctSCCeffect.bgadj = 100 * (SCC - SCComitingMGS [[N]]) / SCC.bgadj,
                  Distinct_KOs_in_MGS = as.vector (KO_types_per_MGS [[N]][names (SCComitingMGS [[N]])]), row.names = names (SCComitingMGS [[N]])
      ) [rev (order (pctSCCeffect = 100 * (SCC - SCComitingMGS[[N]]) / SCC)),] 
    })
    
    names (DeltaSCCperMGS) <- names (KOsets)
    return (DeltaSCCperMGS)
    
  } 
  
  else {
    
    print ("'z' is provided, calculating partial spearman correlation, adjusting for 'z'")
    
    ### Code for partial spearman correlation, i.e. adjusting for a confounding factor
    y = y [sample.list,]
    z = z [sample.list,]
    
    ### correlating the mean KO signal (calculated with 'KOprofile_FUN') to 
    ### HOMA-IR, and then taking the median SCC for the module
    print ("calculating partial spearman correlation, using all genes")
    partialSCCallMGS <- sapply (KOsets, function (KOset) {	
      KOprofiles <- t (sapply (KOset, KOprofile_FUN))
      return (median (apply (KOprofiles [rowSums (KOprofiles) > 0,], 1, function (x) { pcor.test (x = x [sample.list], y = y, z = z, method = "spearman")$estimate })))
    })
    
    ### correlating the mean KO signal while excluding contribution from 'MGS' 
    ### (calculated with 'LOOKOprofile') to HOMA-IR, and then taking the median 
    ### SCC for the module.
    ### Repeating for all MGS in 'MGS.list'
    print ("calculating partial spearman correlation, leaving-one-MGS-out")
    partialSCComitingMGS <- lapply (KOsets, function (KOset) {
      sapply (intersect (unique (unlist (KO2MGS [KOset])), MGS.list), function (MGS) { # looping over all MGGs in 'MGS.list'
        KOprofiles <- t (sapply (KOset, function (KO) {LOOKOprofile (KO, MGS) }))
        return (median (apply (KOprofiles [rowSums (KOprofiles) > 0,], 1, function (x) { pcor.test (x = x [sample.list], y = y, z = z, method="spearman")$estimate })))
      })
    })
    
    ### Summarizing outout
    partialDeltaSCCperMGS<-lapply (names (partialSCComitingMGS), function (N){
      SCC = partialSCCallMGS[[N]]
      SCC.bgadj =  (median (na.omit (KO_cor_HOMAadjBMI[names (KO_cor_HOMAadjBMI) %in% KOsets[[N]]]), na.rm=T) 
                    - median (na.omit (KO_cor_HOMAadjBMI[! (names (KO_cor_HOMAadjBMI) %in% KOsets[[N]])]), na.rm=T)) 
      ### summary(KO_cor_HOMA); adjust % SCC effect for background distribution (i.e. the fact that median(SCC for HOMA vs modules) is negative and not = 0) 
      data.frame (SCC =SCC, SCC.bgadj = SCC.bgadj, SCC_omiting_MGS = partialSCComitingMGS [[N]], DeltaMGS_SCC = partialSCCallMGS [[N]] - partialSCComitingMGS [[N]], 
                  pctSCCeffect = 100 * (SCC - partialSCComitingMGS [[N]]) / SCC, pctSCCeffect.bgadj = 100 * (SCC - partialSCComitingMGS [[N]]) / SCC.bgadj,
                  Distinct_KOs_in_MGS = as.vector (KO_types_per_MGS [[N]][names (partialSCComitingMGS [[N]])]), 
                  row.names = names (partialSCComitingMGS [[N]])
      ) [rev (order (pctSCCeffect = 100 * (SCC - partialSCComitingMGS [[N]]) / SCC)),] 
    })
    
    names (partialDeltaSCCperMGS) <- names (KOsets)
    return (partialDeltaSCCperMGS)
    
  }  
  
}
### End functions

### Next, run these functions to compute the delta SCC values for the test in
### question - change here if partial correlation to account for a confounder is
### desired.

### Performing the Leave-one-MGS-out for HOMA-IR
DeltaSCCperMGS <- LeaveOneMGSOut_FUN (KOsets = KOsets, MGS.list = MGS.list,
                                      sample.list = ctrl.no.na.4MGS, ### the subset of individuals with complete information for phenotype and metagenomic data.
                                      y = phenotypes [,"Homa.IR", drop = F]) 
save (DeltaSCCperMGS, file = "results/delta_SCC_per_MGS.RData")

### Performing the Leave-one-MGS-out for HOMA-IR.bmi.adjusted
# partialDeltaSCCperMGS = LeaveOneMGSOut_FUN (KOsets = KOsets, MGS.list = MGS.list,
#	sample.list = ctrl.no.na.4MGS, y = phenotypes [,"Homa.IR", drop = F], z = phenotypes [,"BMI.kg.m2", drop = F])
# save (partialDeltaSCCperMGS, file = "results/partial_delta_SCC_per_MGS.RData")
