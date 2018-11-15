###
### This is step13_associate.metabolite.clusters.with.KOs.R, part of the Nature Protocol: 
### 'A computational framework for conducting a three-pronged association study
### to integrate high-throughput ‘-omics’ datasets for the identification of
### potential mechanistic links'.
### 26/7-2018
###

###
### Step 13 - Correlate metabolite clusters to functional metagenomic potentials
### (KOs)
###
### The set of metabolite clusters associated with the phenotype of interest
### should next be tested for association with the set of functional metagenome
### features likewise so associated (identified in Step 12). Here once again the
### complex nature of KEGG modules (consisting of multiple KOs) must be taken
### into account.
###
### In this step, the correlations between each metabolite/lipid cluster and
### each KO are determined.
###

ctrl.4KOs = intersect (ctrl, rownames (ko_abundance))

KO_MetLip_cor = matrix (NA, nrow = ncol (ko_abundance), ncol = ncol (MEsMetLip))
rownames (KO_MetLip_cor) = colnames (ko_abundance)
colnames (KO_MetLip_cor) = colnames (MEsMetLip)

for (m in colnames (MEsMetLip)) {
  KO_MetLip_cor [ , m] = apply (ko_abundance [ctrl.4KOs, ], MARGIN = 2, FUN = function (x) 
    cor (x, MEsMetLip [ctrl.4KOs, m],
         method = cor_method, use = "pairwise.complete.obs"))
}
