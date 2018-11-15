###
### This is step15_plot.metabolome.microbiome.functional.analysis.R, part of the Nature Protocol: 
### 'A computational framework for conducting a three-pronged association study
### to integrate high-throughput ‘-omics’ datasets for the identification of
### potential mechanistic links'.
### 26/7-2018
###

###
### Step 15 - Plot metabolome-microbiome functional analysis results
###
### Create visual representation of the generated results.
###

### load modified heatmap function from the gplots R-package, which allows for
### both multi-column row-sidebar and showing text within cells (here
### significance). Note this function (heatmap.3) is different from the
### 'heatmap3' R package and the heatmap.3 function in GMD package The latest
### version can be obtained from https://gist.github.com/amcdavid/5439787 We
### have tested the script with the version loaded below.
source("Rscripts/heatmap_3.R")

### BH adjust p-values for associations between KEGG modules and metabolite
### modules
kegg_metlip_p_adj = p.adjust (kegg_metlip_p, method = "BH")
dim (kegg_metlip_p_adj) = dim (kegg_metlip_p)
rownames (kegg_metlip_p_adj) = rownames (kegg_metlip_p)
colnames (kegg_metlip_p_adj) = colnames (kegg_metlip_p)

### Select only significant associations to show in heatmap (to make it visually
### comprehensible) I.e. exclude KEGG modules and metabolite modules with not at
### least one significant association. In this example, all rows/columns are
### kept.
tmp = kegg_metlip_p_adj < final.fdr.cutoffs
issig = rowSums (tmp, na.rm = T)
kegg = na.omit (names (issig [issig > 0]))
length (kegg)
issig = colSums (tmp, na.rm=T)
metlip = na.omit (names (issig [issig > 0]))
length (metlip)

### create matrix for heatmap
plotmat = kegg_metlip_est [kegg, metlip]

### create matrix with significance stars
plotmat_p = kegg_metlip_p_adj [kegg, metlip] ### p-values to make stars for heatmap
stars = matrix ("", ncol = ncol (plotmat_p), nrow = nrow (plotmat_p))
rownames (stars) = rownames (plotmat_p)
colnames (stars) = colnames (plotmat_p)
for (z in 1:ncol (stars)) {
  for (j in 1:nrow (stars)) {
    if (plotmat_p [j, z] < 0.1) {
      stars [j, z] = "+"
    }
    if (plotmat_p [j, z] < 0.01) {
      stars [j, z] = "*"
    }
    if (plotmat_p [j, z] < 0.001) {
      stars [j, z] = "**"
    }
  }
}
rm (plotmat_p)

### make sidebar with phenotype associations
pheno = c ("HOMA.ir", "HOMA.ir_BMI.adj.partial")
phenobar = matrix (NA, ncol = length (pheno), nrow = length (kegg))
colnames (phenobar) = pheno
rownames (phenobar) = kegg

for (j in pheno) {
  tmp = p.adjust (cor_HOMA.IR$keggmodules [ , j, "p.value"], method="BH") ### adjust column wise on phenotypes
  for (k in kegg) {
    if (tmp [k] >= 0.1) {
      phenobar [k, j] = "grey"
    } else {
      if (cor_HOMA.IR$keggmodules [k , j, "estimate"] > 0) {
        phenobar [k, j] = "darkblue"
      } else if (cor_HOMA.IR$keggmodules [k , j, "estimate"] < 0){
        phenobar[k, j] = "darkred"
      } else {
        phenobar[k, j] = "white"
      }
    }
  }
}

### Rename column names of phenobar
colnames (phenobar) = c("HOMA-IR", "HOMA-IR.BMI.adj")

### Note, normally one would cluster the rows and columns in the heatmap as
### shown below but for the paper (Pedersen et al, 2016) we needed to group the
### KEGG modules by biological similarity to make a higher-level annotation for
### the figure and thus made a manual arrangement.

### cluster rows of matrix 
d = dist (plotmat)
dend = as.dendrogram (hclust (d, method = "average"))
rm (d)

### cluster columns of matrix 
d2 = dist (t (plotmat))
dend2 = as.dendrogram (hclust (d2, method = "average"))
rm (d2)

### make heatmap with all the selected KEGG modules
pdf ("results/heatmap_KEGG_vs_metabolite_clusters.pdf", height = 10, width = 10)
heatmap.3 (plotmat,
           Rowv = dend, 
           Colv = dend2, 
           dendrogram = "column",
           col = rev (redblue (100)),
           symbreaks = T,
           key = T,
           symkey = T,
           keysize = 1,
           KeyValueName = "SCC.bg.adj",
           lhei = c(0.6, 4),
           cellnote = stars, 
           notecol = "black",
           notecex = 1.1,
           trace = "none",
           labRow = module_mapping_clean [rownames (plotmat)],
           labCol = cluster_mapping_file [match( colnames (plotmat), cluster_mapping_file$New_Name), "label"],
           RowSideColors = t (phenobar [rownames (plotmat), ]), 
           side.height.fraction = 0.35,    
           NumColSideColors = dim (phenobar) [2], 
           margins = c (23, 33),
           cexRow = 1,
           cexCol = 1
)
dev.off () 
