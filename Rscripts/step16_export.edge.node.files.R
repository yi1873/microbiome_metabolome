###
### This is step16_export.edge.node.files.R, part of the Nature Protocol: 
### 'A computational framework for conducting a three-pronged association study
### to integrate high-throughput ‘-omics’ datasets for the identification of
### potential mechanistic links'.
### 26/7-2018
###

###
### Step 16 (optional) - Export metabolome-microbiome associations for network
### analysis
###
### Further exploration of such high dimensional association data can be
### performed using software for network analysis, e.g. Cytoscape
### (www.cytoscape.org) or igraph (www.igraph.org). Here we provide code for
### exporting an edge-file with pairwise association scores and FDR-values
### between metabolite clusters and KEGG modules as well as a corresponding node
### attribute file; both in .txt-format. Several tutorials for importing,
### visualizing and analyzing networks in Cytoscape can be found here:
### https://github.com/cytoscape/cytoscape-tutorials/wiki.
###

### Note, for simplicity, the node and edge files are made for a network
### representing the heatmap plotted in Step 15 (i.e. KEGG modules in
### "keggmodules2test" and metabolite clusters in "metlip2test")

### Make edge file with the following four columns: 
### - KEGG_module, 
### - Metabolite_cluster, 
### - association estimate (SCC.bg.adj), 
### - FDR adjusted p-values
tmp.edge.file.p_adj = melt (kegg_metlip_p_adj [keggmodules2test, metlip2test])
tmp.edge.file.est   = melt (kegg_metlip_est [keggmodules2test, metlip2test])
colnames (tmp.edge.file.p_adj) = c ("KEGG_module", "Metabolite_cluster", "p.adjust")
colnames (tmp.edge.file.est)   = c ("KEGG_module", "Metabolite_cluster", "SCC.bg.adj")
edge.file = cbind (tmp.edge.file.est, "p.adjust" = tmp.edge.file.p_adj$p.adjust)
rm (tmp.edge.file.p_adj, tmp.edge.file.est)

### Save the edge-file as a tab-seperated txt.file
write.table (edge.file, file = "results/edge_file.txt", sep = "\t", row.names = F, col.names = T, quote = F)

### Make node annotation file with the following three columns: 
### - Node_name (corresponds to names in the edge-file), 
### - Description (where existing), 
### - Examples (example metabolites, only for metabolite clusters)
### Get metabolite annotation
node.file = cluster_mapping_file [cluster_mapping_file$New_Name %in% metlip2test, c ("New_Name", "Description", "Examples")]
colnames (node.file) = c ("Node_name", "Description", "Examples")
### Add KEGG module annotation
node.file = rbind (node.file, data.frame ("Node_name" = keggmodules2test, "Description" = module_mapping_clean [keggmodules2test], "Examples" = rep ("", time = length (keggmodules2test))))

### Save the node-file as a tab-seperated txt.file
write.table (node.file, file = "results/node_file.txt", sep = "\t", row.names = F, col.names = T, quote = F)
