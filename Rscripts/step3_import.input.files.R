###
### This is step3_import.input.files.R, part of the Nature Protocol: 
### 'A computational framework for conducting a three-pronged association study
### to integrate high-throughput ‘-omics’ datasets for the identification of
### potential mechanistic links'.
### 26/7-2018
###

###
### Step 3 - Import input files
### 
### The following input data is assumed to exist within the top/data subdirectory:
###
### Data-files: 
### -	phenotypes.tab: 
###   File with clinical phenotypes (columns) per
###   individual (rows). Used to test for associations with, or for confounder
###   analysis. Individuals are labeled ‘idv’ followed by a number, e.g. ‘idv001’.
###
### -	metabolomic.tab / lipidomic.tab: 
###   Input data matrix for abundance of 325
###   polar metabolites or 876 molecular lipids per individual. Note that no
###   additional normalization is done in this script, so data is assumed to be
###   comparable in these regards. Such different data types are eventually merged
###   into a single set of metabolite cluster abundances. Individual
###   metabolite/lipids are named M or L (for specifying a polar metabolite or
###   molecular lipid, respectively), followed by a number and lastly the
###   annotation or ‘unknown’ in case of unannotated metabolites/lipids, e.g.
###   ‘M_20_Valine’.
###
### -	MGS_abundance.tab: 
###   File with the abundance (e.g. median gene abundance) of
###   MGSs (columns) per individual (rows). These are assumed to have been
###   rarefied to comparable depth or otherwise normalized. For historical
###   reasons, the MGSs are labeled ‘T2DCAG’ followed by a number, e.g.
###   ‘T2DCAG00001’.
###
### -	KO_abundance.tab: 
###   File with the abundance of each KO (columns) per
###   individual (rows). The data is assumed to be rarefied to comparable depth or
###   otherwise normalized. For this, the software tool rtk58 can be used.
###
### - gene_abundance_sub.tab: 
###   File with the abundance of each catalog gene
###   (subset-version) in each individual assumed to be rarefied to comparable
###   depth or otherwise normalized.
### 
### Annotation-files:
### -	cluster_mapping_file.tab: 
###   Input file with annotation for metabolite
###   clusters, as available from curation of data in the specific dataset. The
###   WGCNA clustering algorithm names the generated clusters with color codes.
###   This mapping file simply facilitate renaming to more meaningful cluster
###   descriptions. Here, the serum polar metabolite and serum molecular lipid
###   clusters are labelled M01–M35 and L01–L39, respectively, and collectively
###   termed metabolite clusters.
###
### -	MGS_taxonomy.tab: 
###   File with taxonomic annotation of the MGSs (rows) used
###   in the analysis. Each row contains the following information:
###   species_taxonomy, species_pct, genus_taxonomy, genus_pct, family_taxonomy,
###   family_pct, order_taxonomy, order_pct, phylum_taxonomy, phylum_pct, where
###   x_pct is the percentage of the MGS genes that can be annotated (by sequence
###   similarity) to the taxonomy of the MGS. If no taxonomy can be assigned to
###   the MGS, the value will be NA.
###
### -	KEGG_modules.tab: 
###   File containing definition of KEGG gene functional
###   modules (rows) specifying which KO gene groups constitute each KEGG module.
###   The first two columns contain KEGG module entry (number) and name, the third
###   column list all KOs separated by semicolon. Any other functional annotation
###   used analogously could be swapped in instead of the name. This file can be
###   obtained by downloading KEGG modules from
###   www.genome.jp/kegg-bin/get_htext?ko00002.keg (www.genome.jp/kegg/ --> KEGG
###   MODULE --> KEGG modules); download htext and then running the provided
###   script ‘parse_kegg.pl’ (after changing input filename).
###
### -	KO_to_MGS.tab: 
###   File listing for each KO (rows) the MGSs (space-separated)
###   it is a member of. This is based on the KO annotation of the gene catalog
###   (gene_to_KO.tab) and the information, what gene from the gene catalog is
###   within each MGS as given in the list MGS_to_gene.tab.
###
### -	gene_to_KO.tab: 
###   File containing the KO annotation (if any) per gene (rows)
###   in the gene catalog (constituting 7,328,469 genes). Genes are labeled
###   ‘RefCat620’ followed by a number (1…7328469), e.g. ‘RefCat620.1’. Used to
###   create KO_to_MGS.tab. Note, for the purpose of this protocol and to
###   considerably reduce size of input data, only the subset of catalogue genes
###   with KO annotation (n = 2,205,769) are provided in files with gene abundance
###   or annotation (gene_to_KO.tab, MGS_to_gene.tab and gene_abundance_sub.tab)
###   as only those genes are used in the driver-species analysis.
###
### -	MGS_to_gene.tab: 
###   List of genes binned into a given MGS (rows). Used to
###   create KO_to_MGS.tab. Rather than gene names the file contains the index
###   value of the gene. i.e. position of the gene in the gene catalogue
###   (subset-version, thus 1…2205769).
###
### All demonstration files are tab-delimited text files, but other formats
### would equally work after modifying the respective file-import commands in
### the R-script. 
###
### For all demonstration data, pseudonymised sample names were re-randomised to
### generate anonymised data.
###

options (stringsAsFactors = FALSE)

### - phenotypes.tab

phenotypes = read.table (file = "data/phenotypes.tab", row.names = 1, header = T, sep = "\t")
ctrl = rownames (subset (phenotypes, Diabetes == "nonDiabetic")) ### set of control samples for analyses restricted to non-diabetic individuals.

### - metabolomic.tab

metabolomic = read.table (file = "data/metabolomic.tab", row.names = 1, header = T, sep = "\t")

### - lipidomic.tab

lipidomic = read.table (file = "data/lipidomic.tab", row.names = 1, header = T, sep = "\t")

### - cluster_mapping_file.tab

cluster_mapping_file = read.table (file = "data/cluster_mapping_file.tab", header = T, sep = "\t", row.names = 1)
cluster_mapping_file$label =  sapply (rownames (cluster_mapping_file),  function (x) paste (cluster_mapping_file [x, "New_Name"], cluster_mapping_file [x, "Description"], sep = ": "))

### - MGS_abundance.tab

mgs_abundance = read.table (file = "data/MGS_abundance.tab", sep = "\t", row.names = 1, header = T)

### - MGS_taxonomy.tab

mgs_taxonomy  = read.table (file = "data/MGS_taxonomy.tab", sep = "\t", row.names = 1, header = T)

### - KEGG_modules.tab

tmp = read.table ("data/KEGG_modules.tab", sep = "\t")
koann = strsplit (tmp[,3], split = ";")
names (koann) = tmp[,1]

module_mapping = tmp[,2] ### description of Kegg modules
names (module_mapping) = tmp[,1] ; rm (tmp)

### remove KEGG references in square brackets for more clean names for plotting
module_mapping_clean = sapply(module_mapping, function(x) strsplit(x, " \\[")[[1]][1])

### - KO_abundance.tab

ko_abundance = read.table (file = "data/KO_abundance.tab", sep = "\t", row.names = 1, header = T)

### - KO_to_MGS.tab

tmp = read.table ("data/KO_to_MGS.tab", sep = "\t", strip.white = T)
KO2MGS = strsplit (tmp[,2], split = " ")
names (KO2MGS) = tmp[,1] ; rm (tmp)

### - gene_to_KO.tab

gene2KO = read.table ("data/gene_to_KO.tab", sep = "\t", row.names = 1, header = T)
KO2gene = tapply (1:nrow (gene2KO), gene2KO[,1], c) ; rm (gene2KO)

### - MGS_to_gene.tab

tmp = read.table ("data/MGS_to_gene.tab", sep = "\t", strip.white = T)
MGS2gene = strsplit (tmp[,2], split = " ")
names (MGS2gene) = tmp[,1] ; rm (tmp)

### - gene_abundance_sub.tab

gene_abundance_sub = data.frame( fread ("data/gene_abundance_sub.tab", sep = "\t", header = T), row.names = 1)
