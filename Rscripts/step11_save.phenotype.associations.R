###
### This is step11_save.phenotype.associations.R, part of the Nature Protocol: 
### 'A computational framework for conducting a three-pronged association study
### to integrate high-throughput ‘-omics’ datasets for the identification of
### potential mechanistic links'.
### 26/7-2018
###

###
### Step 11 - Save phenotype associations
###
### Save the (BMI corrected) HOMA-IR association of metabolite clusters, MGSs
### and KEGG modules calculated in Step 8-10.
###

tmp.excel.file = "results/HOMA.IR_associations.xlsx"
write.xlsx (paste ("Associations with HOMA-IR and HOMA-IR adjusted for BMI"), 
            sheetName = "info", file = tmp.excel.file, row.names = F, col.names = F)

for (tmp_name in names (cor_HOMA.IR)) {
  tmpMat = cor_HOMA.IR [[paste (tmp_name)]]
  out = cbind (tmpMat [, "HOMA.ir", c("estimate", "p.value")], "p.adjust" = p.adjust (tmpMat [, "HOMA.ir","p.value"], method = "BH"), 
               tmpMat [, "HOMA.ir_BMI.adj.partial", c ("estimate", "p.value")], "p.adjust" = p.adjust (tmpMat [, "HOMA.ir_BMI.adj.partial", "p.value"], method = "BH"))
  colnames (out) = paste (rep (c ("HOMA.IR", "HOMA.IR_BMI.adj.partial"), each = 3), colnames (out), sep = "_")
  if (tmp_name == "metlip") { rownames (out) = cluster_mapping_file [match( rownames (out), cluster_mapping_file$New_Name), "label"] }
  write.xlsx (x = out, file = tmp.excel.file, append = T, sheetName = paste (tmp_name, sep = "_"), row.names = T, col.names = T)            
}
