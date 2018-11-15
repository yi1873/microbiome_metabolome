###
### This is step19_plot.leave.one.MGS.out.results.R, part of the Nature Protocol: 
### 'A computational framework for conducting a three-pronged association study
### to integrate high-throughput ‘-omics’ datasets for the identification of
### potential mechanistic links'.
### 26/7-2018
###

###
### Step 19 - Plotting leave-one-MGS-out results
###
### Having computed these results, we plot the top driver taxa for the
### phenotype-associated gene functions (Figure 4).
### These include sub-plots for:
### - Distribution/density plots of correlations for KOs in a KEGG module vs 
###   all other KOs
### - Distribution of correlations when leave-one-MGS-out. 
### - Distribution of correlations when leave-one-MGS-out. bg.adj. SCC 
###   (i.e. what was shown in Figure 3c+d in Pedersen et al., 2016)
### Median SCC for KO within a module (red) and all other remaining KOs (green)
### are indicated in the first two plots.
###

### Specify pdf-file for output.
pdf (file = "results/density_plot_SCC_HOMA.IR.pdf", width = 6 * 1, height = 3 * 3)

for (m in names (KOsets)) {
  
  ### Distribution of correlations for KOs in modules vs all other KOs
  incat =     na.omit (KO_cor_HOMA [names (KO_cor_HOMA) %in% KOsets [[m]] ]) ### select all correlations between HOMA-IR and KOs in the KEGG module
  incat2 =    data.frame ("KO_cor_HOMA" = incat, "cat" = "KOs in module")
  notincat =  na.omit (KO_cor_HOMA [! (names (KO_cor_HOMA) %in% KOsets [[m]])]) ### select all correlations between HOMA-IR and KOs NOT in the KEGG module   
  notincat2 = data.frame ("KO_cor_HOMA" = notincat, "cat" = "KOs not in module")
  tmp.df = rbind (incat2, notincat2)
  cdf <- ddply (tmp.df, "cat", summarise, rating.median = median (KO_cor_HOMA)) ### find median for each group
  g1 <- ggplot (tmp.df, aes (x = KO_cor_HOMA, fill = cat)) + 
    geom_density (alpha = 0.4) +
    geom_vline (data = cdf, aes (xintercept = rating.median, colour = cat), linetype = "dashed", size = 0.8) +
    ggtitle (paste0 ("KEGG module: ", m, 
                     "\n", strsplit (module_mapping [m], split = "\\[|,")[[1]][1], 
                     "\n", length (KOsets [[m]]), " KOs in module vs all remaining ", (length (KO_cor_HOMA) - length (KOsets[[m]])), " KOs", "\n")) +
    xlab ("SCC for KOs and HOMA-IR") +
    theme_classic () + theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ())
  
  ### Plot Distribution of correlations when leave-one-MGS-out.   
  if (nrow (DeltaSCCperMGS [[m]]) > 1) {
    g2 <- ggplot (DeltaSCCperMGS [[m]], aes (x = SCC_omiting_MGS)) + 
      geom_density (alpha = 1, fill = "grey", aes (y = ..scaled..)) + 
      geom_segment (aes (y = -0.1, yend = -0.02, x = SCC_omiting_MGS, xend = SCC_omiting_MGS)) +
      geom_vline (data = cdf, aes (xintercept = rating.median, colour = cat), linetype = "dashed", size = 0.8) +
      ggtitle (paste0 ("KEGG module: ", m, 
                       "\n", strsplit (module_mapping [m], split = "\\[|,")[[1]][1], 
                       "\nSCC for leave-1-MGS-out",
                       "\nNumber of MGSs = ", nrow (DeltaSCCperMGS [[m]]))) +
      xlab ("SCC for KOs and HOMA-IR") +
      theme_classic () + theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ())
    
  } else {
    
    ### Special circumstance when there is only one MGS (then one cannot make a 
    ### density plot, instead a black line is plotted at the value when that 
    ### one MGS is left out).
    
    g2 <- ggplot () +
      scale_x_continuous (limits = range (c (DeltaSCCperMGS [[m]][, "SCC_omiting_MGS"], cdf$rating.median))) +
      scale_y_continuous (name = "", limits = c (0, 1)) +
      geom_vline (data = cdf, aes (xintercept = rating.median, colour = cat), linetype = "dashed", size = 0.8) +
      geom_vline (data = DeltaSCCperMGS [[m]], aes (xintercept = SCC_omiting_MGS), linetype = "longdash", color = "black",size = 1) + 
      ggtitle (paste0 ("KEGG module: ", m, 
                       "\n", strsplit (module_mapping [m], split = "\\[|,") [[1]][1], 
                       "\nSCC for leave-1-MGS-out",
                       "\nNumber of MGSs = ", nrow (DeltaSCCperMGS [[m]]))) +
      xlab ("SCC for KOs and HOMA-IR") +
      theme_classic () + theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ()) 
    
  }
  
  ### Plot distribution of correlations when leave-one-MGS-out. Bg.adjust SCC 
  ### (i.e. what we show in Figure 3c+d)
  tmp.DeltaSCCperMGS = DeltaSCCperMGS [[m]]
  ### Calculate delta-SCC in relation to the background adjusted SCC. 
  ### i.e the 'plottedvalue.s.m' shown as the second last equation in the 
  ### methods section in Pedersen et al., 2016.
  tmp.DeltaSCCperMGS$DeltaMGS_SCC.bgadj = tmp.DeltaSCCperMGS$SCC.bgadj - tmp.DeltaSCCperMGS$DeltaMGS_SCC
  ### specify the range of the x-axis to include 0, i.e. [min:0] for negative 
  ### correlations and [0:max] for positive correlations
  x.range = c (min (0, tmp.DeltaSCCperMGS$DeltaMGS_SCC.bgadj), max (0, tmp.DeltaSCCperMGS$DeltaMGS_SCC.bgadj))
  
  if (nrow (tmp.DeltaSCCperMGS) > 1 ) {
    
    g3 <- ggplot (tmp.DeltaSCCperMGS, aes (x = DeltaMGS_SCC.bgadj)) + 
      geom_density (alpha = 1, fill = "grey", aes (y = ..scaled..)) + 
      geom_segment (aes (y = -0.1, yend = -0.02, x = DeltaMGS_SCC.bgadj, xend = DeltaMGS_SCC.bgadj)) +
      ggtitle (paste0 ("KEGG module: ", m, 
                       "\n", strsplit (module_mapping [m], split = "\\[|,")[[1]][1], 
                       "\nbg.adj.SCC for leave-1-MGS-out",
                       "\nNumber of MGSs = ", nrow (tmp.DeltaSCCperMGS))) +
      xlab ("bg.adj.SCC for KOs and HOMA-IR") + xlim (x.range) +
      theme_classic () + theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ())
    
  } else { 
    
    ### Special circumstance when there is only one MGS (then one cannot make a 
    ### density plot)
    
    g3 <- ggplot () +
      ggtitle (paste0 (m, 
                       "\n", strsplit (module_mapping [m], split="\\[|,")[[1]][1], 
                       "\nbg.adj.SCC for leave-1-MGS-out",
                       "\nNumber of MGSs = ", nrow (tmp.DeltaSCCperMGS), " - consequently not showing a density/rug plot")) 
    
  }
  
  g = plot_grid (g1, g2, g3, ncol = 1, nrow = 3, align = "v", labels = c ("a", "b", "c"), axis = "rl")
  # g = plot_grid (g1, g3, ncol = 1, nrow = 2, align = "v", labels = c ("a", "b", "c"), axis = "rl") ### use this if deleting g2
  plot (g)
  
}

dev.off ()  
