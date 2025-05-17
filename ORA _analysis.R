library(clusterProfiler)
library(readxl)
library(writexl)
library(org.Hs.eg.db)
library(ReactomePA)
library(dplyr)
library(openxlsx)
library(ggplot2)


set.seed(123)
setwd("path to folder")

gene_entrez_<-read_xlsx("gene_list",sheet = 1)%>%as.data.frame()
gene_entrez <- gene_entrez_$Entrez.ID.


go_bp_up <- enrichGO(gene = gene_entrez,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID",
                     ont = "BP", 
                     pvalueCutoff =1,
                     pAdjustMethod = "BH",
                     minGSSize = 3,
                     maxGSSize = 900,
                     readable = TRUE)

go_mf_up <- enrichGO(gene = gene_entrez,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID",
                     ont = "MF",  
                     pAdjustMethod = "BH",
                     pvalueCutoff = 1,
                     minGSSize = 3,
                     maxGSSize = 900,
                     readable = TRUE)

go_cc_up <- enrichGO(gene = gene_entrez,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID",
                     ont = "CC",  
                     pAdjustMethod = "BH",
                     pvalueCutoff = 1,
                     minGSSize = 3,
                     maxGSSize = 900,
                     readable = TRUE)

reactome_enrich_up <- enrichPathway(gene = gene_entrez, 
                                    organism = "human", 
                                    pvalueCutoff = 1,
                                    pAdjustMethod = "BH",
                                    minGSSize = 3,
                                    maxGSSize = 900,
                                    readable = TRUE)
reactome_enrich_up <- setReadable(reactome_enrich_up, "org.Hs.eg.db", "ENTREZID")


sheets <- list(
  "GO_BP_UP" = subset(as.data.frame(go_bp_up@result), pvalue <= 0.05),
  "GO_MF_UP" = subset(as.data.frame(go_mf_up@result), pvalue <= 0.05),
  "GO_CC_UP" = subset(as.data.frame(go_cc_up@result), pvalue <= 0.05),
  "Reactome_UP" = subset(as.data.frame(reactome_enrich_up@result), pvalue <= 0.05)
)

write_xlsx(sheets, "enrched_results.xlsx")
#############################

library(DOSE)

dot_plot_bp<-dotplot(go_bp_up, showCategory=10,orderBy = "pvalue",color = "pvalue", title = "GO:BP") 
dot_plot_mf<-dotplot(go_mf_up, showCategory=10,orderBy = "pvalue",color = "pvalue", title = "GO:MF") 
dot_plot_cc<-dotplot(go_cc_up, showCategory=10,orderBy = "pvalue",color = "pvalue", title = "GO:CC")
dot_plot_r<-dotplot(reactome_enrich_up, showCategory=10,orderBy = "pvalue",color = "pvalue", title= "Reactome") 

library(gridExtra)
library(grid)

combined_plot <- grid.arrange(dot_plot_bp,dot_plot_cc, dot_plot_mf, dot_plot_r, 
                              top = textGrob("Functional Enrichment Analysis  of Ugregulated genes(Group A vs.  Group B)",gp=gpar(fontsize=20)),
                              ncol = 2) 

ggsave(filename = "Dotplot_upregulated_genes", plot = combined_plot, 
       device = "pdf", 
       width = 18, height = 18)

