library(clusterProfiler)
library(org.Hs.eg.db)
library(readxl)
library(writexl)
library(openxlsx)
library(ReactomePA)

df = read_xlsx("DEG.xlsx",sheet = 1)
gene_conversion <- bitr(df$Gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
df_converted <- merge(df, gene_conversion, by.x = "Gene", by.y = "SYMBOL")


original_gene_list <- -log10(df_converted$pvalue) * sign(df_converted$FC)
# name the vector
names(original_gene_list) <- df_converted$ENTREZID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

###GSEA
Gsea_BP<- gseGO(geneList = gene_list,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",  
                     nPerm = 1000,
                     minGSSize = 10,
                     maxGSSize = 500,
                     pvalueCutoff = 1,
                     verbose = TRUE)

Gsea_MF <- gseGO(geneList = gene_list,
                     OrgDb = org.Hs.eg.db,
                     ont = "MF",  
                     nPerm = 1000,
                     minGSSize = 10,
                     maxGSSize = 500,
                     pvalueCutoff = 1,
                     verbose = TRUE)

Gsea_reactome <- gsePathway(geneList = gene_list,
                            organism = "human", 
                            nPerm = 1000,
                            minGSSize = 10,
                            maxGSSize = 500,
                            pvalueCutoff = 1,
                            pAdjustMethod = "BH",
                            verbose = TRUE)


bp_results <- as.data.frame(Gsea_BP) %>% dplyr::filter(pvalue < 0.05)
mf_results <- as.data.frame(Gsea_MF) %>% dplyr::filter(pvalue < 0.05)
reactome_results <- as.data.frame(Gsea_reactome) %>% dplyr::filter(pvalue < 0.05)

# Create a new workbook
wb <- createWorkbook()

# Add sheets for each result
addWorksheet(wb, "BP")
addWorksheet(wb, "MF")
addWorksheet(wb, "Reactome")

# Write the results to each sheet
writeData(wb, sheet = "BP", bp_results)
writeData(wb, sheet = "MF", mf_results)
writeData(wb, sheet = "Reactome", reactome_results)

# Save the workbook
saveWorkbook(wb, "GSEA_results.xlsx", overwrite = TRUE)

#plot

library(DOSE)
dot_plot_bp<-dotplot(Gsea_BP, showCategory=20, split=".sign",color = "pvalue") + facet_grid(.~.sign)
dot_plot_mf<-dotplot(Gsea_MF, showCategory=20, split=".sign",color = "pvalue") + facet_grid(.~.sign)
dot_plot_r<-dotplot(Gsea_reactome, showCategory=20, split=".sign",color = "pvalue") + facet_grid(.~.sign)


ggsave(filename = "dotplot_reactome.tiff", plot = dot_plot_r, 
       device = "tiff", 
       width = 22, height = 22, units = "in", 
       dpi = 300)
ggsave(filename = "dotplot_BP.tiff", plot = dot_plot_bp, 
       device = "tiff", 
       width = 22, height = 22, units = "in", 
       dpi = 300)
ggsave(filename = "dotplot_MF.tiff", plot = dot_plot_mf, 
       device = "tiff", 
       width = 22, height = 22, units = "in", 
       dpi = 300)

