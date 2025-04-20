#load library

library(readxl)
library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyverse)
library(writexl)
library(openxlsx)
library(tidyverse)

setwd("path to folder")
data <- read_excel("DEG List.xlsx",sheet = 1)

original_gene_list <- -log10(data$Pvalue) * sign(data$FC)
names(original_gene_list) <-data$Gene

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# custom gmt files upload
msigdb_bp <- read.gmt("m5.go.bp.v2024.1.Mm.symbols.gmt")
msigdb_cc <- read.gmt("m5.go.cc.v2024.1.Mm.symbols.gmt")
msigdb_mf <- read.gmt("m5.go.mf.v2024.1.Mm.symbols.gmt")
msigdb_reactome <-read.gmt("m2.cp.reactome.v2024.1.Mm.symbols.gmt")

set.seed(123)

#perform GSEA
gsea_bp <- GSEA(gene_list, TERM2GENE = msigdb_bp, pvalueCutoff = 1,
                     minGSSize = 5,
                     maxGSSize = 900)

gsea_cc <- GSEA(gene_list, TERM2GENE = msigdb_cc, pvalueCutoff = 1,
                minGSSize = 5,
                maxGSSize = 900)

gsea_mf <- GSEA(gene_list, TERM2GENE = msigdb_mf, pvalueCutoff = 1,
                minGSSize = 5,
                maxGSSize = 900)

gsea_reactome <- GSEA(gene_list, TERM2GENE = msigdb_reactome, pvalueCutoff = 1,
                minGSSize = 5,
                maxGSSize = 900)

wb <- createWorkbook()

# Add worksheets to the workbook for each result
addWorksheet(wb, "GSEA_BP")
addWorksheet(wb, "GSEA_CC")
addWorksheet(wb, "GSEA_MF")
addWorksheet(wb, "GSEA_Reactome")

# Write data to each sheet
writeData(wb, sheet = "GSEA_BP", subset(as.data.frame(gsea_bp), pvalue < 0.05))
writeData(wb, sheet = "GSEA_CC", subset(as.data.frame(gsea_cc), pvalue < 0.05))
writeData(wb, sheet = "GSEA_MF", subset(as.data.frame(gsea_mf), pvalue < 0.05))
writeData(wb, sheet = "GSEA_Reactome", subset(as.data.frame(gsea_reactome), pvalue < 0.05))

# Save the workbook
saveWorkbook(wb, "GSEA_Enrichment_Results.xlsx", overwrite = TRUE)


#plot

pathways_of_interest <- gsea_reactome[c("REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT" ,
                                  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT"),] 
foldchange_list <- data$Estimate
names(foldchange_list) <- data$Gene

cnetplot(gsea_reactome,
         showCategory = pathways_of_interest$Description, 
         categorySize="pvalue",foldChange = foldchange_list
)+
  scale_color_gradientn(
    colors = c("darkgreen","white", "darkred"),
    name = "Fold Change",
    limits = c(min(foldchange_list), max(foldchange_list))
  )

