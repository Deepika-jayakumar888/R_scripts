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
go_bp_up <- setReadable(go_bp_up, "org.Hs.eg.db", "ENTREZID")


go_mf_up <- enrichGO(gene = gene_entrez,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID",
                     ont = "MF",  
                     pAdjustMethod = "BH",
                     pvalueCutoff = 1,
                     minGSSize = 3,
                     maxGSSize = 900,
                     readable = TRUE)
go_mf_up <- setReadable(go_mf_up, "org.Hs.eg.db", "ENTREZID")


go_cc_up <- enrichGO(gene = gene_entrez,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID",
                     ont = "CC",  
                     pAdjustMethod = "BH",
                     pvalueCutoff = 1,
                     minGSSize = 3,
                     maxGSSize = 900,
                     readable = TRUE)
go_cc_up <- setReadable(go_cc_up, "org.Hs.eg.db", "ENTREZID")

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


go_bp_df <- as.data.frame(go_bp_up@result)
go_mf_df <- as.data.frame(go_mf_up@result)
go_cc_df <- as.data.frame(go_cc_up@result)
reactome_df <- as.data.frame(reactome_enrich_up@result)

go_bp_df$Regulation <- "GO:BP"
go_mf_df$Regulation <- "GO:MF"
go_cc_df$Regulation <- "GO:CC"
reactome_df$Regulation <- "Reactome"

reactome_df$FoldEnrichment <- 0
reactome_df$zScore <- 0
reactome_df$RichFactor <- 0



combined_df <- rbind(go_bp_df,go_cc_df, reactome_df,go_mf_df)
unique(combined_df$Regulation)

top10_combined_df <- combined_df %>%
  filter(pvalue <= 0.05) %>% 
  group_by(Regulation) %>%
  arrange(pvalue) %>%
  slice_head(n = 10) %>%
  ungroup()


pdf("dotplot.pdf", width = 14, height = 8)
dotplot_combined <- ggplot(top10_combined_df, aes(x = Count, y = reorder(Description, Count), color = pvalue)) +
  geom_point(aes(size = Count)) +
  scale_color_gradient(low = "red", high = "blue") +
  facet_wrap(~ Regulation, scales = "free_y",ncol = 2, nrow = 2) +
  theme_minimal() +
  labs(title = "Functional Enrichment Analysis of  ",
       x = "Gene Count",
       y = "GO Term")+
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"), 
    legend.title = element_text(face = "bold")
  )

print(dotplot_combined)

dev.off()

