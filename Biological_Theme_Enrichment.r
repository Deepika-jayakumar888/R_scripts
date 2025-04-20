
#################
compareCluster,helps to  calculate enriched functional profiles of each gene clusters and aggregate the results into a single object. 

website link : https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-comparecluster.html

#################

#load library

library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(writexl)
library(readxl)
library(openxlsx)
library(enrichplot)
library(DOSE)
library(ReactomePA)



setwd("path to gene list file")

g_id<- read_xlsx("Gene_entrezid_list.xlsx",sheet = 1)%>%as.data.frame()
row.names(g_id)<- g_id[,1]
g_id<-g_id[,-1,drop = F]



data<- read.csv("gene_cluster.csv",header = T)%>%as.data.frame()


cluster_list <- list()

for (cluster_num in sort(unique(data$belongs.to.cluster))) {
  data_cluster <- data[data$belongs.to.cluster == cluster_num, ]
  vector <- data_cluster$X  # Extract gene names
  
  # Extract only gene IDs as a character vector
  id_vector <- as.character(g_id[rownames(g_id) %in% vector, , drop = TRUE])
  
  cluster_list[[paste0("Cluster_", cluster_num)]] <- id_vector
}



#perform compare cluster
enrich_bp <- compareCluster(cluster_list,
                         fun = "enrichGO",
                         OrgDb = org.Hs.eg.db,
                         ont = "BP",
                         keyType="ENTREZID",
                         pvalueCutoff=1,
                         pAdjustMethod = "BH",
                         qvalueCutoff = 1)
enrich_bp <- setReadable(enrich_bp, "org.Hs.eg.db", "ENTREZID")


enrich_mf <- compareCluster(cluster_list,
                            fun = "enrichGO",
                            OrgDb = org.Hs.eg.db,
                            ont = "MF",
                            keyType="ENTREZID",
                            pvalueCutoff=1,
                            pAdjustMethod = "BH",
                            qvalueCutoff = 1)
enrich_mf <- setReadable(enrich_mf, "org.Hs.eg.db", "ENTREZID")


enrich_reactome <- compareCluster(cluster_list,
                            fun = "enrichPathway",
                            organism = "human",
                            pvalueCutoff=1,
                            pAdjustMethod = "BH",
                            qvalueCutoff = 1,
                            minGSSize = 5,
                            maxGSSize = 900)
enrich_reactome <- setReadable(enrich_reactome, "org.Hs.eg.db", "ENTREZID")

##BP
enrich_bp_df <-as.data.frame(enrich_bp)
enrich_bp_df<-enrich_bp_df[enrich_bp_df$pvalue <= 0.05, ]
enrich_bp_df <-as.data.frame(enrich_bp_df)
write_xlsx(enrich_bp_df,"biological_process.xlsx")


enrich_bp@compareClusterResult <- enrich_bp@compareClusterResult[enrich_bp@compareClusterResult$pvalue <= 0.05, ]
pdf("Dotplot_BP_cluster.pdf",height = 48,width=14)
dotplot(enrich_bp,color = "pvalue",showCategory = 10)
dev.off()

##MF

enrich_mf_df <-as.data.frame(enrich_mf)
enrich_mf_df<-enrich_mf_df[enrich_mf_df$pvalue <= 0.05, ]
enrich_mf_df <-as.data.frame(enrich_mf_df)
write_xlsx(enrich_mf_df,"MF_enrichment_analysis.xlsx")


enrich_mf@compareClusterResult <- enrich_mf@compareClusterResult[enrich_mf@compareClusterResult$pvalue <= 0.05, ]
pdf("Dotplot_MF_cluster.pdf",height = 16,width=14)
dotplot(enrich_mf,color = "pvalue")
dev.off()

##reactome


enrich_reactome_df <-as.data.frame(enrich_reactome)
enrich_reactome_df<-enrich_reactome_df[enrich_reactome_df$pvalue <= 0.05, ]
enrich_reactome_df <-as.data.frame(enrich_reactome_df)
write_xlsx(enrich_reactome_df,"reactome_enrichment_analysis.xlsx")


enrich_reactome@compareClusterResult <- enrich_reactome@compareClusterResult[enrich_reactome@compareClusterResult$pvalue <= 0.05, ]
pdf("Dotplot_reactome_cluster.pdf",height = 18,width=14)
dotplot(enrich_reactome,color = "pvalue")
dev.off()


