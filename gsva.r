#load library

library(dplyr)
library(GSVA)
library(GSEABase)
library(tidyverse)
library(readxl)
library(writexl)
library(openxlsx)
library(clusterProfiler)

setwd("path to folder")
set.seed(123)
gene_expression<-read_xlsx("gene_expression.xlsx",sheet = 1)%>%as.data.frame()
row.names(gene_expression)<-gene_expression[,1]
gene_expression <- gene_expression[,-1]

gene_expression <- gene_expression %>% mutate(across(everything(), as.numeric))
gene_expression <- log2(gene_expression+1)
gene_expression <- as.matrix(gene_expression)


gene_sets <- read.gmt("Immune_Response_Type.gmt")
class(gene_sets)
gene_sets_1 <- as.list(split(gene_sets$gene, gene_sets$term))



zscore_params <- zscoreParam(
  exprData = gene_expression,
  geneSets = gene_sets_1,
  minSize = 1,
  maxSize = 600
  )
gsva_results <- gsva(zscore_params,
                     verbose = FALSE)

gsva_results <-as.data.frame(gsva_results)
gsva_results$Cell_Type <- rownames(gsva_results)
gsva_results <-gsva_results[,c(3,1,2)]
write_xlsx(gsva_results,"Cell_Type_GeneSet.xlsx")

