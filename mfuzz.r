######################################################################################

Mfuzz is an R package used for clustering genes based on time-series expression data. 
It performs soft clustering using the fuzzy c-means algorithm, allowing genes to belong to multiple clusters based on their expression patterns. 
This method helps identify co-expressed and potentially co-regulated genes by assigning a membership value (ranging from 0–1) to each gene. 
The analysis requires normalized gene expression data as input.

######################################################################################

#install the package
BiocManager::install("Mfuzz")

#load the package
library(Mfuzz)
library(readxl)
library(Biobase)
library(dplyr)
library(limma)


setwd("../Deepika/mfuzz")
set.seed(123)

#preparing the data for soft clustering
#load a gene expression matrix
data <- read_xlsx("Normalised_gene expression data.xlsx",sheet = 1)%>%as.data.frame()
row.names(data)<-data$Genes
data<-data[,-1]
data <- data%>%mutate(across(everything(),as.numeric))


#load metadata
metadata <- read_xlsx("Annotations.xlsx",sheet = 1)%>%as.data.frame()
metadata <-metadata[order(metadata$Diagnosis),]
row.names(metadata)<-metadata$Sample_ID
metadata <- metadata[,-1,drop = F]

unique(metadata$Diagnosis)
identical(rownames(metadata),colnames(data))

data <- t(data)%>%as.data.frame()
identical(rownames(metadata),rownames(data))
data<-as.matrix(data)

#take the average of gene expression based on sample metadata
average_expression <- avereps(data, ID = metadata$Diagnosis)
average_expression <- as.data.frame(average_expression)
average_expression <- t(average_expression)%>%as.data.frame()
#write.csv(average_expression,"average_expression_expression_data.csv")


# mfuzz
average_expression <-as.matrix(average_expression)          
eset <- new("ExpressionSet", exprs = average_expression)
eset

# preprocessing steps
# filter genes which has missing value above the threshold
eset_s <- filter.NA(eset, thres=0.25)

# fill missing values by the average values expression
eset_s <- fill.NA(eset_s,mode="mean")

# exclude  genes with small changes in expression or variation
eset_s <- filter.std(eset_s,min.std=0,visu = T)# this step can be excluded as mfuzz perform soft clustering and it is robust to noise.However, if the number of
genes with small expression changes is large, such pre-filtering may be necessary.

#the expression values of genes were standardised to have a mean value of zero and a standard deviation of one
eset_s <- standardise(eset)

# Estimate for optimal fuzzifier m
m_val <- mestimate(eset_s)
m_val

#estimation of optimised number of clusters
Dmin(eset_s, m = m_val,crange=seq(2, 20, 1),repeats = 100, visu = TRUE)

#Function for soft clustering
cl2 <- mfuzz(eset_s,c= 9, m_val)


pdf("mfuzz_clustered_genes.pdf") 
mfuzz.plot(eset_s, cl = cl2, mfrow = c(1, 1),
           time.labels = colnames(eset_s), new.window  = FALSE)
dev.off()


write.csv(cl2$membership,file = "memship_value.csv")
write.csv(cl2$cluster,file = "mfuzz_cluster.csv")
write.csv(cl2$size,file = "mfuzz_size.csv")

#based on the membership value you can select the cutoff to extract the genes that belong to each cluster.

######################################################
