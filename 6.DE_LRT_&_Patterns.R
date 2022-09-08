### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
###   Differential expression analysis across Braak lewy body stages ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


### Packages required:
library("DESeq2")
library("tibble")
library("DEGreport")
library("ggplot2")
library("ggpubr")
library("dplyr")
library("tidyr")


# load dds
load("./dds_84_qSvs.Rda")
dim(dds)

coldata <- colData(dds)
coldata <- as.data.frame(coldata)


# Check the design formula
design(dds)


### export counts as txt file to upload on GEO
dds_norm <- estimateSizeFactors(dds)
counts <- counts(dds_norm, normalized=TRUE)

write.table(counts, file="./counts.txt")


#########################
# Likelihood ratio test #
#########################

dds_LRT <- DESeq(dds, test="LRT", reduced= ~ Sex + Age_death + qSV1 + qSV2 + qSV3 + qSV4 + qSV5)

### Save results
res_LRT <- results(dds_LRT)
summary(res_LRT)
# Number of DE genes FDR < 0.1
sum(res_LRT$padj < 0.1, na.rm = TRUE) # 799

### subset significant genes 
res_LRT_0.05 <- results(dds_LRT, alpha=0.05)
summary(res_LRT_0.05)
# Number of DE genes FDR < 0.05
sum(res_LRT_0.05$padj < 0.05, na.rm = TRUE) # 266


### even though there are fold changes present they are not directly associated with the actual hypothesis test.
### We are unable to set a fold change criteria here since the statistic is not generated from any one pairwise comparison. 
### This list includes genes that can be changing in any number of combinations across the four factor levels. 


### Create a tibbl  e for LRT results
res_LRT_tb <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()


### Save file with all genes
write.csv(as.data.frame(res_LRT_tb), file="./LRT_genes.csv",  row.names = FALSE)

### Subset to return genes with padj < 0.1
sigLRT_genes_0.1 <- res_LRT_tb %>% 
  dplyr::filter(padj < 0.1)

### Create a tibbl  e for LRT results
res_LRT_0.05_tb <- res_LRT_0.05 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

### Subset to return genes with padj < 0.05
sigLRT_genes_0.05 <- res_LRT_0.05_tb %>% 
  dplyr::filter(padj < 0.05)

### Order results by pvalue
sigLRT_genes_0.05 <- sigLRT_genes_0.05[order(sigLRT_genes_0.05$pvalue),]
sigLRT_genes_0.1 <- sigLRT_genes_0.1[order(sigLRT_genes_0.1$pvalue),]

### Save file with significant genes
write.csv(as.data.frame(sigLRT_genes_0.05), file="./sigLRT_genes_0.05.csv",  row.names = FALSE)
write.csv(as.data.frame(sigLRT_genes_0.1), file="./sigLRT_genes_0.1.csv",  row.names = FALSE)



################################
# Expression pattern analysis  #
################################


### Get name of significant genes
sigLRT_genes_name_0.05 <- sigLRT_genes_0.05 %>% 
  pull(gene)

### Obtain rlog values for those significant genes

### Log transform for plotting of patterns
rld <- rlogTransformation(dds_LRT)
rld_mat <- assay(rld)

### Obtain rlog values for the significant genes
cluster_rlog_0.05 <- rld_mat[sigLRT_genes_name_0.05, ]

dds_LRT$Braak_aSyn_groups_Plot <- as.factor(dds_LRT$Braak_aSyn_groups_Plot)
dds_LRT$Braak_aSyn_groups_Plot_names <- as.character(dds_LRT$Braak_aSyn_groups_Plot)

clusters_0.05 <- degPatterns(cluster_rlog_0.05, 
                             metadata = colData(dds_LRT), 
                             summarize = "Braak_aSyn_groups", 
                             time = "Braak_aSyn_groups_Plot", 
                             col= NULL , 
                             minc = 1, # 15 default
                             consensusCluster = FALSE)

### Figure 5
cluster_0.5_plot <- clusters_0.05$plot

cluster_0.5_plot <- cluster_0.5_plot +
  scale_color_manual(values = c("Black")) +
  theme_minimal() +
  labs(title="Clusters of DE genes (FDR < 0.05)", col="Braak Stage") +
  theme(legend.position = "none", plot.title = element_text(size=16, margin=margin(0,0,0,0)))


ggsave("./Clusters_0.05_5.svg", plot = cluster_0.5_plot, device = "svg", width = 10, height = 10, units = 'in', dpi =300)

# Extract the genes in each group
cluster_groups <- clusters_0.05$df
group1 <- clusters_0.05$df %>%
  dplyr::filter(cluster == 1)

group2 <- clusters_0.05$df %>%
  dplyr::filter(cluster == 2)

group3 <- clusters_0.05$df %>%
  dplyr::filter(cluster == 3)

group4 <- clusters_0.05$df %>%
  dplyr::filter(cluster == 4)

group5 <- clusters_0.05$df %>%
  dplyr::filter(cluster == 5)

group6 <- clusters_0.05$df %>%
  dplyr::filter(cluster == 6)

group7 <- clusters_0.05$df %>%
  dplyr::filter(cluster == 7)

group8 <- clusters_0.05$df %>%
  dplyr::filter(cluster == 8)



### Subset genes present in each group
sigLRT_genes_0.05 <- as.data.frame(sigLRT_genes_0.05)
rownames(sigLRT_genes_0.05) <- sigLRT_genes_0.05$gene

sigLRT_genes_0.05_1 <- sigLRT_genes_0.05[group1$genes,]
sigLRT_genes_0.05_2 <- sigLRT_genes_0.05[group2$genes,]
sigLRT_genes_0.05_3 <- sigLRT_genes_0.05[group3$genes,]
sigLRT_genes_0.05_4 <- sigLRT_genes_0.05[group4$genes,]
sigLRT_genes_0.05_5 <- sigLRT_genes_0.05[group5$genes,]
sigLRT_genes_0.05_6 <- sigLRT_genes_0.05[group6$genes,]
sigLRT_genes_0.05_7 <- sigLRT_genes_0.05[group7$genes,]
sigLRT_genes_0.05_8 <- sigLRT_genes_0.05[group8$genes,]


write.csv(as.data.frame(sigLRT_genes_0.05_1), file="./sigLRT_genes_0.05_1_F.csv",  row.names = FALSE)
write.csv(as.data.frame(sigLRT_genes_0.05_2), file="./sigLRT_genes_0.05_2_F.csv",  row.names = FALSE)
write.csv(as.data.frame(sigLRT_genes_0.05_3), file="./sigLRT_genes_0.05_3_F.csv",  row.names = FALSE)
write.csv(as.data.frame(sigLRT_genes_0.05_4), file="./sigLRT_genes_0.05_4_F.csv",  row.names = FALSE)
write.csv(as.data.frame(sigLRT_genes_0.05_5), file="./sigLRT_genes_0.05_5_F.csv",  row.names = FALSE) 
write.csv(as.data.frame(sigLRT_genes_0.05_6), file="./sigLRT_genes_0.05_6_F.csv",  row.names = FALSE)
write.csv(as.data.frame(sigLRT_genes_0.05_7), file="./sigLRT_genes_0.05_7_F.csv",  row.names = FALSE)
write.csv(as.data.frame(sigLRT_genes_0.05_8), file="./sigLRT_genes_0.05_8_F.csv",  row.names = FALSE) 

