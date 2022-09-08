### ### ### ### ### ### ### ### ### ### ### ### ### 
###   Differential pairwise expression analysis  ###
### ### ### ### ### ### ### ### ### ### ### ### ### 

### Packages required:
library("DESeq2")
library("tibble")
library("ggplot2")
library("ggpubr")
library("dplyr")
library("tidyr")
library("ensembldb")


### Load dds
load("./dds_84_qSvs.Rda")
coldata <- colData(dds)
coldata <- as.data.frame(coldata)

# Check the desing formula
design(dds)


#####################################################################################################################
### Deferential expression analysis with SEX, AGE and 5 qSVs as covariates
#####################################################################################################################

#########################
# Wald test #
#########################

dds_W <- DESeq(dds)


### Save results
### Extract results for specific comparisons
res_0vs14 <- results(dds_W, contrast =c("Braak_aSyn_groups","1","0")) # Braak Lewy body stage 1-4 vs Braak Lewy body stage 0
summary(res_0vs14)

res_0vs5 <- results(dds_W, contrast =c("Braak_aSyn_groups","2","0")) # Braak Lewy body stage 5 vs Braak Lewy body stage 0
summary(res_0vs5)

res_0vs6 <- results(dds_W, contrast =c("Braak_aSyn_groups","3","0")) # Braak Lewy body stage 6 vs Braak Lewy body stage 0
summary(res_0vs6)

### subset significant genes FDR < 0.05
# Braak Lewy body stage 1-4 vs Braak Lewy body stage 0
res_0vs14_0.05 <- results(dds_W, contrast =c("Braak_aSyn_groups","1","0"), alpha=0.05)
summary(res_0vs14_0.05)
# Number of DE genes
sum(res_0vs14_0.05$padj < 0.05, na.rm = TRUE) # 0

# Braak Lewy body stage 5 vs Braak Lewy body stage 0
res_0vs5_0.05 <- results(dds_W, contrast =c("Braak_aSyn_groups","2","0"), alpha=0.05)
summary(res_0vs5_0.05)
# Number of DE genes
sum(res_0vs5_0.05$padj < 0.05, na.rm = TRUE) # 979

# Braak Lewy body stage 6 vs Braak Lewy body stage 0
res_0vs6_0.05 <- results(dds_W, contrast =c("Braak_aSyn_groups","3","0"), alpha=0.05)
summary(res_0vs6_0.05)
# Number of DE genes
sum(res_0vs6_0.05$padj < 0.05, na.rm = TRUE) # 38



############ 0 vs 1-4

### Create a tibbl  e for LRT results
res_0vs14_tb <- res_0vs14 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

### Save file with all genes
write.csv(as.data.frame(res_0vs14_tb), file="./res_0vs14_genes.csv",  row.names = FALSE)


############ 0 vs 5

### Create a tibbl  e for LRT results
res_0vs5_tb <- res_0vs5 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

### Save file with all genes
write.csv(as.data.frame(res_0vs5_tb), file="./res_0vs5_genes.csv",  row.names = FALSE)

### Subset to return genes with padj < 0.1
res_0vs5_0.1 <- res_0vs5_tb %>% 
  dplyr::filter(padj < 0.1)

### Create a tibbl  e for LRT results
res_0vs5_0.05_tb <- res_0vs5_0.05 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

### Subset to return genes with padj < 0.05
res_0vs5_0.05 <- res_0vs5_0.05_tb %>% 
  dplyr::filter(padj < 0.05)

### Order genes by pvalue
res_0vs5_0.05 <- res_0vs5_0.05[order(res_0vs5_0.05$pvalue),]
res_0vs5_0.1 <- res_0vs5_0.1[order(res_0vs5_0.1$pvalue),]

### Save file with significant genes
write.csv(as.data.frame(res_0vs5_0.05), file="./res_0vs5_0.05.csv",  row.names = FALSE)
write.csv(as.data.frame(res_0vs5_0.1), file="./res_0vs5_0.1.csv",  row.names = FALSE)


############ 0 vs 6
 
### Create a tibbl  e for LRT results
 res_0vs6_tb <- res_0vs6 %>%
   data.frame() %>%
   rownames_to_column(var="gene") %>% 
   as_tibble()
 
### Save file with all genes
write.csv(as.data.frame(res_0vs6_tb), file="./res_0vs6_genes.csv",  row.names = FALSE)
 
### Subset to return genes with padj < 0.1
res_0vs6_0.1 <- res_0vs6_tb %>% 
   dplyr::filter(padj < 0.1)
 
### Create a tibbl  e for LRT results
res_0vs6_0.05_tb <- res_0vs6_0.05 %>%
   data.frame() %>%
   rownames_to_column(var="gene") %>% 
   as_tibble()
 
### Subset to return genes with padj < 0.05
res_0vs6_0.05 <- res_0vs6_0.05_tb %>% 
   dplyr::filter(padj < 0.05)

### Order genes by pvalue
res_0vs6_0.05 <- res_0vs6_0.05[order(res_0vs6_0.05$pvalue),]
res_0vs6_0.1 <- res_0vs6_0.1[order(res_0vs6_0.1$pvalue),]
 
### Save file with significant genes
write.csv(as.data.frame(res_0vs6_0.05), file="./res_0vs6_0.05.csv",  row.names = FALSE)
write.csv(as.data.frame(res_0vs6_0.1), file="./res_0vs6_0.1.csv",  row.names = FALSE)

 
