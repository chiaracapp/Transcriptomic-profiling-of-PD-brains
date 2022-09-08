### ### ### ### ### ### ### 
###   Selected eQTLs    ###
### ### ### ### ### ### ### 


### Packages required:
library("DESeq2")
library("dplyr")
library("rstatix")
library("corrplot")
library("circlize")
library("limma")
library("ggplot2")
library("r2symbols")
library("ggpubr")

### Load dds
load("./dds_84_qSvs.Rda")
coldata <- colData(dds)
coldata <- as.data.frame(coldata)

### Load PD risk genes
Risk_genes <- read.csv2("./Nalls_2019_RiskSNP_Extended_CC.csv", sep = ",")

counts <- counts(dds)
genes_all <- rownames(counts)

# Delete duplicated genes
Risk_genes <- Risk_genes[!duplicated(Risk_genes$Nearest.Gene),]

# Delete NA
Risk_genes <- Risk_genes[!is.na(Risk_genes$Ensembl.ID),] 

# Keep only genes detected with RNA-Seq
Risk_genes_T <- Risk_genes[Risk_genes$Ensembl.ID %in% genes_all,] 

### Load significant genes
genes_0vs5_0.05 <- read.csv2("./res_0vs5_0.05_ID.csv", sep = ",")

### Is any of the risk genes among the DE genes?
Risk_genes_0vs5 <- Risk_genes_T[Risk_genes_T$Nearest.Gene %in% genes_0vs5_0.05$hgnc_symbol,]

# Save
write.csv(as.data.frame(Risk_genes_0vs5), file="./Risk_genes_0vs5_0.05.csv",  row.names = FALSE)

###### eQTL analysis for PD risk genes DE at Braak Lewy body stage 5

## genotype file as YG IDs, add them to coldata
YG_IDs <- read.csv2("./Sample_ID_New.csv", sep = ",")
YG_IDs <- YG_IDs[YG_IDs$X %in% rownames(coldata),]

coldata <- merge(coldata, YG_IDs, by.x = "row.names", by.y = "X")
rownames(coldata) <- coldata$Row.names
coldata$Row.names <- NULL


### Add position ID (Xchr.position.A1.A2_RiskA) to rsID SNPs
SNPs_ID <- c("X1.154898185.G.C_C", "X2.102396963.T.A_A", "X4.77110365.C.G_G", "X5.134199105.C.A_A", "X16.52969426.G.A_A")
Risk_genes_0vs5$SNPs_ID <-  SNPs_ID

### Load genotypes
Genotypes <- read.csv2("./Chiara_samples_PDsnps_new.txt", sep = "")

### Keep only 84 samples
Genotypes <- Genotypes[Genotypes$FID %in% coldata$Sample_Name,] # 1 samples is missing YG498
rownames(Genotypes) <- Genotypes$FID

### Keep only SNPs of interest in genotype
Genotypes <- Genotypes[,colnames(Genotypes) %in% SNPs_ID] 


# Transform counts
vsd <- vst(dds, blind = FALSE) 

# Adjust counts for covariates
covariates <- coldata[,c(14:18)] # 5 qSVs
assay(vsd) <- limma::removeBatchEffect(assay(vsd),batch = vsd$Sex, batch2 = scale(vsd$Age_death), covariates =  as.matrix(covariates),  model.matrix(~ vsd$Braak_aSyn_groups))
assay(vsd)

Adjusted_counts <- assay(vsd)

### Keep only 5 genes
Adjusted_counts <- Adjusted_counts[rownames(Adjusted_counts) %in% Risk_genes_0vs5$Ensembl.ID,]

### Divide samples in 4 groups
Adjusted_counts <- t(Adjusted_counts)
Adjusted_counts <- as.data.frame(Adjusted_counts)
Adjusted_counts$Group <- coldata$Braak_aSyn_groups_Plot
Adjusted_counts$YGID <- coldata$Sample_Name
Adjusted_counts <- Adjusted_counts[!(Adjusted_counts$YGID == "YG498"),] # remove sample without genotype YG498

Adjusted_counts <- merge(Adjusted_counts, Genotypes, by.x = "YGID", by.y = "row.names")

### Plot gene by corresponding SNP + Linear regression
Adjusted_counts$X1.154898185.G.C_C <- factor(Adjusted_counts$X1.154898185.G.C_C)
plot(Adjusted_counts$X1.154898185.G.C_C, Adjusted_counts$ENSG00000163344,)
lm_1 <- lm(ENSG00000163344 ~ X1.154898185.G.C_C, data=Adjusted_counts)
summary(lm_1) # p-value: 0.738 (no GG present)

Adjusted_counts$X2.102396963.T.A_A <- factor(Adjusted_counts$X2.102396963.T.A_A)
plot(Adjusted_counts$X2.102396963.T.A_A, Adjusted_counts$ENSG00000071054,)
lm_2 <- lm(ENSG00000071054 ~ X2.102396963.T.A_A, data=Adjusted_counts)
summary(lm_2) # p-value: 0.005071 ######

Adjusted_counts$X4.77110365.C.G_G <- factor(Adjusted_counts$X4.77110365.C.G_G)
plot(Adjusted_counts$X4.77110365.C.G_G, Adjusted_counts$ENSG00000138760,)
lm_5 <- lm(ENSG00000138760 ~ X4.77110365.C.G_G, data=Adjusted_counts)
summary(lm_5) # p-value: 0.8156

Adjusted_counts$X5.134199105.C.A_A <- factor(Adjusted_counts$X5.134199105.C.A_A)
plot(Adjusted_counts$X5.134199105.C.A_A, Adjusted_counts$ENSG00000181904,)
lm_9 <- lm(ENSG00000181904 ~ X5.134199105.C.A_A, data=Adjusted_counts)
summary(lm_9) # p-value: 0.8825 (no CC)

Adjusted_counts$X16.52969426.G.A_A <- factor(Adjusted_counts$X16.52969426.G.A_A)
plot(Adjusted_counts$X16.52969426.G.A_A, Adjusted_counts$ENSG00000177200,)
lm_8 <- lm(ENSG00000177200 ~ X16.52969426.G.A_A, data=Adjusted_counts)
summary(lm_8) # p-value:  0.3063


 p_values <- c(0.738, 0.005071, 0.8156,  0.8825, 0.3063 )
 p_adj <- p.adjust(p_values, method = "BH", n = length(p_values))
 p_adj # 0.882500 0.025355 0.882500 0.882500 0.765750
 


#### Figure 7 a
#### MAP4K4
summary(lm_2)

MAP4K4_All <- ggplot(Adjusted_counts, aes(x=as.factor(X2.102396963.T.A_A), y=ENSG00000071054, fill=as.factor(X2.102396963.T.A_A))) + 
  geom_boxplot() +
  geom_point(size = 2) +
  scale_fill_manual(values = c("snow3","snow3","snow3")) + 
  #geom_abline(intercept =13.58042, slope = -0.19015 , color="#9E0142" )+
  scale_x_discrete(labels=c("TT","TA","AA")) +
  labs(title= expression(paste(italic("MAP4K4"), "")),  
       y= "Adjusted Counts", 
       x= "rs11683001") +
  theme_minimal() +
  theme(legend.position = "none", plot.title = element_text(size=12, margin=margin(0,0,0,0)), text = element_text(size = 10))

# Save
ggsave("./MAP4K4_All.svg", plot = MAP4K4_All, device = "svg", width = 4, height = 4, units = 'in', dpi =600)


