###   Covariate Selection


### Packages required:
library("DESeq2")
library("sva")
library("dplyr")
library("rstatix")
library("corrplot")
library("circlize")
library("limma")
library("ggplot2")
library("ggpubr")
library("ComplexHeatmap")

### Load dds
load("./dds_84.Rda")
coldata <- colData(dds)

#####################################################################################################################
### qSVA: correct for RNA quality
#####################################################################################################################

# Import degradation matrix
degMat <- read.delim("./degradation.matrix_84.txt") 
SampleNames <- degMat$X

# Save path to degradation matrix
covFile <- "./degradation.matrix_84.txt"  

# Import table with number of reads per sample and create a vector
Mreads <- read.csv("./MReads_84.csv")
Reads <- Mreads$M.Seqs

# Generate degradation matrix
degCovAdj2 = read.degradation.matrix(
  covFiles = covFile, # coverage file(s) for degradation regions 		-> path to matrix saved as a txt
  sampleNames = SampleNames, # sample names; creates column names of degradation matrix -> rownames of degradation matrix
  readLength = 150, # read length in base pairs
  totalMapped = Reads, # how many reads per sample (library size normalization)  -> you can get the total number of reads either from the first line of the flagstats files or with the following command: samtools view -c /Volumes/Data/RNA-Seq/Results/HISAT2/LINUX/sample_X.bam
  normFactor = 80,  # common library size to normalize to, 80M as default
  type="region_matrix_all") # whether input are individual 'bwtool' output, 'region_matrix' run on individual samples, or 'region_matrix' run on all samples together

### Estimate qSVs
qSVs <- qsva(degCovAdj2)
coldata <- merge(coldata, qSVs, by.x = "row.names", by.y = "row.names")


### Correlation between qSVs and Group, Age , Sex, PMD, RIN, Cells, 5 qSVs
## prepare dataframe
coldata <- coldata %>% mutate(Sex =case_when(Sex == "F" ~ "1", 
                                              Sex == "M" ~ "2"))

coldata <- coldata %>% mutate(Neuro =case_when(Neuropath_diagnosis == "CONTR" ~ "0", 
                                                Neuropath_diagnosis == "iLBD" ~ "1",
                                                Neuropath_diagnosis == "PD" ~ "2",
                                                Neuropath_diagnosis == "PDD" ~ "3"))
## Load cell composition
load("/Volumes/Data/RNA-Seq/Results/Final_09.2021/Braak_Stage_Groups/R_datasets/4_groups_cell.Rda")
coldata$Neurons <- cell$Neurons
coldata$Astrocytes <- cell$Astrocytes
coldata$VLMC <- cell$VLMC
coldata$Tanycytes <- cell$Tanycytes
coldata$OPC <- cell$OPC
coldata$Oligodendrocytes <- cell$Oligodendrocytes
coldata$Ependymal <- cell$Ependymal
coldata$Unknown <- cell$Unknown
coldata$Endothelial <- cell$Endothelial
coldata$Microglia <- cell$Microglia
coldata$NFO <- cell$NFO
coldata <- as.data.frame(coldata)
coldata <- coldata%>% 
  rename(Age = Age_death,
         PMD = PMD_min,
         Group = Braak_aSyn_groups,
         qSV1 = PC1,
         qSV2 = PC2,
         qSV3 = PC3,
         qSV4 = PC4,
         qSV5 = PC5)
coldata_cor$Group <- as.numeric(coldata_cor$Group)
coldata_cor$Sex <- as.numeric(coldata_cor$Sex)



coldata_cor <- coldata[, c("Group","Sex", "Age","PMD","RIN",
                            "Neurons","Astrocytes"  , "VLMC" ,  "Tanycytes" , "OPC" ,"Oligodendrocytes","Ependymal",   "Unknown" ,  "Endothelial" , "Microglia", "NFO",
                            "qSV1" ,  "qSV2" , "qSV3" , "qSV4", "qSV5" )]


## Pearson correlation
Pearson_matrix_1 <-  cor_mat(coldata_cor, method = "pearson")

# get the p.values
Pearson_matrix_p_1 <- Pearson_matrix_1 %>% cor_get_pval()

Pearson_matrix_1 <- as.data.frame(Pearson_matrix_1)
rownames(Pearson_matrix_1) <- Pearson_matrix_1$rowname
Pearson_matrix_1$rowname <- NULL
Pearson_matrix_1 <- mutate_all(Pearson_matrix_1, function(x) as.numeric(as.character(x)))
Pearson_matrix_1 <- as.matrix(Pearson_matrix_1)

Pearson_matrix_p_1 <- as.data.frame(Pearson_matrix_p_1)
rownames(Pearson_matrix_p_1) <- Pearson_matrix_p_1$rowname
Pearson_matrix_p_1$rowname <- NULL
Pearson_matrix_p_1 <- mutate_all(Pearson_matrix_p_1, function(x) as.numeric(as.character(x)))
Pearson_matrix_p_1 <- as.matrix(Pearson_matrix_p_1)


### Plot tau coefficients
# Initialize file path

cor_plot_1 <- { 
  corrplot(Pearson_matrix_1, p.mat = Pearson_matrix_p_1, sig.level = 0.05, type = "lower", diag = FALSE, 
           tl.col = "black", tl.srt = 45, insig='blank', tl.cex = 1)
  recordPlot()
}

cor_plot_1

### Plot p-values
Pearson_matrix_p_log_1 <- -log10(Pearson_matrix_p_1)

### Figure 4a
Pearson_matrix_p_log_plot_1 <- Pearson_matrix_p_log_1[c(5,4,6:16),17:21]

col_fun = colorRamp2(c(0, 1.3, 12), c("grey", "white", "red2"))

png("./P-value_cor_qSVs.png", width=5,height= 6,units="in",res=600)
Hm <- Heatmap(Pearson_matrix_p_log_plot_1, col = col_fun, name = " ", cluster_columns = FALSE, cluster_rows = FALSE, column_names_rot = 0, column_names_centered = 1, row_names_side = "left", rect_gp = gpar(col = "white", lwd = 1), width = unit(3, "in"), 
              height = unit(5, "in"),
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(Pearson_matrix_p_log_plot_1[i, j] > 1.29)
                  grid.text(sprintf("%.1f", Pearson_matrix_p_log_plot_1[i, j]), x, y, gp = gpar(fontsize = 12))
              })
draw(Hm)
Hm_plot <- recordPlot()
dev.off()

print(Hm_plot)
ggsave("./P-value_cor_qSVs.svg", plot = Hm_plot, device = "svg", width = 15, height = 7, units = 'in', dpi =300)
print(Plot)


# save dds with qSVs in col data
dds$Braak_aSyn_groups <- as.factor(dds$Braak_aSyn_groups)
dds$qSV1 <- coldata$qSV1
dds$qSV2 <- coldata$qSV2
dds$qSV3 <- coldata$qSV3
dds$qSV4 <- coldata$qSV4
dds$qSV5 <- coldata$qSV5

design(dds) <- ~ Sex + Age_death + qSV1 + qSV2 + qSV3 + qSV4 + qSV5 + Braak_aSyn_groups
save(dds, file="./dds_84_qSvs.Rda")


#####################################################################################################################
### PCA plots
#####################################################################################################################

coldata <- colData(dds)
coldata <- as.data.frame(coldata)

### variance stabilizing transformation (VST)
vsd <- vst(dds, blind = FALSE) # differences between the variables in the design will not contribute to the expected variance-mean trend of the experiment.

### Perform PCA to visualize data
my_colors_ID <- c("0" = "#C6DBEF", "1-4" = "#6BAED6", "5" = "#08519C", "6" = "#08306B")

  
Groups <- plotPCA(vsd, intgroup = c( "Braak_aSyn_groups_Plot"))
Groups<- Groups  + 
            scale_color_manual(values = my_colors_ID) +
            stat_ellipse(type="norm") +
            theme_minimal() +
            labs(title="Uncorrected", col="Braak Stage") +
            geom_point(shape = 1,size = 3, colour = "black") +
            theme(text = element_text(size = 10),
                  plot.title = element_text(size=12, margin=margin(0,0,0,0)))

print(Groups)


### PCA after adjusting for RIN using removeBatchEffect() form limma package
covariates <- coldata[,c(14:18)]

### Sex + Age + qSvs
vsd_1 <- vsd
assay(vsd_1) <- limma::removeBatchEffect(assay(vsd_1),batch = vsd_1$Sex, batch2 = scale(vsd_1$Age_death), covariates =  as.matrix(covariates),  model.matrix(~ vsd_1$Braak_aSyn_groups))

Plot1 <- plotPCA(vsd_1, intgroup = c( "Braak_aSyn_groups_Plot"))
Plot1 <- Plot1 +  
            scale_color_manual(values = my_colors_ID) +
            stat_ellipse(type="norm", linetype=1) +
            theme_minimal() +
            labs(title="Corrected by Sex, Age at Death and 5 qSVs", col="Braak Stage") +
            geom_point(shape = 1,size = 3, colour = "black") +
            theme(text = element_text(size = 10),
              plot.title = element_text(size=12, margin=margin(0,0,0,0)))

print(Plot1)

### Figure 4b
figure1 <- ggarrange(Groups, Plot1,  
                     labels = c("c","d"),
                     ncol = 1,
                     nrow = 2,
                     common.legend = TRUE, legend = "bottom")
### Figure 4
figure <- ggarrange(Hm_plot, figure1,  
                     ncol = 2,
                     nrow = 1)








