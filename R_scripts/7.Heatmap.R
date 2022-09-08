### ### ### ### ### ### ### ### ### ### ###
###   Heatmap  Supplementary figure 3   ###
### ### ### ### ### ### ### ### ### ### ###


### Packages required:
library("DESeq2")
library("limma")
library("ggplot2")
library("ComplexHeatmap")
library("dtwclust") 
library("circlize")


### Load dds
load("./dds_84_qSvs.Rda")
coldata <- colData(dds)

### Transform counts
vsd <- vst(dds, blind = FALSE) 

### adjusting for sex, age at death and 5 qSVs using removeBatchEffect() form limma package
covariates <- coldata[,c(14:18)]

### Sex + Age + Covariates
vsd_adj <- vsd
assay(vsd_adj) <- limma::removeBatchEffect(assay(vsd_adj),batch = vsd_adj$Sex, batch2 = scale(vsd_adj$Age_death), covariates =  as.matrix(covariates),  model.matrix(~ vsd_adj$Braak_aSyn_groups))

##### HEATMAP with significant genes
genes_0.05 <- read.csv2("./sigLRT_genes_0.05.csv", sep = ",")

### Make matrix with only significant genes
matrix_adj <- assay(vsd_adj)
matrix_adj <- matrix_adj[rownames(matrix_adj) %in% genes_0.05$gene,]

# normalize to Z scores 
matrix_adj_Z <- zscore(matrix_adj) # this function normalizes the data by row automatically 
matrix_adj_Z <- t(matrix_adj_Z) # transpose the matrix to plot MS/Ctrl as rows

### Add row names 
rownames(matrix_adj_Z) <- vsd_adj$Braak_aSyn_groups_Plot
matrix_adj_Z <- t(matrix_adj_Z)

### Colors for heatmap
col<- colorRampPalette(colors)(256)
my_colors_ID <- c("0" = "#C6DBEF", "1-4" = "#6BAED6", "5" = "#08519C", "6" = "#08306B")

### Annotation
col_ha <- list(Braak_Stage = my_colors_ID)
column_ha = HeatmapAnnotation(Braak_Stage = vsd_adj$Braak_aSyn_groups_Plot, col = col_ha,
                       annotation_legend_param = list(title= "Braak Stage"),
                       show_annotation_name = FALSE)

png("./Heatmap.png",width=8.5,height=6,units="in",res=600)
Heatmap(matrix_adj_Z, name= "heatmap", col = col,
        bottom_annotation = column_ha,
        clustering_method_rows = "complete",
        show_row_names = FALSE,
        show_column_names = FALSE,
        column_km = 3,
        #column_gap = unit(c(0.5,0.25), "cm"),
        width = unit(6, "in"), 
        height = unit(5, "in"),
        column_dend_height = unit(1.5, "cm"),
        row_dend_width = unit(1.5, "cm"),
        row_gap = unit(0.5, "mm"),
        column_gap = unit(0.5, "mm"),
        column_title = NULL,
        row_title = NULL,
        heatmap_legend_param = list(title = "Z-Score",legend_height = unit(4, "cm")))

dev.off()
 

