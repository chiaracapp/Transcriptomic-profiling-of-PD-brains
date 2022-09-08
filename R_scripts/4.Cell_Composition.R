###   Cell Composition

### Packages required:
library("rstatix")
library("dplyr")
library("reshape2") 
library("DESeq2")
library("ggplot2")
library("biomaRt")
library("markerGeneProfile")
library("homologene")
library("ggpubr")
library("corrplot")

### Load dds
load("./dds_84.Rda")
coldata <- colData(dds)


#####################################################################################################################
### Cell deconvolution using Scaden 
#####################################################################################################################

### Steps in Scaden:
###   1. Process: Scaden performs pre-processing on your training data, making sure it has the same genes as your prediction (bulk) data and performing some data transformations to make the data suitable for machine learning.
###               For this step you need the training data downloaded (mouse brain) and the bulk data in the form of counts normalized by library size.
###               Scaden will create a new file for training which only contains the intersection of genes between the training and the prediction data. Furthermore, the training data will be log2-transformed and scaled to the range [0,1].
###   2. Train:  Scaden consists of three deep neural network models. By default, each of them will be trained for 5,000 steps, which is the recommended number of training steps.
###   3. Predict: Scaden predicts the cell composition


## Normalize counts to library size using the median ratio method described in equation 5 in Anders and Huber (2010)
dds_Scaden <- estimateSizeFactors(dds)

## Export counts from dds for cell estimation with scaden
counts_Scaden <- counts(dds_Scaden, normalized = TRUE)
counts_Scaden <- as.data.frame(counts_Scaden)

# Add column with gene symbols corresponding to ensembl ID
listMarts()
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

genes.table <- getBM(filters= "ensembl_gene_id", 
                    attributes= c("ensembl_gene_id","hgnc_symbol", "description","gene_biotype"),
                    values= rownames(counts_Scaden), mart= ensembl)

# merge gene symbols to counts
counts_gene <- merge(counts_Scaden, genes.table, by.x = 'row.names', by.y = 'ensembl_gene_id')

# organize columns, annotations in first columns 
counts_gene <- counts_gene[,c(1,86:88, 2:85)]

# delete unwanted columns (Keep only gene symbols)
counts_gene <- counts_gene[,-c(1,3,4)]

# convert to mouse names
human_genes <- counts_gene[,1]

# Basic function to convert human to mouse gene names
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

genesV2 = getLDS(attributes = c("hgnc_symbol"), 
                 filters = "hgnc_symbol", 
                 values = human_genes , 
                 mart = human, 
                 attributesL = c("mgi_symbol"), 
                 martL = mouse, 
                 uniqueRows=T) 

genesV3 <- genesV2[!duplicated(genesV2[,c(1)]),]

counts_mouse <- merge(counts_gene, genesV3, by.x = "hgnc_symbol", by.y = "HGNC.symbol") # 14.603 genes
counts_mouse <- counts_mouse[,-c(1)] # remove column with human symbols
counts_mouse <- counts_mouse[,c(85, 1:84)] # move column with mouse symbols in front

write.table(counts_mouse, "./Counts_mouse_84.txt",  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

### Run Scaden

#      ____                _            
#     / ___|  ___ __ _  __| | ___ _ __  
#     \___ \ / __/ _` |/ _` |/ _ \ '_ \ 
#      ___) | (_| (_| | (_| |  __/ | | |
#     |____/ \___\__,_|\__,_|\___|_| |_|

# scaden, version 1.1.2

# (base) chiaraca@L25087 84_samples % scaden process mouse_brain.h5ad Counts_mouse_84.txt

# (base) chiaraca@L25087 84_samples %  scaden train processed.h5ad  --steps 5000 --model_dir model                                                                                                                                                                                                                                 train.py:98

# (base) chiaraca@L25087 84_samples % scaden predict --model_dir model Counts_mouse_84.txt

# a file called scaden_predictions.txt with the percentage of each cell type in each sample is obtained


### Comparison in cell composition between groups

cell <- read.delim("./scaden_predictions.txt", header = TRUE, sep = "\t", dec = ".")

# add columns with braak stage and neuropath
cell <- merge(cell, coldata, by.x = "X", by.y = "row.names")
cell <- as.data.frame(cell)
# remove unwanted columns
cell <- cell[,-c(17:20,22:23,25)]
rownames(cell) <- cell[,"X"]
cell[,"X"] <- NULL


cell <- cell %>% mutate(Sex = case_when(Sex == "F" ~ "1", 
                                        Sex == "M" ~ "2"))

cell <- cell %>% mutate(Diagnosis = case_when(Neuropath_diagnosis == "CONTR" ~ "0", 
                                              Neuropath_diagnosis == "iLBD" ~ "1",
                                              Neuropath_diagnosis == "PD" ~ "2",
                                              Neuropath_diagnosis == "PDD" ~ "3"))

cell$Neuropath_diagnosis <- NULL
cell <- mutate_all(cell, function(x) as.numeric(as.character(x)))

cell <- cell %>% 
      rename(Age = Age_death,
             PMD = PMD_min,
             Group = Braak_aSyn_groups)

cell <- cell[, c("Group","Diagnosis", "Sex", "Age", "PMD","RIN", 
                       "Neurons","Astrocytes"  , "VLMC" ,  "Tanycytes" , "OPC" ,"Oligodendrocytes","Ependymal",   "Unknown" ,  "Endothelial" , "Microglia", "NFO")]

cell <- mutate_all(cell, function(x) as.numeric(as.character(x)))

save(cell, file="./4_groups_cell.Rda")

### Pearson correlation
Pearson_matrix <-  cor_mat(cell, method = "pearson")

# get the p.values
Pearson_matrix_p <- Pearson_matrix %>% cor_get_pval()

# Arrange for plotting
Pearson_matrix <- as.data.frame(Pearson_matrix)
rownames(Pearson_matrix) <- Pearson_matrix$rowname
Pearson_matrix$rowname <- NULL
Pearson_matrix <- mutate_all(Pearson_matrix, function(x) as.numeric(as.character(x)))
Pearson_matrix <- as.matrix(Pearson_matrix)

Pearson_matrix_p <- as.data.frame(Pearson_matrix_p)
rownames(Pearson_matrix_p) <- Pearson_matrix_p$rowname
Pearson_matrix_p$rowname <- NULL
Pearson_matrix_p <- mutate_all(Pearson_matrix_p, function(x) as.numeric(as.character(x)))
Pearson_matrix_p <- as.matrix(Pearson_matrix_p)

### Figure 3a
### Plot tau coefficients
# Initialize file path
png("./Cell_pearson.png",width=4.5,height=6,units="in",res=600)
  corrplot(Pearson_matrix, p.mat = Pearson_matrix_p, sig.level = 0.05, type = "lower", diag = FALSE, 
           tl.col = "black", tl.srt = 45, insig='blank', tl.cex = 0.7, mar = c(0.5,0.5,0.5,0), cl.cex = 0.7)
Plot <- recordPlot()
print(Plot)
dev.off()
  

#### Prepare dataset for ancova

cell_2 <- cell
cell_2$Group <- as.factor(cell_2$Group)
cell_2$Sex <- as.factor(cell_2$Sex)

######## ANCOVA #############

res.Neurons <- cell_2 %>% anova_test(Neurons ~ Sex + Age + PMD + RIN + Group)
get_anova_table(res.Neurons)

res.Astrocytes <- cell_2 %>% anova_test(Astrocytes ~ Sex + Age + PMD + RIN + Group)
get_anova_table(res.Astrocytes)

res.VLMC <- cell_2 %>% anova_test(VLMC ~ Sex + Age + PMD + RIN + Group)
get_anova_table(res.VLMC)

res.Tanycytes <- cell_2 %>% anova_test(Tanycytes ~ Sex + Age + PMD + RIN + Group)
get_anova_table(res.Tanycytes)

res.OPC <- cell_2 %>% anova_test(OPC ~ Sex + Age + PMD + RIN + Group)
get_anova_table(res.OPC)

res.Oligodendrocytes <- cell_2 %>% anova_test(Oligodendrocytes ~ Sex + Age + PMD + RIN + Group)
get_anova_table(res.Oligodendrocytes)

res.Ependymal <- cell_2 %>% anova_test(Ependymal ~ Sex + Age + PMD + RIN + Group)
get_anova_table(res.Ependymal)

res.Unknown <- cell_2 %>% anova_test(Unknown ~ Sex + Age + PMD + RIN + Group)
get_anova_table(res.Unknown)

res.Endothelial <- cell_2 %>% anova_test(Endothelial ~ Sex + Age + PMD + RIN + Group)
get_anova_table(res.Endothelial)

res.Microglia <- cell_2 %>% anova_test(Microglia ~ Sex + Age + PMD + RIN + Group)
get_anova_table(res.Microglia)

### Plot cell composition Figure 3b
# reshape the dataset
cell_3 <- cell_2
cell_3$Sex <- NULL
cell_3$Age <- NULL
cell_3$PMD <- NULL
cell_3$RIN <- NULL
cell_3$Diagnosis <- NULL


cell_plot <- aggregate(. ~ Group, cell_3, function(x) list(mean = round(mean(x),3), sd =round(sd(x),3)))
cell_plot <- do.call("data.frame", cell_plot) # flatten

Neurons <- cell_plot[,1:3]
Neurons$cell_type <- c(rep("Neurons", 4))
names(Neurons) <- c("Braak_aSyn_groups", "Mean", "SD", "Cell_Type")
Astrocytes <- cell_plot[,c(1,4,5)]
Astrocytes$cell_type <- c(rep("Astrocytes", 4))
names(Astrocytes) <- c("Braak_aSyn_groups", "Mean", "SD", "Cell_Type")
VLMC <- cell_plot[,c(1,6,7)]
VLMC$cell_type <- c(rep("VLMC", 4))
names(VLMC) <- c("Braak_aSyn_groups", "Mean", "SD", "Cell_Type")
Tanycytes <- cell_plot[,c(1,8,9)]
Tanycytes$cell_type <- c(rep("Tanycytes", 4))
names(Tanycytes) <- c("Braak_aSyn_groups", "Mean", "SD", "Cell_Type")
OPC <- cell_plot[,c(1,10,11)]
OPC$cell_type <- c(rep("OPC", 4))
names(OPC) <- c("Braak_aSyn_groups", "Mean", "SD", "Cell_Type")
Oligodendrocytes <- cell_plot[,c(1,12,13)]
Oligodendrocytes$cell_type <- c(rep("Oligodendrocytes", 4))
names(Oligodendrocytes) <- c("Braak_aSyn_groups", "Mean", "SD", "Cell_Type")
Ependymal <- cell_plot[,c(1,14,15)]
Ependymal$cell_type <- c(rep("Ependymal", 4))
names(Ependymal) <- c("Braak_aSyn_groups", "Mean", "SD", "Cell_Type")
Unknown <- cell_plot[,c(1,16,17)]
Unknown$cell_type <- c(rep("Unknown", 4))
names(Unknown) <- c("Braak_aSyn_groups", "Mean", "SD", "Cell_Type")
Endothelial <- cell_plot[,c(1,18,19)]
Endothelial$cell_type <- c(rep("Endothelial", 4))
names(Endothelial) <- c("Braak_aSyn_groups", "Mean", "SD", "Cell_Type")
Microglia <- cell_plot[,c(1,20,21)]
Microglia$cell_type <- c(rep("Microglia", 4))
names(Microglia) <- c("Braak_aSyn_groups", "Mean", "SD", "Cell_Type")
NFO <- cell_plot[,c(1,20,21)]
NFO$cell_type <- c(rep("NFO", 4))
names(NFO) <- c("Braak_aSyn_groups", "Mean", "SD", "Cell_Type")


cellplot <- rbind(Neurons,Astrocytes,VLMC,Tanycytes,OPC,Oligodendrocytes,Ependymal,Unknown,Endothelial,Microglia,NFO)

cellplot <- cellplot %>% mutate(Braak_aSyn_group_Plot =
                                case_when(Braak_aSyn_groups == "0" ~ "0", 
                                          Braak_aSyn_groups == "1" ~ "1-4", 
                                          Braak_aSyn_groups == "2" ~ "5",
                                          Braak_aSyn_groups == "3" ~ "6" ))

cellplot$Mean <- as.numeric(cellplot$Mean)
cellplot$SD <- as.numeric(cellplot$SD)

my_colors_ID <- c("0" = "#C6DBEF", "1-4" = "#6BAED6", "5" = "#08519C", "6" = "#08306B")

# make stacked bar chart
Plot1 <- ggplot(cellplot, aes(x= reorder(Cell_Type, Mean),y=Mean,fill=factor(Braak_aSyn_group_Plot)))+
            geom_bar(stat="identity",position="dodge",color = "black") +
            theme_minimal() +
            scale_fill_manual(values = my_colors_ID) + 
            scale_x_discrete(guide = guide_axis(n.dodge = 1)) +
            geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2,position=position_dodge(.9)) +
            coord_flip() +
            labs( x = NULL, y = "Cell-type %", fill = "Braak Stage") +
            theme(plot.title = element_text(size=12, margin=margin(0,0,0,0)))


# Figure 3
figure <- ggarrange(Plot, Plot1, 
                     labels = c("a","b"),
                     ncol = 2,
                     nrow = 1)

ggsave("./Cells.svg", plot = figure, device = "svg", width = 8.5, height = 8, units = 'in', dpi =600)
