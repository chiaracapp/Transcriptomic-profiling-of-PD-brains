###   Collapsing transcript level-quantification obtained with Salmon onto gene-level quantification using the tximport R package and pre-filtering  ###

### Packages required:
library("tximport") 
library("biomaRt")
library("ggpubr") 
library("rstatix")
library("dplyr")
library("DESeq2")
library("ggplot2")
library("ensembldb")
library("AnnotationHub")


### Create a named vector files pointing to the quantification files
Files <- file.path("/Volumes/Data/RNA-Seq/Results/Salmon2", list.files("/Volumes/Data/RNA-Seq/Results/Salmon2"), "quant.sf")
names(Files) <- list.files("/Volumes/Data/RNA-Seq/Results/Salmon2")
Files

### Create tx2gene dataframe using Ensembl release v101 transcriptome which corresponds to Gencode v35 transcriptome

### Prepare library package for v101
## Load the annotation resource.
ah <- AnnotationHub()
ahEdb <- query(ah, pattern = c("Homo Sapiens", "EnsDb", 101))
ahEdb <- ahEdb[[1]]

## Create tx2gene dataframe
tx2gene <- transcripts(ahEdb, columns = c("tx_id", "gene_id"),
                      return.type = "DataFrame")

## Save dataframe with tx_id=enst and gene_id=ensg as csv file
write.csv(tx2gene, "./gene_map101.csv")


### Run tximport
count_data <- tximport(files = Files, 
                       type = "salmon", 
                       tx2gene = tx2gene, 
                       ignoreTxVersion = TRUE)          



### Open file with samples´ information and organize it
coldata <- read.csv("./Coldata_Final_84.csv", row.names=1, sep= ",", stringsAsFactors=FALSE)
coldata$Sample_Name <- NULL
coldata$Brain_surgery <- NULL
col_fact <- c("Sex","Neuropath_diagnosis","Braak_aSyn_stage","Braak_NFT_stage","Instrument","Run.Number","Flow.cell.ID" )
coldata[col_fact] <- lapply(coldata[col_fact], factor)
coldata <- coldata[, c("Sex","Age_death","RIN","PMD_min", "Neuropath_diagnosis","Braak_aSyn_stage","Braak_NFT_stage","Instrument","Run.Number","Flow.cell.ID")]
coldata <- as.data.frame(coldata)

## Add column with neuropathological groups based on Braak Lewy body stage
coldata <- coldata %>% mutate(Braak_aSyn_groups =
                                case_when(Braak_aSyn_stage == "0" ~ "0", 
                                          Braak_aSyn_stage == "1" ~ "1", 
                                          Braak_aSyn_stage == "2" ~ "1", 
                                          Braak_aSyn_stage == "3" ~ "1", 
                                          Braak_aSyn_stage == "4" ~ "1", 
                                          Braak_aSyn_stage == "5" ~ "2",
                                          Braak_aSyn_stage == "6" ~ "3" ))

coldata$Braak_aSyn_groups <- as.factor(coldata$Braak_aSyn_groups)

coldata <- coldata %>% mutate(Braak_aSyn_groups_Plot =
                                case_when(Braak_aSyn_groups == "0" ~ "0", 
                                          Braak_aSyn_groups == "1" ~ "1-4", 
                                          Braak_aSyn_groups == "2" ~ "5",
                                          Braak_aSyn_groups == "3" ~ "6" ))

coldata$Braak_aSyn_groups_Plot <- as.factor(coldata$Braak_aSyn_groups_Plot)


### Construct the DESeqDataSet object from the matrix of counts and the sample information table 
dds <- DESeqDataSetFromTximport(txi = count_data,
                                colData = coldata,
                                design = ~  Braak_aSyn_groups) 

dim(dds) # 60237 genes and 84 samples


### Filter out transcripts encoded by the mitochondrial genome

## Get mitochondrial genes
mito <- genes(ahEdb, filter = ~ seq_name == "MT")
mitoIDs <- names(mito) # 37 mitochondrial transcripts

## Obtain the indices of only desired genes
mitoIDs <- which(!rownames(dds) %in% mitoIDs)

## Cut your desired genes in the DESeq object
dds <- dds[mitoIDs, ]

dim(dds) # 60200 genes and 84 samples -> 37 mitocondrial genes have been removed

### Pre-filtering
## Perform a minimal pre-filtering to keep only rows that have more than 0 counts in all samples (coutns>0)
keep <- counts(dds) > 0 
keep1 <- keep[apply(keep,1, function(x) all(x!="FALSE")),]
keep2 <- rownames(keep1)
dds <- dds[keep2,]

dim(dds) # 21.615 genes and 84 samples -> 38.570 genes have been removed because they had at least 1 sample with 0 reads

### save dds
save(dds, file="./dds_84.Rda")

 
  