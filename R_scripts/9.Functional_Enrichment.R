### ### ### ### ### ### ### ### ### ### 
###   Functional enrichemnt analysis  ###
### ### ### ### ### ### ### ### ### ### 


### Packages required:
library("ggplot2")
library("dplyr")
library("ermineR")
library("colorspace")
library("qdapTools")
library("irr")
library("grid")
library("RColorBrewer")
library("kableExtra")
library("treemapify")


#### Needed functions obtained from Nido, G. S., Dick, F., Toker, L., Petersen, K., Alves, G., Tysnes, O. B., Jonassen, I., Haugarvoll, K., & Tzoulis, C. (2020). 
#### Common gene expression signatures in Parkinson's disease are driven by changes in cell composition. Acta Neuropathol Commun, 8(1), 55. https://doi.org/10.1186/s40478-020-00932-7 

## 
GetAnnoFiles <- function(platform){
  platform <- as.character(platform)
  if(length(list.files(pattern=platform)) == 0){
    download.file(paste0("https://gemma.msl.ubc.ca/annots/", platform, "_noParents.an.txt.gz"), destfile=paste0(platform, ".gz"))
  }
  if(length(list.files(pattern=platform)) > 1){
    print("Multiple annotation files matching the platform exist")
  } else {
    warning("Using existing annotation file, consider updating")
  }
  Anno_file <- read.table(paste0(platform, ".gz"), comment="#", header=T, quote='"', sep="\t")
  return(Anno_file)
}

##
gene_pathway_matrix <- function(input){
  m <- table(input[[2]], input[[1]]) == 1
  m[,match(unique(input[[1]]), colnames(m))]
}

##
overlap_matrix <- function (gmmat){
  #require(parallel)
  Overlap <- sapply(1:ncol(gmmat), function(p1) {
    sapply(1:ncol(gmmat), function(p2) {
      sum(gmmat[,p1] & gmmat[,p2]) / min(sum(gmmat[,p1]), sum(gmmat[,p2]))
    })
  })
  rownames(Overlap) <- colnames(gmmat)
  colnames(Overlap) <- colnames(gmmat)
  return(Overlap)
}

##
update_simmat <- function(simmat, gmmat, failed, passed){
  diag(simmat) <- 1
  #gene_pathway_mat <- gene_pathway_matrix(pathways)
  genes_in_cluster <- rep(0, length(rownames(gmmat)))
  names(genes_in_cluster) <- rownames(gmmat)
  # genes in the union of the two pathways
  merged <- unique(c(rownames(gmmat)[gmmat[,failed]],
                     rownames(gmmat)[gmmat[,passed]]))
  genes_in_cluster <- ifelse(names(genes_in_cluster) %in% merged, 1, 0)
  names(genes_in_cluster) <- rownames(gmmat)
  new_kappa_col <- sapply(1:nrow(simmat), function(x){
    irr::kappa2(data.frame(gmmat[,x],genes_in_cluster))$value
  })
  simmat[,passed] <- new_kappa_col
  simmat <- simmat[-which(rownames(simmat)==failed), -which(rownames(simmat)==failed)]
  diag(simmat) <- NA
  return(simmat)
}

##
cluster_pathways <- function(pathways, threshold=0.4, subsetsize=length(unique(pathways[[1]])),
                             lowThreshold=50, highThreshold=1000, method=c("kappa", "jaccard", "overlap")) {
  # check pathways input is data.frame
  stopifnot("data.frame" %in% class(pathways))
  # max number of pathways allowed are 1000, safety measure
  totPaths <- length(unique(pathways[[1]]))
  if (subsetsize > totPaths) {
    warning(paste0("Subset size chosen (", subsetsize, ") greater than the number of pathways (", totPaths, "), setting subsetsize <- ", totPaths))
    subsetsize <- totPaths
  }
  if (subsetsize > 1000){
    warning("Maximum subset size allowed is 1000, setting subsetsize <- 1000")
    subsetsize <- 1000
  }
  
  method <- method[1]
  
  # Filter to only top pathways
  pathways[,1] <- as.character(pathways[,1])
  pathways[,2] <- as.character(pathways[,2])
  pathways <- pathways[pathways[[1]] %in% head(unique(pathways[[1]]), n=subsetsize),]
  
  # membership of genes to pathways and pathway sizes
  membership_matrix <- gene_pathway_matrix(pathways)
  pSizes <- colSums(membership_matrix)
  
  # similarity of pathways
  if (method == "kappa") {
    simmat <- sim_matrix(membership_matrix)
  } else if (method == "jaccard") {
    simmat <- jaccard_matrix(membership_matrix)
  } else if (method == "overlap") {
    simmat <- overlap_matrix(membership_matrix)
  } else {
    stop(paste0("Method \"", method, "\" unknown"))
  }
  diag(simmat) <- NA
  
  pathwayNames <- colnames(simmat)   # pathway names
  if (any(pathwayNames != unique(pathways[[1]]))) stop("HORRIBLE ERROR")
  
  cluster <- matrix(rep(0, subsetsize*subsetsize), nrow=subsetsize)
  diag(cluster) <- 1
  colnames(cluster) <- pathwayNames
  rownames(cluster) <- pathwayNames
  
  iter <- 1
  while(max(simmat[upper.tri(simmat, diag=FALSE)]) >= threshold){
    message(paste0("Iteration ", iter, "..."))
    index <- which(simmat == max(simmat[upper.tri(simmat, diag = FALSE)]), arr.ind=TRUE)
    t_i = index[1,1] # first row index
    t_j = index[1,2] # first col index
    t_both <- c(t_i, t_j)
    pair_size <- c(pSizes[rownames(simmat)[t_i]], pSizes[rownames(simmat)[t_j]])
    pair_names <- c(rownames(simmat)[t_i], rownames(simmat)[t_j])
    if ( (pair_size[1] <= lowThreshold || pair_size[1] >= highThreshold) && (pair_size[2] <= lowThreshold || pair_size[2] >= highThreshold) ) {
      failed <- rownames(simmat)[max(t_both)]
      passed <- rownames(simmat)[min(t_both)]
      simmat <- update_simmat(simmat, membership_matrix, failed, passed)
      cluster[failed,passed] <- 1
      cluster[failed,failed] <- -1
    } else if ( pair_size[1] <= lowThreshold || pair_size[1] >= highThreshold ) {
      failed <- rownames(simmat)[t_i]
      passed <- rownames(simmat)[t_j]
      simmat <- update_simmat(simmat, membership_matrix, failed, passed)
      cluster[failed,passed] <- 1
      cluster[failed,failed] <- -1
    } else if ( pair_size[2] <= lowThreshold || pair_size[2] >= highThreshold ) {
      failed <- rownames(simmat)[t_j]
      passed <- rownames(simmat)[t_i]
      simmat <- update_simmat(simmat, membership_matrix, failed, passed)
      cluster[failed,passed] <- 1
      cluster[failed,failed] <- -1
    } else {
      failed <- rownames(simmat)[max(t_both)]
      passed <- rownames(simmat)[min(t_both)]
      simmat <- update_simmat(simmat, membership_matrix, failed, passed)
      cluster[failed,passed] <- 1
      cluster[failed,failed] <- -1
    }
    iter <- iter + 1
  }
  
  # remodel "cluster" structure to extract the titles
  cluster_idx <- which(diag(cluster)==1)
  clusters <- apply(cluster[,cluster_idx], 2, function(x) names(which(x==1)))
  data.frame(title=rep(names(clusters), sapply(clusters, length)),
             member=unlist(clusters,  FALSE),
             stringsAsFactors=FALSE,
             row.names=NULL)
}



##
draw_treemap_F <- function(dtf, title="", vp=NULL) {
  treemap::treemap(
    dtf=dtf,
    index=c("Name"), 
    vSize="NumGenes", 
    vColor="log10P",
    sortID="color",
    range = c(-40,40),
    title=title,
    type="value",
    palette=rev(brewer.pal(10, "RdBu")),
    border.col=c("white"),
    border.lwds=c(2),
    overlap.labels=0,
    bg.labels=0,
    fontsize.labels = 8,
    fontfamily.title = "sans",
    inflate.labels = FALSE,
    force.print.labels=TRUE,
    fontsize.title = 12,
    position.legend = "none",
    vp=NULL)
}



### Load genes
all_genes_0vs5 <- read.csv2("./res_0vs5_genes_ID.csv", sep = ",")
all_genes_0vs6 <- read.csv2("./res_0vs6_genes_ID.csv", sep = ",")
all_genes_0vs14 <- read.csv2("./res_0vs14_genes_ID.csv", sep = ",")

## Get annotation file linking every gene to a GO ID
GenericHumanAnno <- GetAnnoFiles("Generic_human")
## Save Annotation
save(GenericHumanAnno, file="./GenericHumanAnno.Rda")

## If re running script, use annotation file saved above
load("./GenericHumanAnno.Rda")


######## We perform enrichment analyses after transforming the p-values to take into account the direction of change in expression


##### 0 vs 1-4
#################

# Prepare dataset
str(all_genes_0vs14)
all_genes_0vs14$log2FoldChange <- as.numeric(all_genes_0vs14$log2FoldChange)
all_genes_0vs14$pvalue <- as.numeric(all_genes_0vs14$pvalue)

# Transform p-values
all_genes_0vs14$DownPval <- apply(all_genes_0vs14 %>% dplyr::select(log2FoldChange, pvalue), 1, function(x){
  if(x[1] < 0){
    x[2]/2
  } else {
    1-x[2]/2
  }
})
all_genes_0vs14$DownPvalAdj <- p.adjust(all_genes_0vs14$DownPval, "BH")

all_genes_0vs14$UpPval <- apply(all_genes_0vs14 %>% dplyr::select(log2FoldChange, pvalue), 1, function(x){
  if(x[1] > 0){
    x[2]/2
  } else {
    1-x[2]/2
  }
})
all_genes_0vs14$UpPvalAdj <- p.adjust(all_genes_0vs14$UpPval, "BH")


### Find pathways for genes down-regulated from Braak Lewy body stage 0 to Braak Lewy body stages 1-4

# scores = A data.frame. Rownames have to be gene identifiers (eg. probes, must be unique), followed by any number of columns. 
all_genes_0vs141 <- all_genes_0vs14[!(is.na(all_genes_0vs14$hgnc_symbol) | all_genes_0vs14$hgnc_symbol==""), ] 
all_genes_0vs141 <- all_genes_0vs141[!duplicated(all_genes_0vs141$hgnc_symbol), ] 
rownames(all_genes_0vs141) <- all_genes_0vs141$hgnc_symbol


EnrichListPDdown_0vs14 <- gsr(scores=all_genes_0vs141, scoreColumn="DownPval",
                        bigIsBetter=FALSE, logTrans=TRUE, annotation=GenericHumanAnno, 
                        aspects=c("Biological Process"),
                        iterations=200000)

EnrichListPDdown_Res_0vs14 <- EnrichListPDdown_0vs14$results %>%
  mutate(Direction="DOWN") %>%
  arrange(Pval)


EnrichListPDup_0vs14 <-gsr(scores=all_genes_0vs141, scoreColumn="UpPval",
                     bigIsBetter=FALSE, logTrans=TRUE, annotation=GenericHumanAnno, 
                     aspects=c("Biological Process"),
                     iterations=200000)

EnrichListPDup_Res_0vs14 <- EnrichListPDup_0vs14$results %>%
  mutate(Direction="UP") %>%
  arrange(Pval)

res_0vs14 <- bind_rows(EnrichListPDdown_Res_0vs14, EnrichListPDup_Res_0vs14) %>%
  arrange(CorrectedPvalue)

# Save 
write.csv(as.data.frame(res_0vs14), file="./ermine_all_0vs14.csv",  row.names = FALSE)

### Subset only significant pathways
res_0vs14_0.05 <- res_0vs14 %>%
  dplyr::filter(CorrectedPvalue < 0.05) %>%
  arrange(CorrectedPvalue) 

### Count pathways
sum(res_0vs14_0.05$Direction == "UP", na.rm = TRUE)
sum(res_0vs14_0.05$Direction == "DOWN", na.rm = TRUE)

write.csv(as.data.frame(res_0vs14_0.05), file="./ermine_all_0vs14_0.05.csv",  row.names = FALSE)


# SIMPLIFYING PATHWAYS

maxP <- 200
clusters_0vs14 <- list()

# UP
obj_0vs14 <- EnrichListPDup_Res_0vs14 %>%
  dplyr::filter(CorrectedPvalue < 0.05) %>%
  arrange(CorrectedPvalue) %>%
  dplyr::select(Name, GeneMembers, CorrectedPvalue, NumGenes)

pathways_0vs14 <- apply(obj_0vs14, 1, function(x) unlist(strsplit(x[[2]], '\\|')))
names(pathways_0vs14) <- obj_0vs14$Name
pathways_0vs14 <- qdapTools::list2df(pathways_0vs14, col1="gene", col2="Name")[,c(2,1)]
clustersUP_0vs14 <- cluster_pathways(pathways_0vs14, method="overlap", subsetsize=maxP)

names(clustersUP_0vs14) <- c("Name", "members")
clustersUP_0vs14 <- left_join(clustersUP_0vs14, obj_0vs14 %>% dplyr::select(-GeneMembers, members=Name), by="members") %>%
  mutate(log10P=-log10(CorrectedPvalue)) %>%
  arrange(CorrectedPvalue)
clusters_0vs14[["UP"]] <- clustersUP_0vs14

# DOWN
obj_D_0vs14 <- EnrichListPDdown_Res_0vs14 %>%
  dplyr::filter(CorrectedPvalue < 0.05) %>%
  arrange(CorrectedPvalue) %>%
  dplyr::select(Name, GeneMembers, CorrectedPvalue, NumGenes)
pathways_D_0vs14 <- apply(obj_D_0vs14, 1, function(x) unlist(strsplit(x[[2]], '\\|')))
names(pathways_D_0vs14) <- obj_D_0vs14$Name
pathways_D_0vs14 <- qdapTools::list2df(pathways_D_0vs14, col1="gene", col2="Name")[,c(2,1)]

clustersDOWN_0vs14 <- cluster_pathways(pathways_D_0vs14, method="overlap", subsetsize=maxP)

names(clustersDOWN_0vs14) <- c("Name", "members")
clustersDOWN_0vs14 <- left_join(clustersDOWN_0vs14, obj_D_0vs14 %>% dplyr::select(-GeneMembers, members=Name), by="members") %>%
  mutate(log10P=-log10(CorrectedPvalue)) %>%
  arrange(CorrectedPvalue)
clusters_0vs14[["DOWN"]] <- clustersDOWN_0vs14

## Drwa Treemap: Figure 6 a
maxD <- 100
svglite("./Treemap_0vs14_F.svg", width = 8, height = 3.5, pointsize = 10)
draw_treemap_F(dtf=bind_rows(clusters_0vs14[["UP"]],
                             clusters_0vs14[["DOWN"]] %>% 
                             mutate(log10P=-log10P)) %>%
                             arrange(desc(abs(log10P))) %>%
                             head(n=maxD),
                             title = "a Pathways enriched at Braak stage 1-4",
                             vp=NULL)
dev.off()


# Combine UP and DOWN including the simplified pathway name. Save file and count number of simplified pathways
clustersDOWN_0vs14$Direction <- rep("DOWN", length(clustersDOWN_0vs14$Name))
clustersUP_0vs14$Direction <- rep("UP", length(clustersUP_0vs14$Name))

res_0vs14_0.05_simp <- bind_rows(clustersDOWN_0vs14, clustersUP_0vs14)

# Save
write.csv(as.data.frame(res_0vs14_0.05_simp), file="./ermine_all_0vs14_0.05_simp.csv",  row.names = FALSE)

# Count number of simplified pathways
res_0vs14_0.05_simp$Name <- as.factor(res_0vs14_0.05_simp$Name)
nlevels(res_0vs14_0.05_simp$Name) 


##### 0 vs 5
#################

# Prepare dataset
str(all_genes_0vs5)
all_genes_0vs5$log2FoldChange <- as.numeric(all_genes_0vs5$log2FoldChange)
all_genes_0vs5$pvalue <- as.numeric(all_genes_0vs5$pvalue)

# Transform p-values
all_genes_0vs5$DownPval <- apply(all_genes_0vs5 %>% dplyr::select(log2FoldChange, pvalue), 1, function(x){
  if(x[1] < 0){
    x[2]/2
  } else {
    1-x[2]/2
  }
})
all_genes_0vs5$DownPvalAdj <- p.adjust(all_genes_0vs5$DownPval, "BH")

all_genes_0vs5$UpPval <- apply(all_genes_0vs5 %>% dplyr::select(log2FoldChange, pvalue), 1, function(x){
  if(x[1] > 0){
    x[2]/2
  } else {
    1-x[2]/2
  }
})
all_genes_0vs5$UpPvalAdj <- p.adjust(all_genes_0vs5$UpPval, "BH")


### Find pathways for genes down-regulated from Braak Lewy body stage 0 to Braak Lewy body stages 5

# scores = A data.frame. Rownames have to be gene identifiers (eg. probes, must be unique), followed by any number of columns. 
all_genes_0vs51 <- all_genes_0vs5[!(is.na(all_genes_0vs5$hgnc_symbol) | all_genes_0vs5$hgnc_symbol==""), ] 
all_genes_0vs51 <- all_genes_0vs51[!duplicated(all_genes_0vs51$hgnc_symbol), ] 
rownames(all_genes_0vs51) <- all_genes_0vs51$hgnc_symbol


EnrichListPDdown_0vs5 <- gsr(scores=all_genes_0vs51, scoreColumn="DownPval",
      bigIsBetter=FALSE, logTrans=TRUE, annotation=GenericHumanAnno, 
      aspects=c("Biological Process"),
      iterations=200000)

EnrichListPDdown_Res_0vs5 <- EnrichListPDdown_0vs5$results %>%
    mutate(Direction="DOWN") %>%
    arrange(Pval)


EnrichListPDup_0vs5 <-gsr(scores=all_genes_0vs51, scoreColumn="UpPval",
                     bigIsBetter=FALSE, logTrans=TRUE, annotation=GenericHumanAnno, 
                     aspects=c("Biological Process"),
                     iterations=200000)

EnrichListPDup_Res_0vs5 <- EnrichListPDup_0vs5$results %>%
  mutate(Direction="UP") %>%
  arrange(Pval)

res_0vs5 <- bind_rows(EnrichListPDdown_Res_0vs5, EnrichListPDup_Res_0vs5) %>%
              arrange(CorrectedPvalue)

# Save
write.csv(as.data.frame(res_0vs5), file="./ermine_all_0vs5.csv",  row.names = FALSE)

### Subset only significant pathways
res_0vs5_0.05 <- res_0vs5 %>%
  dplyr::filter(CorrectedPvalue < 0.05) %>%
  arrange(CorrectedPvalue) 

### Count pathways
sum(res_0vs5_0.05$Direction == "UP", na.rm = TRUE) 
sum(res_0vs5_0.05$Direction == "DOWN", na.rm = TRUE) 

# Save
write.csv(as.data.frame(res_0vs5_0.05), file="./ermine_all_0vs5_0.05.csv",  row.names = FALSE)


# SIMPLIFYING PATHWAYS

clusters_0vs5 <- list()

# UP
obj_0vs5 <- EnrichListPDup_Res_0vs5 %>%
  dplyr::filter(CorrectedPvalue < 0.05) %>%
  arrange(CorrectedPvalue) %>%
  dplyr::select(Name, GeneMembers, CorrectedPvalue, NumGenes)

pathways_0vs5 <- apply(obj_0vs5, 1, function(x) unlist(strsplit(x[[2]], '\\|')))
names(pathways_0vs5) <- obj_0vs5$Name
pathways_0vs5 <- qdapTools::list2df(pathways_0vs5, col1="gene", col2="Name")[,c(2,1)]
clustersUP_0vs5 <- cluster_pathways(pathways_0vs5, method="overlap", subsetsize=maxP)

names(clustersUP_0vs5) <- c("Name", "members")
clustersUP_0vs5 <- left_join(clustersUP_0vs5, obj_0vs5 %>% dplyr::select(-GeneMembers, members=Name), by="members") %>%
  mutate(log10P=-log10(CorrectedPvalue)) %>%
  arrange(CorrectedPvalue)
clusters_0vs5[["UP"]] <- clustersUP_0vs5

# DOWN
obj_D_0vs5 <- EnrichListPDdown_Res_0vs5 %>%
  dplyr::filter(CorrectedPvalue < 0.05) %>%
  arrange(CorrectedPvalue) %>%
  dplyr::select(Name, GeneMembers, CorrectedPvalue, NumGenes)
pathways_D_0vs5 <- apply(obj_D_0vs5, 1, function(x) unlist(strsplit(x[[2]], '\\|')))
names(pathways_D_0vs5) <- obj_D_0vs5$Name
pathways_D_0vs5 <- qdapTools::list2df(pathways_D_0vs5, col1="gene", col2="Name")[,c(2,1)]

clustersDOWN_0vs5 <- cluster_pathways(pathways_D_0vs5, method="overlap", subsetsize=maxP)

names(clustersDOWN_0vs5) <- c("Name", "members")
clustersDOWN_0vs5 <- left_join(clustersDOWN_0vs5, obj_D_0vs5 %>% dplyr::select(-GeneMembers, members=Name), by="members") %>%
  mutate(log10P=-log10(CorrectedPvalue)) %>%
  arrange(CorrectedPvalue)
clusters_0vs5[["DOWN"]] <- clustersDOWN_0vs5

#Draw Treemap: Figure 6 b
maxD <- 100
svglite("./Treemap_0vs5_F.svg", width = 8, height = 3.5, pointsize = 10)
draw_treemap_F(dtf=bind_rows(clusters_0vs5[["UP"]],
                             clusters_0vs5[["DOWN"]] %>% 
                             mutate(log10P=-log10P)) %>%
                             arrange(desc(abs(log10P))) %>%
                             head(n=maxD),
                             title= "b Pathways enriched at Braak stage 5",
                             vp=NULL)
dev.off()


# Combine UP and DOWN including the simplified pathway name. Save file and count number of simplified pathways
clustersDOWN_0vs5$Direction <- rep("DOWN", length(clustersDOWN_0vs5$Name))
clustersUP_0vs5$Direction <- rep("UP", length(clustersUP_0vs5$Name))

res_0vs5_0.05_simp <- bind_rows(clustersDOWN_0vs5, clustersUP_0vs5)

# Save
write.csv(as.data.frame(res_0vs5_0.05_simp), file="./ermine_all_0vs5_0.05_simp.csv",  row.names = FALSE)

# Count simplified pathways
res_0vs5_0.05_simp$Name <- as.factor(res_0vs5_0.05_simp$Name)
nlevels(res_0vs5_0.05_simp$Name) 


##### 0 vs 6
#################

# Prepare dataset
str(all_genes_0vs6)
all_genes_0vs6$log2FoldChange <- as.numeric(all_genes_0vs6$log2FoldChange)
all_genes_0vs6$pvalue <- as.numeric(all_genes_0vs6$pvalue)

# Transfrom p-values
all_genes_0vs6$DownPval <- apply(all_genes_0vs6 %>% dplyr::select(log2FoldChange, pvalue), 1, function(x){
  if(x[1] < 0){
    x[2]/2
  } else {
    1-x[2]/2
  }
})
all_genes_0vs6$DownPvalAdj <- p.adjust(all_genes_0vs6$DownPval, "BH")

all_genes_0vs6$UpPval <- apply(all_genes_0vs6 %>% dplyr::select(log2FoldChange, pvalue), 1, function(x){
  if(x[1] > 0){
    x[2]/2
  } else {
    1-x[2]/2
  }
})
all_genes_0vs6$UpPvalAdj <- p.adjust(all_genes_0vs6$UpPval, "BH")

### Find pathways for genes down-regulated from Braak Lewy body stage 0 to Braak Lewy body stages 6

#scores = A data.frame. Rownames have to be gene identifiers (eg. probes, must be unique), followed by any number of columns. 
all_genes_0vs61 <- all_genes_0vs6[!(is.na(all_genes_0vs6$hgnc_symbol) | all_genes_0vs6$hgnc_symbol==""), ]
all_genes_0vs61 <- all_genes_0vs61[!duplicated(all_genes_0vs61$hgnc_symbol), ]
rownames(all_genes_0vs61) <- all_genes_0vs61$hgnc_symbol


EnrichListPDdown_0vs6 <- gsr(scores=all_genes_0vs61, scoreColumn="DownPval",
                        bigIsBetter=FALSE, logTrans=TRUE, annotation=GenericHumanAnno, 
                        aspects=c("Biological Process"),
                        iterations=200000)

EnrichListPDdown_Res_0vs6 <- EnrichListPDdown_0vs6$results %>%
  mutate(Direction="DOWN") %>%
  arrange(Pval)


EnrichListPDup_0vs6 <-gsr(scores=all_genes_0vs61, scoreColumn="UpPval",
                     bigIsBetter=FALSE, logTrans=TRUE, annotation=GenericHumanAnno, 
                     aspects=c("Biological Process"),
                     iterations=200000)

EnrichListPDup_Res_0vs6 <- EnrichListPDup_0vs6$results %>%
  mutate(Direction="UP") %>%
  arrange(Pval)

res_0vs6 <- bind_rows(EnrichListPDdown_Res_0vs6, EnrichListPDup_Res_0vs6) %>%
  arrange(CorrectedPvalue)

# Save
write.csv(as.data.frame(res_0vs6), file="./ermine_all_0vs6.csv",  row.names = FALSE)

### Subset only significant pathways
res_0vs6_0.05 <- res_0vs6 %>%
  dplyr::filter(CorrectedPvalue < 0.05) %>%
  arrange(CorrectedPvalue) 

### Count significant pathways
sum(res_0vs6_0.05$Direction == "UP", na.rm = TRUE)
sum(res_0vs6_0.05$Direction == "DOWN", na.rm = TRUE) 

# Save
write.csv(as.data.frame(res_0vs6_0.05), file="./ermine_all_0vs6_0.05.csv",  row.names = FALSE)


# SIMPLIFYING PATHWAYS

clusters_0vs6 <- list()

# UP
obj_0vs6 <- EnrichListPDup_Res_0vs6 %>%
  dplyr::filter(CorrectedPvalue < 0.05) %>%
  arrange(CorrectedPvalue) %>%
  dplyr::select(Name, GeneMembers, CorrectedPvalue, NumGenes)

pathways_0vs6 <- apply(obj_0vs6, 1, function(x) unlist(strsplit(x[[2]], '\\|')))
names(pathways_0vs6) <- obj_0vs6$Name
pathways_0vs6 <- qdapTools::list2df(pathways_0vs6, col1="gene", col2="Name")[,c(2,1)]
clustersUP_0vs6 <- cluster_pathways(pathways_0vs6, method="overlap", subsetsize=maxP)

names(clustersUP_0vs6) <- c("Name", "members")
clustersUP_0vs6 <- left_join(clustersUP_0vs6, obj_0vs6 %>% dplyr::select(-GeneMembers, members=Name), by="members") %>%
  mutate(log10P=-log10(CorrectedPvalue)) %>%
  arrange(CorrectedPvalue)
clusters_0vs6[["UP"]] <- clustersUP_0vs6

# DOWN
obj_D_0vs6 <- EnrichListPDdown_Res_0vs6 %>%
  dplyr::filter(CorrectedPvalue < 0.05) %>%
  arrange(CorrectedPvalue) %>%
  dplyr::select(Name, GeneMembers, CorrectedPvalue, NumGenes)
pathways_D_0vs6 <- apply(obj_D_0vs6, 1, function(x) unlist(strsplit(x[[2]], '\\|')))
names(pathways_D_0vs6) <- obj_D_0vs6$Name
pathways_D_0vs6 <- qdapTools::list2df(pathways_D_0vs6, col1="gene", col2="Name")[,c(2,1)]

clustersDOWN_0vs6 <- cluster_pathways(pathways_D_0vs6, method="overlap", subsetsize=maxP)

names(clustersDOWN_0vs6) <- c("Name", "members")
clustersDOWN_0vs6 <- left_join(clustersDOWN_0vs6, obj_D_0vs6 %>% dplyr::select(-GeneMembers, members=Name), by="members") %>%
  mutate(log10P=-log10(CorrectedPvalue)) %>%
  arrange(CorrectedPvalue)
clusters_0vs6[["DOWN"]] <- clustersDOWN_0vs6

# Draw Treemap: Figure 6 c
maxD <- 100
svglite("./Treemap_0vs6_F.svg", width = 8, height = 3.5, pointsize = 10)
draw_treemap_F(dtf=bind_rows(clusters_0vs6[["UP"]],
                             clusters_0vs14[["DOWN"]] %>% mutate(log10P=-log10P)) %>%
                 arrange(desc(abs(log10P))) %>%
                 head(n=maxD),
               title= "c Pathways enriched at Braak stage 6",
               vp=NULL)
dev.off()



# Combine UP and DOWN including the simplified pathway name. Save file and count number of simplified pathways
clustersDOWN_0vs6$Direction <- rep("DOWN", length(clustersDOWN_0vs6$Name))
clustersUP_0vs6$Direction <- rep("UP", length(clustersUP_0vs6$Name))

res_0vs6_0.05_simp <- bind_rows(clustersDOWN_0vs6, clustersUP_0vs6)

# Save 
write.csv(as.data.frame(res_0vs6_0.05_simp), file="./ermine_all_0vs6_0.05_simp.csv",  row.names = FALSE)

# Count simplified pathways
res_0vs6_0.05_simp$Name <- as.factor(res_0vs6_0.05_simp$Name)
nlevels(res_0vs6_0.05_simp$Name) #56






# Analysis of enrichment signal

########## 1-4 vs 5
#  up
sigPs_1 <- unique(c(res_0vs14 %>% dplyr::filter(CorrectedPvalue<0.05, Direction=="UP") %>% pull(Name),
                    res_0vs5 %>% dplyr::filter(CorrectedPvalue<0.05, Direction=="UP") %>% pull(Name)))

# In this list we have on top, what loses significance in the upregulated list
# and in the bottom, what becomes significantly upregulated
up_1 <- inner_join(res_0vs14 %>% dplyr::filter(Direction=="UP", Name %in% sigPs_1) %>% dplyr::select(Name, padj_14=CorrectedPvalue),
                   res_0vs5 %>%  dplyr::filter(Direction=="UP", Name %in% sigPs_1) %>% dplyr::select(Name, padj_5=CorrectedPvalue),
                   by=c("Name")) %>%
  mutate(Delta=-log10(padj_5)+log10(padj_14)) %>% 
  dplyr::filter(!(padj_14 < 0.05 & padj_5 < 0.05)) %>%
  arrange(Delta)


# down
sigPs_D_1 <- unique(c(res_0vs14 %>% dplyr::filter(CorrectedPvalue<0.05, Direction=="DOWN") %>% pull(Name),
                      res_0vs5 %>% dplyr::filter(CorrectedPvalue<0.05, Direction=="DOWN") %>% pull(Name)))

# In this list we have on top, what loses significance in the downregulated list
# and in the bottom, what becomes significantly downregulated
down_1 <- inner_join(res_0vs14 %>% dplyr::filter(Direction=="DOWN", Name %in% sigPs_D_1) %>% dplyr::select(Name, padj_14=CorrectedPvalue),
                     res_0vs5 %>% dplyr::filter(Direction=="DOWN", Name %in% sigPs_D_1) %>% dplyr::select(Name, padj_5=CorrectedPvalue),
                     by=c("Name")) %>%
  mutate(Delta=-log10(padj_5)+log10(padj_14)) %>% 
  dplyr::filter(!(padj_14 < 0.05 & padj_5 < 0.05)) %>%
  arrange(Delta)

# Save
write.csv(up_1, "./pvalue_delta_enrichment_UP_14_5.csv")
write.csv(down_1, "./pvalue_delta_enrichment_DOWN_14_5.csv")



########## 5 vs 6
#  up
sigPs_2 <- unique(c(res_0vs5 %>% dplyr::filter(CorrectedPvalue<0.05, Direction=="UP") %>% pull(Name),
                  res_0vs6 %>% dplyr::filter(CorrectedPvalue<0.05, Direction=="UP") %>% pull(Name)))

# In this list we have on top, what loses significance in the upregulated list
# and in the bottom, what becomes significantly upregulated
up_2 <- inner_join(res_0vs5 %>% dplyr::filter(Direction=="UP", Name %in% sigPs_2) %>% dplyr::select(Name, padj_5=CorrectedPvalue),
                    res_0vs6 %>%  dplyr::filter(Direction=="UP", Name %in% sigPs_2) %>% dplyr::select(Name, padj_6=CorrectedPvalue),
                    by=c("Name")) %>%
  mutate(Delta=-log10(padj_6)+log10(padj_5)) %>% 
  dplyr::filter(!(padj_5 < 0.05 & padj_6 < 0.05)) %>%
  arrange(Delta)

# down
sigPs_D_2 <- unique(c(res_0vs5 %>% dplyr::filter(CorrectedPvalue<0.05, Direction=="DOWN") %>% pull(Name),
                  res_0vs6 %>% dplyr::filter(CorrectedPvalue<0.05, Direction=="DOWN") %>% pull(Name)))

# In this list we have on top, what loses significance in the downregulated list
# and in the bottom, what becomes significantly downregulated
down_2 <- inner_join(res_0vs5 %>% dplyr::filter(Direction=="DOWN", Name %in% sigPs_D_2) %>% dplyr::select(Name, padj_5=CorrectedPvalue),
                      res_0vs6 %>% dplyr::filter(Direction=="DOWN", Name %in% sigPs_D_2) %>% dplyr::select(Name, padj_6=CorrectedPvalue),
                      by=c("Name")) %>%
  mutate(Delta=-log10(padj_6)+log10(padj_5)) %>% 
  dplyr::filter(!(padj_5 < 0.05 & padj_6 < 0.05)) %>%
  arrange(Delta)

# Save
write.csv(up_2, "./pvalue_delta_enrichment_UP_5_6.csv")
write.csv(down_2, "./pvalue_delta_enrichment_DOWN_5_6.csv")


### Simplify pathways gaining or losing significance from one Braak Lewy body stage to the following one

#### SIMPLIFY UP
sig_14_up <- res_0vs14 %>% dplyr::filter(CorrectedPvalue<0.05, Direction=="UP")
sig_5_up <-   res_0vs5 %>% dplyr::filter(CorrectedPvalue<0.05, Direction=="UP")
sig_6_up <-   res_0vs6 %>% dplyr::filter(CorrectedPvalue<0.05, Direction=="UP")

### 14 - 5 up
sig_14_5_up <- rbind(sig_14_up, sig_5_up)
sig_14_5_up <- unique(sig_14_5_up)
sig_14_5_up <- sig_14_5_up %>% select(Name, GeneMembers,CorrectedPvalue, NumGenes)
sig_14_5_up <- as_tibble(sig_14_5_up)

pathways_up_14_5 <- apply(sig_14_5_up, 1, function(x) unlist(strsplit(x[[2]], '\\|')))
names(pathways_up_14_5) <- sig_14_5_up$Name
pathways_up_14_5 <- qdapTools::list2df(pathways_up_14_5, col1="gene", col2="Name")[,c(2,1)]
clustersUP_up_14_5 <- cluster_pathways(pathways_up_14_5, method="overlap", subsetsize=maxP)

simpUP_up_14_5 <- merge(up_1, clustersUP_up_14_5, by.x = "Name", by.y = "member")
simpUP_up_14_5 <- simpUP_up_14_5 %>% arrange(Delta)

# Save
write.csv(simpUP_up_14_5 , "./pvalue_delta_enrichment_UP_14_5_simp.csv")


### 5- 6 up
sig_5_6_up <- rbind(sig_5_up, sig_6_up)
sig_5_6_up <- sig_5_6_up[!duplicated(sig_5_6_up[ , "Name"]),]
sig_5_6_up <- sig_5_6_up %>% select(Name, GeneMembers,CorrectedPvalue, NumGenes)
sig_5_6_up <- as_tibble(sig_5_6_up)

pathways_up_5_6 <- apply(sig_5_6_up, 1, function(x) unlist(strsplit(x[[2]], '\\|')))
names(pathways_up_5_6) <- sig_5_6_up$Name
pathways_up_5_6 <- qdapTools::list2df(pathways_up_5_6, col1="gene", col2="Name")[,c(2,1)]
clustersUP_up_5_6 <- cluster_pathways(pathways_up_5_6, method="overlap", subsetsize=maxP)

simpUP_up_5_6 <- merge(up_2, clustersUP_up_5_6, by.x = "Name", by.y = "member")
simpUP_up_5_6 <- simpUP_up_5_6 %>% arrange(Delta)

# Save
write.csv(simpUP_up_5_6, "./pvalue_delta_enrichment_UP_5_6_simp.csv")


#### SIMPLIFY DOWN
sig_14_down <- res_0vs14 %>% dplyr::filter(CorrectedPvalue<0.05, Direction=="DOWN")
sig_5_down <-   res_0vs5 %>% dplyr::filter(CorrectedPvalue<0.05, Direction=="DOWN")
sig_6_down <-   res_0vs6 %>% dplyr::filter(CorrectedPvalue<0.05, Direction=="DOWN")

### 14 - 5 down
sig_14_5_down <- rbind(sig_14_down, sig_5_down)
sig_14_5_down <- sig_14_5_down[!duplicated(sig_14_5_down[ , "Name"]),]
sig_14_5_down <- sig_14_5_down %>% select(Name, GeneMembers,CorrectedPvalue, NumGenes)
sig_14_5_down <- as_tibble(sig_14_5_down)

pathways_down_14_5 <- apply(sig_14_5_down, 1, function(x) unlist(strsplit(x[[2]], '\\|')))
names(pathways_down_14_5) <- sig_14_5_down$Name
pathways_down_14_5 <- qdapTools::list2df(pathways_down_14_5, col1="gene", col2="Name")[,c(2,1)]
clustersDOWN_down_14_5 <- cluster_pathways(pathways_down_14_5, method="overlap", subsetsize=maxP)

simpDOWN_14_5 <- merge(down_1, clustersDOWN_down_14_5, by.x = "Name", by.y = "member")
simpDOWN_14_5 <- simpDOWN_14_5 %>% arrange(Delta)

# Save
write.csv(simpDOWN_14_5, "./pvalue_delta_enrichment_DOWN_14_5_simp.csv")

### 5 - 6 down
sig_5_6_down <- rbind(sig_5_down, sig_6_down)
sig_5_6_down <- sig_5_6_down[!duplicated(sig_5_6_down[ , "Name"]),]
sig_5_6_down <- sig_5_6_down %>% select(Name, GeneMembers,CorrectedPvalue, NumGenes)
sig_5_6_down <- as_tibble(sig_5_6_down)

pathways_down_5_6 <- apply(sig_5_6_down, 1, function(x) unlist(strsplit(x[[2]], '\\|')))
names(pathways_down_5_6) <- sig_5_6_down$Name
pathways_down_5_6 <- qdapTools::list2df(pathways_down_5_6, col1="gene", col2="Name")[,c(2,1)]
clustersDOWN_down_5_6 <- cluster_pathways(pathways_down_5_6, method="overlap", subsetsize=maxP)

simpDOWN_5_6 <- merge(down_2, clustersDOWN_down_5_6, by.x = "Name", by.y = "member")
simpDOWN_5_6 <- simpDOWN_5_6 %>% arrange(Delta)

# Save
write.csv(simpDOWN_5_6, "./pvalue_delta_enrichment_DOWN_5_6_simp.csv")

### Combine up and down regulated

### 1-4 to 5
## add column with direction
simpUP_up_14_5$Direction <- rep("UP", n=length(simpUP_up_14_5))
simpDOWN_14_5$Direction <- rep("DOWN", n=length(simpDOWN_14_5))
simpUP_up_5_6$Direction <- rep("UP", n=length(simpUP_up_5_6))
simpDOWN_5_6$Direction <- rep("DOWN", n=length(simpDOWN_5_6))

simp_all_14_5 <- rbind(simpUP_up_14_5,simpDOWN_14_5) %>%
  arrange(Delta)
simp_all_14_5 <- simp_all_14_5[, c("title", "Name", "padj_14", "padj_5", "Delta", "Direction")]

sum(simp_all_14_5$Direction == "UP" & simp_all_14_5$Delta > 0, na.rm = TRUE) # pathway up regulated that gained significance
sum(simp_all_14_5$Direction == "DOWN" & simp_all_14_5$Delta > 0, na.rm = TRUE) # pathway down regulated that gained significance
sum(simp_all_14_5$Direction == "UP" & simp_all_14_5$Delta < 0, na.rm = TRUE) # pathway up regulated that lost significance
sum(simp_all_14_5$Direction == "DOWN" & simp_all_14_5$Delta < 0, na.rm = TRUE) #18 pathway down regulated that lost significance

# Save
write.csv(simp_all_14_5, "./pvalue_delta_enrichment_all_14_5_simp.csv")

### 5 to 6
simp_all_5_6 <- rbind(simpUP_up_5_6,simpDOWN_5_6) %>%
  arrange(Delta)
simp_all_5_6 <- simp_all_5_6[, c("title", "Name", "padj_5", "padj_6", "Delta", "Direction")]

sum(simp_all_5_6$Direction == "UP" & simp_all_5_6$Delta > 0, na.rm = TRUE) # pathway up regulated that gained significance
sum(simp_all_5_6$Direction == "DOWN" & simp_all_5_6$Delta > 0, na.rm = TRUE) # pathway down regulated that gained significance
sum(simp_all_5_6$Direction == "UP" & simp_all_5_6$Delta < 0, na.rm = TRUE) # pathway up regulated that lost significance
sum(simp_all_5_6$Direction == "DOWN" & simp_all_5_6$Delta < 0, na.rm = TRUE) # pathway down regulated that lost significance

# Save
write.csv(simp_all_5_6, "./pvalue_delta_enrichment_all_5_6_simp.csv")


#### Concordant pathways

###Upregulated

## 1-4 and 5
cs <- c("0vs14", "0vs5")
concordantUP_14_5 <- inner_join(
  EnrichListPDup_Res_0vs14 %>%
    arrange(CorrectedPvalue) %>%
    distinct(ID, .keep_all=TRUE) %>%
    select(Name, GeneMembers, ID, RawScore, Pval, CorrectedPvalue, NumGenes) %>%
    filter(CorrectedPvalue < 0.05),
  EnrichListPDup_Res_0vs5 %>%
    arrange(CorrectedPvalue) %>%
    distinct(ID, .keep_all=TRUE) %>%
    select(Name, GeneMembers, ID, RawScore, Pval, CorrectedPvalue, NumGenes) %>%
    filter(CorrectedPvalue < 0.05),
  by=c("Name", "ID"), suffix=paste0("_", cs)) %>%
  mutate(GeneMembers=GeneMembers_0vs14) %>%
  select(Name, GeneMembers, CorrectedPvalue_0vs14, CorrectedPvalue_0vs5, NumGenes_0vs14, NumGenes_0vs5)
# UP: no concordant pathway between 1-4 and 5

## 5 and 6
cs1 <- c("0vs5", "0vs6")
concordantUP_5_6 <- inner_join(
  EnrichListPDup_Res_0vs5 %>%
    arrange(CorrectedPvalue) %>%
    distinct(ID, .keep_all=TRUE) %>%
    select(Name, GeneMembers, ID, RawScore, Pval, CorrectedPvalue, NumGenes) %>%
    filter(CorrectedPvalue < 0.05),
  EnrichListPDup_Res_0vs6 %>%
    arrange(CorrectedPvalue) %>%
    distinct(ID, .keep_all=TRUE) %>%
    select(Name, GeneMembers, ID, RawScore, Pval, CorrectedPvalue, NumGenes) %>%
    filter(CorrectedPvalue < 0.05),
  by=c("Name", "ID"), suffix=paste0("_", cs1)) %>%
  mutate(GeneMembers=GeneMembers_0vs5) %>%
  select(Name, GeneMembers, CorrectedPvalue_0vs5, CorrectedPvalue_0vs6, NumGenes_0vs5, NumGenes_0vs6)

# Save
write.csv(concordantUP_5_6, "./concordantUP_5_6.csv")

## 1-4 and 6
cs2 <- c("0vs14", "0vs6")
concordantUP_14_6 <- inner_join(
  EnrichListPDup_Res_0vs14 %>%
    arrange(CorrectedPvalue) %>%
    distinct(ID, .keep_all=TRUE) %>%
    select(Name, GeneMembers, ID, RawScore, Pval, CorrectedPvalue, NumGenes) %>%
    filter(CorrectedPvalue < 0.05),
  EnrichListPDup_Res_0vs6 %>%
    arrange(CorrectedPvalue) %>%
    distinct(ID, .keep_all=TRUE) %>%
    select(Name, GeneMembers, ID, RawScore, Pval, CorrectedPvalue, NumGenes) %>%
    filter(CorrectedPvalue < 0.05),
  by=c("Name", "ID"), suffix=paste0("_", cs2)) %>%
  mutate(GeneMembers=GeneMembers_0vs14) %>%
  select(Name, GeneMembers, CorrectedPvalue_0vs14, CorrectedPvalue_0vs6, NumGenes_0vs14, NumGenes_0vs6)
# UP: 0 concordant pathway between 1-4 and 6


#### Downregulated

## 1-4 and 5
concordantDOWN_14_5 <- inner_join(
  EnrichListPDdown_Res_0vs14 %>%
    arrange(CorrectedPvalue) %>%
    distinct(ID, .keep_all=TRUE) %>%
    select(Name, GeneMembers, ID, RawScore, Pval, CorrectedPvalue, NumGenes) %>%
    filter(CorrectedPvalue < 0.05),
  EnrichListPDdown_Res_0vs5 %>%
    arrange(CorrectedPvalue) %>%
    distinct(ID, .keep_all=TRUE) %>%
    select(Name, GeneMembers, ID, RawScore, Pval, CorrectedPvalue, NumGenes) %>%
    filter(CorrectedPvalue < 0.05),
  by=c("Name", "ID"), suffix=paste0("_", cs)) %>%
  mutate(GeneMembers=GeneMembers_0vs14) %>%
  select(Name, GeneMembers, CorrectedPvalue_0vs14, CorrectedPvalue_0vs5, NumGenes_0vs14, NumGenes_0vs5)

# Simplify pathways
pathways_D_14_5 <- apply(concordantDOWN_14_5, 1, function(x) unlist(strsplit(x[[2]], '\\|')))
names(pathways_D_14_5) <- concordantDOWN_14_5$Name
pathways_D_14_5 <- qdapTools::list2df(pathways_D_14_5, col1="gene", col2="Name")[,c(2,1)]

clustersDOWN_14_5 <- cluster_pathways(pathways_D_14_5, method="overlap", subsetsize=maxP, threshold=0.5)

names(clustersDOWN_14_5) <- c("Name", "members")
clustersDOWN_14_5 <- left_join(clustersDOWN_14_5, concordantDOWN_14_5 %>% select(-GeneMembers, members=Name), by="members") %>%
  mutate(log10P_0vs14=-log10(CorrectedPvalue_0vs14), log10P_0vs5=-log10(CorrectedPvalue_0vs5),
         log10P=sqrt(log10P_0vs14*log10P_0vs5), NumGenes=(NumGenes_0vs14+NumGenes_0vs5)/2) %>%
  arrange(CorrectedPvalue_0vs14)

# Save
write.csv(clustersDOWN_14_5, "./concordantDOWN_14_5.csv")

## 5 and 6
concordantDOWN_5_6 <- inner_join(
  EnrichListPDdown_Res_0vs5 %>%
    arrange(CorrectedPvalue) %>%
    distinct(ID, .keep_all=TRUE) %>%
    select(Name, GeneMembers, ID, RawScore, Pval, CorrectedPvalue, NumGenes) %>%
    filter(CorrectedPvalue < 0.05),
  EnrichListPDdown_Res_0vs6 %>%
    arrange(CorrectedPvalue) %>%
    distinct(ID, .keep_all=TRUE) %>%
    select(Name, GeneMembers, ID, RawScore, Pval, CorrectedPvalue, NumGenes) %>%
    filter(CorrectedPvalue < 0.05),
  by=c("Name", "ID"), suffix=paste0("_", cs1)) %>%
  mutate(GeneMembers=GeneMembers_0vs5) %>%
  select(Name, GeneMembers, CorrectedPvalue_0vs5, CorrectedPvalue_0vs6, NumGenes_0vs5, NumGenes_0vs6)

# Simplify pathways
pathways_D_5_6 <- apply(concordantDOWN_5_6, 1, function(x) unlist(strsplit(x[[2]], '\\|')))
names(pathways_D_5_6) <- concordantDOWN_5_6$Name
pathways_D_5_6 <- qdapTools::list2df(pathways_D_5_6, col1="gene", col2="Name")[,c(2,1)]

clustersDOWN_5_6 <- cluster_pathways(pathways_D_5_6, method="overlap", subsetsize=maxP, threshold=0.5)

names(clustersDOWN_5_6) <- c("Name", "members")
clustersDOWN_5_6 <- left_join(clustersDOWN_5_6, concordantDOWN_5_6 %>% select(-GeneMembers, members=Name), by="members") %>%
  mutate(log10P_0vs5=-log10(CorrectedPvalue_0vs5), log10P_0vs6=-log10(CorrectedPvalue_0vs6),
         log10P=sqrt(log10P_0vs5*log10P_0vs6), NumGenes=(NumGenes_0vs5+NumGenes_0vs6)/2) %>%
  arrange(CorrectedPvalue_0vs5)

# Save
write.csv(clustersDOWN_5_6, "./concordantDOWN_5_6.csv")

## 1-4 and 6
concordantDOWN_14_6 <- inner_join(
  EnrichListPDdown_Res_0vs14 %>%
    arrange(CorrectedPvalue) %>%
    distinct(ID, .keep_all=TRUE) %>%
    select(Name, GeneMembers, ID, RawScore, Pval, CorrectedPvalue, NumGenes) %>%
    filter(CorrectedPvalue < 0.05),
  EnrichListPDdown_Res_0vs6 %>%
    arrange(CorrectedPvalue) %>%
    distinct(ID, .keep_all=TRUE) %>%
    select(Name, GeneMembers, ID, RawScore, Pval, CorrectedPvalue, NumGenes) %>%
    filter(CorrectedPvalue < 0.05),
  by=c("Name", "ID"), suffix=paste0("_", cs2)) %>%
  mutate(GeneMembers=GeneMembers_0vs14) %>%
  select(Name, GeneMembers, CorrectedPvalue_0vs14, CorrectedPvalue_0vs6, NumGenes_0vs14, NumGenes_0vs6)

# Simplify pathways
pathways_D_14_6 <- apply(concordantDOWN_14_6, 1, function(x) unlist(strsplit(x[[2]], '\\|')))
names(pathways_D_14_6) <- concordantDOWN_14_6$Name
pathways_D_14_6 <- qdapTools::list2df(pathways_D_14_6, col1="gene", col2="Name")[,c(2,1)]

clustersDOWN_14_6 <- cluster_pathways(pathways_D_14_6, method="overlap", subsetsize=maxP, threshold=0.5)

names(clustersDOWN_14_6) <- c("Name", "members")
clustersDOWN_14_6 <- left_join(clustersDOWN_14_6, concordantDOWN_14_6 %>% select(-GeneMembers, members=Name), by="members") %>%
  mutate(log10P_0vs14=-log10(CorrectedPvalue_0vs14), log10P_0vs6=-log10(CorrectedPvalue_0vs6),
         log10P=sqrt(log10P_0vs14*log10P_0vs6), NumGenes=(NumGenes_0vs14+NumGenes_0vs6)/2) %>%
  arrange(CorrectedPvalue_0vs14)

# Save
write.csv(clustersDOWN_14_6, "./concordantDOWN_14_6.csv")

####Find pathways common to all 3 stages

cs3 <- c("1", "2")
concordantDOWN_common <- inner_join(
  concordantDOWN_14_5,
  concordantDOWN_5_6,
  by=c("Name"), suffix=paste0("_", cs3))

concordantDOWN_common$GeneMembers_2 <- NULL
concordantDOWN_common$CorrectedPvalue_0vs5_2 <- NULL
concordantDOWN_common$NumGenes_0vs5_2 <- NULL
concordantDOWN_common <- concordantDOWN_common[, c("Name","GeneMembers_1","CorrectedPvalue_0vs14","CorrectedPvalue_0vs5_1", "CorrectedPvalue_0vs6",
                                                   "NumGenes_0vs14" , "NumGenes_0vs5_1"  ,  "NumGenes_0vs6")]

# Save
write.csv(concordantDOWN_common, "./concordantDOWN_common.csv")

# Simplify pathways
pathways_D_common <- apply(concordantDOWN_common, 1, function(x) unlist(strsplit(x[[2]], '\\|')))
names(pathways_D_common) <- concordantDOWN_common$Name
pathways_D_common <- qdapTools::list2df(pathways_D_common, col1="gene", col2="Name")[,c(2,1)]

clustersDOWN_common <- cluster_pathways(pathways_D_common, method="overlap", subsetsize=maxP, threshold=0.5)

names(clustersDOWN_common) <- c("Name", "members")
clustersDOWN_common <- left_join(clustersDOWN_common, concordantDOWN_common %>% select(-GeneMembers_1, members=Name), by="members") 

# Save
write.csv(clustersDOWN_common, "./concordantDOWN_common_simp.csv")



