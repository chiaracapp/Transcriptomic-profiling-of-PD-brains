### Demographics

### Packages required:
library("DESeq2")
library("dplyr")
library("ggplot2")
library("rstatix")

### Demographics

### Load dds
load("./dds_84.Rda")

coldata <- colData(dds)
coldata <- as.data.frame(coldata)

dem <- read.csv2("./ClinicopathVUmc_forChiara.csv", header = TRUE, sep = ",")
ID <- read.csv2("./Sample_ID.csv", header = TRUE, sep = ",")

ID <- ID[ID$X %in% rownames(coldata),]
dem <- dem[dem$Case_ID %in% ID$Sample_Name,]

dem <- merge(dem, ID, by.x = "Case_ID", by.y = "Sample_Name")
dem[,c(1:5,7,11,12,13,16)] <- NULL

coldata <- merge(coldata, dem, by.x = "row.names", by.y="X")
coldata$pH <- as.numeric(coldata$pH)

### Save file to make Supplementary Table 1
write.csv(as.data.frame(coldata), file="/Supplementary_Table_1.csv",  row.names = FALSE)


### Chi square test for categorical individuals’ demographics
## SEX
sex <- table(coldata$Braak_aSyn_groups_Plot, coldata$Sex)
sex_res <- chisq.test(sex)

## NFT stage
NFT <- table(coldata$Braak_aSyn_groups_Plot, coldata$Braak_NFT_stage)
chisq.test(NFT)

## Brain surgery
brain <- table(coldata$Braak_aSyn_groups_Plot, coldata$Brain_surgery)
chisq.test(brain)

## CERAD stage
cerad <- table(coldata$Braak_aSyn_groups_Plot, coldata$CERAD_stage)
chisq.test(cerad)

## AD level
ad <- table(coldata$Braak_aSyn_groups_Plot, coldata$AD.level)
chisq.test(ad)

## CAA type
caa <- table(coldata$Braak_aSyn_groups_Plot, coldata$CAA_type)
chisq.test(caa)


### Kruskal-Wallis Test for continuous individuals’ demographics
## RIN 
stat.test_RIN <- coldata %>%
  kruskal_test(RIN ~ Braak_aSyn_groups_Plot) %>%
  add_significance()

## Age at death 
stat.test_Age <- coldata %>%
  kruskal_test(Age_death ~ Braak_aSyn_groups_Plot) %>%
  add_significance()

## PMD
stat.test_PMD <- coldata %>%
  kruskal_test(PMD_min ~ Braak_aSyn_groups_Plot) %>%
  add_significance()

## pH
stat.test_pH <- coldata %>%
  kruskal_test(pH ~ Braak_aSyn_groups_Plot) %>%
  add_significance()


## dementia duration
coldata <- coldata %>% mutate(Dementia_duration_num =
                                case_when(Dementia_duration == "no dementia" ~ 0, 
                                          Dementia_duration == "1" ~  1, 
                                          Dementia_duration == "2" ~ 2,
                                          Dementia_duration == "3" ~ 3,
                                          Dementia_duration == "5" ~ 5,
                                          Dementia_duration == "6" ~ 6,
                                          Dementia_duration == "7" ~ 7,
                                          Dementia_duration == "8" ~ 8))

stat.test_demd <- coldata %>%
  kruskal_test(Dementia_duration_num ~ Braak_aSyn_groups_Plot) %>%
  add_significance()

## disease duraiton
coldata_case <- dplyr::filter(coldata, Braak_aSyn_groups_Plot != "0")

stat.test_dd <- coldata_case %>%
  kruskal_test(Disease_duration ~ Braak_aSyn_groups_Plot) %>%
  add_significance()


### Plot demographics and make Figure 2

## Sex
plot_sex <- ggplot(coldata, aes(Braak_aSyn_groups_Plot, fill = Sex)) +
  geom_bar(stat="count", width = 0.7, color = "black") +
  theme_minimal() +
  labs(title="Sex", x = NULL, y = "Count", fill = "Sex") +
  scale_fill_grey(start = 0.9, end = 0) +
  theme(legend.title = element_blank(),
        text = element_text(size = 16),
        plot.title = element_text(size=20, hjust = 0.5))

## NFT
plot_NFT <- ggplot(coldata, aes(Braak_aSyn_groups_Plot, fill = Braak_NFT_stage)) +
  geom_bar(stat="count", width = 0.7, color = "black") +
  theme_minimal() +
  labs(title="Braak NFT stage", x = NULL, y = "Count", fill = "Braak_NFT_stage") +
  scale_fill_grey(start = 0.9, end = 0) +
  theme(legend.title = element_blank(),
        text = element_text(size = 16),
        plot.title = element_text(size=20, hjust = 0.5))

## Brain surgery
plot_brain <- ggplot(coldata, aes(Braak_aSyn_groups_Plot, fill = Brain_surgery)) +
  geom_bar(stat="count", width = 0.7, color = "black") +
  theme_minimal() +
  labs(title="Brain surgery", x = NULL, y = "Count", fill = "Brain_surgery") +
  scale_fill_grey(start = 0.9, end = 0) +
  theme(legend.title = element_blank(),
        text = element_text(size = 16),
        plot.title = element_text(size=20, hjust = 0.5))

## CERAD stage
coldata$CERAD_stage <- factor(coldata$CERAD_stage, levels = c("frequent", "moderate", "sparse", "none"))

plot_cerad <- ggplot(coldata, aes(Braak_aSyn_groups_Plot, fill = CERAD_stage)) +
  geom_bar(stat="count", width = 0.7, color = "black") +
  theme_minimal() +
  labs(title="CERAD stage", x = NULL, y = "Count", fill = "CERAD_stage") +
  scale_fill_grey(start = 0.9, end = 0) +
  theme(legend.title = element_blank(),
        text = element_text(size = 16),
        plot.title = element_text(size=20, hjust = 0.5))

## Ad level
plot_ad <- ggplot(coldata, aes(Braak_aSyn_groups_Plot, fill = AD.level)) +
  geom_bar(stat="count", width = 0.7, color = "black") +
  theme_minimal() +
  labs(title="AD level", x = NULL, y = "Count", fill = "AD.level") +
  scale_fill_grey(start = 0.9, end = 0) +
  theme(legend.title = element_blank(),
        text = element_text(size = 16),
        plot.title = element_text(size=20, hjust = 0.5))


## CAA type
coldata$CAA_type <- factor(coldata$CAA_type, levels = c("2", "1", "no", "NA"))

plot_caa <- ggplot(coldata, aes(Braak_aSyn_groups_Plot, fill = CAA_type)) +
  geom_bar(stat="count", width = 0.7, color = "black") +
  theme_minimal() +
  labs(title="CAA type", x = NULL, y = "Count", fill = "CAA_type") +
  scale_fill_grey(start = 0.9, end = 0) +
  theme(legend.title = element_blank(),
        text = element_text(size = 16),
        plot.title = element_text(size=20, hjust = 0.5))


## Set colors
my_colors_ID <- c("0" = "#C6DBEF", "1-4" = "#6BAED6", "5" = "#08519C", "6" = "#08306B")

## RIN
plot_rin <- ggplot(coldata, aes(Braak_aSyn_groups_Plot, RIN, fill = Braak_aSyn_groups_Plot)) +
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, fill= "black") +
  theme_minimal() +
  scale_fill_manual(values = my_colors_ID) + 
  ylim(5,10.5) +
  labs(title="RIN", x = NULL, y = "RIN values") +
  theme(legend.position = "none",
        text = element_text(size = 16),
        plot.title = element_text(size=20, hjust = 0.5))

## Age at death 
plot_age <- ggplot(coldata, aes(Braak_aSyn_groups_Plot, Age_death, fill = Braak_aSyn_groups_Plot)) +
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, fill= "black") +
  theme_minimal() +
  scale_fill_manual(values = my_colors_ID) + 
  ylim(60,100) +
  labs(title="Age at death", x = NULL, y = "Years") +
  theme(legend.position = "none",
        text = element_text(size = 16),
        plot.title = element_text(size=20, hjust = 0.5))

## PMD
plot_PMD <- ggplot(coldata, aes(Braak_aSyn_groups_Plot, PMD_min, fill = Braak_aSyn_groups_Plot)) +
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, fill= "black") +
  theme_minimal() +
  scale_fill_manual(values = my_colors_ID) + 
  ylim(60,1300) +
  labs(title="PMD", x = NULL, y = "Minutes") +
  theme(legend.position = "none",
        text = element_text(size = 16),
        plot.title = element_text(size=20, hjust = 0.5))

## pH
plot_pH <- ggplot(coldata, aes(Braak_aSyn_groups_Plot, pH, fill = Braak_aSyn_groups_Plot)) +
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, fill= "black") +
  theme_minimal() +
  scale_fill_manual(values = my_colors_ID) + 
  ylim(5.8,7.8) +
  labs(title="pH", x = NULL, y = "pH") +
  theme(legend.position = "none",
        text = element_text(size = 16),
        plot.title = element_text(size=20, hjust = 0.5))

## dementia duration
plot_demd <- ggplot(coldata, aes(Braak_aSyn_groups_Plot, Dementia_duration_num, fill = Braak_aSyn_groups_Plot)) +
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, fill= "black") +
  theme_minimal() +
  scale_fill_manual(values = my_colors_ID) + 
  ylim(0,8) +
  labs(title="Dementia duration", x = NULL, y = "Years") +
  theme(legend.position = "none",
        text = element_text(size = 16),
        plot.title = element_text(size=20, hjust = 0.5))


## disease duraiton
my_colors_ID_1 <- c( "1-4" = "#6BAED6", "5" = "#08519C", "6" = "#08306B")

plot_dd <- ggplot(coldata_case, aes(Braak_aSyn_groups_Plot, Disease_duration, fill = Braak_aSyn_groups_Plot)) +
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, fill= "black") +
  theme_minimal() +
  scale_fill_manual(values = my_colors_ID_1) + 
  ylim(0,45) +
  labs(title="Disease duration", x = NULL, y = "Years") +
  theme(legend.position = "none",
        text = element_text(size = 16),
        plot.title = element_text(size=20, hjust = 0.5))


## Combine figures
continuous <- ggarrange(plot_age, plot_PMD, plot_demd, plot_dd, plot_rin, plot_pH,
                    ncol = 3,
                    nrow = 2)

categorical <- ggarrange(plot_sex, plot_NFT, plot_brain, plot_cerad, plot_ad, plot_caa,
                        ncol = 3,
                        nrow = 2)

figure <- ggarrange(continuous, categorical,
                         labels = c("a","b"),
                         ncol = 1,
                         nrow = 2)



ggsave("./Demographics.png", plot = figure, device = "png", width = 15, height = 15, units = 'in', dpi =300)
ggsave("./Demographics.svg", plot = figure, device = "svg", width = 10, height = 10, units = 'in', dpi =600)

