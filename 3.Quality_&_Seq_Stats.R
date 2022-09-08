###  Samples´ quality and sequencing stats

### Packages required:
library("DESeq2")
library("dplyr")
library("ggplot2")


### Load dds
load("./dds_84.Rda")

## oranize samples´ data
coldata <- colData(dds)
coldata <- as.data.frame(coldata)
coldata <- coldata %>% mutate(Braak_aSyn_groups_Plot =
                                                                   case_when(Braak_aSyn_groups == "0" ~ "0", 
                                                                             Braak_aSyn_groups == "1" ~ "1-4", 
                                                                             Braak_aSyn_groups == "2" ~ "5",
                                                                             Braak_aSyn_groups == "3" ~ "6" ))

## Import sequencing stats
seq_stats <- read.csv2("./Seq_stats.csv", header = TRUE, sep = ",")

## combine dataframes
coldata <- merge(coldata, seq_stats, by.x= "row.names", by.y = "Sample")
rownames(coldata) <- coldata$Row.names
coldata$Row.names <- NULL
col_num <- c("Salmon_percent_mapped","percent_gc","Hisat_mapped_passed_pct")
coldata[col_num] <- lapply(coldata[col_mun], numeric)
coldata <- coldata %>% mutate(Neuropath_diagnosis =
                                      case_when(Neuropath_diagnosis == "CONTR" ~ "Ctrl", 
                                                Neuropath_diagnosis == "iLBD" ~ "iLBD", 
                                                Neuropath_diagnosis == "PD" ~ "PD", 
                                                Neuropath_diagnosis == "PDD" ~ "PDD"))

### Set colors for figures
my_colors_ID <- c("0" = "#C6DBEF", "1-4" = "#6BAED6", "5" = "#08519C", "6" = "#08306B")
my_colors_neuropath <- c("Ctrl" = "#C7E9C0", "iLBD" = "#74C476", "PD" = "#238B45", "PDD" = "#00441B")


### Supplementaru Fugure 1
plot1.1 <- ggplot(coldata, aes(Braak_aSyn_groups_Plot, fill = Braak_aSyn_stage)) +
  geom_bar(stat="count", width = 0.7, color = "black") +
  theme_minimal() +
  scale_fill_manual(values = my_colors) + 
  labs(x = NULL, y = "Count", fill = "Braak Stage") +
  theme(plot.title = element_text(size=12, margin=margin(0,0,20,0)))

plot1.2 <- ggplot(coldata, aes(Braak_aSyn_groups_Plot, fill = Neuropath_diagnosis)) +
  geom_bar(stat="count", width = 0.7, color = "black") +
  theme_minimal() +
  scale_fill_manual(values = my_colors_neuropath)+ 
  labs( x = NULL, y = "Count", fill = "Diagnosis") +
  theme(plot.title = element_text(size=12, margin=margin(0,0,20,0)))

### Combine figures 
figure <- ggarrange(plot1.2, plot1.1, 
                        labels = c("a","b"),
                        ncol = 2,
                        nrow = 1)

ggsave("./Groups.svg", plot = figure, device = "svg", width = 8, height = 4, units = 'in', dpi =600)


### RIN 
## Kruskal_Wallis Test
stat.test_RIN <- coldata %>%
  kruskal_test(RIN ~ Braak_aSyn_groups) %>%
  add_significance()

plot2 <- ggplot(coldata, aes(Braak_aSyn_groups_Plot, RIN, fill = Braak_aSyn_groups_Plot)) +
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, fill= "black") +
  theme_minimal() +
  scale_fill_manual(values = my_colors_ID) + 
  ylim(5,10.5) +
  labs(title="RIN values", x = NULL, y = "RIN values") +
  theme(legend.position = "none",
        plot.title = element_text(size=16, margin=margin(0,0,20,0)))

# Mean RIN
mean(coldata$RIN) # 7.8
sd(coldata$RIN) # 0.679
min(coldata$RIN) # 5.4
max(coldata$RIN) # 9


## %GC
## Kruskal_Wallis Test
stat.test_GC <- coldata %>%
  kruskal_test(percent_gc ~ Braak_aSyn_groups) %>%
  add_significance()

# Pairwise comparisons by Dunn’s test
pwc_GC <- coldata %>% 
  dunn_test(percent_gc ~ Braak_aSyn_groups, p.adjust.method = "bonferroni") 

plot3 <- ggplot(coldata, aes(Braak_aSyn_groups_Plot, percent_gc, fill = Braak_aSyn_groups_Plot)) +
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, fill= "black") +
  theme_minimal() +
  scale_fill_manual(values = my_colors_ID) + 
  geom_segment(aes(x = "0", y = 54.375, xend = "5", yend = 54.375)) +
  geom_segment(aes(x = "0", y = 54.375, xend = "0", yend = 54.095)) +
  geom_segment(aes(x = "5", y = 54.375, xend = "5", yend = 54.095)) +
  geom_text(x = 2, y = 54.625, label = expression(paste("Dunn test, ", italic(p.adjust), " = 0.00242"))) +
  geom_segment(aes(x = "1-4", y = 53.375, xend = "5", yend = 53.375)) +
  geom_segment(aes(x = "1-4", y = 53.375, xend = "1-4", yend = 53.095)) +
  geom_segment(aes(x = "5", y = 53.375, xend = "5", yend = 53.095)) +
  geom_text(x = 2.5, y = 53.625, label = expression(paste("Dunn test, ", italic(p.adjust), " = 0.0202"))) +
  ylim(45,55) +
  labs(title="GC content (%)", x = NULL, y = "% GC", subtitle = expression(paste("Kruskal-Wallis, ",italic(X^2), " (3) = 15.5, ", italic(p), " = 0.00141, ", italic(n)," = 84")) ) +
  theme(legend.position = "none",
        plot.title = element_text(size=16, margin=margin(0,0,20,0)))

ggsave("./Braak_Stage_Groups_GC.png", plot = plot3, device = "png", width = 5, height = 5, units = 'in', dpi =300)


### Supplementary Figure 2

## Mapped % Salmon 
## Kruskal_Wallis Test
stat.test_Salmon <- coldata %>%
  kruskal_test(Salmon_percent_mapped ~ Braak_aSyn_groups) %>%
  add_significance()

plot4 <- ggplot(coldata, aes(Braak_aSyn_groups_Plot, Salmon_percent_mapped, fill = Braak_aSyn_groups_Plot)) +
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, fill= "black") +
  theme_minimal() +
  scale_fill_manual(values = my_colors_ID) + 
  ylim(40,100) +
  labs(title="Reads mapped (%) with Salmon", x = NULL, y = "% mapped reads", fill = "Braak Stage") +
  theme(
        plot.title = element_text(size=12, margin=margin(0,0,20,0)))


## Mapped % Hisat
## Kruskal_Wallis Test
stat.test_Hisat <- coldata %>%
  kruskal_test(Hisat_mapped_passed_pct ~ Braak_aSyn_groups) %>%
  add_significance()

plot5 <- ggplot(coldata, aes(Braak_aSyn_groups_Plot, Hisat_mapped_passed_pct, fill = Braak_aSyn_groups_Plot)) +
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, fill= "black") +
  theme_minimal() +
  scale_fill_manual(values = my_colors_ID) + 
  ylim(94,100) +
  labs(title="Reads mapped (%) with HISAT2", x = NULL, y = "% mapped reads", fill = "Braak stage") +
  theme(
        plot.title = element_text(size=12, margin=margin(0,0,20,0)))


### Combine figures
figure1 <- ggarrange(plot4, plot5, 
                    labels = c("a","b"),
                    ncol = 2,
                    nrow = 1,
                    common.legend = TRUE, legend = "bottom")

ggsave("./Mapped_reads.svg", plot = figure1, device = "svg", width = 8, height = 4, units = 'in', dpi =600)





