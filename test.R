data <- read_csv("bf591-finalproject-mproberts99/data/GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.csv")
BiocManager::install("fgsea")
library(fgsea)
library(tidyverse)

res2 <- data %>% 
  dplyr::select(symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(stat))

ranks <- deframe(res2)

pathways.hallmark <- gmtPathways("/usr4/bf527/monicapr/bf591-finalproject-mproberts99/data/h.all.v7.5.1.symbols.gmt")

fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks)

write_csv(fgseaRes, file = 'fgsea_results_GSE64810.csv')

filtered <- fgseaRes %>% arrange(padj) %>% slice_head(n=50) %>% mutate(status= ifelse(NES > 0, 'positive', 'negative'))
filtered$pathway <- gsub("\\_", " ", filtered$pathway)
filtered <- filtered %>% arrange(status)
pathways <- factor(filtered$pathway)
filtered$pathway <- factor(filtered$pathway, levels=unique(pathways))
status <- factor(filtered$status)
filtered$status <- factor(filtered$status, levels=unique(status))
filtered %>% ggplot(aes(x=stringr::str_wrap(pathway, 40), y=NES, fill=status)) +
  geom_col() +
  coord_flip() +
  theme_light() +
  theme(axis.text = element_text(size = 5)) +
  xlab("")

samples <- read_csv('/usr4/bf527/monicapr/bf591-finalproject-mproberts99/data/SraRunTable.csv')

means <- sapply(samples, mean)
sd <- sapply(samples, sd)


