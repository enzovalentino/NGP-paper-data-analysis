setwd("/path/to/files")
# With this command I process the files obtained by eggNOG mapper and KEGGaNOG (https://github.com/iliapopov17/KEGGaNOG) on the metagenomes

### PANEL B

# importing and reading archive
library(tidyverse)
library(pheatmap)
library(archive)
library(ggpubr)

file <- "kegganog_tables.tar.gz"
archive(file)

lista <- list()
for (i in 1:dim(archive(file))[1]) {
lista[[i]] <- read.table(archive_read(file, file = i), header = T, sep = "\t", dec =".", check.names = F)
names(lista)[i] <- gsub("_path.*", "", archive(file)[i,"path"])
}

pathways <- Reduce(
  function(x, y, ...) dplyr::full_join(x, y, ...),
  lista)
colnames(pathways)[1] <- "Campione"

rm(lista)
# 178 pathways

names_pathways <- colnames(pathways)[2:ncol(pathways)]

# loading metadata
metadata <- readxl::read_xlsx("metadata.xlsx", sheet = 1)
metadata[,c(11:14)] <- str_split_fixed(metadata$Campione, "_", n = 4)
metadata$Campione <- paste0(metadata$V2, "_", metadata$V1, "_", metadata$V3, "_", metadata$V4)
metadata <- metadata[,-c(11:14)]
metadata$Colon_part <- gsub("_.*", "", metadata$Campione)
metadata$Colon_part2 <- ifelse(grepl("^AC", metadata$Colon_part), "Ascending", ifelse(grepl("^DC", metadata$Colon_part), "Descending", "Transverse"))
metadata$Periodo <- gsub("Controllo", "Stabilization", metadata$Periodo)
metadata$Periodo <- gsub("Trattamento", "Treatment", metadata$Periodo)
metadata$Condizione <- gsub("Controllo", "Fiber", metadata$Condizione)
metadata$Condizione <- gsub("NGP", "Fiber + NGP", metadata$Condizione)
metadata$Campione <- gsub("M_0", "M_", metadata$Campione)
metadata$Campione <- gsub("L_0", "L_", metadata$Campione)

# merging
pathways <- dplyr::full_join(metadata, pathways)
pathways <- as.data.frame(pathways)
rownames(pathways) <- pathways$Campione

# clustering
pheatmap(t(dplyr::select(pathways, all_of(names_pathways))), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         clustering_distance_rows = "binary", 
         clustering_distance_cols = "binary", 
         clustering_method = "ward.D", 
         annotation_col=pathways[,c(3,4,5)], 
         cutree_cols = 2, 
         border_color = NA, 
         fontsize = 7)

## Lumen ====
### clustering only lumen and ascending ====
tmp <- dplyr::filter(pathways, Source == "Lumen" & Colon_part2 == "Ascending") %>% dplyr::select(all_of(names_pathways))
torm <- c(colnames(tmp)[which(colSums(tmp) == 0)], colnames(tmp)[which(colSums(tmp) == nrow(tmp))])
pheatmap(t(dplyr::filter(pathways, Source == "Lumen" & Colon_part2 == "Ascending") %>% dplyr::select(all_of(names_pathways)) %>% dplyr::select(-all_of(torm))), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         clustering_distance_rows = "manhattan", 
         clustering_distance_cols = "manhattan", 
         clustering_method = "ward.D", 
         annotation_color = list(Condizione = c("Fiber" = "#F8766D", "Fiber + NGP" = "#00BFC4"), Periodo = c("Stabilization" = "#00a1d5", "Treatment" = "#df8f44", "Washout" = "#6a6599")),
         annotation_col=dplyr::filter(pathways, Source == "Lumen")[,c(3,4,5)], 
         cutree_cols = 2, 
         show_colnames = F,
         border_color = NA, 
         fontsize = 9)


### PANEL A
ggplot(dplyr::filter(pathways, Source == "Lumen" & Colon_part2 == "Ascending"), aes(x = Periodo, y = `mevalonate pathway`, color = Periodo))+
  geom_boxplot(alpha = 0, size = 0.8)+
  geom_jitter()+
  ggpubr::geom_pwc()+
  facet_wrap(~Condizione)+
  labs(x = "", y = "mevalonate pathway\ncompleteness")+
  theme_classic2()+
  theme(strip.text = element_text(face = "bold", size = 11))+
  scale_x_discrete(labels = c("Stabilization", "Treatment", "Washout (1+2)"))+
  scale_color_manual(values = c("Stabilization" = "#00a1d5", "Treatment" = "#df8f44", "Washout" = "#6a6599"))+
  guides(color = 'none')



### PANEL C
tmp <- dplyr::filter(pathways, Source == "Lumen" & Colon_part2 == "Ascending") %>% dplyr::select(all_of(names_pathways))
torm <- c(colnames(tmp)[which(colSums(tmp) == 0)], colnames(tmp)[which(colSums(tmp) == nrow(tmp))])

pcoa1 <- cmdscale(vegan::vegdist(dplyr::filter(pathways, Source == "Lumen" & Colon_part2 == "Ascending") %>% dplyr::select(all_of(names_pathways)) %>% dplyr::select(-all_of(torm)), method = "bray"), k = 3, eig = T)
ggplot(as.data.frame(pcoa1$points), aes(x = V1, y = V2, color = dplyr::filter(pathways, Source == "Lumen" & Colon_part2 == "Ascending")$Periodo, shape = dplyr::filter(pathways, Source == "Lumen" & Colon_part2 == "Ascending")$Condizione))+
  geom_point(size = 4)+
  scale_color_manual(values = c("Stabilization" = "#00a1d5", "Treatment" = "#df8f44", "Washout" = "#6a6599"))+
  labs(color = "", shape ="", 
       x = "", 
       y = "")
