# Starting ====
setwd("path")

# importing libraries
library(tidyverse)
library(tidygraph)
library(vegan)
library(patchwork)
library(corrr)
library(ggraph)

# Taxonomic analysis ====
species <- read.table("analysis/ngp_mpa_species.txt", header = T, sep = "\t", dec = ".")
rownames(species) <- species[,1]
species <- as.data.frame(t(species[,-1]))
rownames(species) %in% metadata$Campione
rownames(species) <- gsub("_4_", "_04_", rownames(species))
rownames(species) <- gsub("_6_", "_06_", rownames(species)) 
rownames(species) <- gsub("_7_", "_07_", rownames(species)) 
rownames(species) <- gsub("_8_", "_08_", rownames(species)) 
rownames(species) <- gsub("_9_", "_09_", rownames(species)) 
rownames(species)[which(!(rownames(species) %in% metadata$Campione))] 

species <- species[order(rownames(species)),]
metadata <- metadata[order(metadata$Campione),]

names_species <- colnames(species)
names_meta <- colnames(metadata)

species$Campione <- rownames(species)
species <- merge(species, metadata, by = "Campione", all.x = T)

species_ngps$Colon_part4 <- ifelse(species_ngps$Colon_part2 == "Ascending", "Ascending", "Transverse + Descending")
#### B. salyersiae ====
ggplot(filter(species_ngps, taxon == "Bacteroides_salyersiae"), aes(x = Periodo, y = Relative_abundance, color = Condizione))+
  geom_boxplot(alpha = 0)+
  geom_jitter(position = position_dodge2(width = .75))+
  ggpubr::geom_pwc(y.position = .6)+
  scale_y_sqrt()+
  facet_wrap(Source~factor(Colon_part4, levels = c("Ascending", "Transverse + Descending")))+
  labs(color = "", x = "", y = "Relative abundance", title = "Bacteroides salyersiae")+
  theme_classic()+
  theme(plot.title = element_text(size = 14, face = "italic", hjust = .5), 
        strip.text = element_text(size = 10, face = "bold"))

#### B. uniformis ====
ggplot(filter(species_ngps, taxon == "Bacteroides_uniformis"), aes(x = Periodo, y = Relative_abundance, color = Condizione))+
  geom_boxplot(alpha = 0)+
  geom_jitter(position = position_dodge2(width = .75))+
  ggpubr::geom_pwc(y.position = 3.1)+
  scale_y_sqrt()+
  facet_wrap(Source~factor(Colon_part4, levels = c("Ascending", "Transverse + Descending")))+
  labs(color = "", x = "", y = "Relative abundance", title = "Bacteroides uniformis")+
  theme_classic()+
  theme(plot.title = element_text(size = 14, face = "italic", hjust = .5), 
        strip.text = element_text(size = 10, face = "bold"))

#### B. thetaiotaomicron ====
ggplot(filter(species_ngps, taxon == "Bacteroides_thetaiotaomicron"), aes(x = Periodo, y = Relative_abundance, color = Condizione))+
  geom_boxplot(alpha = 0)+
  geom_jitter(position = position_dodge2(width = .75))+
  ggpubr::geom_pwc(y.position = 4.2)+
  scale_y_sqrt()+
  facet_wrap(Source~factor(Colon_part4, levels = c("Ascending", "Transverse + Descending")))+
  labs(color = "", x = "", y = "Relative abundance", title = "Bacteroides thetaiotaomicron")+
  theme_classic()+
  theme(plot.title = element_text(size = 14, face = "italic", hjust = .5), 
        strip.text = element_text(size = 10, face = "bold"))
