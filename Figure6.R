setwd("/path/to/folder")

# for each metagenome, from the relative abundance species-level composition of metagenomes computed by MetaPhlAn3 with the mpa_v30_CHOCOPhlAn_201901 database, I calculated the Gut Microbiome Wellness Index 2 (accoridng to https://www.nature.com/articles/s41467-024-51651-9)

# GMWI2 command:
# wget https://raw.githubusercontent.com/danielchang2002/GMWI2/refs/heads/main/src/gmwi2_metaphlan_output.py
# wget https://raw.githubusercontent.com/danielchang2002/GMWI2/refs/heads/main/src/GMWI2/GMWI2_databases/GMWI2_model.joblib
# python3 gmwi2_metaphlan_output.py metagenome_metaphlan_output.txt GMWI2_model.joblib output_prefix

library(tidyverse)
a <- read.table("allsamples_gmwi2.txt", header = F, sep = "\t", dec =".", col.names = c("Campione", "GMWI"))
metadata <- readxl::read_xlsx("metadata.xlsx", sheet = 1)
metadata[,c(11:14)] <- str_split_fixed(metadata$Campione, "_", n = 4)
metadata <- metadata[,-c(11:14)]
metadata$Colon_part <- gsub("_.*", "", metadata$Campione)
metadata$Colon_part2 <- ifelse(grepl("^AC", metadata$Colon_part), "Ascending", ifelse(grepl("^DC", metadata$Colon_part), "Descending", "Transverse"))
metadata$Periodo <- gsub("Controllo", "Stabilization", metadata$Periodo)
metadata$Periodo <- gsub("Trattamento", "Treatment", metadata$Periodo)
metadata$Condizione <- gsub("Controllo", "Fiber", metadata$Condizione)
metadata$Condizione <- gsub("NGP", "Fiber + NGP", metadata$Condizione)
metadata$Periodo <- ifelse((grepl("solo con fibre", metadata$Note) & grepl("Washout", metadata$Periodo)), "Washout1", ifelse((grepl("feed std", metadata$Note) & grepl("Washout", metadata$Periodo)), "Washout2", metadata$Periodo))
a$Campione %in% metadata$Campione

a <- dplyr::left_join(a, metadata, by = "Campione")

a$Colon_part4 <- ifelse(grepl("_AC", a$Campione), "Ascending", ifelse(grepl("_TC", a$Campione), "Transverse", "Descending")) 
a$Periodo2 <- ifelse(grepl("Washout", a$Periodo), "Washout", a$Periodo)

ggplot(a, aes(x = Periodo, y = GMWI, color = Condizione))+
  geom_boxplot(alpha = 0)+
  geom_jitter(position = position_dodge2(width = .75))+
  ggpubr::geom_pwc()+
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = .7, alpha = .7)+
  facet_wrap(~Source)+
  labs(color = "", x = "", y = "Gut Microbiome Wellnes Index", title = "")+
  theme_classic()+
  theme(plot.title = element_text(size = 14, face = "italic", hjust = .5), 
        strip.text = element_text(size = 10, face = "bold"))
