setwd("/path/to/folder/with/files")
## PANEL A
x <- read.table("~/Library/Mobile Documents/com~apple~CloudDocs/PhD/Progetti/Genomi_Alessia/polyphenols/anaerobe_genomes_polyphenols.txt")

#library(tidyverse)
library(reshape2)
library(pheatmap)
library(vegan)
#library(Rtsne)
library(ggplot2)
library(patchwork)
library(dplyr)

# calcolo coverage (length - gaps / slen) 
x$V14 <- ((x$V4 - x$V5) / x$V7)*100

# grouping by geneid, sorting by coverage and pident, extracting the best matches and filtering out genes with pident and coverage < 70
x_filt <- x %>% group_by(V1) %>% arrange(desc(V14), desc(V3)) %>% filter(row_number() == 1) %>% filter(V3 >= 30 && V14 >= 50)
rm(x)

# NO GENES!

# per rileggere x_filt, per ricominciare l'analisi:
#x_filt <- read.table("~/Library/Mobile Documents/com~apple~CloudDocs/PhD/Progetti/Genomi_Alessia/polyphenols/anaerobe_genomes_polyphenols.txt", header = T, sep = "\t", dec = ".", quote = "")

x_filt$V15 <- gsub("_gene_id.*", "", x_filt$V1) # sampleID
colnames(x_filt) <- c("qseqid", "sseqid", "pident", "length", "gaps", "qlen", "slen","qstart", "qend", "sstart", "send", "evalue", "bitscore", "coverage", "GenomeID")

# salvo l'oggetto x_filt per sicurezza
write.table(x_filt, "~/Library/Mobile Documents/com~apple~CloudDocs/PhD/Progetti/Genomi_Alessia/polyphenols/anaerobe_genomes_polyphenols_matches_filt_sorted.txt", row.names = F, col.names = T, sep = "\t", dec = ".", quote = F)

# ora collasso i geni x i soggetti

pres_abs <- x_filt %>% group_by(sseqid, GenomeID) %>% summarise(n = n())
pres_abs_m <- reshape2::dcast(pres_abs, V2 ~ V15)
pres_abs_m[is.na(pres_abs_m)] <- 0
rownames(pres_abs_m) <- pres_abs_m[,1]
pres_abs_m <- as.data.frame(pres_abs_m[,-1])

# binary
pres_abs_m_b <- ifelse(pres_abs_m > 0, 1, 0)





## PANEL B
library(archive)
library(tidyverse)

# for each genomes, aminoacidic sequences of genes predicted with Prokka were aligned to the CAZy database (version 07262023) throug diamond
# all the diamond output files from all the genomes were compressed in a tar.gz archive, that was imported and processed in R

file <- "ngp_shime_mags_cazy.tar.gz"
dimens <- archive(file)

# uploading cazy matches
cazy_matches <- list()
for (i in 1:dim(dimens)[1]) {
 cazy_matches[[i]] <- read.table(archive_read(file, file = i), header = F, sep = "\t", dec =".", col.names = c("qseqid", "sseqid", "pident", "length", "gaps", "qlen", "slen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
 names(cazy_matches)[i] <- gsub("_cz.txt", "", dimens[i,1])
}

# colnames of single data frames
cazy_matches <- Reduce( 
 function(x, y, ...) rbind(x, y, ...),
 cazy_matches
)

cazy_matches$SampleID <- gsub("_gene_.*", "", cazy_matches$qseqid)

# filtering matches
cazy_matches <- cazy_matches %>% mutate(coverage = (length - gaps)/slen*100) # to calculate coverage
cazy_matches_filt <- cazy_matches %>% group_by(qseqid) %>% arrange(desc(coverage), desc(pident)) %>% filter(row_number()==1) %>% filter(pident >= 70 & coverage >= 70)
cazy_matches_filt <- read.table("cazy_ngp_genomes_filtmatches.txt", header = T, sep = "\t", dec = ".", quote = "")
cazy_matches_filt$cazyme <- gsub('^.*?\\|', '', cazy_matches_filt$sseqid)
cazy_matches_filt$bin_id <- gsub('_L_[^_]+$', '', cazy_matches_filt$qseqid)

cazy_matches_filt_matrix <- cazy_matches_filt %>% group_by(bin_id, cazyme) %>% count() %>% pivot_wider(names_from = cazyme, values_from = n)
cazy_matches_filt_matrix[is.na(cazy_matches_filt_matrix)] <- 0
names_cazy <- colnames(cazy_matches_filt_matrix)[-1]

metadata <- readxl::read_xlsx("NGP_SHIME_MAGs_metadata.xlsx", sheet = 1)
 
cazy_matches_filt_matrix <- left_join(cazy_matches_filt_matrix, metadata, by = "bin_id")
cazy_matches_filt_matrix <- as.data.frame(cazy_matches_filt_matrix)
rownames(cazy_matches_filt_matrix) <- cazy_matches_filt_matrix$bin_id
pheatmap::pheatmap(t(cazy_matches_filt_matrix %>% dplyr::select(., all_of(names_cazy)) %>% mutate(across(all_of(names_cazy), ~ifelse(.x > 0, 1, 0)))), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         clustering_distance_rows = "manhattan", 
         clustering_distance_cols = "manhattan", 
         clustering_method = "ward.D", 
         show_colnames = F, 
         show_rownames = F,
         drop_levels = T, 
         cutree_cols = 2, 
         border_color = NA, 
         fontsize = 13)
