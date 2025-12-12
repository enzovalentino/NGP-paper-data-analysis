setwd("/path/to/folder") 

### PANEL A
scfa <- readxl::read_xlsx("scfa.xlsx", sheet = 1)
library(tidyverse)
scfa <- drop_na(scfa, Sample)

library(pheatmap)
scfa <- as.data.frame(scfa)
rownames(scfa) <- scfa$Sample
pheatmap(scfa[,c(2:4)], 
         cluster_rows = FALSE,
         cluster_cols = FALSE, 
         clustering_distance_rows = "canberra", 
         clustering_distance_cols = "correlation", 
         clustering_method = "ward.D", 
         gaps_row = c(5, 8, 10, 12, 14),
         show_colnames = TRUE, 
         labels_row = gsub("\\+.*", "", scfa$Sample), 
         display_numbers = T, 
         number_color = "black", 
         drop_levels = T, 
         border_color = "black",
         fontsize = 12)

# saved as pdf named scfa.pdf, cairo disabled, 6 x 7.5 inches, potrait mode

### PANEL B
