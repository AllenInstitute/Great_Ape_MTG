#### Figure S6B ####
#### Upset plot for microglia - great apes ####

library(stringr)
library(ggplot2)
library(UpSetR)
library(here)

# pairwise comparisons:

# human vs chimp:
h_c <- read.table(here("data/microglia", "Micro-PVM_human_vs_chimp_sig_genes.csv"), sep=",", header=TRUE) 
h_c <- h_c[h_c$padj<0.01 & (h_c$log2FoldChange<(-0.5) | h_c$log2FoldChange>(0.5)),] 
HC<-h_c[(1)]
HC$HC<-1

# human vs gorilla:
h_g <- read.table(here("data/microglia", "Micro-PVM_human_vs_gorilla_sig_genes.csv"), sep=",", header=TRUE) 
h_g <- h_g[h_g$padj<0.01 & (h_g$log2FoldChange<(-0.5) | h_g$log2FoldChange>(0.5)),] 
HG<-h_g[(1)]
HG$HG<-1

# chimp vs gorilla:
c_g <- read.table(here("data/microglia", "Micro-PVM_chimp_vs_gorilla_sig_genes.csv"), sep=",", header=TRUE) 
c_g <- c_g[c_g$padj<0.01 & (c_g$log2FoldChange<(-0.5) | c_g$log2FoldChange>(0.5)),] 
CG<-c_g[(1)]
CG$CG<-1

# matrix:
merge_matrix <- merge(HC, HG, by="gene", all=TRUE)
merge_matrix <- merge(merge_matrix, CG, by="gene", all=TRUE)
merge_matrix[is.na(merge_matrix)] <- 0

# upset plot:
upset(merge_matrix, keep.order = TRUE, order.by = "freq",
      mainbar.y.label = "Number of DEGs", sets.x.label = "DEGs per pairwise comparison",
      set_size.show=TRUE, empty.intersections=TRUE, point.size=1.9, line.size=0.7,
      sets = c("CG", "HG", "HC"))
ggsave("micro_upset_fig_s6b.pdf")
