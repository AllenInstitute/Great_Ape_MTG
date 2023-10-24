#### Figure S6A ####
#### Upset plot for microglia - primates ####

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

# human vs rhesus:
h_r <- read.table(here("data/microglia", "Micro-PVM_human_vs_rhesus_sig_genes.csv"), sep=",", header=TRUE) 
h_r <- h_r[h_r$padj<0.01 & (h_r$log2FoldChange<(-0.5) | h_r$log2FoldChange>(0.5)),] 
HR<-h_r[(1)]
HR$HR<-1

# chimp vs rhesus:
c_r <- read.table(here("data/microglia", "Micro-PVM_chimp_vs_rhesus_sig_genes.csv"), sep=",", header=TRUE) 
c_r <- c_r[c_r$padj<0.01 & (c_r$log2FoldChange<(-0.5) | c_r$log2FoldChange>(0.5)),] 
CR<-c_r[(1)]
CR$CR<-1

# gorilla vs rhesus:
g_r <- read.table(here("data/microglia", "Micro-PVM_gorilla_vs_rhesus_sig_genes.csv"), sep=",", header=TRUE) 
g_r <- g_r[g_r$padj<0.01 & (g_r$log2FoldChange<(-0.5) | g_r$log2FoldChange>(0.5)),] 
GR<-g_r[(1)]
GR$GR<-1

# human vs marmoset:
h_m <- read.table(here("data/microglia", "Micro-PVM_human_vs_marmoset_sig_genes.csv"), sep=",", header=TRUE) 
h_m <- h_m[h_m$padj<0.01 & (h_m$log2FoldChange<(-0.5) | h_m$log2FoldChange>(0.5)),] 
HM<-h_m[(1)]
HM$HM<-1

# chimp vs marmoset:
c_m <- read.table(here("data/microglia", "Micro-PVM_chimp_vs_marmoset_sig_genes.csv"), sep=",", header=TRUE) 
c_m <- c_m[c_m$padj<0.01 & (c_m$log2FoldChange<(-0.5) | c_m$log2FoldChange>(0.5)),] 
CM<-c_m[(1)]
CM$CM<-1

# gorilla vs marmoset:
g_m <- read.table(here("data/microglia", "Micro-PVM_gorilla_vs_marmoset_sig_genes.csv"), sep=",", header=TRUE) 
g_m <- g_m[g_m$padj<0.01 & (g_m$log2FoldChange<(-0.5) | g_m$log2FoldChange>(0.5)),] 
GM<-g_m[(1)]
GM$GM<-1

# rhesus vs marmoset:
r_m <- read.table(here("data/microglia", "Micro-PVM_rhesus_vs_marmoset_sig_genes.csv"), sep=",", header=TRUE) 
r_m <- r_m[r_m$padj<0.01 & (r_m$log2FoldChange<(-0.5) | r_m$log2FoldChange>(0.5)),] 
RM<-r_m[(1)]
RM$RM<-1

# matrix:
merge_matrix <- merge(HC, HG, by="gene", all=TRUE)
merge_matrix <- merge(merge_matrix, CG, by="gene", all=TRUE)
merge_matrix <- merge(merge_matrix, HR, by="gene", all=TRUE)
merge_matrix <- merge(merge_matrix, CR, by="gene", all=TRUE)
merge_matrix <- merge(merge_matrix, GR, by="gene", all=TRUE)
merge_matrix <- merge(merge_matrix, HM, by="gene", all=TRUE)
merge_matrix <- merge(merge_matrix, CM, by="gene", all=TRUE)
merge_matrix <- merge(merge_matrix, GM, by="gene", all=TRUE)
merge_matrix <- merge(merge_matrix, RM, by="gene", all=TRUE)
merge_matrix[is.na(merge_matrix)] <- 0

# upset plot:
upset(merge_matrix,  nsets = 10, nintersects = 60, number.angles=0,
      keep.order = TRUE, order.by = "freq",
      mainbar.y.label = "Number of DEGs", sets.x.label = "DEGs per pairwise comparison",
      set_size.show=TRUE, empty.intersections=TRUE,
      point.size=1.9, line.size=0.7,
      mb.ratio = c(0.50,0.50))
ggsave("micro_upset_fig_s6a.pdf")
