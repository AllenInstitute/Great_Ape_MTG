#### Figure S6F ####
#### Heatmap of human microglia divergent genes in SynGO database ####

library(matrixStats)
library(gplots)
library(stringr)
library(here)

# human astrocyte DEGs:
h_c <- read.table(here("data/microglia", "Micro-PVM_human_vs_chimp_sig_genes.csv"), sep=",", header=TRUE) 
h_c <- h_c[h_c$padj<0.01 & (h_c$log2FoldChange<(-0.5) | h_c$log2FoldChange>(0.5)),] 
h_c<-h_c[-c(2:7)]
h_c$h_c <- "TRUE"
h_g <- read.table(here("data/microglia", "Micro-PVM_human_vs_gorilla_sig_genes.csv"), sep=",", header=TRUE) 
h_g <- h_g[h_g$padj<0.01 & (h_g$log2FoldChange<(-0.5) | h_g$log2FoldChange>(0.5)),] 
h_g<-h_g[-c(2:7)]
h_g$h_g <- "TRUE"
union_hc_hg <- merge(h_c,h_g,by="gene",all=TRUE) 
union_hc_hg[is.na(union_hc_hg)] <- "FALSE"

# human microglia DEGs in SynGO:
syngo_terms <- read.table(here("data/syngo", "syngo_level0_terms_id.txt"), sep="\t", header=TRUE) 
syngo_genes <- read.table(here("data/syngo", "syngo_level0_ontologies.txt"), sep="\t", header=TRUE) 
syngo_terms_to_test <- syngo_terms$syngo_term_name
syngo <- unlist(str_split(syngo_genes$genes...hgnc_symbol[syngo_genes$GO.term.name == syngo_terms_to_test], ";"))
union_hc_hg$SynGO <- ifelse(union_hc_hg$gene %in% syngo, TRUE, FALSE)
union_hc_hg <- union_hc_hg[which(union_hc_hg$SynGO==TRUE),]
hdeg_syngo <- union_hc_hg$gene

# gene expression:
data <- readRDS(here("data/microglia","MicroPVM1_human_vs_all_dds.RDS"))
norm.data <- DESeq2::counts(data, normalized=TRUE)
norm.dataframe <- as.data.frame(norm.data)
gene.counts <- norm.dataframe[rownames(norm.dataframe) %in% hdeg_syngo, ]
mean.count <- gene.counts
mean.count$human_mean <- rowMeans(gene.counts[grepl("human", colnames(gene.counts))])
mean.count$chimp_mean <- rowMeans(gene.counts[grepl("chimp", colnames(gene.counts))])
mean.count$gorilla_mean <- rowMeans(gene.counts[grepl("gorilla", colnames(gene.counts))])
mean.count <- mean.count[,23:25]
# fold change of gene expression:
fold.change <- mean.count
fold.change$g_fc <- log2(mean.count$gorilla_mean/rowMeans(fold.change[c(1,2)]))
fold.change$c_fc <- log2(mean.count$chimp_mean/rowMeans(fold.change[c(1,3)]))
fold.change$h_fc <- log2(mean.count$human_mean/rowMeans(fold.change[c(2,3)]))
fold.change <- fold.change[which(fold.change$human_mean>5 | fold.change$chimp_mean>5 | fold.change$gorilla_mean>5), ] 
fold.change <- fold.change[which(abs(fold.change$h_fc)>1.3), ] 
# genes with low reads in 10X dataset
low_reads <- c("DISC1", "RPL17", "ATP6V0C", "PPP3R1", "PPP2CA", "RPL26", "VAMP2", "BLOC1S6")
fold.change <- fold.change[!(rownames(fold.change) %in% low_reads), ] 
# Inf values:
fold.change$h_fc[fold.change$h_fc == "-Inf"] <- -6.189618 

# plot heatmap:
heatmap<-fold.change[,4:6]
heatmap<-heatmap[order(-heatmap$h_fc),] 
heatmap<-data.matrix(heatmap)
break_heatmap <-c(seq(-6.2,6.2,by=0.2))
color_heatmap <-colorRampPalette(c("grey25","white","#4876ff"))(62)
heatmap.2(heatmap, Colv=FALSE, Rowv=FALSE, dendrogram='none',
          col=color_heatmap, breaks=break_heatmap, scale="none", trace="none",
          key=TRUE, key.title=NA, key.ylab=NA, density.info=c("none"),
          srtCol=45, cexCol=0.6, cexRow=0.6,
          colsep=1:ncol(heatmap),rowsep=1:nrow(heatmap),
          sepcolor="black")
