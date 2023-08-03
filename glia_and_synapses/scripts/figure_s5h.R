#### Figure S5H - Heatmap ####
#### Expression changes of neurotransmitter receptors and transporters in astrocytes ####

library(matrixStats)
library(gplots)
library(stringr)
library(here)

# human astrocyte DEGs:
h_c <- read.table(here("data/astrocytes", "Astro_human_vs_chimp_sig_genes.csv"), sep=",", header=TRUE) 
h_c <- h_c[h_c$padj<0.01 & (h_c$log2FoldChange<(-0.5) | h_c$log2FoldChange>(0.5)),] 
h_c<-h_c[-c(2:7)]
h_c$h_c <- "TRUE"
h_g <- read.table(here("data/astrocytes", "Astro_human_vs_gorilla_sig_genes.csv"), sep=",", header=TRUE) 
h_g <- h_g[h_g$padj<0.01 & (h_g$log2FoldChange<(-0.5) | h_g$log2FoldChange>(0.5)),] 
h_g<-h_g[-c(2:7)]
h_g$h_g <- "TRUE"
union_hc_hg <- merge(h_c,h_g,by="gene",all=TRUE) 
union_hc_hg[is.na(union_hc_hg)] <- "FALSE"

# human astrocyte DEGs among neurotransmiter receptors and transporters:
receptors <- read.table(here("data/astrocytes", "Astro_receptors_transporters.txt"), sep="\t", header=TRUE) 
union_hc_hg$receptors <- ifelse(union_hc_hg$gene %in% receptors$gene_name, TRUE, FALSE)
union_hc_hg <- union_hc_hg[which(union_hc_hg$receptors==TRUE),]
hdeg_receptors <- union_hc_hg$gene

# gene expression:
data <- readRDS(here("data/astrocytes","Astro1_human_vs_all_dds.RDS"))
norm.data <- DESeq2::counts(data, normalized=TRUE)
norm.dataframe <- as.data.frame(norm.data)
gene.counts <- norm.dataframe[rownames(norm.dataframe) %in% hdeg_receptors, ]
mean.count <- gene.counts
mean.count$human_mean <- rowMeans(gene.counts[grepl("human", colnames(gene.counts))])
mean.count$chimp_mean <- rowMeans(gene.counts[grepl("chimp", colnames(gene.counts))])
mean.count$gorilla_mean <- rowMeans(gene.counts[grepl("gorilla", colnames(gene.counts))])
mean.count <- mean.count[,22:24]
# fold change of gene expression:
fold.change <- mean.count
fold.change$h_fc <- log2(mean.count$human_mean/rowMeans(fold.change[c(2,3)]))
fold.change$c_fc <- log2(mean.count$chimp_mean/rowMeans(fold.change[c(1,3)]))
fold.change$g_fc <- log2(mean.count$gorilla_mean/rowMeans(fold.change[c(1,2)]))
fold.change <- fold.change[which(fold.change$human_mean>5 | fold.change$chimp_mean>5 | fold.change$gorilla_mean>5), ] 

# plot heatmap:
heatmap<-fold.change[,4:6]
heatmap<-data.matrix(heatmap)
rownames(heatmap) <- factor(rownames(heatmap), levels = c("GRIN1","GRIN2A","GRIN2B","GRIA1","GRIA2","GRIA4","GRIK4","GRIK5","GRM5","GRM7",
                                                          "GABRA2","GABRB2","GABRB3","GABRG3","GABBR1","GABBR2",
                                                          "P2RY14","ADORA1","ADORA2B",
                                                          "SLC17A7","SLC6A1"))
break_heatmap <-c(seq(-3.8,3.8,by=0.2))
color_heatmap <-colorRampPalette(c("grey25","white","#4876ff"))(38)
heatmap.2(heatmap, Colv=FALSE, Rowv=FALSE, dendrogram='none',
          col=color_heatmap, breaks=break_heatmap, scale="none", trace="none",
          key=TRUE, key.title=NA, key.ylab=NA, density.info=c("none"),
          srtCol=45, cexCol=0.6, cexRow=0.6,
          colsep=1:ncol(heatmap),rowsep=1:nrow(heatmap),
          sepcolor="black")
