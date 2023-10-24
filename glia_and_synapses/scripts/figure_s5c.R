#### Figure S5C #####
#### Expression of genes associated with perisynaptic astrocytic processes ####

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
union_hc_hg$human_divergence <- ifelse((union_hc_hg$h_c == 'TRUE' | union_hc_hg$h_g == 'TRUE'), TRUE, FALSE)
union_hc_hg <- union_hc_hg[,c(1,4)]

# human astrocyte DEGs and genes associated with perisynaptic astrocytic processes:
astro_cleft <- read.table(here("data/astrocytes", "Takano_2020_astrocyte_cleft_proteome.txt"), sep="\t", header=TRUE) 

# gene expression:
data <- readRDS(here("data/astrocytes","Astro1_human_vs_all_dds.RDS"))
norm.data <- DESeq2::counts(data, normalized=TRUE)
norm.dataframe <- as.data.frame(norm.data)
gene.counts <- norm.dataframe[rownames(norm.dataframe) %in% astro_cleft$gene, ]
mean.count <- gene.counts
mean.count$human_mean <- rowMeans(gene.counts[grepl("human", colnames(gene.counts))])
mean.count$chimp_mean <- rowMeans(gene.counts[grepl("chimp", colnames(gene.counts))])
mean.count$gorilla_mean <- rowMeans(gene.counts[grepl("gorilla", colnames(gene.counts))])
mean.count <- mean.count[,22:24]
# fold change of gene expression:
fold.change <- mean.count
fold.change$g_fc <- log2(mean.count$gorilla_mean/rowMeans(fold.change[c(1,2)]))
fold.change$c_fc <- log2(mean.count$chimp_mean/rowMeans(fold.change[c(1,3)]))
fold.change$h_fc <- log2(mean.count$human_mean/rowMeans(fold.change[c(2,3)]))
fold.change <- fold.change[which(fold.change$human_mean>5 | fold.change$chimp_mean>5 | fold.change$gorilla_mean>5), ] 

# plot heatmap of human divergent genes:
gene_set_astro_cleft_dirvegent <- merge(fold.change,union_hc_hg,by.x="row.names",by.y="gene",all.x=TRUE)
row.names(gene_set_astro_cleft_dirvegent) <- gene_set_astro_cleft_dirvegent$Row.names
gene_set_astro_cleft_dirvegent <- gene_set_astro_cleft_dirvegent[which(gene_set_astro_cleft_dirvegent$human_divergence == "TRUE"), ] 
hm_divergent<-gene_set_astro_cleft_dirvegent[,5:7]
hm_divergent<-hm_divergent[order(-hm_divergent$h_fc),] 
hm_divergent<-data.matrix(hm_divergent)
break_heatmap <-c(seq(-4.2,4.2,by=0.2))
color_heatmap <-colorRampPalette(c("grey25","white","#4876ff"))(42)
heatmap.2(hm_divergent, Colv=FALSE, Rowv=FALSE, dendrogram='none',
          col=color_heatmap, breaks=break_heatmap, scale="none", trace="none",
          key=TRUE, key.title=NA, key.ylab=NA, density.info=c("none"),
          srtCol=45, cexCol=0.6, cexRow=0.6,
          colsep=1:ncol(hm_divergent),rowsep=1:nrow(hm_divergent),
          sepcolor="black")

# plot heatmap of genes with shared expression:
gene_set_astro_cleft_shared <- merge(fold.change,union_hc_hg,by.x="row.names",by.y="gene",all.x=TRUE)
row.names(gene_set_astro_cleft_shared) <- gene_set_astro_cleft_shared$Row.names
gene_set_astro_cleft_shared<-gene_set_astro_cleft_shared[is.na(gene_set_astro_cleft_shared$human_divergence),]
hm_shared<-gene_set_astro_cleft_shared[,5:7]
hm_shared<-hm_shared[order(-hm_shared$h_fc),] 
hm_shared<-data.matrix(hm_shared)
break_heatmap <-c(seq(-4.2,4.2,by=0.2))
color_heatmap <-colorRampPalette(c("grey25","white","#4876ff"))(42)
heatmap.2(hm_shared, Colv=FALSE, Rowv=FALSE, dendrogram='none',
          col=color_heatmap, breaks=break_heatmap, scale="none", trace="none",
          key=TRUE, key.title=NA, key.ylab=NA, density.info=c("none"),
          srtCol=45, cexCol=0.6, cexRow=0.6,
          colsep=1:ncol(hm_shared),rowsep=1:nrow(hm_shared),
          sepcolor="black")
