#### Figure S6H #####
#### Bar graphs - gene count per protein family ####

library(ggplot2)
library(reshape2)
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

# human astrocyte DEGs among neurotransmiter receptors and transporters:
receptors <- read.table(here("data/astrocytes", "Astro_receptors_transporters.txt"), sep="\t", header=TRUE) 

# gene expression:
data <- readRDS(here("data/astrocytes","Astro1_human_vs_all_dds.RDS"))
norm.data <- DESeq2::counts(data, normalized=TRUE)
norm.dataframe <- as.data.frame(norm.data)
gene.counts <- norm.dataframe[rownames(norm.dataframe) %in% receptors$gene_name, ]
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

# gene expression of receptors and transporters:
data.table <- merge(receptors, fold.change, by.x = "gene_name", by.y = 'row.names', all.x = TRUE)
data.table <- merge(data.table, h_c, by.x = "gene_name", by.y = 'gene', all.x = TRUE)
data.table <- merge(data.table, h_g, by.x = "gene_name", by.y = 'gene', all.x = TRUE)
data.table$h_c[is.na(data.table$h_c)] <- FALSE
data.table$h_g[is.na(data.table$h_g)] <- FALSE

# summary function:
summary_function <- function(fam, data.table) {
  return(c(fam, 
           nrow(data.table[which(data.table$Family_name == fam
                                 & data.table$h_fc > 0 & (data.table$h_c == 'TRUE' | data.table$h_g == 'TRUE')
                                 & (data.table$human_mean > 5 | data.table$chimp_mean > 5 | data.table$gorilla_mean > 5)), ]),
           nrow(data.table[which(data.table$Family_name == fam
                                 & data.table$h_fc < 0 & (data.table$h_c == 'TRUE' | data.table$h_g == 'TRUE')
                                 & (data.table$human_mean > 5 | data.table$chimp_mean > 5 | data.table$gorilla_mean > 5)), ]),
           as.numeric (nrow(data.table[which(data.table$Family_name == fam
                                 & ((data.table$h_c == 'FALSE' & data.table$h_g == 'FALSE')
                                    | (data.table$human_mean < 5 & data.table$chimp_mean < 5 & data.table$gorilla_mean < 5))), ]))
  ))
}

# initialize output table:
col.names = c('protein_family', 'up', 'down', 'no')
colClasses = c(rep("character",1), rep("numeric",3))
output_table <- read.table(text = "",
                           colClasses = colClasses,
                           col.names = col.names)
# list:
families_list <- c("Glutamate ionotropic receptor AMPA", "Glutamate ionotropic receptor NMDA", "Glutamate ionotropic receptor kainate",
                   "Glutamate metabotropic receptor", "GABA type A receptor", "GABA type B receptor",
                   "Adenosine receptor", "Purinergic receptor", "Glutamate transporter",
                   "GABA transporter")
# loop:
for(fam in families_list) {
  output_table <- rbind(output_table,summary_function(fam, data.table))
}
colnames(output_table) <- col.names

# data table:
plot_data <- melt(output_table, id.vars=c("protein_family"))
colnames(plot_data) <- c('protein_family', 'expression_change', 'gene_count')
plot_data$expression_change <- factor(plot_data$expression_change, levels = c("no","up","down"))
order_protein_families <- unique(families_list)
plot_data$protein_family <- factor(plot_data$protein_family, levels = order_protein_families)
plot_data[,3] <- as.numeric(plot_data[,3])

# plot:
ggplot(data=plot_data, aes(y=gene_count, x=protein_family, fill=expression_change))+
  geom_bar(position="stack",stat="identity",color="#000000",size=0.25)+
  scale_fill_manual(values=c("#cccccc","#4876ff","#a6c4ff"))+
  scale_y_continuous(expand=c(0,0),limits=c(0,22),breaks=c(0,5,10,15,20))+
  theme_bw()+
  theme(axis.text.x=element_text(color="#000000",family="Helvetica",angle=90,vjust=0.5,hjust=1),
        axis.text.y=element_text(color="#000000",family="Helvetica"),
        axis.ticks = element_line(color="#000000", size=0.25), axis.line = element_line(color="#000000", size=0.25),
        panel.background = element_blank(), panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none")

