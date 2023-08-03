#### Figure panels: S6L ####
#### Gene expression in oligodendrocytes ####

library(ggplot2)
library(here)

# gene expression:
data <- readRDS(here("data/oligodendrocytes","Oligo1_human_vs_all_dds.RDS"))
norm.data <- DESeq2::counts(data, normalized=TRUE)
gene = "CNTNAP2" ##gene of interest
gene.counts <- data.frame(norm.data[rownames(norm.data) == gene, ])
colnames(gene.counts) <- "expr"
gene.counts$species[grepl("human", rownames(gene.counts))] <- "H"
gene.counts$species[grepl("chimp", rownames(gene.counts))] <- "C"
gene.counts$species[grepl("gorilla", rownames(gene.counts))] <- "G"
gene.counts$species[grepl("rhesus", rownames(gene.counts))] <- "R"
gene.counts$species[grepl("marmoset", rownames(gene.counts))] <- "M"
species <- c("H", "C", "G", "R", "M")
gene.counts$species <- factor(gene.counts$species, levels=species)

# boxplot:
ggplot(gene.counts, aes(x=species, y=expr, fill=species))+
  geom_boxplot(outlier.shape = NA, color="#000000",size=0.25)+
  scale_fill_manual(values=c("#4876ff", "#387d7a", "#ee7942", "#006400", "#8b1c62"))+
  theme_bw()+
  theme(axis.text = element_text(color="#000000",family="Helvetica"),
        axis.ticks = element_line(color="#000000", size=0.25), axis.line = element_line(color="#000000", size=0.25), 
        panel.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none")+
  labs(title = gene)
ggsave("oligo_gene_expression.pdf")
