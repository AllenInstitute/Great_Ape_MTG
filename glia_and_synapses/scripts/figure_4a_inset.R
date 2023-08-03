#### Figure 4A - inset plot ####
#### Highly divergent genes in astrocytes ####

library(stringr)
library(dplyr)
library(here)
library(ggplot2)

# astrocyte DEG data:
astro_h_c <- read.table(here("data/Astrocytes", "Astro_human_vs_chimp_sig_genes.csv"), sep=",", header=TRUE) 
astro_h_c <- astro_h_c[astro_h_c$padj<0.01 & (astro_h_c$log2FoldChange >= 0.5 | astro_h_c$log2FoldChange <= -0.5),] 
astro_h_g <- read.table(here("data/Astrocytes", "Astro_human_vs_gorilla_sig_genes.csv"), sep=",", header=TRUE) 
astro_h_g <- astro_h_g[astro_h_g$padj<0.01 & (astro_h_g$log2FoldChange >= 0.5 | astro_h_g$log2FoldChange <= -0.5),] 
astro_c_g <- read.table(here("data/Astrocytes", "Astro_chimp_vs_gorilla_sig_genes.csv"), sep=",", header=TRUE) 
astro_c_g <- astro_c_g[astro_c_g$padj<0.01 & (astro_c_g$log2FoldChange >= 0.5 | astro_c_g$log2FoldChange <= -0.5),] 

# human-specific intersection:
astro_intersect_h <- merge(astro_h_c,astro_h_g,by="gene",all=FALSE)
astro_intersect_h <- astro_intersect_h[!(astro_intersect_h$gene %in% astro_c_g$gene),]

# chimp-specific intersection:
astro_intersect_c <- merge(astro_h_c,astro_c_g,by="gene",all=FALSE)
astro_intersect_c <- astro_intersect_c[!(astro_intersect_c$gene %in% astro_h_g$gene),]

# gorilla-specific intersection:
astro_intersect_g <- merge(astro_h_g,astro_c_g,by="gene",all=FALSE)
astro_intersect_g <- astro_intersect_g[!(astro_intersect_g$gene %in% astro_h_c$gene),]

# number of species-specific highly divergent genes:
df <- data.frame(matrix(ncol = 3, nrow = 6))
colnames(df) <- c('Specificity', 'DEG_direction', 'Number_of_genes')
df[c(1)] <- rep(c('H','C','G'),each=2)
df[c(2)] <- rep(c('up','down'),times=3)
df[1,3] <- nrow(astro_intersect_h[which(astro_intersect_h$log2FoldChange.x > log2(10) & astro_intersect_h$log2FoldChange.y > log2(10)), ])
df[2,3] <- nrow(astro_intersect_h[which(astro_intersect_h$log2FoldChange.x < -log2(10) & astro_intersect_h$log2FoldChange.y < -log2(10)), ])
df[3,3] <- nrow(astro_intersect_c[which(astro_intersect_c$log2FoldChange.x < -log2(10) & astro_intersect_c$log2FoldChange.y > log2(10)), ])
df[4,3] <- nrow(astro_intersect_c[which(astro_intersect_c$log2FoldChange.x > log2(10) & astro_intersect_c$log2FoldChange.y < -log2(10)), ])
df[5,3] <- nrow(astro_intersect_g[which(astro_intersect_g$log2FoldChange.x < -log2(10) & astro_intersect_g$log2FoldChange.y < -log2(10)), ])
df[6,3] <- nrow(astro_intersect_g[which(astro_intersect_g$log2FoldChange.x > log2(10) & astro_intersect_g$log2FoldChange.y > log2(10)), ])

# plot
df$Specificity <- factor(df$Specificity, levels = c("H","C","G"))
df$DEG_direction <- factor(df$DEG_direction, levels = c("up","down"))
plot_color = c("#4876ff","#387d7a","#ee7942")
ggplot(data=df, aes(y=Number_of_genes, x=Specificity, alpha=DEG_direction, fill=Specificity))+
  geom_bar(position="stack",stat="identity",color="#000000",size=0.25)+
  scale_fill_manual(values=plot_color)+
  scale_alpha_manual(values=c(1,0.5))+
  scale_y_continuous(expand=c(0,0),limits=c(0,345),breaks=c(0,50,100,150,200,250,300))+
  theme_bw()+
  theme(axis.text.x=element_text(color="#000000",family="Helvetica"),
        axis.text.y=element_text(color="#000000",family="Helvetica"),
        axis.ticks = element_line(color="#000000", size=0.25),
        axis.line = element_line(color="#000000", size=0.25),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none")
ggsave("astro_highly_divergent_genes.pdf")
