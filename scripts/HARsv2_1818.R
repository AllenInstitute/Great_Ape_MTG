# plots for Fig. 6E, S24

library(universalmotif)
library(ggseqlogo)
library(motifStack)
library(cowplot)
library(GenomicRanges)
library(karyoploteR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyverse)
library(here)

# defaults ----
HAR.ID <- "HARsv2_1818"
font <- "Helvetica"
human_color <- "#4876FF"
chimp_color <- "#387D7A"
gorilla_color <- "#EE7942"
rhesus_color <- "#006400"
marmoset_color <- "#8B1C62"
species_colors <- c(human_color, chimp_color, gorilla_color, rhesus_color, marmoset_color)
EN_color <- "#8B0000"
IN_color <- "#00008B"
NN_color <- "#006400"
ct_colors <- c("EN" = EN_color, "IN" = IN_color, "NN" = NN_color)
cluster_color_table <- read.table(here("data/HARs_hCONDELs_HAQERs", "colors_for_figure.txt"), sep="\t", header=FALSE)
cluster_colors <- cluster_color_table$V3
names(cluster_colors) <- cluster_color_table$V1
cluster_names <- cluster_color_table$V2
names(cluster_names) <- cluster_color_table$V1
cluster_order_table <- read.table(here("data/HARs_hCONDELs_HAQERs", "cross_species_clusters_order.txt"), 
                                  sep="\t", header=FALSE)
cluster_order_table <- merge(cluster_order_table, cluster_color_table, by.x="V2", by.y="V2")
cluster_order <- cluster_order_table$V1.y[order(cluster_order_table$V1.x)]
species_short_names <- c("H", "C", "G", "R", "M")
names(species_short_names) <- c("Human", "Chimp", "Gorilla", "Rhesus", "Marmoset")
inh_neurons <- cluster_order_table$V1.y[cluster_order_table$V1.x >=20 & cluster_order_table$V1.x <= 50]

# plot Girskis et al., Neuron 2021 MPRA data ----
HAR.MPRA.data = read.table(here("data/HARs_hCONDELs_HAQERs/HARsv2_1818_Girskis2021_MPRA_SHSY5Y_data.txt"), 
                           sep = "\t", header = TRUE) %>%
  mutate(expr = as.numeric(expr), 
         species = factor(species, levels=c("Human", "Chimp")))
ggplot(HAR.MPRA.data, aes(x=species, y=expr, fill=species)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c(human_color, chimp_color), guide="none") +
  geom_jitter(width=0.25) + 
  theme_classic(base_size=20) +
  scale_x_discrete(labels=species_short_names) +
  theme(legend.position="none", axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        text = element_text(family = font)) + 
  labs(x="", y=paste0(HAR.ID, " enhancer activity"))
ggsave("HARsv2_1818-Girskis2021-MPRA-SHSY5Y.pdf")

# plot PTPRG data from pseudobulk counts ----
PTPRG_data = read.table(here("data/HARs_hCONDELs_HAQERs", "PTPRG_norm_expr.txt"), 
                        sep = "\t", header = TRUE)
ct_of_interest = c("L5ET2", "MicroPVM1", "Vip2", "Vip6")

ggplot(PTPRG_data[PTPRG_data$cell_type == "L5ET2", ] %>%
         mutate(species = factor(species, levels = names(species_short_names))), 
       aes(x=species, y=expr, fill=species)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.25) +
  theme_classic(base_size=20) +
  labs(x="", y=paste0(gene, " expression")) +
  scale_fill_manual(values=species_colors) +
  scale_x_discrete(labels=species_short_names) +
  theme(legend.position = "none", text = element_text(family = font),
        axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=20))
ggsave("PTPRG_L5ET2_byspecies.pdf")

# Fig. S24A
ggplot(PTPRG_data[PTPRG_data$species == "Human", ] %>%
         mutate(cell_type = factor(cell_type, levels=cluster_order)), 
       aes(x=cell_type, y=expr, fill=cell_type)) +
  geom_boxplot() +
  theme_classic(base_size=22) +
  labs(x="", y=paste0(gene, " expression")) +
  scale_fill_manual(values=cluster_colors) +
  scale_x_discrete(labels=cluster_names) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=13), legend.position = "none", 
        text = element_text(family = font))
ggsave("PTPRG_human.pdf")

# Fig. S24B
ggplot(PTPRG_data[PTPRG_data$cell_type %in% ct_of_interest, ], 
       aes(x=species, y=expr, fill=species)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.25) +
  theme_bw(base_size=22) +
  labs(x="", y=paste0(gene, " expression")) +
  scale_fill_manual(values=species_colors) + 
  theme(legend.position = "none", text = element_text(family = font), 
        axis.text.x = element_text(size=13)) +
  facet_wrap(~cell_type, scales="free_y", labeller = labeller(cell_type = cluster_names))
ggsave("PTPRG_L5ET2_MicroPVM1_Vip2_Vip6_byspecies.pdf")

# plot TWIST1 data from pseudobulk counts ----
TWIST1_data = read.table(here("data/HARs_hCONDELs_HAQERs", "TWIST1_norm_expr.txt"), 
                        sep = "\t", header = TRUE)

chimp_TWIST1_data <- TWIST1_data %>%
  mutate(cell_type = factor(cell_type, levels=cluster_order)) %>%
  filter(species == "Chimp")
chimp_TWIST1_data.sub = chimp_TWIST1_data %>%
  filter(cell_type %in% c("L5ET2", inh_neurons, "MicroPVM1")) %>%
  mutate(cell_type2 = if_else(!(cell_type %in% c("L5ET2", "MicroPVM1")), "IN", cell_type), 
         cell_type2 = factor(cell_type2, levels = c("L5ET2", "IN", "MicroPVM1")))
ggplot(chimp_TWIST1_data.sub, aes(x=cell_type2, y=expr, fill=cell_type2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.25) +
  scale_fill_manual(values = c(cluster_colors[names(cluster_colors) == "L5ET2"],
                               "IN" = IN_color,
                               cluster_colors[names(cluster_colors) == "MicroPVM1"])) +
  theme_classic(base_size=20) +
  labs(x="", y=paste0(g, " expression")) +
  scale_x_discrete(labels = c("L5ET2" = "L5\nET_2", "IN" = "IN", "MicroPVM1" = "Micro-\nPVM_1")) +
  theme(legend.position = "none", axis.text.y = element_text(size=20),
        text = element_text(family = font),
        axis.text.x = element_text(size=20))
ggsave("TWIST1_L5ET2_IN_MicroPVM1.pdf")

# Fig. S24C
ggplot(chimp_TWIST1_data, aes(x=cell_type, y=expr, fill=cell_type)) +
  geom_boxplot() +
  theme_classic(base_size=22) +
  labs(x="", y=paste0(g, " expression")) +
  scale_fill_manual(values=cluster_colors) +
  scale_x_discrete(labels=cluster_names) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=13), legend.position = "none",
        text = element_text(family = font))
ggsave("TWIST1_chimp.pdf")

# plot Song et al., Nature 2020 PLACseq data ----
PTPRG_color="#cc9900"
HAR_color="#cc276e"

pp <- getDefaultPlotParams(plot.type=1)
pp$ideogramheight <- 0
pp$bottommargin <- 50
pp$data1inmargin <- 5

# get genes
PTPRG.region <- toGRanges("chr3:61,561,569-62,297,609")
PTPRG.kp <- plotKaryotype(zoom = PTPRG.region)
PTPRG.gene <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg38.knownGene,
                                    karyoplot=PTPRG.kp,
                                    plot.transcripts = TRUE, 
                                    plot.transcripts.structure = TRUE)
PTPRG.gene <- addGeneNames(PTPRG.gene)
PTPRG.gene <- mergeTranscripts(PTPRG.gene)

FHIT.region <- toGRanges("chr3:59,747,277-61,251,452")
FHIT.kp <- plotKaryotype(zoom = FHIT.region)
FHIT.gene <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg38.knownGene,
                                    karyoplot=FHIT.kp,
                                    plot.transcripts = TRUE, 
                                    plot.transcripts.structure = TRUE)
FHIT.gene <- addGeneNames(FHIT.gene)
FHIT.gene <- mergeTranscripts(FHIT.gene)

pdf("region_FHIT-to-PTPRG.pdf", width=8, height=3)
plot.region <- toGRanges("chr3:61210000-61650000")
kp <- plotKaryotype(zoom = plot.region, plot.params=pp, cex=1.4)
kpAddBaseNumbers(kp, tick.dist=100000, minor.tick.dist=50000, cex=1.2)
kpPlotGenes(kp, data=PTPRG.gene, r0=0, r1=0.35, add.gene.names=FALSE, gene.name.position = "left", gene.name.cex = 1.25, col=PTPRG_color)
kpText(kp, chr="chr3", x=61610000, y=0.16, labels="PTPRG", cex=1.4, col=PTPRG_color)
kpPlotGenes(kp, data=FHIT.gene, r0=0, r1=0.35, add.gene.names=FALSE, gene.name.position = "right", gene.name.cex = 1)
kpText(kp, chr="chr3", x=61230000, y=0.16, labels="FHIT", cex=1.2, col="black")

# plot HARsv2_1818
kpPlotRegions(kp, r0=0, r1=0.14, data=data.frame(chr=c("chr3"), start=c(61283266), end=c(61283416)), col=HAR_color)
kpText(kp, chr="chr3", x=61340000, y=0.07, labels="HARsv2_1818", cex=1.4, col=HAR_color)

# plot enhancer-promoter contact (from Song2020)
start.contact <- toGRanges(data.frame("chr3", 61280000, 61285000))
end.contact <- toGRanges(data.frame("chr3", 61560000, 61565000))
kpPlotLinks(kp, data=start.contact, data2=end.contact, col="black", r0=0.17, r1=0.5)

dev.off()

# plot TWIST1 motif ----
# MA1123.1 = TWIST1 (pfm from JASPAR, downloaded on 220223)
TWIST1 <- read_jaspar(here("data/HARs_hCONDELs_HAQERs", "MA1123.1.jaspar"))
TWIST1 <- convert_motifs(TWIST1, class = "motifStack-pfm")
TWIST1.RC <- motif_rc(TWIST1)
p1 <- ggseqlogo(TWIST1.RC@mat) + theme(axis.text.x = element_blank()) + 
  theme(plot.margin = unit(c(0, 1, 0, 1), "cm"))
aln = data.frame(letter=strsplit("GTACAACTGTAGTGTACATCTGTAGTGTACATCTGTAGTGTACATCTGTGGTGTACATCAGCAGT", "")[[1]],
                 species = rep(c("H", "C", "G", "R", "M"), each = 13), 
                 x = rep(1:13, 5))
aln$species <- factor(aln$species, levels=rev(c("H", "C", "G", "R", "M")))
aln$mut <- 0
aln$mut[6] <- 1
aln$mut[c(19, 32, 45, 58)] <- 2
aln$mut <- factor(aln$mut, levels=c(0,1,2))
p2 <- ggplot(aln, aes(x, species)) +
  geom_text(aes(label=letter, color=mut)) + 
  scale_x_continuous(breaks=1:10, expand = c(0.105, 0)) + 
  xlab('') + 
  scale_color_manual(values=c('black', "#038542", "#d72639")) + 
  theme_logo() + 
  theme(legend.position = 'none', axis.text.x = element_blank(), 
        plot.margin = unit(c(0, 1, 0, 1), "cm"), text = element_text(family = font))
plot_grid(p1, p2, ncol=1, align="v")
ggsave("motif_TWIST1_withseq.pdf", width=4.5, height=2.25)
