library(tidyverse)
library(here)

# defaults ----
cluster_color_table <- read.table(here("data/HARs_hCONDELs_HAQERs", "colors_for_figure.txt"), sep="\t", header=FALSE)
cluster_colors <- cluster_color_table$V3
names(cluster_colors) <- cluster_color_table$V1
cluster_names <- cluster_color_table$V2
names(cluster_names) <- cluster_color_table$V1
font <- "Helvetica"
EN_color <- "#8B0000"
IN_color <- "#00008B"
NN_color <- "#006400"
ct_colors <- c("EN" = EN_color, "IN" = IN_color, "NN" = NN_color)
cluster_order_table <- read.table(here("data/HARs_hCONDELs_HAQERs", "cross_species_clusters_order.txt"), 
                            sep="\t", header=FALSE)
cluster_order_table <- merge(cluster_order_table, cluster_color_table, by.x="V2", by.y="V2")
cluster_order <- cluster_order_table$V1.y[order(cluster_order_table$V1.x)]
ct_info <- read.table(here("data/HARs_hCONDELs_HAQERs", "celltypeInfo.txt"), header=TRUE, sep="\t")
cell_types = read.table(here("data/HARs_hCONDELs_HAQERs", "cell_types.txt"), header = TRUE, sep = "\t")

# regulatory domain size of human DEGs (Fig. S16D) ----
regdom_by_ct = read.table(here("data/HARs_hCONDELs_HAQERs", "human_regdom_by_ct.txt"), 
                          sep = "\t", header = TRUE)
ggplot(regdom_by_ct %>%
         mutate(cell_type = factor(cell_type, levels=cluster_order)), 
       aes(x=cell_type, y=length, fill=cell_type)) + 
  geom_boxplot(size=0.2) + 
  scale_y_log10() + 
  theme_classic(base_size=22) + 
  scale_fill_manual(values = cluster_colors) + 
  scale_x_discrete(labels=cluster_names) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12), 
        legend.position = "none", text = element_text(family = font)) + 
  labs(x="", y="size of regulatory domain (bp)")
ggsave("boxplot_human_DEGs_regdom.pdf")

# proportion of HARs/hCONDELs/HAQERs near human DEGs ----
# Fig. 6A
human_DEGs_HARs_hCONDELs = read.table(here("data/HARs_hCONDELs_HAQERs", "human_DEGs_HARs_hCONDELs.txt"), 
                          sep="\t", header=TRUE)
human_DEGs_HARs_hCONDELs <- merge(human_DEGs_HARs_hCONDELs, ct_info, by="cell_type") 
ggplot(human_DEGs_HARs_hCONDELs %>%
         mutate(prop_HARs_hCONDELs = as.numeric(prop_HARs_hCONDELs), 
                sig = ifelse(padj_genes_near_HARs_hCONDELs < 0.05 & 
                                padj_HARs_hCONDELs < 0.05, "*", NA)) %>%
         arrange(desc(prop_HARs_hCONDELs)) %>%
         mutate(cell_type = factor(cell_type, levels=cell_type)), 
       aes(x=cell_type, y=prop_HARs_hCONDELs, fill=cat)) + 
  geom_col() + 
  geom_text(aes(label=sig), vjust=0.2, size=6) +
  scale_fill_manual(values=ct_colors) +
  labs(x="", y=paste0("Proportion of ", species_name, " DEGs\nnear HARs or hCONDELs"), fill="") +
  theme_classic(base_size=22) + 
  scale_y_continuous(expand = c(0, 0), limits=c(0,0.46)) +
  scale_x_discrete(labels=cluster_names) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12), 
        axis.title.y = element_text(size=20, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        text = element_text(family=font))
ggsave("barplot_human_DEGs_nearHARsorhCONDELs.pdf")

# Fig. S16A
ggplot(human_DEGs_HARs_hCONDELs %>%
         mutate(prop_HARs = as.numeric(prop_HARs), 
                sig = ifelse(padj_genes_near_HARs < 0.05 & 
                               padj_HARs < 0.05, "*", NA)) %>%
         arrange(desc(prop_HARs)) %>%
         mutate(cell_type = factor(cell_type, levels=cell_type)), 
       aes(x=cell_type, y=prop_HARs, fill=cell_type)) + 
  geom_col() + 
  geom_text(aes(label=sig), vjust=0.2, size=6) +
  scale_fill_manual(values=cluster_colors) +
  labs(x="", y=paste0("Proportion of ", species_name, "\nDEGs near HARs"), fill="") +
  theme_classic(base_size=22) + 
  scale_y_continuous(expand = c(0, 0), limits=c(0,0.46)) +
  scale_x_discrete(labels=cluster_names) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12), 
        axis.title.y = element_text(size=20, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        text = element_text(family=font), legend.position="none")
ggsave("barplot_human_DEGs_nearHARs.pdf")

# Fig. S16B
ggplot(human_DEGs_HARs_hCONDELs %>%
         mutate(prop_hCONDELs = as.numeric(prop_hCONDELs), 
                sig = ifelse(padj_genes_near_hCONDELs < 0.05 & 
                               padj_hCONDELs < 0.05, "*", NA)) %>%
         arrange(desc(prop_hCONDELs)) %>%
         mutate(cell_type = factor(cell_type, levels=cell_type)), 
       aes(x=cell_type, y=prop_hCONDELs, fill=cell_type)) + 
  geom_col() + 
  geom_text(aes(label=sig), vjust=0.2, size=6) +
  scale_fill_manual(values=cluster_colors) +
  labs(x="", y=paste0("Proportion of ", species_name, "\nDEGs near hCONDELs"), fill="") +
  theme_classic(base_size=22) + 
  scale_y_continuous(expand = c(0, 0), limits=c(0,0.23)) +
  scale_x_discrete(labels=cluster_names) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12), 
        axis.title.y = element_text(size=20, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        text = element_text(family=font), legend.position="none")
ggsave("barplot_human_DEGs_nearhCONDELs.pdf")

# Fig. S16C
human_DEGs_HAQERs = read.table(here("data/HARs_hCONDELs_HAQERs", "human_DEGs_HAQERs.txt"), 
            sep="\t", header = TRUE)
human_DEGs_HAQERs <- merge(human_DEGs_HAQERs, ct_info, by="cell_type")
ggplot(human_DEGs_HAQERs %>%
         mutate(prop_HAQERs = as.numeric(prop_HAQERs), 
                sig = ifelse(padj_genes_near_HAQERs < 0.05 & 
                               padj_HAQERs < 0.05, "*", NA)) %>%
         arrange(desc(prop_HAQERs)) %>%
         mutate(cell_type = factor(cell_type, levels=cell_type)),
       aes(x=cell_type, y=prop_HAQERs, fill=cell_type)) + 
  geom_col() + 
  geom_text(aes(label=sig), vjust=0.2, size=6) +
  scale_fill_manual(values=cluster_colors) +
  labs(x="", y=paste0("Proportion of ", species_name, "\nDEGs near HAQERs"), fill="") +
  theme_classic(base_size=22) + 
  scale_y_continuous(expand = c(0, 0), limits=c(0,0.23)) +
  scale_x_discrete(labels=cluster_names) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12), 
        axis.title.y = element_text(size=20, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        text = element_text(family=font), legend.position = "none")
ggsave("barplot_human_DEGs_nearHAQERs.pdf")

# proportion of HARs/hCONDELs/HAQERs near DEGs between humans and chimps (Fig. S17) ----
hvc_DEGs_HARs_hCONDELs_HAQERs = read.table(here("data/HARs_hCONDELs_HAQERs", 
                                                "human_hvc_DEGs_HARs_hCONDELs_HAQERs.txt"), 
                                           sep = "\t", header = TRUE)
hvc_DEGs_HARs_hCONDELs_HAQERs <- merge(hvc_DEGs_HARs_hCONDELs_HAQERs, cell_types, by="cell_type")

# HARs
ggplot(hvc_DEGs_HARs_hCONDELs_HAQERs %>%
         mutate(prop_HARs = as.numeric(prop_HARs), 
                sig = ifelse(padj_genes_near_HARs < 0.05 & 
                               padj_HARs < 0.05, "*", NA)) %>%
         arrange(desc(prop_HARs)) %>%
         mutate(cell_type = factor(cell_type, levels=cell_type)), 
       aes(x=cell_type, y=prop_HARs, fill=cat)) + 
  geom_col() + 
  geom_text(aes(label=sig), vjust=0.2, size=6) +
  scale_fill_manual(values=ct_colors) +
  labs(x="", y=paste0("Proportion of ", species_name, "\nDEGs near HARs"), fill="") +
  theme_classic(base_size=22) + 
  scale_y_continuous(expand = c(0, 0), limits=c(0,0.24)) +
  scale_x_discrete(labels=cluster_names) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12), 
        axis.title.y = element_text(size=20, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        text = element_text(family=font))
ggsave("barplot_hvc_DEGs_nearHARs.pdf")

# hCONDELs
ggplot(hvc_DEGs_HARs_hCONDELs_HAQERs %>%
         mutate(prop_hCONDELs = as.numeric(prop_hCONDELs), 
                sig = ifelse(padj_genes_near_hCONDELs < 0.05 & 
                               padj_hCONDELs < 0.05, "*", NA)) %>%
         arrange(desc(prop_hCONDELs)) %>%
         mutate(cell_type = factor(cell_type, levels=cell_type)), 
       aes(x=cell_type, y=prop_hCONDELs, fill=cat)) + 
  geom_col() + 
  geom_text(aes(label=sig), vjust=0.2, size=6) +
  scale_fill_manual(values=ct_colors) +
  labs(x="", y=paste0("Proportion of ", species_name, "\nDEGs near hCONDELs"), fill="") +
  theme_classic(base_size=22) + 
  scale_y_continuous(expand = c(0, 0), limits=c(0,0.14)) +
  scale_x_discrete(labels=cluster_names) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12), 
        axis.title.y = element_text(size=20, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        text = element_text(family=font))
ggsave("barplot_hvc_DEGs_nearhCONDELs.pdf")

# HAQERs
ggplot(hvc_DEGs_HARs_hCONDELs_HAQERs %>%
         mutate(prop_HAQERs = as.numeric(prop_HAQERs), 
                sig = ifelse(padj_genes_near_HAQERs < 0.05 & 
                               padj_HAQERs < 0.05, "*", NA)) %>%
         arrange(desc(prop_HAQERs)) %>%
         mutate(cell_type = factor(cell_type, levels=cell_type)), 
       aes(x=cell_type, y=prop_HAQERs, fill=cat)) + 
  geom_col() + 
  geom_text(aes(label=sig), vjust=0.2, size=6) +
  scale_fill_manual(values=ct_colors) +
  labs(x="", y=paste0("Proportion of ", species_name, "\nDEGs near HAQERs"), fill="") +
  theme_classic(base_size=22) + 
  scale_y_continuous(expand = c(0, 0), limits=c(0,0.09)) +
  scale_x_discrete(labels=cluster_names) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12), 
        axis.title.y = element_text(size=20, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        text = element_text(family=font))
ggsave("barplot_hvc_DEGs_nearHAQERs.pdf")

# HARs/hCONDELs near human, chimp, gorilla, and rhesus (Fig. S18) ----
merge_data = rbind(read.table(here("data/HARs_hCONDELs_HAQERs", "human_DEGs_HARs_hCONDELs.txt"), 
                              sep="\t", header=TRUE) %>% 
                     mutate(species="H"), 
                   read.table(here("data/HARs_hCONDELs_HAQERs", "chimp_DEGs_HARs_hCONDELs.txt"), 
                              sep="\t", header=TRUE) %>%
                     mutate(species="C"), 
                   read.table(here("data/HARs_hCONDELs_HAQERs", "gorilla_DEGs_HARs_hCONDELs.txt"), 
                              sep="\t", header=TRUE) %>%
                     mutate(species="G"), 
                   read.table(here("data/HARs_hCONDELs_HAQERs", "rhesus_DEGs_HARs_hCONDELs.txt"), 
                              sep="\t", header=TRUE) %>%
                     mutate(species="R")) %>% 
  select(species, cell_type, prop_HARs, prop_hCONDELs, prop_HARs_hCONDELs, 
         padj_genes_near_HARs, padj_genes_near_hCONDELs, 
         padj_genes_near_HARs_hCONDELs, padj_HARs, padj_hCONDELs, padj_HARs_hCONDELs) %>%
  mutate(cell_type = factor(cell_type, levels=cluster_order), 
         species = factor(species, levels=c("R", "G", "C", "H")))

num_sig <- data.frame(rbind(c("H", nrow(human_data[human_data$padj_genes_near_HARs_hCONDELs < 0.05 & 
                                                     human_data$padj_HARs_hCONDELs < 0.05, ]), 
                              nrow(human_data[human_data$padj_genes_near_HARs < 0.05 & 
                                                human_data$padj_HARs < 0.05, ]), 
                              nrow(human_data[human_data$padj_genes_near_hCONDELs < 0.05 & 
                                                human_data$padj_hCONDELs < 0.05, ])), 
                            c("C", nrow(chimp_data[chimp_data$padj_genes_near_HARs_hCONDELs < 0.05 & 
                                                     chimp_data$padj_HARs_hCONDELs < 0.05, ]), 
                              nrow(chimp_data[chimp_data$padj_genes_near_HARs < 0.05 & 
                                                chimp_data$padj_HARs < 0.05, ]), 
                              nrow(chimp_data[chimp_data$padj_genes_near_hCONDELs < 0.05 & 
                                                chimp_data$padj_hCONDELs < 0.05, ])),
                            c("G", nrow(gorilla_data[gorilla_data$padj_genes_near_HARs_hCONDELs < 0.05 & 
                                                       gorilla_data$padj_HARs_hCONDELs < 0.05, ]), 
                              nrow(gorilla_data[gorilla_data$padj_genes_near_HARs < 0.05 & 
                                                  gorilla_data$padj_HARs < 0.05, ]), 
                              nrow(gorilla_data[gorilla_data$padj_genes_near_hCONDELs < 0.05 & 
                                                  gorilla_data$padj_hCONDELs < 0.05, ])),
                            c("R", nrow(rhesus_data[rhesus_data$padj_genes_near_HARs_hCONDELs < 0.05 & 
                                                      rhesus_data$padj_HARs_hCONDELs < 0.05, ]), 
                              nrow(rhesus_data[rhesus_data$padj_genes_near_HARs < 0.05 & 
                                                 rhesus_data$padj_HARs < 0.05, ]), 
                              nrow(rhesus_data[rhesus_data$padj_genes_near_hCONDELs < 0.05 & 
                                                 rhesus_data$padj_hCONDELs < 0.05, ]))))
colnames(num_sig) <- c("species", "HARs_hCONDELs", "HARs", "hCONDELs")
num_sig$HARs_hCONDELs <- as.numeric(num_sig$HARs_hCONDELs)
num_sig$HARs <- as.numeric(num_sig$HARs)
num_sig$hCONDELs <- as.numeric(num_sig$hCONDELs)
num_sig$species <- factor(num_sig$species, levels = c("R", "G", "C", "H"))

# plot max p-value of two tests
# HARs and hCONDELs
merge_data$padj <- pmax(merge_data$padj_genes_near_HARs_hCONDELs, merge_data$padj_HARs_hCONDELs)
plot_HARs_and_hCONDELs <- 
  ggplot(merge_data %>%
           mutate(padj=-1*log10(padj)), 
         aes(x=cell_type, y=species, color=padj, size=prop_HARs_hCONDELs)) + 
  geom_point() + 
  scale_color_gradientn(colors = c("gray", "gray", "pink", "darkred"),
                        values = scales::rescale(c(0, 1.3, 1.31, 2, 16)), 
                        limits=c(0,16.1)) +
  scale_size_continuous(range = c(1, 5), limits=c(0,0.43)) +
  theme_bw(base_size=16) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12), 
        text = element_text(family = font)) +
  scale_x_discrete(labels=cluster_names) +
  labs(y="", x="", color="-log10(padj)", size="proportion")
ggsave(plot_HARs_and_hCONDELs, file="dotplot_HARs_hCONDELs.pdf")

ggplot(num_sig, aes(x=species, y=HARs_hCONDELs, fill=species)) + 
  geom_col() + 
  theme_pubr(base_size=16) + 
  scale_y_continuous(expand=c(0,0), position="right") + 
  labs(x="", y="Number of cell types enriched\nfor genes near HARs/hCONDELs") + 
  theme(legend.position = "none") + 
  scale_fill_manual(values=species_colors) + 
  coord_flip()
ggsave("barplot_num_sig_cell_types_HARs_hCONDELs.pdf")

prop = merge(merge_data %>% filter(species != "H") %>% 
               mutate(sig_other = if_else(padj < 0.05, 1, 0)) %>% 
               select(species, cell_type, prop_HARs_hCONDELs, sig_other) %>% 
               rename(prop = prop_HARs_hCONDELs, comp = species), 
             merge_data %>% mutate(sig = if_else(padj < 0.05, 1, 0)) %>% 
               filter(species == "H") %>% 
               select(cell_type, prop_HARs_hCONDELs, sig) %>% 
               rename(H = prop_HARs_hCONDELs), 
             by="cell_type") %>% 
  mutate(sig_join = if_else(sig == 0 & sig_other == 0, "not sig", 
                            if_else(sig == 1 & sig_other == 1, "both", 
                                    if_else(sig == 1, "H", "other")))) %>%
  select(cell_type, comp, H, prop, sig_join)
prop$comp <- factor(prop$comp, levels = c("C", "G", "R"))
ggplot(prop, aes(x=prop, y=H, color=sig_join)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dotted") +
  theme_pubr(base_size=16) +
  scale_x_continuous(expand=c(0,0), limits=c(0.09, 0.45)) + 
  scale_y_continuous(expand=c(0,0), limits = c(0.09, 0.45)) + 
  scale_color_manual(values = c("both" = "#ff0000", "H" = human_color, 
                                "C" = chimp_color, "G" = gorilla_color, 
                                "R" = rhesus_color, "not sig" = "#bbbbbb", 
                                "other" = "#000000")) +
  labs(x = "Proportion of DEGs near HARs/hCONDELs", 
       y = "Proportion of human DEGs\nnear HARs/hCONDELs", 
       color = "adj. pval < 0.05") +
  facet_wrap(~comp)
ggsave("scatterplot_HARs_hCONDELs.pdf")

# HARs
merge_data$padj <- pmax(merge_data$padj_genes_near_HARs, merge_data$padj_HARs)
plot_HARs <- 
  ggplot(merge_data %>%
           mutate(padj=-1*log10(padj)), 
         aes(x=cell_type, y=species, color=padj, size=prop_HARs)) + 
  geom_point() + 
  scale_color_gradientn(colors = c("gray", "gray", "pink", "darkred"),
                        values = scales::rescale(c(0, 1.3, 1.31, 2, 16)), 
                        limits=c(0,16.1)) +
  scale_size_continuous(range = c(1, 5), limits=c(0,0.43)) +
  theme_bw(base_size=16) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12), 
        text = element_text(family = font)) +
  scale_x_discrete(labels=cluster_names) +
  labs(y="", x="", color="-log10(padj)", size="proportion")
ggsave(plot_HARs, file="dotplot_HARs.pdf")

ggplot(num_sig, aes(x=species, y=HARs, fill=species)) + 
  geom_col() + 
  theme_pubr(base_size=16) + 
  scale_y_continuous(expand=c(0,0), position="right") + 
  labs(x="", y="Number of cell types enriched\nfor genes near HARs") + 
  theme(legend.position = "none") + 
  scale_fill_manual(values=species_colors) + 
  coord_flip()
ggsave("barplot_num_sig_cell_types_HARs.pdf")

# hCONDELs 
merge_data$padj <- pmax(merge_data$padj_genes_near_hCONDELs, merge_data$padj_hCONDELs)
plot_hCONDELs <- 
  ggplot(merge_data %>%
           mutate(padj=-1*log10(padj)), 
         aes(x=cell_type, y=species, color=padj, size=prop_hCONDELs)) + 
  geom_point() + 
  scale_color_gradientn(colors = c("gray", "gray", "pink", "darkred"),
                        values = scales::rescale(c(0, 1.3, 1.31, 2, 16)), 
                        limits=c(0,16.1)) +
  scale_size_continuous(range = c(1, 5), limits=c(0,0.43)) +
  theme_bw(base_size=16) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12), 
        text = element_text(family = font)) +
  scale_x_discrete(labels=cluster_names) +
  labs(y="", x="", color="-log10(padj)", size="proportion")
ggsave(plot_hCONDELs, file="dotplot_hCONDELs.pdf")

ggplot(num_sig, aes(x=species, y=hCONDELs, fill=species)) + 
  geom_col() + 
  theme_pubr(base_size=16) + 
  scale_y_continuous(expand=c(0,0), position="right") + 
  labs(x="", y="Number of cell types enriched\nfor genes near hCONDELs") + 
  theme(legend.position = "none") + 
  scale_fill_manual(values=species_colors) + 
  coord_flip()
ggsave("barplot_num_sig_cell_types_hCONDELs.pdf")

# Enrichment of hDEGs (near HARs/hCONDELs) in SynGO (Fig. 6B) ----
SynGO_data <- read.table(here("data/HARs_hCONDELs_HAQERs", "SynGO_hDEGs_HARs_hCONDELs.txt"), 
                         header=TRUE, sep="\t") %>%
  select(gene_group, OR_hDEG_enrich, padj_hDEG_enrich, OR_hDEG_near_HAR_hCONDEL_enrich, 
         padj_hDEG_near_HAR_hCONDEL_enrich)

cats_to_plot <- c("SynGO", "trans-synaptic signaling", "synapse assembly", 
                  "presynaptic membrane", "postsynaptic membrane", "metabolism", 
                  "transport")
names(cats_to_plot) <- c("SynGO", "SynGO-L3-BP-trans-synaptic signaling", 
                         "SynGO-L3-BP-synapse assembly", 
                         "SynGO-L3-CC-presynaptic membrane", 
                         "SynGO-L3-CC-postsynaptic membrane", 
                         "SynGO-L2-BP-metabolism", "SynGO-L2-BP-transport")
SynGO_data.sub <- SynGO_data[SynGO_data$gene_group %in% names(cats_to_plot), ]
SynGO_data.sub <- rbind(SynGO_data.sub %>% 
                          select(gene_group, OR_hDEG_enrich, padj_hDEG_enrich) %>%
                          rename(OR = OR_hDEG_enrich, 
                                 padj = padj_hDEG_enrich) %>%
                          mutate(stat = "hDEG_enrich"), 
                        SynGO_data.sub %>% 
                          select(gene_group, OR_hDEG_near_HAR_hCONDEL_enrich, 
                                 padj_hDEG_near_HAR_hCONDEL_enrich) %>%
                          rename(OR = OR_hDEG_near_HAR_hCONDEL_enrich, 
                                 padj = padj_hDEG_near_HAR_hCONDEL_enrich) %>%
                          mutate(stat = "hDEG_near_HAR_hCONDEL_enrich")) %>%
  mutate(padj = -1*log10(padj), 
         gene_group = factor(gene_group, levels=rev(names(cats_to_plot))))

ggplot(SynGO_data.sub, aes(x=stat, y=gene_group, color=padj, size=OR)) + 
  geom_point() + 
  scale_color_gradientn(colors = c("gray", "gray", "pink", "darkred"), 
                        values = scales::rescale(c(0, 1.3, 1.31, 1.32, 5, 16))) +
  scale_size_continuous(range = c(1, 9)) +
  theme_bw(base_size=16) + 
  scale_y_discrete(labels=cats_to_plot) + 
  scale_x_discrete(labels=c("hDEG_enrich" = "hDEGs", 
                            "hDEG_near_HAR_hCONDEL_enrich" = "hDEGs near\nHAR/hCONDEL")) +
  labs(y="", x="", color="-log10(padj)")
ggsave("dotplot_SynGO_enrichments.pdf")
