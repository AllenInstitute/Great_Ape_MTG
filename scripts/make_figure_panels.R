library(tidyverse)
library(magrittr)
library(RColorBrewer)
library(rstatix)
library(ggbeeswarm)
library(GGally)
library(ggVennDiagram)


theme_set(theme_bw())



#### Figure 2 ####
to_plot <- readRDS("data/downsampled_subclass_spearman_corr.RDS")


# Figure 2D. Subclass expression correlations
to_plot %>%
  filter(species_1 != species_2) %>% 
  ggplot() +
  geom_tile(aes(x = species_1, y = species_2 %>% fct_rev(), fill = spearman_corr)) +
  scale_fill_viridis_c(option = "B", na.value = "gray") +
  labs(x = "", y = "", fill = "Spearman corr.") +
  facet_wrap(~celltype, ncol = 9) +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2))


# Add evo dist
evo_distm <- matrix(c(0, 6, 7, 27, 38,
                      6, 0, 7, 27, 38,
                      7, 7, 0, 27, 38,
                      27, 27, 27, 0, 38,
                      38, 38, 38, 38, 0), 
                    nrow = 5, ncol = 5, byrow = TRUE,
                    dimnames = list(c("human", "chimp", "gorilla", "rhesus", "marmoset"), 
                                    c("human", "chimp", "gorilla", "rhesus", "marmoset")))

to_plot$evo_dist <- NA
for (i in 1:nrow(to_plot)) {
  to_plot$evo_dist[i] <- evo_distm[to_plot$species_1[i], to_plot$species_2[i]]
}


to_plot$species_1_2 <- paste0(to_plot$species_1, "_", to_plot$species_2)

to_plot$cellsuperclass <- ifelse(to_plot$celltype %in% c("OPC", "Astro", "Oligo", "Endo", "Micro-PVM", "VLMC"), 
                                 "Non-neuronal", "Neuronal")

to_plot$cellclass <- NA
to_plot$cellclass[to_plot$celltype %in% c("OPC", "Astro", "Oligo")] <- "Glial"
to_plot$cellclass[to_plot$celltype %in% c("Endo", "Micro-PVM", "VLMC")] <- "Non-neural"
to_plot$cellclass[to_plot$celltype %in% c("L2/3 IT", "L4 IT", "L5 IT", "L6 IT", "L6 IT Car3",
                                          "L5 ET", "L5/6 NP", "L6 CT", "L6b")] <- "Excitatory"
to_plot$cellclass[to_plot$celltype %in% c("Lamp5_Lhx6", "Lamp5", "Sncg", "Vip", "Pax6", 
                                          "Chandelier", "Pvalb", "Sst", "Sst Chodl")] <- "Inhibitory"

to_plot_nodupes <- to_plot[! duplicated(to_plot$spearman_corr), ]
to_plot_nodupes <- na.omit(to_plot_nodupes)


# Remove low sample size correlations
to_plot <- droplevels(subset(to_plot, ! celltype %in% c("Sst Chodl", "Endo", "VLMC")))


# Plot human vs. other sp.
type_colors <- setNames(c(rev(grey.colors(nlevels(to_plot$celltype) - 4)), 
                          brewer.pal(4, "Set2")), 
                        levels(to_plot$celltype))


# Figure 2E. Species evo trend
p1 <- to_plot %>% 
  filter(species_1 == "human") %>%
  # filter(cellsuperclass == "Neuronal") %>% 
  ggplot(aes(x = evo_dist, y = spearman_corr, color = celltype)) +
  geom_line(size = 2) +
  scale_color_manual(values = type_colors) +
  scale_x_continuous(breaks = c(0, 6, 7, 27, 38)) +
  xlab("Evo. distance (mya)") +
  ylab("Spearman correlation") +
  ggtitle(species1) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())
plot(p1)

ggsave(p1, filename = "output/downsampled_expr_cor_evo_trend.pdf", height = 4, width = 7)


# Figure 2F. Human-accelerated evo vs. great apes
greatape_cor_set <- c("human_chimp", "human_gorilla", "chimp_gorilla")

p1 <- to_plot %>% 
  filter(species_1_2 %in% greatape_cor_set & cellsuperclass == "Neuronal") %>%
  ggplot(aes(x = fct_relevel(species_1_2, greatape_cor_set), y = spearman_corr)) +
  geom_boxplot(width = 0.5) +
  geom_jitter(width = 0.1, height = 0) +
  facet_wrap(~cellsuperclass, scales = "free_y", nrow = 2) +
  xlab("") +
  ylab("Spearman correlation") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
plot(p1)

ggsave(p1, filename = "output/subclass_cor_boxplots_greatapes.pdf", width = 2.5, height = 4)



#### Figure S4 ####

species_cor_set <- c(
  "human_chimp", 
  "human_gorilla", "chimp_gorilla",
  "human_rhesus", "chimp_rhesus", "gorilla_rhesus", 
  "human_marmoset", "chimp_marmoset", "gorilla_marmoset", "rhesus_marmoset")

# Figure S4D. All species cor
p1 <- to_plot %>% 
  filter(species_1_2 %in% species_cor_set) %>%
  ggplot(aes(x = fct_relevel(species_1_2, species_cor_set), y = spearman_corr)) +
  geom_boxplot(width = 0.5) +
  geom_jitter(width = 0.1, height = 0) +
  facet_wrap(~cellsuperclass, scales = "free_y", nrow = 2) +
  xlab("") +
  ylab("Spearman correlation") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
plot(p1)

ggsave(p1, filename = "output/subclass_cor_boxplots.pdf", width = 5, height = 10)



# Compare between to within-species variation
to_plot_within <- to_plot %>% 
  filter(species_1 == species_2)

within_cor_mean <- NULL
for(i in 1:nrow(to_plot)) {
  cor1 <- to_plot_within %>% 
    filter(species_1 == to_plot$species_1[i] & celltype == to_plot$celltype[i]) %>% 
    pull(spearman_corr)
  
  cor2 <- to_plot_within %>% 
    filter(species_1 == to_plot$species_2[i] & celltype == to_plot$celltype[i]) %>% 
    pull(spearman_corr)
  
  cor_mean <- mean(c(cor1, cor2))
  
  within_cor_mean <- c(within_cor_mean, cor_mean)
}

to_plot$within_cor_mean <- within_cor_mean


# Figure S4C. Within vs. between species variation - SI
p1 <- to_plot %>% 
  filter(species_1_2 %in% species_cor_set) %>%
  # filter(cellsuperclass == "Non-neuronal") %>%
  arrange(celltype, -evo_dist) %>% 
  ggplot(aes(x = within_cor_mean, y = spearman_corr, 
             color = celltype, size = evo_dist)) +
  # geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  scale_size_continuous(range = c(1, 5)) +
  scale_color_manual(values = type_colors) +
  guides(color = "none") +
  xlab("Average within-species similarity (r)") +
  ylab("Between-species similarity (r)")
plot(p1)

ggsave(p1, file = "output/between_vs_within_cor_scatter.pdf", height = 3, width = 4)


# Figure S4C. Ratio of between:within species variation - SI
p1 <- to_plot %>% 
  filter(species_1_2 %in% species_cor_set) %>%
  mutate(celltype2 = ifelse(cellsuperclass == "Non-neuronal", 
                            as.character(celltype), cellsuperclass), 
         cor_ratio = spearman_corr/within_cor_mean) %>% 
  ggplot(aes(x = fct_reorder(celltype2, cor_ratio, mean), y = cor_ratio,  
             color = celltype)) +
  # geom_jitter(width = 0.2) +
  geom_quasirandom() +
  # scale_size_continuous(range = c(1, 5)) +
  scale_color_manual(values = type_colors) +
  guides(color = "none") +
  labs(x = "", 
       y = "Between:within-species correlation") +
  guides(x = guide_axis(angle = 45))
plot(p1)

ggsave(p1, file = "output/between_within_cor_ratio.pdf", height = 4, width = 3)


# Stats
to_plot %>% 
  filter(species_1_2 %in% species_cor_set) %>%
  mutate(celltype2 = ifelse(cellsuperclass == "Non-neuronal", 
                            as.character(celltype), 
                            as.character(cellsuperclass)), 
         cor_ratio = spearman_corr/within_cor_mean) %>% 
  t_test(cor_ratio ~ celltype2, ref.group = "Neuronal",
         var.equal = TRUE, alternative = "greater") %>% 
  add_significance()



# Gene-gene correlations
species <- c("human", "chimp", "gorilla", "rhesus", "marmoset")
expr <- list()

for (species1 in species) {
  f_in <- paste0("data/subclass_log2cpm_means_orthologous_", 
                 species1, "_symbols.RDS")
  expr1 <- readRDS(f_in)
  expr1$species <- species1
  
  expr[[species1]] <- as_tibble(expr1)
}

expr_all <- bind_rows(expr)


# Corrs
cor_fun <- function(df, cor.method = "pearson") {
  cor_all <- cor(df[, 2:6])[1, -1]
  return(cor_all)
}

max_fun <- function(df) {
  max_expr <- max(df$human)
  return(max_expr)
}

expr_nest <- expr_all %>% 
  pivot_longer(cols = Astro:VLMC) %>% 
  pivot_wider(names_from = species, values_from = value) %>% 
  group_by(genes) %>% 
  nest()

expr_cor <- expr_nest %>%
  mutate(h_c = map_dbl(data, ~ cor(.$human, .$chimp))) %>% 
  mutate(h_g = map_dbl(data, ~ cor(.$human, .$gorilla))) %>% 
  mutate(h_r = map_dbl(data, ~ cor(.$human, .$rhesus))) %>% 
  mutate(h_m = map_dbl(data, ~ cor(.$human, .$marmoset))) %>% 
  mutate(r_m = map_dbl(data, ~ cor(.$rhesus, .$marmoset))) %>% 
  mutate(h_cl = map_chr(data, ~ .$name[which.max(.$human)])) %>% 
  mutate(h_max = map_dbl(data, ~ max(.$human))) %>% 
  select(-data) 

expr_corl <- expr_cor %>%
  pivot_longer(cols = h_c:h_m) %>% 
  mutate(name = fct_relevel(name, c("h_c", "h_g", "h_r", "h_m"))) %>% 
  group_by(genes) %>% 
  mutate(cor_mean = mean(value, na.rm = TRUE))


# Figure S4H. Correlation boxplots
expr_corl %>% 
  filter(h_max > 2) %>%
  ggplot(aes(x = name, y = value)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.25, color = "blue") +
  xlab("Species comparison") +
  ylab("Pearson correlation")


# Figure S4I. Div gene Venn
expr_div <- expr_corl %>% 
  filter(h_max > 2) %>%
  filter(value < 0.25)

div_genes_list <- list()

for (pw1 in c("h_c", "h_g", "h_r", "h_m")) {
  div_genes1 <- expr_div %>% 
    filter(name == pw1) %>% 
    pull(genes)
  
  div_genes_list[[pw1]] <- div_genes1
}

ggVennDiagram(div_genes_list, 
              category.names = c("Chimp", "Gorilla", "Rhesus", "Marmoset"),  
              label_alpha = 0, edge_size = 0) +
  scale_fill_gradientn(colors = c("white", rep("light blue", 4))) +
  guides(fill = "none")


# Figure S4J. Plot divergent genes
to_plot <- expr_all %>% 
  pivot_longer(cols = Astro:VLMC) %>% 
  mutate(species = fct_relevel(species, c("human", "chimp", "gorilla", "rhesus", "marmoset")), 
         name = fct_relevel(name, c("OPC", "Astro", "Oligo", 
                                    "VLMC", "Endo", "Micro-PVM", 
                                    "Lamp5_Lhx6", "Lamp5", "Sncg", "Pax6", "Vip",  
                                    "Chandelier", "Pvalb", "Sst", "Sst Chodl", 
                                    "L2/3 IT", "L4 IT", "L5 IT", "L6 IT", "L6 IT Car3",
                                    "L5 ET", "L5/6 NP", "L6 CT", "L6b")))

to_plot$cellclass <- NA
to_plot$cellclass[to_plot$name %in% c("Astro", "OPC", "Oligo", 
                                      "Micro-PVM", "Endo", "VLMC")] <- "Non-neuronal"
to_plot$cellclass[to_plot$name %in% c("L2/3 IT", "L4 IT", "L5 IT", "L6 IT", "L6 IT Car3",
                                      "L5 ET", "L5/6 NP", "L6 CT", "L6b")] <- "Excitatory"
to_plot$cellclass[to_plot$name %in% c("Lamp5_Lhx6", "Lamp5", "Sncg", "Vip", "Pax6", 
                                      "Chandelier", "Pvalb", "Sst", "Sst Chodl")] <- "Inhibitory"

div_genes <- c("FAM177B", "MEPE", "PRLR")  # NPNT  # Representative human/ape-specific genes

# Dot plot
to_plot %>% 
  filter(genes %in% div_genes) %>% 
  ggplot(aes(x = species, y = fct_rev(name), size = value, color = cellclass)) +
  geom_point() +
  facet_wrap(~ genes, nrow = 1) +
  scale_size_area() +
  xlab("Species") +
  ylab("Cell subclass") +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))




#### Figure S14 ####

# Sst proportions
species <- c("Human", "Chimp", "Gorilla", "Rhesus", "Marmoset")

meta_sp <- list()
for (species1 in species) {
  meta_fn <-
    paste0(
      "data/", species1, "_Master_metadata_for_plots_and_sharing_12_16_21.RDS"
    )
  meta1 <- readRDS(meta_fn)
  
  if (! "layer" %in% colnames(meta1)) {
    meta1 <- cbind(meta1[, 1:10], 
                   layer = rep(NA, nrow(meta1)), 
                   meta1[, 11:18])
  }
  
  meta_sp[[species1]] <- meta1
}

meta <- do.call("rbind", meta_sp) %>% 
  mutate(species = as_factor(species)) %>% 
  mutate(species = fct_relevel(species, c("human", "chimp", "gorilla", "rhesus", "marmoset")))



# Figure S14K. Props
meta %>% 
  filter(subclass == "Sst" & tech == "10Xv3") %>% 
  add_count(species) %>% 
  add_count(species, cross_species_cluster) %>% 
  mutate(prop = nn / n) %>% 
  select(species, cross_species_cluster, prop) %>% 
  unique() %>% 
  mutate(cl_target = if_else(cross_species_cluster == "Sst_1", "yes", "no")) %>% 
  ggplot(aes(x = cross_species_cluster, y = prop, shape = species, color = cl_target)) +
  geom_jitter(width = 0.1) +
  scale_y_log10() +
  scale_color_manual(values = c("black", "red")) +
  scale_shape_manual(values = c(15, 16, 17, 3, 4)) +
  guides(color = FALSE) +
  ylab("Proportion") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank())



# Figure S14O. Dopamine receptor expression dot plot
# First run scripts for Figure S4 above
div_genes <- c("DRD1", "DRD2", "DRD3")  # SST TH receptors

to_plot %>% 
  filter(genes %in% div_genes) %>% 
  ggplot(aes(x = species, y = fct_rev(name), size = value, color = cellclass)) +
  geom_point() +
  facet_wrap(~ genes, nrow = 1) +
  scale_size_area() +
  xlab("Species") +
  ylab("Cell subclass") +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


# Figure S14N. Monkey cor
# Correlations
drd_cor <- to_plot %>% 
  filter(genes %in% div_genes) %>% 
  pivot_wider(names_from = species, values_from = value) %>% 
  select(-name, -cellclass) %>% 
  group_by(genes) %>% 
  nest() %>% 
  mutate(pw_cor = map(data, ~ cor(., method = "spearman")))

apes <- species[1:3]

drd_cor_subset <- tibble(
  species1 = rep(factor(apes, levels = apes), 2), 
  species2 = factor(c(rep("rhesus", 3), rep("marmoset", 3)), 
                    levels = c("rhesus", "marmoset")), 
  DRD1 = as.vector(drd_cor$pw_cor[[1]][apes, c("rhesus", "marmoset")]),
  DRD2 = as.vector(drd_cor$pw_cor[[2]][apes, c("rhesus", "marmoset")]),
  DRD3 = as.vector(drd_cor$pw_cor[[3]][apes, c("rhesus", "marmoset")])
) %>%
  pivot_longer(cols = c("DRD1", "DRD2", "DRD3")) %>% 
  group_by(name, species1) %>% 
  mutate(mean_cor = mean(value))

drd_cor_subset %>% 
  ggplot(aes(x = species1, y = value, shape = species2)) +
  geom_col(aes(y = mean_cor/2)) +
  geom_point() +
  facet_wrap(~ name) +
  scale_shape_manual(values = c(3, 4)) +
  xlab("") +
  ylab("Correlation") +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))




#### Figure 5, S15 ####
deg <- readRDS(file = "data/aggregated_DEG_lists.RDS")

deg <- deg %>%
  rename(gene = genes)

# Save csv
if (! file.exists("output/species_degs.csv")) {
  deg %>% 
    select(species, cluster, gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>% 
    filter(! is.na(padj) & padj < 0.05) %>% 
    write_csv(file = "output/species_degs.csv")
}

species_order <- c("human", "chimp", "gorilla", "rhesus", "marmoset")
lineage <- data.frame(species = species_order, mya = c(6, 6, 9, 25, 40 + 15))

deg_sig <- deg %>%
  left_join(lineage, by = "species") %>% 
  mutate(species = factor(species, levels = species_order), 
         sig_diff = if_else(padj <  0.01 & abs(log2FoldChange) > 0.5, sign(log2FoldChange), 0), 
         sig_diff = replace(sig_diff, log2FoldChange == 0, NA)) %>%
  filter(sig_diff %in% c(-1, 1)) %>%
  group_by(species, gene) %>%
  mutate(deg_types = n_distinct(cluster)) %>%
  group_by(species, deg_types) %>%
  mutate(deg_bin = n_distinct(gene)) %>%
  group_by(species) %>%
  mutate(deg_bin_prop = deg_bin / n_distinct(gene)) %>%
  group_by(gene) %>%
  mutate(num_lineage = n_distinct(species)) %>%
  ungroup()


# Add cluster meta
gene_cnt <- read_csv(file = "data/consensus_type_gene_cnt.csv")

gene_cnt_summary <- gene_cnt %>% 
  filter(species == "human") %>% 
  group_by(cross_species_cluster) %>% 
  summarize(median_genes = median(nFeature_RNA)) 

cluster_meta <- read_csv(file = "data/cross_species_cluster_counts.csv", 
                         show_col_types = FALSE) %>% 
  left_join(gene_cnt_summary, by = c("cluster_full_name" = "cross_species_cluster"))

deg_sig <- deg_sig %>% 
  left_join(cluster_meta, by = "cluster") 


# Add cluster expression count
gene_spec <- read_csv(file = "data/gene_spec_tau.csv")

deg_sig <- deg_sig %>% 
  left_join(gene_spec, by = c("species", "gene")) %>% 
  mutate(species = factor(species, levels = species_order))


# Figure S15B. DEG count vs. log2FC
p1 <- deg_sig %>%
  filter(species == "human") %>%
  group_by(species, cluster, sig_diff) %>% 
  mutate(log2fc_median = median(log2FoldChange), 
         deg_count = n_distinct(gene), 
         deg_prop = deg_count/median_genes) %>% 
  select(species, cluster, class, ds_count, full_count, sig_diff, log2fc_median, deg_count, deg_prop) %>% 
  unique() %>%
  ggplot(aes(x = deg_prop, y = abs(log2fc_median), 
             color = class, size = deg_count)) +
  geom_point(alpha = 0.5) +
  labs(x = "DEG proportion", y = "Abs. log2(fold-change)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot(p1)

ggsave(p1, filename = "output/deg_prop_vs_fc.pdf", width = 4, height = 3)


# Figure S15C. Many lineage-specific changes only occur in one or a few cell types
p1 <- deg_sig %>% 
  select(species, sig_diff, deg_types, deg_bin, deg_bin_prop) %>% 
  unique() %>% 
  ggplot(aes(x = deg_types, y = deg_bin, color = species)) +
  # geom_histogram(bins = 57) +
  geom_point() +
  geom_line() +
  xlab("No. of consensus types") +
  ylab("No. of DEGs") +
  ggtitle("DEG cell type specificity")
plot(p1)

ggsave(p1, filename = "output/deg_specificity.pdf", width = 4, height = 3)


# Figure S15D. DEG specificity
p1 <- deg_sig %>%
  # filter(species == "human") %>% 
  select(species, gene, sig_diff, gene_spec, marker) %>% unique() %>% 
  ggplot(aes(x = gene_spec, color = species)) +
  geom_density(adjust = 1.5) +
  labs(y = "Density", 
       x = "Cell type specificity")
plot(p1)

ggsave(p1, filename = "output/deg_specificity_types_species.pdf", width = 3.5, height = 3)


# Figure S15E. Human DEG specificity vs. log2FC
p1 <- deg_sig %>% 
  filter(species == "human") %>%
  # filter(abs(log2FoldChange) < 20) %>% 
  group_by(gene) %>% 
  mutate(deg_types = n()) %>%
  ggplot(aes(x = as_factor(deg_types), y = abs(log2FoldChange))) +
  geom_boxplot() +
  xlab("DEG specificity (number of consensus types)") +
  ylab("DEG magnitude (log2FC)") +
  ggtitle("Human") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
plot(p1)

ggsave(p1, filename = "output/human_deg_spec_log2FC.pdf", width = 10, height = 5)


# Figure S15F. Num. species with DEG (10x change)
p1 <- deg_sig %>% 
  group_by(gene) %>%
  mutate(num_lineage_10x = case_when(abs(log2FoldChange) > log2(10) ~ n_distinct(species))) %>% 
  select(gene, num_lineage, num_lineage_10x) %>%  
  pivot_longer(cols = c(num_lineage, num_lineage_10x), names_to = "fc", values_to = "species") %>% 
  drop_na() %>%
  unique() %>% 
  ggplot(aes(x = as_factor(species), fill = fc)) + 
  geom_bar(position = "identity") +
  scale_fill_brewer(palette = "Blues", name = "Fold-change") +
  labs(x = "No. of species",
       y = "No. of DEGs") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot(p1)

ggsave(p1, filename = "output/deg_number_of_species.pdf", width = 3.5, height = 3)




# GO analysis
go_class <- "cc"
combined_results <- readRDS(str_c("data/cross_species_cluster_GO", go_class, ".RDS"))

# Add cluster meta
cluster_meta <- read_csv(file = "data/cross_species_cluster_counts.csv", show_col_types = FALSE)
combined_results <- combined_results %>% 
  left_join(cluster_meta, by = "cluster")

#plot data
to_plot <- combined_results %>%
  filter(p.adjust < 0.005) %>%
  
  group_by(species, Description) %>%
  mutate(species_term_counts = n(), 
         species_term_mean_degs = mean(Count)) %>%
  ungroup() %>%
  
  arrange(species, desc(species_term_counts)) %>%
  
  distinct(ID, Description, species, species_term_counts, species_term_mean_degs) %>%
  
  group_by(ID) %>%
  mutate(term_count = n(), 
         term_human = if_else(any(species == "human"), "yes", "no")) %>%
  ungroup() %>%
  
  # filter(term_count > 1 | species == "human") %>%
  filter(term_human == "yes" & term_count > 0) %>%
  
  mutate(species = species %>% as_factor() %>% fct_relevel(c("human", "chimp", "gorilla", "rhesus", "marmoset")))


# GO semantic similarity
go_semsim <- read.table(str_c("data/pw_deg/20220328/go_", go_class, "_sem_sim.txt"), row.names = 1)
keep_goid <- match(unique(to_plot$ID), row.names(go_semsim))
go_semdist <- as.dist(1 - go_semsim[keep_goid, keep_goid])
hc1 <- hclust(go_semdist)

num_cl <- ifelse(go_class == "cc", 4, 6)
cl1 <- cutree(hc1, num_cl)
go_cl <- data.frame(ID = names(cl1), go_cluster = cl1)

go_id_order <- hc1$labels[hc1$order]
go_order <- to_plot$Description[match(go_id_order, to_plot$ID)]

to_plot <- to_plot %>%
  left_join(go_cl, by = "ID") %>% 
  arrange(species, desc(term_count)) %>% 
  mutate(Description = Description %>% as_factor() %>% fct_relevel(go_order))


# Figure 5D. DEG GO enrichment
p1 <- to_plot %>%
  ggplot() +
  geom_point(aes(x = species, y = Description, size = species_term_mean_degs, 
                 color = as_factor(go_cluster))) +
  labs(size = "Average No. of DEGs",
       color = "GO cluster", 
       x = "Species enriched gene list", 
       y = "", 
       title = go_class) +
  theme(aspect.ratio = 5,
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
plot(p1)

pdf(file = str_c(output_dir, "dotplot_GO_", go_class, ".pdf"),
    useDingbats = FALSE,
    height = 5,
    width = 12, 
    onefile = TRUE)
plot(as.dendrogram(hc1))
print(p1)
dev.off()



# Cluster counts for gene categories
combined_results %>% 
  filter(pvalue < 0.05 & grepl("ribosom", Description) & species == "human") %>% 
  group_by(cluster) %>% 
  count() %>% 
  arrange(-n)


# GO pvalue by cell type
go_keep <- data.frame(rbind(cbind("1_Ribosome", c("GO:0022626", "GO:0022625")), 
                            cbind("2_ECM", "GO:0062023"), 
                            cbind("3_Axon", c("GO:0030175", "GO:0044309", "GO:0030424")), 
                            cbind("4_Synapse", c("GO:0005911", "GO:0005912", "GO:0043197", "GO:0099572", 
                                                 "GO:0045211", "GO:0014069", "GO:0032279", "GO:0098978", 
                                                 "GO:0098984", "GO:0098794", "GO:0097060"))))
colnames(go_keep) <- c("go_set", "ID")

to_plot <- combined_results %>% 
  mutate(cluster = factor(cluster, levels = sort(unique(cluster))), 
         class = factor(class), 
         species = factor(species, levels = c("human", "chimp", "gorilla", "rhesus", "marmoset"))) %>%
  filter(ID %in% go_keep$ID) %>% 
  filter(pvalue < 0.05) %>%
  left_join(go_keep, by = "ID") %>% 
  group_by(species, cluster, go_set) %>% 
  mutate(min_p = min(pvalue))


# Figure S15H. Cell types with GO enrichment
p1 <- to_plot %>% 
  ggplot(aes(x = fct_reorder(fct_reorder(cluster, -as.numeric(cluster)), as.numeric(class)), 
             y = species, size = -log10(min_p), color = class)) +
  geom_point() +
  coord_flip() +
  facet_wrap(~go_set, nrow = 1) +
  guides(color = "none") +
  labs(x = "", 
       y = "") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
plot(p1)

ggsave(p1, filename = "output/GO_cc_sig0.05_cell_types.pdf", height = 8, width = 6)


# Figure 5E. Number of cell types per GO category
p2 <- to_plot %>% 
  group_by(species, go_set, class) %>% 
  filter(species == "human") %>%
  summarize(cnt = n_distinct(cluster)) %>% 
  ggplot(aes(x = fct_rev(go_set), y = cnt, fill = class)) +
  # facet_wrap(~ species) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(y = "Number of consensus types", 
       x = "")
plot(p2)

ggsave(p2, filename = "output/go_enrich_types.pdf", height = 3, width = 4)
