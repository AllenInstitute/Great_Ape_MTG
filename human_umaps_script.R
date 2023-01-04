library(tidyverse)
library(scattermore)
library(Seurat)
library(here)

#load and wrangle data
meta <- readRDS(here("data", "human_meta.RDS"))
reductions <- readRDS("data", "human_all_umap_coords_2d.RDS"))

to_plot <- reductions %>%
    as_tibble(rownames = "sample_id") %>%
    filter(sample_id %in% meta$sample_id) %>%
    left_join(meta, by = "sample_id")

#plot 10x nuclei by cluster
tmp <- meta %>% distinct(cluster, cluster_color)
cluster_colors <- tmp$cluster_color
names(cluster_colors) <- tmp$cluster

p1 <- to_plot %>%
    filter(species_tech %>% str_detect("10x")) %>%
    sample_n(n()) %>%
    ggplot() +
    geom_scattermore(aes(x = UMAP_1, y = UMAP_2, color = cluster), pointsize = 1) +
    scale_color_manual(values = cluster_colors) +
    
    labs(x = "",
         y = "") +
    
    theme(aspect.ratio = 1,
          legend.position = "none",
          panel.background = element_rect(color = "black", fill = "white"),
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
p1

#plot SS nuclei by cluster
p2 <- to_plot %>%
    filter(species_tech %>% str_detect("ss")) %>%
    sample_n(n()) %>%
    ggplot() +
    geom_scattermore(aes(x = UMAP_1, y = UMAP_2, color = cluster), pointsize = 1) +
    scale_color_manual(values = cluster_colors) +
    
    labs(x = "",
         y = "") +
    
    theme(aspect.ratio = 1,
          legend.position = "none",
          panel.background = element_rect(color = "black", fill = "white"),
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
p2

#plot all nuclei by donor
tmp <- meta %>% distinct(donor, donor_color)
donor_colors <- tmp$donor_color
names(donor_colors) <- tmp$donor

p3 <- to_plot %>%

    sample_n(n()) %>%
    ggplot() +
    geom_scattermore(aes(x = UMAP_1, y = UMAP_2, color = donor), pointsize = 1) +
    scale_color_manual(values = donor_colors) +
    
    labs(x = "",
         y = "") +
    
    theme(aspect.ratio = 1,
          legend.position = "none",
          panel.background = element_rect(color = "black", fill = "white"),
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
p3

#plot SS nuclei by layer
p4 <- to_plot %>%
    filter(species_tech %>% str_detect("_ss")) %>%
    sample_n(n()) %>%
    
    ggplot() +
    geom_scattermore(aes(x = UMAP_1, y = UMAP_2, color = layer), pointsize = 1) +
    scale_color_brewer(palette = "RdBu", direction = -1) +
    
    labs(x = "",
         y = "",
         color = "Layer") +
    
    theme(aspect.ratio = 1,
          
          panel.background = element_rect(color = "black", fill = "white"),
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
p4
