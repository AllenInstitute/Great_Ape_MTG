library(Seurat)
library(tidyverse)
library(here)

# load in metadata and dendrogram
meta <- readRDS(here("data", "human_meta.RDS"))
dend <- readRDS(here("data", "human_dend.RDS"))

dend_order <- labels(dend)    

#Plot 1: layer heatmap

#check all SS clusters are represented
ss_cl <- meta %>%
  filter(species_tech %>% str_detect("_ss")) %>%
  distinct(cluster) %>% pull() 
dend_order[which(dend_order %in% ss_cl == FALSE)]

#check layer info
meta %>% distinct(layer) 

#make plotting object
to_plot <- meta %>%
  filter(species_tech %>% str_detect("_ss")) %>%
  
  group_by(cluster) %>%
  mutate(cluster_n = n()) %>%
  ungroup() %>%
  
  group_by(cluster, layer) %>%
  mutate(cluster_layer_n = n()) %>%
  ungroup() %>%
  
  mutate(layer_prop = cluster_layer_n / cluster_n) %>%
  
  distinct(cluster, layer, layer_prop) %>%
  
  #factorize
  mutate(layer = layer %>% as_factor() %>% fct_relevel(c("1", "2", "3", "4", "5", "6")) %>% fct_rev(),
         cluster = cluster %>% as_factor() %>% fct_relevel(dend_order))

p1 <- to_plot %>%
  ggplot() +
  geom_tile(aes(x = cluster, y = layer, fill = layer_prop)) +
  
  scale_fill_gradient(low = "white", high = "black") +
  
  labs(fill = "Proportion", 
       y = "Layer",
       x = "") +
  
  theme(aspect.ratio = 6/length(dend_order),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid = element_blank())

p1

#Plot 2: proportions by class
to_plot <- meta %>%
  filter(species_tech %>% str_detect("10x"),
         dataset != "AIBS_L5") %>%
  
  group_by(donor, class) %>%
  mutate(donor_n = n()) %>%
  ungroup() %>%
  
  group_by(donor, cluster) %>%
  mutate(donor_cluster_n = n()) %>%
  ungroup() %>%
  
  mutate(cluster_prop = donor_cluster_n / donor_n * 100) %>%
  
  group_by(cluster) %>%
  summarise(mean_prop = mean(cluster_prop),
            sd_prop = sd(cluster_prop)) %>%
  ungroup() %>%
  
  mutate(cluster = cluster %>% as_factor() %>% fct_relevel(dend_order)) 

tmp <- meta %>% distinct(cluster, cluster_color)
colors_use <- tmp$cluster_color
names(colors_use) <- tmp$cluster

p2 <- to_plot %>%
  ggplot() +
  geom_errorbar(aes(x = cluster, ymin = mean_prop - sd_prop, ymax = mean_prop + sd_prop)) +
  geom_point(aes(x = cluster, y = mean_prop, color = cluster), size = 2, alpha = 0.6) +
  
  scale_color_manual(values = colors_use) +
  scale_y_log10() +
  
  labs(y = "Percent of class",
       x = "") +  
  
  theme(aspect.ratio = 12/length(dend_order),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
        legend.position = "none",
        panel.background = element_rect(color = "black", fill = "white"),
        panel.grid.major.y = element_line(color = "lightgrey", size = 1, linetype = 2),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

p2

#Plot 3: Stacked barplot dataset
to_plot <- meta %>% 
  
  mutate(
    dataset_label = case_when(
      dataset == "AIBS_L5" ~ "Cv3 layer 5 only",
      species_tech %>% str_detect("ss") ~ "SSv4 layer dissected",
      TRUE ~ "Cv3 all layers"),
    
    dataset_label  = dataset_label %>% as_factor() %>% fct_relevel(c("Cv3 layer 5 only", "Cv3 all layers", "SSv4 layer dissected")),
    
    cluster = cluster %>% as_factor() %>% fct_relevel(dend_order))


p3 <- to_plot %>%
  ggplot() +
  geom_bar(aes(x = cluster, fill = dataset_label), position = "fill") +
  
  scale_fill_manual(values = c("#236D8E", "#FFA500", "#18965A")) +
  
  labs(x = "",
       y = "") +
  
  theme(aspect.ratio = 6/length(dend_order),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
        legend.position = "none",
        panel.background = element_rect(color = "black", fill = "white"),
       )

p3


#Plot 4: Stacked barplot donor
to_plot <- meta %>% 
  
  mutate(cluster = cluster %>% as_factor() %>% fct_relevel(dend_order),
         donor = donor %>% as_factor() %>% fct_relevel(c("H200.1030", "H200.1025", "H19.30.002", "H19.30.001", "H18.30.002", "H200.1023", "H18.30.001")))

tmp <- meta %>% distinct(donor, donor_color)
colors_use <- tmp$donor_color
names(colors_use) <- tmp$donor

p4 <- to_plot %>%
  ggplot() +
  geom_bar(aes(x = cluster, fill = donor), position = "fill") +
  
  scale_fill_manual(values = colors_use) +
  
  labs(x = "",
       y = "") +
  
  theme(aspect.ratio = 6/length(dend_order),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
        panel.background = element_rect(color = "black", fill = "white"),
       )

p4


#Plot 5: Violin plot genes detected
to_plot <- meta %>% 
  filter(species_tech %>% str_detect("10x")) %>%
  mutate(cluster = cluster %>% as_factor() %>% fct_relevel(dend_order))


tmp <- meta %>% distinct(cluster, cluster_color)
colors_use <- tmp$cluster_color
names(colors_use) <- tmp$cluster


p5 <- to_plot %>%
  ggplot() +
  geom_violin(aes(x = cluster, y = nFeature_RNA, fill = cluster), scale = "width") +
  geom_abline(slope = 0, intercept = 3000, linetype = 2, color = "lightgrey") +
  geom_abline(slope = 0, intercept = 6000, linetype = 2, color = "lightgrey") +
  geom_abline(slope = 0, intercept = 9000, linetype = 2, color = "lightgrey") +
  
  
  scale_fill_manual(values = colors_use) +
  
  labs(x = "",
       y = "Genes detected (10x)") +
  
  theme(aspect.ratio = 6/length(dend_order),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
        panel.grid = element_blank(),
        legend.position = "none",
        panel.background = element_rect(color = "black", fill = "white"),
  )

p5
