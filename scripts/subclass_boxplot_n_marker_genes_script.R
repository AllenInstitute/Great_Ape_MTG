library(tidyverse)
library(here)

#read in species subclass marker gene lists
human <- readRDS(here("DEG_lists", "human_subclass_wilcox_markers.RDS")) 
chimp <- readRDS(here("DEG_lists", "chimp_subclass_wilcox_markers.RDS"))
gorilla <- readRDS(here("DEG_lists", "gorilla_subclass_wilcox_markers.RDS"))
rhesus <- readRDS(here("DEG_lists", "rhesus_subclass_wilcox_markers.RDS"))
marmoset <- readRDS(here("DEG_lists", "marmoset_subclass_wilcox_markers.RDS"))

#format gene lists
format_deg_lists <- function(meta){
    meta %>% as_tibble() %>%
        filter(p_val_adj < 0.01, 
               avg_log2FC > 1) %>%
        distinct(cluster, gene) %>%
        
        group_by(cluster) %>%
        summarise(gene_count = n()) %>%
        ungroup()
}  

to_plot <- bind_rows(
    format_deg_lists(human) %>% mutate(species = "human"),
    format_deg_lists(chimp) %>% mutate(species = "chimp"),
    format_deg_lists(gorilla) %>% mutate(species = "gorilla"),
    format_deg_lists(rhesus) %>% mutate(species = "rhesus"),
    format_deg_lists(marmoset) %>% mutate(species = "marmoset")
)

#find number of conserved subclass marker genes to report as separate point on boxplot
orthologous_genes <- readRDS(here("data", "orthologous_genes.RDS"))

combined_genes <- bind_rows(human %>%
                                as_tibble() %>%
                                filter(gene %in% orthologous_genes$human_symbol) %>%
                                mutate(species = "human"),
                            
                            chimp %>%
                                as_tibble() %>%
                                filter(gene %in% orthologous_genes$chimp_symbol) %>%
                                mutate(species = "chimp") %>%
                                left_join(
                                    orthologous_genes %>% as_tibble() %>%
                                        select(human_symbol, chimp_symbol) %>%
                                        set_names("new_symbol", "gene"),
                                    by = "gene"
                                ) %>%
                                select(-gene) %>%
                                rename("gene" = "new_symbol"),
                            
                            gorilla %>%
                                as_tibble() %>%
                                filter(gene %in% orthologous_genes$gorilla_symbol) %>%
                                mutate(species = "gorilla") %>%
                                left_join(
                                    orthologous_genes %>% as_tibble() %>%
                                        select(human_symbol, gorilla_symbol) %>%
                                        set_names("new_symbol", "gene"),
                                    by = "gene"
                                ) %>%
                                select(-gene) %>%
                                rename("gene" = "new_symbol"),
                            
                            rhesus %>%
                                as_tibble() %>%
                                filter(gene %in% orthologous_genes$rhesus_symbol) %>%
                                mutate(species = "rhesus") %>%
                                left_join(
                                    orthologous_genes %>% as_tibble() %>%
                                        select(human_symbol, rhesus_symbol) %>%
                                        set_names("new_symbol", "gene"),
                                    by = "gene"
                                ) %>%
                                select(-gene) %>%
                                rename("gene" = "new_symbol"),
                            
                            marmoset %>%
                                as_tibble() %>%
                                filter(gene %in% orthologous_genes$marmoset_symbol) %>%
                                mutate(species = "marmoset") %>%
                                left_join(
                                    orthologous_genes %>% as_tibble() %>%
                                        select(human_symbol, marmoset_symbol) %>%
                                        set_names("new_symbol", "gene"),
                                    by = "gene"
                                ) %>%
                                select(-gene) %>%
                                rename("gene" = "new_symbol")
)

n_conserved_genes <- combined_genes %>%
    filter(avg_log2FC > 1,
           p_val_adj < 0.01) %>%
    
    group_by(cluster, gene) %>%
    mutate(gene_count = n()) %>%
    ungroup() %>%
    
    filter(gene_count == 5) %>%
    distinct(cluster, gene) %>%
    
    group_by(cluster) %>%
    summarise(n_conserved = n())

#plot boxplot
subclass_order <- c("L5/6 NP", "L6 IT Car3", "L6 CT", "L6b", "L5 ET",
                    "L2/3 IT", "L6 IT", "L5 IT", "L4 IT",
                    "Lamp5_Lhx6", "Lamp5", "Sncg", "Vip", "Pax6", 
                    "Chandelier", "Pvalb", "Sst", "Sst Chodl",
                    "OPC", "Astro", "Oligo", "VLMC", "Endo", "Micro-PVM")

n_conserved_genes <- n_conserved_genes %>%
    mutate(cluster = cluster %>% as_factor() %>% fct_relevel(subclass_order))

p1 <- to_plot %>% 
    mutate(species = species %>% as_factor() %>% fct_relevel(c("human", "chimp", "gorilla", "rhesus", "marmoset")),
           cluster = cluster %>% as_factor() %>% fct_relevel(subclass_order)) %>%
    
    ggplot() +
    geom_boxplot(aes(x = cluster, y = gene_count), outlier.shape = NA) +
    geom_jitter(aes(x = cluster, y = gene_count, color = species), width = 0.2, size = 3, shape = 16) +
    geom_point(data = n_conserved_genes, aes(x = cluster, y = n_conserved), size = 3, shape = 17) +
    
    scale_color_manual(values = c("royalblue1", "#387D7A", "sienna2", "darkgreen", "maroon4")) + 
    scale_y_log10(breaks = c(1, 10, 100, 1000),
                  labels = c(1, 10, 100, 1000),
                  limits = c(1,1200)) +
    
    labs(y = "Number of subclass\nmarker genes",
         x = "",
         color = "") +
    
    theme(
        aspect.ratio = 1/2,
        panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "lightgrey"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )
p1

#plot density plot
p2 <- to_plot %>%
    mutate(class = case_when(
        cluster %in% subclass_order[1:5] ~ "Exc: non-IT",
        cluster %in% subclass_order[6:9] ~ "Exc: IT",
        cluster %in% subclass_order[10:18] ~ "Inh",
        cluster %in% subclass_order[19:24] ~ "NN"
    )) %>%
    
    ggplot() +
    geom_density(aes(y = gene_count, fill = class), alpha = 0.3) +
    geom_point(aes(x = 0, y = gene_count, color = class)) +
    
    scale_fill_manual(values = c("#4D9EC3", "#8CBB25", "#B46C00", "#574E3C")) + 
    scale_color_manual(values = c("#4D9EC3", "#8CBB25", "#B46C00", "#574E3C")) + 
    
    scale_y_log10(breaks = c(1, 10, 100, 1000),
                  labels = c(1, 10, 100, 1000),
                  limits = c(1,1200)) +
    
    labs(x = "", 
         y = "density") +
    
    theme(
        aspect.ratio = 3,
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )

p2
