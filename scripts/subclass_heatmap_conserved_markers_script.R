library(tidyverse)
library(scrattch.hicat)
library(here)

#read in data and calculate expression means
orthologous_genes <- readRDS(here("data", "orthologous_genes.RDS"))

human_mat <- readRDS(here("data", "human_mat.RDS"))
human_meta <- readRDS(here("data", "human_meta.RDS"))
human_meta <- human_meta %>% filter(species_tech %>% str_detect("10x"))
human_mat <- human_mat[orthologous_genes$human_symbol, human_meta$sample_id]

chimp_mat <- readRDS(here("data", "chimp_mat.RDS"))
chimp_meta <- readRDS(here("data", "chimp_meta.RDS"))
chimp_meta <- chimp_meta %>% filter(species_tech %>% str_detect("10x"))
chimp_mat <- chimp_mat[orthologous_genes$chimp_symbol, chimp_meta$sample_id]
rownames(chimp_mat) <- orthologous_genes$human_symbol

gorilla_mat <- readRDS(here("data", "gorilla_mat.RDS"))
gorilla_meta <- readRDS(here("data", "gorilla_meta.RDS"))
gorilla_meta <- gorilla_meta %>% filter(species_tech %>% str_detect("10x"))
gorilla_mat <- gorilla_mat[orthologous_genes$gorilla_symbol, gorilla_meta$sample_id]
rownames(gorilla_mat) <- orthologous_genes$human_symbol

rhesus_mat <- readRDS(here("data", "rhesus_mat.RDS"))
rhesus_meta <- readRDS(here("data", "rhesus_meta.RDS"))
rhesus_meta <- rhesus_meta %>% filter(species_tech %>% str_detect("10x"))
rhesus_mat <- rhesus_mat[orthologous_genes$rhesus_symbol, rhesus_meta$sample_id]
rownames(rhesus_mat) <- orthologous_genes$human_symbol

marmoset_mat <- readRDS(here("data", "marmoset_mat.RDS"))
marmoset_meta <- readRDS(here("data", "marmoset_meta.RDS"))
marmoset_meta <- marmoset_meta %>% filter(species_tech %>% str_detect("10x"))
marmoset_mat <- marmoset_mat[orthologous_genes$marmoset_symbol, marmoset_meta$sample_id]
rownames(marmoset_mat) <- orthologous_genes$human_symbol

#calculate subclass means
cl <- human_meta$subclass
names(cl) <- human_meta$sample_id
human_mat <- cpm(human_mat)
human_mat@x <- log2(human_mat@x + 1)
human_means <- get_cl_means(human_mat, cl)

cl <- chimp_meta$subclass
names(cl) <- chimp_meta$sample_id
chimp_mat <- cpm(chimp_mat)
chimp_mat@x <- log2(chimp_mat@x + 1)
chimp_means <- get_cl_means(chimp_mat, cl)

cl <- gorilla_meta$subclass
names(cl) <- gorilla_meta$sample_id
gorilla_mat <- cpm(gorilla_mat)
gorilla_mat@x <- log2(gorilla_mat@x + 1)
gorilla_means <- get_cl_means(gorilla_mat, cl)

cl <- rhesus_meta$subclass
names(cl) <- rhesus_meta$sample_id
rhesus_mat <- cpm(rhesus_mat)
rhesus_mat@x <- log2(rhesus_mat@x + 1)
rhesus_means <- get_cl_means(rhesus_mat, cl)

cl <- marmoset_meta$subclass
names(cl) <- marmoset_meta$sample_id
marmoset_mat <- cpm(marmoset_mat)
marmoset_mat@x <- log2(marmoset_mat@x + 1)
marmoset_means <- get_cl_means(marmoset_mat, cl)

human_means <- human_means %>% as_tibble(rownames = "genes")
chimp_means <- chimp_means %>% as_tibble(rownames = "genes")
gorilla_means <- gorilla_means %>% as_tibble(rownames = "genes")
rhesus_means <- rhesus_means %>% as_tibble(rownames = "genes")
marmoset_means <- marmoset_means %>% as_tibble(rownames = "genes")

#read in species subclass marker gene lists

human <- readRDS(here("DEG_lists", "human_subclass_wilcox_markers.RDS")) 
chimp <- readRDS(here("DEG_lists", "chimp_subclass_wilcox_markers.RDS"))
gorilla <- readRDS(here("DEG_lists", "gorilla_subclass_wilcox_markers.RDS"))
rhesus <- readRDS(here("DEG_lists", "rhesus_subclass_wilcox_markers.RDS"))
marmoset <- readRDS(here("DEG_lists", "marmoset_subclass_wilcox_markers.RDS"))

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

#filter and format orthologous DEGs
genes_to_plot <- combined_genes %>%
    
    #filter sign genes
    filter(p_val_adj < 0.01,
           avg_log2FC > 1) %>%
    
    
    #filter DEGs in all 5 species
    group_by(cluster, gene) %>%
    mutate(gene_count = n()) %>%
    ungroup() %>%
    
    filter(gene_count == 5) %>%
    distinct(cluster, gene) %>%
    rename("marker_of" = "cluster")

#add expression values to marker list
combined_data <- genes_to_plot %>%
    left_join(
        human_means %>%
            set_names(str_c("humanxxx", colnames(human_means))) %>%
            rename("gene" = "humanxxxgenes"),
        by = "gene"
    ) %>%
    
    left_join(
        chimp_means %>%
            set_names(str_c("chimpxxx", colnames(chimp_means))) %>%
            rename("gene" = "chimpxxxgenes"),
        by = "gene"
    ) %>%
    
    left_join(
        gorilla_means %>%
            set_names(str_c("gorillaxxx", colnames(gorilla_means))) %>%
            rename("gene" = "gorillaxxxgenes"),
        by = "gene"
    ) %>%
    
    left_join(
        rhesus_means %>%
            set_names(str_c("rhesusxxx", colnames(rhesus_means))) %>%
            rename("gene" = "rhesusxxxgenes"),
        by = "gene"
    ) %>%
    
    left_join(
        marmoset_means %>%
            set_names(str_c("marmosetxxx", colnames(marmoset_means))) %>%
            rename("gene" = "marmosetxxxgenes"),
        by = "gene"
    ) 

to_plot <- combined_data %>%
    mutate(plot_label = str_c(marker_of, "_", gene)) %>%
    gather(key = "species_cluster", value = "expression", -c(marker_of, gene, plot_label)) %>%
    separate(species_cluster, into = c("species", "cluster"), sep = "xxx", remove = FALSE)


#format subclass means
subclass_order <- c("L5/6 NP", "L6 IT Car3", "L6 CT", "L6b", "L5 ET",
                    "L2/3 IT", "L6 IT", "L5 IT", "L4 IT",
                    "Lamp5_Lhx6", "Lamp5", "Sncg", "Vip", "Pax6", 
                    "Chandelier", "Pvalb", "Sst", "Sst Chodl",
                    "OPC", "Astro", "Oligo", "VLMC", "Endo", "Micro-PVM")

#wrangle data and order factors
to_plot <- to_plot %>%
    filter(cluster %in% subclass_order) %>%
    
    mutate(marker_of = marker_of %>% as_factor() %>% fct_relevel(subclass_order),
           cluster = cluster %>% as_factor() %>% fct_relevel(subclass_order),
           species = species %>% as_factor() %>% fct_relevel(c("human", "chimp", "gorilla", "rhesus", "marmoset"))) %>%
    
    arrange(marker_of, cluster, species) %>%
    
    mutate(plot_label = plot_label %>% as_factor(),
           species_cluster = species_cluster %>% as_factor()) %>%
    
    #scale data 
    group_by(species, gene) %>%
    mutate(min = min(expression),
           max = max(expression)) %>%
    ungroup() %>%
    
    mutate(scaled_expression = (expression - min) / (max - min) )

#make plot
p1 <- to_plot %>%
    ggplot() +
    geom_tile(aes(x = plot_label, y = species_cluster, fill = scaled_expression)) +
    
    scale_fill_gradient2(low = "lightblue", mid = "white", high = "darkred", midpoint = 0.5) +
    
    labs(x = str_c(to_plot %>% distinct(plot_label) %>% nrow(), " DEGs") ) +
    theme(
        aspect.ratio = 1/2,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 6))

p1
