library(tidyverse)
library(scrattch.hicat)
library(scattermore)
library(here)

#load data and calculate cl means/medians
orthologous_genes <- readRDS(here("data", "orthologous_genes.RDS"))

#human
human_mat <- readRDS(here("data", "human_mat.RDS"))
human_meta <- readRDS(here("data", "human_meta.RDS"))

human_meta <- human_meta %>% 
    filter(species_tech %>% str_detect("10x"))

human_mat <- human_mat[ , human_meta$sample_id]
human_mat <- cpm(human_mat)
human_mat@x <- log2(human_mat@x + 1)
human_mat <- human_mat[orthologous_genes$human_symbol, ]

cl <- human_meta$subclass
names(cl) <- human_meta$sample_id

human_medians <- get_cl_medians(human_mat, cl)
human_means <- get_cl_means(human_mat, cl)
rm(human_mat)
gc()

#chimp
chimp_mat <- readRDS(here("data", "chimp_mat.RDS"))
chimp_meta <- readRDS(here("data", "chimp_meta.RDS"))

chimp_meta <- chimp_meta %>% 
    filter(species_tech %>% str_detect("10x"))

chimp_mat <- chimp_mat[ , chimp_meta$sample_id]
chimp_mat <- cpm(chimp_mat)
chimp_mat@x <- log2(chimp_mat@x + 1)
chimp_mat <- chimp_mat[orthologous_genes$chimp_symbol, ]
rownames(chimp_mat) <- orthologous_genes$human_symbol

cl <- chimp_meta$subclass
names(cl) <- chimp_meta$sample_id

chimp_medians <- get_cl_medians(chimp_mat, cl)
chimp_means <- get_cl_means(chimp_mat, cl)
rm(chimp_mat)
gc()

#gorilla
gorilla_mat <- readRDS(here("data", "gorilla_mat.RDS"))
gorilla_meta <- readRDS(here("data", "gorilla_meta.RDS"))

gorilla_meta <- gorilla_meta %>% 
    filter(species_tech %>% str_detect("10x"))

gorilla_mat <- gorilla_mat[ , gorilla_meta$sample_id]
gorilla_mat <- cpm(gorilla_mat)
gorilla_mat@x <- log2(gorilla_mat@x + 1)
gorilla_mat <- gorilla_mat[orthologous_genes$gorilla_symbol, ]
rownames(gorilla_mat) <- orthologous_genes$human_symbol

cl <- gorilla_meta$subclass
names(cl) <- gorilla_meta$sample_id

gorilla_medians <- get_cl_medians(gorilla_mat, cl)
gorilla_means <- get_cl_means(gorilla_mat, cl)
rm(gorilla_mat)
gc()

#rhesus
rhesus_mat <- readRDS(here("data", "rhesus_mat.RDS"))
rhesus_meta <- readRDS(here("data", "rhesus_meta.RDS"))

rhesus_meta <- rhesus_meta %>% 
    filter(species_tech %>% str_detect("10x"))

rhesus_mat <- rhesus_mat[ , rhesus_meta$sample_id]
rhesus_mat <- cpm(rhesus_mat)
rhesus_mat@x <- log2(rhesus_mat@x + 1)
rhesus_mat <- rhesus_mat[orthologous_genes$rhesus_symbol, ]
rownames(rhesus_mat) <- orthologous_genes$human_symbol

cl <- rhesus_meta$subclass
names(cl) <- rhesus_meta$sample_id

rhesus_medians <- get_cl_medians(rhesus_mat, cl)
rhesus_means <- get_cl_means(rhesus_mat, cl)
rm(rhesus_mat)
gc()

#marmoset
marmoset_mat <- readRDS(here("data", "marmoset_mat.RDS"))
marmoset_meta <- readRDS(here("data", "marmoset_meta.RDS"))

marmoset_meta <- marmoset_meta %>% 
    filter(species_tech %>% str_detect("10x"))

marmoset_mat <- marmoset_mat[ , marmoset_meta$sample_id]
marmoset_mat <- cpm(marmoset_mat)
marmoset_mat@x <- log2(marmoset_mat@x + 1)
marmoset_mat <- marmoset_mat[orthologous_genes$marmoset_symbol, ]
rownames(marmoset_mat) <- orthologous_genes$human_symbol

cl <- marmoset_meta$subclass
names(cl) <- marmoset_meta$sample_id

marmoset_medians <- get_cl_medians(marmoset_mat, cl)
marmoset_means <- get_cl_means(marmoset_mat, cl)
rm(marmoset_mat)
gc()



#combined data
combine_means_and_medians <- function(means, medians, species_name){
    
    means %>%
        as_tibble(rownames = "genes") %>%
        
        gather(key = "subclass",
               value = "mean_expression",
               -genes) %>%
        
        mutate(
            species = species_name,
            link = str_c(genes, subclass)) %>%
        
        left_join(
            medians %>%
                as_tibble(rownames = "genes") %>%
                
                gather(key = "subclass",
                       value = "median_expression",
                       -genes) %>%
                
                mutate(link = str_c(genes, subclass)) %>%
                select(-c(genes, subclass)),
            
            by = "link"
        )
    
}

human_combined <- combine_means_and_medians(human_means, human_medians, "human")    
chimp_combined <- combine_means_and_medians(chimp_means, chimp_medians, "chimp")    
gorilla_combined <- combine_means_and_medians(gorilla_means, gorilla_medians, "gorilla")    
rhesus_combined <- combine_means_and_medians(rhesus_means, rhesus_medians, "rhesus")    
marmoset_combined <- combine_means_and_medians(marmoset_means, marmoset_medians, "marmoset")    

combined_data <- bind_rows(human_combined, chimp_combined, gorilla_combined, rhesus_combined, marmoset_combined)

#find pairwise correlations
correlation_function <- function(species_1, species_2, celltype){
    #filter to only genes with cl medians >0 in both species for spearman corr
    tmp <- combined_data %>%
        filter(species %in% c(species_1, species_2),
               subclass == celltype, 
               median_expression != 0) %>%
        
        group_by(genes) %>%
        mutate(gene_counts = n()) %>%
        ungroup() %>%
        filter(gene_counts == 2)
    
    data_1 <- tmp %>%
        filter(species == species_1,
               subclass == celltype) %>%
        arrange(genes)
    
    data_2 <- tmp %>%
        filter(species == species_2,
               subclass == celltype) %>%
        arrange(genes)
    
    cor(data_1$median_expression, data_2$median_expression, method = "spearman")
}

correlation_results <- tibble(species_1 = rep(rep(c("human", "chimp", "gorilla", "rhesus", "marmoset"), 5), 24),
                              species_2 = rep(c(rep("human", 5), rep("chimp", 5), rep("gorilla", 5), rep("rhesus", 5), rep("marmoset", 5)), 24),
                              celltype = sort(rep(unique(human_combined$subclass), 25))) %>%
    mutate(spearman_corr = 999)

for(i in 1:nrow(correlation_results)){
    correlation_results$spearman_corr[i] <- correlation_function(correlation_results$species_1[i], 
                                                                 correlation_results$species_2[i], 
                                                                 correlation_results$celltype[i])
}

#plot results
p1 <- correlation_results %>%
    mutate(species_1 = species_1 %>% as_factor() %>% fct_relevel(c("human", "chimp", "gorilla", "rhesus", "marmoset")),
           species_2 = species_2 %>% as_factor() %>% fct_relevel(c("human", "chimp", "gorilla", "rhesus", "marmoset")),
           celltype = celltype %>% as_factor() %>% fct_relevel(c("L5/6 NP", "L6 IT Car3", "L6 CT", "L6b", "L5 ET",
                                                                 "L2/3 IT", "L6 IT", "L5 IT", "L4 IT",
                                                                 "Lamp5_Lhx6", "Lamp5", "Sncg", "Vip", "Pax6", 
                                                                 "Chandelier", "Pvalb", "Sst", "Sst Chodl", 
                                                                 "OPC", "Astro", "Oligo", "VLMC", "Endo", "Micro-PVM"))
    ) %>%
    
    filter(!(celltype %in% c("VLMC", "Endo"))) %>%
    
    ggplot() +
    geom_tile(aes(x = species_1, y = species_2 %>% fct_rev(), fill = spearman_corr)) +
    
    scale_fill_viridis_c(option = "B") +
    
    labs(x = "",
         y = "", 
         fill = "Spearman corr.",
         caption = "Spearman correlation between species subclass medians.\nOnly genes with >0 median value in both species were used.") +
    
    facet_wrap(~celltype, ncol = 9) +
    
    theme(aspect.ratio = 1,
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2))

p1
