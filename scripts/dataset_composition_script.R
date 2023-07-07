library(tidyverse)
library(here)

meta_human <- readRDS(here("data", "human_meta.RDS"))
meta_chimp <- readRDS(here("data", "chimp_meta.RDS"))
meta_gorilla <- readRDS(here("data", "gorilla_meta.RDS"))
meta_rhesus <- readRDS(here("data", "rhesus_meta.RDS"))
meta_marmoset <- readRDS(here("data", "marmoset_meta.RDS"))


meta <- bind_rows(meta_human, meta_chimp, meta_gorilla, meta_rhesus, meta_marmoset)

#plot dataset barplot
to_plot <- meta %>%
    mutate(species = species %>% as_factor()) %>%
    
    mutate(dataset_label = case_when(
        tech == "SSv4" ~ "SSv4 layer dissected",
        dataset == "AIBS_L5" ~ "Cv3 layer 5 only",
        TRUE ~ "Cv3 all layers"
    )) %>%
    
    mutate(dataset_label = dataset_label %>% 
               as_factor() %>% 
               fct_relevel(c("Cv3 all layers", "Cv3 layer 5 only", "SSv4 layer dissected")) %>% 
               fct_rev()) %>%
    
    group_by(species, dataset_label) %>%
    summarise(label_count = n()) %>%
    ungroup() %>%
    
    mutate(label_count_label = str_c(round(signif(label_count, digits = 3) / 1000, 1), "K"))

p1 <- to_plot %>%
    ggplot() +
        geom_col(aes(x = label_count, y = dataset_label, fill = dataset_label), 
                 color = "black", width = 0.9) +
        geom_text(aes(x = label_count, y = dataset_label, label = label_count_label)) +
    
        labs(x = "", 
             y = "") +
        
        scale_fill_manual(values = c("#18965A", "#236D8E", "#FFA500")) +
        scale_x_continuous(breaks = c(0, 50000, 100000, 150000),
                           labels = c("0", "50k", "100k", "150k"), 
                           limits = c(0, 150000)) +
    
        facet_wrap(~species, ncol = 1) +
    
    theme_light() +
    
    theme(
        axis.title.y = element_blank(),
        strip.background.x = element_rect(fill = "#2c3e50"),
        strip.text.x = element_text(face = "bold")
    )
p1

#plot dataset barplot
to_plot <- meta %>%
    mutate(species = species %>% as_factor()) %>%
    
    mutate(dataset_label = case_when(
        tech == "SSv4" ~ "SSv4 layer dissected",
        dataset == "AIBS_L5" ~ "Cv3 layer 5 only",
        TRUE ~ "Cv3 all layers"
    )) %>%
    
    mutate(dataset_label = dataset_label %>% 
               as_factor() %>% 
               fct_relevel(c("Cv3 all layers", "Cv3 layer 5 only", "SSv4 layer dissected")) %>% 
               fct_rev()) %>%
    
    distinct(dataset_label, species, donor, sex) %>%
    
    group_by(dataset_label, species, sex) %>%
    summarise(sex_count = n()) %>%
    ungroup()

p2 <- to_plot %>%
    ggplot() +
        geom_col(aes(x = sex_count, y = dataset_label, fill = sex), 
                 color = "black", width = 0.9) +
    
        labs(x = "", 
             y = "") +
        scale_fill_manual(values = c("#CD2990", "#3A5FCD")) +
    
        facet_wrap(~species, ncol = 1) +
    
    theme_light() +
    
    theme(
        axis.title.y = element_blank(),
        strip.background.x = element_rect(fill = "#2c3e50"),
        strip.text.x = element_text(face = "bold")
    )
