library(tidyverse)
library(Seurat)
library(here)


#load data
human <- readRDS(here("data", "human_meta.RDS"))
chimp <- readRDS(here("data", "chimp_meta.RDS"))
gorilla <- readRDS(here("data", "gorilla_meta.RDS"))
rhesus <- readRDS(here("data", "rhesus_meta.RDS"))
marmoset <- readRDS(here("data", "marmoset_meta.RDS"))

subtype_labels <- readRDS(here("data", "cux2_annotations.RDS"))

#plot stacked barplot
meta <- bind_rows(human, chimp, gorilla, rhesus, marmoset)


to_plot <- subtype_labels %>%
    mutate(species = species %>% as_factor() %>% fct_relevel(c("human", "chimp", "gorilla", "rhesus", "marmoset"))) %>%
    
    left_join(
        meta %>% select(sample_id, donor),
        by = "sample_id"
    )  %>%
    group_by(donor) %>%
    mutate(donor_count = n()) %>%
    ungroup() %>%
    
    group_by(donor, cux2_status) %>%
    mutate(donor_status_count = n()) %>%
    ungroup() %>%
    
    mutate(donor_status_prop = donor_status_count / donor_count) %>%
    
    distinct(cux2_status, species, donor, donor_status_prop) %>%
    
    group_by(species, cux2_status) %>%
    summarise(mean_prop = mean(donor_status_prop),
              sd_prop = sd(donor_status_prop)) %>%
    ungroup() %>%
    
    mutate(x_label = str_c(cux2_status, "_", species) %>% as_factor())
    
    
p1 <- to_plot %>%
    ggplot() +
    geom_col(aes(x = x_label, y = mean_prop, fill = cux2_status)) +
    geom_errorbar(aes(x = x_label, ymin = mean_prop, ymax = mean_prop + sd_prop)) +
    geom_abline(slope = 0, intercept = 0.5, color = "grey", linetype = 2) +
    geom_text(aes(x = x_label, y = mean_prop + sd_prop + 0.05, label = round(mean_prop, 2))) +
    
    labs(x = "") +
    
    facet_wrap(~species, nrow = 1, scales = "free_x") +
    
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = NA, color = "black"))


p1
