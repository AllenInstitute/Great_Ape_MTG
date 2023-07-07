library(tidyverse)
library(here)

# load in metadata 

human <- readRDS(here("data", "human_meta.RDS"))
chimp <- readRDS(here("data", "chimp_meta.RDS"))
gorilla <- readRDS(here("data", "gorilla_meta.RDS"))
rhesus <- readRDS(here("data", "rhesus_meta.RDS"))
marmoset <- readRDS(here("data", "marmoset_meta.RDS"))

meta <- bind_rows(human, chimp, gorilla, rhesus, marmoset)

#Plot 2: proportions by class
to_plot <- meta %>%
    filter(species_tech %>% str_detect("10x"),
           dataset != "AIBS_L5") %>%
    
    group_by(donor, class) %>%
    mutate(donor_n = n()) %>%
    ungroup() %>%
    
    group_by(donor, subclass) %>%
    mutate(donor_subclass_n = n()) %>%
    ungroup() %>%
    
    mutate(subclass_prop = donor_subclass_n / donor_n * 100) %>%
    
    group_by(species, subclass) %>%
    summarise(mean_prop = mean(subclass_prop),
              sd_prop = sd(subclass_prop)) %>%
    ungroup() %>%
    
    left_join(meta %>% distinct(subclass, class), by = "subclass") %>%
    
    mutate(subclass = subclass %>% as_factor() %>% fct_relevel(
        c("L5/6 NP", "L6 IT Car3", "L6 CT", "L6b", "L5 ET",
          "L2/3 IT", "L6 IT", "L5 IT", "L4 IT",
          "Lamp5_Lhx6", "Lamp5", "Sncg", "Vip", "Pax6", 
          "Chandelier", "Pvalb", "Sst", "Sst Chodl", 
          "OPC", "Astro", "Oligo", "VLMC", "Endo", "Micro-PVM")
    )) %>% 
    
    mutate(class = class %>% as_factor() %>% fct_relevel(c("exc", "inh", "glia")),
           species = species %>% as_factor() %>% fct_relevel(c("human", "chimp", "gorilla", "rhesus", "marmoset"))) %>%
    arrange(subclass, species) %>%
    mutate(species_subclass = str_c(species, "_", subclass) %>% as_factor())

tmp <- meta %>% distinct(subclass, subclass_color)
colors_use <- tmp$subclass_color
names(colors_use) <- tmp$subclass

p1 <- to_plot %>%
    ggplot() +
    geom_errorbar(aes(x = species_subclass, ymin = mean_prop - sd_prop, ymax = mean_prop + sd_prop)) +
    geom_point(aes(x = species_subclass, y = mean_prop, color = subclass), size = 2, alpha = 0.6, shape = 16) +
    
    scale_color_manual(values = colors_use) +
    scale_y_log10(breaks = c(0, 0.1, 1, 10, 100),
                  limits = c(0.1, 100)) +
    
    labs(y = "Percent of class",
         x = "") +  
    
    facet_wrap(~class, scales = "free", ncol = 1) +
    theme(aspect.ratio = 1/3,
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
          legend.position = "none",
          panel.background = element_rect(color = "black", fill = "white"),
          panel.grid.major.y = element_line(color = "lightgrey", size = 1, linetype = 2),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())

p1


#run statistics
#find donor, species, subclass average proportion of class
group_avg <- meta %>%
    filter(species_tech %>% str_detect("10x"),
           dataset != "AIBS_L5") %>%
    
    group_by(donor, class) %>%
    mutate(donor_n = n()) %>%
    ungroup() %>%
    
    group_by(donor, subclass) %>%
    mutate(donor_subclass_n = n()) %>%
    ungroup() %>%
    
    mutate(subclass_prop = donor_subclass_n / donor_n * 100) %>%
    distinct(species, donor, subclass, subclass_prop) %>%
    arrange(subclass, species)
group_avg


#Run one-way ANOVA
model <- aov(data = group_avg, subclass_prop ~ species * subclass)
summary(model)

results <- tibble(p_val = 99.9, comparison = "blank", subclass = "blank")
subclasses <- group_avg %>% distinct(subclass) 

for(i in 1:nrow(subclasses)){
    subclass_oi <- subclasses$subclass[i]
    tmp <- group_avg %>% filter(subclass == subclass_oi)
    tmp <- pairwise.t.test(tmp$subclass_prop, tmp$species, p.adj = "none")
    
    tmp <- tmp$p.value %>% 
        as_tibble(rownames = "group1") %>%
        gather(key = "group2", value = "p_val", - group1) %>%
        mutate(comparison = str_c(group1, group2)) %>%
        filter(comparison %>% str_detect("human"),
               p_val < 0.05) %>%
        mutate(comparison = comparison %>% str_remove_all("human")) %>%
        select(-c(group1, group2)) %>%
        mutate(subclass = subclass_oi)
     results <- bind_rows(results, tmp)
    
}
results <- results %>% filter(comparison != "blank")
results


results$p_val_adj_fdr <-  p.adjust(p = results$p_val, method = "fdr", n = 4 * 24) #correct for comparisons to human- 4 species 24 subclasses
results

results %>%
    filter(p_val_adj_fdr < 0.05)
