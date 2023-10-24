#### Figure S18 ####
#### hDEGs near HARs and hCONDELs in SynGO dataset ####

library(stringr)
library(reshape2)
library(ggplot2)
library(here)

# syngo:
syngo_terms <- read.table(here("data/syngo_analysis", "syngo_terms_id.txt"), sep="\t", header=TRUE) 
syngo_genes <- read.table(here("data/syngo_analysis", "syngo_ontologies.txt"), sep="\t", header=TRUE) 
syngo_terms_to_plot <- syngo_terms[which(syngo_terms$syngo_level == "CC_Level_2"),2] # Here, choose a SynGO Level and Category.

# cell type list:
cell_types <- unlist(read.table(here("data/syngo_analysis", "cell_type_list.txt"), sep="\t", header=TRUE))

# summary function:
summary_function <- function(ct, SynGO_term, all_genes, genes_in_SynGO_term, syngo_terms) {
  return(data.frame(cell_type = ct,
                    SynGO_term_name = syngo_terms$syngo_term_name[syngo_terms$syngo_ref_id == SynGO_term],
                    SynGO_Ref_ID = SynGO_term,
                    SynGO_Level = syngo_terms$syngo_level[syngo_terms$syngo_ref_id == SynGO_term],
                    num_SynGO_term =  length(genes_in_SynGO_term), 
                    num_all_genes = length(all_genes$genes),
                    num_DEGs = sum(all_genes$DEG),
                    num_DEGs_in_SynGO = sum(all_genes$DEG & all_genes$SynGO),
                    num_DEGs_HARs_hCONDELs_in_SynGO = sum(all_genes$DEG_and_HARorhCONDEL & all_genes$SynGO),
                    ratio_DEGs_per_SynGO = sum(all_genes$DEG & all_genes$SynGO)/length(genes_in_SynGO_term),
                    ratio_DEGs_HARs_hCONDELs_per_DEGs_per_SynGO = sum(all_genes$DEG_and_HARorhCONDEL & all_genes$SynGO)/length(genes_in_SynGO_term),
                    names_DEGs_in_SynGO = str_flatten(all_genes$genes[all_genes$DEG == TRUE & all_genes$SynGO == TRUE], collapse=";"),
                    names_DEGs_HARs_hCONDELs_in_SynGO = str_flatten(all_genes$genes[all_genes$DEG_and_HARorhCONDEL == TRUE & all_genes$SynGO == TRUE], collapse=";")
  ))
}

# initialize output table
col.names = c("cell_type", 
              "SynGO_term_name", 
              "SynGO_Ref_ID", 
              "SynGO_Level",
              "num_SynGO_term", 
              "num_all_genes",
              "num_DEGs", 
              "num_DEGs_in_SynGO", 
              "num_DEGs_HARs_hCONDELs_in_SynGO", 
              "ratio_DEGs_per_SynGO",
              "ratio_DEGs_HARs_hCONDELs_per_DEGs_per_SynGO",
              "names_DEGs_in_SynGO",
              "names_DEGs_HARs_hCONDELs_in_SynGO")
colClasses = c(rep("character",4), rep("integer", 5), rep("numeric",2), rep("character",2))
output_table <- read.table(text = "", colClasses = colClasses, col.names = col.names)

# for loop cell types:
for(ct in cell_types) {
  sig_genes <- read.table(here("data/human_differential_expression", paste0(ct, "_DEGs_HARs_hCONDELs.txt")), sep="\t", header=TRUE)
  sig_genes$HARs_count<-0
  idx<-which(!is.na(sig_genes$HARs.FlankingandIntron))
  valx<-unlist(lapply(strsplit(sig_genes$HARs.FlankingandIntron[idx], split=","),length))
  sig_genes$HARs_count[idx]<-valx
  sig_genes$hCONDELs_count<-0
  idy<-which(!is.na(sig_genes$hCONDELs.FlankingandIntron))
  valy<-unlist(lapply(strsplit(sig_genes$hCONDELs.FlankingandIntron[idy], split=","),length))
  sig_genes$hCONDELs_count[idy]<-valy
  sig_genes$sum_HARs_hCONDELs<-sig_genes$HARs_count+sig_genes$hCONDELs_count
  all_genes <- read.table(here("data/human_differential_expression",paste0(ct,"_human_vs_all_markers.txt")), sep="\t", header=TRUE)
  all_genes <- na.omit(all_genes)
  all_genes$DEG <- ifelse(all_genes$genes %in% sig_genes$genes, TRUE, FALSE)
  all_genes$DEG_and_HARorhCONDEL <- ifelse(all_genes$genes %in% sig_genes$genes[sig_genes$HARorhCONDEL == TRUE], TRUE, FALSE)
  # initialize ct output table:
  ct_output_table <- read.table(text = "", colClasses = colClasses, col.names = col.names)
  for(SynGO_term in syngo_terms_to_plot) {
    genes_in_SynGO_term <- unlist(str_split(syngo_genes$genes...hgnc_symbol[syngo_genes$user.interface.reference.code == SynGO_term], ";"))
    all_genes$SynGO <- ifelse(all_genes$genes %in% genes_in_SynGO_term, TRUE, FALSE)
    sig_genes$SynGO <- ifelse(sig_genes$genes %in% genes_in_SynGO_term, TRUE, FALSE)
    ct_output_table <- rbind(ct_output_table, 
                             summary_function(ct, SynGO_term, all_genes, genes_in_SynGO_term, syngo_terms))
  }
  output_table <- rbind(output_table, ct_output_table)
}

# plot data:
order_celltypes <- rev(cell_types)
output_table$cell_type <- factor(output_table$cell_type, levels = order_celltypes)
output_table$SynGO_Ref_ID <- factor(output_table$SynGO_Ref_ID, levels = syngo_terms_to_plot)
ggplot()+
  geom_point(data=output_table, aes(x=SynGO_Ref_ID,y=cell_type,size=num_DEGs_in_SynGO,color=num_DEGs_HARs_hCONDELs_in_SynGO),shape=16,alpha=1)+
  scale_color_gradientn(colours=c("#455bc7","#646bee","#916dc2","#bf6e97","#ee7370","#f08e6b","#f3aa67","#f7c666"), breaks=seq(0,14,2), limits = c(0, 13))+
  scale_size_continuous(range = c(-1, 3.5), breaks = c(5,10,15,20,25,30), limits = c(0, 30))+
  theme_bw()+
  theme(axis.text.x=element_text(color="#000000",family="Helvetica",angle=90,hjust=1,vjust=0.5),
        axis.text.y=element_text(color="#000000",family="Helvetica"),
        axis.ticks = element_line(color="#000000", size=0.25),
        panel.grid.major = element_line(size = 0.25))+
  labs(color = "Near HARs/hCONDELs", size = "No. hDEGs")
ggsave("hdegs_har_hcondel_syngo.pdf")
