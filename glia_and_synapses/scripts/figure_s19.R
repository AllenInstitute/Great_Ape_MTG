#### Figure S19 ####
#### Human gene expression divergence of synaptic compartments and processes across consensus cell types ####

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
summary_function <- function(genes, ct, SynGO_term_id) {
  return(c(ct, SynGO_term_id,
           length(genes$genes[genes$log2FoldChange > 0 & genes$SynGO==TRUE]),
           length(genes$genes[genes$log2FoldChange < 0 & genes$SynGO==TRUE]),
           median(genes$log2FoldChange[genes$log2FoldChange > 0 & genes$SynGO==TRUE]),
           median(abs(genes$log2FoldChange[genes$log2FoldChange < 0 & genes$SynGO==TRUE]))
  ))
}

# initialize output table:
col.names = c("Cell_type", "SynGO_term_id", "up_Number_hDEGs", "down_Number_hDEGs", "up_Median_abs_log2FC", "down_Median_abs_log2FC")
colClasses = c(rep("character", 2), rep("numeric", 4))
output_table <- read.table(text = "",
                           colClasses = colClasses,
                           col.names = col.names)
output_table <- rbind(output_table,col.names)
colnames(output_table) <- output_table[1,]

# for loop cell types:
for(ct in cell_types){
  genes <- read.table(here("data/human_differential_expression",paste0(ct,"_human_vs_all_markers.txt")), sep="\t", header=TRUE)
  genes$DEG <- ifelse(abs(genes$log2FoldChange) > 0.5 & genes$padj < 0.05, "TRUE", "FALSE")
  genes$DEG[is.na(genes$DEG)] <- FALSE
  genes <- genes[genes$DEG==TRUE,]
  for(SynGO_term_id in syngo_terms_to_plot){
    genes_in_SynGO_term <- unlist(str_split(syngo_genes$genes...hgnc_symbol[syngo_genes$user.interface.reference.code == SynGO_term_id], ";"))
    genes$SynGO <- ifelse(genes$genes %in% genes_in_SynGO_term, TRUE, FALSE)
    output_table <- rbind(output_table, summary_function(genes, ct, SynGO_term_id))
  }
}
output_table[is.na(output_table)] = 0
output_table <- output_table[-1,]

# data table:
#gene number
table_numbers <- output_table[c(1,2,3,4)]
table_numbers <- melt(table_numbers, id.vars=c('Cell_type', 'SynGO_term_id'),var='gene_group')
table_numbers$value <- as.numeric(table_numbers$value)
colnames(table_numbers)[4] <- "Number_hDEGs"
table_numbers$DEG_direction[grepl("up", table_numbers$gene_group)] <- "Up"
table_numbers$DEG_direction[grepl("down", table_numbers$gene_group)] <- "Down"
table_numbers <- table_numbers[-3]
#fold change
table_log2FC <- output_table[c(1,2,5,6)]
table_log2FC <- melt(table_log2FC, id.vars=c('Cell_type', 'SynGO_term_id'),var='gene_group')
table_log2FC$value <- as.numeric(table_log2FC$value)
colnames(table_log2FC)[4] <- "log2FC"
table_log2FC$DEG_direction[grepl("up", table_log2FC$gene_group)] <- "Up"
table_log2FC$DEG_direction[grepl("down", table_log2FC$gene_group)] <- "Down"
table_log2FC <- table_log2FC[-3]
#merge
table_hDEGs <- merge(table_numbers,table_log2FC)

# plot:
order_celltypes <- rev(cell_types)
table_hDEGs$Cell_type <- factor(table_hDEGs$Cell_type, levels = order_celltypes)
table_hDEGs$SynGO_term_id <- factor(table_hDEGs$SynGO_term_id, levels = syngo_terms_to_plot)
ggplot()+
  geom_point(data=table_hDEGs, aes(y=Cell_type, x=log2FC, size=Number_hDEGs, color=DEG_direction), shape=16, alpha=0.8)+
  facet_grid(. ~ SynGO_term_id)+
  scale_size_continuous(range = c(-0.1, 4.5), breaks = c(1,5,10,15,20,25))+ #, limits=c(0,54)
  scale_x_continuous(expand=c(0,0), limits=c(0,7.26))+
  theme_bw()+
  theme(axis.text=element_text(color="#000000",family="Helvetica"),
        text=element_text(color="#000000",family="Helvetica"),
        axis.ticks = element_line(color="#000000", size=0.25),
        panel.spacing.x = unit(0, "lines"),
        panel.grid.major = element_line(size=0.25),panel.grid.minor = element_blank())
ggsave("hdegs_syngo.pdf")
