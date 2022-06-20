#
# Bubble Plots for the 2-way GLM results
#
library(tidyverse)
library(ggplot2)
library(ggpubr)


# Gene lists
#-Obtaining the subsets that will be used for the comparison-#
# All genes present in Pelayo's initial .vcf file
all_genes <- scan('/local/mmoradiellos/work_stuff/4_germline_variant_freq_gam_matrix/genes_subsets_ALL_CCG_CPGs/genes_sorted_uniq.txt', what = "character")
all_genes_without_semmic <- unique(unlist(lapply(all_genes, function(x) strsplit(x, split = ";")[[1]]))) # Remove misleading entries

# CPGs to compare approaches, from our Gitlab; keep those present in our input file
CPGs_list <- scan(file = '/local/mmoradiellos/work_stuff/4_germline_variant_freq_gam_matrix/genes_subsets_ALL_CCG_CPGs/CPG.csv', what = 'character')
CPGs_list_available <- CPGs_list[CPGs_list %in% all_genes]

# Cancer Consensus Genes from COSMIC https://cancer.sanger.ac.uk/census (GRCh37-COSMIC v95)
ccg_list_tab <- read_tsv('/local/mmoradiellos/work_stuff/4_germline_variant_freq_gam_matrix/genes_subsets_ALL_CCG_CPGs/Census_allMon_Feb_14_11_23_46_2022.tsv')
ccg_list <- unique(filter(ccg_list_tab, is.na(Germline), Somatic == "yes")$`Gene Symbol`)
ccg_list_available <- ccg_list[ccg_list %in% all_genes]
somatic_drivers <- ccg_list_available[!(ccg_list_available %in% CPGs_list_available)]

# List of genes with cpg-like features (includes canonical CPGs, DNA-repair, chromatin, other possible drivers, etc)
new_set <- read_tsv('/local/mmoradiellos/work_stuff/5_gam_loh_merges_model_inputs/snakemake_gam+loh_models_new/input_files/custom_gene_sets/CPG_PPI_input_merging_v2.txt')
cpg_like <- new_set$Gene[!(new_set$Gene %in% CPGs_list_available)]

# Empty list to be filled and later used to compare
sig_genes_compar <- list()
sig_tissue_compar <- list()



##################################
# Runs changes from this point on #

class <- 'allclass'
gene_sets_to_test <- c('all_genes', 'cpgs_plus_others')  
gene_set <- gene_sets_to_test[2]

# Complete FDR tested tables
parent_path <- './merged_outputs_correct_pop_size/'
path_glm_table <- list.files(path = paste0(parent_path, class, '/', class, '_glm_results/', gene_set, '/fdr_complete_tables'), pattern = '*.tsv', full.names = T)

cancer_types <- if (gene_set == 'all_genes') {
  gsub(pattern = paste0('.*[_]([A-Za-z]+)_FDR0[.]2_Germ0[.]1_withFDR_', class, '[.]tsv$'), replacement = '\\1', x = path_glm_table)
} else if(gene_set == 'cpgs_plus_others'){
  gsub(pattern = paste0('.*[/]custom[_]([A-Za-z]+)_FDR0[.]2_Germ0[.]1_withFDR[.]tsv$'), replacement = '\\1', x = path_glm_table)
}

all_files <- setNames(lapply(path_glm_table, function(x) {
  read_tsv(x)
}), cancer_types) # Creates list of GLM outputs

sig_counts_complete_path <- read_tsv( file = list.files(path = paste0(parent_path, class, '/', class, '_glm_results/', gene_set, '/sig_counts'), pattern = '*gene_names*', full.names =  T))

if (class == 'clinvar' | class == 'plof') {
  sig_counts <- sig_counts_complete_path %>% select(Cancer_type, FDR_0.2.GermFreq_0.1) %>% drop_na() # In some cases clinvar and plof had more sig. results for Freq > 0.5 due to having less sample to FDR correct 
} else {
  sig_counts <- sig_counts_complete_path %>% select(Cancer_type, FDR_0.2.GermFreq_0.5) %>% drop_na()
} 

# Check which counts come from tissues with less than 200 patients, removes them of the search list
sig_counts_to_check <-
  sig_counts[unlist(lapply(sig_counts$Cancer_type, function(ct){ return(all_files[[ct]]['Size'][[1]][1] > 200)})), ]

genes_to_extract <- unique(unlist(strsplit(as.vector(sig_counts_to_check %>% pull(colnames(sig_counts_to_check)[2])), split = ',')))

tables <- setNames(lapply(cancer_types, function(c_t){
  
  if (dim(all_files[[c_t]])[1] == 0) {
    return() # Some cancer types did not have any FDR-testeable results
  } else {
    if (all_files[[c_t]]['Size'][[1]][1] < 200){
      return()
    } else {
      all_files[[c_t]]['tissue'] <- rep(c_t, nrow(all_files[[c_t]]))
      
      table <-  all_files[[c_t]] %>%
        relocate(tissue) %>%
        filter(Gene %in% genes_to_extract) %>%
        select(tissue, Gene, interaction_estimate, p_value, twoway_fdr)
      
      table['significant'] <- unlist(lapply(table %>% pull(Gene), function(gene){ 
        if ( dim(sig_counts %>% filter(Cancer_type == c_t))[1] == 0  ){
          return('no')
        } else if(gene %in% unlist(strsplit(as.vector(sig_counts %>% filter(Cancer_type == c_t) %>% pull(colnames(sig_counts)[2])), split = ','))){
          return('yes')
        } else {
          return('no')
        }
      }))
      return(table)
    }
  }
}), cancer_types)

bubble_table <- do.call(rbind, tables)
bubble_table$log_p_value <- -log10(bubble_table$p_value  + 0.000001)
bubble_table$log_twoway_fdr <- -log10(bubble_table$twoway_fdr)
bubble_table$tissue_reorder <- gsub(pattern = 'Pancancer', replacement = 'aPancancer', bubble_table$tissue)

bubbles <- ggplot(bubble_table, aes(x=Gene, y=tissue_reorder, color = significant, size = log_p_value, fill = interaction_estimate))+
  geom_point(shape=21, stroke = 1.3) + 
  scale_color_manual(values = c("grey", "green4")) +
  scale_fill_gradient2(low = "blue", high = "red") +
  theme_bw() +
  labs(color = 'FDR Significant', size = expression(~-Log[10]~(italic(P)-value)), fill = 'Interaction Estimate') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


if (gene_set == 'all_genes') {
  # Might include other possible genes that aren't in any list previously obtained
  bubble_table$gene_type <- unlist(lapply(bubble_table$Gene, function(gene){
    if (gene %in% CPGs_list_available) {
      return( factor('Known CPG', levels = c('Known CPG', 'CPG-like characteristics', 'Known somatic driver', 'Other')) )
    } else if (gene %in% cpg_like) {
      return( factor('CPG-like characteristics', levels = c('Known CPG', 'CPG-like characteristics', 'Known somatic driver', 'Other')) )
    } else if (gene %in% somatic_drivers) {
      return( factor('Known somatic driver', levels = c('Known CPG', 'CPG-like characteristics', 'Known somatic driver', 'Other')) )
    } else {
      return( factor('Other', levels = c('Known CPG', 'CPG-like characteristics', 'Known somatic driver', 'Other')) )
    }
  })) 
} else if(gene_set == 'cpgs_plus_others'){
  # Can only be either CPGs or CPG-like 
  bubble_table$gene_type <- unlist(lapply(bubble_table$Gene, function(gene){
    if (gene %in% CPGs_list_available) {
      return( factor('Known CPG', levels = c('Known CPG', 'CPG-like characteristics')) )
    } else if (gene %in% cpg_like) {
      return( factor('CPG-like characteristics', levels = c('Known CPG', 'CPG-like characteristics')) )
    }
  }))
}

gene_type_df <- distinct(bubble_table %>% select(Gene, gene_type) %>% arrange(Gene))

# Color palettes for each gene type
all_genes_colors <- c('Known CPG' = '#F8766D',
                      'CPG-like characteristics' = '#00BFC4',
                      'Known somatic driver' = '#12664F',
                      'Other' = '#7B8CDE') 
cpgs_plus_colors <- c('Known CPG' = '#F8766D', 'CPG-like characteristics' = '#00BFC4')

# Bar to indicate the type of gene according to our lists
gene_type <- ggplot(gene_type_df, aes(x = Gene, y = 1, fill = gene_type)) +
  geom_tile() +
  scale_fill_manual( values = if (gene_set == 'all_genes') all_genes_colors else (cpgs_plus_colors) ) +
  theme_bw() +
  labs(fill = 'Gene Set') +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size = 10, face = 'bold'),
        axis.title.x = element_text(size = 11),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

#ggarrange(bubbles, gene_type, heights = c(15, 3), ncol = 1, align = 'v') # To check first merge
#https://stackoverflow.com/questions/51848658/ggplot-adding-tracking-colors-below-x-axis


### Size of tested genes for each tissue
yticks_tissue <- gsub('aPancancer', 'Pancancer', sort(unique(bubble_table$tissue_reorder))) # To only show here the tissue types

tested_tissues <- unique(bubble_table$tissue)
reznor <- setNames(lapply(tested_tissues, function(c_t){
  if (all_files[[c_t]]['Size'][[1]][1] < 200){
    return()
  } else { 
    if (class == 'clinvar' | class == 'plof') {
      sub <- all_files[[c_t]] %>% filter(Freq_Germline > 0.1)
    } else{
      sub <- all_files[[c_t]] %>% filter(Freq_Germline > 0.5)
    }
    return( data.frame('tissue' = c_t, 'FDR_test_size' = nrow(sub) ))
  }
}), tested_tissues)
fdr_test_size <- do.call(rbind, reznor)  
fdr_test_size$tissue_reorder <- gsub('Pancancer', 'aPancancer', unique(fdr_test_size$tissue))

fdr_size_plot <- ggplot(fdr_test_size, aes(x = 1, y = tissue_reorder, fill = FDR_test_size)) +
  geom_tile(color = 'black', size = 0.2) +
  scale_fill_continuous(low = '#abdbe3', high = '#063970' ) +
  theme_bw() +
  labs(fill = 'Num. Genes \ntested for FDR') +
  ylab('Cancer Types') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.y = element_text(size = 11),
        axis.text.y = element_text(face = "bold", size = 9, angle = 40),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "left") +  # Like Void but does not remove x axis
  scale_y_discrete(labels = yticks_tissue) 

### Avg. LOH Frequency & Avg. Germline Frequency of tested tissues
mean_freqs <- setNames(lapply(tested_tissues, function(c_t){
  if (all_files[[c_t]]['Size'][[1]][1] < 200){
    return()
  } else { 
    return( data.frame('tissue' = c_t,
                       'avg_LOH_Freq' = mean(all_files[[c_t]]['Freq_LOH'][[1]]),
                       'avg_Germline_Freq' = mean(all_files[[c_t]]['Freq_Germline'][[1]])))
  }
}), tested_tissues)
mean_freqs_table <- do.call(rbind, mean_freqs)
mean_freqs_table$tissue_reorder <- gsub('Pancancer', 'aPancancer', sort(unique(bubble_table$tissue)))

loh_mean_freqs_plot <- ggplot(mean_freqs_table, aes(x = 1, y = tissue_reorder, fill = avg_LOH_Freq)) +
  geom_tile(color = 'black', size = 0.2) +
  scale_fill_gradient2(mid = '#AAD922', high = '#2F4B26', limits = c(0, 100)) +
  theme_bw() +
  labs(fill = 'Avg. LOH\nFreq (%)') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "left")

germ_mean_freqs_plot <- ggplot(mean_freqs_table, aes(x = 1, y = tissue_reorder, fill = avg_Germline_Freq)) +
  geom_tile(color = 'black', size = 0.2) +
  scale_fill_gradient2(low = '#F17F29', high = '#F96900', limits = c(0, 1)) +
  theme_bw() +
  labs(fill = 'Avg. Germline\nFreq (%)') +
  ylab('Cancer Types') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "left")

# Merging legends
legend_1 <- get_legend(fdr_size_plot)
legend_2 <- get_legend(loh_mean_freqs_plot)
legend_3 <- get_legend(germ_mean_freqs_plot)
side_legends <- ggarrange(legend_1, legend_2, legend_3, nrow=3, align = 'v')

# Combining plots
rm_legend <- function(p){p + theme(legend.position = "none")}
side_plots <- ggarrange(rm_legend(fdr_size_plot), rm_legend(loh_mean_freqs_plot), rm_legend(germ_mean_freqs_plot),
                        nrow = 1, widths = c(1, 0.32, 0.32))

side_plots_with_legends <- ggarrange(side_legends, side_plots, widths = c(0.4, 0.6))

print(class)
print(gene_set)
#--------------------------------------#
# |{ Final Combination of all Plots }|
#--------------------------------------#
ggarrange(
  ggarrange(side_plots_with_legends, (ggplot() + theme_void()), heights = c(15, 3), ncol = 1 , align = 'v'),
  ggarrange(bubbles, gene_type, heights = c(15, 2.8), ncol = 1, align = 'v'),
  widths = c(3.7,15)
)

 # 1818 792 for the image

# Add current results to list of comparisons
sig_genes_compar[[paste0(class, '_', gene_set)]] <- gene_type_df$Gene
sig_tissue_compar[[paste0(class, '_', gene_set)]] <- sig_counts_to_check$Cancer_type

#------------------------------#
# |{ Class-specific Results }|
#------------------------------#

# -- |This section must be done after all four pathogenic classes are done|-- #
# Unique_genes_found
sig_genes <- lapply(names(sig_genes_compar), function(class){
  expand.grid(x = class, sig_genes_compar[[class]])
  #expand.grid(gsub(x = class, '_[a-z]+', ''), sig_genes_compar[[class]]) # Just class, not gene set detail
})
class_gene <- do.call(rbind, sig_genes) 
keep_singles <- function(v){
  v[!(v %in% v[duplicated(v)])] 
}
class_gene %>% filter(Var2 %in% keep_singles(class_gene$Var2))


sig_tissues <- lapply(names(sig_tissue_compar), function(class){
  expand.grid(x = class, sig_tissue_compar[[class]])
  #expand.grid(gsub(x = class, '_[a-z]+', ''), sig_genes_compar[[class]]) # Just class, not gene set detail
})
class_tissue <- do.call(rbind, sig_tissues)
class_tissue %>% filter(Var2 %in% keep_singles(class_tissue$Var2))
