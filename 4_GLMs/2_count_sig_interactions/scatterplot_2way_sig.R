#
# Scatterplots for 2-way Sig. Counts #
#

# Takes the tables with all gene pairs that passed quality thresholds 
# and that were used to perform FDR correction; these entries are
# represented in scatterplots and significant co-occurrences are highlighted 

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(patchwork)

#
# Effect Size vs p-value ( /sig. LOH -  /sig. no loh )
#

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

# Graphics Function
scatterplot_patho <- function(class, gene_set, cancer_type = 'Pancancer'){
  parent_path <- '/local/mmoradiellos/work_stuff/5_gam_loh_merges_model_inputs/bck_snakemake_gamloh_models/out_all_genes_test_pop_size_corrected/merged_outputs_correct_pop_size/'
  path_glm_table <- list.files( path = paste0(parent_path, class), pattern = cancer_type, full.names = T)
  
  table_comp <- read_tsv(path_glm_table)
  
  loh_cut <- ifelse( median(table_comp$Freq_LOH) < 5, 5, median(table_comp$Freq_LOH) )
  
  if(class == 'clinvar' | class == 'plof'){ ## Ask Solip if this makes sense, me estoy rayando ya por todo
    germ_cut = 0.1
  } else {
    germ_cut = 0.5
  }
    
  table_sub <- table_comp %>% filter( if (gene_set == 'cpgs_plus_others') ( (Freq_LOH > loh_cut ) & ( Freq_Germline > germ_cut ) & (Gene %in% new_set$Gene) )
                                      else ( (Freq_LOH > loh_cut ) & ( Freq_Germline > germ_cut) ) ) %>% select(Gene, interaction_estimate, p_value, Freq_Germline)
  
  
  if (gene_set == 'all_genes') {
    
    # Might include other possible genes that aren't in any list previously obtained
    table_sub$gene_type <- unlist(lapply(table_sub$Gene, function(gene){
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
    
    shapes <- c(21, 22, 23, 24)
    names(shapes) <- c('Known CPG', 'CPG-like characteristics', 'Known somatic driver', 'Other')
    
  } else if(gene_set == 'cpgs_plus_others'){
    
    # Can only be either CPGs or CPG-like 
    table_sub$gene_type <- unlist(lapply(table_sub$Gene, function(gene){
      if (gene %in% CPGs_list_available) {
        return( factor('Known CPG', levels = c('Known CPG', 'CPG-like characteristics')) )
      } else if (gene %in% cpg_like) {
        return( factor('CPG-like characteristics', levels = c('Known CPG', 'CPG-like characteristics')) )
      }
    }))
    
    shapes <- c(21, 22)
    names(shapes) <- c('Known CPG', 'CPG-like characteristics')
  }

  table_sub$log_p_value <- -log10(table_sub$p_value + 0.000001 )
  table_sub$twoway_fdr <- p.adjust(table_sub$p_value, method = 'BH')
  table_sub$log_10_twoway_fdr <- -log10(table_sub$twoway_fdr + 0.000001 )
  
  scat <- ggplot(table_sub, aes(x = log_p_value, y = interaction_estimate, fill = log_p_value, color = log_p_value, shape = gene_type, size = Freq_Germline)) +
    geom_point(alpha = 0.8) +
    scale_shape_manual(values = shapes) + 
    scale_fill_gradient(low = 'skyblue', high = '#BF211E', limits = c(0, 6)) +
    scale_color_gradient(low = 'skyblue', high = '#BF211E', limits = c(0, 6)) +
    scale_size_continuous( limits = c(0, 5), range = c(2, 7) ) +
    ggrepel::geom_text_repel(aes(label = ifelse(twoway_fdr < 0.2 ,as.character(Gene),'')), size = 6, position = 'identity', box.padding = 0.5)  +
    ggtitle(paste(class, gene_set, cancer_type)) +
    #ylab('Effect Size (enrichment of pathogenic variants in LOH samples)') +
    labs(size = 'Germline Variant \n Frequency',
         shape = 'Gene Type') +
    xlab( expression(~-Log[10]~(italic(P)-value)) ) +
    theme(plot.title = element_text(hjust = 0.5, size = 15)) +
    theme_light() 
  
  return(scat)
}

patho_class <- c('clinvar', 'plof', 'delmis', 'allclass')
gene_sets_to_test <- c('all_genes', 'cpgs_plus_others')  
gene_set_test <- gene_sets_to_test[2]
plot_list <- lapply(patho_class, function(p){scatterplot_patho(p, gene_set = gene_set_test)})
plot = wrap_plots(plot_list, guides = 'collect') 
