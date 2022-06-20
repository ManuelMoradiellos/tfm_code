######################################################
#--{ Count Significant Genetic Interactions 2-way }--#
######################################################

# Various functions to count significant interactions 
# from the GLMs outputs: 
#
# 1) Obtain a matrix of significant counts using different FDR and Germline Mut. Freq. cutoffs 
# 2) Obtain the proportion of significant gene interactions over the total number of genes tested FDR 
# 3) Obtain the entries for the table that are used for the FDR adjusting using the less strict conditions

library(tidyverse)
library(stats)

# |- Functions -| #

significant_interactions_counter <- function(table_list = all_files, c_types = cancer_names,
                                             gene_set = 'all_genes', cpg_list = cpgs, custom_list = '',
                                             fdr_conditions = c(0.1, 0.2),
                                             germline_freq_conditions = c(0.1, 0.5, 1),
                                             count_cpgs_in_all = FALSE, possible_novel = FALSE,
                                             save_table = TRUE,
                                             out_subdir = '',
                                             file_suffix = '_significant_counts.tsv'){
  
  ## Given a list of 2-way GLM outputs, returns a matrix of significant counts
  ## found when filtering with different conditions (Germline Freq. & FDR Cutoffs)
  ## giving the results per cancer type and combination of conditions
  #
  # gene_set : Either 'all_genes', 'cpgs' or 'custom' (Which genes to focus on)
  # cpg_list | custom_list : List of gene names to focus on
  # fdr_conditions : Set of FDR cutoffs to apply (in %)
  # germline_freq_conditions : Set Germline Freq. cutoffs to test (in %)
  # count_cpgs_in_all : Perform the FDR using all genes but focus on CPGs found
  #                     within that set to compare both approaches
  # possible_novel : Performs FDR using all available genes but focus on those
  #                  that aren't in the canonical ~150 CPGs list as they could
  #                  be novel CPG-like genes
  # save_table & out_subdir : If you want to save resulting table and where to do so
  
  counts <- setNames(lapply(c_types, function(cancer_type){ # Loops per cancer type
    fdr_count <- setNames(lapply( fdr_conditions, function(fdr_cut){         # Per FDR Cutoff
      germ_count <- setNames(lapply( germline_freq_conditions, function(germ_cut){  # Per Germline Freq. Cutoff
        
        # Some cancer types (THYM, LAML) have really low LOH freq. so using the median resulted in
        # an overestimation of significant associations. Decided on 5% as it seemed reasonable
        loh_cut <- ifelse( median(table_list[[cancer_type]]$Freq_LOH) < 5, 5, median(table_list[[cancer_type]]$Freq_LOH) )
        
        filt_table <- table_list[[cancer_type]] %>%
          filter( if (gene_set == 'cpgs') ( (Freq_LOH > loh_cut ) & ( Freq_Germline > germ_cut ) & (Gene %in% cpg_list) )
                  else if (gene_set == 'custom') ( (Freq_LOH > loh_cut ) & ( Freq_Germline > germ_cut ) & (Gene %in% custom_list) )
                  else ( (Freq_LOH > loh_cut ) & ( Freq_Germline > germ_cut) ) ) # May filter by CPGs
        filt_table$twoway_fdr <- p.adjust(filt_table$p_value, method = "BH")  # Adjusted FDRs/q-values from filtered table
        
        if ( (gene_set == 'all_genes') & count_cpgs_in_all ) { # AllGenes to apply FDR but now keep just CPGs
          filt_table <- filt_table %>%  filter(Gene %in% cpg_list) 
        } else if ( (gene_set == 'all_genes') & possible_novel ) { # Keeps significant findings in not canonical CPGs
          filt_table <- filt_table %>%  filter( !(Gene %in% cpg_list) )
        }
        # Counts Positive direction Interaction and under FDR Cutoff: Significant Associations
        return( dim(filt_table %>% filter( (interaction_estimate > 0) & (twoway_fdr < fdr_cut) ))[1] ) 
      }), c( paste0('GermFreq_', germline_freq_conditions ) ) )
    }), c( paste0('FDR_', fdr_conditions ) ) )
  }), c_types)
  
  sig_counts_df <- counts %>% map(as.data.frame) %>% bind_rows() # Creates "matrix" of all information
  rownames(sig_counts_df) <- c_types
  if (save_table) {
    if ( (gene_set == 'all_genes') & count_cpgs_in_all ) {
      testing <- 'all_genes_counting_CPGs' # Change name to reference correctly what its being stored
    }
    write_tsv(x = data.frame('Cancer_type' = rownames(sig_counts_df), sig_counts_df),
              file = paste0(dirname(files_to_read[1]), out_subdir, gene_set, file_suffix) )
  }
  return( data.frame('Cancer_type' = rownames(sig_counts_df), sig_counts_df) )
}


proportion_tested_genes <- function(table_list = all_files, sig_counts_matrix,
                                    c_types = cancer_names, gene_set = 'all_genes', 
                                    cpg_list = cpgs, custom_list = '', 
                                    fdr_conditions = c(0.1, 0.2),
                                    germline_freq_conditions = c(0.1, 0.5, 1),
                                    count_cpgs_in_all = FALSE, possible_novel = FALSE,
                                    save_table = TRUE, out_subdir = '/sig_counts_delmis/',
                                    file_suffix_proportion_count = '_proportion_sig_tested_fdr.tsv',
                                    file_suffix_tested_count = '_count_sig_tested_fdr.tsv'){
  
  ## Given a list of 2-way GLM outputs, returns a matrix of the total number of 
  ## genes used for the FDR correction (size of genes tested for significance)
  ## when filtering with different conditions (Germline Freq. & FDR Cutoffs)
  ## giving the results per cancer type and combination of conditions
  ## Returns both the total number of genes tested and the proportion of
  ## {significant_counts/tested_genes}
  #
  # gene_set : Either 'all_genes', 'cpgs' or 'custom' (Which genes to focus on)
  # sig_counts_matrix : Result of 'significant_interactions_counter' function,
  #                     needed for the proportion count
  # cpg_list | custom_list : List of gene names to focus on
  # fdr_conditions : Set of FDR cutoffs to apply (in %)
  # germline_freq_conditions : Set Germline Freq. cutoffs to test (in %)
  # count_cpgs_in_all : Perform the FDR using all genes but focus on CPGs found
  #                     within that set to compare both approaches
  # possible_novel : Performs FDR using all available genes but focus on those
  #                  that aren't in the canonical ~150 CPGs list as they could
  #                  be novel CPG-like genes
  # save_table & out_subdir : If you want to save resulting table and where to do so
  # file_suffix_* : File suffix names to the written tables
  
  counts <- setNames(lapply(c_types, function(cancer_type){ # Loops per cancer type
    fdr_count <- setNames(lapply( fdr_conditions, function(fdr_cut){         # Per FDR Cutoff
      germ_count <- setNames(lapply( germline_freq_conditions, function(germ_cut){  # Per Germline Freq. Cutoff
        
        # Some cancer types (THYM, LAML) have really low LOH freq. so using the median resulted in
        # an overestimation of significant associations. Decided on 5% as it seemed reasonable
        loh_cut <- ifelse( median(table_list[[cancer_type]]$Freq_LOH) < 5, 5, median(table_list[[cancer_type]]$Freq_LOH) )
        
        filt_table <- table_list[[cancer_type]] %>%
          filter( if (gene_set == 'cpgs') ( (Freq_LOH > loh_cut ) & ( Freq_Germline > germ_cut ) & (Gene %in% cpg_list) )
                  else if (gene_set == 'custom') ( (Freq_LOH > loh_cut ) & ( Freq_Germline > germ_cut ) & (Gene %in% custom_list) )
                  else ( (Freq_LOH > loh_cut ) & ( Freq_Germline > germ_cut) ) ) # May filter by CPGs
        filt_table$twoway_fdr <- p.adjust(filt_table$p_value, method = "BH")  # Adjusted FDRs/q-values from filtered table
        
        if ( (gene_set == 'all_genes') & count_cpgs_in_all ) { # AllGenes to apply FDR but now keep just CPGs
          filt_table <- filt_table %>%  filter(Gene %in% cpg_list) 
        } else if ( (gene_set == 'all_genes') & possible_novel ) { # Keeps significant findings in not canonical CPGs
          filt_table <- filt_table %>%  filter( !(Gene %in% cpg_list) )
        }
        return( dim(filt_table)[1] ) # Returns the number of genes used for the FDR correction
      }), c( paste0('GermFreq_', germline_freq_conditions ) ) )
    }), c( paste0('FDR_', fdr_conditions ) ) )
  }), c_types)
  
  size_counts_df <- counts %>% map(as.data.frame) %>% bind_rows()
  sig_counts_df <- sig_counts_matrix[, -1]
  if (save_table) {
    if ( (gene_set == 'all_genes') & count_cpgs_in_all ) {
      testing <- 'all_genes_counting_CPGs' # Change name to reference correctly what its being stored
    }
    write_tsv(x = data.frame('Cancer_type' = rownames(sig_counts_matrix), round(sig_counts_df/size_counts_df, 4) * 100),
              file = paste0(dirname(files_to_read[1]), out_subdir, gene_set, file_suffix_proportion_count) ) # % of Sig. Genes over all used for FDR
    write_tsv(x = data.frame('Cancer_type' = rownames(sig_counts_matrix), size_counts_df),
              file = paste0(dirname(files_to_read[1]), out_subdir, gene_set, file_suffix_tested_count) ) # Total number of Genes used for FDR
  }
  return( data.frame('Cancer_type' = rownames(sig_counts_matrix), size_counts_df) )
}


obtain_complete_fdr_tested_gene_info <- function(table_list = all_files, c_types = cancer_names,
                                                 gene_set = 'all_genes', cpg_list = cpgs, custom_list = '',
                                                 fdr_conditions = c(0.2),
                                                 germline_freq_conditions = c(0.5),
                                                 count_cpgs_in_all = FALSE, possible_novel = FALSE,
                                                 save_table = TRUE,
                                                 out_subdir = '/sig_counts_allclass/',
                                                 files_suffix = '_FDR0.2_Germ0.5_withFDR'){
  # Obtains all rows for genes that pass the less strict set
  # of conditions, those that are used for the FDR adjustment
  
  less_stric_cond_tables_fdr <- setNames(lapply(c_types, function(cancer_type){
    fdr_count <- setNames(lapply( fdr_conditions, function(fdr_cut){
      germ_count <- setNames(lapply( germline_freq_conditions, function(germ_cut){
        
        # Some cancer types (THYM, LAML) have really low LOH freq. so using the median resulted in
        # an overestimation of significant associations. Decided on 5% as it seemed reasonable
        loh_cut <- ifelse( median(table_list[[cancer_type]]$Freq_LOH) < 5, 5, median(table_list[[cancer_type]]$Freq_LOH) )
        
        filt_table <- table_list[[cancer_type]] %>%
          filter( if (gene_set == 'cpgs') ( (Freq_LOH > loh_cut) & ( Freq_Germline > germ_cut) & (Gene %in% cpgs) )
                  else if (gene_set == 'custom') ( (Freq_LOH > loh_cut ) & ( Freq_Germline > germ_cut ) & (Gene %in% custom_list) )
                  else ( (Freq_LOH > loh_cut) & ( Freq_Germline > germ_cut) ) ) # May filter by CPGs
        filt_table$twoway_fdr <- p.adjust(filt_table$p_value, method = "BH") # Adjusted FDRs/q-values from filtered table
        
        if ( (gene_set == 'all_genes') & count_cpgs_in_all ) { # AllGenes to apply FDR but now keep just CPGs
          filt_table <- filt_table %>%  filter(Gene %in% cpg_list) 
        } else if ( (gene_set == 'all_genes') & possible_novel ) { # Keeps significant findings in not canonical CPGs
          filt_table <- filt_table %>%  filter( !(Gene %in% cpg_list) )
        }
        return( filt_table %>% relocate(twoway_fdr, .after = p_value) )
      }), c( paste0('GermFreq_', germline_freq_conditions ) ) )
    }), c( paste0('FDR_', fdr_conditions ) ) )
  }), c_types)
  
  nemuri <- setNames(lapply(c_types, function(cancer_type){
    lingus <- as.data.frame(less_stric_cond_tables_fdr[cancer_type][1][1])
    colnames(lingus) <- gsub(pattern = paste0(cancer_type, '.FDR_0.2.GermFreq_0.5.'), replacement = '', x = colnames(lingus)) # Keep clear colnames, not list of lists tags
    if (save_table) {
      write_tsv(lingus, file = paste0(dirname(files_to_read[1]), out_subdir, gene_set, '_', cancer_type, files_suffix, '.tsv') )
    }
    return(lingus)
    }), c_types)
  return(nemuri)
}


shikari <- function(table_list = all_files, c_types = cancer_names,
                                             gene_set = 'all_genes', cpg_list = cpgs, custom_list = '',
                                             fdr_conditions = c(0.1, 0.2),
                                             germline_freq_conditions = c(0.1, 0.5, 1),
                                             count_cpgs_in_all = FALSE, possible_novel = FALSE,
                                             save_table = TRUE,
                                             out_subdir = '',
                                             file_suffix = '_significant_gene_names.tsv'){
  
  ## Given a list of 2-way GLM outputs, returns a matrix of the names of significant genes
  ## found when filtering with different conditions (Germline Freq. & FDR Cutoffs)
  ## giving the results per cancer type and combination of conditions
  #
  # gene_set : Either 'all_genes', 'cpgs' or 'custom' (Which genes to focus on)
  # cpg_list | custom_list : List of gene names to focus on
  # fdr_conditions : Set of FDR cutoffs to apply (in %)
  # germline_freq_conditions : Set Germline Freq. cutoffs to test (in %)
  # count_cpgs_in_all : Perform the FDR using all genes but focus on CPGs found
  #                     within that set to compare both approaches
  # possible_novel : Performs FDR using all available genes but focus on those
  #                  that aren't in the canonical ~150 CPGs list as they could
  #                  be novel CPG-like genes
  # save_table & out_subdir : If you want to save resulting table and where to do so
  
  sig_gene_names <- setNames(lapply(c_types, function(cancer_type){ # Loops per cancer type
    fdr_count <- setNames(lapply( fdr_conditions, function(fdr_cut){         # Per FDR Cutoff
      germ_count <- setNames(lapply( germline_freq_conditions, function(germ_cut){  # Per Germline Freq. Cutoff
        
        # Some cancer types (THYM, LAML) have really low LOH freq. so using the median resulted in
        # an overestimation of significant associations. Decided on 5% as it seemed reasonable
        loh_cut <- ifelse( median(table_list[[cancer_type]]$Freq_LOH) < 5, 5, median(table_list[[cancer_type]]$Freq_LOH) )
        
        filt_table <- table_list[[cancer_type]] %>%
          filter( if (gene_set == 'cpgs') ( (Freq_LOH > loh_cut ) & ( Freq_Germline > germ_cut ) & (Gene %in% cpg_list) )
                  else if (gene_set == 'custom') ( (Freq_LOH > loh_cut ) & ( Freq_Germline > germ_cut ) & (Gene %in% custom_list) )
                  else ( (Freq_LOH > loh_cut ) & ( Freq_Germline > germ_cut) ) ) # May filter by CPGs
        filt_table$twoway_fdr <- p.adjust(filt_table$p_value, method = "BH")  # Adjusted FDRs/q-values from filtered table
        
        if ( (gene_set == 'all_genes') & count_cpgs_in_all ) { # AllGenes to apply FDR but now keep just CPGs
          filt_table <- filt_table %>%  filter(Gene %in% cpg_list) 
        } else if ( (gene_set == 'all_genes') & possible_novel ) { # Keeps significant findings in not canonical CPGs
          filt_table <- filt_table %>%  filter( !(Gene %in% cpg_list) )
        }
        # Counts Positive direction Interaction and under FDR Cutoff: Significant Associations
        sig <- filt_table %>% filter( (interaction_estimate > 0) & (twoway_fdr < fdr_cut) )
        return( paste(sig$Gene, collapse = ', ') ) 
      }), c( paste0('GermFreq_', germline_freq_conditions ) ) )
    }), c( paste0('FDR_', fdr_conditions ) ) )
  }), c_types)
  
  sig_gene_names_df <- sig_gene_names %>% map(as.data.frame) %>% bind_rows() # Creates "matrix" of all information
  rownames(sig_gene_names_df) <- c_types
  if (save_table) {
    if ( (gene_set == 'all_genes') & count_cpgs_in_all ) {
      testing <- 'all_genes_counting_CPGs' # Change name to reference correctly what its being stored
    }
    write_tsv(x = data.frame('Cancer_type' = rownames(sig_gene_names_df), sig_gene_names_df),
              file = paste0(dirname(files_to_read[1]), out_subdir, gene_set, file_suffix) )
  }
  return( data.frame('Cancer_type' = rownames(sig_gene_names_df), sig_gene_names_df) )
}


# |- Data Imports -| #
files_to_read <- list.files(path = "../../bck_snakemake_gamloh_models/out_all_genes_test_pop_size_corrected/merged_outputs_corrected_pop_size/merged_outputs/allclass", pattern = "\\.tsv$", full.names = T)
# Regexp to extract cancer type within filename, useful if we just want to include some cancer types as it detects them automatically
cancer_names <- gsub(".*[_]([^.]+)[.]tsv$", "\\1", files_to_read) 
all_files <- setNames(lapply(files_to_read,function(x) {
  read_tsv(x)
}), cancer_names) # Creates list of GLM outputs

cpgs <- scan('../input_files/custom_gene_sets/CPG.csv', what = 'character') # 152 Canonical CPGs list

# New set of Genes to test, these were curated by checking
# how many features they had that were related to CPG
# characteristics. ('If its a CPG', 'Next to DNA Repair Genes', 'In the PIK_mTOR route', 
#                  'Related or Next to Cell Cycle Genes', 'Chromatin Related')  
# It contains canonical CPGs and other possible novel ones
new_set <- read_tsv('../input_files/custom_gene_sets/CPG_PPI_input_merging_v2.txt')
genes_with_cpg_features <- new_set$Gene

###############################
#--{ Script Usage/Showcase }--#
###############################

# All Genes
all_genes_sig_counts <- significant_interactions_counter(out_subdir = '/allclass_glm_results/all_genes/sig_counts/',
                                                         file_suffix = '_significant_counts_allclass.tsv')

all_genes_sig_genes <- shikari(out_subdir = '/allclass_glm_results/all_genes/sig_counts/',
                                                         file_suffix = '_significant_gene_names_allclass.tsv')

all_genes_number_of_tested_genes <- proportion_tested_genes(sig_counts_matrix = all_genes_sig_counts,
                                                            out_subdir = '/allclass_glm_results/all_genes/sig_counts/',
                                                            file_suffix_proportion_count = '_proportion_sig_tested_fdr_allclass.tsv',
                                                            file_suffix_tested_count = '_count_sig_tested_fdr_allclass.tsv')

all_genes_fdr_complete <-  obtain_complete_fdr_tested_gene_info(out_subdir = '/allclass_glm_results/all_genes/fdr_complete_tables/',
                                                                
                                                                files_suffix = '_FDR0.2_Germ0.5_withFDR_allclass')

# CPGs
cpgs_sig_counts <- significant_interactions_counter(gene_set = 'cpgs',
                                                    out_subdir = '/delmis_glm_results/cpgs/sig_counts/',
                                                    file_suffix = '_significant_counts_delmis.tsv')

cpgs_sig_genes <- shikari(gene_set = 'cpgs',
                               out_subdir = '/delmis_glm_results/cpgs/sig_counts/',
                               file_suffix = '_significant_gene_names_delmis.tsv')

cpgs_number_of_tested_genes <- proportion_tested_genes(sig_counts_matrix = all_genes_sig_counts, gene_set = 'cpgs',
                                                       out_subdir = '/delmis_glm_results/cpgs/sig_counts/',
                                                       file_suffix_proportion_count = '_proportion_sig_tested_fdr_delmis.tsv',
                                                       file_suffix_tested_count = '_count_sig_tested_fdr_delmis.tsv')

cpgs_fdr_complete <-  obtain_complete_fdr_tested_gene_info(gene_set = 'cpgs',
                                                           out_subdir = '/delmis_glm_results/cpgs/fdr_complete_tables/',
                                                           files_suffix = '_FDR0.2_Germ0.1_withFDR_delmis')

# All Genes but counting CPGs
all_genes_cpgs_sig_counts <- significant_interactions_counter(count_cpgs_in_all = TRUE,
                                                              out_subdir = '/delmis_glm_results/all_genes_counting_cpgs/sig_counts/',
                                                              file_suffix = '_significant_counts_delmis.tsv')

all_genes_cpgs_sig_genes <- shikari(count_cpgs_in_all = TRUE,
                                    out_subdir = '/delmis_glm_results/all_genes_counting_cpgs/sig_counts/',
                                    file_suffix = '_significant_gene_names_delmis.tsv')


all_genes_cpgs_number_of_tested_genes <- proportion_tested_genes(count_cpgs_in_all = TRUE, sig_counts_matrix = all_genes_cpgs_sig_counts,
                                                                 out_subdir = '/delmis_glm_results/all_genes_counting_cpgs/sig_counts/',
                                                                 file_suffix_proportion_count = '_proportion_sig_tested_fdr_delmis.tsv',
                                                                 file_suffix_tested_count = '_count_sig_tested_fdr_delmis.tsv')
                                                                 
all_genes_cpgs_fdr_complete <-  obtain_complete_fdr_tested_gene_info(count_cpgs_in_all = TRUE,
                                                                     out_subdir = '/delmis_glm_results/all_genes_counting_cpgs/fdr_complete_tables/',
                                                                     files_suffix = '_FDR0.2_Germ0.1_withFDR_delmis')

# All Genes but not CPGs
all_genes_not_cpgs_sig_counts <- significant_interactions_counter(possible_novel = TRUE,
                                                                  out_subdir = '/delmis_glm_results/all_genes_excluding_cpgs/sig_counts/',
                                                                  file_suffix = '_not_cpgs_significant_counts_delmis.tsv')

all_genes_not_cpgs_sig_genes <- shikari(possible_novel = TRUE,
                                        out_subdir = '/delmis_glm_results/all_genes_excluding_cpgs/sig_counts/',
                                        file_suffix = '_not_cpgs_significant_gene_names_delmis.tsv')

all_genes_not_cpgs_number_of_tested_genes <- proportion_tested_genes(possible_novel = TRUE, sig_counts_matrix = all_genes_not_cpgs_sig_counts,
                                                                     out_subdir = '/delmis_glm_results/all_genes_excluding_cpgs/sig_counts/',
                                                                     file_suffix_proportion_count = '_not_cpgs_proportion_sig_tested_fdr_delmis.tsv',
                                                                     file_suffix_tested_count = '_not_cpgs_count_sig_tested_fdr_delmis.tsv')

all_genes_not_cpgs_fdr_complete <-  obtain_complete_fdr_tested_gene_info(possible_novel = TRUE,
                                                                         out_subdir = '/delmis_glm_results/all_genes_excluding_cpgs/fdr_complete_tables/',
                                                                         files_suffix = '_not_cpgs_FDR0.2_Germ0.1_withFDR')

# Genes with CPG-like features(including CPGs)
cpg_like_features_sig_counts <- significant_interactions_counter(gene_set = 'custom', custom_list = genes_with_cpg_features,
                                                                 out_subdir = '/delmis_glm_results/genes_features_cpg_like/sig_counts/',
                                                                 file_suffix = '_cpg_like_features_significant_counts_delmis.tsv')

cpg_like_features_sig_genes <- shikari(gene_set = 'custom', custom_list = genes_with_cpg_features,
                                       out_subdir = '/delmis_glm_results/genes_features_cpg_like/sig_counts/',
                                       file_suffix = '_cpg_like_features_significant_gene_names_delmis.tsv')

cpg_like_features_number_of_tested_genes <- proportion_tested_genes(sig_counts_matrix = cpg_like_features_sig_counts,
                                                                    gene_set = 'custom', custom_list = genes_with_cpg_features,
                                                                    out_subdir = '/delmis_glm_results/genes_features_cpg_like/sig_counts/')

cpg_like_features_fdr_complete <-  obtain_complete_fdr_tested_gene_info(gene_set = 'custom', custom_list = genes_with_cpg_features,
                                                                out_subdir = '/delmis_glm_results/genes_features_cpg_like/fdr_complete_tables/')


c( sum(all_genes_sig_counts[, -1]), sum(cpgs_sig_counts[, -1]),
   sum(all_genes_cpgs_sig_counts[, -1]), sum(all_genes_not_cpgs_sig_counts[, -1]),
   sum(cpg_like_features_sig_counts[, -1]))


