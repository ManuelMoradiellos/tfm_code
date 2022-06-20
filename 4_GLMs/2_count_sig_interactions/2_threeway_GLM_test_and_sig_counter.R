##############################################
#--{ Perform 3-way GLMs between two genes }--#
##############################################

# This script does two main things:

# 1) Serves as testing grounds on the 3-way GLM code #    to add to the snakemake
#    pipeline (obtain all combinations from 2-way sig GeneA and somatic GeneB)
#
# 2) Count significant higher-order interactions from the results

library(tidyverse)
library(stats)

# |{ Functions }| #

threway_glm_contingency_inputter <- function(patho_class = 'clinvar', what_to_do = 'read',
                                 dir_path_to_files_to_read = '/local/mmoradiellos/work_stuff/6_3way_two_genes_models/3way_EURonly_across_classes/',
                                 path_to_rds = '/local/mmoradiellos/work_stuff/6_3way_two_genes_models/3way_outs/new_input_tables/3way_geneA_geneB_'){
  # Script that is used to take 3-way info files (based on my
  # results) and turn them into a list of inputs for the GLMs.
  #
  # patho_class : {allclass, clinvar, delmis, plof} Pathogenic class to focus on
  # what_to_do : {write, read} Depending on what you want to do, create the list
  #              of inputs and write them into .RDS files or read previously
  #              created inputs.
  # dir_path_to_files_to_read : The path to the directory that contains the 33 
  #                             tables with data to create the inputs
  # path_to_rds : Where are the list of inputs in .RDS either located or 
  #               wanted to be stored
  if (what_to_do == 'write') {
    # Obtains a vector with all the filenames containing data for the 3-way GLM inputs
    ##files_to_read <- list.files(path = dir_path_to_files_to_read, pattern = paste0('^', patho_class, '_'), full.names = T) # Missing 'Pancancer' ## For older tables
    files_to_read <- list.files(path = paste0(dir_path_to_files_to_read, patho_class), pattern = paste0('^', patho_class, '_'), full.names = T )
    # Regexp to extract cancer type within filename, useful if we just want to include some cancer types as it detects them automatically
    cancer_names <- gsub(paste0('.*[/]', patho_class, '_([A-Za-z]+)_gam.*[.]txt$'), '\\1', files_to_read)
    threeway_inputs <- setNames(lapply(files_to_read,function(x) {
      read_tsv(x)
    }), cancer_names) # Creates list of GLM outputs
    
    # Obtain 3-way GLMs inputs by various transformations to the previous tables
    # Returns inputs into a list of list in which the levels are Cancer_Type -> Gene-Pair 
    threeway_glm_inputs <- setNames(lapply(cancer_names, function(c_t){
      uniq_pairs <- threeway_inputs[[c_t]] %>% pull(Gene) %>% unique() # List of Unique Pairs to extract contingency tables
      
      input_tables <- setNames(lapply(uniq_pairs, function(gene_pair) {
        pair_rows <- filter(threeway_inputs[[c_t]], Gene == gene_pair) # Find both rows for each gene pair
        table_row <- unlist( c(pair_rows[1, ], pair_rows[2, 8:11]) ) # "Pastes" relevant columns from Somatic Mutation Yes or NO from each gene pair into a single vector
        names(table_row) <- unlist(c(colnames(pair_rows)[1:7], paste0(colnames(pair_rows)[8:11], "_target_yes"), paste0(colnames(pair_rows)[8:11], "_target_no")))
        contingency_table <- data.frame("Germline_A" = c(1, 0, 1, 0, 1, 0, 1, 0),
                                       "LOH_A" = c(1, 1, 0, 0, 1, 1, 0, 0),
                                       "Somatic_B" = c(1, 1, 1, 1, 0, 0, 0, 0),
                                       "N" = table_row[8: length(table_row)]) # Manually check expected combination that must be stored
        contingency_table <- contingency_table[order(contingency_table$Germline_A, contingency_table$LOH_A, contingency_table$Somatic_B) ,] 
        contingency_table$N <- as.numeric(contingency_table$N) + 1 # Add pseudocounts
        return(contingency_table)
      }), uniq_pairs) 
      return(input_tables)
    }), cancer_names) #cancer_names
    saveRDS(threeway_glm_inputs, file = paste0(path_to_rds, patho_class, '.RDS'))
    return(threeway_glm_inputs)
  } else if (what_to_do == 'read'){
    return( readRDS(file = paste0(path_to_rds, patho_class, '.RDS')) )
  }
}
threeway_contingency_glm_inputs <- threway_glm_contingency_inputter(patho_class = 'clinvar', what_to_do = 'write') # Calling the fucntion to create our inputs (Contingency GLMs)


threeway_glm_table_inputter <- function(patho_class = 'clinvar', what_to_do = 'write',
                                    dir_path_to_files_to_read = '/local/mmoradiellos/work_stuff/6_3way_two_genes_models/3way_EURonly_across_classes/',
                                    path_to_rds = '/local/mmoradiellos/work_stuff/6_3way_two_genes_models/3way_outs/new_input_tables/3way_geneA_geneB_'){
  # Script that is used to take 3-way info files (based on my
  # results) and turn them into a list of inputs for the GLMs.
  #
  # patho_class : {allclass, clinvar, delmis, plof} Pathogenic class to focus on
  # what_to_do : {write, read} Depending on what you want to do, create the list
  #              of inputs and write them into .RDS files or read previously
  #              created inputs.
  # dir_path_to_files_to_read : The path to the directory that contains the 33 
  #                             tables with data to create the inputs
  # path_to_rds : Where are the list of inputs in .RDS either located or 
  #               wanted to be stored
  if (what_to_do == 'write') {
    # Obtains a vector with all the filenames containing data for the 3-way GLM inputs
    ##files_to_read <- list.files(path = dir_path_to_files_to_read, pattern = paste0('^', patho_class, '_'), full.names = T) # Missing 'Pancancer' ## For older tables
    files_to_read <- list.files(path = paste0(dir_path_to_files_to_read, patho_class), pattern = paste0('^', patho_class, '_'), full.names = T )
    # Regexp to extract cancer type within filename, useful if we just want to include some cancer types as it detects them automatically
    cancer_names <- gsub(paste0('.*[/]', patho_class, '_([A-Za-z]+)_gam.*[.]txt$'), '\\1', files_to_read)
    threeway_inputs <- setNames(lapply(files_to_read,function(x) {
      read_tsv(x)
    }), cancer_names) # Creates list of GLM outputs
    
    # Obtain 3-way GLMs inputs by various transformations to the previous tables
    # Returns inputs into a list of list in which the levels are Cancer_Type -> Gene-Pair 
    complete_tables <- setNames(lapply( names(threeway_inputs), function(cancer_type){
      
      uniq_pairs <- threeway_inputs[[cancer_type]] %>% pull(Gene) %>% unique() # List of Unique Pairs to extract contingency tables
      
      if ( length(uniq_pairs) == 0 ) {
        return( data.frame("Gene" = 'none', "Pop" = 'none', "Size" = 0, "Freq_GermlineA" = 0, 
                           "Freq_LOHA" = 0, "Freq_SomaticB" = 0, "Germline_LOH_NoSomaticB" = 0,
                           "NoGermline_LOH_NoSomaticB" = 0, "Germline_NoLOH_NoSomaticB" = 0,
                           "NoGermline_NoLOH_NoSomaticB" = 0, "Germline_LOH_SomaticB" = 0,
                           "NoGermline_LOH_SomaticB" = 0, "Germline_NoLOH_SomaticB"  = 0,
                           "NoGermline_NoLOH_SomaticB" = 0) )
      } # Empty dataframe with the same format
      
      input_tables <- lapply(uniq_pairs, function(gene_pair){
        pair_rows <- filter(threeway_inputs[[cancer_type]], Gene == gene_pair)
        table_row <- unlist( c(pair_rows[1,], pair_rows[2,8:11]) ) # Obtain both rows for the same genetic pair
        names(table_row)[8:11] <- paste0(names(table_row)[8:11], '_SomaticB') # Add tag to distinguish data
        names(table_row)[12:15] <- paste0(names(table_row)[12:15], '_NoSomaticB')
        table_row <- table_row[-7] # Remove 'target' column
      }) 
      final_df <- as.data.frame(do.call(rbind, input_tables)) # Combine rows into table
      final_df[, 3:14] <- as.numeric(unlist(final_df[, 3:14])) # Adequate type
      final_df[, 7:14] <- final_df[, 7:14] + 1 # Pseudocounts as included for the GLMs
      return(final_df)
    }), names(threeway_inputs))
    saveRDS(complete_tables, file = paste0(path_to_rds, patho_class, '_tablelike.RDS'))
    return(complete_tables)
  } else if (what_to_do == 'read'){
    return( readRDS(file = paste0(path_to_rds, patho_class, '_tablelike.RDS')) )
  }
}


# {These two functions were lovingly stolen from:
# from https://stackoverflow.com/questions/26539441/remove-null-elements-from-list-of-lists}
## First, a function that tests whether an obkect is either NULL or a list of NULLS
is.NullOb <- function(x) is.null(x) | all(sapply(x, is.null))

rmNullObs <- function(x) {
  # Recursively step down into list, removing all such objects 
  x <- Filter(Negate(is.NullOb), x)
  lapply(x, function(x) if (is.list(x)) rmNullObs(x) else x)
}

extract_gene_pairs_to_test <- function(gene_sets_to_test = c('all_genes', 'cpgs_plus_others'),
                                  patho_class = 'clinvar', 
                                  fdr_cut = 0.2, germ_freq_cut = 0.5, somatic_B_freq_cut = 2,
                                  path_to_sig_names_matrix = '/local/mmoradiellos/work_stuff/5_gam_loh_merges_model_inputs/bck_snakemake_gamloh_models/out_all_genes_test_pop_size_corrected/merged_outputs_correct_pop_size/',
                                  dir_path_to_threeway_tables = '/local/mmoradiellos/work_stuff/6_3way_two_genes_models/3way_EURonly_across_classes/'
                                 ){
  
  # Extracts the names of those genes that showed significant association
  # in the 2-way models into a list of cancer_types -> genes_to_test 
  # for the 3-way GLMs between two genes.
  # Genes to test include the significant genes and all somatic mutation
  # paited genes (GeneB) to be test
  
  # From the total number of significant genes found across my tests,
  # we've decided to focus on the following sets (could be changed) 
  # genes_sets_to_test : 'all_genes' & 'cpgs_plus_others'
  # patho_class : Which pathogenic class to focus on (extract data from)
  #
  # Also we've decided to check the significant associations found 
  # at 20% FDR (0.2 due to scale) and 0.5% Germline Mutation Frequency for
  # ALLclass and DelMis, 0.1% for CliVar and pLoF
  #
  # fdr_cut and germ_freq_cut respectively 
  # somatic_B_freq_cut : Threshold for Somatic Mutation Frequency for Gene B
  # path_to_sig_names_matrix : Directory where the matrix with the symbols
  #                            of the significant genes is found
  # dir_path_to_files_to_read : The path to the directory that contains the 33 
  #                             tables with data to create the inputs
  
  # Extract path of significant matrix with the names of the genes to test
  sig_genes_tables_path <- sapply(gene_sets_to_test, function(gene_set){
    list.files(path =paste0(path_to_sig_names_matrix, patho_class, '/',
                            patho_class, '_glm_results/', gene_set, '/sig_counts'),
               pattern = 'gene_names', full.names = T)
  })
  # Reads matrices and extracts just the columns of the desired cutoffs
  sig_tables <- setNames(lapply(sig_genes_tables_path,function(x) {
    tab <- read_tsv(x)
    tab[, c('Cancer_type', paste0('FDR_', fdr_cut, '.GermFreq_', germ_freq_cut))]
  }), gene_sets_to_test) 
  
  ## Extracts 3-way data tables to find all matching GeneB
  ## that are to be tested with the significant GeneA from prev. matrices
  ## Obtains a vector with all the filenames containing data for the 3-way GLM inputs
  threway_tables_to_read <- files_to_read <- list.files(path = paste0('/local/mmoradiellos/work_stuff/6_3way_two_genes_models/3way_EURonly_across_classes/', patho_class), pattern = paste0('^', patho_class, '_'), full.names = T )
  # Regexp to extract cancer type within filename, useful if we just want to include some cancer types as it detects them automatically
  cancer_names <- gsub(paste0('.*[/]', patho_class, '_([A-Za-z]+)_gam.*[.]txt$'), '\\1', threway_tables_to_read)
  threeway_inputs <- setNames(lapply(threway_tables_to_read,function(x) {
    read_tsv(x)
  }), cancer_names) # Creates list of tables per cancer type
  
  # Takes previous column and turns it into a list of cancer type and genes found
  to_test_possible_null <- setNames(lapply(gene_sets_to_test, function(gene_set){
    setNames(lapply(cancer_names, function(c_t){
      # If we have no sig. genes for a cancer type return NULL
      if ( is.na(sig_tables[[gene_set]] %>% filter(Cancer_type == c_t) %>% pull( paste0('FDR_', fdr_cut, '.GermFreq_', germ_freq_cut) ) ) ){
        return()
      } else { # Return vector of all gene symbols found
        return( strsplit(x = (sig_tables[[gene_set]] %>% filter(Cancer_type == c_t) %>% pull( paste0('FDR_', fdr_cut, '.GermFreq_', germ_freq_cut) ) ), split = ',') )
      }
    }), cancer_names)
  }), gene_sets_to_test)
  
  # Previous list contains NULL entries, removes them to keep the ones to test
  to_test <- rmNullObs(to_test_possible_null) 
  
  # Create a list of gene_set -> cancer_type -> vector of genes to test
  pairs_to_test <- setNames(lapply(gene_sets_to_test, function(gene_set){
    test_ct <- setNames( lapply( names(to_test[[gene_set]]), function(cancer_types){
      gene_pairs_found <- setNames(lapply(to_test[[gene_set]][[cancer_types]][[1]], function(gene_to_test){
        
        threeway_inputs[[cancer_types]][grepl(pattern = paste0(gene_to_test, '_'), x = threeway_inputs[[cancer_types]]$Gene) , ] %>% 
          filter(Freq_SomaticB > somatic_B_freq_cut) %>% 
          pull(Gene) %>% unique() # As rows are duplicated
        
      }), to_test[[gene_set]][[cancer_types]][[1]])
      return( unlist(gene_pairs_found) ) # Returns all gene-pairs as vector for easier filtering of GLM inputs
    }), names(to_test[[gene_set]]))
  }), gene_sets_to_test)
  
  pairs_to_test_nonNA <- rmNullObs(pairs_to_test)
  
  return(pairs_to_test)
}


##################################################
#--{ Code for Genetic Interactions 3-way GLMs }--#
##################################################

pathogenic_class <- 'allclass'

# Read the complete GLM inputs for a certain pathogenic class and obtain
# our list of 2-way sig. genes to be test in higher order interaction (GeneA_GeneB)
threeway_contingency_glm_inputs <- threway_glm_contingency_inputter(patho_class = pathogenic_class, what_to_do = 'read')
threeway_table_glm_inputs <- threeway_glm_table_inputter(patho_class = pathogenic_class, what_to_do = 'read')
pairs_to_test <- extract_gene_pairs_to_test(patho_class = pathogenic_class, somatic_B_freq_cut = 2)

# To check gene pairs found vs. the ones that appear in the matrices
# do.call(c, (unique(lapply(str_split(pairs_to_test$all_genes$LGG, pattern = '_'), `[[`, 1))))

# Filter complete 3-way GLM inputs to just keep the necessary ones
inputs_to_test_possempty <- setNames(lapply( names(pairs_to_test) , function(gene_set){
  setNames( lapply( names(pairs_to_test[[gene_set]]), function(cancer_type){
    threeway_contingency_glm_inputs[[cancer_type]][ pairs_to_test[[gene_set]][[cancer_type]] ]
  }), names(pairs_to_test[[gene_set]]) )
}), names(pairs_to_test) )

lengths(inputs_to_test_possempty$all_genes)
lengths(inputs_to_test_possempty$cpgs_plus_others)

# Table to check the number of gene_pairs to test, check in case some are empty (one sig. gene isn't found on Somatic Table)
write_tsv(rbind( data.frame('Cancer_set' = 'all_genes','Cancer_type' = names(inputs_to_test_possempty$all_genes), 'Total' = lengths(inputs_to_test_possempty$all_genes)),
                 data.frame('Cancer_set' = 'cpgs_plus_others', 'Cancer_type' = names(inputs_to_test_possempty$cpgs_plus_others), 'Total' = lengths(inputs_to_test_possempty$cpgs_plus_others)) ),
          file = paste0('/local/mmoradiellos/work_stuff/6_3way_two_genes_models/3way_outs/new_input_tables/', pathogenic_class, '/3way_GLM_out_geneA_geneB_', pathogenic_class, '_tested_pairs_count.tsv') )

inputs_to_test <- lapply(inputs_to_test_possempty, function(x) Filter(length, x)) # Removes empty cases to avoid future errors 
lengths(inputs_to_test$all_genes)
lengths(inputs_to_test$cpgs_plus_others)

# Formula for the models
## Equivalent to "N ~ Germline_A + LOH_A + Somatic_B + Germline_A:LOH_A + Germline_A:Somatic_B + LOH_A:Somatic_B + Germline_A:LOH_A:Somatic_B"
formula_3way_two_genes <- "N ~ Germline_A * LOH_A * Somatic_B" 

# Perform GLMs for all relevant pairs and store results into an .RDS
glm_outputs_3way <- setNames(lapply( names(inputs_to_test), function(gene_set){
  setNames(lapply( names(inputs_to_test[[gene_set]]), function(cancer_type){
    setNames(lapply( names(inputs_to_test[[gene_set]][[cancer_type]]), function(gene_pair){
      
      glm_3way <- glm(formula = formula_3way_two_genes,
                      data = inputs_to_test[[gene_set]][[cancer_type]][[gene_pair]],
                      family = poisson(link = 'log'),
                      control = glm.control(epsilon = 1e-6, maxit=100))
      return(glm_3way)
    }), names(inputs_to_test[[gene_set]][[cancer_type]]) )
  }), names(inputs_to_test[[gene_set]]) )
}), names(inputs_to_test) )

saveRDS(glm_outputs_3way, file = paste0('/local/mmoradiellos/work_stuff/6_3way_two_genes_models/3way_outs/new_input_tables/', pathogenic_class, '/3way_GLM_out_geneA_geneB_', pathogenic_class, '_output.RDS'))

# Obtains coefficients return by the GLM in a row-like manner
glm_tab_rows_3way <- setNames(lapply( names(inputs_to_test), function(gene_set){
  setNames(lapply( names(inputs_to_test[[gene_set]]), function(cancer_type){
    setNames(lapply( names(inputs_to_test[[gene_set]][[cancer_type]]), function(gene_pair){
      
      glm_3way <- glm(formula = formula_3way_two_genes,
                      data = inputs_to_test[[gene_set]][[cancer_type]][[gene_pair]],
                      family = poisson(link = 'log'),
                      control = glm.control(epsilon = 1e-6, maxit=100))
      
      vec <- as.vector(t(summary(glm_3way)$coefficients)) # Each column as row
      names(vec) <- paste0(rep(gsub(pattern = '[()]', replacement = '', x = rownames(summary(glm_3way)$coefficients)), each = 4), c('_estimate', '_std_error', '_z_value', '_p_value'))
      return(vec)
    }), names(inputs_to_test[[gene_set]][[cancer_type]]) )
  }), names(inputs_to_test[[gene_set]]) )
}), names(inputs_to_test) )


######################################################
#--{ Count Significant Genetic Interactions 3-way }--#
######################################################

# Takes the previous rows and collapses them into data frames
# so the results can be a list of gene_set -> cancer_types -> table of all tested pairs
fdr_corr <- 'cancer_type' # cancer_type or gene_a

if (fdr_corr == 'gene_a') {
  threeway_glm_results_table <- setNames( lapply( names(glm_tab_rows_3way), function(gene_set){
    setNames( lapply( names(glm_tab_rows_3way[[gene_set]]), function(cancer_type){
      genesA <- paste0(unique(lapply(str_split(names(glm_tab_rows_3way[[gene_set]][[cancer_type]]), pattern = '_'), `[[`, 1)), '_') # Extract which are the GeneA tested at each point 
      genesA_fdr_rows <- lapply(genesA, function(curr_geneA){
        temp_list <- glm_tab_rows_3way[[gene_set]][[cancer_type]][grepl(names(glm_tab_rows_3way[[gene_set]][[cancer_type]]), pattern = curr_geneA)] # Filter elements corresponding to that GeneA
        temp_df <- do.call(rbind, temp_list)
        df <- data.frame( cbind('Gene_pairs' = rownames(temp_df), temp_df) )
        df$Germline_A.LOH_A.Somatic_B_fdr <- p.adjust(df$Germline_A.LOH_A.Somatic_B_p_value, method = "BH") # FDR Correction at tested GeneA level 
        return( df )
      })
      pre_final_df <- do.call(rbind, genesA_fdr_rows)
      return(  merge(threeway_table_glm_inputs[[cancer_type]], pre_final_df, by.x = 'Gene', by.y = 'Gene_pairs')  ) 
    }), names(glm_tab_rows_3way[[gene_set]]) )
  }) , names(glm_tab_rows_3way) )
  saveRDS(threeway_glm_results_table, file = paste0('../../../6_3way_two_genes_models/3way_outs/new_input_tables/', pathogenic_class, '/3way_GLM_geneA_geneB_', pathogenic_class, 'output_tables.RDS'))
} else if (fdr_corr == 'cancer_type'){
  threeway_glm_results_table <- setNames(lapply( names(glm_tab_rows_3way), function(gene_set){
    setNames(lapply( names(glm_tab_rows_3way[[gene_set]]), function(cancer_type){
      temp_df <- do.call(rbind, glm_tab_rows_3way[[gene_set]][[cancer_type]] )
      df <- data.frame( cbind('Gene_pairs' = rownames(temp_df), temp_df) )
      df$Germline_A.LOH_A.Somatic_B_fdr <- p.adjust(df$Germline_A.LOH_A.Somatic_B_p_value, method = "BH")
      return( merge(threeway_table_glm_inputs[[cancer_type]], df, by.x = 'Gene', by.y = 'Gene_pairs') )
    }), names(glm_tab_rows_3way[[gene_set]]) )
  }), names(glm_tab_rows_3way) )
}

# Includes all results in one single table with columns to differentiate between gene set and cancer types
almost_there_table <- lapply( names(threeway_glm_results_table), function(gene_set){
  almost <- lapply( names(threeway_glm_results_table[[gene_set]]), function(cancer_type){
    return( data.frame( cbind('Cancer_type' = rep(cancer_type, nrow(threeway_glm_results_table[[gene_set]][[cancer_type]])), threeway_glm_results_table[[gene_set]][[cancer_type]]) ) )
  })
  there <- do.call(rbind, almost)
  return( data.frame( cbind('Gene_set' = rep(gene_set, nrow(there))), there )  )
})
final_table <- do.call(rbind, almost_there_table)
final_table %>% filter(Germline_A.LOH_A.Somatic_B_fdr < 0.4)
write_tsv(final_table, file = paste0('../../../6_3way_two_genes_models/3way_outs/new_input_tables/', pathogenic_class, '/3way_GLM_geneA_geneB_', pathogenic_class, '_whole_table.tsv'))

allclass <- read_tsv(file = paste0('../../../6_3way_two_genes_models/3way_outs/new_input_tables/', 'allclass', '/3way_GLM_geneA_geneB_', 'allclass', '_whole_table.tsv'))
allclass %>% filter(Germline_A.LOH_A.Somatic_B_fdr < 0.4) %>% pull(Gene)
