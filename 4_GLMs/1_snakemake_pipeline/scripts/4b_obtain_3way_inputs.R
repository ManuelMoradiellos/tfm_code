####################################################
#--{ Obtain Inputs for 3-way Interaction Models }--#
####################################################

# Using the GermlineMut-LOH table merged by pathogenic class,
# cancer type and all population group for each cancer type
# obtains the 3-way GLMs input by keeping those genes in which
# there is, at least, one account of GermlineMut or LOH event.
#
# Returns .RDS with a list of GLM inputs per gene and a 
# table with the same information in a wide format that also
# includes various columns with the 'diff'/'effect size' per
# available population group and a few extra comparison with
# the baseline group (as requested by Solip)

library(tidyverse)

in_merged_matrix_by_patho_ct <- snakemake@input[['merged_matrix_by_patho_ct']]
  
out_3way_model_input_rds <- snakemake@output[['threeway_model_input_rds']]
out_3way_model_table <- snakemake@output[['threeway_model_table']]
  
metadata_table <- snakemake@params[['metadata_table']]
current_cancer <- snakemake@params[['cancer_type']]
add_pseudocounts <- snakemake@params[['add_pseudocounts']]


# |- Functions -| #
binary_matrix_to_freq_vector <- function(sub_binary_matrix,
                                         add_median = FALSE){
  ## Obtains mean frequency vector of either Mutation 
  ## or LOH across genes of the current cancer type
  
  freq_vector <- colMeans(sub_binary_matrix) * 100 # In % format
  names(freq_vector) <- gsub( x = names(freq_vector), pattern = '*\\.germline_mut|*\\.loh_event', replacement = "")
  
  if (add_median){
    # After event freq. per gene across cancer adds median
    freq_vector['median_freq'] <- median(freq_vector) # Maybe this is not so useful for the inputs
  }
  return(freq_vector)
}

obtain_event_freq_vectors <- function(complete_binary_matrix){
  ## Obtains vectors with event freq. from binary Mut-LOH matrix
  
  # Transforms the input matrix into two 'cancer_type_pancancer x gene'
  # ones to check the frequencies of Mutation and LOH events separately
  mut_subset <- complete_binary_matrix[, c(grepl(pattern = '*\\.germline_mut', x = colnames(complete_binary_matrix)) )]
  loh_subset <- complete_binary_matrix[, c(grepl(pattern = '*\\.loh_event', x = colnames(complete_binary_matrix)) )]
  
  # Obtain gene-event frequency vectors
  mut_freq_vector <- binary_matrix_to_freq_vector(sub_binary_matrix = mut_subset, add_median = FALSE)
  loh_freq_vector <- binary_matrix_to_freq_vector(sub_binary_matrix = loh_subset, add_median = FALSE)
  
  return(list(mut_freqs = mut_freq_vector, loh_freqs = loh_freq_vector))
}

# Data Imports
merged_matrix_patho_ct <- read_tsv(in_merged_matrix_by_patho_ct)
patient_metadata <- read_tsv(metadata_table)

#  ¯¯¯¯¯¯¯¯¯¯¯¯¯¯  #
# |- Main Tasks -| #
#  --------------  #

# Takes complete GermlineMut-LOH matrix and filters
# it just to the patients with the current cancer type
if (current_cancer != 'Pancancer'){
  cancer_patients_with_pop <- patient_metadata %>% filter(cancer_type == current_cancer) %>% select(bcr_patient_barcode, mixed.PAM.Cluster)
  merged_matrix_patho_ct <- merged_matrix_patho_ct %>% filter(Patient.Barcode %in% cancer_patients_with_pop$bcr_patient_barcode)
  merged_matrix_patho_ct <- merge(cancer_patients_with_pop, merged_matrix_patho_ct, by.x = 'bcr_patient_barcode', by.y = 'Patient.Barcode')
} else{ # Unfiltered for Pancancer
  cancer_patients_with_pop <- patient_metadata %>% select(bcr_patient_barcode, mixed.PAM.Cluster)
  merged_matrix_patho_ct <- merge(cancer_patients_with_pop, merged_matrix_patho_ct, by.x = 'bcr_patient_barcode', by.y = 'Patient.Barcode')
}

# |- Obtaining Event Frequencies -| #
#
# Takes binary matrix and obtains Germline Mutation and LOH freqs. (in %)
# to use include in the final 3-way GLM's outputs or to merge them
# with the 2-way ones
freq_vectors_list <- obtain_event_freq_vectors(merged_matrix_patho_ct[ , -c(1, 2)])
mut_loh_freqs_columns_complete <- data.frame(Gene = names(freq_vectors_list$mut_freqs),
                                             Freq_Germline = freq_vectors_list$mut_freqs,
                                             Freq_LOH = freq_vectors_list$loh_freqs)

# Removes genes in which there are no counts for any of the events
mut_freqs_nonzero <- freq_vectors_list$mut_freqs[freq_vectors_list$mut_freqs > 0]
loh_freqs_nonzero <- freq_vectors_list$loh_freqs[freq_vectors_list$loh_freqs > 0]
genes_to_keep <- intersect(names(mut_freqs_nonzero), names(loh_freqs_nonzero))

# |- List of 3-way GLM Inputs -| #
#
# Stores each contingency table per gene as different elements of a list
# to be used by the GLM in future steps as is the expected format
input_list_collapsed <- setNames(lapply(genes_to_keep, function(x) {
  germline_ID <- paste0(x, '.germline_mut') # To search each gene in the matrix
  LOH_ID <- paste0(x, '.loh_event')
  out <- as.data.frame(table(factor(merged_matrix_patho_ct[[germline_ID]], levels = c(0,1)),
                             factor(merged_matrix_patho_ct[[LOH_ID]], levels = c(0,1)),
                             factor(merged_matrix_patho_ct[["mixed.PAM.Cluster"]], levels = c("EUR", "AFR", "AMR", "ASIAN")  ))) # The Baseline Population 
                                                                                                                                 # can be changed by reordering 
  colnames(out) <-  c('germline_mut', 'loh', "pop", 'N') 
  if (add_pseudocounts) {out$N <- out$N + 1} # To avoid having many zeroes, as it returns bigger standard errors in the GLMs
  return(out)
}), nm = genes_to_keep)
saveRDS(input_list_collapsed, file = out_3way_model_input_rds) 

# |- Contingency Tables/Count Tables for 3-way GLMs -| #
#
# Goes through each gene of the target list and returns contingency table
# results as a one-dimensional vector that will be used as an easy way
# to check model results
table_rows <- lapply(names(input_list_collapsed), function(x){
  return( c(x, as.numeric(as.vector(unlist(input_list_collapsed[[x]]['N'])) )))
  })
contigency_table_counts <- as.data.frame(do.call(rbind, table_rows)) # Combines all rows into a df

# Goes through all available population groups of the cancer type (some of them
# aren't present in various populations) and obtains better column names
column_name_comb <- do.call(paste, c(input_list_collapsed[[1]][c('pop', 'germline_mut','loh')], sep = '-'))
column_name_comb <- gsub('-0-0', '_noMut_noLOH',
               gsub('-1-0', '_Mut_noLOH',
               gsub('-0-1', '_noMut_LOH',
               gsub('-1-1', '_Mut_LOH', column_name_comb)))) # Translates all factor codes into comprehensible versions 
colnames(contigency_table_counts) <- c('Gene', column_name_comb)

# Computes all of Solip's recommended values
# These are different 'diff' values for all population groups, these are like proportions of the counts
for (pop_group in unique(cancer_patients_with_pop$mixed.PAM.Cluster[!is.na(cancer_patients_with_pop$mixed.PAM.Cluster)])) {
  ML <- as.numeric( contigency_table_counts[ , paste0(pop_group, '_Mut_LOH')] )
  nML <- as.numeric( contigency_table_counts[ , paste0(pop_group, '_noMut_LOH')] )
  MnL <- as.numeric( contigency_table_counts[ , paste0(pop_group, '_Mut_noLOH')] )
  nMnL <- as.numeric(contigency_table_counts[ , paste0(pop_group, '_noMut_noLOH')] )
  contigency_table_counts[, paste0(pop_group, '_diff')] <- c(ML / (ML + nML)) - (MnL / (MnL + nMnL))
}

# More comparisons that are used to quickly check along 
# the coefficients and estimates returned by the models
contigency_table_counts[, 'Min_diff'] <- sapply(1:nrow(contigency_table_counts), function(n_row){
  min(contigency_table_counts[n_row, grepl(x = colnames(contigency_table_counts), pattern = '*_diff')])
})
contigency_table_counts[, 'Max_diff'] <- sapply(1:nrow(contigency_table_counts), function(n_row){
  max(contigency_table_counts[n_row, grepl(x = colnames(contigency_table_counts), pattern = '*_diff')])
})
contigency_table_counts[, 'Min_EUR'] <- sapply(1:nrow(contigency_table_counts), function(n_row){
  contigency_table_counts[n_row, 'EUR_diff'] - contigency_table_counts[n_row, 'Min_diff'] 
})
contigency_table_counts[, 'Max_EUR'] <- sapply(1:nrow(contigency_table_counts), function(n_row){
  contigency_table_counts[n_row, 'EUR_diff'] - contigency_table_counts[n_row, 'Max_diff'] 
})
write_tsv( contigency_table_counts, file = out_3way_model_table)