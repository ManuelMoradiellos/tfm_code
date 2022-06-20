####################################################
#--{ Obtain Inputs for 2-way Interaction Models }--#
####################################################

# Using the GermlineMut-LOH table merged by pathogenic class,
# cancer type and population group obtains the 2-way GLMs
# input by keeping those genes in which there is, at least,
# one account of GermlineMut or LOH event.
# Returns .RDS with a list of GLM inputs per gene and a 
# table with the same information plus event frequencies to
# create a better-looking final output

library(tidyverse)

merged_matrix_by_patho_ct_pop <- snakemake@input[['merged_matrix_by_patho_ct_pop']]
  
out_2way_model_input_rds <- snakemake@output[['twoway_model_input_rds']]
out_2way_model_table <- snakemake@output[['twoway_model_table']]
  
population <- snakemake@params[['population']]
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

# |- Mandatory Data Imports -|
merged_patho_ct_pop <- read_tsv(merged_matrix_by_patho_ct_pop)

#  ¯¯¯¯¯¯¯¯¯¯¯¯¯¯  #
# |- Main Tasks -| #
#  --------------  # 

# |- Obtaining Event Frequencies -| #
#
# Takes binary matrix and obtains Germline Mutation and LOH freqs. (in %)
# to use include in the final 2-way GLM's outputs
freq_vectors_list <- obtain_event_freq_vectors(merged_patho_ct_pop)
mut_loh_freqs_columns_complete <- data.frame(Gene = names(freq_vectors_list$mut_freqs),
                                             Freq_Germline = freq_vectors_list$mut_freqs,
                                             Freq_LOH = freq_vectors_list$loh_freqs)

# Removes genes in which there are no counts for any of the events
mut_freqs_nonzero <- freq_vectors_list$mut_freqs[freq_vectors_list$mut_freqs > 0]
loh_freqs_nonzero <- freq_vectors_list$loh_freqs[freq_vectors_list$loh_freqs > 0]
genes_to_keep <- intersect(names(mut_freqs_nonzero), names(loh_freqs_nonzero))

# |- List of 2-way GLM Inputs -| #
#
# Stores each contingency table per gene as different elements of a list
# to be used by the GLM in future steps as is the expected format
## Changed from dplyr pipes to various base functions and approaches as it seems to work
glm_models_per_gene_list <- setNames(lapply(genes_to_keep, function(gene) {
  germline_ID <- paste0(gene, '.germline_mut')
  LOH_ID <- paste0(gene, '.loh_event')
  contingency_table <- as.data.frame(table(factor(merged_patho_ct_pop[[germline_ID]], levels = c(0,1)),
                                           factor(merged_patho_ct_pop[[LOH_ID]], levels = c(0,1))))
  colnames(contingency_table) <-  c('germline_mut', 'loh', 'N') # Names to include in GLM's formula
  if (add_pseudocounts) {contingency_table$N <- contingency_table$N + 1} # To avoid having many zeroes, as it returns bigger standard errors in the GLMs
  contingency_table <- contingency_table[order(contingency_table$germline_mut, contingency_table$loh), ] 
  return(contingency_table)
}), nm = genes_to_keep) # setNames outside of lapply to name each element of list
saveRDS(glm_models_per_gene_list, file = out_2way_model_input_rds)

# |- Table of 2-way GLM Inputs -| #
#
# Goes through each gene of the target list and returns contingency table
# results as a one-dimensional vector that will be used as an easy way
# to check model results
glm_model_rows_list <- lapply(genes_to_keep, function(gene) {
  germline_ID <- paste0(gene, '.germline_mut')
  LOH_ID <- paste0(gene, '.loh_event')
  contingency_row <- as.vector(table(factor(merged_patho_ct_pop[[germline_ID]], levels = c(0,1)),
                                     factor(merged_patho_ct_pop[[LOH_ID]], levels = c(0,1))))
  table_row <- c(gene, population, sum(contingency_row), contingency_row) 
  return( table_row )
})
model_input <- as.data.frame(do.call(rbind, glm_model_rows_list)) # Combines all row into df
colnames(model_input) <- c('Gene', 'Pop', 'Size', 'NoGermline_NoLOH', 'Germline_NoLOH', 'NoGermline_LOH', 'Germline_LOH')
if (add_pseudocounts){
  model_input$NoGermline_NoLOH <- as.numeric(model_input$NoGermline_NoLOH) + 1
  model_input$Germline_NoLOH <- as.numeric(model_input$Germline_NoLOH) + 1
  model_input$NoGermline_LOH <- as.numeric(model_input$NoGermline_LOH) + 1
  model_input$Germline_LOH <- as.numeric(model_input$Germline_LOH) + 1
} 

model_input_final <- merge(model_input, mut_loh_freqs_columns_complete, by = 'Gene') # Adds event frequency columns

# Reorders model input to match Solip's format, avoiding population group 
# info when we are considering all populations as a whole
if (population != 'ALLpop'){
  model_input_final <- model_input_final[, c('Gene', 'Pop', 'Size', 'Freq_Germline', 'Freq_LOH', 'Germline_LOH', 'NoGermline_LOH', 'Germline_NoLOH', 'NoGermline_NoLOH')]
} else {
  model_input_final <- model_input_final[, c('Gene', 'Size', 'Freq_Germline', 'Freq_LOH', 'Germline_LOH', 'NoGermline_LOH', 'Germline_NoLOH', 'NoGermline_NoLOH')]
}
write_tsv(model_input_final, file = out_2way_model_table)