library(tidyverse)

# |- Functions -| #
binary_matrix_to_freq_vector <- function(sub_binary_matrix,
                                         add_median = FALSE){
  ## Obtains mean frequency vector of either Mutation 
  ## or LOH across genes of the current cancer type
  
  freq_vector <- colMeans(sub_binary_matrix) * 100 # In % format
  names(freq_vector) <- gsub( x = names(freq_vector), pattern = '*\\.germline_mut|*\\.somatic_mut', replacement = "")
  
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
  germline_subset <- complete_binary_matrix[, c(grepl(pattern = '*\\.germline_mut', x = colnames(complete_binary_matrix)) )]
  somatic_subset <- complete_binary_matrix[, c(grepl(pattern = '*\\.somatic_mut', x = colnames(complete_binary_matrix)) )]
  
  # Obtain gene-event frequency vectors
  germline_freq_vector <- binary_matrix_to_freq_vector(sub_binary_matrix = germline_subset, add_median = FALSE)
  somatic_freq_vector <- binary_matrix_to_freq_vector(sub_binary_matrix = somatic_subset, add_median = FALSE)
  
  return(list(germline_freqs = germline_freq_vector, somatic_freqs = somatic_freq_vector))
}

# |- Mandatory Data Imports -|
in_merged_matrix_by_patho_ct <- snakemake@input[['merged_matrix_by_patho_ct_pop']] # Only EUR
in_somatic_mutation_matrix <- snakemake@input[['somatic_mutation_events']]

out_geneA_geneB_glm_inputs_rds <- snakemake@output[['geneA_geneB_glm_inputs_rds']]
out_geneA_geneB_glm_inputs_table <- snakemake@output[['geneA_geneB_glm_inputs_table']]

add_pseudocounts <- snakemake@params[['add_pseudocounts']]
filt_B_freq <- as.numeric(snakemake@params[['filt_B_freq']]) # 2%


#  ¯¯¯¯¯¯¯¯¯¯¯¯¯¯  #
# |- Main Tasks -| #
#  --------------  # 

merged_matrix_patho_ct <- read_tsv(in_merged_matrix_by_patho_ct)
somatic_mutation_matrix <- read_tsv(in_somatic_mutation_matrix)
germsomloh_matrix_patho_ct = merge(merged_matrix_patho_ct, somatic_mutation_matrix, by = "Patient.Barcode")


# |- Obtaining Event Frequencies -| #
#
# Takes binary matrix and obtains Germline Mutation and LOH freqs. (in %)
# to use include in the final 2-way GLM's outputs
freq_vectors_list <- obtain_event_freq_vectors(germsomloh_matrix_patho_ct)

germline_drivers_to_test <- freq_vectors_list$germline_freqs[freq_vectors_list$germline_freqs > 0]
somatic_drivers_to_test <- freq_vectors_list$somatic_freqs[freq_vectors_list$somatic_freqs >= filt_B_freq]

freqs_columns_complete <- data.frame(expand.grid(Gene_A = names(germline_drivers_to_test), Gene_B = names(somatic_drivers_to_test)),
                                              expand.grid(Germline_Freq = germline_drivers_to_test, Somatic_Freq = somatic_drivers_to_test))
germ_som_freqs_columns <- freqs_columns_complete %>% unite(gene_pair, c(Gene_A, Gene_B))


# Obtain Model Inputs for .RDS
geneA_geneB_2way_glm_input <- setNames(lapply( names(germline_drivers_to_test), function(germ_gene) {
    cont_table <- setNames(lapply( names(somatic_drivers_to_test) , function(somatic_gene) {
        if (germ_gene == somatic_gene | ( germ_gene == 'TP53BP1' & somatic_gene == 'TP53' )) {
            return() # Not interested in cis- interactions
        } else {
            gene_germ = factor(germsomloh_matrix_patho_ct[[paste0(germ_gene, ".germline_mut")]], levels = c(0,1))
            gene_som = factor(germsomloh_matrix_patho_ct[[paste0(somatic_gene, ".somatic_mut")]], levels = c(0,1))
            input = as.data.frame(table(gene_germ, gene_som))
            colnames(input) = c("Germline_A", "Somatic_B", "N")
            if (add_pseudocounts) {input$N = input$N + 1}
            return(input[order(input$Germline_A, input$Somatic_B),])            
        }
    }), paste0(germ_gene, '_', names(somatic_drivers_to_test) ) )
    return(Filter(Negate(is_empty), cont_table))    
}), names(germline_drivers_to_test) )

geneA_geneB_2way_glm_input_flt <- unlist(geneA_geneB_2way_glm_input, recursive = FALSE) # Combines all gene-pairs contingency tables into same level (flattens)
names(geneA_geneB_2way_glm_input_flt) <- gsub(x = names(geneA_geneB_2way_glm_input_flt), pattern = '([A-Z0-9]+)[.]', replacement = '' )
saveRDS(geneA_geneB_2way_glm_input_flt, file = out_geneA_geneB_glm_inputs_rds)


# |- Table of 2-way GLM Inputs -| #
#
# Goes through each gene of the target list and returns contingency table
# results as a one-dimensional vector that will be used as an easy way
# to check model results
geneA_geneB_2way_glm_table <- setNames(lapply( names(germline_drivers_to_test), function(germ_gene) {
    cont_table <- setNames(lapply( names(somatic_drivers_to_test) , function(somatic_gene) {
        if (germ_gene == somatic_gene | ( germ_gene == 'TP53BP1' & somatic_gene == 'TP53' )) {
            return( c("Size" = 0, "NoGermA_NoSomB" = 0, "GermA_NoSomB" = 0, "NoGermA_SomB" = 0, "GermA_SomB" = 0)) 
        } else {
            gene_germ = factor(germsomloh_matrix_patho_ct[[paste0(germ_gene, ".germline_mut")]], levels = c(0,1))
            gene_som = factor(germsomloh_matrix_patho_ct[[paste0(somatic_gene, ".somatic_mut")]], levels = c(0,1))
            input = as.data.frame(table(gene_germ, gene_som))
            colnames(input) = c("Germline_A", "Somatic_B", "N")
            if (add_pseudocounts) {input$N = input$N + 1}
            contingency_row = as.vector(input$N)
            names(contingency_row) = c("NoGermA_NoSomB", "GermA_NoSomB", "NoGermA_SomB", "GermA_SomB")
            row <- c('Size' = sum(contingency_row), contingency_row)
            return(row)
        }
    }), paste0(germ_gene, '_', names(somatic_drivers_to_test) ) )
    missing_genesymbols_cols <- do.call(rbind, cont_table)
    return( cbind( 'gene_pair' = names(cont_table) , missing_genesymbols_cols)  )
}), names(germline_drivers_to_test) )

all_pairs_table <- as.data.frame(do.call(rbind, geneA_geneB_2way_glm_table))
all_pairs_table <- all_pairs_table %>% filter(Size != 0)
glm_table_final <- merge(germ_som_freqs_columns, all_pairs_table, by = 'gene_pair')
          
write_tsv(glm_table_final, file = out_geneA_geneB_glm_inputs_table)