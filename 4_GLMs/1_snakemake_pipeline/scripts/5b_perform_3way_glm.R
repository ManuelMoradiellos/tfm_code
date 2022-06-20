##########################################
#--{ Perform 3-way Interaction Models }--#
##########################################

# Reads the list of inputs (contingency tables per gene) to apply
# 2-way GLMs accounting for the interaction between Germline Mutation:LOH
# Returns list of outputs per gene and a single table with some 
# of the relevant coefficients and statistics from each model

library(tidyverse)
library(stats)

input_3way_model_input_rds <- snakemake@input[['threeway_model_input_rds']]
input_3way_model_table <- snakemake@input[['threeway_model_table']]
  
out_3way_glm_rds <- snakemake@output[['threeway_model_output_rds']]
out_3way_glm_tsv <- snakemake@output[['threeway_model_output_table']]


input_3way_glm <- readRDS(input_3way_model_input_rds)
model_table <- read_tsv(input_3way_model_table)

# Formula for the models
## Equivalent to "N ~ germline_mut + loh + pop + germline_mut:loh + germline_mut:pop + loh:pop + germline_mut:loh:pop"
model_cis_3_way_pop <- "N ~ germline_mut * loh * pop" 

# Loops for each gene and stores GLM output into a list of results per gene
glm_outputs_3way <- setNames(lapply(names(input_3way_glm), function(gene){
  glm_germut_loh_pop <- glm(formula = model_cis_3_way_pop, data = input_3way_glm[[gene]],
                            family = poisson(link = "log"), 
                            control = glm.control(epsilon = 1e-6, maxit=100))
  return(glm_germut_loh_pop)
}), names(input_3way_glm))
saveRDS(glm_outputs_3way, file = out_3way_glm_rds)

# Extract some of the main coefficients/output of the GLM to ease its interpretation
# These are the ones mainly used by Solip
glm_3way_rows <- setNames(lapply(names(input_3way_glm), function(gene){
  glm_germut_loh_pop <- glm(formula = model_cis_3_way_pop, data = input_3way_glm[[gene]],
                            family = poisson(link = "log"), 
                            control = glm.control(epsilon = 1e-6, maxit=100))
  # Search the columns we are interested in
  inter_names <- rownames(summary(glm_germut_loh_pop)$coefficients)[grepl('^germline_mut1:*', rownames(summary(glm_germut_loh_pop)$coefficients))][-1]
  coeffs <- t(summary(glm_germut_loh_pop)$coefficients[inter_names, c(1, 4)]) # Only focus on a few columns
  vec <- as.vector(coeffs)
  names(vec) <- paste0(sort(rep(colnames(coeffs), 2)), c('_estimate', '_p_value')) # Adequate column names following the correct order
  return(vec)
}), names(input_3way_glm))

out_3way_table <- do.call(rbind, glm_3way_rows)
out_3way_table <- data.frame(cbind('Gene' = rownames(out_3way_table), out_3way_table)) # Adds genes for future merges
threeway_glm_final_table <- merge(out_3way_table, model_table, by = 'Gene') # Merges with table of proportions/contigency counts
colnames(threeway_glm_final_table) <- gsub(x = colnames(threeway_glm_final_table), pattern = '\\.', replacement = ':')
write_tsv(threeway_glm_final_table, file = out_3way_glm_tsv)