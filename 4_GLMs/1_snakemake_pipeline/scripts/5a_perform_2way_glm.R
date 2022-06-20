##########################################
#--{ Perform 2-way Interaction Models }--#
##########################################

# Reads the list of inputs (contingency tables per gene) to apply
# 2-way GLMs accounting for the interaction between Germline Mutation:LOH
# Returns list of outputs per gene and a single table with some 
# of the relevant coefficients and statistics from each model

library(stats)
library(tidyverse)

input_2way_model_input_rds <- snakemake@input[['twoway_model_input_rds']]
input_2way_model_table <- snakemake@input[['twoway_model_table']]
  
out_2way_glm_rds <- snakemake@output[['twoway_model_output_rds']]
out_2way_glm_tsv <- snakemake@output[['twoway_model_output_table']]

# |- Data Imports -| #
input_2way_glm <- readRDS(input_2way_model_input_rds)
model_table <- read_tsv(input_2way_model_table)

# Formula for the models
model_cis_2_way_nopop <- "N ~ germline_mut + loh + germline_mut:loh"

# Loops for each gene and stores GLM output object into a list
glm_outputs <- setNames(lapply(names(input_2way_glm), function(gene){
  glm_germut_loh <- glm(formula = model_cis_2_way_nopop, data = input_2way_glm[[gene]],
                        family = poisson(link = "log"), 
                        control = glm.control(epsilon = 1e-6, maxit=100))
  return(glm_germut_loh)
}), names(input_2way_glm))
saveRDS(glm_outputs, file = out_2way_glm_rds)

# Extract some of the main coefficients/output of the GLM to ease its interpretation
glm_results_rows <- lapply(names(input_2way_glm), function(gene){
  glm_germut_loh <- glm(formula = model_cis_2_way_nopop, data = input_2way_glm[[gene]],
                        family = poisson(link = "log"), 
                        control = glm.control(epsilon = 1e-6, maxit=100))
  vec <- c(gene, glm_germut_loh$aic, summary(glm_germut_loh)$coefficients['germline_mut1:loh1', ])
  return(vec)
})

glm_interaction_results <- as.data.frame(do.call(rbind, glm_results_rows)) # Combines all row into df
colnames(glm_interaction_results) <- c('Gene', 'AIC', 'interaction_estimate', 'std_error', 'z_value', 'p_value')
twoway_glm_final_table <- merge(model_table, glm_interaction_results, by = 'Gene') # Merges with Germline and LOH Freqs. as well as contingency
write_tsv(twoway_glm_final_table, file = out_2way_glm_tsv)