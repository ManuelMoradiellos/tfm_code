library(tidyverse)
library(stats)

in_geneA_geneB_glm_inputs_rds <- snakemake@input[['geneA_geneB_glm_inputs_rds']]
in_geneA_geneB_glm_inputs_table <- snakemake@input[['geneA_geneB_glm_inputs_table']]
  
out_geneA_geneB_glm_outputs_rds <- snakemake@output[['geneA_geneB_glm_outputs_rds']]
out_geneA_geneB_glm_outputs_table <- snakemake@output[['geneA_geneB_glm_outputs_table']]

# |- Data Imports -| #
input_2way_glm <- readRDS(in_geneA_geneB_glm_inputs_rds)
model_table <- read_tsv(in_geneA_geneB_glm_inputs_table)

# Formula for the models
formula_trans_2_way <- "N ~ Germline_A * Somatic_B"

# Loops for each gene pair and stores GLM output object into a list
glm_outputs <- setNames(lapply(names(input_2way_glm), function(gene_pair){
  glm_germsoh <- glm(formula = formula_trans_2_way, data = input_2way_glm[[gene_pair]],
                        family = poisson(link = "log"), 
                        control = glm.control(epsilon = 1e-6, maxit=100))
  return(glm_germsoh)
}), names(input_2way_glm))
saveRDS(glm_outputs, file = out_geneA_geneB_glm_outputs_rds)

# Extract some of the main coefficients/output of the GLM to ease its interpretation
glm_results_rows <- lapply(names(input_2way_glm), function(gene_pair){
  glm_germsoh <- glm(formula = formula_trans_2_way, data = input_2way_glm[[gene_pair]],
                        family = poisson(link = "log"), 
                        control = glm.control(epsilon = 1e-6, maxit=100))
  vec <- c(gene_pair, glm_germsoh$aic, summary(glm_germsoh)$coefficients['Germline_A1:Somatic_B1', ])
  return(vec)
})
glm_interaction_results <- as.data.frame(do.call(rbind, glm_results_rows)) # Combines all row into df
colnames(glm_interaction_results) <- c('gene_pair', 'GeneAGeneB_AIC', 'GeneAGeneB_interaction_estimate', 'GeneAGeneB_std_error', 'GeneAGeneB_z_value', 'GeneAGeneB_p_value')
twoway_glm_final_table <- merge(model_table, glm_interaction_results, by = 'gene_pair') # Merges with Germline and Somatic Freqs. as well as contingency
write_tsv(twoway_glm_final_table, file = out_geneA_geneB_glm_outputs_table)