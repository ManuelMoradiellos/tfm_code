##########################################
#--{ Perform 3-way Interaction Models }--#
##########################################

# Reads the list of inputs (contingency tables per gene) to apply
# 2-way GLMs accounting for the interaction between Germline Mutation:LOH
# Returns list of outputs per gene and a single table with some 
# of the relevant coefficients and statistics from each model

library(tidyverse)
library(stats)

input_3way_model_input_rds <- snakemake@input[['threeway_AB_model_input_rds']]

out_3way_glm_rds <- snakemake@output[['threeway_AB_model_output_rds']]
print(snakemake@output[['threeway_AB_model_output_rds']])


input_3way_glm <- readRDS(input_3way_model_input_rds)

# Formula for the models
## Equivalent to "N ~ germline_mut + loh + pop + germline_mut:loh + germline_mut:pop + loh:pop + germline_mut:loh:pop"

message(colnames(input_3way_model_input_rds[[1]][[1]]))

model_cis_3_way_pop <- "N ~ Germline_A * LOH_A * Somatic_B" 

# Loops for each gene and stores GLM output into a list of results per gene
glm_outputs_3way <- setNames(lapply(names(input_3way_glm), function(gene_germloh){
    setNames(lapply(names(input_3way_glm[[gene_germloh]]), function(gene_som){
        glm_germut_loh_pop <- glm(formula = model_cis_3_way_pop, data = input_3way_glm[[gene_germloh]][[gene_som]],
                                    family = poisson(link = "log"), 
                                    control = glm.control(epsilon = 1e-6, maxit=100))
        return(glm_germut_loh_pop)
    }), names(names(input_3way_glm[[gene_germloh]])))
}), names(input_3way_glm))
saveRDS(glm_outputs_3way, file = out_3way_glm_rds)


