#################################################
#--{ Merge EUR 2-way and 3-way output tables }--#
#################################################

# This step just merges the relevant coefficients and contingency 
# table counts from the 2-way GLMs for EUR population (baseline)
# and the 3-way GLMs outputs to showcase a comparison of both
# models in a single table.

library(tidyverse)

out_2way_glm_tsv <- snakemake@input[['twoway_model_output_table']]
out_3way_glm_tsv <- snakemake@input[['threeway_model_output_table']]
  
out_merged_glm_outputs_tsv <- snakemake@output[['merged_glm_outputs_tsv']]

table_2way <- read_tsv(out_2way_glm_tsv)
table_3way <- read_tsv(out_3way_glm_tsv)
colnames(table_3way) <- gsub(x = colnames(table_3way), pattern = '\\.', replacement = ':') # When its coerced to a tibble it loses the formatting of '.' as ':' in colnames 
write_tsv(merge(table_2way, table_3way, by = 'Gene'), file = out_merged_glm_outputs_tsv)