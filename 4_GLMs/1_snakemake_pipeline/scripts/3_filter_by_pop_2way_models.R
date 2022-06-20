##################################################
#--{ Filter Merged Matrix by Population Group }--#
##################################################

# Takes a GermlineMut-LOH Merged Matrix filtered
# by cancer type and further filters by population
# group.
# This matrices are mostly relevant for 2-way 
# GermlineMuy:LOH interaction models as they are
# perfomed at each population level separatedly 
# (AFR, AMR, ASIAN, EUR and ALLpop -previous 
# combined as one single population-)

library(tidyverse)

merged_matrix_by_patho_ct <- snakemake@input[['merged_matrix_by_patho_ct']]
  
out_merged_matrix_by_patho_ct_pop <- snakemake@output[['merged_matrix_by_patho_ct_pop']]
  
metadata_table <- snakemake@params[['metadata_table']]
population <- snakemake@params[['population']]

# |- Data Import -| #
merged_matrix <- read_tsv(merged_matrix_by_patho_ct)

#  ¯¯¯¯¯¯¯¯¯¯¯¯¯¯  #
# |- Main Tasks -| #
#  --------------  #

if (population != 'ALLpop'){ 
  # Reads metadata and extracts patients of that population
  patient_metadata <- read_tsv(metadata_table) 
  patients_vec <- patient_metadata %>%
    filter(mixed.PAM.Cluster == population) %>% 
    pull(bcr_patient_barcode)
  
  # Extract rows for those patients when focusing on one population group
  merged_matrix_pop_filtered <- merged_matrix %>%  
    filter(Patient.Barcode %in% patients_vec)
} else { # ALLpop is a group that combines all populations as one
  merged_matrix_pop_filtered <- merged_matrix
}
write_tsv(merged_matrix_pop_filtered, file = out_merged_matrix_by_patho_ct_pop)