#############################################
#--{ Filter Merged Matrix by Cancer Type }--#
#############################################

# Takes a pathogenic class merged MUT-LOH event merged binary
# matrix and filters it according to the cancer type of its patients

library('tidyverse')

input_merged_matrix <- snakemake@input[['merged_matrix_by_patho']]
  
out_merged_matrix_by_patho_ct <- snakemake@output[['merged_matrix_by_patho_ct']]
  
metadata_table <- snakemake@params[['metadata_table']]
current_cancer <- snakemake@params[['cancer_type']]

# |- Data Import -| #
gamloh_merged <- read_tsv(input_merged_matrix)
patient_metadata <- read_tsv(metadata_table)

#  ¯¯¯¯¯¯¯¯¯¯¯¯¯¯  #
# |- Main Tasks -| #
#  --------------  #

if (current_cancer == 'Pancancer') {
  write_tsv(gamloh_merged, file = out_merged_matrix_by_patho_ct) # Keep all samples for Pancancer level
} else {
  cancer_patients <- patient_metadata %>% filter(cancer_type == current_cancer) %>% pull(bcr_patient_barcode) # Filter by Population Group
  write_tsv(gamloh_merged %>% filter(Patient.Barcode %in% cancer_patients), file = out_merged_matrix_by_patho_ct)
} 