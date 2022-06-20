######################################################
#--{ GAM & LOH Matrices Merge by Pathogenic Class }--#
######################################################

# Takes a GAM (Binary Mutational Event matrix) and an LOH matrix
# (Binary Loss-of-Heterozygosity) and merges them by patient. 
# Each of these GAMs comes from a different set of pathogenic
# variants defined by our pathogenic classes (pLoF, Delmis, 
# Clinvar and Allclass -the previous combined-).

library(tidyverse)

input_gam <- snakemake@input[['gam_matrix']]
input_loh <- snakemake@input[['loh_matrix']]
  
output_merged_matrix_patho <- snakemake@output[['merged_matrix_by_patho']]

# |- Basic Imports -| #
gam_complete_matrix <- read_tsv(input_gam)
loh_complete_matrix <- read_tsv(input_loh)

#  ¯¯¯¯¯¯¯¯¯¯¯¯¯¯  #
# |- Main Tasks -| #
#  --------------  #
loh_complete_matrix$samples <- substr(x = loh_complete_matrix$samples,
                                       start = 1, stop = 12) # From Sample to Patient barcode (TCGA)

merged_matrix <- merge(gam_complete_matrix, loh_complete_matrix,
                       by.x = 'Patient.Barcode', by.y = 'samples') # Merges by rowname, keeping all column names with '.x' and '.y' for duplicates

merged_matrix <- merged_matrix[c(1, order(sub("\\.[xy]$", "", names(merged_matrix)[-1])) + 1)]  # Reorders data frame so columns for the same genes are next to each other

merged_matrix <- merged_matrix %>% select( matches( "\\.Barcode|\\.x|\\.y")  ) # Does not keep columns that were unique to one of the tables
# \\ may not be necessary
colnames(merged_matrix) <- str_replace_all(string = colnames(merged_matrix),
                                           c(".x" =  ".germline_mut",
                                             ".y" = ".loh_event")) # Rename columns according to matrix of origin

write_tsv(merged_matrix, file = output_merged_matrix_patho)