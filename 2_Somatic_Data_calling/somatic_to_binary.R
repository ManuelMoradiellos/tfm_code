library(tidyverse)

ref_data_cell <- read_tsv('/local/mmoradiellos/work_stuff/6_3way_two_genes_models/CELL2018_CancerFitness_All_All.txt') # Need to add reference for this somatic mutation data
germloh_merged_patient_barcodes <- scan('/local/mmoradiellos/work_stuff/6_3way_two_genes_models/gam_loh_merged_patient_barcodes.txt', what = 'character') # All merged GermlineMut + LOH binary matrices have the same patients,

# Keep common patients with the data used for the 2-way models and then all columns that are related to somatic mutation 
somatic_mut_ref <- ref_data_cell %>% select( c('Event', contains(':mutation')) ) %>% filter(Event %in% germloh_merged_patient_barcodes) 
colnames(somatic_mut_ref) <- gsub(pattern = ':mutation', replacement = '.somatic_mut', x = colnames(somatic_mut_ref) )
colnames(somatic_mut_ref)[1] <- 'Patient.Barcode' 
