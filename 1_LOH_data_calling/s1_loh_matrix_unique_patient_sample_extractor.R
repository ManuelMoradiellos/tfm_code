## As some patients are represented multiple times due to having samples from
# different origin (TCGA Codes from 01-06, without 04) these need to be
# 'harmonized' to the same type present in the germline variant calling step.

library('tidyverse')

loh_complete_matrix <- read_tsv('../input_files/gam_loh_complete_matrices/loh/loh_matrix_complete_corrected.tsv')

# The codes present refer to following levels: Primary Solid Tumor, Recurrent Solid Tumor,
# Primary Blood Derived Cancer - Peripheral Blood, Additional - New Primary, Metastatic

# This code extracts all samples and prioritizes those referent to Primary
# Solid Tumor or uses the next one that's available

test <- data.frame('samples' = loh_complete_matrix$patients,
                   'barcode' = substr(x = loh_complete_matrix$patients, start = 1, stop = 12),
                   'sub' = as.numeric(substr(x = loh_complete_matrix$patients, start = 14, stop = 15)))
to_keep <- test %>% group_by(barcode) %>% filter(sub == min(sub)) %>% pull(samples) # Obtain samples in 'ranked' order from 01 to 06
loh_nondup_matrix <- loh_complete_matrix %>% filter(patients %in% to_keep)
colnames(loh_nondup_matrix)[1] <- 'samples'
write_tsv(loh_nondup_matrix, file = '../input_files/gam_loh_complete_matrices/loh/loh_matrix_nodup_corrected.tsv')
