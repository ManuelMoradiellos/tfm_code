library(tidyverse)

## One of the encountered issues on creating the matrix was that there were a lot of
## NAs for some patient-gene entries as the regions they are located within did not
## appeared to be studied in some patients.  In the original study, they used the
## ABSOLUTE algorithm/method for the LOH calls, it seems that they only check some regions
## that may have differing modal absolute copy number of each sample and then double-check
## to see if they're cancerous cell and really have that LOH going on (with various indicators).
## So our criteria is that all genes that appear as NA for a sample weren't suspect of LOH
## to begin with and, thus, we encode them as LOH = 0 (comparing the graphs and freq. tables
## it seems to check out with the references, so it was a good approach )  

## |-{ R CODE }-| ##
mapping <- read.table("2_out_genome_mapping/mapped_samples.tsv", header = FALSE, sep ="\t")
mapping$V1 <- factor(mapping$V1, levels = unique(mapping$V1)) # Turn both into factors for table()
mapping$V3 <- factor(mapping$V3, levels = unique(mapping$V3)) 
mapping_1 <- filter(mapping, V2 == 1) ## Only filter LOH == 1 to lower the size of the table

# Checks whether the Sample-Gene pair appears and writes them as either
# 1 or 0 depending on whether that pair did not appear in the matrix
table_mapping <- as.data.frame.matrix(table(mapping_1$V1, mapping_1$V3))

# Sometimes the algorithm counts duplicates as there can be multiple
# LOH events, so those are shifted down to 1 as I'm interested into binary
table_mapping[table_mapping > 1] <- 1

# Adding a column with the samples barcodes
table_mapping$samples <- rownames(table_mapping)
table_mapping <- table_mapping %>% relocate(samples)

# Some patients are repeated as they have TCGA samples from various origins (From 01 - 06, w.o. 04)
# Primary Solid Tumor, Recurrent Solid Tumor, Primary Blood Derived Cancer - Peripheral Blood, Additional - New Primary, Metastatic
samples_type <- data.frame('samples' = table_mapping$samples,
                           'barcode' = substr(x = table_mapping$samples, start = 1, stop = 12),
                           'sub' = as.numeric(substr(x = table_mapping$samples, start = 14, stop = 15)))

to_keep <- samples_type %>% group_by(barcode) %>% filter(sub == min(sub)) %>% pull(samples)
table_mapping_nodup <- table_mapping %>% filter(samples %in% to_keep)
dim(table_mapping_nodup)  ## [1] 10625 19657  
write.table(table_mapping_nodup, file = "3_out_binary_LOH_matrix_R/loh_matrix_complete.tsv", sep = "\t")