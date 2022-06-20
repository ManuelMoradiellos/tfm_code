#########################
# Perpetrated by Manuel # 
#########################

library('tidyverse')

# This script takes the TCGA barcodes from the chr1 'rawest' file from
# Luis' filtering steps and then creates a common Population Column in PAM 
# naming checking both PAM cluster and PanCancer Supplementary Info to
# fill in the missing gaps.

# We use chr1 as it contains all samples
# comes from 'cut -f139 all.chr.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.SNV.TCGA_AF.header.chr1.tsv' 
# located in cluster's lgjimeno/tcga_pipeline_divide_by_chr_clean_restart_B_complete/out_pcac
chr1_barcodes <- read_tsv('pcac_chr1_samplesIDs.13-01-22.uniq.tsv') # shortened the name
chr1_barcodes$SAMPLE <- gsub(x = chr1_barcodes$SAMPLE, pattern = '-...-...-....-..', replacement = '')
colnames(chr1_barcodes) <- 'bcr_patient_barcode'

# Now reads both annotation files and creates a merged version with the previous one
pam_pop_info <- read_tsv('../old_filtered_tables/population_cluster_DriveFile/TCGApopinfo.cluters.and.subcluster.tsv')
colnames(pam_pop_info)[1] <- 'bcr_patient_barcode'
pam_merge <- merge(chr1_barcodes, pam_pop_info, by = 'bcr_patient_barcode')

supp_pop_info <- read_tsv("../old_filtered_tables//population_pancancer_sup1/TCGA.PanCan_CDR_Sup1.original_clinical_data.tsv",
                          col_select = c(1, 5))
colnames(supp_pop_info)[2] <- 'PanCancer_race'
supp_merge <- merge(chr1_barcodes, supp_pop_info, by = 'bcr_patient_barcode')

# Creates a combined matrix of the two annotations, keeping rows with missing
# data from any column (mainly PCA.PAM) to fill in coming steps
new_annotation <- merge(supp_merge, pam_merge, all.x = TRUE)

similar_clusters <- c('ASIAN', 'BLACK OR AFRICAN AMERICAN', 'WHITE')
 
new_annotation$mixed.PAM.Cluster <- sapply(1:nrow(new_annotation), function(x){
  if (is.na(new_annotation[x, 'PCA.PAM.Ancestry.Cluster'])) {
    
    if (new_annotation[x, 'PanCancer_race'] == 'ASIAN') {
      return('ASIAN')
    }
    if (new_annotation[x, 'PanCancer_race'] == 'BLACK OR AFRICAN AMERICAN') {
      return('AFR')
    }
    if (new_annotation[x, 'PanCancer_race'] == 'WHITE') {
      return('EUR')
    }
    if (!(new_annotation[x, 'PanCancer_race'] %in% similar_clusters)) {
      return(NA)
    }
  } else {return(new_annotation[x, 'PCA.PAM.Ancestry.Cluster'])}
})

# Now I do some checks recommended by Luis

## First, to see the correspondance of the new column created
check_1 <- new_annotation[which(is.na(new_annotation$PCA.PAM.Ancestry.Cluster)), ]
table(check_1$PanCancer_race, check_1$mixed.PAM.Cluster, useNA = 'always')

write_tsv(new_annotation, file = 'TCGA_barcode_population_PanCancerSupp_PopInfoPAM.13-01-22.tsv')

## This is needed for the snakemake pipeline
pop_counts <- table(new_annotation$mixed.PAM.Cluster, useNA = 'always')
saveRDS(pop_counts, file = 'PAM_pop_counts.RDS')


###############
# Seq. Center #
###############
chr1_barcodes <- read_tsv('pcac_chr1_samplesIDs.13-01-22.uniq.tsv') # shortened the name
chr1_barcodes$BARCODE <- gsub(x = chr1_barcodes$SAMPLE, pattern = '-...-...-....-..', replacement = '')
chr1_barcodes$CENTER <- gsub(x = chr1_barcodes$SAMPLE, pattern = '....-..-....-...-...-....-', replacement = '')
colnames(chr1_barcodes) <- c('SAMPLE', 'TCGA.PatientID', 'Center.Code')

center_info <- read_tsv('../../TCGA.Center_Codes.tsv', col_select = c(1,5))
colnames(center_info) <- c('Center.Code', 'Center.ShortName')
patientseqcenter <- merge(chr1_barcodes, center_info, by = 'Center.Code')

patientseqcenter <- patientseqcenter[, 3:4]
write_tsv(patientseqcenter, 'TCGA.PatientBarcode.SeqCenter.13-01-22.tsv')