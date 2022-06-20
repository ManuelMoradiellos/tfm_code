#########################
# Perpetrated by Manuel #  
#########################

# Takes the results from the sample filtering and RGDV prioritization
# step and extracts the relevant columns that will be used for plink
# inside the PCA plotter script (performs PCAs and plots them)

library('tidyverse')

PCAplotter_input_from_luis_input <- function(input_file){
  
  filename <- strsplit(input_file, split = '.tsv')[[1]]
  
  output_dir <- paste0('plink_PCAplotter_inputs_', filename)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  start_tsv <- (read.table(input_file, header = TRUE, sep = "\t", row.names = 1))
  
  ######################
  # Obtain .tsv inputs #
  ######################
  
  PCAplotter_tsv <- start_tsv
  
  # Obtain the first nine columns of the .vcf from our variant IDs
  names(PCAplotter_tsv) <- gsub(x = names(PCAplotter_tsv), pattern = '\\....\\....\\.....\\...',
                           replacement = '')
  names(PCAplotter_tsv) <- gsub(x = names(PCAplotter_tsv), pattern = '\\.', replacement = '-')
  
  # ALL_POP does not need any filtering
  write_tsv(PCAplotter_tsv, file = paste0('./', output_dir, '/', filename, '.ALL_POP_plink_PCAplotter_input2.tsv'),quote = 'none')
  
  # We now read our population metadata (mixing two annotation files to fill in
  # the missing info)
  pop_info <- read_tsv(file = './TCGA_barcode_population_PanCancerSupp_PopInfoPAM.13-01-22.tsv',
                       col_select = c(1, 5))
  colnames(pop_info) <- c('TCGA.PatientID', 'PCA.PAM.Ancestry.Cluster') # For plink-PCA script
  
  for (pop in c('AFR', 'AMR', 'EAS', 'EUR')) {
    # Reads samples from annotation file that match each population cluster
    filtered_samples <- as.vector(pop_info[which(pop_info$PCA.PAM.Ancestry.Cluster == pop), ][[1]])
    filtered_df <- subset(PCAplotter_tsv, select = filtered_samples) 
    write_tsv(filtered_df, file = paste0('./', output_dir, '/', filename, '.', pop, '_plink_PCAplotter_input2.tsv'), quote = 'none')
  }
  
  ######################
  # Obtain .vcf inputs # 
  ######################
  
  PCAplotter_vcf <- start_tsv
  
  # Translate Luis' coding to standard genotype info
  PCAplotter_vcf[PCAplotter_vcf == 0] <- '0/0'
  PCAplotter_vcf[PCAplotter_vcf == 1] <- '0/1'
  PCAplotter_vcf[PCAplotter_vcf == 2] <- '1/1'
  
  # Obtain the first nine columns of the .vcf from our variant IDs
  names(PCAplotter_vcf) <- gsub(x = names(PCAplotter_vcf), pattern = '\\....\\....\\.....\\...',
                           replacement = '')
  names(PCAplotter_vcf) <- gsub(x = names(PCAplotter_vcf), pattern = '\\.', replacement = '-')
  
  split_list <- lapply(rownames(PCAplotter_vcf), function(x){
    sep <- strsplit(x, split = '_')[[1]]
    return(c(sep[1], sep[2], x,  sep[3], sep[4], sep[5], 'PASS', sep[5], 'GT'))
  } )
  # Binds all results of previous function to create a matrix
  split_df <- as.data.frame(do.call(what = rbind, args = split_list), row.names = row.names(PCAplotter_vcf))
  
  # Adds the names of those columns as expected
  vcf_main <- c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                'INFO', 'FORMAT')
  colnames(split_df) <- vcf_main
  final_vcf <- merge(split_df, PCAplotter_vcf, by = 0) # Merged by original rownames
  final_vcf <-  final_vcf[,-1] # Removes rownames 
  
  # ALL_POP does not need any filtering
  write_tsv(final_vcf, file = paste0('./', output_dir, '/', filename, '.ALL_POP.prevcf.tsv'),quote = 'none')
  
  system(command = sprintf('cat ../mock_vcfH_header.txt %s/%s.ALL_POP.prevcf.tsv > %s/%s.ALL_POP_plink_PCAplotter_input1.vcf',
                           output_dir, filename, output_dir, filename))
  
  # We now read our population metadata (mixing two annotation files to fill in
  # the missing info)
  pop_info <- read_tsv(file = './TCGA_barcode_population_PanCancerSupp_PopInfoPAM.13-01-22.tsv',
                       col_select = c(1, 5))
  colnames(pop_info) <- c('TCGA.PatientID', 'PCA.PAM.Ancestry.Cluster') # For plink-PCA script
  
  for (pop in c('AFR', 'AMR', 'EAS', 'EUR')) {
    # Reads samples from annotation file that match each population cluster
    filtered_samples <- as.vector(pop_info[which(pop_info$PCA.PAM.Ancestry.Cluster == pop), ][[1]])
    filtered_df <- subset(final_vcf, select = c(vcf_main, filtered_samples)) # Subsets to keep first nine cols
    # and the samples of each population
    write_tsv(filtered_df, file = paste0('./', output_dir, '/', filename, '.', pop, '.prevcf.tsv'), quote = 'none')
    system(command = sprintf('cat ../mock_vcfH_header.txt %s/%s.prevcf.tsv > %s/%s_plink_PCAplotter_input1.vcf',
                             output_dir, paste0(filename, '.', pop), output_dir,
                             paste0(filename, '.', pop)))
  }
}

PCAplotter_input_from_luis_input('PCA.input.3levels.common_0.05_exonic.tsv') # Done for all tested Common Variants MAFs