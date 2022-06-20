####################################################
#--{ Obtain Inputs for 3-way Interaction Models }--#
####################################################

# Using the GermlineMut-LOH table merged by pathogenic class,
# cancer type and all population group for each cancer type
# obtains the 3-way GLMs input by keeping those genes in which
# there is, at least, one account of GermlineMut or LOH event.
#
# Returns .RDS with a list of GLM inputs per gene and a 
# table with the same information in a wide format that also
# includes various columns with the 'diff'/'effect size' per
# available population group and a few extra comparison with
# the baseline group (as requested by Solip)

library(tidyverse)

in_merged_matrix_by_patho_ct <- snakemake@input[['merged_matrix_by_patho_ct_pop']]
print(in_merged_matrix_by_patho_ct)
merged_matrix_patho_ct <- read.table(in_merged_matrix_by_patho_ct, sep = "\t", header = TRUE)
somatic_mutation_matrix <- read.table(snakemake@input[['somatic_mutation_matrix']],sep = "\t", header = TRUE)

print(colnames(merged_matrix_patho_ct)[1:10])
print(colnames(somatic_mutation_matrix)[1:10])

germsomloh_matrix_patho_ct = merge(merged_matrix_patho_ct, somatic_mutation_matrix, by = "Patient.Barcode")

out_3way_AB_model_input_rds <- snakemake@output[['threeway_AB_model_input_rds']]
  
metadata_table <- snakemake@params[['metadata_table']]
patient_metadata <- read_tsv(metadata_table)
tested_genes_list = readRDS(snakemake@input[['list_of_genes_to_test']])

#EUR_patients = filter(patient_metadata, mixed.PAM.Cluster == "EUR")$brc_patient_barcode
#somatic_mutation_matrix <- read_tsv(snakemake@input[['somatic_mutation_matrix']])
#somatic_mutation_matrix_EUR = filter(somatic_mutation_matrix, Patient.Barcode %in% EUR_patients)

current_cancer <- snakemake@params[['cancer_type']]
patho <- snakemake@params[['patho']]
geneset <- snakemake@params[['geneset']]


add_pseudocounts <- snakemake@params[['add_pseudocounts']]
filt_B_freq <- as.numeric(snakemake@params[['filt_B_freq']])


## get model
germloh_genes = tested_genes_list[[patho]][[geneset]][[current_cancer]][[1]]
som_freqs = colMeans(germsomloh_matrix_patho_ct[, grepl(pattern = "somatic_mut$", colnames(germsomloh_matrix_patho_ct))])
somatic_genes = sapply(names(som_freqs)[som_freqs >= filt_B_freq], function(x) strsplit(x, split = "\\.")[[1]][1])
print(germloh_genes)
print(somatic_genes)

gml_input = lapply(germloh_genes, function(germloh_gene) {
    print(germloh_gene)
    temp1 = lapply(somatic_genes, function(somatic_gene) {
        print(somatic_gene)
        gene_germ = factor(germsomloh_matrix_patho_ct[[paste0(germloh_gene, ".germline_mut")]], levels = c(0,1))
        gene_loh = factor(germsomloh_matrix_patho_ct[[paste0(germloh_gene, ".loh_event")]], levels = c(0,1))
        gene_som = factor(germsomloh_matrix_patho_ct[[paste0(somatic_gene, ".somatic_mut")]], levels = c(0,1))
        input = as.data.frame(table(gene_germ, gene_loh, gene_som))
        print(dim(input))
        print(is(input))
        colnames(input) = c("Germline_A", "LOH_A",  "Somatic_B", "N")
        if (add_pseudocounts) {input$N = input$N + 1}
        return(input)
    })
    names(temp1) = somatic_genes
    return(temp1)    
})

names(gml_input) = germloh_genes

saveRDS(gml_input, out_3way_AB_model_input_rds)


