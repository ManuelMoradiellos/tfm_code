######################################################
#--{ Count Significant Genetic Interactions 2-way }--# # Perpetrated by Manuel Moradiellos
###################################################### 

# Takes the results from the RDGV calling step and compares
# the enrichment results of the four pathogenic classes for
# three gene set using a Wilcoxon Test and boxplots.

library('tidyverse')
library('patchwork')
library('formattable')
library('htmltools')
library('webshot')
library('ggsignif')

#############
# Functions #
#############

gam_obtainer_patientID <- function(patho_table){
  # Returns table of alterations
  patho_table$consensus_gene_symbol <- factor(patho_table$consensus_gene_symbol, levels = sort(unique(all_classes_annotable$consensus_gene_symbol))) # To account for all possibilities
  patho_table$bcr_patient_barcode <- factor(patho_table$bcr_patient_barcode, levels = sort(unique(all_classes_annotable$bcr_patient_barcode)))
  gam <- table(patho_table$consensus_gene_symbol, patho_table$bcr_patient_barcode) ## Maybe including NAs is not right, adds one to each dimension
  gam[gam > 1] <- 1  # We are only interested in a binary table of alterations per gene
  return(gam)
}

gam_obtainer_sampleID <- function(patho_table){
  # Returns table of alterations
  patho_table$consensus_gene_symbol <- factor(patho_table$consensus_gene_symbol, levels = sort(unique(all_classes_annotable$consensus_gene_symbol))) # To account for all possibilities
  patho_table$SAMPLE <- factor(patho_table$SAMPLE, levels = sort(unique(all_classes_annotable$SAMPLE)))
  gam <- table(patho_table$consensus_gene_symbol, patho_table$SAMPLE) ## Maybe including NAs is not right, adds one to each dimension
  gam[gam > 1] <- 1  # We are only interested in a binary table of alterations per gene
  return(gam)
}

mut_freq_across_cancer <- function(gam, patient_cancer_table = patient_cancerabv){
  ## Extract samples and number of alterations of each gene
  blank_matrix <- matrix(NA, nrow = nrow(gam), ncol = length(unique(patient_cancer_table$cancer_type_abbreviation))) # Creates blank matrix of size equal to input's
  colnames(blank_matrix) <- unique(patient_cancer_table$cancer_type_abbreviation)                                    # and with the same rownames (gene symbols)
  rownames(blank_matrix) <- rownames(gam)
  
  # Goes trough each of the cancer types and extracts the patient samples ID
  # that had each of those, to then obtain the necessary subset of genes 
  # alterated in each cancer type to compute their alteration frequency
  # (num. of single alterations for each cancer type divided by size of genes universe)
  for (cancer_type in unique(patient_cancer_table$cancer_type_abbreviation)) {
    patients <- patient_cancer_table %>% 
      filter(cancer_type_abbreviation == cancer_type) %>% 
      pull(patientID)
    check <- colnames(gam) %in% patients # Avoids 'subscript out of bounds' error
    #print(paste0(cancer_type, ' : ', sum(check)))
    blank_matrix[ , cancer_type] <- (rowMeans(gam[, check])) # Obtains the real alteration freq. for all possible 
  }                                                                           # gene and then adds it to the blank matrix under  
  return(blank_matrix)                                                        # each corresponding cancer type 
}

mut_freq_pancancer <- function(mut_freq_table){
  # Converts into final dataframe and computes mutation frequency in Pancancer level
  mut_freq_table <- data.frame(Gene = row.names(mut_freq_table), mut_freq_table) #Converts to data frame
  mut_freq_table$Pancancer <- rowMeans(mut_freq_table[, -1]) # Adds Pancancer freq.
  mut_freq_table <- mut_freq_table %>%
    mutate(max_freq = do.call(pmax, c(select(., c(colnames(mut_freq_table[, -c(1, 34)])))))) # Add max. freq. for each gene across cancer types
  return(mut_freq_table)
}

mut_freq_subset_and_plot <- function(matrix, patho_class = 'pLoF', ycoord_cut = 0.005,
                                 cpgs = CPGs_list_available, ccg = ccg_list_available_not_CPG, rando = random_genes_1K){
  # Subsets each mutation frequency table and creates violin plots
  # for each one of those and then combines them in the same plot
  
  if (!dir.exists( paste0(patho_class, '_outdir') )){
    dir.create( paste0(patho_class, '_outdir') )
  }
  
  #-----------------------#
  #- Matrices Subsetting -#
  #-----------------------#
  mutfreq_matrix_cpgs <- data.frame(matrix) %>% 
    filter(rownames(matrix) %in% CPGs_list_available)
  
  mutfreq_matrix_ccg <- data.frame(matrix) %>% 
    filter(rownames(matrix) %in% ccg_list_available_not_CPG)
  
  mutfreq_matrix_random <- data.frame(matrix) %>% 
    filter(rownames(matrix) %in% random_genes_1K)
  
  #----------------#
  #- Violin Plots -#
  #----------------#
  to_plot <- list('CPGs' = mutfreq_matrix_cpgs, 'Somatic Drivers' = mutfreq_matrix_ccg,
                  '1k Random Gene Set' = mutfreq_matrix_random)
  # Goes through each of the subsets to plot them
  violin_plots <- lapply(names(to_plot), function(x){
    alteration_freq_violin_plotter_patch(matrix = to_plot[[x]], condition = patho_class,
                                         gene_subset = x, ycoord_cut = ycoord_cut)
  })
  violin_plots <- lapply(violin_plots, function(x){ x$labels$title <- NULL ; return(x) }) # Maybe does nothing as there's 
                                                                                          # no title at each plot level
  leg <- cowplot::get_legend(violin_plots[[1]] + theme(legend.position = "top")) # Extracts legend to put in onto the top
  
  ggsave(plot = ((cowplot::plot_grid(leg) / (violin_plots[[1]] / violin_plots[[2]] / violin_plots[[3]])) + 
                   plot_layout(heights = unit(c(1, 17), "cm")) +
                   plot_annotation(title = paste0("Germline Variant Frequency across Cancer Types for ", patho_class, "class")) +
                   plot_layout(guides = "collect", ncol = 1) & theme(legend.position = "none", plot.title = element_text(size = 21))),
         path = paste0('./', patho_class, '_outdir'), filename = paste0('germline_var_freq_', patho_class, '_violin_plots.png'), width = 14.4, height = 9.09, units = "in", dpi = 300, device = png)
}


export_formattable <- function(f, file, width = "100%", height = NULL, 
                               background = "white", delay = 0.2){
  # https://stackoverflow.com/questions/38833219/command-for-exporting-saving-table-made-with-formattable-package-in-r
  # "Export" formattable results as .pngs #
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}


#-------#
#-PLOTS-#
#-------#

alteration_freq_violin_plotter <- function(matrix, condition = 'Clinvar+Delmis+pLoF', ycoord_cut = 0.005){
  df_to_plot <- data.frame(matrix[, order(colMeans(matrix), decreasing = TRUE)]) # Orders column by mean
  df_to_plot %>% gather(col, val) %>% 
    mutate(col = factor(col, unique(col))) %>% 
    ggplot(aes(col, val, colour = col)) + 
    geom_violin(aes(fill = col), lwd = 0.3, alpha = 0.8, show.legend = FALSE) + # Removes colour from legend
    ggtitle(paste0('Germline Variant Frequency across Cancer Types for ', condition , ' class')) +
    theme_update(plot.title = element_text(hjust = 0.5, vjust = 1, size = 18)) +
    #labs(subtitle = paste0('Ordered by mean alteration frequency for ', nrow(df_to_plot), ' genes')) + ## Not so useful anymore, working with genes as factors always returns total
    labs(subtitle = 'Ordered by mean alteration frequency') +
    theme_update(plot.subtitle = element_text(vjust = 1, hjust = 0.5, size = 11)) +
    xlab('TCGA Cancer Type') +
    ylab('Variant Frequency') +
    theme(axis.title.x = element_text(size = 14, vjust = -0.8, margin = margin(t = 0, r = 0, b = 10, l = 0)),
          axis.text.x = element_text(size = 9, angle = 45, vjust = 0.6),
          axis.text.y = element_text(size = 9),
          axis.title.y = element_text(size = 14, vjust = 0.5, margin = margin(t = 0, r = 15, b = 0, l = 2.5))) +
    coord_cartesian(ylim = c(0, ycoord_cut)) + # Change limit of y-axis to better see the detailing WITHOUT ALTERING ANY STAT SUMMARY COMPUTATION
    stat_summary(fun = "mean", geom = "crossbar", aes(color = "Mean"), na.rm = F, fatten = 1.5, width = 0.7, alpha = 0.6) +
    stat_summary(fun ="median", geom = "crossbar", aes(color = "Median"), na.rm = F, fatten = 1.5, width = 0.7, alpha = 0.6) +
    scale_color_manual("Statistics", values = c(Mean = "blue1", Median = "black")) 
}

alteration_freq_violin_plotter_patch <- function(matrix, condition = 'Clinvar+Delmis+pLoF', gene_subset = '', ycoord_cut = 0.005){
  df_to_plot <- data.frame(matrix) # Orders column by mean
  
  # Removes some axis titles depending on the plot to create a more cohesive design
  if (gene_subset == 'CPGs'){
    xaxis_title <- ''
    yaxis_title <- ''
  } else if (gene_subset == 'Cancer Consensus Gene') {
    xaxis_title <- ''
    yaxis_title <- 'Variant Frequency'
  } else if (gene_subset == 'Random Gene Set'){
    xaxis_title <- 'TCGA Cancer Type'
    yaxis_title <- ''
  }
  
  vplot <- df_to_plot %>% gather(col, val) %>% 
    mutate(col = factor(col, unique(col))) %>% 
    ggplot(aes(col, val, colour = col)) + 
    geom_violin(aes(fill = col), lwd = 0.3, alpha = 0.8, show.legend = FALSE) + # Removes colour from legend
    # ggtitle(paste0('Germline Variant Frequency across Cancer Types for ', condition , ' class')) +
    labs(subtitle = gene_subset) + ## Not so useful anymore, working with genes as factors always returns total
    theme_update(plot.subtitle = element_text(vjust = 1, hjust = 0.5, size = 13)) +
    xlab(xaxis_title) +
    ylab(yaxis_title) +
    theme(axis.title.x = element_text(size = 19, vjust = -0.8, margin = margin(t = 0, r = 0, b = 10, l = 0)),
          axis.text.x = element_text(size = 9, angle = 45, vjust = 0.6),
          axis.text.y = element_text(size = 9),
          axis.title.y = element_text(size = 19, vjust = 0.5, margin = margin(t = 0, r = 15, b = 0, l = 2.5)),
          legend.position = 'top') +
    coord_cartesian(ylim = c(0, ycoord_cut)) + # Change limit of y-axis to better see the detailing WITHOUT ALTERING ANY STAT SUMMARY COMPUTATION
    stat_summary(fun = "mean", geom = "crossbar", aes(color = "Mean"), na.rm = F, fatten = 1.5, width = 0.5, alpha = 0.4) +
    stat_summary(fun ="median", geom = "crossbar", aes(color = "Median"), na.rm = F, fatten = 1.5, width = 0.5, alpha = 0.4) +
    scale_color_manual("Statistics", values = c(Mean = "blue1", Median = "black"))
}

#################################################
# {|- Data Imports and basic table handling -|} #
#################################################

# Three main pathogenic groups
######################### NEWEST
clinvar_filtered_annotable <- read_tsv('latest_germline_pipeline_output/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.PEXT.clinvar.2nd.rare_0.001.header.tsv') # Non-overlapping with pLoF
plof_filtered_annotable <- read_tsv('latest_germline_pipeline_output/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.PEXT.plof.2nd.rare_0.001.header.tsv')
delmis_filtered_annotable <- read_tsv('latest_germline_pipeline_output/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.PEXT.delmis.2nd.rare_0.001.header.tsv')

# Combination of all models
exclu_models <- rbind(plof_filtered_annotable, delmis_filtered_annotable)
all_classes_annotable <- unique(rbind(exclu_models, clinvar_filtered_annotable)) # May have repeated variants between Clinvar and pLoF
write_tsv(all_classes_annotable, file = 'input_tables_march_25/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.PEXT.all_classes.2nd.rare_0.001.header.tsv')

## Read PatientBarcode (12 char.) - Cancer Type Metadata, then creates level for Colorrectal
# patient_metadata <- read_tsv('../../1_boolean_table_and_loh_events_count/input_files/Survival_SupplementalTable_S1_20171025_xena_sp_1') # We will use a consensus table instead (combines TCGA+GDC+Xena)
# patient_cancerabv <- patient_metadata[, c(2,3)]
#  patient_cancerabv <- patient_cancerabv %>% 
#     mutate(cancer_type_abbreviation = replace(cancer_type_abbreviation, (cancer_type_abbreviation == 'READ'|cancer_type_abbreviation == 'COAD'), 'COADREAD'))
patient_cancerabv <- read_tsv('clinical_annotation/consensus_cancer_type_table.tsv') # We will use a consensus table instead (combines TCGA+GDC+Xena)
colnames(patient_cancerabv) <-c('patientID', 'cancer_type_abbreviation')

#-|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-#            
#-Obtaining the subsets that will be used for the comparison-#
# All genes present in Pelayo's initial .vcf file
all_genes <- scan('genes_subsets_ALL_CCG_CPGs/genes_sorted_uniq.txt', what = "character")
all_genes_without_semmic <- unique(unlist(lapply(all_genes, function(x) strsplit(x, split = ";")[[1]]))) # Remove misleading entries

# CPGs to compare approaches, from our Gitlab; keep those present in our input file
CPGs_list <- scan(file = 'genes_subsets_ALL_CCG_CPGs/CPG.csv', what = 'character')
CPGs_list_available <- CPGs_list[CPGs_list %in% all_genes]

# Cancer Consensus Genes from COSMIC https://cancer.sanger.ac.uk/census (GRCh37-COSMIC v95)
ccg_list_tab <- read_tsv('genes_subsets_ALL_CCG_CPGs/Census_allMon_Feb_14_11_23_46_2022.tsv')
ccg_list <- unique(filter(ccg_list_tab, is.na(Germline), Somatic == "yes")$`Gene Symbol`)
ccg_list_available <- ccg_list[ccg_list %in% all_genes]
somatic_drivers <- ccg_list_available[!(ccg_list_available %in% CPGs_list_available)] # Truly somatic drivers

# Extracts cancer related genes and obtains a random subset of 1k genes which aren't among those 
all_cancer_related_genes <- unique(c(CPGs_list, ccg_list))
all_genes_without_semmic_and_cancer <- all_genes_without_semmic[!(all_genes_without_semmic %in% all_cancer_related_genes)]
set.seed(666)
random_genes_1K <- sample(all_genes_without_semmic_and_cancer, size = 1000, replace = FALSE)

#########################################
# {|- Main Work and Functions Usage -|} #
#########################################

# Using our Consensus Gene Symbol and TCGA's Patient Barcode we obtain the
# alteration matrix, only count if a gene is altered at least once per patient
gam_allmodels <- gam_obtainer_patientID(all_classes_annotable)
total_genes <- dim(gam_allmodels)[1] 

gam_plof <- gam_obtainer_patientID(plof_filtered_annotable)
gam_delmis <- gam_obtainer_patientID(delmis_filtered_annotable)
gam_clinvar <- gam_obtainer_patientID(clinvar_filtered_annotable)

# For each GAM we obtain a matrix of alteration 
# frequency for each gene at cancer level, 
mutfreq_all <- mut_freq_across_cancer(gam_allmodels)
mutfreq_plof <- mut_freq_across_cancer(gam_plof)
mutfreq_clinvar <- mut_freq_across_cancer(gam_clinvar)
mutfreq_delmis <- mut_freq_across_cancer(gam_delmis)

# Plotting
# alteration_freq_violin_plotter(mutfreq_all, ycoord_cut = 0.01)
# alteration_freq_violin_plotter(mutfreq_plof, ycoord_cut = 0.004, condition = 'pLoF')
# alteration_freq_violin_plotter(mutfreq_clinvar, ycoord_cut = 0.001, condition = 'Clinvar')
# alteration_freq_violin_plotter(mutfreq_delmis, ycoord_cut = 0.004, condition = 'Delmis')

# Converts to data frame and adds Pancancer frequency column
# to each class as well as maximum freq. found for gene
mutfreq_all_pan <- mut_freq_pancancer(mutfreq_all)
mutfreq_plof_pan <- mut_freq_pancancer(mutfreq_plof)
mutfreq_clinvar_pan <- mut_freq_pancancer(mutfreq_clinvar)
mutfreq_delmis_pan <- mut_freq_pancancer(mutfreq_delmis)

df_list <- setNames(list(mutfreq_clinvar_pan, mutfreq_plof_pan, mutfreq_delmis_pan, mutfreq_all_pan),
                    c('ClinVar', 'pLoF', 'DelMis', 'ALLclass'))

pancancer_boxplots <- lapply(names(df_list), function(patho_class){
  
  mutfreq_matrix_cpgs <- data.frame(df_list[[patho_class]]) %>% 
    filter(rownames(df_list[[patho_class]]) %in% CPGs_list_available)
  mutfreq_matrix_cpgs <- mutfreq_matrix_cpgs %>%
    mutate(max_freq = do.call(pmax, c(select(., c(colnames(mutfreq_matrix_cpgs[, -c(1, 34)]))))))
  mutfreq_matrix_cpgs$class <- rep('CPGs', dim(mutfreq_matrix_cpgs)[1])
  mutfreq_matrix_cpgs$class <- factor(mutfreq_matrix_cpgs$class, levels = c('CPGs', 'Somatic Drivers', '1K Random Genes'))

  mutfreq_matrix_ccg <- data.frame(df_list[[patho_class]]) %>% 
    filter(rownames(df_list[[patho_class]]) %in% ccg_list_available_not_CPG)
  mutfreq_matrix_ccg <- mutfreq_matrix_ccg %>%
    mutate(max_freq = do.call(pmax, c(select(., c(colnames(mutfreq_matrix_ccg[, -c(1, 34)]))))))
  mutfreq_matrix_ccg$class <- rep('Somatic Drivers', dim(mutfreq_matrix_ccg)[1])
  mutfreq_matrix_ccg$class <- factor(mutfreq_matrix_ccg$class, levels = c('CPGs', 'Somatic Drivers', '1K Random Genes'))

  mutfreq_matrix_random <- data.frame(df_list[[patho_class]]) %>% 
    filter(rownames(df_list[[patho_class]]) %in% random_genes_1K)
  mutfreq_matrix_random <- mutfreq_matrix_random %>%
    mutate(max_freq = do.call(pmax, c(select(., c(colnames(mutfreq_matrix_random[, -c(1, 34)]))))))
  mutfreq_matrix_random$class <- rep('1K Random Genes', dim(mutfreq_matrix_random)[1])
  mutfreq_matrix_random$class <- factor(mutfreq_matrix_random$class, levels = c('CPGs', 'Somatic Drivers', '1K Random Genes'))
  
  mutfreq_combined <- bind_rows(mutfreq_matrix_cpgs, mutfreq_matrix_ccg, mutfreq_matrix_random)
  write_tsv(x = mutfreq_combined, file = paste0('latest_germline_pipeline_output/germline_pancancer_freq_table_', patho_class, '_all_test_genesets.tsv'))
  
  if (patho_class == 'ALLclass'){
    xaxis_title <- ''
    yaxis_title <- ''
  } else if (patho_class == 'pLoF') {
    xaxis_title <- ''
    yaxis_title <- ''
  } else if (patho_class == 'ClinVar'){
    xaxis_title <- ''
    yaxis_title <- ''
  } else if (patho_class == 'DelMis'){
    xaxis_title <- ''
    yaxis_title <- ''
  }
  
  do_violins <- function(df, type) {
    plot <- ggplot(df, aes(x = class, y = Pancancer, fill = class)) + 
      geom_boxplot(notch = TRUE) +   
      geom_point(position = position_jitter(seed = 1, width = 0.2), alpha = 0.05, show.legend = FALSE) +
      geom_signif(data = df, comparisons = list( c('CPGs', '1K Random Genes'), c('CPGs', 'Somatic Drivers'), c('Somatic Drivers', '1K Random Genes') ), map_signif_level = T) +
      ggtitle(type) + 
      theme_update(plot.title = element_text(hjust = 0.5, size = 14)) +
      xlab(xaxis_title) +
      ylab(yaxis_title) +
      coord_cartesian(ylim = c(0, 0.006)) +
      stat_summary(fun = "mean", geom = "crossbar", aes(color = "Mean"), na.rm = F, fatten = 1.5, width = 0.75, alpha = 0.4) +
      stat_summary(fun ="median", geom = "crossbar", aes(color = "Median"), na.rm = F, fatten = 1.5, width = 0.4, alpha = 0.4) +
      scale_color_manual("Statistics", values = c(Mean = "red", Median = "black")) +
      scale_fill_brewer(palette = 'BuPu') +
      theme_update(axis.title.x = element_text(size = 18, vjust = -0.8, margin = margin(t = 0, r = 0, b = 10, l = 0)),
            axis.text.x = element_text(size = 11, hjust = 0.6, vjust = 0.6),
            axis.text.y = element_text(size = 11),
            axis.title.y = element_text(size = 18, vjust = 0.5, margin = margin(t = 0, r = 15, b = 0, l = 2.5))) 
    return(plot)
  }
  do_violins(mutfreq_combined, patho_class)
})

# leg <- cowplot::get_legend(pancancer_boxplots[[1]] + theme(legend.position = "top"))
# leg_2 <- leg$grobs[[2]]
# boxplots <- cowplot::plot_grid(leg_2)/ (pancancer_boxplots[[1]] + pancancer_boxplots[[2]] + pancancer_boxplots[[3]] + pancancer_boxplots[[4]]) +
#   plot_annotation(title = "Pancancer Germline Variant Frequency across Pathogenic Classes") +
#   plot_layout(heights = unit(c(1, 15), "cm")) +
#   plot_layout(guides = "collect") & theme(legend.position = "none", plot.title = element_text(size = 21))
# 
# pgrid <- cowplot::plot_grid(pancancer_boxplots[[1]], pancancer_boxplots[[2]] , pancancer_boxplots[[3]] , pancancer_boxplots[[4]], ncol = 4)
# p <- cowplot::plot_grid(pgrid, leg_2, ncol = 2, rel_widths = c(1, .1))

ggsave(plot = wrap_plots(pancancer_boxplots, guides = "collect") + plot_layout(ncol = 4),
         #plot_annotation(title = "Pancancer Germline Variant Frequency across Pathogenic Classes"),
       path = '/local/mmoradiellos/TFM_code/3_Germline_data/', filename = 'pancancer_germline_frequency_across_pathogenic_class_666.svg',
       width = 20.4, height = 4.59, units = "in", dpi = 300, device = 'svg')

# Save each resulting table
need <- F
if (need){
  write_tsv(mutfreq_all_pan, file = 'out_dir_input_tables_march_25/germ_var_freq_across_cancer_types_ClinvarDelmisPlof.tsv')
  write_tsv(mutfreq_plof_pan, file = 'out_dir_input_tables_march_25/germ_var_freq_across_cancer_types_Plof.tsv')
  write_tsv(mutfreq_clinvar_pan, file = 'out_dir_input_tables_march_25/germ_var_freq_across_cancer_types_Clinvar.tsv')
  write_tsv(mutfreq_delmis_pan, file = 'out_dir_input_tables_march_25/germ_var_freq_across_cancer_types_Delmis.tsv')
}
  
#################
# Wilcoxon Test #
#################
clinvar_df <- read_tsv('latest_germline_pipeline_output/germline_pancancer_freq_table_ClinVar_all_test_genesets.tsv')
plof_df <- read_tsv(file = 'latest_germline_pipeline_output/germline_pancancer_freq_table_pLoF_all_test_genesets.tsv')
delmis_df <- read_tsv(file = 'latest_germline_pipeline_output/germline_pancancer_freq_table_DelMis_all_test_genesets.tsv')
all_types_df <- read_tsv(file = 'latest_germline_pipeline_output/germline_pancancer_freq_table_ALLclass_all_test_genesets.tsv')


clinvar_CPG_cancer_genes <- wilcox.test(filter(clinvar_df, class == "Somatic Drivers")$Pancancer,
                                        filter(clinvar_df, class == "CPGs")$Pancancer)
clinvar_cancer_gene_random_genes <- wilcox.test(filter(clinvar_df, class == "Somatic Drivers")$Pancancer,
                                                filter(clinvar_df, class == "1K Random Genes")$Pancancer)
clinvar_CPG_random_genes <- wilcox.test(filter(clinvar_df, class == "1K Random Genes")$Pancancer,
                                        filter(clinvar_df, class == "CPGs")$Pancancer)


plof_CPG_cancer_genes <- wilcox.test(filter(plof_df, class == "Somatic Drivers")$Pancancer,
                                     filter(plof_df, class == "CPGs")$Pancancer)
plof_cancer_gene_random_genes <- wilcox.test(filter(plof_df, class == "Somatic Drivers")$Pancancer,
                                             filter(plof_df, class == "1K Random Genes")$Pancancer)
plof_CPG_random_genes <- wilcox.test(filter(plof_df, class == "CPGs")$Pancancer,
                                     filter(plof_df, class == "1K Random Genes")$Pancancer)

delmis_CPG_cancer_genes <- wilcox.test(filter(delmis_df, class == "Somatic Drivers")$Pancancer,
                                       filter(delmis_df, class == "CPGs")$Pancancer)
delmis_cancer_gene_random_genes <- wilcox.test(filter(delmis_df, class == "Somatic Drivers")$Pancancer,
                                               filter(delmis_df, class == "1K Random Genes")$Pancancer)
delmis_CPG_random_genes <- wilcox.test(filter(delmis_df, class == "1K Random Genes")$Pancancer,
                                       filter(delmis_df, class == "CPGs")$Pancancer)

all_types_CPG_cancer_genes <- wilcox.test(filter(all_types_df, class == "Somatic Drivers")$Pancancer,
                                          filter(all_types_df, class == "CPGs")$Pancancer)
all_types_cancer_gene_random_genes <- wilcox.test(filter(all_types_df, class == "Somatic Drivers")$Pancancer,
                                                  filter(all_types_df, class == "1K Random Genes")$Pancancer)
all_types_CPG_random_genes <- wilcox.test(filter(all_types_df, class == "1K Random Genes")$Pancancer,
                                          filter(all_types_df, class == "CPGs")$Pancancer)


p.values_dataframe <- data.frame("ClinVar" = c(clinvar_CPG_cancer_genes$p.value,
                                               clinvar_CPG_random_genes$p.value,
                                               clinvar_cancer_gene_random_genes$p.value),
                                 "pLoF" = c(plof_CPG_cancer_genes$p.value,
                                            plof_CPG_random_genes$p.value,
                                            plof_cancer_gene_random_genes$p.value),
                                 "DelMis" = c(delmis_CPG_cancer_genes$p.value,
                                              delmis_CPG_random_genes$p.value,
                                              delmis_cancer_gene_random_genes$p.value),
                                 "All types" = c(all_types_CPG_cancer_genes$p.value,
                                                 all_types_CPG_random_genes$p.value,
                                                 all_types_cancer_gene_random_genes$p.value))
row.names(p.values_dataframe) <- c("CPGs_vs_cancer_genes", "CPGs_vs_random_genes", "cancer_genes_vs_random_genes")
export_formattable( f = formattable(p.values_dataframe), file = 'latest_germline_pipeline_output//march25_wilcox_666.png')


###########################################
# LOH and GAM combining for Linear Models #
###########################################
# Maybe I can try to create the GAM with samples instead of patients,
# but nonetheless I'll try to do both just to leave them here

#######  All of the upcoming code will be used to combine LOH Binary table with GAM matrices ########
loh_binary_table_possibly_wrong <- read.table('../1_boolean_table_and_loh_events_count/loh_counts_output/all_data_by_genes_whitelisted_cnl_booltable.tsv',
                               sep = "\t", row.names = 1, header = TRUE) # This is the right one, the transposed isn't really useful
loh_binary_table <- read.table('input_tables_march_25/loh_cpgs_duppatients.tsv', sep = "\t", header = TRUE) # This was the one used to check for common patient barcodes and gene symbols; contains duplicates
loh_binary_table_gamsize <- read.table('all_data_by_genes_whitelisted_loh_booltable_shared_gam_dimensions.tsv', sep = "\t", row.names = 1, header = TRUE)

gam_plof_sampleID <- gam_obtainer_sampleID(plof_filtered_annotable)
gam_delmis_sampleID <- gam_obtainer_sampleID(delmis_filtered_annotable)
gam_clinvar_sampleID <- gam_obtainer_sampleID(clinvar_filtered_annotable)
gam_allclasses_sampleID <- gam_obtainer_sampleID(all_classes_annotable)

gam_plof_sampleID_t <- t(gam_plof_sampleID)
gam_delmis_sampleID_t <- t(gam_delmis_sampleID)
gam_clinvar_sampleID_t <- t(gam_clinvar_sampleID)
gam_allclasses_sampleID_t <- t(gam_allclasses_sampleID)

##############
# Patient level for normal LOH table and transposed GAMs
loh_binary_table_patient <- loh_binary_table
rownames(loh_binary_table_patient) <- gsub(x = rownames(loh_binary_table_patient), pattern = "-...-...-....-..", replacement = "")
loh_binary_table_patient_cpgs <- loh_binary_table_patient[, colnames(loh_binary_table_patient) %in% CPGs_list_available]

gam_plof_patientID_t <- as.data.frame.matrix(gam_plof_sampleID_t) # To go from table to df while keeping the structure
rownames(gam_plof_patientID_t) <- gsub(x = rownames(gam_plof_patientID_t), pattern = "-...-...-....-..", replacement = "")
gam_plof_patientID_t_cpgs <- gam_plof_patientID_t[, colnames(gam_plof_patientID_t) %in% CPGs_list_available]

gam_delmis_patientID_t <- as.data.frame.matrix(gam_delmis_sampleID_t)
rownames(gam_delmis_patientID_t) <- gsub(x = rownames(gam_delmis_patientID_t), pattern = "-...-...-....-..", replacement = "")
gam_delmis_patientID_t_cpgs <- gam_delmis_patientID_t[, colnames(gam_delmis_patientID_t) %in% CPGs_list_available]

gam_clinvar_patientID_t <- as.data.frame.matrix(gam_clinvar_sampleID_t)
rownames(gam_clinvar_patientID_t) <- gsub(x = rownames(gam_clinvar_patientID_t), pattern = "-...-...-....-..", replacement = "")
gam_clinvar_patientID_t_cpgs <- gam_clinvar_patientID_t[, colnames(gam_clinvar_patientID_t) %in% CPGs_list_available]

gam_allclasses_patientID_t <- as.data.frame.matrix(gam_allclasses_sampleID_t)
rownames(gam_allclasses_patientID_t) <- gsub(x = rownames(gam_allclasses_patientID_t), pattern = "-...-...-....-..", replacement = "")
gam_allclasses_patientID_t_cpgs <- gam_allclasses_patientID_t[, colnames(gam_allclasses_patientID_t) %in% CPGs_list_available]

write_tsv(data.frame(patients = rownames(loh_binary_table_patient_cpgs), loh_binary_table_patient_cpgs),
          file = 'out_dir_input_tables_march_25/gam+loh_matrices/loh_binary_table_patient_cpgs_gamdimensions.tsv')

loh_binary_table_patient_cpgs <- read_tsv('out_dir_input_tables_march_25/gam+loh_matrices/loh_binary_table_patient_cpgs_gamdimensions.tsv')

##################################################
# Testing the alternated combination of both table

# First we make some subsets of our data
loh_binary_table_patient_sub <- 
  loh_binary_table_patient[c('TCGA-02-0003', 'TCGA-02-0033', 'TCGA-02-0047', 'TCGA-02-0055', 'TCGA-02-2466'),
                         c('A1BG', 'A1CF', 'A2M', 'A2ML1', 'A3GALT2')]

gam_allclasses_patientID_t_sub <- (gam_allclasses_patientID_t[1:5 , c('A1BG', 'A1CF', 'A2M', 'A2ML1', 'A3GALT2')])

# Code to merge LOH and GAM matrices into one for the Linear Models
test_merge <- merge(gam_allclasses_patientID_t_sub, loh_binary_table_patient_sub, by = 0) # Merges by rowname, keeping all column names with '.x' and '.y' for duplicates
test_merge <- test_merge[c(1, order(sub("\\.[xy]$", "", names(test_merge)[-1])) + 1)]  # Reorders data frame so columns for the same genes are next to each other
test_merge <- test_merge %>% select( matches( ".names|\\.x|\\.y")  ) # Does not keep columns that were unique to one of the tables
colnames(test_merge) <- str_replace_all(string = colnames(test_merge), c("Row.names" = "Patient.Barcode" , ".x" =  ".germline_mut", ".y" = ".loh_event")) # Rename columns according to matrix of origin


##################################################
# Retrieving results for GAM+LOH CPGs merges at cluster, compute mut. and loh freqs.
all_gamloh <- read.table('out_dir_input_tables_march_25/gam+loh_matrices/matrices_merges_cluster/clinvar_loh_cpgs_merged.tsv', sep = "\t", row.names = 1, header = TRUE)

loh_colums <- grep(".loh_event", names(all_gamloh))
all_gamloh_justloh <- all_gamloh %>% select(loh_colums)
round(mean(as.matrix(all_gamloh_justloh)) * 100, 2)

germmut_colums <- grep(".germline_mut", names(all_gamloh))
all_gamloh_justgermmut<- all_gamloh %>% select(germmut_colums)
mean(as.matrix(all_gamloh_justgermmut))

#--------------------------#
# Test and Issues to check #
#

# These are various test and issues that I've been checking during the process

#-| Collect Problem Genes |-#
## Obtain the max. freq of each for all cancer types,
## the extract the name of those genes over 0.1 to
## check in the original table by other means
casestudy_vartable <- delmis_filtered_annotable %>% 
  filter(consensus_gene_symbol %in% rownames(data.frame(mutfreq_delmis) %>% 
                                               mutate(gene_max_freq = do.call(pmax, c(select(., colnames(.))))) %>% 
                                               filter(gene_max_freq > 0.1)))

af_issue <- readxl::read_xlsx(path ='./issue_all_pathoclass_variant_annotation_table_check.xlsx')
ggplot(data = af_issue, aes(x = as.numeric(AF_gnomad_WES), y = AF_TCGA)) +
  geom_point() 
  #coord_cartesian(ylim = c(0, 0.002), xlim = c(0, 0.0011))
# Entry with AF > 0.1
movida <- af_issue %>% filter(as.numeric(AF_gnomad_WES) >= 0.1)

# Entries without germline code in SAMPLE
all_classes_annotable %>% filter(SAMPLE %in% grep(pattern = 'TCGA-..-....-01.-*', colnames(gam_delmis_sampleID), value = T))

# Counts of different Variant-IDs and samples with it
af_TCGA_counts <- all_classes_annotable %>%
  count(Otherinfo3, consensus_gene_symbol, AF_TCGA) %>% 
  arrange(consensus_gene_symbol, desc(n))

############################################
# Retrieving genes due to different naming
############################################
loh_symbols <- colnames(loh_binary_table)
gam_symbols_in_loh <- colnames(gam_plof_sampleID_t)[(colnames(gam_plof_sampleID_t) %in% loh_symbols) ]
gam_symbols_not_in_loh <- colnames(gam_plof_sampleID_t)[! ((colnames(gam_plof_sampleID_t) %in% loh_symbols)) ]
loh_symbols_not_in_gam <- loh_symbols[! (loh_symbols %in% colnames(gam_clinvar_sampleID_t)) ]

library(EnsDb.Hsapiens.v79)
dbCon <- org.Hs.eg_dbconn()
# write your SQL query
sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
# execute the query on the database
aliasSymbol <- dbGetQuery(dbCon, sqlQuery)
# subset to get your results
gam_symbols_not_in_loh_aliases <- aliasSymbol[which(aliasSymbol[, 5] %in% gam_symbols_not_in_loh), c(2, 5)]
gam_loh_translated_matches <- dplyr::filter(gam_symbols_not_in_loh_aliases, alias_symbol %in% loh_symbols_not_in_gam)

# Subset of LOH that was missing but now retrieved
loh_to_translate <- loh_binary_table[, gam_loh_translated_matches$alias_symbol]
colnames(loh_to_translate) <- gam_loh_translated_matches$symbol

# Subset of LOH that already had matches with GAM
loh_not_to_translate <- loh_binary_table[, gam_symbols_in_loh]

gazpacho <- cbind(loh_not_to_translate, loh_to_translate)
rownames(gazpacho) <- gsub(x = rownames(gazpacho), pattern = "-...-...-....-..", replacement = "")

gam_patients_in_loh <- rownames(gam_plof_sampleID_t)[(rownames(gam_plof_sampleID_t) %in% rownames(loh_binary_table_patient)) ]
gazpacho_shared_patients <- gazpacho[gam_patients_in_loh, ]
write_tsv(data.frame(patients = rownames(gazpacho_shared_patients), gazpacho_shared_patients),
            file = 'all_data_by_genes_whitelisted_loh_booltable_shared_gam_dimensions.tsv')


test_merge <- merge(gam_plof_sampleID_t, gazpacho_shared_patients, by = 0) # Merges by rowname, keeping all column names with '.x' and '.y' for duplicates
test_merge <- test_merge[c(1, order(sub("\\.[xy]$", "", names(test_merge)[-1])) + 1)]  # Reorders data frame so columns for the same genes are next to each other
test_merge <- test_merge %>% select( matches( ".names|\\.x|\\.y")  ) # Does not keep columns that were unique to one of the tables
colnames(test_merge) <- str_replace_all(string = colnames(test_merge), c("Row.names" = "Patient.Barcode" , ".x" =  ".mut", ".y" = ".loh")) # Rename columns according to matrix of origin

###########
# Issue on duplicated LOH Patients when comparing LOH mean freq. between now and October's task
possibly_loh_filtered <- read.table('../1_boolean_table_and_loh_events_count/loh_counts_output/all_data_by_genes_whitelisted_cnl_booltable.tsv', sep = '\t', row.names = 1, header = TRUE)
complete_loh <- read.table('../1_boolean_table_and_loh_events_count/second_task/loh_counts_output/loh_matrix_complete.tsv', sep = '\t', row.names = 1, header = TRUE)

loh_cpgs_allsamples <- read.table('input_tables_march_25/loh_cpgs_allsamples.tsv', sep = '\t', header = TRUE, row.names = 1) # LOH for CPGs, has all samples (repeated patients)

# Extract which patients have multiple samples
df_samples <- data.frame("sample" = row.names(loh_cpgs_allsamples), "patients" = substr(rownames(loh_cpgs_allsamples), start = 1, stop = 12)) 
dup <- names(table(df_samples$patients))[table(df_samples$patients) >= 2]
df_samples_dup <- filter(df_samples, patients %in% dup) %>% order(~patients, decreasing = TRUE) # Filter and next "group" by patients
df_samples_dup[order(df_samples_dup$patients), ]
mean(grepl(pattern = "-01$", df_samples_dup$sample)) # Check how many samples per "duplicated patient" we may have

# From all repeated patients extracts samples that mostly be related to primary tumor
dup_samples_to_keep <- sapply(unique(df_samples_dup$patients), function(patient){
  return( sort( filter(df_samples_dup, df_samples_dup$patients == patient)$sample)[1] )
})
# Extract those samples that were not duplicates and combine them with the previously selected
no_dup <- names(table(df_samples$patients))[table(df_samples$patients) == 1]
df_samples_no_dup <- filter(df_samples, patients %in% no_dup)
samples_filter <- c(df_samples_no_dup$sample, unname(dup_samples_to_keep))

# Filter the complete LOH matrix with the selected patients (common to the GAMs)
loh_cpgs_nodup_samples <- loh_cpgs_allsamples[rownames(loh_cpgs_allsamples) %in% samples_filter, ]
loh_cpgs_nodup_patients <- loh_cpgs_nodup_samples
rownames(loh_cpgs_nodup_patients) <- substr(rownames(loh_cpgs_nodup_patients), start = 1, stop = 12) # Change to patient barcode

# Check which patients have both matrices in common to create CPGs matrices faster to merge
common_patients_gam_loh.nodup <- rownames(gam_plof_patientID_t)[rownames(gam_plof_patientID_t) %in% rownames(loh_cpgs_nodup_patients)]
loh_cpgs_common_patient_with_gam <- loh_cpgs_nodup_patients[common_patients_gam_loh.nodup, ]
gam_plof_patientID_t_cpgs_common_patient_with_loh <- gam_plof_patientID_t_cpgs[common_patients_gam_loh.nodup, ]
gam_delmis_patientID_t_cpgs_common_patient_with_loh <- gam_delmis_patientID_t_cpgs[common_patients_gam_loh.nodup, ]
gam_clinvar_patientID_t_cpgs_common_patient_with_loh <- gam_clinvar_patientID_t_cpgs[common_patients_gam_loh.nodup, ]
gam_allclasses_patientID_t_cpgs_common_patient_with_loh <- gam_allclasses_patientID_t_cpgs[common_patients_gam_loh.nodup, ]

write_tsv(data.frame(patients = rownames(loh_cpgs_common_patient_with_gam), loh_cpgs_common_patient_with_gam),
          file = 'out_dir_input_tables_march_25/gam+loh_matrices/loh_cpgs_common_patient_with_gam.tsv')
