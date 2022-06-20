#########################
# Perpetrated by Manuel #        22-02-22 (New MergeID Criteria)
#########################

## This script takes the result of a file annotated by VEP browser
## with the predictors' score CADD, LoFtool and SpliceAI and performs 
## various filters and graphical comparisons for us to decide on a  
## suitable threshold for those tools.

# Once the filtering steps are validated, they 
# will added to the final snakemake pipeline

# Loads all R packages that we will need
packages <- c("tidyverse","ggplot2", "RColorBrewer", "patchwork", 
              "formattable", "webshot", "htmltools", "flextable")
lapply(packages, require, character.only = TRUE)


#####################################
# {|-Functions that will be used-|} #
#####################################

#---------------#
#-DATA HANDLING-#
#---------------#

vep_output_relevantcols <- function(tsv_filename){
  df_vep_annot <- read_tsv(file = tsv_filename)
  df_vep_annot_relevantcols <- df_vep_annot[, c('#Uploaded_variation', 'Gene', 'Feature', 
                                                    'LoFtool', 'SpliceAI_pred_DP_AG', 
                                                    'SpliceAI_pred_DP_AL', 'SpliceAI_pred_DP_DG',
                                                    'SpliceAI_pred_DP_DL', 'SpliceAI_pred_DS_AG',
                                                    'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG',
                                                    'SpliceAI_pred_DS_DL', 'SpliceAI_pred_SYMBOL')]
  #colnames(df_vep_annot_relevantcols)[c(1, 4)] <- c('Uploaded_variation', 'CADD_PHRED_VEP') # Remove problematic character
  colnames(df_vep_annot_relevantcols)[1] <- 'Uploaded_variation' # Remove problematic character
  df_vep_annot_relevantcols$Feature <- str_replace(df_vep_annot_relevantcols$Feature, pattern = '.[0-9]+$', replacement = '') # Remove version number, not useful
  return(df_vep_annot_relevantcols %>% distinct()) # Features where returning duplicated rows...
}

last_exon_checker <- function(dataframe){
  # Obtains a subset of the dataframe keeping those rows where we have
  # NAs for the 'EXON' entry (in which exon on the gene is the variant)
  # or when that last bit its TRUE
  ## It's a bit naive as we end up with a lot of NAs but it was easy to do
  dataframe$in_last_exon <- sapply(1:nrow(dataframe), function(x){
    if( dataframe[x, 'EXON'] == '-' ){
      return(NA)
    } else {
      ifelse(length(unique(strsplit(dataframe[x, 'EXON'], split = '/')[[1]])) == 1, TRUE, FALSE)
    }
  })
  return(dataframe) # Maybe we want to keep all entries if they are not in last exon
  #return(dataframe %>% 
  #         filter( is.na(in_last_exon) | (in_last_exon == TRUE)   ) )
}


clinvar_class_sapply <- function(dataframe){
  # Adds a column to the dataframe checking whether the clinical #
  # significance falls under our list of "pathogenic" tags       #
  dataframe$clinvar_patho <- sapply(1:nrow(dataframe), function(x){
    clinvar_vector <- strsplit(dataframe[x, 'CLIN_SIG'], split = ',')[[1]]
    outvalue <- any(clinvar_vector %in% c('likely_pathogenic', 'risk_factor', 'pathogenic', 'pathogenic/likely_pathogenic', '_risk_factor'))
    return(outvalue)
  })
  return(dataframe)
}


spliceai_damag_sapply <- function(dataframe, spliceai_thr = 0.8){
  # First obtains the maximum SpliceAI score from each entry
  # as it is the one to check as recommended by the authors, 
  # then checks if its over the threshold
  dataframe <- dataframe %>% 
    mutate(spliceai_max = do.call(pmax, c(select(., c("SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL")))))
    #mutate(spliceai_max = do.call(pmax, c(select(., c("SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DL"))))) # Mischan's criteria, does not make much sense to me
  dataframe$spliceai_damag <- sapply(1:nrow(dataframe), function(x){
    if ( dataframe[x , 'Func.refGene'] != 'splicing' ) {
      return(FALSE) }
    else if ( is.na(dataframe[x , 'spliceai_max']) ) {
      return(TRUE) }
    else {
      ifelse(dataframe[x , 'spliceai_max'] > spliceai_thr, TRUE, FALSE) } 
  })
  return(dataframe)
}


cadd_damag_sapply <- function(dataframe, cadd_thr = 15){
  # Adds a column which shows whether a variants passes CADD score or not
  dataframe$cadd_damag <- sapply(1:nrow(dataframe), function(x){
    if ( dataframe[x ,'Func.refGene'] != 'exonic' ){
      return(FALSE) }
    else if ( is.na(dataframe[x ,'CADD_phred']) ){
      return(NA) }
    else {
      ifelse(dataframe[x ,'CADD_phred'] >= cadd_thr, TRUE, FALSE)  }
  })
  return(dataframe)
}


loftool_damag_sapply <- function(dataframe, loftool_thr = 0.334){
  # Adds a column which shows whether a variants passes LoFtool test or not
  dataframe$loftool_damag <- sapply(1:nrow(dataframe), function(x){
    if ( dataframe[x ,'Func.refGene'] != 'exonic' ){
      return(FALSE) }
    else if ( is.na(dataframe[x ,'LoFtool']) ){
      return(NA) }
    else {
      ifelse(dataframe[x ,'LoFtool'] < loftool_thr, TRUE, FALSE) }
  })
  return(dataframe)
}


cadd_loftool_damag_sapply <- function(dataframe, cadd_thr = 0.8, loftool_thr = 0.334){
  # Adds a column where it test whether an exonic function
  # passses CADD threshold and if not, LoFtool's 
  dataframe$cadd_loftool_damag <- sapply(1:nrow(dataframe), function(x){
    if ( dataframe[x ,'Func.refGene'] != 'exonic' ){
      return(FALSE) }
    else if ( is.na(dataframe[x ,'CADD_phred']) ){
      return(NA) }
    else if ( dataframe[x ,'CADD_phred'] >= cadd_thr ){
      return(TRUE) }
    else if ( is.na(dataframe[x ,'LoFtool']) ){
      return(NA) }
    else{
      ifelse(dataframe[x ,'LoFtool'] < loftool_thr, TRUE, FALSE) }
  })
  return(dataframe)
}


wholecriteria_damag_sapply <- function(dataframe, spliceai_thr = 0.8, cadd_thr = 15,
                                       loftool_thr = 0.334){
  # Applies our current criteria in the correct order: Checks if variant is 
  # located in terminal exon or if we have NA for that entry,  then checks if
  # it is splicing (if so, SpliceAI) or exonic (if so, CADD_phred or LoFtool)
  #
  # This criteria just removes variants located in the terminal exon that
  # may not truly be damaging, and keeps the rest (NAs for some predictor
  # scores, NAs and FALSE for last terminal)
  dataframe$whole_criteria_keep <- sapply(1:nrow(dataframe), function(x){
    if ( is.na(dataframe[x, 'in_last_exon']) | dataframe[x, 'in_last_exon'] == FALSE ){
      return(TRUE)
    }
    else {
      if ( dataframe[x, 'Func.refGene'] == 'splicing' ){
        if ( is.na(dataframe[x, 'spliceai_max']) ){
          return(TRUE) }  
        else {
          ifelse(dataframe[x, 'spliceai_max'] >= spliceai_thr, TRUE,  FALSE) } }
      else if ( dataframe[x, 'Func.refGene'] == 'exonic' ){
        if ( !(is.na(dataframe[x, 'CADD_phred'])) ){
          ifelse( dataframe[x, 'CADD_phred'] >= cadd_thr, TRUE, FALSE) }
        else if ( !(is.na(dataframe[x, 'LoFtool'])) ){
          ifelse(dataframe[x ,'LoFtool'] < loftool_thr, TRUE, FALSE) }
        else {
          return(TRUE) } }
    }
  })
  return(dataframe)
}


tcgagnomadfilter_paste_sapply <- function(dataframe){
  # Creates a column pasting two results to check for possible biases in scoring
  dataframe$filter_comb <- sapply(1:nrow(dataframe), function(x){
    if ( is.na(dataframe[x, 'FILTER_gnomad']) ){
      return( paste0(dataframe[x, 'FILTER_ID'],'-NULL') )  }
    else {
      return( paste0(dataframe[x, 'FILTER_ID'],'-',dataframe[x, 'FILTER_gnomad']) )   }
  })
  return(dataframe)
}


dataframe_vep_filterer <- function(dataframe, cadd_thr = 15, loftool_thr = 0.334,
                                   spliceai_thr = 0.8, cpgs){
  # Takes a data frame and adds various boolean comparison of tests for # 
  # latter filtering, although some of them will be only for variants   #
  # in the terminal exon                                                #
  
  # 0. # Checks if variants are in last exon or not
  dataframe_filt <- last_exon_checker(dataframe)
  # 1. # Adds a column which shows whether splicing variants passes SpliceAI test or not
  dataframe_filt <- spliceai_damag_sapply(dataframe_filt)
  # 2. # Adds a column which shows whether exonic variants passes CADD test or not
  dataframe_filt <- cadd_damag_sapply(dataframe_filt)
  # 3. # Adds a column which shows whether exonic variants passes LoFtool test or not
  dataframe_filt <- loftool_damag_sapply(dataframe_filt)
  # 4. # Adds a column which shows whether exonic variants passes our criteria or not
  dataframe_filt <- cadd_loftool_damag_sapply(dataframe_filt)
  # 5. # Applies our filtering criteria which tests if any variant
  #      passes our combined criteria of SpliceAI, CADD, LoFtool
  dataframe_filt <- wholecriteria_damag_sapply(dataframe_filt)
  # 6. # Checks whether the gene that has that variant is a CPG
  dataframe_filt$cpg_or_not <- dataframe_filt$consensus_gene_symbol %in% cpgs
  # 7. # Adds a column which shows if clinical significance is pathogenic 
  dataframe_filt <- clinvar_class_sapply(dataframe_filt)
  # 8. # Adds a column with the combination of TCGA and gnomAD filters
  dataframe_filt <- tcgagnomadfilter_paste_sapply(dataframe_filt)
  return(dataframe_filt)
}


dataframe_vep_filterer_snakemake <- function(dataframe){
  # Takes a data frame and performs various boolean  tests for # 
  # filtering, although some of them will be only for variants #
  # in the terminal exon                                       #
  
  cols_to_remove <-c('Uploaded_variation', 'Gene.y', 'Feature.y', 
                     'LoFtool', 'SpliceAI_pred_DP_AG', 'SpliceAI_pred_DP_AL',
                     'SpliceAI_pred_DP_DG', 'SpliceAI_pred_DP_DL', 'SpliceAI_pred_DS_AG',
                     'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL',
                     'SpliceAI_pred_SYMBOL', 'in_last_exon', 'spliceai_max',
                     'whole_criteria_keep', 'pasteID_genetrans')
  
  # 0. # Changes some column types
  dataframe <- important_cols_as_num(dataframe)
  # 1. # Checks if variants are in last exon or not
  dataframe<- last_exon_checker(dataframe)
  # 2. # Obtains the maximum SpliceAI score from each entry as #
  #      it is the one to check as recommended by the authors  #
  dataframe <- dataframe %>% 
    mutate(spliceai_max = do.call(pmax, c(select(., c("SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL")))))
  # 3. # Applies our filtering criteria which tests if any variant in terminal
  #      exon passes our combined criteria of SpliceAI, CADD, LoFtool
  dataframe <- wholecriteria_damag_sapply(dataframe)
  ## I don't understand why the previous function  ##
  ## returns a list and sometimes raises errors:   ##
  #dataframe$whole_criteria_keep <- unlist(dataframe$whole_criteria_keep)
  # 4. # Removes some variant types that are not of interest
  dataframe <- dataframe %>% filter(Func.refGene != 'exonic;splicing')
  # 5. # Filters out variants that do not pass our criteria
  dataframe_filt <- (dataframe %>% filter(whole_criteria_keep != FALSE))
  # 6. # Removes columns added to do the filtering as they are no further useful
  dataframe_filt <- select(dataframe_filt, -cols_to_remove)
  return(dataframe_filt)
}

snakemake_merger <- function(table_to_match, vep_annotated){
  # Work on unique variants as VarID + Gene + Transcript, create uniqueID based on that
  rare_pLOF_unique_vars$pasteID_genetrans <-
    paste(rare_pLOF_unique_vars$Otherinfo3, rare_pLOF_unique_vars$Gene, rare_pLOF_unique_vars$Feature, sep = "_")
  
  # Read VEP-Browser output keeping relevant columns and create uniqueID as before
  plof_vep_annot_nodup <- vep_output_relevantcols(vep_file)
  plof_vep_annot_nodup$pasteID_genetrans <- 
    paste(plof_vep_annot_nodup$Uploaded_variation, plof_vep_annot_nodup$Gene.VEP, plof_vep_annot_nodup$Feature.VEP, sep = "_")
  
  # Match both files and keep all rows of our input that do not match
  rareplof_matched <- merge(rare_pLOF_unique_vars, plof_vep_annot_nodup, by = 'pasteID_genetrans', all.x = TRUE)
  
  # Apply CADD+LoFtool/SpliceAI filters on terminal exon variants and write results
  rareplof_matched_filt <- dataframe_vep_filterer(rareplof_matched)
  write.table(rareplof_matched_filt, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}
  


gene_grouping <- function(dataframe, cpg_list = CPGs_list_available,
                          ccg_list = ccg_list_available_not_CPG, 
                          random_gene_set = random_genes_1K){
  # For some comparison we want to divide each gene into individual groupings
  # This adds a column according to each type they belong to
  dataframe$gene_group <- sapply(1:nrow(dataframe), function(x){
    if ( dataframe[x, 'consensus_gene_symbol'] %in% cpg_list ){
      return(factor('CPG')) }
    else if ( dataframe[x, 'consensus_gene_symbol'] %in% ccg_list ){
      return(factor('CancerConsensusGene')) }
    else if ( dataframe[x, 'consensus_gene_symbol'] %in% random_gene_set ){
      return(factor('1KRandomGeneSet')) }
  })
  return(dataframe)
}


important_cols_as_num <- function(dataframe){
  # Some columns are read as characters due to having '-' as NAs,
  # so we need to change their value types to perform the numeric
  # filters
  dataframe[, c("CADD_phred", "LoFtool", "SpliceAI_pred_DS_AG", 
                         "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", 
                         "SpliceAI_pred_DS_DL")] <- 
    sapply(dataframe[, c("CADD_phred","LoFtool", "SpliceAI_pred_DS_AG", 
                                  "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", 
                                  "SpliceAI_pred_DS_DL")], as.numeric)
  return(dataframe)
}


obtain_freq_table <- function(dataframe){
  # Obtain a table with the number of mutations of each gene per patient,
  # disclosing the type of gene according to three defined groups
  gam_matrix <- as.matrix(table( factor(dataframe$consensus_gene_symbol, levels = c(CPGs_list_available, ccg_list_available_not_CPG, random_genes_1K)),
                                 factor(dataframe$bcr_patient_barcode, levels = all_patients)   ))
  gam_matrix[gam_matrix > 2] <- 1 
  freqtable <- data.frame('Gene.freq' = rowMeans(gam_matrix))
  freqtable$consensus_gene_symbol <- rownames(freqtable)
  return(distinct(merge(freqtable, dataframe[, c('consensus_gene_symbol', 'gene_group')], by = 'consensus_gene_symbol')))
}


export_formattable <- function(f, file, width = "100%", height = NULL, 
                               background = "white", delay = 0.2){
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

violin_plotter_loftool_spliceai_scores <- function(dataframe, scores, title, y_axis_title, thres, thres_legend){
  # Violin plot that compares a predictor score for exonic #
  # and splicing variants returning the mean and median    #
  violin <- ggplot(dataframe, aes(x = Func.refGene, y = get(scores), colour = Func.refGene)) +
    geom_violin(aes(fill= Func.refGene), na.rm = FALSE, lwd = 0.2) +
    scale_fill_manual("Variant Type", values = c(exonic = "orange", splicing = "lightblue")) +
    geom_hline(aes(yintercept = thres, linetype = thres_legend), color = 'red') +
    scale_linetype_manual(name = '', values = 2) +
    stat_summary(fun = "mean", geom = "point", aes(color = "Mean"), na.rm = TRUE) +
    stat_summary(fun ="median", geom = "point", aes(color = "Median"), na.rm = TRUE) +
    scale_color_manual("Statistics", values = c(Mean = "purple", Median = "black")) +
    ggtitle(title) +
    theme_update(plot.title = element_text(hjust = 0.5, size = 15)) +
    xlab('Variant Type (Func.refGene)') +
    ylab(y_axis_title) +
    theme_update(plot.title = element_text(hjust = 0.5, size = 15)) 
  print(violin)
}


violin_mean_mutation_freq <- function(dataframe, plottitle = '', plot_subtitle = '', 
                                  xaxis_title = 'Gene Type', yaxis_title = 'Mean Mutation Frequency per Gene'){
  ggplot(dataframe, aes(x = gene_group, y = Gene.freq, fill = gene_group)) +
    geom_violin(aes(x = gene_group), na.rm = FALSE, lwd = 0.3, alpha = 0.45) +
    geom_jitter(alpha = 0.08) +
    scale_fill_manual("Gene Type", values = c('darkolivegreen2', 'lightblue', 'plum2')) +
    stat_summary(fun = "mean", geom = "crossbar", aes(color = "Mean"), na.rm = TRUE, fatten = 2, width = 0.3) +
    stat_summary(fun ="median", geom = "crossbar", aes(color = "Median"), na.rm = TRUE, fatten = 2, width = 0.3) +
    scale_color_manual("Statistics", values = c(Mean = "blue1", Median = "red1")) +
    ggtitle(plottitle) +
    xlab(xaxis_title) +
    ylab(yaxis_title) + ylim(0, 0.016) +
    theme_update(plot.title = element_text(hjust = 0.5, vjust = -0.3, size = 12)) +
    labs(subtitle = plot_subtitle) +
    theme(axis.title.x = element_text(size = 15, vjust = -1.3),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 11),
          axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          panel.background = element_rect(fill = 'gray92')) +
    guides(fill = guide_legend(override.aes = list(linetype = 0))) # Removes marks from legend
}


barplotter_loftool_notpass <- function(dataframe, y_plot = 'Freq', title, y_axis_title){
  # Barplot that shows either the Total Counts ('Freq') or the     #
  # percentage ('Proportion') of variants that do not pass LoFtool #
  
  # Obtains freq. table 
  ptv_table <- as.data.frame(table(dataframe$ExonicFunc.refGene))
  colnames(ptv_table)[1] <- 'PTV'
  levels(ptv_table$PTV)[1] <- 'Splicing'
  
  barplot <- ggplot(data = ptv_table, aes(x = PTV, y = get(y_plot), fill = PTV)) +
    geom_bar(stat = "Identity") +
    geom_text(aes(label=get(y_plot)), vjust=-0.3, size=4.5) +
    ggtitle(title) +
    xlab('Protein-Truncating Variant Type') +
    ylab(y_axis_title) +
    theme_update(plot.title = element_text(hjust = 0.5, size = 15)) +
    theme(axis.title.x = element_text(size = 12, vjust = 1),
          axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11),
          axis.title.y = element_text(size = 12, vjust = 1),
          legend.position = "none")
  print(barplot)
}


density_plotter <- function(dataframe, x_plot = 'gnomAD_AF', title, 
                            x_axis_title = 'gnomAD Allele Frequency',
                            default_colors = FALSE){
  densiplot <- ggplot(dataframe, aes(x = as.numeric(get(x_plot)), colour = loftool_damag) ) +
    geom_density( alpha= 0.5) +
    ggtitle(title) + 
    xlab(x_axis_title) +
    ylab('Density') +
    theme_update(plot.title = element_text(hjust = 0.5, size = 15)) +
    theme(axis.title.x = element_text(size = 12, vjust = 1),
          axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11),
          axis.title.y = element_text(size = 12, vjust = 1))
  if (default_colors) {
    densiplot = densiplot + 
      scale_colour_manual('LoFtool Threshold',  labels = c ('Pass', 'Not Pass', 'NA'),
                          values = c('#F8766D', '#00BFC4'), na.value = 'green4')}
  else {
    densiplot = densiplot + 
      scale_colour_manual('LoFtool Threshold', labels = c ('Pass', 'Not Pass', 'NA'),
                          values = c('darkmagenta', 'royalblue'), na.value = 'orange2')}
  print(densiplot)
}


stack_plotter <- function(dataframe, default_colors = TRUE){
  #
  # Creates stacked bar plots with the top ten genes with most counts for each
  # category, including the rest as a comparison for possible biases
  #
  gene_pass <- as.data.frame(table(dataframe[, c('consensus_gene_symbol', 'loftool_damag')], useNA = 'always')) # Keep NAs
  loftool_damag_levels <- c('TRUE', 'FALSE', 'NA') # To plot each level
  plot_list <- lapply(1:length(loftool_damag_levels), function(filter_pass){
    if (loftool_damag_levels[filter_pass] == 'NA') {
      # To extract the top ten genes with most counts for each level
      subset_top10 <- gene_pass[order(is.na(gene_pass$loftool_damag), gene_pass$Freq, decreasing = TRUE), ]$consensus_gene_symbol[1:10]
      subtitle_level <-  'NA'
      xaxis <- ' '
      yaxis <- ' '
    }
    else {
      subset_top10 <- gene_pass[order(gene_pass$loftool_damag == loftool_damag_levels[filter_pass], gene_pass$Freq, decreasing = TRUE), ]$consensus_gene_symbol[1:10]
      if (loftool_damag_levels[filter_pass] == 'TRUE'){
        subtitle_level <- 'Pass'
        xaxis <- ' '
        yaxis <- 'Total Count'}
      else {
        subtitle_level <- 'Not Pass' 
        xaxis <- 'Gene'
        yaxis <- ' '}
    }
    to_plot <- gene_pass[gene_pass$consensus_gene_symbol %in% subset_top10, ] # subsets data to plot
    stackplot <- to_plot %>%
      # Creating factors, levels and using mapping to represent  
      # the bar plot in the order with most counts for each level
      mutate(consensus_gene_symbol = factor(to_plot$consensus_gene_symbol, levels = subset_top10, ordered = TRUE)) %>%
      ggplot(to_plot, mapping = aes(x = consensus_gene_symbol, y = Freq, fill = loftool_damag)) +
      geom_bar(position = 'stack', stat = 'identity') +
      xlab(xaxis) +
      ylab(yaxis) + ylim(0, 125) +
      labs(subtitle = paste0('Ordered by most counts in ', subtitle_level)) +
      theme_update(plot.subtitle = element_text(hjust = 0.5, size = 9)) +
      theme(axis.title.x = element_text(size = 14, vjust = 1),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 11),
            axis.title.y = element_text(size = 14, vjust = 1)) 
    if (default_colors) {
      stackplot = stackplot + 
        scale_colour_manual('LoFtool Threshold',  labels = c ('Not Pass', 'Pass', 'NA'),
                            values = c('#F8766D', '#00BFC4'), na.value = 'green4')}
    else {
      stackplot = stackplot + 
        scale_colour_manual('LoFtool Threshold', labels = c ('Not Pass', 'Pass', 'NA'),
                            values = c('darkmagenta', 'royalblue'), na.value = 'orange2')}
    return(stackplot)
  })
  wrap_plots(plot_list, guides = 'collect', tag_level = 'new') +
    plot_annotation(title = 'Top 10 Genes in each LoFtool filter level')
}


#---------------#
#-PRETTY TABLES-#
#---------------#

filters_summary_tables_old <- function(dataframe, saveTABLE = FALSE){
  
  total_entries <- nrow(dataframe) 
  
  ###################################################
  # 1. Percentages of entries with predictors score #
  ###################################################
  
  merged_table_entries_with_loftool <- nrow(dataframe[!(is.na(dataframe$loftool)), ])
  merged_table_entries_with_spliceAI <-nrow(dataframe[!(is.na(dataframe$spliceai_max)), ])
  merged_table_entries_noNAs <- nrow(dataframe[  !(is.na(dataframe$loftool)) & !(is.na(dataframe$spliceai_max)), ])
  merged_table_entries_someNAs <- nrow(dataframe[  !(is.na(dataframe$loftool)) | !(is.na(dataframe$spliceai_max)), ])
  
  total_and_missing_entries <- data.frame(num_vars = c(merged_table_entries_with_loftool,
                                                       merged_table_entries_with_spliceAI,
                                                       merged_table_entries_noNAs,
                                                       merged_table_entries_someNAs),
                                          row.names = c('With LofTool score',
                                                        'With SpliceAI score',
                                                        "With both predictors' scores",
                                                        "With any predictors' score"))
  # This creates row.names as a column with no table to make it bold with ` ` = formatter...
  total_and_missing_entries$dataframe <- rownames(total_and_missing_entries)
  rownames(total_and_missing_entries) <- c()
  total_and_missing_entries <- total_and_missing_entries[, c(2,1)]
  colnames(total_and_missing_entries) <- c(' ', '% of variants')
  
  percentage_with_info <- formattable(total_and_missing_entries, align = c('l', rep('l', ncol(total_and_missing_entries) -1)),
                                      list(` ` = formatter("span", style = ~style(display = "block", font.weight = "bold")),
                                           `% of variants` = formatter("span", x ~ percent(x / total_entries))))
  print(percentage_with_info)
  
  ################################################################
  # 2. Percentages of entries that pass the predictors threshold #
  ################################################################
  
  merged_table_loftool_damaging_nrow <- nrow(dataframe %>% filter( !is.na(loftool_damag) & (loftool_damag == TRUE) ))
  merged_table_spliceAI_damaging_nrow <- nrow(dataframe %>% filter( !is.na(spliceai_damag) & (spliceai_damag == TRUE) ))
  merged_table_loftool_AND_spliceai_damaging_nrow <- nrow(dataframe %>% filter( !is.na(loftool_damag) & (loftool_damag == TRUE) & !is.na(spliceai_damag) & (spliceai_damag == TRUE) ))
  matched_annotation_noNA_loftool_OR_spliceai_damaging_nrow <- nrow(dataframe %>% filter( !is.na(loftool_damag) & !is.na(spliceai_damag) & ( (loftool_damag == TRUE) | (spliceai_damag == TRUE)) ))
  matched_annotation_someNA_loftool_OR_spliceai_damaging_nrow <- nrow(dataframe %>% filter( (loftool_damag == TRUE) | (spliceai_damag == TRUE) ))
  
  pred_compar  <- data.frame(lof = c(merged_table_loftool_damaging_nrow, NA),
                             spai = c(merged_table_spliceAI_damaging_nrow, NA),
                             lofANDspai = c(merged_table_loftool_AND_spliceai_damaging_nrow, NA),
                             lofORspai = c(matched_annotation_noNA_loftool_OR_spliceai_damaging_nrow, 
                                           matched_annotation_someNA_loftool_OR_spliceai_damaging_nrow),
                             row.names = c('Without NAs', 'With NAs in any predictor'))
  #pred_compar[is.na(pred_compar)] = '' # Remove NAs from appearing on table if not percentage
  
  pred_compar$blank <- rownames(pred_compar)
  rownames(pred_compar) <- c()
  pred_compar <- pred_compar[, c(5, c(1:4))]
  colnames(pred_compar) <- c(' ', 'LoFtools Damaging', 'SpliceAI Damaging',
                             'LoFtools AND SpliceAI damaging', 'LoFtools OR SpliceAI damaging')
  filter_passes <- formattable(pred_compar, align = c('l', rep('l', ncol(total_and_missing_entries) -1)),
                               list(` ` = formatter("span", style = ~style(display = "block", font.weight = "bold")),
                                    `LoFtools OR SpliceAI damaging` = formatter("span", x ~ percent( x / total_entries)),
                                    `LoFtools Damaging` =  formatter("span", x ~ percent( x / total_entries)),
                                    `SpliceAI Damaging` =  formatter("span", x ~ percent( x / total_entries)),
                                    `LoFtools AND SpliceAI damaging` =  formatter("span", x ~ percent( x / total_entries))))
  print(filter_passes)
  
  cat(paste0('Number of variants for ', deparse(substitute(dataframe)) ,' that pass LoFtool: ', nrow(dataframe %>% filter(loftool_damag == TRUE)), '\n'))
  
}

filters_summary_tables <- function(dataframe, saveTABLE = FALSE){
  ## Updated version that includes CADD 
  ## score and a different table design
  
  dataframe_last_exon <- dataframe %>% filter(in_last_exon, TRUE)
  dataframe_not_keep <- dataframe %>% filter(whole_criteria_keep == FALSE)
  total_entries_merged <- nrow(dataframe) 
  
  ###################################################
  # 1. Percentages of entries with predictors score #
  ###################################################
  
  # Metrics for overall (from the total of entries)
  merged_table_exonic_entries_with_loftool <- nrow(dataframe[ (!(is.na(dataframe$LoFtool)) & dataframe$Func.refGene == 'exonic'), ])
  merged_table_exonic_entries_with_cadd <- nrow(dataframe[ (!(is.na(dataframe$CADD_phred)) & dataframe$Func.refGene == 'exonic'), ])
  merged_table_exonic_entries_with_both <- nrow(dataframe[ (!(is.na(dataframe$CADD_phred)) & !(is.na(dataframe$LoFtool)) & dataframe$Func.refGene == 'exonic'), ])
  merged_table_splicing_entries_with_spliceai <- nrow(dataframe[ (!(is.na(dataframe$spliceai_max)) & dataframe$Func.refGene == 'splicing'), ])
  
  # Metrics for variants in terminal exons (those in which we apply the criteria)
  merged_table_exonic_entries_with_loftool_le <- nrow(dataframe_last_exon[ (!(is.na(dataframe_last_exon$LoFtool)) & dataframe_last_exon$Func.refGene == 'exonic'), ])
  merged_table_exonic_entries_with_cadd_le <- nrow(dataframe_last_exon[ (!(is.na(dataframe_last_exon$CADD_phred)) & dataframe_last_exon$Func.refGene == 'exonic'), ])
  merged_table_exonic_entries_with_both_le <- nrow(dataframe_last_exon[ (!(is.na(dataframe_last_exon$CADD_phred)) & !(is.na(dataframe_last_exon$LoFtool)) & dataframe_last_exon$Func.refGene == 'exonic'), ])
  merged_table_splicing_entries_with_spliceai_le <- nrow(dataframe_last_exon[ (!(is.na(dataframe_last_exon$spliceai_max)) & dataframe_last_exon$Func.refGene == 'splicing'), ])
  
  
  total_and_missing_entries <- data.frame(last_exon_counts = c(merged_table_exonic_entries_with_cadd_le, 
                                                               merged_table_exonic_entries_with_loftool_le, 
                                                               merged_table_exonic_entries_with_both_le, 
                                                               merged_table_splicing_entries_with_spliceai_le),
                                          row.names = c('Exonic with CADD score',
                                                        'Exonic with LoFtool score',
                                                        'Exonic with CADD and LoFtool scores',
                                                        'Splicing with SpliceAI score'))
  # This creates row.names as a column with no table to make it bold with ` ` = formatter...
  total_and_missing_entries$blank <- rownames(total_and_missing_entries)
  rownames(total_and_missing_entries) <- c()
  total_and_missing_entries <- total_and_missing_entries[, c(2, 1)]
  colnames(total_and_missing_entries) <- c(' ', 'Num. (%) of variants from terminal exons')
  
  percentage_with_info <- formattable(total_and_missing_entries, align = c('l', rep('c', ncol(total_and_missing_entries) -1)),
                                      list(` ` = formatter("span", style = ~formattable:::style(display = "block", font.weight = "bold")),
                                           `Num. (%) of variants from terminal exons` = formatter("span", x ~ paste( x, ' (', percent(x / nrow(dataframe_last_exon)), ')') )))
  print(percentage_with_info)
  
  ###############################
  # 2. Basic Statistics Summary #
  ###############################
  
  filter_compar <- data.frame(total_counts = c(total_entries_merged, nrow(dataframe_last_exon), nrow(dataframe_not_keep)),
                              row.names = c('Entries Merged VEP Annot.', 'Entries in Terminal Exons', 'Removed Variants') )
  #filter_compar[is.na(filter_compar)] = '' # Remove NAs from appearing on table if not percentage

  filter_compar$blank <- rownames(filter_compar)
  rownames(filter_compar) <- c()
  filter_compar <- filter_compar[, c(2, 1)]
  colnames(filter_compar) <- c(' ', 'Num. (%) of variants regarding merged entries')
  filter_summ <- formattable(filter_compar, align = c('l', rep('c', ncol(total_and_missing_entries) -1)),
                               list(` ` = formatter("span", style = ~formattable:::style(display = "block", font.weight = "bold")),
                                    `Num. (%) of variants regarding merged entries` = formatter("span", x ~ paste(x, ' (', percent( x / total_entries_merged), ')' ) )))
  print(filter_summ)
}

predictors_summary_tables <- function(plof_dataframe, delmis_dataframe, saveTABLE = FALSE){
  
  #--|--|--|--|--|--|--|--|--|--|--|--|--#            
  # pLoF class - Merged with VEP browser #
  #--|--|--|--|--|--|--|--|--|--|--|--|--#
  #
  # BASIC
  total_entries_plof <- nrow(plof_dataframe) 
  total_last_exon_plof <- nrow(plof_dataframe %>% filter(in_last_exon, TRUE))
  total_removed_plof <- nrow(plof_dataframe %>% filter(whole_criteria_keep == FALSE))
  exonic_entries_plof <- plof_dataframe %>% filter(Func.refGene == 'exonic')
  splicing_entries_plof <- plof_dataframe %>% filter(Func.refGene == 'splicing')
  #
  # CASES
  plof_cadd <- nrow( exonic_entries_plof %>% filter( !(is.na(CADD_phred)) ) )
  plof_cadd_na <- nrow( exonic_entries_plof %>% filter( is.na(CADD_phred)) )
  plof_loftool <- nrow( exonic_entries_plof %>% filter( !(is.na(LoFtool)) ) )
  plof_loftool_na <- nrow( exonic_entries_plof %>% filter( is.na(LoFtool)) )
  plof_bothCADDlof <- nrow( exonic_entries_plof %>% filter( !(is.na(CADD_phred)) & !(is.na(LoFtool)) ) )
  plof_bothCADDlof_na <- nrow( exonic_entries_plof %>% filter( is.na(CADD_phred) & is.na(LoFtool)) )
  plof_spliceai <- nrow( splicing_entries_plof %>% filter( !(is.na(spliceai_max)) ) )
  plof_spliceai_na <- nrow( splicing_entries_plof %>% filter( is.na(spliceai_max) ) )
  
  #-|--|--|--|--|--|--|--|--|--|--|--|--|--#            
  # DelMis class - Merged with VEP browser #
  #-|--|--|--|--|--|--|--|--|--|--|--|--|--#
  #
  # BASIC
  total_entries_delmis <- nrow(delmis_dataframe) 
  total_last_exon_delmis <- nrow(delmis_dataframe %>% filter(in_last_exon, TRUE))
  total_removed_delmis <- nrow(delmis_dataframe %>% filter(whole_criteria_keep == FALSE))
  exonic_entries_delmis <- delmis_dataframe %>% filter(Func.refGene == 'exonic')
  splicing_entries_delmis <- delmis_dataframe %>% filter(Func.refGene == 'splicing')
  #
  # CASES
  delmis_cadd <- nrow( exonic_entries_delmis %>% filter( !(is.na(CADD_phred)) ) )
  delmis_cadd_na <- nrow( exonic_entries_delmis %>% filter( is.na(CADD_phred)) )
  delmis_loftool <- nrow( exonic_entries_delmis %>% filter( !(is.na(LoFtool)) ) )
  delmis_loftool_na <- nrow( exonic_entries_delmis %>% filter( is.na(LoFtool)) )
  delmis_bothCADDlof <- nrow( exonic_entries_delmis %>% filter( !(is.na(CADD_phred)) & !(is.na(LoFtool)) ) )
  delmis_bothCADDlof_na <- nrow( exonic_entries_delmis %>% filter( is.na(CADD_phred) & is.na(LoFtool)) )
  delmis_spliceai <- nrow( splicing_entries_delmis %>% filter( !(is.na(spliceai_max)) ) )
  delmis_spliceai_na <- nrow( splicing_entries_delmis %>% filter( is.na(spliceai_max) ) )
  
  ###################################################
  # 1. Percentages of entries with predictors score #
  ###################################################
  
  total_and_missing_entries <- data.frame(plof_scores = c(plof_cadd, plof_loftool, plof_bothCADDlof, plof_spliceai, total_entries_plof),
                                          plof_nas = c(plof_cadd_na, plof_loftool_na, plof_bothCADDlof_na, plof_spliceai_na, total_entries_plof),
                                          missense_scores = c(delmis_cadd, delmis_loftool, delmis_bothCADDlof, delmis_spliceai, total_entries_delmis),
                                          missense_nas = c(delmis_cadd_na, delmis_loftool_na, delmis_bothCADDlof_na, delmis_spliceai_na, total_entries_delmis),
                                          row.names = c('CADD annovar (Exonic Variants)',
                                                        'LoFtool (Exonic Variants)',
                                                        'CADD and LoFtool (Exonic Variants)',
                                                        'SpliceAI (Splicing Variants)',
                                                        'Total Number of Variants (in merged file)'))
  # This creates row.names as a column with no table to make it bold with ` ` = formatter...
  total_and_missing_entries$blank <- rownames(total_and_missing_entries)
  rownames(total_and_missing_entries) <- c()
  total_and_missing_entries <- total_and_missing_entries[, c(5, 1:4)]
  colnames(total_and_missing_entries) <- c(' ', 'Num. (%) variants with scores LOF',
                                           'Num. (%) variants with NAs LOF',
                                           'Num. (%) variants with scores Missense',
                                           'Num. (%) variants with NAs Missense')
  
  percentage_with_info <- formattable(total_and_missing_entries, align = c('l', rep('c', ncol(total_and_missing_entries) -1)),
                                      list(` ` = formatter("span", style = ~ formattable:::style(display = "block", font.weight = "bold")),
                                           `Num. (%) variants with scores LOF` = formatter("span", x ~ paste( x, ' (', percent(x / total_entries_plof), ')')),
                                           `Num. (%) variants with NAs LOF` = formatter("span", x ~ paste( x, ' (', percent(x / total_entries_plof), ')' )),
                                           `Num. (%) variants with scores Missense` = formatter("span", x ~ paste( x, ' (', percent(x / total_entries_delmis), ')')),
                                           `Num. (%) variants with NAs Missense`  = formatter("span", x ~ paste(x, ' (', percent(x / total_entries_delmis), ')'))
                                           ))
  print(percentage_with_info)
  
  
  ###############################
  # 2. Basic Statistics Summary #
  ###############################
  
  last_exon_removed_counts <- data.frame(total_counts = c(total_last_exon_plof, total_removed_plof, total_entries_plof),
                              row.names = c('# Variants in terminal exon', '# of Removed Variants', 'Total Number of Variants in merged file') )
  #filter_compar[is.na(filter_compar)] = '' # Remove NAs from appearing on table if not percentage
  
  last_exon_removed_counts$blank <- rownames(last_exon_removed_counts)
  rownames(last_exon_removed_counts) <- c()
  last_exon_removed_counts <- last_exon_removed_counts[, c(2, 1)]
  colnames(last_exon_removed_counts) <- c(' ', 'Num. variants (% respect to total merged file) LOF')
  filter_summ <- formattable(last_exon_removed_counts, align = c('l', rep('c', ncol(last_exon_removed_counts) -1)),
                             list(` ` = formatter("span", style = ~formattable:::style(display = "block", font.weight = "bold")),
                                  `Num. variants (% respect to total merged file) LOF` = formatter("span", x ~ paste(x, ' (', percent( x / total_entries_plof), ')' ) )
                                  ))
  print(filter_summ)
}

check_loftool_pass <- function(dataframe){
  cat(paste0('Number of variants for ', deparse(substitute(dataframe)) ,' that pass LoFtool: ', nrow(dataframe %>% filter(loftool_damag == TRUE)), '\n'))
}


#########################################
# {|-Data Imports and Basic Handling-|} #
#########################################

# Our table, to merge VEP's annotation while keeping the format that we desire
plof_table_to_match <- read_tsv(file = 'filtered_tables_postNAissue_march/germline.TCGAqc_filt.GQ.DP.NSY.VAF.AD.BCD.PASS.TCGA_AF.FUNC.plof.rare_0.001.header.tsv')
delmiss_table_to_match <- read_tsv(file = 'filtered_tables_postNAissue_feb//delmis/germline.TCGAqc_filt.GQ.DP.NSY.VAF.AD.BCD.PASS.TCGA_AF.FUNC.delmis.rare_0.001.header.tsv')


# We take the VEPbrowser-annotated file keeping the
# relevant columns and removing possible duplicates 
plof_vep_annot_relevantcols_nodup <- vep_output_relevantcols('filtered_tables_postNAissue_march/germline.TCGAqc_filt.GQ.DP.NSY.VAF.AD.BCD.PASS.TCGA_AF.FUNC.sorted.plof.rare_0.001.header.vep_annotated.txt')
delmiss_vep_annot_relevantcols_nodup <- vep_output_relevantcols('filtered_tables_postNAissue_feb//delmis/germline.TCGAqc_filt.GQ.DP.NSY.VAF.AD.BCD.PASS.TCGA_AF.FUNC.sorted.delmis.rare_0.001.header.vep_loftool_spliceai_masked_annotated.uniq.vcf')

# To focus on unique variants (some exist due to annotation or samples).
# I filter the same way Luis does and thus, we will have more representative
# proportions in our various distributions plots but I also keep one where
# there may be duplicates at sample level to perform mean gene mutation
# frequency where we need total population info
#
#table_to_match_unique_vars_old <- table_to_match[!(duplicated(table_to_match[, c('Otherinfo3', 'Gene')])), ]
#table_to_match_unique_vars_freqs_old <- table_to_match[!(duplicated(table_to_match[, c('Otherinfo3', 'Gene', 'SAMPLE')])), ]

# Here I define Unique variants as VarID + Gene + Transcript to obtain a better
# merge with our input table, as even though the resulting table is smaller than
# before (with VarID + Gene) I'm more confident in the merge that it is performed
# as now the rows really do correspond to the entries

not_snamake <- FALSE
if(not_snakemake){
  plof_table_to_match_unique_vars <- plof_table_to_match[!(duplicated(plof_table_to_match[, c('Otherinfo3', 'Gene', 'Feature')])), ]
} else{
  plof_table_to_match_unique_vars <- plof_table_to_match
}
#table_to_match_unique_vars_freqs <- table_to_match[!(duplicated(table_to_match[, c('Otherinfo3', 'Gene', 'Feature', 'SAMPLE')])), ]
delmiss_table_to_match_unique_vars <- delmiss_table_to_match[!(duplicated(delmiss_table_to_match[, c('Otherinfo3', 'Gene', 'Feature')])), ]

# All genes present in Pelayo's annotated .vcf
all_genes <- scan("genes_sorted_uniq.txt", what = "character")
all_genes_without_semmic <- unique(unlist(lapply(all_genes, function(x) strsplit(x, split = ";")[[1]]))) # Remove misleading entries

# All patients' barcodes, needed for mutation frequency 
all_patients <- sapply(scan("./samples_ID.txt", what = "character"), function(x) substr(x, 1, 12))

# CPGs to compare approaches, from our Gitlab; keep those present in our input file
CPGs_list <- scan(file = 'CPG.csv', what = 'character')
CPGs_list_available <- CPGs_list[CPGs_list %in% all_genes]

# Cancer Consensus Genes from COSMIC https://cancer.sanger.ac.uk/census (GRCh37-COSMIC v95)
ccg_list_tab <- read_tsv('Census_allMon_Feb_14_11_23_46_2022.tsv')
ccg_list <- unique(filter(ccg_list_tab, is.na(Germline), Somatic == "yes")$`Gene Symbol`)
ccg_list_available <- ccg_list[ccg_list %in% all_genes]
ccg_list_available_not_CPG <- ccg_list_available[!(ccg_list_available %in% CPGs_list_available)]

# Extract cancer related genes and obtain a random 
# subset of 1k genes which aren't among those 
all_cancer_related_genes <- unique(c(CPGs_list, ccg_list))
all_genes_without_semmic_and_cancer <- all_genes_without_semmic[!(all_genes_without_semmic %in% all_cancer_related_genes)]
set.seed(1)
random_genes_1K <- sample(all_genes_without_semmic_and_cancer, size = 1000, replace = FALSE)


#######################
# {|-Work Per Se-|} #
#######################

# Subset to only keep the entries that match either of our vectors of interest
# for the old mutation freqs. table and the current unique variant approach
plof_table_to_match_unique_vars_subs <- plof_table_to_match_unique_vars[plof_table_to_match_unique_vars$consensus_gene_symbol %in% c(CPGs_list_available, ccg_list_available_not_CPG, random_genes_1K),]
#table_to_match_unique_vars_freqs_subs <- table_to_match_unique_vars_freqs[table_to_match_unique_vars_freqs$consensus_gene_symbol %in% c(CPGs_list_available, ccg_list_available_not_CPG, random_genes_1K),]
delmiss_table_to_match_unique_vars_subs <- delmiss_table_to_match_unique_vars[delmiss_table_to_match_unique_vars$consensus_gene_symbol %in% c(CPGs_list_available, ccg_list_available_not_CPG, random_genes_1K),]

# Add a column with factor for each of the groupings
plof_table_to_match_unique_vars_subs <- gene_grouping(plof_table_to_match_unique_vars_subs)
#table_to_match_unique_vars_freqs_subs <- gene_grouping(table_to_match_unique_vars_freqs_subs)
delmiss_table_to_match_unique_vars_subs <- gene_grouping(delmiss_table_to_match_unique_vars_subs)

old_merge_id <- FALSE # To not add extra unneeded columns
if(old_merge_id){
  # Now we add the uniqueID based on VariantID-Gene combination to both to-merge tables
  plof_vep_annot_relevantcols_nodup$pasteID_geneonly <- paste(plof_vep_annot_relevantcols_nodup$Uploaded_variation, plof_vep_annot_relevantcols_nodup$Gene, sep = '_')
  table_to_match_unique_vars_subs$pasteID_geneonly <- paste(table_to_match_unique_vars_subs$Otherinfo3, table_to_match_unique_vars_subs$Gene, sep = "_")
  table_to_match_unique_vars_freqs_subs$pasteID_geneonly <- paste(table_to_match_unique_vars_freqs_subs$Otherinfo3, table_to_match_unique_vars_freqs_subs$Gene, sep = "_")
}

## CURRENT APPROACH
# Different ID for non-unique variants (repeated if they are on different transcripts, better criteria for the next merging process) 
plof_vep_annot_relevantcols_nodup$pasteID_genetrans <- paste(plof_vep_annot_relevantcols_nodup$Uploaded_variation, plof_vep_annot_relevantcols_nodup$Gene, plof_vep_annot_relevantcols_nodup$Feature, sep = "_")
plof_table_to_match_unique_vars_subs$pasteID_genetrans <- paste(plof_table_to_match_unique_vars_subs$Otherinfo3, plof_table_to_match_unique_vars_subs$Gene, plof_table_to_match_unique_vars_subs$Feature, sep = "_")

delmiss_vep_annot_relevantcols_nodup$pasteID_genetrans <- paste(delmiss_vep_annot_relevantcols_nodup$Uploaded_variation, delmiss_vep_annot_relevantcols_nodup$Gene, delmiss_vep_annot_relevantcols_nodup$Feature, sep = "_")
delmiss_table_to_match_unique_vars_subs$pasteID_genetrans <- paste(delmiss_table_to_match_unique_vars_subs$Otherinfo3, delmiss_table_to_match_unique_vars_subs$Gene, delmiss_table_to_match_unique_vars_subs$Feature, sep = "_")

# Combine both tables and change some column to filter them correctly, keep all that do not 
#matched_annotation_uniqvars_old <- merge(table_to_match_unique_vars_subs, plof_vep_annot_relevantcols_nodup, by = 'pasteID_geneonly', all.x = TRUE) ## Old approach
plof_matched_annotation_uniqvars <- merge(plof_table_to_match_unique_vars_subs, plof_vep_annot_relevantcols_nodup, by = 'pasteID_genetrans', all.x = TRUE)
#matched_annotation_for_mutfreqs <- merge(table_to_match_unique_vars_freqs_subs, plof_vep_annot_relevantcols_nodup, by = 'pasteID_genetrans', all.x = TRUE)
delmiss_matched_annotation_uniqvars <- merge(delmiss_table_to_match_unique_vars_subs, delmiss_vep_annot_relevantcols_nodup, by = 'pasteID_genetrans', allx.true = TRUE)

# Change all predictors score of interest to numeric type,
# the NAs that introduces are relevant to us
plof_matched_annotation_uniqvars <- important_cols_as_num(plof_matched_annotation_uniqvars)
#matched_annotation_for_mutfreqs <- important_cols_as_num(matched_annotation_for_mutfreqs)
delmiss_matched_annotation_uniqvars <- important_cols_as_num(delmiss_matched_annotation_uniqvars)

####################################################################
# Obtain Different Subsets for Filtering Criteria Summary Analysis #
####################################################################
# In this part of the script we will focus on obtaining the tables 
# that summarize the basic statistics of how many variants do we keep
# with each combination of filters

# We now obtain some columns with detailed information and filter in various ways
plof_matched_annotation_uniqvars_moreinfo <- dataframe_vep_filterer(plof_matched_annotation_uniqvars,
                                                                 cpgs = CPGs_list, spliceai_thr = 0.8,
                                                                 cadd_thr = 15, loftool_thr = 0.334)
plof_matched_annotation_uniqvars_moreinfo$whole_criteria_keep_unlist <- unlist(plof_matched_annotation_uniqvars_moreinfo$whole_criteria_keep)
plof_matched_annotation_uniqvars_moreinfo <- plof_matched_annotation_uniqvars_moreinfo %>% filter(Func.refGene != 'exonic;splicing')

delmiss_matched_annotation_uniqvars_moreinfo <- dataframe_vep_filterer(delmiss_matched_annotation_uniqvars,
                                                                    cpgs = CPGs_list, spliceai_thr = 0.8,
                                                                    cadd_thr = 15, loftool_thr = 0.334)
delmiss_matched_annotation_uniqvars_moreinfo <- delmiss_matched_annotation_uniqvars_moreinfo %>% filter(Func.refGene != 'exonic;splicing')

# To obtain the summary, need to do it for CPGs and non-CPGs
# CPGs #
CPG_plof_matched_annotation_uniqvars_moreinfo <- plof_matched_annotation_uniqvars_moreinfo %>% filter(cpg_or_not == TRUE)
CPG_delmiss_matched_annotation_uniqvars_moreinfo <- delmiss_matched_annotation_uniqvars_moreinfo %>% filter(cpg_or_not == TRUE)

predictors_summary_tables(plof_dataframe = CPG_plof_matched_annotation_uniqvars_moreinfo, 
                          delmis_dataframe = CPG_delmiss_matched_annotation_uniqvars_moreinfo)

# non-CPG genes #
nonCPG_plof_matched_annotation_uniqvars_moreinfo <- plof_matched_annotation_uniqvars_moreinfo %>% filter(cpg_or_not == FALSE)
nonCPG_delmiss_matched_annotation_uniqvars_moreinfo <- delmiss_matched_annotation_uniqvars_moreinfo %>% filter(cpg_or_not == FALSE)

predictors_summary_tables(plof_dataframe = nonCPG_plof_matched_annotation_uniqvars_moreinfo,
                          delmis_dataframe = nonCPG_delmiss_matched_annotation_uniqvars_moreinfo)


## To extract the list of variants that are removed by this filters
mergeIDs_to_remove <- (plof_matched_annotation_uniqvars_moreinfo %>%
    filter(whole_criteria_keep == FALSE))$'pasteID_genetrans'


###### Total dataset
plof_march_table_to_match <- plof_table_to_match_unique_vars
plof_march_table_to_match$pasteID_genetrans <- paste(plof_march_table_to_match$Otherinfo3, plof_march_table_to_match$Gene, plof_march_table_to_match$Feature, sep = "_")
plof_march_matched_annotation_uniqvars <- merge(plof_march_table_to_match, plof_vep_annot_relevantcols_nodup, by = 'pasteID_genetrans', all.x = TRUE)
plof_march_matched_annotation_uniqvars <- important_cols_as_num(plof_march_matched_annotation_uniqvars)

plof_march_matched_annotation_uniqvars <- last_exon_checker(plof_march_matched_annotation_uniqvars)
plof_march_matched_annotation_uniqvars <- spliceai_damag_sapply(plof_march_matched_annotation_uniqvars)
plof_march_matched_annotation_uniqvars <- cadd_damag_sapply(plof_march_matched_annotation_uniqvars)
plof_march_matched_annotation_uniqvars <- loftool_damag_sapply(plof_march_matched_annotation_uniqvars)
plof_march_matched_annotation_uniqvars <- cadd_loftool_damag_sapply(plof_march_matched_annotation_uniqvars)
plof_march_matched_annotation_uniqvars <- wholecriteria_damag_sapply(plof_march_matched_annotation_uniqvars)
plof_march_matched_annotation_uniqvars$whole_criteria_keep_unlist <- unlist(plof_march_matched_annotation_uniqvars$whole_criteria_keep)


# Other tests #
quick_manual_check <- FALSE
if (quick_manual_check) {
  ## Manual check of those that will be removed to see if there
  ## has been any mistake when applying the filtering criteria
  plof_matched_annotation_uniqvars_moreinfo[, c('in_last_exon', 'Func.refGene', 'gene_group', 'CADD_phred', 'LoFtool', 'spliceai_max', 'whole_criteria_keep')] %>%
    filter(in_last_exon == TRUE & whole_criteria_keep == FALSE)

  # To see the differences between those who have NAs for the Exon information (x/out_of_total) and for nonNAs
  # Even though there is some noise, makes sense that most of the NAs for that field are coming from variants
  # that are related to regulatory regions, downstream and upstream, splicing sites, etc...
  table(filter(plof_matched_annotation_uniqvars_moreinfo, is.na(in_last_exon) )$Consequence) %>% as.data.frame() %>% arrange(desc(Freq))
  table(filter(plof_matched_annotation_uniqvars_moreinfo, !is.na(in_last_exon) )$Consequence) %>% as.data.frame() %>% arrange(desc(Freq))
}

test_possible_loftool_threshold <- FALSE
if (test_possible_loftool_threshold){
  # This was done to check a suitable theshold for LoFtool
  salmorejo_cpg_clinvar <- matched_annotation_moreinfo %>% filter( (cpg_or_not == TRUE)  & (clinvar_patho == TRUE) ) 
  salmorejo_cpg_nonclinvar <- matched_annotation_moreinfo %>% filter( (cpg_or_not == TRUE)  & (clinvar_patho == FALSE) ) 
  salmorejo_noncpg_clinvar <- matched_annotation_moreinfo %>% filter( (cpg_or_not == FALSE)  & (clinvar_patho == TRUE) ) 
  salmorejo_noncpg_nonclinvar <- matched_annotation_moreinfo %>% filter( (cpg_or_not == FALSE)  & (clinvar_patho == FALSE) ) 
  filters_summary_tables_old(salmorejo_cpg_clinvar)
}

# Correspondence of Annovar CADD Scores (the ones in the input file)
# and the VEP annotated ones (from the browser version)
plot(matched_annotation_uniqvars$CADD_phred, matched_annotation_uniqvars$CADD_PHRED)