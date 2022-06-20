#########################
# Perpetrated by Manuel #  Updated Version 18/01/2022 
#########################

# Loads all R packages that we will need
x <- c("tidyverse","ggplot2", "RColorBrewer", "plinkQC", "patchwork")
lapply(x, require, character.only = TRUE)

# YOU MUST HAVE plink INSTALLED IN YOUR WORKSTATION!
# Stores plink's path for 'plinkQC' steps
path_plink <- gsub("plink: ", "", system('whereis plink', intern = TRUE))

# Colour-blindness palette if I want to use it
color_blind <- c("#009E73","#0072B2", "#D55E00", "#CC79A7")


pcaineitor_pca <- function(start_dir, graph_output_dir = 'png_PCA_graphs', vcf.loc,
                           tsv.loc, pop.loc, center.loc, high_ld_region,
                           opt_prf = '', further_prun = FALSE, ld_par = '50 5 0.2',
                           hwe = '5e-324', pca.dim = 8000, pve_dim = 25,
                           pca_dim_to_plot = 4, 
                           filter_outliers = FALSE, outliers_PC = 'PC3',
                           outliers_PC_thr = -0.07, outliers_PC_method = 'le',
                           save_img = TRUE, save_ann_pca = TRUE,
                           plot_title = 'Common Variants across All Population',
                           plot_subtitle = 'No IQR Filter applied - Outliers Not Removed'){
  
#############################
# Description and arguments #
#############################
  #
  # Takes a .vcf file and uses various plink commands to perform some QC steps
  # before performing PCA and, lastly, plots some PCs combinations 
  #
  # {Args}:
  #   start_dir: Directory where you want to store all the jobs from plink,
  #              and where all the other inputs files' path will be referenced
  #   graph_out_dir: Name of the directory to store all plots, creates it
  #                  if it doesn't exist
  #   vcf.loc : The location path of the input .vcf referenced as if you
  #             you were currently in 'start_dir'
  #   tsv.loc : The location path of the .tsv returned by Luis' pipeline. This 
  #             file is similar to the .vcf although it does not have the nine
  #             canonical .vcf headers and the genotypes only reference to the
  #             variants count per sample (0 - 0/0, 1 - 0/1, 2 - 1/1); needed 
  #             for some number of variants per sample plots
  #   pop.loc : The location path of the population cluster - barcode .tsv, used
  #             to paint PCAs with population and subpopulation info.
  #             This was obtained from another script of mine that checks the main
  #             population metadata files that we used to identify each sample's
  #             ancestry; that script is named 'obtain_missing_PAM_popclusters.R'
  #   center.loc : The location path of the sequencing center that's the origin
  #                of the samples.
  #                This was obtained from the last two digits of each barcode
  #                with R code by hand; anytime we obtain a new input we should
  #                update this file by mapping those digits to the codes from:
  #                https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/center-codes  
  #   high_ld_region : Path to high LD regions file obtained from here
  #                    https://github.com/meyer-lab-cshl/plinkQC/blob/master/inst/extdata/high-LD-regions-hg19-GRCh37.bed 
  #   opt_prf : An optional prefix for the output files, just to return a more
  #             detailed name to distinguish results for different tests
  #   further_prun : TRUE if you want to perform LD pruning within your
  #                   sample alongside removing those located in reference 
  #                   high LD regions (for hg19), FALSE if just the latter
  #   ld_par : Parameters for plink's --indep-pairwise for LD pruning as in
  #            <windows size> <step size> <r^2 cutoff>             
  #            https://www.cog-genomics.org/plink/1.9/ld#indep
  #   hwe : Filters out variants which have a Hardy-Weinberg equilibrium
  #         exact test p-value below this threshold. The default was chose
  #         by counting our genotypes and inputing them in this calculator:
  #         https://www.cog-genomics.org/software/stats 
  #   pca.dim : By defaults plink extracts the top 20 dimensions and it has
  #             an upper limit of 8000, hence the default value of this argument.
  #             It will return some warnings if your input has less that 8k 
  #             dimensions/samples, so you can tune it up to your experiment
  #             if you don't want to see them but thse warnings don't have any 
  #             negative impact.
  #             To obtain the percentage of explained variance its better to 
  #             extract all possible dimensions, not plink's 20 default.
  #   pve_dim : Number of PCs to show in percent variance explained barplot
  #   pca_dim_to_plot : Number of components/dimensions to plot, it will 
  #                     generate plots for each possible combination (I recommend
  #                     4 as it is manageable, 5 if you want to explore more)
  #   filter_outliers : During some testing, we saw some outliers at some PC.
  #                     First run the script with this argument as FALSE, identify
  #                     in which PC (for now only able to work on one) and at what
  #                     value/threshold (above or below) you do see outliers so 
  #                     they can be removed in the next runs of the analysis by
  #                     creating a tab delimited file that plink asks for.
  #   outliers_PC : Selects PC where the outliers where identified from first
  #                 script run; currently, this is implemented to work on just
  #                 single PC at a time.
  #   outliers_PC_thr : Indicate a positive or negative value seen in the PC 
  #                     combinations where outliers where found; samples below
  #                     this threshold at the chosen PC will be removed.
  #   outliers_PC_method : 'le' for less or equal than threshold, 'ge' for
  #                         greater or equal than
  #   save_img : If you want to store all plots at 'graph_out_dir'
  #   save_ann_pca : To store a .tsv of the PCA results for PCs 1-4
  #   plot_title : Title of the main PCA plots
  #   plot_subtitle : Subtitle of the main PCA plots (to include QC steps info,
  #                   or whatever specification that you see fit)
  #
  # {Returns} :
  #   plink returns a lot of files, but probably the more useful outputs from
  #   this script will be
  #     {plink.fileprefix}.5_PCA.eigenval : Eigenvalues of the different PCs,        
  #       remember that at most it will return 8000 as it is plink's limit 
  #     {plink.fileprefix}.5_PCA.eigenvec : The values of the PCs themselves
  #     {plink.fileprefix}._cluster_PCA1-4.tsv : Subset of the whole PCA with  
  #                                              PCs 1 through 4   
  #   PCA_graphs/ : Directory for all the plots created, including the 
  #                 percentage variance explained per PC and the various
  #                 combinations of PCs and metadata files
  #
  
  # Changes directory and creates folder for the plots
  setwd(start_dir)
  if (!dir.exists(graph_output_dir)) {
    dir.create(graph_output_dir)
  }
  
  # Creates the prefix for all files taking into account your parent folder
  # (should be the population cluster of analysis) and its parent folder
  # opt_prf if you want to be more specific of your analysis
  plink.fileprefix <- paste0(opt_prf, '_', basename(dirname(getwd())), '_', basename(getwd()))
  
  # Checks if the most important input file path is right
  if (!file.exists(vcf.loc)) {
    stop('Error: Could not find input .vcf!\n')
  }

  ####################
  # QC Filter Steps  #
  ####################
  
  ##|{ .vfc to .bed }|## 
  #
  # This writes the .vcf into plink's binary format as it is more efficient to work with
  system(command = sprintf('plink --vcf %s --make-bed --out %s.0_init',
                           vcf.loc, plink.fileprefix))
  
  ##|{ Outliers Handling }|##
  #
  # Once you have ran this whole function the first time you can choose which PC
  # has outliers to obtain their IDs for future extraction of those for better PCAs
  if (filter_outliers){
    pca_14 <- read_tsv(file = paste0(plink.fileprefix, '.5_PCA.PCs1-4.outliers_not_removed.tsv'))
    if (outliers_PC_method == 'le') { # Lower or equal than threshold are removed
      outids <- pull(select(filter(pca_14, get(outliers_PC) <= outliers_PC_thr), 'TCGA.PatientID'))
    }else if (outliers_PC_method == 'ge') { # Greater or equal than threshold are removed
      outids <- pull(select(filter(pca_14, get(outliers_PC) >= outliers_PC_thr), 'TCGA.PatientID'))
    }
    write.table(data.frame(outids, outids), file = paste0(plink.fileprefix, '.outliers.tsv'), sep = " ",
                row.names = FALSE, col.names = FALSE, quote = FALSE) # Saves outliers IDs
    
    # Remove outliers (samples) identified prior to any performing any QC steps on the data
    system(command = sprintf('plink --bfile %s.0_init --remove %s.outliers.tsv --make-bed --out %s.1_init',
                             plink.fileprefix, plink.fileprefix, plink.fileprefix))
  } else {
    warning("Could not find Outliers Barcode file.
            \nNo sample will be removed beforehand, but you should check some PCs
            \nfor possible outliers!\n")
    # This is the laziest way to rename the file I could think of
    system(command = sprintf('plink --bfile %s.0_init --make-bed --out %s.1_init',
                             plink.fileprefix, plink.fileprefix))
  }

  ##|{ Heterozygosity Rate Calling }|##
  #
  # By default and in some QC references they use a heterozygosity rate +-3
  # standard deviations away from the mean ratio it computes
  fail_het_imiss <- check_het_and_miss(indir = '.', name = paste0(plink.fileprefix, '.1_init'),
                                       hetTh = 3, interactive = FALSE,
                                       label_fail = FALSE, path2plink = path_plink,
                                       showPlinkOutput = FALSE)
  
  # Write the IDs that failed the test to extract them from our plink binary
  # file (.bed) in further steps
  write.table(data.frame(fail_het_imiss$fail_het$FID, fail_het_imiss$fail_het$IID),
              file = paste0(plink.fileprefix, '.2_heterotest.fail.IDs'), sep = " ",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Naming something as 'posthetero' may be my single greatest achievement
  system(command = sprintf('plink --bfile %s.1_init --remove %s.2_heterotest.fail.IDs --make-bed --out %s.2_posthetero',
                           plink.fileprefix, plink.fileprefix, plink.fileprefix))
  
  ##|{ LD Computation within our variants + High LD Region }|##
  #
  # Just an optional QC step to compare the results, as I haven't found a way
  # to estimate a good r² cutoff and, in our current dataset, it removes
  # almost half of our variants but returns much better intrapopulation clustering!
  system(command = sprintf('plink -bfile %s.2_posthetero --exclude %s --range --indep-pairwise %s --make-bed --out %s.3_postld',
                           plink.fileprefix, high_ld_region, ld_par, plink.fileprefix))
  
  ##|{ Autosomes and Hardy-Weinberg Equilibrium  }|##
  #
  # First, removes those variants that where pruned out by the Linkage 
  # Disequilibrium test, then excludes all non-autosomal
  # variants to avoid some variant genotyping issues due to different
  # copies in sexual chromosomes.
  # Lastly, it performs a Hardy-Weinberg equilibrium test to filter out
  # some more variants below the test p-value threshold 'hwe'
  system(command = sprintf('plink --bfile %s.3_postld  --extract %s.3_postld.prune.in --autosome --hwe %s --make-bed --out %s.4_hweq',
                           plink.fileprefix, plink.fileprefix, hwe, plink.fileprefix))
  
  ##################
  # PCA with plink #
  ##################
  
  # Performs the PCA with final QC-filtered file
  # It also asks for the covariance matrix
  system(command = sprintf('plink --bfile %s.4_hweq --make-rel square0 --pca %s --out %s.5_PCA',
                           plink.fileprefix, pca.dim, plink.fileprefix))
  
  # plink has a limitation of 8000 components/samples to extract,
  # if we have more than that limit we overcome this issue by 
  # obtaining the diagonal of the covariance matrix, as it is 
  # always returned with all the eigenvalues regardless of
  # plink's limitations 
  #
  # Thus, this awk command goes printing all elements of the diagonal
  system(command = sprintf("awk '{print $NR}' %s.5_PCA.rel > %s.5_PCA.diagonal.covariance", plink.fileprefix, plink.fileprefix))
  
  
  ################
  # PCA Plotting #
  ################
  
  # We obtain the eigenvectors of our PCA
  pca <- read_table(dir(pattern = '*.5_PCA.eigenvec'), col_names = FALSE)
  pca <- pca[, -2] # Plink duplicates ID column, not useful to us 
  names(pca)[1] <- "TCGA.PatientID"
  names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
  
  ##|{ Metadata for PCAs }|##
  #
  # Number of variants per sample
  if (!file.exists(tsv.loc)) { 
    stop("Error: Could not find input .tsv or it is not specified!\n")
  } else {
    var_count_tsv <- read_tsv(tsv.loc)
    nvariants <- data.frame('TCGA.PatientID' = colnames(var_count_tsv),
                            'Num.Variants' = colSums(var_count_tsv))
  }
  # Population cluster info
  if (!file.exists(pop.loc)) {
    stop('Error: Could not find population metadata .tsv')
  } else {
    pop_info <- read_tsv(file = pop.loc, col_select = c(1, 4, 5))
    colnames(pop_info) <- c('TCGA.PatientID', 'PCA.PAM.Ancestry.SubCluster', 'PCA.PAM.Ancestry.Cluster') 
  }
  # Sequencing center info
  if (!file.exists(center.loc)) {
    stop('Error: Could not find seq. center metadata .tsv!')
  } else {
    center_info <- read_tsv(center.loc)
  }
  # Now we include all metadata in our pca matrix
  pca <- merge(pca, nvariants, by = 'TCGA.PatientID') 
  pca <- merge(pca, pop_info, by = 'TCGA.PatientID') 
  pca <- merge(pca, center_info, by = 'TCGA.PatientID')
  
  # Stores four main PCs and metadata to share it to the group if they want to
  # produce their own plots
  # ALSO, this is needed in order to filter out outliers present in one of the
  # first components, as after a first run of the script you can choose PC and
  # threshold to remove those IDs
  if (save_ann_pca) {
    if (filter_outliers){
      sub_pca <-pca[, c(1:5, ncol(pca)-2, ncol(pca)-1, ncol(pca))]
      write.table(sub_pca, file = paste0(plink.fileprefix, '.5_PCA.PCs1-4.tsv'),
                  sep = "\t", quote = FALSE, row.names = FALSE)
    }else{
      sub_pca <-pca[, c(1:5, ncol(pca)-2, ncol(pca)-1, ncol(pca))]
      write.table(sub_pca, file = paste0(plink.fileprefix, '.5_PCA.PCs1-4.outliers_not_removed.tsv'),
                  sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }
  
  # Reads eigenvalues to plot '% explained variance' per component
  eigenval <- scan(dir(pattern = '*.5_PCA.eigenval'))
  if (length(eigenval) == 8000) { # Eigenvalues might be over the limit of plink
    covarmatrix_diag <- scan(dir(pattern = '*.5_PCA.diagonal.covariance')) # To get the total sum of eigenvalues
    pve <- data.frame(PC = 1:length(eigenval), p.v.e = eigenval/sum(covarmatrix_diag)*100)
  } else {
    pve <- data.frame(PC = 1:length(eigenval), p.v.e = eigenval/sum(eigenval)*100)
  }
  
  # Limits the max. PCs to make a more readable barplot
  pve_sub <- pve[1:pve_dim, ] 
  # Plotting pve for each component
  pve_plot <- ggplot(pve_sub, aes(PC, p.v.e)) + geom_bar(stat = 'identity', fill = '#7E52A0') + 
    geom_text(aes(label = round(p.v.e, 2)), position = position_dodge(.9), vjust = -0.5, size = 2.8) +
    ylab("Percentage variance explained") + ggtitle("Percentage Variance Explained of each Principal Component") + 
    theme_update(plot.title = element_text(hjust = 0.5)) +
    theme(axis.title.y =  element_text(size = 14, vjust = 1.5)) +
    scale_x_continuous(breaks = 1:nrow(pve_sub))
  if (save_img) {
    ggsave(plot = pve_plot, path = paste0('./', graph_output_dir), height = 5, width = 8,
           filename = paste0(plink.fileprefix, '_population_pve.png'), device = png)
  }
  
  # We mostly focused on PC1 trough 4, but for some tests we did also checked more
  # Obtains all possible PCs combinations and loops through each metadata level
  components <- combn(seq(1, pca_dim_to_plot), m = 2) 
  for (cluster in c('Center.ShortName', 'PCA.PAM.Ancestry.Cluster', 'PCA.PAM.Ancestry.SubCluster', 'Num.Variants')) {
    if (cluster == 'Center.ShortName') {
      legend <- 'Center'
      color <- 'Paired'
    } else if (cluster == 'PCA.PAM.Ancestry.Cluster'){
      legend <- 'Population'
      color <- 'Dark2'
    } else if (cluster == 'PCA.PAM.Ancestry.SubCluster') {
      legend <- 'Subpopulation'
      color <- 'Dark2'
    } else if (cluster == 'Num.Variants'){
      legend <- 'Number of variants'
    }
    
    # Goes for each PCs combination to plot all PCAs, will return some warnings
    # due to the interaction of dplyr and the variable calling needed to specify
    # aes(), but I wasn't able to find a solution to those
    plot_list <- lapply(1:ncol(components), function(comb) {
      # Extract the two dimensions to compare at each step
      a <- components[, comb][1]
      b <- components[, comb][2]
      
      # To use one of the color palette defined before as it fits the data
      if (cluster != 'Num.Variants') {
        pop_plot <- ggplot(data = pca, aes(pca[[a +1]], pca[[b +1]], col = eval(as.name(paste(cluster))) )) + 
          geom_point(shape = 10, size = 3) +  scale_color_brewer(palette = color, na.value = 'black') +
          #geom_point(shape = 10, size = 3) +  scale_fill_manual(palette = color_blind, na.value = 'black') +
          xlab(paste0('PC', a, ' (', signif(pve$p.v.e[a], 3), '%)')) +
          ylab(paste0('PC', b, ' (', signif(pve$p.v.e[b], 3), '%)')) +
          labs(color = legend, subtitle = plot_subtitle) + ggtitle(plot_title) + 
          theme_update(plot.title = element_text(hjust = 0.5, size = 15)) +
          theme(axis.title.x = element_text(size = 12, vjust = 1),
                axis.title.y = element_text(size = 12, vjust = 1),
                plot.subtitle = element_text(hjust = 0.5),
                panel.background = element_rect(fill = "#F9F9F9",
                                                colour = "#F9F9F9",
                                                size = 0.5, linetype = "solid"))
        print(pop_plot)
        
        # For the number of variants we use a gradient color scheme
      } else {
        pop_plot <- ggplot(data = pca, aes(x = pca[[a +1]], y = pca[[b +1]], col = eval(as.name(paste(cluster))) )) +
          geom_point(shape = 10, size = 3) + scale_color_gradient(low = 'skyblue', high = 'navy' ) +
          xlab(paste0('PC', a, ' (', signif(pve$p.v.e[a], 3), '%)')) +
          ylab(paste0('PC', b, ' (', signif(pve$p.v.e[b], 3), '%)')) +
          labs(color= legend, subtitle = plot_subtitle) + ggtitle(plot_title) +
          theme_update(plot.title = element_text(hjust = 0.5, size = 15)) +
          theme(axis.title.x = element_text(size = 12, vjust = 1),
                axis.title.y = element_text(size = 12, vjust = 1),
                plot.subtitle = element_text(hjust = 0.5),
                panel.background = element_rect(fill = "#F9F9F9",
                                                colour = "#F9F9F9",
                                                size = 0.5, linetype = "solid"))
        print(pop_plot)
      }
      
      # Stores each individual PCA with a dynamic name
      if (save_img) {
        ggsave(plot = pop_plot, path = paste0('./', graph_output_dir),
               height = 1784, width = 3171, unit = 'px', dpi = 300,
               filename = paste0(plink.fileprefix, '_', cluster, '_PC', a, '_PC', b, '.png'),
               device = 'png')
      }
      # Returns each plot stored in a list thanks to lapply()
      return(pop_plot) 
      })
    
    # Remove all individual titles and subtitles to have a better looking 
    # multiplot
    plot_list <- lapply(plot_list, function(x){ x$labels$title <- NULL ; return(x) })
    plot_list <- lapply(plot_list, function(x){ x$labels$subtitle <- NULL ; return(x) })
    
    # Saves all PCs combinations in a single file (ncol and nrow left auto to 'facet_wrap')
    ggsave(plot = wrap_plots(plot_list, guides = 'collect') +
             plot_annotation(title = plot_title, subtitle = plot_subtitle,
                             theme = theme(plot.subtitle = element_text(hjust = 0.5))),
           path = paste0('./', graph_output_dir), height = 2800, width = 5600,
           unit = 'px', dpi = 300,
           filename = paste0(plink.fileprefix, '_', cluster, '_combined_plots.png'),
           device = png)
  }
}

##### Callings of function

## FIRST run the script with 'filter_outliers' = FALSE and then identify in the 
## plots if there's a PC with outliers and under which value. Then run it again
## defining 'outliers_PC', 'outliers_PC_thr' and 'outliers_PC_method' with those 
pcaineitor_pca(start_dir = './plink_PCAplotter_output/0.05/pcas_4d/all_pop/',
               graph_output_dir = 'png_PCA_graphs_postliers', pca_dim_to_plot = 4,
               further_prun = TRUE, opt_prf = '0.05', filter_outliers = TRUE,
               outliers_PC = 'PC3', outliers_PC_thr = 0.025, outliers_PC_method = 'ge',
               plot_title = 'Common Variants VAF 0.05 across all TCGA Populations',
               plot_subtitle = '(Only Exonic Variants - "Outliers" not removed)',
               vcf.loc = '../../../../plink_PCAplotter_inputs_PCA.input.3levels.common_0.05_exonic/PCA.input.3levels.common_0.05_exonic.ALL_POP_plink_PCAplotter_input1.vcf',
               tsv.loc = '../../../../plink_PCAplotter_inputs_PCA.input.3levels.common_0.05_exonic/PCA.input.3levels.common_0.05_exonic.ALL_POP_plink_PCAplotter_input2.tsv',
               pop.loc = '../../../../TCGA_barcode_population_PanCancerSupp_PopInfoPAM.13-01-22.tsv',
               center.loc = '../../../../TCGA.PatientBarcode.SeqCenter.13-01-22.tsv',
               high_ld_region = '../../../../high-LD-regions-hg19-GRCh37.bed')
setwd('../../../../')
