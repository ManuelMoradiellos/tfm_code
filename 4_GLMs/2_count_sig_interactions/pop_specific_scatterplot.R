library(tidyverse)
library(ggplot2)
library(ggpubr)
library(patchwork)

class = 'allclass'
tissue = 'Pancancer'
complete_table <- read.delim("/local/mmoradiellos/work_stuff/5_gam_loh_merges_model_inputs/bck_snakemake_gamloh_models/out_all_genes_test_pop_size_corrected/merged_outputs_correct_pop_size/allclass/allclass_glm_results/all_genes/fdr_complete_tables/all_genes_Pancancer_FDR0.2_Germ0.1_withFDR_allclass.tsv")
sig_counts <- read_tsv('/local/mmoradiellos/work_stuff/5_gam_loh_merges_model_inputs/bck_snakemake_gamloh_models/out_all_genes_test_pop_size_corrected/merged_outputs_correct_pop_size/allclass/allclass_glm_results/all_genes/sig_counts/all_genes_significant_gene_names_allclass.tsv')

# Extract which genes were significant in 2-way
if (class == 'clinvar' | class == 'plof') {
  sig_counts <- sig_counts %>% select(Cancer_type, FDR_0.2.GermFreq_0.1) %>% drop_na() # In some cases clinvar and plof had more sig. results for Freq > 0.5 due to having less sample to FDR correct 
} else {
  sig_counts <- sig_counts %>% select(Cancer_type, FDR_0.2.GermFreq_0.5) %>% drop_na()
} 
genes_to_extract <- unique(unlist(strsplit(as.vector(sig_counts %>% filter(Cancer_type == tissue) %>% pull(colnames(sig_counts)[2])), split = ',')))


# AFR #
complete_table$fdr_afr <-  p.adjust(complete_table$germline_mut1.loh1.popAFR_p_value, method = 'BH')
complete_table$twoway_sig <- ifelse(complete_table$Gene %in% genes_to_extract, 'yes', 'no')

AFR_vs_EUR <- ggplot(complete_table, aes(x =germline_mut1.loh1_estimate, y=germline_mut1.loh1.popAFR_estimate, size = abs(Freq_Germline),
                                    fill = -log(fdr_afr+0.00001, 10))) +
  geom_point(shape=21, alpha = 0.6) +
  theme_light() +
  scale_fill_gradient2(low = "grey", mid = "white", high = "red", limits = c(0, 6)) +
  geom_abline(slope = 1, lty=2, color = 'grey') +
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 3)) +
  labs(fill = '-log10FDR', size = 'freq') +
  ggrepel::geom_text_repel(aes(label = ifelse(twoway_sig == 'yes', as.character(Gene),'')),
                           size = 6, box.padding = 0.5, max.overlaps = Inf,
                           nudge_x = .15,
                           nudge_y = .5,
                           segment.curvature = -1e-20,
                           arrow = arrow(length = unit(0.015, "npc")))  
# AMR
complete_table$fdr_amr <-  p.adjust(complete_table$germline_mut1.loh1.popAMR_p_value, method = 'BH')

AMR_vs_EUR <- ggplot(complete_table, aes(x =germline_mut1.loh1_estimate, y=germline_mut1.loh1.popAMR_estimate, size = abs(Freq_Germline),
                                      fill = -log(fdr_amr+0.00001, 10))) +
  geom_point(shape=21, alpha = 0.6) +
  theme_light() +
  scale_fill_gradient2(low = "grey", mid = "white", high = "red", limits = c(0, 6)) +
  geom_abline(slope = 1, lty=2, color = 'grey') +
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 3)) +
  labs(fill = '-log10FDR', size = 'freq') +
  ggrepel::geom_text_repel(aes(label = ifelse(twoway_sig == 'yes', as.character(Gene),'')),
                           size = 6, box.padding = 0.5, max.overlaps = Inf,
                           nudge_x = .15,
                           nudge_y = .5,
                           segment.curvature = -1e-20,
                           arrow = arrow(length = unit(0.015, "npc"))) 

# EAS
complete_table$fdr_eas <-  p.adjust(complete_table$germline_mut1.loh1.popASIAN_p_value, method = 'BH')
EAS_vs_EUR <- ggplot(complete_table, aes(x =germline_mut1.loh1_estimate, y=germline_mut1.loh1.popASIAN_estimate, size = abs(Freq_Germline),
                                      fill = -log(fdr_eas+0.00001, 10))) +
  geom_point(shape=21, alpha = 0.6) +
  theme_light() +
  scale_fill_gradient2(low = "grey", mid = "white", high = "red", limits = c(0, 6)) +
  geom_abline(slope = 1, lty=2, color = 'grey') +
  labs(fill = '-log10FDR', size = 'freq') +
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 3)) +
  ggrepel::geom_text_repel(aes(label = ifelse(twoway_sig == 'yes', as.character(Gene),'')),
                           size = 6, box.padding = 0.5, max.overlaps = Inf,
                           nudge_x = .15,
                           nudge_y = .5,
                           segment.curvature = -1e-20,
                           arrow = arrow(length = unit(0.015, "npc"))) 

#####
afr_high <- sort(-log(complete_table$fdr_afr+0.00001, 10), decreasing = T)[2]
AFR_scatter <- ggplot(complete_table, aes(x=germline_mut1.loh1.popAFR_estimate,
                                 y=-log(germline_mut1.loh1.popAFR_p_value+0.00001, 10),
                                 size = abs(Freq_Germline), fill = -log(fdr_afr+0.00001, 10)))+
  geom_point(shape=21) +
  theme_light() +
  scale_fill_gradient2(low = "grey", mid = "white", high = "red", limits = c(0, 6))+
  labs(fill = '-log10FDR', size = 'freq') +
  geom_hline(yintercept=0.69897, linetype="dashed", color = "grey") +
  ggrepel::geom_text_repel(aes(label = ifelse( (-log(fdr_afr+0.00001, 10) >= afr_high) , as.character(Gene),'')),
                           size = 6, box.padding = 0.5, max.overlaps = Inf,
                           nudge_x = .15,
                           nudge_y = .5,
                           arrow = arrow(length = unit(0.015, "npc"))) 
  

amr_high <- sort(-log(complete_table$fdr_amr+0.00001, 10), decreasing = T)[2]
AMR_scatter <- ggplot(complete_table, aes(x=germline_mut1.loh1.popAMR_estimate,
                                          y=-log(germline_mut1.loh1.popAMR_p_value+0.00001, 10),
                                          size = abs(Freq_Germline), fill = -log(fdr_amr+0.00001, 10)))+
  geom_point(shape=21) +
  theme_light() +
  scale_fill_gradient2(low = "grey", mid = "white", high = "red", limits = c(0, 6))+
  labs(fill = '-log10FDR', size = 'freq') +
  geom_hline(yintercept=0.69897, linetype="dashed", color = "grey") +
  ggrepel::geom_text_repel(aes(label = ifelse( (-log(fdr_amr+0.00001, 10) >= amr_high) , as.character(Gene),'')),
                           size = 6, box.padding = 0.5, max.overlaps = Inf,
                           nudge_x = .15,
                           nudge_y = .5,
                           arrow = arrow(length = unit(0.015, "npc"))) 

eas_high <- sort(-log(complete_table$fdr_eas+0.00001, 10), decreasing = T)[2]
EAS_scatter <- ggplot(complete_table, aes(x=germline_mut1.loh1.popASIAN_estimate,
                                          y=-log(germline_mut1.loh1.popASIAN_p_value+0.00001, 10),
                                          size = abs(Freq_Germline), fill = -log(fdr_eas+0.00001, 10)))+
  geom_point(shape=21) +
  theme_light() +
  scale_fill_gradient2(low = "grey", mid = "white", high = "red", limits = c(0, 6))+
  labs(fill = '-log10FDR', size = 'freq') +
  geom_hline(yintercept=0.69897, linetype="dashed", color = "grey") +
  ggrepel::geom_text_repel(aes(label = ifelse( (-log(fdr_eas+0.00001, 10) >= eas_high) , as.character(Gene),'')),
                           size = 6, box.padding = 0.5, max.overlaps = Inf,
                           nudge_x = .15,
                           nudge_y = .5,
                           arrow = arrow(length = unit(0.015, "npc"))) 


AFR_vs_EUR + AMR_vs_EUR + EAS_vs_EUR + AFR_scatter + AMR_scatter + EAS_scatter + plot_layout(ncol = 3, guides = 'collect')



FDR =p.adjust(complete_table$germline_mut1.loh1.popASIAN_p_value, method = 'BH')
EAS1 <- ggplot(complete_table, aes(x =germline_mut1.loh1_estimate, y=germline_mut1.loh1.popASIAN_estimate,
                          size = abs(Freq_Germline), fill = -log(FDR+0.00001, 10)))+
  geom_point(shape=21) + theme_bw()+
  scale_fill_gradient2(low = "grey", mid = "white", high = "red")+
  geom_abline(slope = 1, lty=2, color = 'grey')