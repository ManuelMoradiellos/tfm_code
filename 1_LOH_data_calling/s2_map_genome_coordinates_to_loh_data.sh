#!/bin/bash

# Script that takes the annotated genome and extracts the 
# names of all protein coding genes and their chr. coordinates


# First a python script was created to take the LOH calling using
# ABSOLUTE output (from reference) to extract the coordinates to
# be mapped.
# Those need to be ordered it as a .BED file
awk '{first = $1; $1=""; print $0, first}' ./1_LOH_calling_which_genes_to_map/to_map_nonNA.tsv > ./1_LOH_calling_which_genes_to_map/to_map_nonNA.bed

# Next filtering to only keep protein coding genes
awk -F '[\t]' '$10~ /protein_coding/ {print}' ./0_input_files_genome_annotation/gencode.v38.gene.bed > ./0_input_files_genome_annotation/gencode.v38.protein_coding.gene.bed 

# Now we just need the coords and gene name
awk -F '[\t;]' '{print $1,$2,$3,$12}' ./0_input_files_genome_annotation/gencode.v38.protein_coding.gene.bed > ./0_input_files_genome_annotation/prot_cod_geneCleaned.bed 
# Clean the file up - there are spaces as field separators
awk -F '[ ]' '{print $1,$2,$3,$6}' ./0_input_files_genome_annotation/prot_cod_geneCleaned.bed > ./0_input_files_genome_annotation/prot_cod_geneMinimal.bed

# Remove quotes from gene names
sed 's/"//g' ./0_input_files_genome_annotation/prot_cod_geneMinimal.bed > ./0_input_files_genome_annotation/prot_cod_final.bed
sed -i '1 i\#Chromosome\t#Start\t#End\t#Gene' ./0_input_files_genome_annotation/prot_cod_final.bed # Add header
sed -i 's/ /\t/g' ./0_input_files_genome_annotation/prot_cod_final.bed # Tab delimiter

# I do the mapping of all genes in the chr regions using bedtools intersect
bedtools intersect -a to_map_nonNA.bed -b ./0_input_files_genome_annotation/prot_cod_final.bed -wa -wb | awk 'BEGIN {FS=OFS="\t"} {print $5, $9, $4}' > ./0_input_files_genome_annotation/mapped_samples.tsv