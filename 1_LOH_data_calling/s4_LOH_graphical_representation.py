import os
from sqlite3 import paramstyle
import sys
import csv
import math
from xml.etree.ElementTree import tostringlist
import pandas as pd
import matplotlib.pyplot as plt
from pandas.core.frame import DataFrame
from pandas.core.groupby.generic import AggScalar
import seaborn as sns
from collections import defaultdict

#############
#-- Setup --#
#############

# Please change accordingly to your needs, so everything is tidy
input_dir_loh_metadata = "0_input_files_loh_data/"
input_dir_loh_binary = "3_out_binary_LOH_matrix_R/"
output_dir = "./4_out_LOH_graphical_representation/"
os.makedirs(os.path.dirname(output_dir), exist_ok = True) 


###################
#-- Input Files --#
###################

# Essential Genes List from genome-wide gene essentiality data from
# loss-of-function and CCLE (Cancer Cell Line Encyclopedia) lines,
# Supp. Data 1 from https://www.nature.com/articles/s41467-020-16399-y 

ess_genes = pd.read_excel(os.path.join(input_dir_loh_metadata, '41467_2020_16399_MOESM3_ESM.xlsx'), \
    sheet_name = 0, skiprows = range(5), index_col = 0)  # First five rows are pretty headers, unuseful
ess_genes = ess_genes[ess_genes.columns[-1]]  # Last column is boolean result of essenciality
ess_genes_list = ess_genes[ess_genes == 1].index.tolist() # Len of 1482 instead of 1481 as paper?
                                                          # Apparently 1482 as checked in xlsx


# Survival and clinical information for TCGA samples obtained from TCGA; file obtained from
# https://gitlab.com/ccg-cnio/resources/-/blob/master/TCGA_GAM/resources/MAF_procedure/combined_study_clinical_data.tsv 
patient_info = pd.read_csv(os.path.join(input_dir_loh_metadata, 'Survival_SupplementalTable_S1_20171025_xena_sp_1'), \
    sep = '\t', header = 0, index_col = 0)
cancer_type_per_sample_1 = patient_info[['cancer type abbreviation']] # Patient's cancer type
                                                                      # useful for clustering

patient_info_2 = pd.read_csv(os.path.join(input_dir_loh_metadata, 'combined_study_clinical_data.tsv'), \
    sep = '\t', header = 0, index_col = 2)
cancer_type_per_sample_2 = patient_info_2[['TCGA PanCanAtlas Cancer Type Acronym']] 
cancer_type_per_sample_2.rename(columns=\
    {'TCGA PanCanAtlas Cancer Type Acronym':'cancer type abbreviation'}, inplace = True)    

# We merge both tables and remove some duplicates to get the most out of them
cancer_type_per_sample = pd.concat([cancer_type_per_sample_1,cancer_type_per_sample_2],\
    ignore_index = False, sort = True).reset_index().drop_duplicates(subset = 'index').\
        set_index('index')                                                                
cancer_type_per_sample.index.names = ['samples']

# Solip recommended to merge COAD + READ into one class, similar to 'colorrectal' in
cancer_type_per_sample['cancer type abbreviation'].\
    replace({'COAD':'COADREAD', 'READ':'COADREAD'}, inplace = True)


# This is a somatic drivers obtained from Cancer Gene Census from COSMIC
res_dir = ('./0_input_files_loh_data/')

input_csv = pd.read_csv(os.path.join(res_dir, 'Census_allMon_Feb_14_11_23_46_2022.tsv'),\
    sep = '\t', header = 0, index_col = None)  

somatic_drivers = input_csv.query('Somatic=="yes"')['Gene Symbol'].tolist() 
germline_drivers = input_csv.query('Germline=="yes"')['Gene Symbol'].tolist()

# List of 152 CPGs used in further parts of the TFM, included to compare the 
# LOH frequency with this gene set to check for possible differences between
# the original LOH matrix (this one) and the merged GAM+LOH matrices for 
# the linear models; taking from COSMIC too
cpgs_list = list(pd.read_csv(os.path.join(res_dir, 'CPG.csv'), sep = '\t', header = None, names = ['gene'])['gene'])
cpgs_list = [gene.strip(' ') for gene in cpgs_list] # Remove blank cases just in case

# When filtering with both of the lists later on, I found out that we end up missing 24
# genes from our drivers list.
# 
# Checking by hand, some of the genes were mapped with other names, so we'll CHANGE
# THE CURRENT GENECARDS NAMING OF THOSE TO THE ONES FROM OUR LIST
diff_gene_names = {'H3C2':'HIST1H3B', 'H4C9':'HIST1H4I', 'SEPTIN6':'SEPT6',\
    'SEPTIN9':'SEPT9', 'PIGK':'IGK', 'PTCRA' : 'TRA', 'H3-3B':'H3F3B', 'SEPTIN5' : 'SEPT5',
    'H3-3A':'H3F3A'}

# Why some of the missing genes from the lists are not included? 
#
# HMGN2P46 is a pseudogene
# IGL seems to include a family of immunoglobulin-related pseudogenes and prot. coding
#     the protein coding ones are Immunoglobulin Lamba Like [Num] (IGLL-)
#     IGLL1 and IGLL5 are protein coding ones that were mapped and could be included
# IGH seems to be a similar case, some protein coding genes that we currently have
#     mapped are IGHMBP2 and IGHV3OR16-17
# TRD or other gene symbols do not seem to appear in our mapping
# MTCP1 doesn't seem to have mapped in any sample
# MDS2 does not appear as it is a lncRNA
# PTK6 doesn't seem to have mapped in any sample
# TRP are pseudogenes
# ASPSCR1 doesn't seem to have mapped in any sample
# CARS seems to be non-exonic 
# DUX4L1 are a family of pseudogenes
# MALAT1 are lncRNA
# CRLF2 doesn't seem to have mapped in any sample (chr. X)
# P2RY8 doesn't seem to have mapped in any sample (chr. X)



######################
#-- Main scripting --#
######################

matrix_df = pd.read_csv(os.path.join(input_dir_loh_binary,'loh_matrix_complete.tsv'),\
     sep = '\t', header = 0, index_col = 1)                  
matrix_df = matrix_df.drop(matrix_df.columns[0], axis = 1).astype(int)

# We change the name of the genes mentioned a few section before
matrix_df.rename(columns = diff_gene_names, inplace = True)

# We can now choose the filter to use on the matrix:
def matrix_filter(df = matrix_df, filt = 'essential', write_table = False):
    if filt == 'essential':
        # List of 1481 genes from Nichols et al. (2020)
        matrix_filtered_df = df.filter(items = ess_genes_list, axis = 1)

    elif filt == 'somatic':
        # From 'Cancer_gene_census_3sep_2021_94.tsv'
        matrix_filtered_df = df.filter(items = somatic_drivers, axis = 1)

    elif filt == 'germline':
        # From 'Cancer_gene_census_3sep_2021_94.tsv'
        matrix_filtered_df = df.filter(items = germline_drivers, axis = 1)

    elif filt == 'all_drivers':
        # List of 723 genes from 'Cancer_gene_census_3sep_2021_94.tsv'
        drivers_list = list(set(sorted(somatic_drivers + germline_drivers))) 
        matrix_filtered_df = df.filter(items = drivers_list, axis = 1) 

    elif filt == 'essential+drivers':
        all_list = list(set(sorted(somatic_drivers + germline_drivers + ess_genes_list)))
        matrix_filtered_df = df.filter(items = all_list, axis = 1)

    elif filt == 'cpgs':
        # List of 152 CPGs provided by the lab; used to compare with 
        # further steps of the project that compute LOH freq. for this set
        matrix_filtered_df = df.filter(items = cpgs_list, axis = 1)
    
    elif filt == 'all_genes':
        # Keeps all available genes 
        matrix_filtered_df = df
                     
    if write_table:         
        matrix_filtered_df.to_csv(os.path.join(output_dir, \
            '{}_filtered_loh_booltable.tsv'.format(filt)), sep = '\t') 

    return matrix_filtered_df, filt # Returns 'filt' for file suffixes 


##############################################
#-- Graphical representation of LOH events --#
##############################################

#########   ----------------------------------------------------- 
#-{ 1 }-#  |{ LOH mean|median events in samples per cancer type }|
#########   -----------------------------------------------------

# We filter the matrix to get similar results to those of the reference
matrix_essenfiltered_df, filt = matrix_filter(matrix_df)
matrix_allgenes_df, filt = matrix_filter(matrix_df, filt = 'all_genes')  
matrix_driversfiltered_df, filt = matrix_filter(matrix_df, filt = 'all_drivers') 
 
# Now we count all loh events per patient and add their cancer type for further analysis

def event_counter_per_sample(dfs, what = 'proportion'):  

    counts_df = pd.DataFrame(index=dfs.index) 

    if what == 'count':
        counts_df['total_loh_events'] = dfs.sum(axis = 1) 

    elif what == 'proportion':
        counts_df['proportion_loh'] = dfs.sum(axis = 1).div(dfs.shape[1]).mul(100).round(2) # Divides LOH by total genes

    counts_df = counts_df.join(cancer_type_per_sample)

    return counts_df

loh_proportion_all = event_counter_per_sample(matrix_allgenes_df)
loh_proportion_drivers = event_counter_per_sample(matrix_driversfiltered_df)

def graphical_repr_samples(df = loh_proportion_drivers, filt = 'all_drivers', what = 'proportion', stat = 'mean',\
    saveimage = False, legend_show = False, filesufx = '{}_filtered'.format(filt)):
    
    if what == 'count':
        column = 'total_loh_events'
    elif what == 'proportion':
        column = 'proportion_loh'

    stat = stat.lower()  # Can choose between mean or median to order the results 
    if stat == 'mean':
        df_graph = df.groupby('cancer type abbreviation', sort = False) \
            [column].mean().sort_values(ascending = False)
    elif stat == 'median':
        df_graph = df.groupby('cancer type abbreviation', sort = False) \
            [column].median().sort_values(ascending = False)

    print(df_graph)
    
    plt.figure(figsize=(15,8))
    ax = sns.stripplot(x = 'cancer type abbreviation', y = column,\
        data = df, jitter = 0.1, order = df_graph.index)

    ax.set_title('LOH events in {} genes per sample by TCGA tumor type'.format(filesufx))
    ax.set_xlabel('TCGA tumor type', labelpad = 15)
    ax.set_ylabel('# of LOH events', labelpad = 15)
    plt.xticks(rotation = 90)

    # Adds mean o median bar to each tumor type
    _ = [ax.hlines(y, xmin= i-.25, xmax = i+.25, zorder=3, colors = 'black', label = stat) \
        for i, y in df_graph.reset_index()[column].items()]

    # For mean o median to appear in a legend
    if legend_show:
        handles,l = ax.get_legend_handles_labels()
        labels = [stat + 's']
        plt.legend(handles = handles[0:1], labels = labels)

    if saveimage: # To store the results
        plt.subplots_adjust(left = 0.085, right = 0.95, top = 0.915, bottom = 0.14)
        filename = (os.path.join(output_dir, \
            'loh_events_{}_{}_tsv.svg'.format(stat, filesufx)))
        plt.savefig(fname = filename, format = 'svg', dpi = 300)

    plt.show()

# Saving plots for all_data_by_genes_whitelisted.tsv LOH events
graphical_repr_samples(stat='mean')
graphical_repr_samples(stat='median', saveimage = True) # Really similar to the reference

graphical_repr_samples(df = loh_proportion_all,  filt = 'all_genes', stat = 'median', saveimage = False)



#########   ------------------------------------------- 
#-{ 2 }-#  |{ LOH frequency in a gene per cancer type }|
#########   -------------------------------------------

# For the example we're going to use the somatic & germline cancer drivers list,
matrix_driversfiltered_df, filt = matrix_filter(matrix_df, filt = 'all_drivers')
matrix_allgenes_df, filt = matrix_filter(matrix_df, filt = 'all_genes')  

# Now we count all loh events per gene and add their cancer type for further analysis

def event_counter_per_gene(dfs, count = 'total', filename = 'drivergenes',\
    write_table = False, round = 3):  
    # Takes aout the information of the samples, and aggregates and sums LOH events
    # in each level for all cancer types
    
    tmp_df = dfs.join(cancer_type_per_sample).reset_index().iloc[:,1:]\
        .set_index('cancer type abbreviation')
    # We lose some of the samples info as we don't have info on their cancer type

    counts_df = tmp_df.groupby('cancer type abbreviation').sum()

    if count == 'freq':        
        # Total of patients with each cancer type
        total_cancer = tmp_df.index.value_counts(sort = False)
        counts_df = round(counts_df.div(total_cancer, axis = 0) * 100, round)
        
    if write_table:
       counts_df.to_csv(os.path.join(output_dir, \
            'loh_{}_events_{}_.tsv'.format(count, filename)), sep = '\t')
    return counts_df

loh_freq_cancertype = event_counter_per_gene(dfs = matrix_allgenes_df, count = 'freq', filename = 'all_genes')


def pancancer_freq(dfs, decimal = 3):  
    # To give an average freq. in pan-cancer takes out the information of the samples
    # and does not include the cancer type
    
    total_loh_gene = dfs.join(cancer_type_per_sample).reset_index().iloc[:,1:].sum()
    num_samples = dfs.shape[0] - 1
    freq_gene = total_loh_gene.div(num_samples, axis = 0)
    total_freq = freq_gene.sum()

    return round( (total_freq/(dfs.shape[1])) *100, decimal)  

print(f'{pancancer_freq(matrix_driversfiltered_df)}%') # Should return ~15.xx%