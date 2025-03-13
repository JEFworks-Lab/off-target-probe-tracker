

############################################################################################################

# Compare Visium and Xenium data


############################################################################################################

# import libraries
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
import numpy as np
import matplotlib.patches as mpatches
import pandas as pd
import scanpy as sc
from PIL import Image
import pickle
import scipy
from collections import defaultdict, Counter
import seaborn as sns
import os
import scipy.sparse as sp
import scipy.stats as stats
import ast
from utils import rerasterize_patches, plotRaster, rasterizeGeneExpression_topatches_basedoncenters
from adjustText import adjust_text


############################################################################################################


### read in xenium and visium data ###

# combined data
adata_xenium = sc.read_h5ad('/home/caleb/Desktop/improvedgenepred/data/breastcancer_sample1_rep1_aligned_toxeniumimage/xenium_data_full.h5ad')
adata_visium = sc.read_h5ad('/home/caleb/Desktop/improvedgenepred/data/breastcancer_sample1_rep1_aligned_tovisiumimage/visium_data_full.h5ad')

# make .X a csr matrix
adata_xenium.X = scipy.sparse.csr_matrix(adata_xenium.X)
adata_visium.X = scipy.sparse.csr_matrix(adata_visium.X)

# add array for gene expression
adata_xenium.X_array = pd.DataFrame(adata_xenium.X.toarray(), index=adata_xenium.obs.index)
adata_visium.X_array = pd.DataFrame(adata_visium.X.toarray(), index=adata_visium.obs.index)

# patches
with open('/home/caleb/Desktop/improvedgenepred/data/breastcancer_sample1_rep1_aligned_tovisiumimage/visium_patches_full.pkl', 'rb') as f:
    aligned_visium_dictionary = pickle.load(f)

with open('/home/caleb/Desktop/improvedgenepred/data/breastcancer_sample1_rep1_aligned_toxeniumimage/xenium_patches_full.pkl', 'rb') as f:
    aligned_xenium_dictionary = pickle.load(f)


# scale the data
scaling_factor = 1
for i in aligned_visium_dictionary:
    aligned_visium_dictionary[i].X = sc.pp.log1p(np.round(aligned_visium_dictionary[i].X * scaling_factor))

# scale the data
scaling_factor = 1
for i in aligned_xenium_dictionary:
    aligned_xenium_dictionary[i].X = sc.pp.log1p(np.round(aligned_xenium_dictionary[i].X * scaling_factor))


# get full gene list
full_gene_list = list(adata_visium.var_names)

# subset the data to intersecting genes
intersecting_genes = list(set(adata_xenium.var_names) & set(adata_visium.var_names))
print(len(intersecting_genes))
# subset the data
adata_xenium_sub = adata_xenium[:, intersecting_genes]
adata_visium_sub = adata_visium[:, intersecting_genes]

# get the total counts of genes
adata_visium_sub.var['total_counts_sub'] = list(adata_visium_sub.X.sum(axis=0).A1)
adata_xenium_sub.var['total_counts_sub'] = list(adata_xenium_sub.X.sum(axis=0).A1)

# log the total counts
adata_visium_sub.var['total_counts_sub_log'] = np.log(adata_visium_sub.var['total_counts_sub']+1)
adata_xenium_sub.var['total_counts_sub_log'] = np.log(adata_xenium_sub.var['total_counts_sub']+1)



############################################################################################################


### read in off-target data ###

# read in text file

# read in genes of interest
bp_diff = "0" # NOTE: set this to 0, 5, or 10 in string form
offtarget_genes = pd.read_csv(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/gene_summary_out_{bp_diff}bp_nucmer_final.csv")
offtarget_genes

# get the list of genes in xenium probes
xenium_probes = pd.read_csv("/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/xenium_probes_list.csv")
list(xenium_probes['Gene'])


# adding diff to adata_visium_sub.var
gene_summary_tmp = offtarget_genes.copy()
gene_summary_tmp = gene_summary_tmp.set_index('Gene')
# merge df
merged_df = adata_visium_sub.var.join(gene_summary_tmp, how="left")
# add another column for if genes are in xenium_probes
merged_df['dowehaveprobes'] = merged_df.index.isin(list(xenium_probes['Gene']))

# make adata.var the merged
adata_visium_sub.var = merged_df


# getting what genes are in visium
colors = []
for gene in adata_visium_sub.var.index:

    # gene = "ZEB2"
    aligned_to_list = (
        [] if pd.isna(adata_visium_sub.var.loc[gene, "aligned_to"])
        else adata_visium_sub.var.loc[gene, "aligned_to"].split(',')
    )
    
    genes_in_vis = adata_visium_sub.var.loc[gene, "genes_in_vis"]

    # more than one gene in `aligned_to` OR a single gene that is different from the index
    is_offtarget = len(aligned_to_list) > 1 or (len(aligned_to_list) == 1 and aligned_to_list[0] != gene)

    # off-targets not visium
    is_invisium = is_offtarget and gene != genes_in_vis

    # check if gene is in xenium_probes
    is_aprobe = adata_visium_sub.var.loc[gene, "dowehaveprobes"]

    # assign colors
    if is_invisium:
        colors.append('darkgreen')
    elif is_offtarget:
        colors.append('darkblue')
    elif not is_aprobe:
        colors.append('purple')
    else:
        colors.append('lightgray') 

# double check
print(pd.Series(colors).value_counts())




# get data from your AnnData objects
x_vis = adata_visium_sub.var['total_counts_sub_log']
y_xen = adata_xenium_sub.var['total_counts_sub_log']

# Convert colors list to an array
colors_array = np.array(colors)

# plot
plt.figure(figsize=(10, 6))

# Plot reference lines
plt.plot([min(y_xen), max(y_xen)],
         [min(y_xen), max(y_xen)],
         color='darkred', linewidth=2, linestyle='--')

x_vals = np.linspace(min(x_vis), max(x_vis) + 0.5, 100)
y_vals = np.log(10) + x_vals
# plt.plot(x_vals, y_vals, color='green', linewidth=2)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# Define the desired plotting order
plot_order = ['lightgray', 'purple', 'darkblue', 'darkgreen']

# Plot each group in order
for col in plot_order:
    idx = np.array(x_vis.index)[np.where(colors_array == col)]
    if col in ['lightgray']:
        plt.scatter(x_vis[idx], y_xen[idx], c=col, marker="o", facecolors='none')
    else:
        if col == 'darkblue':
            plt.scatter(x_vis[idx], y_xen[idx], c=col, marker="o", facecolors='none')
        elif col == 'purple':
            plt.scatter(x_vis[idx], y_xen[idx], c=col, marker="o", facecolors='none')
        elif col == 'darkgreen':
            plt.scatter(x_vis[idx], y_xen[idx], c=col, marker="o", facecolors='none')


# Label points that are darkgreen 
texts = []
for i, gene in enumerate(adata_visium_sub.var_names):
    if colors[i] == 'darkgreen':
        txt = plt.text(x_vis[i], y_xen[i], gene, fontsize=8)
        texts.append(txt)

# Also add HDC 
hdc_idx = np.where(adata_visium_sub.var_names == 'HDC')[0]
if hdc_idx.size > 0:
    txt = plt.text(x_vis[hdc_idx][0], y_xen[hdc_idx][0], 'HDC', fontsize=8)
    texts.append(txt)

# adjust text
adjust_text(texts, arrowprops=dict(arrowstyle="-", lw=0.5), force_points=0.5)

# Custom legend
red_patch = mpatches.Patch(color='darkgreen', label='Off-target probes in Visium')
blue_patch = mpatches.Patch(color='darkblue', label='Off-target probes not in Visium')
green_patch = mpatches.Patch(color='purple', label='Custom Probes')
gray_patch = mpatches.Patch(color='lightgray', label='No off-targets')
plt.legend(handles=[red_patch, blue_patch, green_patch, gray_patch])

sns.despine()
# plt.axis("off")
plt.savefig("/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/final_figures/figure1/visium_xenium_comparison.svg", dpi=300, bbox_inches='tight')


############################################################################################################


### plot the gene expression of the off-targets ###

# combine gene expression
bp_diff = "6" # NOTE: set this to 0, 5, or 10 in string form
offtarget_genes = pd.read_csv(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/gene_summary_out_{bp_diff}bp_nucmer_new.csv")
offtarget_genes

# log transform the data
sc.pp.log1p(adata_xenium)
sc.pp.log1p(adata_visium)


# New resolution that doesn;t cause issues
new_resolution_xenium = 275

# Rerasterize the patches
aligned_visium_dictionary_rerastered = rerasterize_patches(aligned_visium_dictionary, new_resolution_xenium)
aligned_xenium_dictionary_rerastered = rerasterize_patches(aligned_xenium_dictionary, new_resolution_xenium)


# plotRaster(adata_xenium.uns["spatial"], aligned_xenium_dictionary_rerastered, color_by='gene_expression', gene_name='GFAP', if_vis=False)
plotRaster(adata_visium.uns["spatial"], aligned_visium_dictionary_rerastered, color_by='gene_expression', gene_name='GFAP', if_vis=False)


# turn off interactive mode so figures can be saved without displaying
plt.ioff()

# drfine your output folder
output_base = "/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/gene_plots_" + bp_diff + "bp"
os.makedirs(output_base, exist_ok=True)

# define the gene lists
gene_list_visium = list(adata_visium.var.index)
gene_list_xenium = list(adata_xenium.var.index)

# subset to get only few genes, all new genes for 5bp but only two in dataset
# offtarget_genes[offtarget_genes['Gene'].isin(['ACTG2', 'APOC1', 'CD83', 'CXCL5', 'CEACAM8', 'ENAH', 'FLNB'])]

# offtarget_genes[offtarget_genes['Gene'].isin(['CXCL5', 'CEACAM8'])]

# offtarget_genes.iloc[4:,:]

# Loop over each row of dataframe
for idx, row in offtarget_genes[offtarget_genes['Gene'].isin(['CEACAM8'])].iterrows():
    # grab main gene
    main_gene = row["Gene"]
    print(main_gene)
    # skip if gene is not in visium gene list
    if main_gene not in gene_list_visium:
        print(f"Skipping {main_gene} because it is not in the Visium gene list.")
        continue
    genes_in_vis = row["genes_in_vis"]  # This is a dictionary (e.g., {'HOXD9': 7, 'HOXB9': 1})
    
    # Create a folder for this gene inside the output folder.
    gene_folder = os.path.join(output_base, main_gene)
    os.makedirs(gene_folder, exist_ok=True)
    
    # Plot Xenium data for the main gene
    fig = plotRaster(adata_xenium.uns["spatial"], aligned_xenium_dictionary_rerastered, color_by='gene_expression', gene_name=main_gene, if_vis=False)
    # Save the figure in the gene folder
    xenium_filename = os.path.join(gene_folder, f"{main_gene}_xenium.png")
    fig.savefig(xenium_filename, bbox_inches='tight', dpi=150)
    plt.close(fig)
    
    # Plot Visium data for the main gene.
    fig = plotRaster(adata_visium.uns["spatial"], aligned_visium_dictionary_rerastered, color_by='gene_expression', gene_name=main_gene, if_vis=True)
    visium_main_filename = os.path.join(gene_folder, f"{main_gene}_visium.png")
    fig.savefig(visium_main_filename, bbox_inches='tight', dpi=150)
    plt.close(fig)
    
    # For each off-target gene
    # Plot Visium data for the off-target
    if ',' in genes_in_vis:
        # for off_target in list(ast.literal_eval(genes_in_vis).keys()):
        for off_target in genes_in_vis.split(','):
            if off_target != main_gene:
                fig = plotRaster(adata_visium.uns["spatial"], aligned_visium_dictionary_rerastered, color_by='gene_expression', gene_name=off_target, if_vis=True)
                off_target_filename = os.path.join(gene_folder, f"{main_gene}_visium_offtarget_{off_target}.png")
                fig.savefig(off_target_filename, bbox_inches='tight', dpi=150)
                plt.close(fig)

# print
print("Plots saved for all genes.")
# turn on interactive mode
plt.ion()


############################################################################################################


### plot aggregated gene expression ###

# combine gene expression
bp_diff = "6" # NOTE: set this to 0, 5, or 10 in string form
offtarget_genes = pd.read_csv(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/gene_summary_out_{bp_diff}bp_nucmer_new.csv")
offtarget_genes


# need to rerasterize the patches

# combined data
adata_xenium = sc.read_h5ad('/home/caleb/Desktop/improvedgenepred/data/breastcancer_sample1_rep1_aligned_tovisiumimage/xenium_data_full.h5ad')
adata_visium = sc.read_h5ad('/home/caleb/Desktop/improvedgenepred/data/breastcancer_sample1_rep1_aligned_tovisiumimage/visium_data_full.h5ad')

# make .X a csr matrix
adata_xenium.X = scipy.sparse.csr_matrix(adata_xenium.X)
adata_visium.X = scipy.sparse.csr_matrix(adata_visium.X)

# add array for gene expression
adata_xenium.X_array = pd.DataFrame(adata_xenium.X.toarray(), index=adata_xenium.obs.index)
adata_visium.X_array = pd.DataFrame(adata_visium.X.toarray(), index=adata_visium.obs.index)

# combine gene expression
for i in range(0, len(offtarget_genes)):
    gene = offtarget_genes['Gene'][i]
    # off-targets
    off_targets = offtarget_genes.loc[i, 'genes_in_vis'].split(',')
    # check if gene is in both adata
    if gene not in adata_visium.var_names or gene not in adata_xenium.var_names:
        print(f"Skipping {gene} because it is not in both datasets.")
        continue
    # check if off-targets are in both adata, if it isnt remove that gene from list
    valid_off_targets = [g for g in off_targets if g in adata_visium.var_names]
    removed_genes = [g for g in off_targets if g not in valid_off_targets]
    if removed_genes:
        print(f"Removed off-target genes for {gene}: {removed_genes}")
    if not valid_off_targets:
        continue
    off_targets = valid_off_targets
    # replace gene with aggregated gene
    # get the sum of the off-targets
    adata_visium[:, gene].X = sp.csr_matrix(adata_visium[:, off_targets].X.sum(axis=1))  # add to gene list
    adata_visium.var_names_make_unique()
    

# redo to get the patches for visium
adata_patches_visium_agg = rasterizeGeneExpression_topatches_basedoncenters(adata_visium.uns['spatial'], adata_visium, adata_visium.obsm['spatial'], patch_size=250, aggregation='sum', visium=False)
len(adata_patches_visium_agg)

# scale the data
scaling_factor = 1
for i in adata_patches_visium_agg:
    adata_patches_visium_agg[i].X = sc.pp.log1p(np.round(adata_patches_visium_agg[i].X * scaling_factor))


# log transform the data
sc.pp.log1p(adata_visium)


# New resolution that doesn;t cause issues
new_resolution_xenium = 275

# Rerasterize the patches
aligned_visium_dictionary_rerastered_agg = rerasterize_patches(adata_patches_visium_agg, new_resolution_xenium)


# turn off interactive mode so figures can be saved without displaying
plt.ioff()

# drfine your output folder
output_base = "/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/gene_plots_" + bp_diff + "bp"
os.makedirs(output_base, exist_ok=True)

# define the gene lists
gene_list_visium = list(adata_visium.var.index)
gene_list_xenium = list(adata_xenium.var.index)

# offtarget_genes.iloc[17:,:]

# Loop over each row of dataframe
for idx, row in offtarget_genes[offtarget_genes['Gene'].isin(['CEACAM8'])].iterrows():
    # grab main gene
    main_gene = row["Gene"]
    print(main_gene)
    # skip if gene is not in visium gene list
    if main_gene not in gene_list_visium:
        print(f"Skipping {main_gene} because it is not in the Visium gene list.")
        continue
    genes_in_vis = row["genes_in_vis"]  # This is a dictionary (e.g., {'HOXD9': 7, 'HOXB9': 1})
    
    # Create a folder for this gene inside the output folder
    gene_folder = os.path.join(output_base, main_gene)
    os.makedirs(gene_folder, exist_ok=True)
    
    # Plot Visium agg data for the main gene
    fig = plotRaster(adata_visium.uns["spatial"], aligned_visium_dictionary_rerastered_agg, color_by='gene_expression', gene_name=main_gene, if_vis=True)
    visium_main_filename = os.path.join(gene_folder, f"{main_gene}_visium_agg.png")
    fig.savefig(visium_main_filename, bbox_inches='tight', dpi=150)
    plt.close(fig)

# print
print("Plots saved for all genes.")
# turn on interactive mode
plt.ion()


############################################################################################################


### Finding Pearson Correlation ###

# NOTE: re load data so it isnt normalized
# NOTE: can change from rep1 to rep2

# combine gene expression
bp_diff = "10" # NOTE: set this to 0, 5, or 10 in string form
offtarget_genes = pd.read_csv(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/gene_summary_filtered_{bp_diff}bp.csv")
offtarget_genes

# # combined data
# adata_xenium = sc.read_h5ad('/home/caleb/Desktop/improvedgenepred/data/breastcancer_sample1_rep1_aligned_tovisiumimage/xenium_data_full.h5ad')
# adata_visium = sc.read_h5ad('/home/caleb/Desktop/improvedgenepred/data/breastcancer_sample1_rep1_aligned_tovisiumimage/visium_data_full.h5ad')
adata_xenium = sc.read_h5ad('/home/caleb/Desktop/improvedgenepred/data/breastcancer_sample1_rep2_aligned/xenium_data_full.h5ad')
adata_visium = sc.read_h5ad('/home/caleb/Desktop/improvedgenepred/data/breastcancer_sample1_rep2_aligned/visium_data_full.h5ad')




# make .X a csr matrix
adata_xenium.X = scipy.sparse.csr_matrix(adata_xenium.X)
adata_visium.X = scipy.sparse.csr_matrix(adata_visium.X)

# add array for gene expression
adata_xenium.X_array = pd.DataFrame(adata_xenium.X.toarray(), index=adata_xenium.obs.index)
adata_visium.X_array = pd.DataFrame(adata_visium.X.toarray(), index=adata_visium.obs.index)

# find correlation between xenium and visium
# calculate pearson correlation between xenium and visium per gene
pearson_corr = []
for gene in offtarget_genes['Gene']:
    # check if gene is in both adata
    if gene not in adata_visium.var_names or gene not in adata_xenium.var_names:
        print(f"Skipping {gene} because it is not in both datasets.")
        continue
    corr = stats.pearsonr(np.log(adata_visium[:, gene].X.toarray().flatten()+1), np.log(adata_xenium[:, gene].X.toarray().flatten()+1))[0]
    pearson_corr.append(corr)


# now do the same but combine the gene expression of off-targets
pearson_corr_agg = []
for i in range(0, len(offtarget_genes)):
    print(i)
    gene = offtarget_genes['Gene'][i]
    # off-targets
    off_targets = list(ast.literal_eval(offtarget_genes.loc[i, 'genes_in_vis']).keys())
    # check if gene is in both adata
    if gene not in adata_visium.var_names or gene not in adata_xenium.var_names:
        print(f"Skipping {gene} because it is not in both datasets.")
        continue
    # check if off-targets are in both adata, if it isnt remove that gene from list
    valid_off_targets = [g for g in off_targets if g in adata_visium.var_names]
    removed_genes = [g for g in off_targets if g not in valid_off_targets]
    if removed_genes:
        print(f"Removed off-target genes for {gene}: {removed_genes}")
    if not valid_off_targets:
        continue
    off_targets = valid_off_targets
    # get corr
    corr = stats.pearsonr(np.log(adata_visium[:, off_targets].X.toarray().sum(axis=1).flatten()+1), np.log(adata_xenium[:, gene].X.toarray().flatten()+1))[0]

    pearson_corr_agg.append(corr)

pearson_corr_agg



common_genes = offtarget_genes['Gene']
# remove genes that are not in both datasets
common_genes = [g for g in common_genes if g in adata_visium.var_names and g in adata_xenium.var_names]


# plot correlation
plt.figure(figsize=(10, 5))
plt.bar(common_genes, pearson_corr, color='#0406EA', alpha=0.5, label='Single Genes')
plt.bar(common_genes, pearson_corr_agg, color='#F2080A', alpha=0.5, label='Aggregated Genes')
plt.xticks(rotation=45)
plt.ylabel('Pearson Correlation')
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.legend()
plt.title('Pearson Correlation between Visium Data and Xenium Data')
sns.despine()
plt.show()



############################################################################################################


### plot other genes ###


# log transform the data
sc.pp.log1p(adata_xenium)
sc.pp.log1p(adata_visium)


# New resolution that doesn;t cause issues
new_resolution_xenium = 275

# Rerasterize the patches
aligned_visium_dictionary_rerastered = rerasterize_patches(aligned_visium_dictionary, new_resolution_xenium)
aligned_xenium_dictionary_rerastered = rerasterize_patches(aligned_xenium_dictionary, new_resolution_xenium)


# drfine your output folder
output_base = "/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/gene_plots_notofftarget"
os.makedirs(output_base, exist_ok=True)

# define the gene lists
gene_list_visium = list(adata_visium.var.index)
gene_list_xenium = list(adata_xenium.var.index)


genes2plot = ['TACSTD2']

# Loop over each row of dataframe
for main_gene in genes2plot:
    # skip if gene is not in visium gene list
    if main_gene not in gene_list_visium:
        print(f"Skipping {main_gene} because it is not in the Visium gene list.")
        continue
    
    # Create a folder for this gene inside the output folder.
    gene_folder = os.path.join(output_base, main_gene)
    os.makedirs(gene_folder, exist_ok=True)
    
    # Plot Xenium data for the main gene
    fig = plotRaster(adata_xenium.uns["spatial"], aligned_xenium_dictionary_rerastered, color_by='gene_expression', gene_name=main_gene, if_vis=False)
    # Save the figure in the gene folder
    xenium_filename = os.path.join(gene_folder, f"{main_gene}_xenium.png")
    fig.savefig(xenium_filename, bbox_inches='tight', dpi=150)
    plt.close(fig)
    
    # Plot Visium data for the main gene.
    fig = plotRaster(adata_visium.uns["spatial"], aligned_visium_dictionary_rerastered, color_by='gene_expression', gene_name=main_gene, if_vis=True)
    visium_main_filename = os.path.join(gene_folder, f"{main_gene}_visium.png")
    fig.savefig(visium_main_filename, bbox_inches='tight', dpi=150)
    plt.close(fig)
    