############################################################################################################

# Plot two images side by side with gene expression colored patches


############################################################################################################

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
import pandas as pd
import scanpy as sc
from PIL import Image
import pickle
import scipy
from collections import defaultdict, Counter
from scipy.stats import linregress
import seaborn as sns


# # file name
# file_name = "breastcancer_xenium_sample1_rep2"
# # resolution
# resolution = 250
# # read in the data
# adata_rep2 = sc.read_10x_h5('/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_xenium_sample1_rep2/cell_feature_matrix.h5')

# # Load the full-resolution spatial data
# # cell_centers = pd.read_csv("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_xenium_sample1_rep2/breastcancer_xenium_sample1_rep2_visium_high_res_STalign.csv.gz", index_col=0)
# # cell_centers = pd.read_csv("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/janesick_nature_comms_2023_companion/xenium_cell_centroids_visium_high_res.csv")
# cell_centers = pd.read_csv("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_xenium_sample1_rep2/breastcancer_xenium_sample1_rep2_fullresolution_STalign.csv.gz", index_col=0)
# cell_centers

# # Load the full-resolution image
# Image.MAX_IMAGE_PIXELS = None
# img_name = "Xenium_FFPE_Human_Breast_Cancer_Rep2_he_image"
# img = np.array(Image.open("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/" + file_name + "/" + img_name + ".tif"))
# # img = np.load("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/janesick_nature_comms_2023_companion/visium_high_res_image.npy")
# # plt.imshow(img)

# # add .obs
# adata_rep2.obs = cell_centers
# # add .obsm
# adata_rep2.obsm["spatial"] = adata_rep2.obs[["x_centroid", "y_centroid"]].to_numpy().astype(int)
# # add image
# adata_rep2.uns['spatial'] = img
# # need to add this for subsetting
# adata_rep2.obs.index = adata_rep2.obs.index.astype(str)


# # get rid of genes that aren't in visium
# gene_list = pd.read_csv("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_xenium_sample1_rep1/rastGexp_df.csv", index_col=0)
# gene_list = [gene for gene in gene_list.index if "BLANK" not in gene and "Neg" not in gene and  "antisense" not in gene]
# gene_list = [gene for gene in gene_list if gene not in ['AKR1C1', 'ANGPT2', 'APOBEC3B', 'BTNL9', 'CD8B', 'POLR2J3', 'TPSAB1']]
# # subset the data
# adata_rep2 = adata_rep2[:, gene_list]

# # make an array of the gene expression data
# adata_rep2.X_array = pd.DataFrame(adata_rep2.X.toarray(), index=adata_rep2.obs.index)

# # need to subset bc there are negative values
# adata_rep2 = adata_rep2[adata_rep2.obs["y_centroid"] > 0]
# adata_rep2 = adata_rep2[adata_rep2.obs["x_centroid"] > 0]

# # Extract patches using `extract_patches_from_centers` function
# adata_patches = rasterizeGeneExpression_topatches(img, adata_rep2, patch_size=resolution, aggregation='sum')
# len(adata_patches)

# # scale the data
# scaling_factor = 1
# for i in adata_patches:
#     # aligned_visium_dictionary[i].X_array = sc.pp.log1p(aligned_visium_dictionary[i].X_array * scaling_factor)
#     adata_patches[i].X = sc.pp.log1p(np.round(adata_patches[i].X * scaling_factor))

# # combine the adata patches
# combined_adata = combine_adata_patches(adata_patches, img)

# # Example call to plot patches based on a specific obs column
# plotRaster(img, adata_patches, color_by='total_expression')
# plotRaster(img, adata_patches, color_by='gene_expression', gene_name='LPL')



############################################################################################################



### read in xenium data ###


# combined data
adata_xenium = sc.read_h5ad('/home/caleb/Desktop/improvedgenepred/data/breastcancer_sample1_rep1_aligned_tovisiumimage/xenium_data_full.h5ad')
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

with open('/home/caleb/Desktop/improvedgenepred/data/breastcancer_sample1_rep1_aligned_tovisiumimage/xenium_patches_full.pkl', 'rb') as f:
    aligned_xenium_dictionary = pickle.load(f)

# plt.imshow(adata_visium.uns['spatial'])

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


# plot the total counts of visium vs xenium as a scatter plot
plt.figure(figsize=(10, 5))
plt.scatter(adata_visium_sub.var['total_counts_sub_log'], adata_xenium_sub.var['total_counts_sub_log'])
plt.xlabel("Visium")
plt.ylabel("Xenium")
plt.title("Total counts of genes in visium vs xenium")
plt.xlim([min(adata_visium_sub.var['total_counts_sub_log'])-.5, max(adata_visium_sub.var['total_counts_sub_log'])+.5])
plt.ylim([min(adata_xenium_sub.var['total_counts_sub_log'])-.5, max(adata_xenium_sub.var['total_counts_sub_log'])+.5])
# plot a line
plt.plot([min(adata_visium_sub.var['total_counts_sub_log'])-.5, max(adata_visium_sub.var['total_counts_sub_log'])+.5], [min(adata_visium_sub.var['total_counts_sub_log'])-.5, max(adata_visium_sub.var['total_counts_sub_log'])+.5], color='red')
plt.show()



# plot the total counts of visium vs xenium as a scatter plot with gene names
plt.figure(figsize=(10, 5))
plt.scatter(adata_visium_sub.var['total_counts_sub_log'], adata_xenium_sub.var['total_counts_sub_log'])

# add gene names
for i, gene in enumerate(adata_visium_sub.var_names):
    plt.text(adata_visium_sub.var['total_counts_sub_log'][i], adata_xenium_sub.var['total_counts_sub_log'][i], gene, fontsize=8)

plt.xlabel("Visium")
plt.ylabel("Xenium")
plt.title("Total counts of genes in visium vs xenium")
plt.xlim([min(adata_visium_sub.var['total_counts_sub_log'])-.5, max(adata_visium_sub.var['total_counts_sub_log'])+.5])
plt.ylim([min(adata_xenium_sub.var['total_counts_sub_log'])-.5, max(adata_xenium_sub.var['total_counts_sub_log'])+.5])
plt.show()






# read in text file
# gene_summary = pd.read_csv("/home/caleb/Desktop/off-target-probe-checker/hisat2/gene_summary_transcriptome.txt", sep="\t")
# gene_summary = pd.read_csv("/home/caleb/Desktop/off-target-probe-checker/caleb/bowtie2/gene_summary_transcriptome.txt", sep="\t")
# gene_summary = pd.read_csv("/home/caleb/Desktop/off-target-probe-checker/hisat2_genome/gene_summary_genome.txt", sep="\t")
# gene_summary = pd.read_csv("/home/caleb/Desktop/off-target-probe-checker/hayden/probe2targets.named.tsv", sep="\t")

gene_summary = pd.read_csv("/home/caleb/Desktop/off-target-probe-checker/caleb/bowtie2/probe2targets.named.5bp.tsv", sep="\t", header=None)
gene_summary.columns = ["probe_id", "n_targets", "target_names"]

# get gene names
gene_summary["probe_gene"] = gene_summary["probe_id"].str.split("|").str[1]
# change column names
gene_summary = gene_summary.rename(columns={"target_names": "probe_hit"})


# create a defaultdict that will hold counter for each probe_gene
gene_dict = defaultdict(Counter)

# iterate over each row in the gene_summary dataframe
for idx, row in gene_summary.iterrows():
    # key is target gene
    key = row['probe_gene']
    
    # get the probe_hit string and remove square brackets
    hits_str = row['probe_hit'].strip()  # remove leading/trailing whitespace
    if hits_str.startswith('[') and hits_str.endswith(']'):
        hits_str = hits_str[1:-1]
    
    # split the string on commas and strip any whitespace.
    hits = [hit.strip() for hit in hits_str.split(',') if hit.strip()]
    
    # add to counter for this probe_gene with each hit
    for hit in hits:
        gene_dict[key][hit] += 1

# convert the defaultdict to a normal dict
gene_dict = dict(gene_dict)
print(gene_dict)

# # check
# gene_dict['ACTG2']
# gene_dict['CEACAM8']

# # genes in paper
# gene_dict['MS4A1']
# gene_dict['ERBB2']
# gene_dict['KRT14'] # only one with off-targets
# gene_dict['SFRP4']
# gene_dict['APOC1']
# gene_dict['TACSTD2']

# make dict a df
gene_summary = pd.DataFrame(gene_dict.items(), columns=["Gene", "Matches"])
print(gene_summary)


# function to compute the difference: (total counts) - (self-match count)
def compute_diff(row):
    matches = row['Matches']
    total = sum(matches.values())
    # Get the count for the self-match if present; otherwise use 0
    self_count = matches.get(row['Gene'], 0)
    return total - self_count

# compute difference
gene_summary['diff'] = gene_summary.apply(compute_diff, axis=1)


gene_summary[gene_summary['Gene'] == "AQP4"]



# count total number of rows where diff != 0
len(gene_summary[gene_summary['diff'] != 0])

# Plot a histogram of the differences
plt.figure(figsize=(8,6))
# Use integer bins from 0 to the maximum diff plus one
bins = range(0, max(gene_summary['diff']) + 2)
plt.hist(gene_summary['diff'], bins=bins, edgecolor='black')
plt.xlabel('Difference')
plt.ylabel('Frequency')
plt.title('Histogram of Differences in Match Counts')
plt.show()






gene_summary_tmp = gene_summary.copy()
gene_summary_tmp = gene_summary_tmp.set_index('Gene')
# merge df
merged_df = adata_visium_sub.var.join(gene_summary_tmp, how="left")
# make adata.var the merged
adata_visium_sub.var = merged_df



import seaborn as sns


# plot the scatterplot but color by binary diff (0 or not 0)
colors = ['gray' if pd.isna(diff) else ('red' if diff != 0 else 'blue') 
          for diff in adata_visium_sub.var["diff"]]
plt.figure(figsize=(8, 8))
plt.scatter(adata_visium_sub.var['total_counts_sub_log'], adata_xenium_sub.var['total_counts_sub_log'], c=colors)
plt.xlabel("Log Total Counts Visium", fontsize=16)
plt.ylabel("Log Total Counts Xenium", fontsize=16)
plt.title("Total counts of genes in visium vs xenium")
plt.xlim([min(adata_visium_sub.var['total_counts_sub_log'])-.5, max(adata_visium_sub.var['total_counts_sub_log'])+.5])
plt.ylim([min(adata_xenium_sub.var['total_counts_sub_log'])-.5, max(adata_xenium_sub.var['total_counts_sub_log'])+.5])
# plot a line
plt.plot([min(adata_visium_sub.var['total_counts_sub_log'])-.5, max(adata_visium_sub.var['total_counts_sub_log'])+.5], [min(adata_visium_sub.var['total_counts_sub_log'])-.5, max(adata_visium_sub.var['total_counts_sub_log'])+.5], color='red', linewidth=2, linestyle='--')
# plot a line for 10X = Y
x_vals = np.linspace(min(adata_visium_sub.var['total_counts_sub_log'])-.5, max(adata_visium_sub.var['total_counts_sub_log'])+.5, 100)
y_vals = np.log(10) + x_vals
plt.plot(x_vals, y_vals, color='green', linewidth=2)
# increase font size
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# Custom legend
import matplotlib.patches as mpatches
red_patch = mpatches.Patch(color='red', label='Off Target Probes')
blue_patch = mpatches.Patch(color='blue', label='No Off Target Probes')
gray_patch = mpatches.Patch(color='gray', label='Custom Probes')
plt.legend(handles=[red_patch, blue_patch, gray_patch])
sns.despine()
plt.show()




############################################################################################################



# make new df that has just genes with off-targets
gene_summary_filtered = gene_summary[gene_summary["diff"] != 0]

# not get only matches that are within adata_visium.var.index


# new column with only genes present in allowed_genes
gene_summary_filtered['Filtered_Matches'] = gene_summary_filtered['Matches'].apply(lambda d: {gene: count for gene, count in d.items() if gene in full_gene_list})

# reset index
gene_summary_filtered = gene_summary_filtered.reset_index(drop=True)

# save to csv
# gene_summary_filtered.to_csv("/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/gene_summary_filtered.csv", index=False)



# get rid or rows with empty lists in Misamatches column
gene_summary_filtered = gene_summary_filtered[gene_summary_filtered["Filtered_Matches"].apply(lambda x: len(x) > 0)]

# get rid of rows where Gene column is the only gene in the Filtered_Matches column
gene_summary_filtered = gene_summary_filtered[
    ~gene_summary_filtered.apply(
        lambda row: len(row['Filtered_Matches']) == 1 and list(row['Filtered_Matches'].keys())[0] == row['Gene'],
        axis=1
    )
]

# reset index
gene_summary_filtered = gene_summary_filtered.reset_index(drop=True)
gene_summary_filtered


def rerasterize_patches(adata_patches, patch_size):
    """
    Adjusts the coordinates of patches to align them on a uniform grid.
    If two patches are merged into one, their expression values are averaged.
    
    Parameters:
    - adata_patches: Dictionary of AnnData objects representing the patches.
    - patch_size: The size of each patch.
    
    Returns:
    - new_adata_patches: Dictionary of AnnData objects with adjusted coordinates.
    """
    new_adata_patches = {}

    # Define grid step size
    step_size = patch_size

    # Create a set of unique x and y coordinates and sort them
    x_coords = sorted({adata_patch.uns['patch_coords'][0] for adata_patch in adata_patches.values()})
    y_coords = sorted({adata_patch.uns['patch_coords'][2] for adata_patch in adata_patches.values()})

    # Snap x and y coordinates to a uniform grid
    x_grid = np.arange(min(x_coords), max(x_coords) + step_size, step_size)
    y_grid = np.arange(min(y_coords), max(y_coords) + step_size, step_size)

    # A dictionary to keep track of merged patches and counts for averaging
    merged_patches = {}
    patch_counts = {}

    for idx, adata_patch in adata_patches.items():
        # print(idx)
        x_start, x_end, y_start, y_end = adata_patch.uns['patch_coords']
        
        # Find the nearest x and y positions on the grid
        x_center = (x_start + x_end) / 2
        y_center = (y_start + y_end) / 2

        new_x_start = x_grid[np.abs(x_grid - x_center).argmin()] - step_size / 2
        new_y_start = y_grid[np.abs(y_grid - y_center).argmin()] - step_size / 2
        new_x_end = new_x_start + patch_size
        new_y_end = new_y_start + patch_size

        # Create a unique key for the grid position
        grid_key = (new_x_start, new_y_start)

        if grid_key in merged_patches:
            # Merge the expression values by taking the mean
            existing_patch = merged_patches[grid_key]
            existing_patch.X = (existing_patch.X * patch_counts[grid_key] + adata_patch.X) / (patch_counts[grid_key] + 1)
            patch_counts[grid_key] += 1
        else:
            # Add the patch to the merged dictionary
            new_adata_patch = adata_patch.copy()
            new_adata_patch.uns['patch_coords'] = [new_x_start, new_x_end, new_y_start, new_y_end]
            merged_patches[grid_key] = new_adata_patch
            patch_counts[grid_key] = 1

    # Convert the merged patches to the final dictionary
    for idx, (grid_key, adata_patch) in enumerate(merged_patches.items()):
        new_adata_patches[idx] = adata_patch

    return new_adata_patches


# log transform the data
sc.pp.log1p(adata_xenium)
sc.pp.log1p(adata_visium)


# New resolution that doesn;t cause issues
new_resolution_xenium = 275

aligned_visium_dictionary_rerastered = rerasterize_patches(aligned_visium_dictionary, new_resolution_xenium)
aligned_xenium_dictionary_rerastered = rerasterize_patches(aligned_xenium_dictionary, new_resolution_xenium)




import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def plotRaster(image, adata_patches, color_by='gene_expression', gene_name=None, if_vis=True):
    """
    Plots patches on the original image, colored by either gene expression or total expression.
    
    Parameters:
    - image: The original image array.
    - adata_patches: Dictionary of AnnData objects representing the patches.
    - color_by: How to color the patches ('gene_expression' or 'total_expression').
    - gene_name: The name of the gene to use if color_by is 'gene_expression'.
    - if_vis: Boolean flag to set the title as "Visium" (True) or "Xenium" (False).
    
    Returns:
    - fig: Matplotlib Figure object.
    """
    # Check inputs
    if color_by == 'gene_expression' and gene_name is None:
        raise ValueError("You must specify a gene_name when color_by='gene_expression'.")
    
    # Collect all values for normalization
    values = []
    for adata_patch in adata_patches.values():
        if color_by == 'gene_expression':
            # Sum the expression for the specified gene
            expression = adata_patch.X[:, adata_patch.var_names.get_loc(gene_name)].sum()
            values.append(expression)
        elif color_by == 'total_expression':
            total_expression = adata_patch.X.sum()
            values.append(total_expression)
    
    values = np.array(values)
    min_value, max_value = values.min(), values.max()
    
    # Plot the original image
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(image)
    
    # Plot each patch with the appropriate color
    for adata_patch in adata_patches.values():
        # Expecting patch_coords as (x_start, x_end, y_start, y_end)
        x_start, x_end, y_start, y_end = adata_patch.uns['patch_coords']
        
        if color_by == 'gene_expression':
            expression = adata_patch.X[:, adata_patch.var_names.get_loc(gene_name)].sum()
            normalized_value = (expression - min_value) / (max_value - min_value) if max_value > min_value else 0
            color = plt.cm.viridis(normalized_value)
        elif color_by == 'total_expression':
            total_expression = adata_patch.X.sum()
            normalized_value = (total_expression - min_value) / (max_value - min_value) if max_value > min_value else 0
            color = plt.cm.viridis(normalized_value)
        
        # Draw a rectangle (square patch)
        rect = patches.Rectangle((x_start, y_start), x_end - x_start, y_end - y_start,
                                 linewidth=1, edgecolor='none', facecolor=color, alpha=1)
        ax.add_patch(rect)
    
    # Create a color bar
    norm = plt.Normalize(min_value, max_value)
    sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation='vertical', fraction=0.03, pad=0.04)
    cbar.set_label(f'{gene_name} Expression' if color_by == 'gene_expression' else "Total Expression")
    
    # Set the plot title
    title_prefix = "Visium" if if_vis else "Xenium"
    ax.set_title(f"{title_prefix} Expression of {gene_name}", fontsize=16)
    
    # Remove axis ticks
    ax.set_xticks([])
    ax.set_yticks([])
    
    return fig






import os
import matplotlib.pyplot as plt
import scipy.sparse as sp


# Turn off interactive mode so figures can be saved without displaying.
plt.ioff()

# Define your output folder for plots.
output_base = "/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/gene_plots"
os.makedirs(output_base, exist_ok=True)

# For convenience, define the gene lists.
gene_list_visium = list(adata_visium.var.index)
gene_list_xenium = list(adata_xenium.var.index)


# Loop over each row of your dataframe that has off-targets.
# gene_summary_filtered should have at least these columns: "Gene" and "Filtered_Matches"
# for idx, row in gene_summary_filtered.iterrows():
for idx, row in gene_summary_filtered.iterrows():
    main_gene = row["Gene"]
    if main_gene not in gene_list_visium:
        print(f"Skipping {main_gene} because it is not in the Visium gene list.")
        continue
    filtered_matches = row["Filtered_Matches"]  # This is a dictionary (e.g., {'HOXD9': 7, 'HOXB9': 1})
    
    # Create a folder for this gene inside the output folder.
    gene_folder = os.path.join(output_base, main_gene)
    os.makedirs(gene_folder, exist_ok=True)
    
    # -------------------------------
    # 1. Plot Xenium data for the main gene.
    # -------------------------------
    # (Assuming plot_xenium_with_centers returns a matplotlib Figure object.)
    # fig = plot_xenium_with_centers(adata_xenium, gene_list_xenium, main_gene, adata_xenium.obsm['spatial'], patch_size=300, if_vis=False)
    fig = plotRaster(adata_xenium.uns["spatial"], aligned_xenium_dictionary_rerastered, color_by='gene_expression', gene_name=main_gene, if_vis=False)
    # Save the figure in the gene folder.
    xenium_filename = os.path.join(gene_folder, f"{main_gene}_xenium.png")
    fig.savefig(xenium_filename, bbox_inches='tight', dpi=150)
    plt.close(fig)
    
    # -------------------------------
    # 2. Plot Visium data for the main gene.
    # -------------------------------
    # fig = plot_xenium_with_centers(adata_visium, gene_list_visium, main_gene, adata_visium.obsm['spatial'], patch_size=300, if_vis=True)
    fig = plotRaster(adata_visium.uns["spatial"], aligned_visium_dictionary_rerastered, color_by='gene_expression', gene_name=main_gene, if_vis=True)

    visium_main_filename = os.path.join(gene_folder, f"{main_gene}_visium.png")
    fig.savefig(visium_main_filename, bbox_inches='tight', dpi=150)
    plt.close(fig)
    
    # -------------------------------
    # 3. For each off-target gene (i.e. any key in Filtered_Matches that is not main_gene)
    # Plot Visium data for the off-target.
    # -------------------------------
    for off_target in filtered_matches.keys():
        if off_target != main_gene:
            # fig = plot_xenium_with_centers(adata_visium, gene_list_visium, off_target, adata_visium.obsm['spatial'], patch_size=300, if_vis=True)
            fig = plotRaster(adata_visium.uns["spatial"], aligned_visium_dictionary_rerastered, color_by='gene_expression', gene_name=off_target, if_vis=True)
            off_target_filename = os.path.join(gene_folder, f"{main_gene}_visium_offtarget_{off_target}.png")
            fig.savefig(off_target_filename, bbox_inches='tight', dpi=150)
            plt.close(fig)

    # # -------------------------------------------------------------------------
    # # (4) Plot an aggregated gene in Visium (sum of all off-targets + main_gene)
    # # -------------------------------------------------------------------------
    # # Build a list of all relevant genes (main + off-targets)
    # all_genes = list(filtered_matches.keys())  # includes main_gene plus others
    # # Filter to only those that exist in adata_visium
    # valid_genes = [g for g in all_genes if g in adata_visium.var_names]
    # if not valid_genes:
    #     # No valid genes to sum => skip
    #     continue

    # # aggregator gene name
    # agg_gene_name = main_gene + "_agg"

    # # list of relevant genes (main + off-targets)
    # all_genes = list(filtered_matches.keys())  
    # if main_gene not in all_genes:
    #     all_genes.append(main_gene)

    # # filter to only those present in var_names
    # valid_genes = [g for g in all_genes if g in adata_visium.var_names]
    # if not valid_genes:
    #     continue  # skip if none present

    # # subset the data to valid genes
    # subset_adata = adata_visium[:, valid_genes]

    # # (1) Convert from log scale back to raw counts.
    # #     Assuming adata_visium.X = np.log1p(raw_counts)
    # if sp.issparse(subset_adata.X):
    #     log_vals = subset_adata.X.toarray()  # dense for demonstration
    # else:
    #     log_vals = subset_adata.X

    # raw_vals = np.expm1(log_vals)  # exp(x) - 1 => revert from log1p

    # # (2) Sum across genes (columns). shape => (#spots, #genes) => sum across axis=1 to get (#spots, )
    # summed_raw = raw_vals.sum(axis=1)  

    # # (3) Re-log the sum => aggregator in log1p scale
    # summed_log = np.log1p(summed_raw)

    # # (4) Insert aggregator gene into adata_visium if not already there
    # if agg_gene_name not in adata_visium.var_names:
    #     # add row in var
    #     from math import nan
    #     new_row = [nan] * len(adata_visium.var.columns)
    #     adata_visium.var.loc[agg_gene_name] = new_row

    #     # now add a new column to .X
    #     if sp.issparse(adata_visium.X):
    #         # existing expression is sparse => hstack
    #         new_col = sp.csr_matrix(summed_log.reshape(-1, 1))
    #         adata_visium._X = sp.hstack([adata_visium.X, new_col])
    #         adata_visium._X = sp.csr_matrix(adata_visium._X)
    #     else:
    #         new_col = summed_log.reshape(-1, 1)
    #         adata_visium._X = np.hstack([adata_visium.X, new_col])
        
    #     # update var_names
    #     adata_visium.var_names = list(adata_visium.var_names) + [agg_gene_name]

    # # (5) Now plot the aggregator with your existing function, e.g. plotRaster()
    # fig = plotRaster(
    #     adata_visium.uns["spatial"],
    #     aligned_visium_dictionary_rerastered,
    #     color_by='gene_expression',
    #     gene_name=agg_gene_name,  # aggregator gene
    #     if_vis=True
    # )

    # # NOTE: need to get new patches for the aggregator gene

    # agg_filename = os.path.join(gene_folder, f"{main_gene}_visium_aggregated.png")
    # fig.savefig(agg_filename, bbox_inches='tight', dpi=150)
    # plt.close(fig)

print("Plots saved for all genes.")






############################################################################################################





### read in xenium data ###


# combined data
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


# combine gene expression
# read in genes of interest
offtarget_genes = pd.read_csv("/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/gene_summary_filtered_5bp.csv")
offtarget_genes

offtarget_genes[['Gene', 'Filtered_Matches']]

import scipy.stats as stats
import ast

# find correlation between xenium and visium
# calculate pearson correlation between xenium and visium per gene
pearson_corr = []
for gene in offtarget_genes['Gene']:
    # check if gene is in both adata
    if gene not in adata_visium.var_names or gene not in adata_xenium.var_names:
        print(f"Skipping {gene} because it is not in both datasets.")
        continue
    corr = stats.pearsonr(np.log10(adata_visium[:, gene].X.toarray().flatten()+1), np.log10(adata_xenium[:, gene].X.toarray().flatten()+1))[0]
    pearson_corr.append(corr)


# now do the same but combine the gene expression of off-targets
pearson_corr_agg = []
for i in range(0, len(offtarget_genes)):
    print(i)
    gene = offtarget_genes['Gene'][i]
    # off-targets
    off_targets = list(ast.literal_eval(offtarget_genes.loc[i, 'Filtered_Matches']).keys())
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
    corr = stats.pearsonr(np.log10(adata_visium[:, off_targets].X.toarray().sum(axis=1).flatten()+1), np.log10(adata_xenium[:, gene].X.toarray().flatten()+1))[0]

    pearson_corr_agg.append(corr)

pearson_corr_agg



common_genes = offtarget_genes['Gene']
# remove genes that are not in both datasets
common_genes = [g for g in common_genes if g in adata_visium.var_names and g in adata_xenium.var_names]

# plot correlation
plt.figure(figsize=(10, 5))
plt.bar(common_genes, pearson_corr, color='blue', alpha=0.5, label='Single Genes')
plt.bar(common_genes, pearson_corr_agg, color='red', alpha=0.5, label='Aggregated Genes')
plt.xticks(rotation=90)
plt.ylabel('Pearson Correlation')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.show()





# import scipy.sparse as sp

# import numpy as np
# import pandas as pd
# import scipy.sparse as sp

# for index, row in offtarget_genes.iterrows():
#     # Use the main gene name + "_agg" as the new aggregated gene
#     gene_name_agg = row['Gene'] + '_agg'
#     matched_genes = row['Filtered_Matches']
    
#     # 1. Determine which of the matched genes are actually present in the AnnData
#     valid_genes = [g for g in matched_genes.keys() if g in adata_visium.var_names]
#     if not valid_genes:
#         continue  # Skip if no valid genes are found in var_names
    
#     # 2. Sum the expression of the valid genes across all cells (obs)
#     if sp.issparse(adata_visium.X):
#         # Convert the sum to a dense 1D array
#         aggregated_expr = np.array(adata_visium[:, valid_genes].X.sum(axis=1)).ravel()
#     else:
#         aggregated_expr = adata_visium[:, valid_genes].X.sum(axis=1)
    

#     # (a) Add a new row in var
#     adata_visium.var.loc[gene_name_agg] = [None] * len(adata_visium.var.columns)
    
#     # (b) Expand the X matrix to add one more column
#     if sp.issparse(adata_visium.X):
#         # Add a new sparse column
#         new_col = sp.csr_matrix(aggregated_expr.reshape(-1, 1))
#         adata_visium._X = sp.hstack([adata_visium.X, new_col])
#         # Convert back to CSR (optional, but often convenient)
#         adata_visium._X = sp.csr_matrix(adata_visium._X)
#     else:
#         # Dense case
#         new_col = aggregated_expr.reshape(-1, 1)
#         adata_visium._X = np.hstack([adata_visium.X, new_col])
    



# # calculate pearson correlation
# from scipy.stats import pearsonr
# # calcilar spearman correlation
# from scipy.stats import spearmanr

# # Calculate the Pearson correlation between the aggregated expression of the main gene and the off-targets.

# # g = "TUBB2B"

# # pearsonr(adata_visium[:, g].X.toarray().flatten(), adata_xenium[:, g].X.toarray().flatten())
# # pearsonr(adata_visium[:, g + '_agg'].X.toarray().flatten(), adata_xenium[:, g].X.toarray().flatten())


# pearson_diff = []
# genes_used = []
# # iterate through each gene in gene_summary_filtered and calculate the pearson correlation
# for index, row in gene_summary_filtered.iterrows():
#     gene_name = row['Gene']

#     if gene_name not in adata_visium.var_names:
#         continue
#     # NOTE: getting rid of a few genes here...
#     genes_used.append(gene_name)

#     # Calculate the Pearson correlation between the aggregated expression of the main gene and the off-targets.
#     corr_visium = pearsonr(adata_visium[:, gene_name].X.toarray().flatten(), adata_xenium[:, gene_name].X.toarray().flatten())[0]
#     corr_visium_agg = pearsonr(adata_visium[:, gene_name + '_agg'].X.toarray().flatten(), adata_xenium[:, gene_name].X.toarray().flatten())[0]

#     pearson_diff.append(corr_visium_agg - corr_visium)

# pearson_diff = pd.DataFrame({'Gene': genes_used, 'Pearson Difference': pearson_diff})

# # merge with gene_summary_filtered with the pearson_diff
# gene_summary_filtered = gene_summary_filtered.merge(pearson_diff, on='Gene', how='left')



# sns.scatterplot(data=gene_summary_filtered, x='diff', y='Pearson Difference')

# # Add gene names to each point
# for i, gene in enumerate(gene_summary_filtered['Gene']):
#     plt.text(gene_summary_filtered['diff'][i], gene_summary_filtered['Pearson Difference'][i], gene, fontsize=8)

# # plot regression line
# sns.regplot(data=gene_summary_filtered, x='diff', y='Pearson Difference', scatter=False)
# plt.show()




# plt.hist(gene_summary_filtered['Pearson Difference'])
# plt.show()



# NOTE: maybe a comparison of total expression per patch?
# NOTE: look into scRNA-seq data to compare as well
# NOTE: look at off target probes in visium!
# NOTE: make alignment more senstive














# # Create the dictionary
# gene_dict = {}

# for i in range(0, len(gene_summary['Mismatches'])):
#     gene_summary['Mismatches'][i]
#     if pd.isna(gene_summary['Mismatches'][i]):
#         gene_dict[gene_summary['Gene'][i]] = []
#     else:
#         gene_list = list(set(gene_summary['Mismatches'][i].split(",")))  # Split by ',' and get unique genes
#         gene_dict[gene_summary['Gene'][i]] = gene_list

# # Print 
# print(gene_dict)

# # get gene list
# # read data
# adata_visium_full = sc.read_visium("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_visium/")
# gene_list = list(adata_visium_full.var_names)

# gene_dict_filtered = {}

# for i in range(len(gene_summary['Mismatches'])):
#     gene_name = gene_summary['Gene'][i]  # Key
#     mismatches = gene_summary['Mismatches'][i]  # Value

#     if pd.isna(mismatches):
#         gene_dict_filtered[gene_name] = []
#     else:
#         # Split, get unique genes, and filter by `gene_list`
#         gene_set = set(mismatches.split(","))
#         filtered_genes = [gene for gene in gene_set if gene in gene_list]

#         gene_dict_filtered[gene_name] = filtered_genes  # Store filtered mismatches

# # Print or use the dictionary
# print(gene_dict_filtered)

# # make dict a df
# gene_summary_filtered = pd.DataFrame(gene_dict_filtered.items(), columns=["Gene", "Mismatches"])
# # get rid or rows with empty lists in Misamatches column
# gene_summary_filtered = gene_summary_filtered[gene_summary_filtered["Mismatches"].apply(lambda x: len(x) > 0)]
# gene_summary_filtered
# len(gene_summary_filtered)


# sort genes
# gene_summary = gene_summary.sort_values('Gene', ascending=True)
# # get difference between columns 1 and 2
# gene_summary["diff"] = gene_summary["Total_Count"] - gene_summary["Matched_Count"]
# # set gene as index
# gene_summary = gene_summary.set_index('Gene')
# # merge df
# merged_df = adata_visium.var.join(gene_summary, how="left")

# # make adata.var the merged
# adata_visium.var = merged_df



# # plot the scatterplot but color by diff
# plt.figure(figsize=(10, 5))
# plt.scatter(np.log(adata_visium.var['total_counts_sub']), np.log(adata_xenium.var['total_counts_sub']), c=np.log(adata_visium.var["diff"]+1))
# plt.xlabel("Visium")
# plt.ylabel("Xenium")
# plt.title("Total counts of genes in visium vs xenium")
# plt.xlim([min(np.log(adata_visium.var['total_counts_sub']))-.5, max(np.log(adata_visium.var['total_counts_sub']))+.5])
# plt.ylim([min(np.log(adata_xenium.var['total_counts_sub']))-.5, max(np.log(adata_xenium.var['total_counts_sub']))+.5])
# # plot a line
# plt.plot([min(np.log(adata_visium.var['total_counts_sub']))-.5, max(np.log(adata_visium.var['total_counts_sub']))+.5], [min(np.log(adata_visium.var['total_counts_sub']))-.5, max(np.log(adata_visium.var['total_counts_sub']))+.5], color='red')
# plt.colorbar()
# plt.show()




# # plot the scatterplot but color by binary diff (0 or not 0)
# colors = ['gray' if pd.isna(diff) else ('red' if diff != 0 else 'blue') 
#           for diff in adata_visium.var["diff"]]
# plt.figure(figsize=(10, 5))
# plt.scatter(np.log(adata_visium.var['total_counts_sub']), np.log(adata_xenium.var['total_counts_sub']), c=colors)
# plt.xlabel("Visium")
# plt.ylabel("Xenium")
# plt.title("Total counts of genes in visium vs xenium")
# plt.xlim([min(np.log(adata_visium.var['total_counts_sub']))-.5, max(np.log(adata_visium.var['total_counts_sub']))+.5])
# plt.ylim([min(np.log(adata_xenium.var['total_counts_sub']))-.5, max(np.log(adata_xenium.var['total_counts_sub']))+.5])
# # plot a line
# plt.plot([min(np.log(adata_visium.var['total_counts_sub']))-.5, max(np.log(adata_visium.var['total_counts_sub']))+.5], [min(np.log(adata_visium.var['total_counts_sub']))-.5, max(np.log(adata_visium.var['total_counts_sub']))+.5], color='red')
# plt.show()



# # plot the scatterplot but color by total alignments
# plt.figure(figsize=(10, 5))
# plt.scatter(np.log(adata_visium.var['total_counts_sub']), np.log(adata_xenium.var['total_counts_sub']), c=np.log(adata_visium.var['Total_Count']))
# plt.xlabel("Visium")
# plt.ylabel("Xenium")
# plt.title("Total counts of genes in visium vs xenium")
# plt.xlim([min(np.log(adata_visium.var['total_counts_sub']))-.5, max(np.log(adata_visium.var['total_counts_sub']))+.5])
# plt.ylim([min(np.log(adata_xenium.var['total_counts_sub']))-.5, max(np.log(adata_xenium.var['total_counts_sub']))+.5])
# # plot a line
# plt.plot([min(np.log(adata_visium.var['total_counts_sub']))-.5, max(np.log(adata_visium.var['total_counts_sub']))+.5], [min(np.log(adata_visium.var['total_counts_sub']))-.5, max(np.log(adata_visium.var['total_counts_sub']))+.5], color='red')
# plt.colorbar()
# plt.show()
