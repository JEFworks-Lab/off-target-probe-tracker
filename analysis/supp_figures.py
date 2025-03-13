

# Import libraries
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
import seaborn as sns
import os
from PIL import Image
import matplotlib.patches as patches
import sys
sys.path.append('..')
from utils import *
import scipy
import pickle
import cv2



############################################################################################################



# should be the name of image data in adata
tissue_section = "CytAssist_FFPE_Human_Breast_Cancer"

# file path where outs data is located
file_path = "/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_visium/"

# read in svg results
gene_list = pd.read_csv("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_xenium_sample1_rep1/rastGexp_df.csv", index_col=0)
gene_list = [gene for gene in gene_list.index if "BLANK" not in gene and "Neg" not in gene and  "antisense" not in gene]
# these were not in the data
gene_list = [gene for gene in gene_list if gene not in ['AKR1C1', 'ANGPT2', 'APOBEC3B', 'BTNL9', 'CD8B', 'POLR2J3', 'TPSAB1']]
len(gene_list)

### Read in adata ###

# read data
adata_visium = sc.read_visium(file_path)
# make unique
adata_visium.var_names_make_unique()
# get mitochondrial gene expression info
adata_visium.var["mt"] = adata_visium.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata_visium, qc_vars=["mt"], inplace=True)

# make spatial position str to integer
# https://discourse.scverse.org/t/data-fomr-new-spatial-transcriptomics-from-10x/1107/6
adata_visium.obsm['spatial'] = adata_visium.obsm['spatial'].astype(int)


# # get new cell centers for high rez image
# READ IN ALIGNED DATA
aligned_visium_points = np.load("/home/caleb/Desktop/improvedgenepred/data/breastcancer_xenium_sample1_rep1/aligned_visium_points_to_xenium_image.npy").astype(int)
adata_visium.obsm['spatial'] = aligned_visium_points

# subet gene list
adata_visium = adata_visium[:, gene_list]


### read in xenium data ### 


# file name
file_name = "breastcancer_xenium_sample1_rep1"
# resolution
resolution = 12
# read in the data
adata_xenium = sc.read_10x_h5('/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_xenium_sample1_rep1/cell_feature_matrix.h5')

# Load the full-resolution spatial data
cell_centers = pd.read_csv(f"/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/{file_name}/{file_name}_fullresolution_STalign.csv.gz", index_col=0)
# cell_centers = pd.read_csv("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/janesick_nature_comms_2023_companion/xenium_cell_centroids_visium_high_res.csv")
# cell_centers = pd.read_csv("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_visium/scaled_spots_for_xenium_image.csv", index_col=0)
# cell_centers.columns = ["x_centroid", "y_centroid"]

# Load the full-resolution image
Image.MAX_IMAGE_PIXELS = None
img_name = "Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image"
img = np.array(Image.open("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/" + file_name + "/" + img_name + ".tif"))
# img = np.load("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/janesick_nature_comms_2023_companion/visium_high_res_image.npy")
# plt.imshow(img)

# add .obs
adata_xenium.obs = cell_centers
# add .obsm
adata_xenium.obsm["spatial"] = adata_xenium.obs[["x_centroid", "y_centroid"]].to_numpy().astype(int)
# add image
adata_xenium.uns['spatial'] = img
# need to add this for subsetting
adata_xenium.obs.index = adata_xenium.obs.index.astype(str)

# subset the data
adata_xenium = adata_xenium[:, gene_list]

# make an array of the gene expression data
adata_xenium.X_array = pd.DataFrame(adata_xenium.X.toarray(), index=adata_xenium.obs.index)



############################################################################################################


### Supplemental Figure 1 ###


# plot visium spots and xenium spots on xenium image


import shapely
from shapely.geometry import MultiPoint
import matplotlib.patches as patches


Image.MAX_IMAGE_PIXELS = None
img_name = "Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image"
img = np.array(Image.open("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/" + file_name + "/" + img_name + ".tif"))

fig, ax = plt.subplots()

# Show your background image
ax.imshow(img)

# Plot Xenium points
x_xenium = adata_xenium.obsm["spatial"][:,0]
y_xenium = adata_xenium.obsm["spatial"][:,1]
ax.scatter(x_xenium, y_xenium, s=1, c="yellow", marker=".", edgecolor='none')

# Use Shapely to get a rotated bounding rectangle for Xenium points
points_xenium = MultiPoint(list(zip(x_xenium, y_xenium)))
rot_rect_xenium = points_xenium.minimum_rotated_rectangle

# Convert Shapely polygon to a Matplotlib patch
x_coords, y_coords = rot_rect_xenium.exterior.xy
polygon_xenium = patches.Polygon(
    xy=list(zip(x_coords, y_coords)),
    fill=False, edgecolor='C0', linewidth=2
)
ax.add_patch(polygon_xenium)

# Plot Visium points
x_visium = adata_visium.obsm["spatial"][:,0]
y_visium = adata_visium.obsm["spatial"][:,1]
ax.scatter(x_visium, y_visium, edgecolor="black", facecolors='none',
           marker="o", linewidths=1, s=10)

# Use Shapely to get a rotated bounding rectangle for Visium points
points_visium = MultiPoint(list(zip(x_visium, y_visium)))
rot_rect_visium = points_visium.minimum_rotated_rectangle

# Convert Shapely polygon to a Matplotlib patch
x_coords_v, y_coords_v = rot_rect_visium.exterior.xy
polygon_visium = patches.Polygon(
    xy=list(zip(x_coords_v, y_coords_v)),
    fill=False, edgecolor='C1', linewidth=2
)
ax.add_patch(polygon_visium)

ax.axis("off")
# plt.show()
plt.savefig("/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/final_figures_supplementary/supp_figure1/supp_fig1a1.png", dpi=300, bbox_inches='tight')
plt.close()

# read in visium image
img = np.load("/home/caleb/Desktop/improvedgenepred/data/breastcancer_visium/spatial/tissue_alignedtoxenium.npy")


fig, ax = plt.subplots()

# Show your background image
ax.imshow(img)

# Plot Xenium points
x_xenium = adata_xenium.obsm["spatial"][:,0]
y_xenium = adata_xenium.obsm["spatial"][:,1]
ax.scatter(x_xenium, y_xenium, s=1, c="yellow", marker=".", edgecolor='none')

# Use Shapely to get a rotated bounding rectangle for Xenium points
points_xenium = MultiPoint(list(zip(x_xenium, y_xenium)))
rot_rect_xenium = points_xenium.minimum_rotated_rectangle

# Convert Shapely polygon to a Matplotlib patch
x_coords, y_coords = rot_rect_xenium.exterior.xy
polygon_xenium = patches.Polygon(
    xy=list(zip(x_coords, y_coords)),
    fill=False, edgecolor='C0', linewidth=2
)
ax.add_patch(polygon_xenium)

# Plot Visium points
x_visium = adata_visium.obsm["spatial"][:,0]
y_visium = adata_visium.obsm["spatial"][:,1]
ax.scatter(x_visium, y_visium, edgecolor="black", facecolors='none',
           marker="o", linewidths=1, s=10)

# Use Shapely to get a rotated bounding rectangle for Visium points
points_visium = MultiPoint(list(zip(x_visium, y_visium)))
rot_rect_visium = points_visium.minimum_rotated_rectangle

# Convert Shapely polygon to a Matplotlib patch
x_coords_v, y_coords_v = rot_rect_visium.exterior.xy
polygon_visium = patches.Polygon(
    xy=list(zip(x_coords_v, y_coords_v)),
    fill=False, edgecolor='C1', linewidth=2
)
ax.add_patch(polygon_visium)

ax.axis("off")
# plt.show()
plt.savefig("/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/final_figures_supplementary/supp_figure1/supp_fig1a2.png", dpi=300, bbox_inches='tight')
plt.close()




### read in aligned data ###

# combined data
adata_xenium = sc.read_h5ad('/home/caleb/Desktop/improvedgenepred/data/breastcancer_sample1_rep1_aligned_toxeniumimage/xeniumdata_xeniumimage_data.h5ad')
adata_visium = sc.read_h5ad('/home/caleb/Desktop/improvedgenepred/data/breastcancer_sample1_rep1_aligned_tovisiumimage/visiumdata_visiumimage_data.h5ad')

# make .X a csr matrix
adata_xenium.X = scipy.sparse.csr_matrix(adata_xenium.X)
adata_visium.X = scipy.sparse.csr_matrix(adata_visium.X)

# add array for gene expression
adata_xenium.X_array = pd.DataFrame(adata_xenium.X.toarray(), index=adata_xenium.obs.index)
adata_visium.X_array = pd.DataFrame(adata_visium.X.toarray(), index=adata_visium.obs.index)

# patches
with open('/home/caleb/Desktop/improvedgenepred/data/breastcancer_sample1_rep1_aligned_tovisiumimage/visiumdata_visiumimage_patches.pkl', 'rb') as f:
    aligned_visium_dictionary = pickle.load(f)

with open('/home/caleb/Desktop/improvedgenepred/data/breastcancer_sample1_rep1_aligned_toxeniumimage/xeniumdata_xeniumimage_patches.pkl', 'rb') as f:
    aligned_xenium_dictionary = pickle.load(f)



# Rerasterize the patches
aligned_visium_dictionary_pred = rerasterize_patches(aligned_visium_dictionary, 276)
aligned_xenium_dictionary_pred = rerasterize_patches(aligned_xenium_dictionary, 276)


# scale the data
scaling_factor = 1
for i in aligned_visium_dictionary_pred:
    # aligned_visium_dictionary[i].X_array = sc.pp.log1p(aligned_visium_dictionary[i].X_array * scaling_factor)
    aligned_visium_dictionary_pred[i].X = sc.pp.log1p(np.round(aligned_visium_dictionary_pred[i].X * scaling_factor))
    # aligned_visium_dictionary[i].X = sc.pp.scale(np.round(aligned_visium_dictionary[i].X * scaling_factor))

# adata_visium.X_array = adata_visium.X_array * scaling_factor

# scale the data
scaling_factor = 1
for i in aligned_xenium_dictionary_pred:
    # aligned_visium_dictionary[i].X_array = sc.pp.log1p(aligned_visium_dictionary[i].X_array * scaling_factor)
    aligned_xenium_dictionary_pred[i].X = sc.pp.log1p(np.round(aligned_xenium_dictionary_pred[i].X * scaling_factor))
    # aligned_xenium_dictionary[i].X = sc.pp.scale(np.round(aligned_xenium_dictionary[i].X * scaling_factor))



# plot
plotRaster(adata_xenium.uns['spatial'], aligned_xenium_dictionary_pred, color_by='total_expression')
plt.savefig("/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/final_figures_supplementary/supp_figure1/supp_fig1b1.svg", dpi=300, bbox_inches='tight')
plotRaster(adata_visium.uns['spatial'], aligned_visium_dictionary_pred, color_by='total_expression')
plt.savefig("/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/final_figures_supplementary/supp_figure1/supp_fig1b2.svg", dpi=300, bbox_inches='tight')




############################################################################################################

# import libraries
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
import scipy
import seaborn as sns
import os
import ast
import anndata
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
import scipy.sparse as sp
# find correlation between xenium and single-cell data for each cluster
import scipy.stats as stats
from sklearn.cluster import KMeans
from adjustText import adjust_text
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.collections as mcoll

# set seed
np.random.seed(0)


### read in xenium data ###

# read in the data
adata_xenium = sc.read_10x_h5('/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_xenium_sample1_rep1/cell_feature_matrix.h5')
# add array for gene expression
adata_xenium.X_array = pd.DataFrame(adata_xenium.X.toarray(), index=adata_xenium.obs.index)
# label genes
adata_xenium.X_array.columns = adata_xenium.var.index


### read in single-cell data ###


# read in single-cell data
file = "/home/caleb/Desktop/off-target-probe-checker/data/SC3pv3_GEX_Breast_Cancer_DTC_Aggr_count_filtered_feature_bc_matrix.h5"
sc_data = sc.read_10x_h5(file)


### compare total expression ###

# get shared genes
shared_genes = np.intersect1d(list(adata_xenium.var.index), list(sc_data.var.index))
len(shared_genes)

# get total expression
xenium_total = adata_xenium.X_array[shared_genes].sum(axis=0)
# make array
sc_data.X_array = pd.DataFrame(sc_data.X.toarray(), index=sc_data.obs.index)
# add gene names
sc_data.X_array.columns = sc_data.var.index

sc_total = sc_data.X_array[shared_genes].sum(axis=0)


### do harmonized clustering ###
## https://portals.broadinstitute.org/harmony/mudan.html

# get counts for shared genes
# Get counts for shared genes
xenium_counts = adata_xenium.X_array[shared_genes]
sc_counts = sc_data.X_array[shared_genes]


# Merge based on columns
cd = pd.concat([xenium_counts, sc_counts], axis=0)

cd.shape

meta = pd.Series(['xenium'] * xenium_counts.shape[0] + ['singlecell'] * sc_counts.shape[0])

# filter out single cells without any counts
vi = cd.sum(axis=1) > 1
cd = cd.loc[vi, :]
meta = np.array(meta)[vi]
pd.Series(meta).value_counts()

# NOTE: MERINGUE does log10 psuedocounts https://github.com/JEFworks-Lab/MERINGUE/blob/master/R/process.R
# do CPM normalization
cd = cd.div(cd.sum(axis=1), axis=0) * 1e6
# log transform
cd = np.log(cd + 1)

# make adata object
adata = ad.AnnData(cd)
adata.obs['batch'] = meta

# do PCA
sc.pp.pca(adata, n_comps=30)

# plot scree plot
sc.pl.pca_variance_ratio(adata, log=True)

# plot
sc.pl.pca(adata, color='batch')

# plot umap
sc.pp.neighbors(adata, use_rep='X_pca', n_pcs=30)
sc.tl.umap(adata)
# sc.pl.umap(adata, color='batch')
fig = sc.pl.umap(adata, color='batch', return_fig=True, title="", legend_loc='none')
# If there is more than one axis, remove the last one (assumed to be the colorbar)
if len(fig.axes) > 1:
    fig.axes[-1].remove()
# Iterate over all axes in the figure and turn them off
for ax in fig.get_axes():
    ax.set_axis_off()
fig.savefig("/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/final_figures_supplementary/supp_figure2/supp_fig2a_data.png", dpi=300, bbox_inches='tight')

fig = sc.pl.umap(adata, color='batch', return_fig=True, title="")
    # Iterate over all axes and remove scatter plot data (PathCollection objects)
for ax in fig.get_axes():
    for artist in ax.get_children():
        if isinstance(artist, mcoll.PathCollection):
            artist.remove()
fig.savefig("/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/final_figures_supplementary/supp_figure2/supp_fig2a_outline.svg", dpi=300, bbox_inches='tight')



# activate numpy2ri
rpy2.robjects.numpy2ri.activate()

# Import the R packages
harmony = importr('harmony')
meringue = importr('MERINGUE')
base = importr('base')

# set seed
ro.r('set.seed(0)')

# Convert to an R matrix 
mat = np.array(cd)
print(mat.shape)
cd_r = ro.r.matrix(mat.ravel(), nrow=mat.shape[0], ncol=mat.shape[1])

# Assign to R's global environment
ro.r.assign("cd", cd_r)

# 'cd' as if it were any R matrix
ro.r('dim(cd)  # Check dimensions in R')

# store your PCA results from Python in R, then call HarmonyMatrix
ro.r.assign("pc", adata.obsm['X_pca'])
ro.r.assign("meta", meta)

# run Harmony
ro.r('''
harmonized <- HarmonyMatrix(pc, meta, do_pca=FALSE, verbose=FALSE, theta=8)
dim(harmonized)
''')

# Get the harmonized matrix
harmonized_r = ro.r('harmonized')
# Convert to numpy array:
harmonized_np = np.array(harmonized_r)
print(harmonized_np.shape)

# Assign to AnnData object
adata.obsm['X_pca_harmony'] = harmonized_np

# do harmony. NOTE: harmonypy wasnt working as well as R version, so I switched
# sc.external.pp.harmony_integrate(adata, key='batch', basis='X_pca', adjusted_basis='X_pca_harmony', theta = 8, lamb=1, max_iter_harmony=50)

# plot Harmony-corrected PCA
sc.pl.embedding(adata, basis='X_pca_harmony', color='batch', components=['1,2'], frameon=False)

# plot harmony on umap
sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_pcs=30)
sc.tl.umap(adata)

fig = sc.pl.umap(adata, color='batch', return_fig=True, title="", legend_loc='none')
# If there is more than one axis, remove the last one (assumed to be the colorbar)
if len(fig.axes) > 1:
    fig.axes[-1].remove()
# Iterate over all axes in the figure and turn them off
for ax in fig.get_axes():
    ax.set_axis_off()
fig.savefig("/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/final_figures_supplementary/supp_figure2/supp_fig2b_data.png", dpi=300, bbox_inches='tight')

fig = sc.pl.umap(adata, color='batch', return_fig=True, title="")
    # Iterate over all axes and remove scatter plot data (PathCollection objects)
for ax in fig.get_axes():
    for artist in ax.get_children():
        if isinstance(artist, mcoll.PathCollection):
            artist.remove()
fig.savefig("/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/final_figures_supplementary/supp_figure2/supp_fig2b_outline.svg", dpi=300, bbox_inches='tight')



############################################################################################################


# perform clustering on harmony-corrected data
sc.tl.leiden(adata, resolution=1)

# do kmeans clustering
kmeans = KMeans(n_clusters=50, random_state=0).fit(adata.obsm['X_pca_harmony'])
adata.obs['kmeans'] = kmeans.labels_.astype(str)

adata.obs['kmeans'].value_counts()

# plot 
sc.pl.umap(adata, color='kmeans')


### get data without cpm normalization ###


# separate out the two datasets
adata_xen = adata[adata.obs['batch'] == 'xenium', :]
adata_sc = adata[adata.obs['batch'] == 'singlecell', :]

# subset sc data to only include cells with more than one count
vi2 = sc_counts.sum(axis=1) > 1
adata_sc_nocpm = sc_data[vi2, :]
# add umap
adata_sc_nocpm.obsm['X_umap'] = adata_sc.obsm['X_umap']
# do cpm normalization and log transform on X
# adata_sc_nocpm.X = adata_sc_nocpm.X / adata_sc_nocpm.X.sum(axis=1) * 1e6
adata_sc_nocpm.X = adata_sc_nocpm.X.toarray()
adata_sc_nocpm.X = np.log(adata_sc_nocpm.X + 1)
adata_sc_nocpm.X = scipy.sparse.csr_matrix(adata_sc_nocpm.X)
# add leiden clustering
adata_sc_nocpm.obs['leiden'] = adata_sc.obs['leiden']
adata_sc_nocpm.obs['kmeans'] = adata_sc.obs['kmeans']
adata_sc_nocpm.var_names_make_unique()
adata_sc_nocpm.uns['neighbors'] = adata_sc.uns['neighbors']


# copy xenium data
xenium_counts_copy = xenium_counts.copy()
# filter out single cells without any counts
vi5 = xenium_counts_copy.sum(axis=1) > 1
xenium_counts_copy = xenium_counts_copy.loc[vi, :]
# log transform
xenium_counts_copy = np.log(xenium_counts_copy + 1)
# make adata object
adata_xen_nocpm = ad.AnnData(xenium_counts_copy)
adata_xen_nocpm.obs['batch'] =  pd.Series(['xenium'] * adata_xen_nocpm.shape[0])
adata_xen_nocpm.obsm['X_umap'] = adata_xen.obsm['X_umap']
# add clustering
adata_xen_nocpm.obs['leiden'] = adata_xen.obs['leiden']
adata_xen_nocpm.obs['kmeans'] = adata_xen.obs['kmeans']
adata_xen_nocpm.var_names_make_unique()
adata_xen_nocpm.uns['neighbors'] = adata_xen.uns['neighbors']



### get correlation between xenium and single-cell data ###


# combine gene expression
bp_diff = "0" # NOTE: set this to 0, 5, or 10 in string form
offtarget_genes = pd.read_csv(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/gene_summary_out_{bp_diff}bp_nucmer_final.csv")
offtarget_genes



# Convert entire X matrix to a (cells x genes) DataFrame
# This can be memory-heavy for large data.
xen_expr = adata_xen_nocpm.to_df()  # rows=cells, columns=genes
xen_expr['kmeans'] = adata_xen_nocpm.obs['kmeans']  # attach cluster labels

# Similarly for single-cell data
sc_expr = adata_sc_nocpm.to_df()
sc_expr['kmeans'] = adata_sc_nocpm.obs['kmeans']
sc_expr.shape

# Build cluster-by-gene matrices
xen_cluster_expr = xen_expr.groupby('kmeans').mean()   # shape (#clusters, #genes)
sc_cluster_expr  = sc_expr.groupby('kmeans').mean()    # shape (#clusters, #genes)

# Align them by columns (genes)
common_genes = offtarget_genes['Gene']
xen_mat = xen_cluster_expr[common_genes].sort_index()
sc_mat  = sc_cluster_expr[common_genes].sort_index()

# Compute Pearson correlation for each gene
# save pearson corr
pearson_corr = []
for gene in common_genes:
    x_vec = xen_mat[gene]
    y_vec = sc_mat[gene]
    # create mask
    mask = ~np.isnan(x_vec) & ~np.isnan(y_vec)
    r, p = stats.pearsonr(x_vec[mask], y_vec[mask])
    pearson_corr.append(r)
    print(f"Gene {gene} cluster-profile correlation: r={r:.3f}, p={p}")



# get counts for all genes, this is raw data
sc_expr_agg = sc_data.X_array
sc_expr_agg.shape

# filter out single cells without any counts
vi3 = sc_counts.sum(axis=1) > 1
sc_expr_agg = sc_expr_agg.loc[vi3, :]
sc_expr_agg.shape



genes_agged = []
# now want to combine gene expression of off-target genes
for i in range(0, len(offtarget_genes)):
    # gene of interest
    gene = offtarget_genes.loc[i, 'Gene']
    # get list of off-targets NOTE this includes gene of interest
    off_targets = offtarget_genes.loc[i, 'genes_in_sc'].split(',')
    # make a new gene in the data
    genes_agged.append(gene)
    # replace gene with off-targets
    sc_expr_agg[gene] = sc_expr_agg[off_targets].sum(axis=1)


# log transform
sc_expr_agg = np.log(sc_expr_agg + 1)
sc_expr_agg.shape

# Similarly for single-cell data
sc_expr_agg['kmeans'] = adata_sc_nocpm.obs['kmeans']


# Build cluster-by-gene matrices
xen_cluster_expr = xen_expr.groupby('kmeans').mean()   # shape (#clusters, #genes)
sc_cluster_expr  = sc_expr_agg.groupby('kmeans').mean()    # shape (#clusters, #genes)

# Align them by columns (genes)
common_genes = offtarget_genes['Gene']
xen_mat = xen_cluster_expr[common_genes].sort_index()
sc_mat  = sc_cluster_expr[common_genes].sort_index()

# Compute Pearson correlation for each gene
# save pearson corr
pearson_corr_agg = []
for gene in common_genes:
    x_vec = xen_mat[gene]
    y_vec = sc_mat[gene]
    # create mask
    mask = ~np.isnan(x_vec) & ~np.isnan(y_vec)
    r, p = stats.pearsonr(x_vec[mask], y_vec[mask])
    pearson_corr_agg.append(r)
    print(f"Gene {gene} cluster-profile correlation: r={r:.3f}, p={p}")



### plot gene expression on umap and dot plots ###

# NOTE: need to set bp_diff above to to make sure oyu get agg correct

# combine gene expression
bp_diff = "0" # NOTE: set this to 0, 5, or 10 in string form
offtarget_genes = pd.read_csv(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/gene_summary_out_{bp_diff}bp_nucmer_final.csv")
offtarget_genes


# Create a custom colormap from gray to red
import matplotlib.colors as mcolors
gray_red = mcolors.LinearSegmentedColormap.from_list(
    'gray_red',
    ['#bfbfbf', '#ff0000']  # light gray to red
)

# make new adata with aggregated genes
adata_sc_nocpm_agg = adata_sc_nocpm.copy()
adata_sc_nocpm_agg.X = scipy.sparse.csr_matrix(sc_expr_agg.iloc[:, :-1])

# iterate over genes
for i in range(0, len(offtarget_genes)):
    # gene of interest
    gene = offtarget_genes.loc[i, 'Gene']

    # just get one gene
    if gene == 'CEACAM6':

        # get list of off-targets
        off_targets = offtarget_genes.loc[i, 'genes_in_sc'].split(',')
        # make a new directory for each gene
        os.makedirs(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/final_figures_supplementary/supp_figure_scplots/umaps_{bp_diff}bp/{gene}", exist_ok=True)

        ### xenium ###

        # plot umap
        fig = sc.pl.umap(adata_xen_nocpm, color=gene, color_map=gray_red, size=5, title = "", return_fig=True)
        # If there is more than one axis, remove the last one (assumed to be the colorbar)
        if len(fig.axes) > 1:
            fig.axes[-1].remove()
        # Iterate over all axes in the figure and turn them off
        for ax in fig.get_axes():
            ax.set_axis_off()
        plt.savefig(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/final_figures_supplementary/supp_figure_scplots/umaps_{bp_diff}bp/{gene}/{gene}_xen_umap_data.png", dpi=300, bbox_inches='tight')
        plt.close()

        # Create the UMAP plot and get the figure object
        fig = sc.pl.umap(adata_xen_nocpm, color=gene, color_map=gray_red, size=5, title = "Xenium Expression of " + gene, return_fig=True)

        # Iterate over all axes and remove scatter plot data (PathCollection objects)
        for ax in fig.get_axes():
            for artist in ax.get_children():
                if isinstance(artist, mcoll.PathCollection):
                    artist.remove()
        # save
        plt.savefig(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/final_figures_supplementary/supp_figure_scplots/umaps_{bp_diff}bp/{gene}/{gene}_xen_umap_outline.svg", dpi=300)
        plt.close()

        ### aggregated single-cell ###

        # plot umap of aggregated gene
        fig = sc.pl.umap(adata_sc_nocpm_agg, color=gene, color_map=gray_red, size=15, title = "", return_fig=True)
        # If there is more than one axis, remove the last one (assumed to be the colorbar)
        if len(fig.axes) > 1:
            fig.axes[-1].remove()
        # Iterate over all axes in the figure and turn them off
        for ax in fig.get_axes():
            ax.set_axis_off()
        plt.savefig(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/final_figures_supplementary/supp_figure_scplots/umaps_{bp_diff}bp/{gene}/{gene}_sc_agg_umap_data.png", dpi=300, bbox_inches='tight')
        plt.close()

        # plot umap of aggregated gene
        fig = sc.pl.umap(adata_sc_nocpm_agg, color=gene, color_map=gray_red, size=15, title = "Aggregated scRNA-seq Expression of " + gene, return_fig=True)
        # Iterate over all axes and remove scatter plot data (PathCollection objects)
        for ax in fig.get_axes():
            for artist in ax.get_children():
                if isinstance(artist, mcoll.PathCollection):
                    artist.remove()
        # save
        plt.savefig(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/final_figures_supplementary/supp_figure_scplots/umaps_{bp_diff}bp/{gene}/{gene}_sc_agg_umap_outline.svg", dpi=300)
        plt.close()
    
        ### single-cell ###

        # plot umap for each off-target
        for off_target in off_targets:
            fig = sc.pl.umap(adata_sc_nocpm, color=off_target, color_map=gray_red, size=15, title = "", return_fig=True)
            # If there is more than one axis, remove the last one (assumed to be the colorbar)
            if len(fig.axes) > 1:
                fig.axes[-1].remove()
            # Iterate over all axes in the figure and turn them off
            for ax in fig.get_axes():
                ax.set_axis_off()
            plt.savefig(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/final_figures_supplementary/supp_figure_scplots/umaps_{bp_diff}bp/{gene}/{off_target}_sc_umap_data.png", dpi=300, bbox_inches='tight')
            plt.close()

            # 
            fig = sc.pl.umap(adata_sc_nocpm, color=off_target, color_map=gray_red, size=15, title =  "scRNA-seq Expression of " + off_target, return_fig=True)
            # Iterate over all axes and remove scatter plot data (PathCollection objects)
            for ax in fig.get_axes():
                for artist in ax.get_children():
                    if isinstance(artist, mcoll.PathCollection):
                        artist.remove()
            
            plt.savefig(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/final_figures_supplementary/supp_figure_scplots/umaps_{bp_diff}bp/{gene}/{off_target}_sc_umap_outline.svg", dpi=300)
            plt.close()


############################################################################################################


### plot gene expression on umap and dot plots of other genes ###

# Create a custom colormap from gray to red
import matplotlib.colors as mcolors
gray_red = mcolors.LinearSegmentedColormap.from_list(
    'gray_red',
    ['#bfbfbf', '#ff0000']  # light gray to red
)

genes2plot = ['HDC']

# iterate over genes
for i in range(0, len(genes2plot)):
    print(i)
    # gene of interest
    gene = genes2plot[i]
    # make a new directory for each gene
    os.makedirs(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/final_figures_supplementary/supp_figure_scplots/umaps_notofftargets/{gene}", exist_ok=True)


    # xenium


    # plot umap
    fig = sc.pl.umap(adata_xen_nocpm, color=gene, color_map=gray_red, size=5, title="", return_fig=True)
    # If there is more than one axis, remove the last one (assumed to be the colorbar)
    if len(fig.axes) > 1:
        fig.axes[-1].remove()
    # Iterate over all axes in the figure and turn them off
    for ax in fig.get_axes():
        ax.set_axis_off()
    plt.savefig(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/final_figures_supplementary/supp_figure_scplots/umaps_notofftargets/{gene}/{gene}_xen_umap_data.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Create the UMAP plot and get the figure object
    fig = sc.pl.umap(adata_xen_nocpm, color=gene, color_map=gray_red, size=5, title="Xenium Expression of " + gene, return_fig=True)

    # Iterate over all axes and remove scatter plot data (PathCollection objects)
    for ax in fig.get_axes():
        for artist in ax.get_children():
            if isinstance(artist, mcoll.PathCollection):
                artist.remove()
    # save
    plt.savefig(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/final_figures_supplementary/supp_figure_scplots/umaps_notofftargets/{gene}/{gene}_xen_umap_outline.svg", dpi=300)
    plt.close()

    # scRNA-seq


    fig = sc.pl.umap(adata_sc_nocpm, color=gene, color_map=gray_red, size=15, title = "", return_fig=True)
    # If there is more than one axis, remove the last one (assumed to be the colorbar)
    if len(fig.axes) > 1:
        fig.axes[-1].remove()
    # Iterate over all axes in the figure and turn them off
    for ax in fig.get_axes():
        ax.set_axis_off()
    plt.savefig(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/final_figures_supplementary/supp_figure_scplots/umaps_notofftargets/{gene}/{gene}_sc_umap_data.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Create the UMAP plot and get the figure object
    fig = sc.pl.umap(adata_sc_nocpm, color=gene, color_map=gray_red, size=15, title = "scRNA-seq Expression of " + gene, return_fig=True)

        # Iterate over all axes and remove scatter plot data (PathCollection objects)
    for ax in fig.get_axes():
        for artist in ax.get_children():
            if isinstance(artist, mcoll.PathCollection):
                artist.remove()
    # save
    plt.savefig(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/final_figures_supplementary/supp_figure_scplots/umaps_notofftargets/{gene}/{gene}_sc_umap_outline.svg", dpi=300)
    plt.close()



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



### plot the gene expression of the off-targets ###

# log transform the data
sc.pp.log1p(adata_xenium)
sc.pp.log1p(adata_visium)


# New resolution that doesn;t cause issues
new_resolution_xenium = 275

# Rerasterize the patches
aligned_visium_dictionary_rerastered = rerasterize_patches(aligned_visium_dictionary, new_resolution_xenium)
aligned_xenium_dictionary_rerastered = rerasterize_patches(aligned_xenium_dictionary, new_resolution_xenium)


# combine gene expression
bp_diff = "0" # NOTE: set this to 0, 5, or 10 in string form
offtarget_genes = pd.read_csv(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting//gene_summary_out_{bp_diff}bp_nucmer_final.csv")
offtarget_genes


# turn off interactive mode so figures can be saved without displaying
plt.ioff()

# drfine your output folder
output_base = "/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/final_figures_supplementary/supp_figure_visplots/APOBEC3B/gene_plots_" + bp_diff + "bp"
os.makedirs(output_base, exist_ok=True)

# define the gene lists
gene_list_visium = list(adata_visium.var.index)
gene_list_xenium = list(adata_xenium.var.index)


# Loop over each row of dataframe
for idx, row in offtarget_genes[offtarget_genes['Gene'].isin(['APOBEC3B'])].iterrows():
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
    xenium_filename = os.path.join(gene_folder, f"{main_gene}_xenium.svg")
    fig.savefig(xenium_filename, bbox_inches='tight', dpi=300)
    plt.close(fig)
    
    # Plot Visium data for the main gene.
    fig = plotRaster(adata_visium.uns["spatial"], aligned_visium_dictionary_rerastered, color_by='gene_expression', gene_name=main_gene, if_vis=True)
    visium_main_filename = os.path.join(gene_folder, f"{main_gene}_visium.svg")
    fig.savefig(visium_main_filename, bbox_inches='tight', dpi=300)
    plt.close(fig)
    
    # For each off-target gene
    # Plot Visium data for the off-target
    if ',' in genes_in_vis:
        # for off_target in list(ast.literal_eval(genes_in_vis).keys()):
        for off_target in genes_in_vis.split(','):
            if off_target != main_gene:
                fig = plotRaster(adata_visium.uns["spatial"], aligned_visium_dictionary_rerastered, color_by='gene_expression', gene_name=off_target, if_vis=True)
                off_target_filename = os.path.join(gene_folder, f"{main_gene}_visium_offtarget_{off_target}.svg")
                fig.savefig(off_target_filename, bbox_inches='tight', dpi=300)
                plt.close(fig)

# print
print("Plots saved for all genes.")
# turn on interactive mode
plt.ion()





### plot aggregated gene expression ###

# combine gene expression
bp_diff = "0" # NOTE: set this to 0, 5, or 10 in string form
offtarget_genes = pd.read_csv(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/gene_summary_out_{bp_diff}bp_nucmer_final.csv")
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
output_base = "/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/final_figures_supplementary/supp_figure_visplots/APOBEC3B/gene_plots_" + bp_diff + "bp"
os.makedirs(output_base, exist_ok=True)

# define the gene lists
gene_list_visium = list(adata_visium.var.index)
gene_list_xenium = list(adata_xenium.var.index)

# offtarget_genes.iloc[17:,:]

# Loop over each row of dataframe
for idx, row in offtarget_genes[offtarget_genes['Gene'].isin(['APOBEC3B'])].iterrows():
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
    visium_main_filename = os.path.join(gene_folder, f"{main_gene}_visium_agg.svg")
    fig.savefig(visium_main_filename, bbox_inches='tight', dpi=300)
    plt.close(fig)

# print
print("Plots saved for all genes.")
# turn on interactive mode
plt.ion()


############################################################################################################


### plot other genes ###

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



# log transform the data
sc.pp.log1p(adata_xenium)
sc.pp.log1p(adata_visium)


# New resolution that doesn;t cause issues
new_resolution_xenium = 275

# Rerasterize the patches
aligned_visium_dictionary_rerastered = rerasterize_patches(aligned_visium_dictionary, new_resolution_xenium)
aligned_xenium_dictionary_rerastered = rerasterize_patches(aligned_xenium_dictionary, new_resolution_xenium)


# drfine your output folder
output_base = "/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/supp_figure_visplots/HDC/gene_plots_notofftarget"
os.makedirs(output_base, exist_ok=True)

# define the gene lists
gene_list_visium = list(adata_visium.var.index)
gene_list_xenium = list(adata_xenium.var.index)


genes2plot = ['HDC']

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
    xenium_filename = os.path.join(gene_folder, f"{main_gene}_xenium.svg")
    fig.savefig(xenium_filename, bbox_inches='tight', dpi=300)
    plt.close(fig)
    
    # Plot Visium data for the main gene.
    fig = plotRaster(adata_visium.uns["spatial"], aligned_visium_dictionary_rerastered, color_by='gene_expression', gene_name=main_gene, if_vis=True)
    visium_main_filename = os.path.join(gene_folder, f"{main_gene}_visium.svg")
    fig.savefig(visium_main_filename, bbox_inches='tight', dpi=300)
    plt.close(fig)
    