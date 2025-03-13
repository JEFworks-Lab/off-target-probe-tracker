
############################################################################################################

# Compare sc-RNAseq and Xenium data


############################################################################################################

# import libraries
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
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

# set seed
np.random.seed(0)


############################################################################################################


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


############################################################################################################


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

# plot
plt.scatter(sc_total, xenium_total)
plt.xscale('log')
plt.yscale('log')
plt.ylabel('Xenium total expression')
plt.xlabel('Single-cell total expression')
plt.show()

# plot names of genes
plt.figure(figsize=(10, 10))
ax = sns.scatterplot(x=sc_total, y=xenium_total)
ax.set(xscale='log', yscale='log')
for i, txt in enumerate(shared_genes):
    ax.annotate(txt, (sc_total[i], xenium_total[i]))
plt.show()



############################################################################################################



### read in off-target data ###

# read in text file

# read in genes of interest
bp_diff = "0" # NOTE: set this to 0, 5, or 10 in string form
offtarget_genes = pd.read_csv(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/gene_summary_out_{bp_diff}bp_nucmer.csv")
offtarget_genes

# get probes
xenium_probes = pd.read_csv("/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/xenium_probes_list.csv")
list(xenium_probes['Gene'])

# get all genes
all_genes = pd.DataFrame(list(adata_xenium.var.index), columns=['Gene'])


df = all_genes.merge(offtarget_genes, on='Gene', how='left')

# add another column for if genes are in xenium_probes
df['dowehaveprobes'] = df['Gene'].isin(list(xenium_probes['Gene']))

df.index = df['Gene']


colors = []
for gene in list(df['Gene']):

    # gene = "ZEB2"
    aligned_to_list = (
        [] if pd.isna(df.loc[gene, "aligned_to"])
        else df.loc[gene, "aligned_to"].split(',')
    )
    
    genes_in_vis = df.loc[gene, "genes_in_vis"]

    # Condition for red: More than one gene in `aligned_to` OR a single gene that is different from the index
    is_offtarget = len(aligned_to_list) > 1 or (len(aligned_to_list) == 1 and aligned_to_list[0] != gene)

    # Blue **only within red genes**: If it meets red conditions *and* has `genes_in_vis`
    is_invisium = is_offtarget and gene != genes_in_vis

    # check if gene is in xenium_probes
    is_aprobe = df.loc[gene, "dowehaveprobes"]

    if is_invisium:
        colors.append('red')
    elif is_offtarget:
        colors.append('blue')
    elif not is_aprobe:
        colors.append('purple')
    else:
        colors.append('gray') 

# double check
print(pd.Series(colors).value_counts())


# plot
plt.figure(figsize=(10, 6))
# plot
plt.scatter(np.log(sc_total+1), np.log(xenium_total+1), c=colors)
plt.xlabel("Log Total Counts scRNA-seq", fontsize=16)
plt.ylabel("Log Total Counts Xenium", fontsize=16)
plt.title("Total counts of genes in scRNA-seq vs Xenium", fontsize=16)

plt.xlim([min(np.log(sc_total+1)) - 0.5, max(np.log(sc_total+1)) + 0.5])
plt.ylim([min(np.log(xenium_total+1)) - 0.5, max(np.log(xenium_total+1)) + 0.5])

# Plot reference lines, 
plt.plot(
    [min(np.log(sc_total+1)) - 0.5, max(np.log(sc_total+1)) + 0.5],
    [min(np.log(sc_total+1)) - 0.5, max(np.log(sc_total+1)) + 0.5],
    color='green', linewidth=2, linestyle='--'
)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# Label points that are red or blue
texts = []  # Store text elements for adjustText
for i, gene in enumerate(adata_xenium.var_names):
    if colors[i] in ['red']:
        txt = plt.text(
            np.log(sc_total+1)[i],
           np.log(xenium_total+1)[i],
            gene, fontsize=8
        )
        texts.append(txt)

# also add HDC in texts
txt = plt.text(
    np.log(sc_total+1)[np.log(sc_total+1).index == "HDC"],
    np.log(xenium_total+1)[np.log(xenium_total+1).index == "HDC"],
    'HDC', fontsize=8
)

# Adjust text to avoid overlap
adjust_text(texts, arrowprops=dict(arrowstyle="-", lw=0.5), force_points=0.5)

# Custom legend
red_patch = mpatches.Patch(color='red', label='Off Target Probes in Visium')
blue_patch = mpatches.Patch(color='blue', label='Off Target Probes no Visium')
green_patch = mpatches.Patch(color='purple', label='Custom Probes')
gray_patch = mpatches.Patch(color='gray', label='No issues')
plt.legend(handles=[red_patch, blue_patch, green_patch, gray_patch])

sns.despine()
plt.show()











############################################################################################################


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
sc.pl.umap(adata, color='batch')


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
sc.pl.umap(adata, color='batch')


############################################################################################################


# perform clustering on harmony-corrected data
sc.tl.leiden(adata, resolution=1)

# do kmeans clustering
kmeans = KMeans(n_clusters=50, random_state=0).fit(adata.obsm['X_pca_harmony'])
adata.obs['kmeans'] = kmeans.labels_.astype(str)

adata.obs['kmeans'].value_counts()

# plot 
sc.pl.umap(adata, color='kmeans')


############################################################################################################


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


################################################################################################################


### get correlation between xenium and single-cell data ###


# combine gene expression
bp_diff = "6" # NOTE: set this to 0, 5, or 10 in string form
offtarget_genes = pd.read_csv(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/gene_summary_out_{bp_diff}bp_nucmer_new.csv")
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




# plot correlation
plt.figure(figsize=(10, 5))
plt.bar(common_genes, pearson_corr, color='#0406EA', alpha=0.5, label='Single Genes')
plt.bar(common_genes, pearson_corr_agg, color='#F2080A', alpha=0.5, label='Aggregated Genes')
plt.xticks(rotation=45)
plt.ylabel('Pearson Correlation')
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.legend()
plt.title('Pearson Correlation between SC Data and Xenium Data')
sns.despine()
plt.show()





############################################################################################################


### plot gene expression on umap and dot plots ###

# NOTE: need to set bp_diff above to to make sure oyu get agg correct

# combine gene expression
bp_diff = "6" # NOTE: set this to 0, 5, or 10 in string form
offtarget_genes = pd.read_csv(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/gene_summary_out_{bp_diff}bp_nucmer_new.csv")
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
    print(i)
    # gene of interest
    gene = offtarget_genes.loc[i, 'Gene']
    # get list of off-targets
    off_targets = offtarget_genes.loc[i, 'genes_in_sc'].split(',')
    # make a new directory for each gene
    os.makedirs(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/umaps_dotplots_{bp_diff}bp/{gene}", exist_ok=True)
    # plot umap
    sc.pl.umap(adata_xen_nocpm, color=gene, color_map=gray_red, size=5, title = "Xenium: " + gene, return_fig=True)
    plt.savefig(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/umaps_dotplots_{bp_diff}bp/{gene}/{gene}_xen_umap.png")
    plt.close()

    # plot umap of aggregated gene
    sc.pl.umap(adata_sc_nocpm_agg, color=gene, color_map=gray_red, size=10, title = "Single Cell Aggregated: " + gene, return_fig=True)
    plt.savefig(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/umaps_dotplots_{bp_diff}bp/{gene}/{gene}_sc_agg_umap.png")
    plt.close()

    # plot umap for each off-target
    for off_target in off_targets:
        sc.pl.umap(adata_sc_nocpm, color=off_target, color_map=gray_red, size=10, title = "Single Cell: " + off_target, return_fig=True)
        plt.savefig(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/umaps_dotplots_{bp_diff}bp/{gene}/{off_target}_sc_umap.png")
        plt.close()


############################################################################################################



### plot gene expression on umap and dot plots of other genes ###

# Create a custom colormap from gray to red
import matplotlib.colors as mcolors
gray_red = mcolors.LinearSegmentedColormap.from_list(
    'gray_red',
    ['#bfbfbf', '#ff0000']  # light gray to red
)

genes2plot = ['TACSTD2']

# iterate over genes
for i in range(0, len(genes2plot)):
    print(i)
    # gene of interest
    gene = genes2plot[i]
    # make a new directory for each gene
    os.makedirs(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/umaps_dotplots_notofftargets/{gene}", exist_ok=True)
    # plot umap
    sc.pl.umap(adata_xen_nocpm, color=gene, color_map=gray_red, size=5, title = "Xenium: " + gene, return_fig=True)
    plt.savefig(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/umaps_dotplots_notofftargets/{gene}/{gene}_xen_umap.png")
    plt.close()

    sc.pl.umap(adata_sc_nocpm, color=gene, color_map=gray_red, size=10, title = "Single Cell: " + gene, return_fig=True)
    plt.savefig(f"/home/caleb/Desktop/off-target-probe-checker/caleb/plotting/umaps_dotplots_notofftargets/{gene}/{gene}_sc_umap.png")
    plt.close()


############################################################################################################

