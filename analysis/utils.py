'''
Misc functions for the pipeline
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from anndata import AnnData
import time
from scipy.sparse import csr_matrix
import anndata as ad




####################################################################################################


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


####################################################################################################


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


####################################################################################################


# Function to subset and aggregate AnnData by bounding box coordinates
def subset_and_aggregate_patch_basedoncenters(adata, image, x_start, x_end, y_start, y_end, used_cells, aggregation='mean', visium=False):
    """Subset an AnnData object based on a spatial range and aggregate the data, ensuring cells are only included in the first patch they appear in."""
    # Extract spatial coordinates
    spatial_coords = adata.obsm["spatial"]

    # filter spots within the bounding box and not already used
    mask = (
        (spatial_coords[:, 0] >= x_start) & (spatial_coords[:, 0] < x_end) &
        (spatial_coords[:, 1] >= y_start) & (spatial_coords[:, 1] < y_end)
    )
    
    # Remove cells that have already been used
    mask = mask & (~adata.obs.index.isin(used_cells))

    # Subset the AnnData object based on the mask
    adata_patch = adata[mask, :]

    # Return None if there are no cells in the patch
    if adata_patch.shape[0] == 0:
        return None

    # Add these cells to the set of used cells
    used_cells.update(adata_patch.obs.index)

    # Aggregate the data within the patch
    if aggregation == 'sum':
        aggregated_data = adata_patch.X.sum(axis=0)
    elif aggregation == 'mean':
        aggregated_data = adata_patch.X.mean(axis=0)
    else:
        raise ValueError("Invalid aggregation method. Use 'sum' or 'mean'.")

    # Create a new AnnData object with aggregated data
    aggregated_data = aggregated_data if isinstance(aggregated_data, csr_matrix) else csr_matrix(aggregated_data)
    new_adata = ad.AnnData(X=aggregated_data)
    
    # Add image patch
    new_adata.uns['spatial'] = image[y_start:y_end, x_start:x_end]
    # Add patch coordinates
    new_adata.uns['patch_coords'] = [x_start, x_end, y_start, y_end]
    
    # Add centroid of new patch
    new_adata.obs['x_centroid'] = (x_start + x_end) / 2
    new_adata.obs['y_centroid'] = (y_start + y_end) / 2

    if visium:
        for field in ['in_tissue', 'array_row', 'array_col']:
            new_adata.obs[field] = adata_patch.obs[field].iloc[0]

    # Add spatial coordinates
    new_adata.obsm["spatial"] = new_adata.obs[["x_centroid", "y_centroid"]].to_numpy().astype(int)

    # Add variables and gene names
    new_adata.var = adata.var
    new_adata.var_names = adata.var_names
    # make sure X is a sparse matrix
    new_adata.X = csr_matrix(new_adata.X)

    return new_adata

# Function to extract patches and aggregate data from an image and AnnData object based on supplied center coordinates
def rasterizeGeneExpression_topatches_basedoncenters(image, adata, center_coords, patch_size=100, aggregation='mean', visium=False):
    """Extract patches centered around supplied coordinates from an image and aggregate AnnData data accordingly."""

    # Initialize variables
    adata_sub_dict = {}
    img_height, img_width, _ = image.shape
    used_cells = set()

    # Loop through each center coordinate
    for patch_index, (x_center, y_center) in enumerate(center_coords):
        # Calculate bounding box around the center coordinate
        x_start = max(0, x_center - patch_size // 2)
        x_end = min(img_width, x_center + patch_size // 2)
        y_start = max(0, y_center - patch_size // 2)
        y_end = min(img_height, y_center + patch_size // 2)

        # Subset and aggregate the AnnData object
        adata_patch = subset_and_aggregate_patch_basedoncenters(adata, image, x_start, x_end, y_start, y_end, used_cells, aggregation, visium)
        
        # Filter out empty patches
        if adata_patch is not None:
            if adata_patch.uns['spatial'].shape == (patch_size, patch_size, 3):
                patch_name = f"patch_{patch_index}"
                adata_sub_dict[patch_name] = adata_patch

    # return the dictionary of patches
    return adata_sub_dict


####################################################################################################