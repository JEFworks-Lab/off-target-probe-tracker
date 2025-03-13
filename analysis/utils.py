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
