
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# Function to plot patches on the original image
def plotRaster(image, adata_patches, color_by='gene_expression', gene_name=None):
    """
    Plots patches on the original image, colored by either gene expression or a column in adata_patches.obs.

    Parameters:
    - image: The original image array.
    - adata_patches: Dictionary of AnnData objects representing the patches.
    - color_by: How to color the patches ('gene_expression' or 'total_expression').
    - gene_name: The name of the gene to use if color_by is 'gene_expression'.
    """
    # Check inputs
    if color_by == 'gene_expression' and gene_name is None:
        raise ValueError("You must specify a gene_name when color_by='gene_expression'.")

    # Collect all values for normalization
    values = []
    for adata_patch in adata_patches.values():
        if color_by == 'gene_expression':
            expression = adata_patch.X[:, adata_patch.var_names.get_loc(gene_name)].sum()
            values.append(expression)
        elif color_by == 'total_expression':
            total_expression = adata_patch.X.sum()
            values.append(total_expression)
    
    # get min and max values
    values = np.array(values)
    min_value, max_value = values.min(), values.max()

    # Plot the original image
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(image)

    # Plot each patch with the appropriate color
    for adata_patch in adata_patches.values():
        x_start, x_end, y_start, y_end = adata_patch.uns['patch_coords']
        
        if color_by == 'gene_expression':
            expression = adata_patch.X[:, adata_patch.var_names.get_loc(gene_name)].sum()
            normalized_value = (expression - min_value) / (max_value - min_value)
            color = plt.cm.viridis(normalized_value)
        elif color_by == 'total_expression':
            total_expression = adata_patch.X.sum()
            normalized_value = (total_expression - min_value) / (max_value - min_value)
            color = plt.cm.viridis(normalized_value)
        
        # Draw a rectangle for the patch
        rect = mpatches.Rectangle((x_start, y_start), x_end - x_start, y_end - y_start,
                                  linewidth=1, edgecolor='none', facecolor=color, alpha=1)
        ax.add_patch(rect)

    # Create a color bar
    norm = plt.Normalize(min_value, max_value)
    sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation='horizontal', fraction=0.1, pad=0.02, shrink = .75)
    # cbar.set_label(f'{gene_name} Expression' if color_by == 'gene_expression' else "total_expression")

    plt.axis('off')


def split_adata_patches(combined_adata, adata_patches_for_image):
    # Dictionary to store the split AnnData objects
    adata_patches = {}
    
    # Get the unique keys (patch identifiers) from the obs index
    unique_keys = combined_adata.obs.index.unique()
    
    # Iterate over the unique keys to split the data
    for key in unique_keys:
        # Subset the combined AnnData for each key
        adata_patch = combined_adata[combined_adata.obs.index == key].copy()
        # Reset the obs index to default to avoid confusion
        adata_patch.obs.index = adata_patch.obs.index.str.split('-').str[-1]

        adata_patch.uns['spatial'] = adata_patches_for_image[key].uns['spatial']
        adata_patch.uns['patch_coords'] = adata_patches_for_image[key].uns['patch_coords']
        # Store the individual AnnData object in the dictionary
        adata_patches[key] = adata_patch
    
    return adata_patches



import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
import anndata

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
        print(idx)
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



# Function to plot patches on two images side by side with a shared color bar
def plotRasterSideBySide(image1, adata_patches1, image2, adata_patches2, color_by='gene_expression', gene_name=None):
    """
    Plots patches on two images side by side, colored by either gene expression or a column in adata_patches.obs.
    A single shared heatmap legend is used.

    Parameters:
    - image1: The first original image array.
    - adata_patches1: Dictionary of AnnData objects representing the patches for the first image.
    - image2: The second original image array.
    - adata_patches2: Dictionary of AnnData objects representing the patches for the second image.
    - color_by: How to color the patches ('gene_expression' or 'total_expression').
    - gene_name: The name of the gene to use if color_by is 'gene_expression'.
    """
    if color_by == 'gene_expression' and gene_name is None:
        raise ValueError("You must specify a gene_name when color_by='gene_expression'.")

    # Collect values across both adata_patches1 and adata_patches2 for normalization
    values = []
    for adata_patch in list(adata_patches1.values()) + list(adata_patches2.values()):
        if color_by == 'gene_expression':
            expression = adata_patch.X[:, adata_patch.var_names.get_loc(gene_name)].sum()
            values.append(expression)
        elif color_by == 'total_expression':
            total_expression = adata_patch.X.sum()
            values.append(total_expression)
    
    # Determine color normalization range
    values = np.array(values)
    min_value, max_value = values.min(), values.max()

    # Set up subplots for side-by-side images
    fig, axes = plt.subplots(1, 2, figsize=(20, 10))
    
    # Helper function to plot patches on each axis
    def plot_patches_on_image(ax, image, adata_patches, title=''):
        ax.imshow(image)
        for adata_patch in adata_patches.values():
            x_start, x_end, y_start, y_end = adata_patch.uns['patch_coords']
            if color_by == 'gene_expression':
                expression = adata_patch.X[:, adata_patch.var_names.get_loc(gene_name)].sum()
                normalized_value = (expression - min_value) / (max_value - min_value)
                color = plt.cm.viridis(normalized_value)
            elif color_by == 'total_expression':
                total_expression = adata_patch.X.sum()
                normalized_value = (total_expression - min_value) / (max_value - min_value)
                color = plt.cm.viridis(normalized_value)
            
            # Draw rectangle for patch
            rect = mpatches.Rectangle((x_start, y_start), x_end - x_start, y_end - y_start,
                                      linewidth=1, edgecolor='none', facecolor=color, alpha=1)
            ax.add_patch(rect)
        # add title
        ax.set_title(title)
        ax.axis('off')

    # Plot patches on the first image
    plot_patches_on_image(axes[0], image1, adata_patches1, title='Ground Truth')
    
    # Plot patches on the second image
    plot_patches_on_image(axes[1], image2, adata_patches2, title='Predicted')
    
    # Create a single color bar for both images
    norm = plt.Normalize(min_value, max_value)
    sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, orientation='horizontal', fraction=0.03, pad=0.04)
    cbar.set_label(f'{gene_name} Expression' if color_by == 'gene_expression' else "Total Expression")

    plt.show()

    

def transform_coordinates(x, y, image_width, image_height, rotation_k=0, flip_axis=None):
    """
    Transforms (x, y) coordinates to match image transformations using np.flip first, then np.rot90.
    
    Parameters:
    - x: Array of original x-coordinates.
    - y: Array of original y-coordinates.
    - image_width: Width of the image.
    - image_height: Height of the image.
    - rotation_k: Number of 90-degree rotations counterclockwise (0, 1, 2, or 3).
    - flip_axis: Axis to flip (None, 0 for vertical, 1 for horizontal).
    
    Returns:
    - x_new: Transformed x-coordinates (array of integers).
    - y_new: Transformed y-coordinates (array of integers).
    """
    # Ensure x and y are numpy arrays
    x = np.asarray(x)
    y = np.asarray(y)

    # Step 1: Apply flipping using np.flip
    if flip_axis == 0:  # Vertical flip
        y = image_height - 1 - y
    elif flip_axis == 1:  # Horizontal flip
        x = image_width - 1 - x

    # Step 2: Apply rotation using np.rot90
    if rotation_k % 4 == 1:  # 90 degrees counterclockwise
        x_new = image_height - 1 - y
        y_new = x
    elif rotation_k % 4 == 2:  # 180 degrees
        x_new = image_width - 1 - x
        y_new = image_height - 1 - y
    elif rotation_k % 4 == 3:  # 270 degrees counterclockwise (90 degrees clockwise)
        x_new = y
        y_new = image_width - 1 - x
    else:  # rotation_k % 4 == 0, no rotation
        x_new, y_new = x, y

    # Ensure the final coordinates are integers
    x_new = np.round(x_new).astype(int)
    y_new = np.round(y_new).astype(int)

    return x_new, y_new

