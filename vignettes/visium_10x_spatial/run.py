import voyagerpy as vp

# Data handling and plotting
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# File management and downloading
import os
import tarfile
import requests
import pathlib
import json
import gget  # mapping ensembl gene IDs to gene symbols

# Graph/matrix handling
import networkx as nx  # graph generation
from scipy.sparse import csr_matrix, issparse  # matrix handling
import scanpy as sc  # PCA and KNN

# Image handling
from scipy.ndimage import binary_opening, binary_closing, binary_fill_holes
from skimage.measure import label, regionprops
from skimage import color, morphology
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from shapely.affinity import scale
from matplotlib.lines import Line2D
from shapely.geometry import Point
from shapely.predicates import contains
from scipy.stats import gaussian_kde

from utils import Checkpoint
checkpoint = Checkpoint()
# Set the relevant directory. This depends on the location of the downloaded data.
outs_dir = pathlib.Path(checkpoint.root_dir / "outs")
# Use this directory as input to the VoyagerPy read* function.

plt.rcParams['axes.grid'] = False
plt.rcParams['axes.edgecolor'] = 'none'  # Hides the axes border
# plt.rcParams['xtick.color'] = 'none'  # Hides x-axis ticks
# plt.rcParams['ytick.color'] = 'none'  # Hides y-axis ticks

plt.rcParams['figure.dpi'] = 120

# Turn on matplotlib interactive mode so we don't need to explicitly call plt.show()
_ = plt.ion()

adata_image_res = 'hires'

root_dir = "/workspace/vignettes/visium_10x_spatial"  #! pathlib.Path(os.getcwd())
data_dir = pathlib.Path(f'{root_dir}/data')
outs_dir = pathlib.Path(f'{root_dir}/outs')

adata = vp.read_10x_visium(
    outs_dir,
    datatype = 'mtx',
    raw = True,
    prefix = None,
    symbol_as_index=False,
    dtype=np.float64,
    res=adata_image_res
)

is_mt = adata.var['symbol'].str.contains('^mt-').values
vp.utils.add_per_cell_qcmetrics(adata, subsets={'mito': is_mt})

visium_spots = vp.spatial.get_visium_spots(adata, with_radius=True)

_ = vp.spatial.set_geometry(adata, geom="spot_poly", values=visium_spots)

qc_features = ["detected", "subsets_mito_percent", "sum"]

image_path = f'{outs_dir}/spatial/tissue_hires_image.png'
im = plt.imread(image_path)

json_path = outs_dir / 'spatial' / 'scalefactors_json.json'

img2 = color.rgb2gray(im)

mask = im[:,:,2] < 0.87

kern = np.ones((3, 3), dtype=np.uint8)

mask_open = binary_opening(mask, structure=kern)

mask_close = binary_closing(mask_open, kern)

mask_label = label(mask_close, connectivity=2)
properties = regionprops(mask_label)

def radii_from_centroid(region):
    centroid = np.array(region.centroid)
    coords = np.array(region.coords)
    distances = np.linalg.norm(coords - centroid, axis=1)
    return distances

# Print area and perimeter for each region
for prop in properties[:6]:
    radii = radii_from_centroid(prop)

max_area = 0
max_label = None
for prop in properties:
    if prop.area > max_area:
        max_area = prop.area
        max_label = prop.label

if max_label is not None:
    inds = np.argwhere(mask_label == max_label)

region = next((prop for prop in properties if prop.label == max_label), None)

if region is not None:
    # Extract the bounding box of the region with the max area
    minr, minc, maxr, maxc = region.bbox

    # Compute indices to zero out
    # Rows outside the bounding box
    row_inds = np.concatenate((np.arange(minr), np.arange(maxr + 1, mask_label.shape[0])))
    # Columns outside the bounding box
    col_inds = np.concatenate((np.arange(minc), np.arange(maxc + 1, mask_label.shape[1])))

    # Set these regions in the mask to zero
    mask_label[row_inds, :] = 0
    mask_label[:, col_inds] = 0

np.unique(mask_label)

unique_labels = np.unique(mask_label)[1:]
filtered_props = [prop for prop in properties if prop.label in unique_labels]

# Create a DataFrame from the filtered properties
data = {'label': [prop.label for prop in filtered_props], 'area': [prop.area for prop in filtered_props]}

df = pd.DataFrame(data)

# Sort by area
df_sorted = df.sort_values(by='area', ascending=False)

mask_label[mask_label == 793] = 0

for unique_label in unique_labels:
    if int(df_sorted['area'][df_sorted['label'] == unique_label]) < 100:
        mask_label[mask_label == unique_label] = 0

mask_label = binary_fill_holes(mask_label)

# Create outlines by subtracting eroded version from the mask
eroded = morphology.erosion(mask_label)
outlines = np.logical_xor(mask_label, eroded)

# Create an RGBA image based on the grayscale image
rgba_img = color.gray2rgba(img2)

# Color map: Modify this part as per the number of labels and required colors
colors = {
    1: (1, 1, 0, 0.3),  # Yellow fill with 30% opacity
    2: (1, 0, 0, 1)     # Red for outlines with full opacity
}

# Apply colors: fill and outlines
for label, color_value in colors.items():
    if label == 1:
        mask = mask_label > 0
        rgba_img[mask] = rgba_img[mask] * (1 - color_value[3]) + np.array(color_value[:3] + (0,)) * color_value[3]
    elif label == 2:
        rgba_img[outlines > 0] = color_value

adata.uns["spatial"]["geom"] = {}
adata.uns["spatial"]["geom"]["tissueBoundary_hires"] = vp.spatial.get_tissue_boundary(adata, size = adata_image_res)

with json_path.open() as file:
    data = json.load(file)
    lowres_boundary_scale_factor = data['tissue_lowres_scalef']
    hires_boundary_scale_factor = data['tissue_hires_scalef']
    spot_scale_factor = data['spot_diameter_fullres']

adata.uns["spatial"]["geom"]["tissueBoundary"] = scale(adata.uns["spatial"]["geom"]["tissueBoundary_hires"], xfact=1/hires_boundary_scale_factor, yfact = 1/hires_boundary_scale_factor, origin = (0,0))
adata.uns["spatial"]["geom"]["tissueBoundary_lowres"] = scale(adata.uns["spatial"]["geom"]["tissueBoundary"], xfact=lowres_boundary_scale_factor, yfact = lowres_boundary_scale_factor, origin = (0,0))
adata.uns["spatial"]["geom"]["tissueBoundary_hires_exterior"] = adata.uns["spatial"]["geom"]["tissueBoundary_hires"].exterior

adata.uns["spatial"]["img"][adata_image_res] = np.swapaxes(adata.uns["spatial"]["img"][adata_image_res], 0, 1)
transposed_im = np.transpose(im, (1, 0, 2))

# adata.uns["spatial"]["img"][adata_image_res] = vp.spatial.rotate_img90(adata, k = 2)[adata_image_res]
# transposed_im = np.rot90(transposed_im, 2)

adata.uns["spatial"]["geom"]["tissueBoundary_hires"] = vp.spatial.get_tissue_boundary(adata, size = adata_image_res)
adata.uns["spatial"]["geom"]["tissueBoundary"] = scale(adata.uns["spatial"]["geom"]["tissueBoundary_hires"], xfact=1/hires_boundary_scale_factor, yfact = 1/hires_boundary_scale_factor, origin = (0,0))
adata.uns["spatial"]["geom"]["tissueBoundary_lowres"] = scale(adata.uns["spatial"]["geom"]["tissueBoundary"], xfact=lowres_boundary_scale_factor, yfact = lowres_boundary_scale_factor, origin = (0,0))
adata.uns["spatial"]["geom"]["tissueBoundary_hires_exterior"] = adata.uns["spatial"]["geom"]["tissueBoundary_hires"].exterior

spatial_coords_df = adata.obs[['pxl_row_in_fullres', 'pxl_col_in_fullres']]

adata.obs['pxl_col_in_fullres_transposed'] = adata.obs.apply(
    lambda row: row['pxl_row_in_fullres'],
    axis=1
)
adata.obs['pxl_row_in_fullres_transposed'] = adata.obs.apply(
    lambda row: row['pxl_col_in_fullres'],
    axis=1
)

adata.obs.rename(columns={
    'pxl_row_in_fullres': 'pxl_row_in_fullres_upright',
    'pxl_col_in_fullres': 'pxl_col_in_fullres_upright',
    'pxl_col_in_fullres_transposed': 'pxl_col_in_fullres',
    'pxl_row_in_fullres_transposed': 'pxl_row_in_fullres'
}, inplace=True)

visium_spots = vp.spatial.get_visium_spots(adata, with_radius=True)
_ = vp.spatial.set_geometry(adata, geom="spot_poly", values=visium_spots)

def check_intersection(row):
    # Create a circle at the given coordinates with the specified radius
    circle = Point(row['pxl_col_in_fullres'], row['pxl_row_in_fullres']).buffer(0.5*spot_scale_factor)

    # Check if the circle intersects any polygon
    found_intersection = circle.intersects(adata.uns["spatial"]["geom"]["tissueBoundary"])

    return found_intersection

def check_contain(row):
    # Create a circle at the given coordinates with the specified radius
    circle = Point(row['pxl_col_in_fullres'], row['pxl_row_in_fullres']).buffer(0.5*spot_scale_factor)

    # Check if the circle intersects any polygon
    found_contain = contains(adata.uns["spatial"]["geom"]["tissueBoundary"], circle)

    return found_contain

adata.obs['intersects_boundary'] = adata.obs.apply(check_intersection, axis=1)
adata.obs['contained_in_boundary'] = adata.obs.apply(check_contain, axis=1)

# CHECKPOINT: intersection_and_contain_df
intersection_and_contain_df = adata.obs[['intersects_boundary', 'contained_in_boundary']]
checkpoint.add("intersection_and_contain_df.csv", intersection_and_contain_df)
# intersection_and_contain_df = checkpoint.sync("intersection_and_contain_df.csv", intersection_and_contain_df, "csv")

def evaluate_conditions(row):
    if bool(row['in_tissue']) and not bool(row['intersects_boundary']):
        return "Space Ranger"
    elif bool(row['intersects_boundary']) and not bool(row['in_tissue']):
        return "Segmentation"
    else:
        return "Same"
    
adata.obs['diff_sr'] = adata.obs.apply(evaluate_conditions, axis=1)

adata.obs['diff_sr'] = adata.obs['diff_sr'].astype('category')

def intersect_without_contain(row):
    if bool(row['intersects_boundary']) and not bool(row['contained_in_boundary']):
        return True
    else:
        return False
    
adata.obs['intersect_without_contain'] = adata.obs.apply(intersect_without_contain, axis=1)

adata.obs['intersect_without_contain'] = adata.obs['intersect_without_contain'].astype('category')


def check_fraction_visium_spot_in_tissue(row):
    # Create a circle at the given coordinates with the specified radius
    circle = Point(row['pxl_col_in_fullres'], row['pxl_row_in_fullres']).buffer(0.5*spot_scale_factor)

    # Check if the circle intersects any polygon
    intersection = circle.intersection(adata.uns["spatial"]["geom"]["tissueBoundary"])

    fraction_in_tissue = intersection.area / circle.area

    return fraction_in_tissue

adata.obs['fraction_in_tissue'] = adata.obs.apply(check_fraction_visium_spot_in_tissue, axis=1)

adata_raw = adata.copy()

mask = adata.obs['fraction_in_tissue'] > 0

# Subset the AnnData object using the mask
adata = adata[mask].copy()

visium_spots = vp.spatial.get_visium_spots(adata, with_radius=True)
_ = vp.spatial.set_geometry(adata, geom="spot_poly", values=visium_spots)

nn_graph = 'visium_graph'
nn_graph_object = vp.spatial.find_visium_graph(
    adata,
    subset=None,
    geom=None,
    graph_key=None,
    inplace=False,  # This is important as we don't want to modify the original graph
)

sparse_matrix = csr_matrix(nx.to_scipy_sparse_array(nn_graph_object, nodelist=sorted(nn_graph_object.nodes())))

adata.obsp[nn_graph] = sparse_matrix

# CHECKPOINT: visium_graph.mtx
# checkpoint.add("visium_graph.mtx", adata.obsp[nn_graph])
# adata.obsp[nn_graph] = checkpoint.sync("visium_graph.mtx", adata.obsp[nn_graph], "mtx")

visium_graph = vp.plotting.draw_graph(adata, nn_graph)
vp.spatial.set_default_graph(adata, nn_graph)
vp.spatial.to_spatial_weights(adata, graph_name=nn_graph)

vp.spatial.moran(adata, feature = qc_features, graph_name = nn_graph, permutations = 200)
adata.uns['spatial']['moran'][nn_graph].loc[qc_features, ["I"]]

# CHECKPOINT: qc_moran
checkpoint.add("qc_moran.csv", adata.uns['spatial']['moran'][nn_graph].loc[qc_features, ["I"]])

vp.spatial.compute_higher_order_neighbors(adata, order=8)
vp.spatial.compute_correlogram(adata, qc_features, order=8)

# CHECKPOINT: qc_correlogram
checkpoint.add("qc_correlogram.csv", adata.uns["spatial"]["moran"]["correlogram"][nn_graph])
# adata.uns["spatial"]["moran"]["correlogram"][nn_graph] = checkpoint.sync("visium_correlogram.csv", adata.uns["spatial"]["moran"]["correlogram"][nn_graph], "csv")

_ = vp.spatial.local_moran(adata, feature = qc_features)

# CHECKPOINT: qc_local_moran
checkpoint.add("qc_local_moran.csv", adata.obsm['local_moran'])

_ = vp.spatial.losh(adata, feature = qc_features)

# CHECKPOINT: qc_losh
checkpoint.add("qc_losh.csv", adata.obsm['losh'])
# adata.obsm['losh'] = checkpoint.sync("visium_losh.csv", adata.obsm['losh'], "csv")

adata.layers['counts'] = adata.X.copy()
vp.utils.log_norm_counts(adata, inplace=True)
adata.layers['logcounts'] = adata.X.copy()

# CHECKPOINT: logcounts.mtx
checkpoint.add("logcounts.mtx", adata.layers['logcounts'])
# adata.layers['logcounts'] = checkpoint.sync("logcounts.mtx", adata.layers['logcounts'], "mtx")

_ = vp.spatial.set_geometry(adata, geom="spot_poly", values=visium_spots)

gene_var = vp.utils.model_gene_var(adata.layers['logcounts'], gene_names=adata.var_names)

# CHECKPOINT: gene_var
checkpoint.add("gene_var.txt", gene_var.index.tolist())
# gene_var = checkpoint.sync("gene_var.csv", gene_var, "csv")

hvgs = vp.utils.get_top_hvgs(gene_var)

# CHECKPOINT AND SYNC: hvgs
checkpoint.add("hvgs.txt", hvgs)
hvgs = checkpoint.sync("hvgs.txt", hvgs, "txt")

adata.var['highly_variable'] = False
adata.var.loc[hvgs, 'highly_variable'] = True

vp.spatial.moran(adata, feature=hvgs, dim='var', graph_name=nn_graph)

# CHECKPOINT: hvgs_moran
checkpoint.add("hvg_moran.csv", adata.uns['spatial']['moran'][nn_graph].loc[hvgs, ["I"]])

hvgs_moransI = adata.uns['spatial']['moran'][nn_graph].loc[hvgs, 'I']
adata.var.loc[hvgs, "moran"] = hvgs_moransI

adata.var.loc[hvgs, ["symbol", "moran"]].sort_values(by='moran')

top_moran = adata.var[adata.var['highly_variable']]['moran'] \
                    .dropna() \
                    .sort_values(ascending=False) \
                    .head(9) \
                    .index.tolist()

# CHECKPOINT: top_moran
checkpoint.add("top_moran.txt", top_moran)

neg_moran = adata.var[adata.var['highly_variable']]['moran'] \
                    .dropna() \
                    .sort_values(ascending=True) \
                    .head(9) \
                    .index.tolist()

# CHECKPOINT: neg_moran
checkpoint.add("neg_moran.txt", neg_moran)

vp.spatial.moran(adata, feature = neg_moran, graph_name = nn_graph, dim = 'var', permutations = 200)

raw_counts = adata.layers['counts']

if issparse(raw_counts):
    # For sparse matrices, compute the mean across the rows (genes)
    mean_expression_across_spots = np.array(raw_counts.mean(axis=0)).flatten()  # mean along columns
else:
    # For dense matrices, you can directly use numpy's mean function
    mean_expression_across_spots = np.mean(raw_counts, axis=0)  # mean along columns

adata.var['mean_expression_across_spots'] = mean_expression_across_spots

adata.X = vp.utils.scale(adata.X, center=True)

sc.tl.pca(adata, use_highly_variable=True, n_comps=30, random_state=1337)

# CHECKPOINT: pca.mtx
checkpoint.add("pca_embedding.mtx", adata.obsm["X_pca"])
checkpoint.add("pca_vec.mtx", adata[:, hvgs].varm["PCs"])

adata.X = adata.layers['logcounts'].copy()

sc.pp.neighbors(
    adata,
    n_pcs=10,
    use_rep='X_pca',
    knn=True,
    n_neighbors=10,
    key_added='knn'
)

dist = adata.obsp['knn_distances'].copy()
dist.data = 1 / dist.data

# row normalize the matrix, this makes the matrix dense.
dist /= dist.sum(axis=1)

# convert dist back to sparse matrix
adata.obsp["knn_weights"] = csr_matrix(dist)

del dist

knn_graph = "knn_weights"

adata.obsp["knn_connectivities"] = (adata.obsp[knn_graph] > 0).astype(int)
vp.spatial.set_default_graph(adata, knn_graph)
vp.spatial.to_spatial_weights(adata, graph_name=knn_graph)

morans = vp.spatial.moran(adata, hvgs, dim = 'var', graph_name=knn_graph)

knn_moran_df = pd.DataFrame(adata.uns['spatial']['moran'][knn_graph])
adata.var['moran_knn'] = adata.var.index.map(knn_moran_df['I'])

top_moran2 = adata.var[adata.var['highly_variable']]['moran_knn'] \
                    .dropna() \
                    .sort_values(ascending=False) \
                    .head(9) \
                    .index.tolist()

# CHECKPOINT: top_moran2
checkpoint.add("top_moran2.txt", top_moran2)

