setwd("/workspace/vignettes/visium_10x_spatial")

suppressMessages({
    library(Voyager)
    library(SpatialFeatureExperiment)
    library(SingleCellExperiment)
    library(ggplot2)
    library(scater)
    library(scuttle)
    library(scran)
    library(stringr)
    library(patchwork)
    library(bluster)
    library(rjson)
    library(EBImage)
    library(terra)
    library(rlang)
    library(sf)
    library(rmapshaper)
    library(dplyr)
    library(BiocParallel)
    library(BiocNeighbors)
    library(reticulate)
})

PY_PATH <- "/home/rstudio/miniconda/envs/voyagerpy/bin/python"   # system("which python", intern = TRUE)
use_python(PY_PATH)
py_config()

# system("pip3 install gget")
gget <- import("gget")


# root_path <- normalizePath(Sys.getenv("RPATH"))
utils_path <- file.path("/workspace/utils.R")   # file.path(root_path, "utils.R")
source(utils_path)
init.dirs()



(sfe <- read10xVisiumSFE(samples = ".", type = "sparse", data = "raw"))

is_mt <- str_detect(rowData(sfe)$symbol, "^mt-")

sfe <- addPerCellQCMetrics(sfe, subsets = list(mito = is_mt))

img <- readImage("outs/spatial/tissue_hires_image.png")

img2 <- img
colorMode(img2) <- Grayscale

mask <- img2[,,3] < 0.87

kern <- makeBrush(3, shape='disc')
mask_open <- opening(mask, kern)

mask_close <- closing(mask_open, kern)

mask_label <- bwlabel(mask_close)
fts <- computeFeatures.shape(mask_label)

summary(fts[,"s.area"])

max_ind <- which.max(fts[,"s.area"])
inds <- which(as.array(mask_label) == max_ind, arr.ind = TRUE)

row_inds <- c(seq_len(min(inds[,1])-1), seq(max(inds[,1])+1, nrow(mask_label), by = 1))
col_inds <- c(seq_len(min(inds[,2])-1), seq(max(inds[,2])+1, nrow(mask_label), by = 1))
mask_label[row_inds, ] <- 0
mask_label[,col_inds] <- 0

unique(as.vector(mask_label))

fts2 <- fts[unique(as.vector(mask_label))[-1],]
fts2 <- fts2[order(fts2[,"s.area"], decreasing = TRUE),]

plot(fts2[,1][-1], type = "l", ylab = "Area")

mask_label[mask_label %in% c(797, as.numeric(rownames(fts2)[fts2[,1] < 100]))] <- 0

mask_label <- fillHull(mask_label)

raster2polygon <- function(seg, keep = 0.2) {
    r <- rast(as.array(seg), extent = ext(0, nrow(seg), 0, ncol(seg))) |> 
        trans() |> flip()
    r[r < 1] <- NA
    contours <- st_as_sf(as.polygons(r, dissolve = TRUE))
    simplified <- ms_simplify(contours, keep = keep)
    list(full = contours,
         simplified = simplified)
}

tb <- raster2polygon(mask_label)

scale_factors <- fromJSON(file = "outs/spatial/scalefactors_json.json")

tb$simplified$geometry <- tb$simplified$geometry / scale_factors$tissue_hires_scalef

tissueBoundary(sfe) <- tb$simplified

sfe <- SpatialFeatureExperiment::transpose(sfe)

sfe$int_tissue <- annotPred(sfe, colGeometryName = "spotPoly", 
                            annotGeometryName = "tissueBoundary",
                            pred = st_intersects)
sfe$cov_tissue <- annotPred(sfe, colGeometryName = "spotPoly", 
                            annotGeometryName = "tissueBoundary",
                            pred = st_covered_by)

# CHECKPOINT: intersection_and_contain_df
intersection_and_contain_df <- colData(sfe)[, c(
	"int_tissue", 
	"cov_tissue"
)] %>% as.data.frame
checkpoint.csv("intersection_and_contain_df.csv", intersection_and_contain_df)
# checkpoint.csv("intersection_and_contain_df.csv", intersection_and_contain_df, sync=TRUE)


sfe$diff_sr <- case_when(sfe$in_tissue == sfe$int_tissue ~ "same",
                         sfe$in_tissue & !sfe$int_tissue ~ "Space Ranger",
                         sfe$int_tissue & !sfe$in_tissue ~ "segmentation") |> 
    factor(levels = c("Space Ranger", "same", "segmentation"))
plotSpatialFeature(sfe, "diff_sr", 
                   annotGeometryName = "tissueBoundary", 
                   annot_fixed = list(fill = NA, size = 0.5, color = "black")) +
    scale_fill_brewer(type = "div", palette = 4)

sfe$diff_int_cov <- sfe$int_tissue != sfe$cov_tissue

spot_ints <- annotOp(sfe, colGeometryName = "spotPoly", 
                     annotGeometryName = "tissueBoundary", op = st_intersection)
sfe$pct_tissue <- st_area(spot_ints) / st_area(spotPoly(sfe)) * 100

sfe_tissue <- sfe[,sfe$int_tissue]

colGraph(sfe_tissue, "visium") <- findVisiumGraph(sfe_tissue)

# CHECKPOINT: visium_graph.mtx
# checkpoint.mtx("visium_graph.mtx", colGraph(sfe_tissue, "visium"))
# checkpoint.mtx("visium_graph.mtx", colGraph(sfe_tissue, "visium"), sync=TRUE)


qc_features <- c("sum", "detected", "subsets_mito_percent")

sfe_tissue <- colDataUnivariate(sfe_tissue, "moran.mc", qc_features, nsim = 200)

# CHECKPOINT: qc_moran
qc_moran <- colFeatureData(sfe_tissue)[qc_features, 'moran.mc_statistic_sample01'] %>% as.data.frame
checkpoint.csv("qc_moran.csv", qc_moran)

sfe_tissue <- colDataUnivariate(sfe_tissue, "sp.correlogram", qc_features,
                                order = 8)

correlogram_list <- colFeatureData(sfe_tissue)[qc_features, 'sp.correlogram_I_sample01']

I_values_list <- list()

for (i in seq_along(correlogram_list)) {
  I_values <- correlogram_list[[i]][, "I"]
  I_values_list[[i]] <- I_values
}

I_values_df <- do.call(rbind, I_values_list)

# Set appropriate row names and column names
row.names(I_values_df) <- qc_features
colnames(I_values_df) <- seq_len(ncol(I_values_df))

# CHECKPOINT: qc_correlogram
checkpoint.csv("qc_correlogram.csv", I_values_df)

sfe_tissue <- colDataUnivariate(sfe_tissue, "localmoran", qc_features)

# Extract the list of data frames
local_moran_list <- localResults(sfe_tissue)[["localmoran"]]

# Extract the 'Ii' column for each feature in qc_features
Ii_values_list <- lapply(qc_features, function(feature) {
  if (feature %in% names(local_moran_list)) {
    return(local_moran_list[[feature]][, "Ii"])
  } else {
    return(NULL)
  }
})

# Combine the Ii values into a data frame
localmoran_values_df <- do.call(cbind, Ii_values_list)

# Set row and column names
row.names(localmoran_values_df) <- row.names(local_moran_list[[1]])
colnames(localmoran_values_df) <- qc_features

# CHECKPOINT: qc_local_moran
checkpoint.csv("qc_local_moran.csv", localmoran_values_df)

sfe_tissue <- colDataUnivariate(sfe_tissue, "LOSH", qc_features)

# Extract the list of data frames
losh_list <- localResults(sfe_tissue)[["LOSH"]]

# Extract the 'Ii' column for each feature in qc_features
Ii_values_list <- lapply(qc_features, function(feature) {
  if (feature %in% names(losh_list)) {
    return(losh_list[[feature]][, "Hi"])
  } else {
    return(NULL)
  }
})

# Combine the Ii values into a data frame
losh_values_df <- do.call(cbind, Ii_values_list)

# Set row and column names
row.names(losh_values_df) <- row.names(losh_list[[1]])
colnames(losh_values_df) <- qc_features

# Display the resulting data frame
print(losh_values_df)

# CHECKPOINT: qc_losh
checkpoint.csv("qc_losh.csv", losh_values_df)

sfe_tissue <- colDataUnivariate(sfe_tissue, "moran.plot", qc_features)

sfe_tissue <- logNormCounts(sfe_tissue)

# CHECKPOINT: logcounts.mtx
checkpoint.mtx("logcounts.mtx", logcounts(sfe_tissue))

dec <- modelGeneVar(sfe_tissue)

# CHECKPOINT: gene_var
checkpoint.txt("gene_var.txt", rownames(dec))

hvgs <- getTopHVGs(dec, n = 2000)

# CHECKPOINT AND SYNC: hvgs
checkpoint.txt("hvgs.txt", hvgs)
checkpoint.txt("hvgs.txt", hvgs, sync=TRUE)

sfe_tissue <- runMoransI(sfe_tissue, features = hvgs, BPPARAM = MulticoreParam(2))

# CHECKPOINT: hvgs_moran
hvgs_moran <- rowData(sfe_tissue)[hvgs, 'moran_sample01'] %>% as.data.frame
rownames(hvgs_moran) <- hvgs
checkpoint.csv("hvgs_moran.csv", hvgs_moran)

top_moran <- rownames(sfe_tissue)[order(rowData(sfe_tissue)$moran_sample01, 
                                        decreasing = TRUE)[1:9]]

# CHECKPOINT: top_moran
checkpoint.txt("top_moran.txt", top_moran)

# gget_info <- gget$info(top_moran)
# 
# rownames(gget_info) <- gget_info$primary_gene_name
# select(gget_info, ncbi_description)

neg_moran <- rownames(sfe_tissue)[order(rowData(sfe_tissue)$moran_sample01, 
                                        decreasing = FALSE)[1:9]]

# CHECKPOINT: neg_moran
checkpoint.txt("neg_moran.txt", neg_moran)

# gget_info_neg <- gget$info(neg_moran)
# 
# rownames(gget_info_neg) <- gget_info_neg$primary_gene_name
# select(gget_info_neg, ncbi_description)

sfe_tissue <- runUnivariate(sfe_tissue, "moran.mc", neg_moran, 
                            colGraphName = "visium", nsim = 200, alternative = "less")

rowData(sfe_tissue)[neg_moran, c("moran_sample01", "moran.mc_p.value_sample01")]

sfe_tissue <- addPerFeatureQCMetrics(sfe_tissue)
names(rowData(sfe_tissue))

plotRowData(sfe_tissue, x = "mean", y = "moran_sample01") +
    scale_x_log10() +
    annotation_logticks(sides = "b") +
    geom_density2d()

sfe_tissue <- runPCA(sfe_tissue, ncomponents = 30, subset_row = hvgs,
                     scale = TRUE) # scale as in Seurat

# CHECKPOINT: pca.mtx
pca.mat <- reducedDim(sfe_tissue, "PCA") %>% Matrix::Matrix(sparse=TRUE)
. <- checkpoint.mtx("pca_embedding.mtx", pca.mat)
pca.rot <- attr(reducedDim(sfe_tissue, "PCA"), "rotation") %>% Matrix::Matrix(sparse=TRUE)
. <- checkpoint.mtx("pca_vec.mtx", pca.rot)

foo <- findKNN(reducedDim(sfe_tissue, "PCA")[,1:10], k=10, BNPARAM=AnnoyParam())
# Split by row
foo_nb <- asplit(foo$index, 1)
dmat <- 1/foo$distance
# Row normalize the weights
dmat <- sweep(dmat, 1, rowSums(dmat), FUN = "/")
glist <- asplit(dmat, 1)
# Sort based on index
ord <- lapply(foo_nb, order)
foo_nb <- lapply(seq_along(foo_nb), function(i) foo_nb[[i]][ord[[i]]])
class(foo_nb) <- "nb"
glist <- lapply(seq_along(glist), function(i) glist[[i]][ord[[i]]])

listw <- list(style = "W",
              neighbours = foo_nb,
              weights = glist)
class(listw) <- "listw"
attr(listw, "region.id") <- colnames(sfe_tissue)

colGraph(sfe_tissue, "knn10") <- listw

sfe_tissue <- runMoransI(sfe_tissue, features = hvgs, BPPARAM = MulticoreParam(2),
                         colGraphName = "knn10", name = "moran_ns")

top_moran2 <- rownames(sfe_tissue)[order(rowData(sfe_tissue)$moran_ns_sample01, 
                                         decreasing = TRUE)[1:9]]

# CHECKPOINT: top_moran2
checkpoint.txt("top_moran2.txt", top_moran2)

# gget_info2 <- gget$info(top_moran2)
# 
# rownames(gget_info2) <- gget_info2$primary_gene_name
# select(gget_info2, ncbi_description)