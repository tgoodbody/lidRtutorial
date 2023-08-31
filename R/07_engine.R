#' ---
#' title: "LAScatalog"
#' ---
#' 
#' ## Relevant resources:
#' 
#' -   [Code](https://github.com/tgoodbody/lidRtutorial/blob/main/R/07_engine.R)
#' -   [lidRbook section](https://r-lidar.github.io/lidRbook/engine.html)
#' 
#' ## Overview
#' 
#' This code performs various operations on LiDAR data using LAScatalog functionality. We visualize and inspect the data, validate the files, clip the data based on specific coordinates, generate a Canopy Height Model (CHM), compute above ground biomass, detect treetops, specify processing options, and use parallel computing.
#' 
#' ## Environment
#' 

# Clear environment
rm(list = ls(globalenv()))

# Load packages
library(lidR)
library(sf)

#' 
#' ## Basic Usage
#' 
#' In this section, we will cover the basic usage of the `lidR` package, including reading LiDAR data, visualization, and inspecting metadata.
#' 
#' ### Read catalog from directory of files
#' 
#' We begin by creating a LAS catalog (`ctg`) from a folder containing multiple LAS files using the `readLAScatalog()` function.
#' 

# Read catalog and drop withheld
ctg <- readLAScatalog(folder = "data/Farm_A/",filter = "-drop_withheld")

#' 
#' ### Inspect catalog
#' 
#' We can inspect the contents of the catalog using standard R functions.
#' 

ctg

#' 
#' ### Visualize catalog
#' 
#' We visualize the catalog, showing the spatial coverage of the LiDAR data *header* extents. The map can be interactive if we use `map = TRUE`. Try clicking on a tile to see its header information.
#' 

plot(ctg)

# Interactive
plot(ctg, map = TRUE)

#' 
#' ### Check coordinate system and extent info
#' 
#' We examine the coordinate system and extent information of the catalog.
#' 

# coordinate system
crs(ctg)
projection(ctg)
st_crs(ctg)

# spatial extents
extent(ctg)
bbox(ctg)
st_bbox(ctg)

#' 
#' ## Validate files in catalog
#' 
#' We validate the LAS files within the catalog using the `las_check` function. It works the same way as it would on a regular `LAS` file.
#' 

# Perform a check on the catalog
las_check(las = ctg)

#' 
#' ## File indexing
#' 
#' We explore indexing of `LAScatalog` objects for efficient processing.
#' 
#' ::: callout-note
#' ## Indexing
#' 
#' **The `lidR` policy has always been: use `LAStools` and `lasindex` for spatial indexing**. If you really don't want, or can't use `LAStools`, then there is a hidden function in `lidR` that users can use (`lidR:::catalog_laxindex()`).
#' :::
#' 

# Instructions for cleaning up any existing .lax files
# (Note: Please replace 'path' with the appropriate path)
path <- "data/Farm_A/"
file_list <- list.files(path)
delete_lax <- file_list[grep("\\.lax$", file_list)]
file.remove(file.path(path, delete_lax))

# check if files have .lax
is.indexed(ctg)

# generate index files
lidR:::catalog_laxindex(ctg)

# check if files have .lax
is.indexed(ctg)

#' 
#' ## Generate CHM
#' 
#' We create a CHM by rasterizing the point cloud data from the catalog.
#' 

# Generate CHM
chm <- rasterize_canopy(las = ctg, res = 0.5, algorithm = p2r(subcircle = 0.15))
plot(chm, col = height.colors(50))

#' 
#' We encounter issues and warnings while generating the CHM. Let's figure out how to fix the warnings and get decent outputs.
#' 

# Check for warnings
warnings()

#' 
#' ## Catalog processing options
#' 
#' We explore and manipulate catalog options.
#' 

# Setting options and re-rasterizing the CHM
opt_filter(ctg) <- "-drop_z_below 0 -drop_z_above 40"
opt_select(ctg) <- "xyz"
chm <- rasterize_canopy(las = ctg, res = 0.5, algorithm = p2r(subcircle = 0.15))
plot(chm, col = height.colors(50))

#' 
#' ## Area-based approach on catalog
#' 
#' In this section, we generate Above Ground Biomass (ABA) estimates using the LAScatalog.
#' 
#' ### Generate ABA output and visualize
#' 
#' We calculate ABA using the `pixel_metrics` function and visualize the results.
#' 

# Generate area-based metrics
model <- pixel_metrics(las = ctg, func = ~max(Z), res = 20)
plot(model, col = height.colors(50))

#' 
#' ## First returns only
#' 
#' We adjust the catalog options to calculate ABA based on first returns only.
#' 

opt_filter(ctg) <- "-drop_z_below 0 -drop_z_above 40 -keep_first"
model <- pixel_metrics(las = ctg, func = ~max(Z), res = 20)
plot(model, col = height.colors(50))

#' 
#' ## Clip a catalog
#' 
#' We clip the `LAS` data in the catalog using specified coordinate groups.
#' 

# Set coordinate groups
x <- c(207846, 208131, 208010, 207852, 207400)
y <- c(7357315, 7357537, 7357372, 7357548, 7357900)

# Visualize coordinate groups
plot(ctg)
points(x, y)

# Clip plots
rois <- clip_circle(las = ctg, xcenter = x, ycenter = y, radius = 30)

# Visualize the LiDAR data with a default color palette
plot(rois[[1]], bg = "white")

# Visualize the LiDAR data with a default color palette
plot(rois[[3]], bg = "white")

#' 
#' ## Validate clipped data
#' 
#' We validate the clipped LAS data using the `las_check` function.
#' 

las_check(rois[[1]])
las_check(rois[[3]])

#' 
#' ## Independent files (e.g. plots) as catalogs
#' 
#' We read an individual LAS file as a catalog and perform operations on it.
#' 

# Instructions for cleaning up any existing .lax files
# (Note: Please replace 'path' with the appropriate path)
path <- paste0(tempdir())
file_list <- list.files(path)
delete_tif <- file_list[grep("\\.tif$", file_list)]
delete_las <- file_list[grep("\\.laz$", file_list)]
file.remove(file.path(path, delete_tif))
file.remove(file.path(path, delete_las))

#' 

# Read single file as catalog
ctg_non_norm <- readLAScatalog(folder = "data/MixedEucaNat.laz")

# Set options for output files
opt_output_files(ctg_non_norm) <- paste0(tempdir(),"/{XCENTER}_{XCENTER}")

# Write file as .laz
opt_laz_compression(ctg_non_norm) <- TRUE

# Get random plot locations and clip
x <- runif(n = 4, min = ctg_non_norm$Min.X, max = ctg_non_norm$Max.X)
y <- runif(n = 4, min = ctg_non_norm$Min.Y, max = ctg_non_norm$Max.Y)
rois <- clip_circle(las = ctg_non_norm, xcenter = x, ycenter = y, radius = 10)

#' 

# Read catalog of plots
ctg_plots <- readLAScatalog(tempdir())

# Set independent files option
opt_independent_files(ctg_plots) <- TRUE
opt_output_files(ctg_plots) <- paste0(tempdir(),"/{XCENTER}_{XCENTER}")

# Generate plot-level terrain models
rasterize_terrain(las = ctg_plots, res = 1, algorithm = tin())

# Check files
path <- paste0(tempdir())
file_list <- list.files(path, full.names = TRUE)
file <- file_list[grep("\\.tif$", file_list)][[1]]

# plot dtm
plot(terra::rast(file))

#' 
#' ## ITD using `LAScatalog`
#' 
#' In this section, we explore Individual Tree Detection (ITD) using the `LAScatalog`. We first configure catalog options for ITD.
#' 

# Set catalog options
opt_filter(ctg) <- "-drop_withheld -drop_z_below 0 -drop_z_above 40"

#' 
#' ### Detect treetops and visualize
#' 
#' We detect treetops and visualize the results.
#' 

# Detect tree tops and plot
ttops <- locate_trees(las = ctg, algorithm = lmf(ws = 3, hmin = 5))
plot(chm, col = height.colors(50))
plot(ttops, add = TRUE, cex = 0.1, col = "black")

#' 
#' ### Specify catalog options
#' 
#' We specify additional catalog options for ITD.
#' 

# Specify more options
opt_select(ctg) <- "xyz"
opt_chunk_size(ctg) <- 300
opt_chunk_buffer(ctg) <- 10

# Detect treetops and plot
ttops <- locate_trees(las = ctg, algorithm = lmf(ws = 3, hmin = 5))
plot(chm, col = height.colors(50))
plot(ttops, add = TRUE, cex = 0.1, col = "black")

#' 
#' ### Parallel computing
#' 
#' In this section, we explore parallel computing using the `lidR` package.
#' 
#' ## Load `future` library
#' 
#' We load the `future` library to enable parallel processing.
#' 

library(future)

#' 
#' ## Specify catalog options
#' 
#' We specify catalog options for parallel processing.
#' 

# Specify options
opt_select(ctg) <- "xyz"
opt_chunk_size(ctg) <- 300
opt_chunk_buffer(ctg) <- 10

# Visualize and summarize the catalog chunks
plot(ctg, chunk = TRUE)
summary(ctg)

#' 
#' ## Single core processing
#' 
#' We perform tree detection using a single core.
#' 

# Process on single core
future::plan(sequential)

# Detect trees
ttops <- locate_trees(las = ctg, algorithm = lmf(ws = 3, hmin = 5))

#' 
#' ## Parallel processing
#' 
#' We perform tree detection using multiple cores in parallel.
#' 

# Process multi-core
future::plan(multisession)

# Detect trees
ttops <- locate_trees(las = ctg, algorithm = lmf(ws = 3, hmin = 5))

#' 
#' ## Revert to single core
#' 
#' We revert to single core processing using `future::plan(sequential)`.
#' 

# Back to single core
future::plan(sequential)

#' 
#' This concludes the tutorial on basic usage, catalog validation, indexing, CHM generation, ABA estimation, data clipping, ITD using catalog, and parallel computing using the `lidR` package in R.
#' 
#' ## Exercises and Questions
#' 
#' ::: callout-tip
#' This exercise is complex because it involves options not yet described. Be sure to use the [lidRbook](https://r-lidar.github.io/lidRbook/) and [package documentation](https://cran.r-project.org/web/packages/lidR/lidR.pdf).
#' :::
#' 
#' Using:

ctg <- readLAScatalog(folder = "data/Farm_A/")

#' #### E1.
#' 
#' Generate a raster of point density for the provided catalog. Hint: Look through the documentation for a function that will do this!
#' 
#' #### E2.
#' 
#' Modify the catalog to have a point density of 10 pts/m2 using the `decimate_points()` function. If you get an error make sure to read the documentation for `decimate_points()` and try: using `opt_output_file()` to write files to a temporary directory.
#' 
#' #### E3.
#' 
#' Generate a raster of point density for this new decimated dataset.
#' 
#' #### E4.
#' 
#' Read the whole decimated catalog as a single `las` file. The catalog isn't very big - not recommended for larger datasets!
#' 
#' #### E5.
#' 
#' Read documentation for the `catalog_retile()` function and merge the decimated catalog into larger tiles.
