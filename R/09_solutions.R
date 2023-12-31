#' ---
#' title: "Excercise Solutions"
#' ---
#' 
#' ## Resources
#' 
#' -   [Code](https://github.com/tgoodbody/lidRtutorial/tree/main/R/09_solutions.R)
#' 
#' ## 1-LAS
#' 

# Load packages
library(lidR)
library(sf)
library(terra)

#' 
#' #### E1.
#' 
#' What are withheld points? Where are they in our pointcloud?
#' 

# According to ASPRS LAS specification http://www.asprs.org/wp-content/uploads/2019/07/LAS_1_4_r15.pdf page 18 "a point # that should not be included in processing (synonymous with Deleted)"

# They are on the edges. It looks like they correspond to a buffer. LAStools makes use of the withheld bit to flag some # points. Without more information on former processing step it is hard to say.

#' 
#' #### E2.
#' 
#' Read the file dropping the withheld points.
#' 

las <- readLAS("data/MixedEucaNat_normalized.laz", filter = "-drop_withheld")
plot(las)

#' 
#' #### E3.
#' 
#' The withheld points seem to be legitimate points that we want to keep.
#' 
#' Try to load the file including the withheld points but get rid of the warning (without using `suppressWarnings()`). Hint: Check available `-set_withheld` filters in `readLAS(filter = "-h")`
#' 

las <- readLAS("data/MixedEucaNat_normalized.laz", filter = "-set_withheld_flag 0")
plot(las, color = "Withheld_flag")

#' 
#' #### E4.
#' 
#' Load only the ground points and plot the point-cloud coloured by the return number of the point. Do it loading the strict minimal amount of memory (4.7 Mb). Hint: use `?lidR::readLAS` and see what `select` options might help.
#' 

las <- readLAS("data/MixedEucaNat_normalized.laz", filter = "-keep_class 2 -set_withheld_flag 0", select = "r")
plot(las, color = "ReturnNumber", legend = T)
format(object.size(las), "Mb")

#' 
#' ## 2-ROI
#' 

plots <- st_read("data/shapefiles/MixedEucaNatPlot.shp")
plot(las@header, map = FALSE)
plot(plots, add = TRUE)

#' 
#' #### E1.
#' 
#' Clip the 5 plots with a radius of 11.3 m,
#' 

inventory <- clip_roi(las, plots, radius = 11.3)
plot(inventory[[2]])

#' 
#' #### E2.
#' 
#' Clip a transect from A `c(203850, 7358950)` to B `c(203950, 7959000)`.
#' 

tr <- clip_transect(las, c(203850, 7358950), c(203950, 7359000), width = 5)
plot(tr, axis = T)

#' 
#' #### E3.
#' 
#' Clip a transect from A `c(203850, 7358950)` to B `c(203950, 7959000)` but reorient it so it is no longer on the XY diagonal. Hint = `?clip_transect`
#' 

## ptr <- clip_transect(las, c(203850, 7358950), c(203950, 7359000), width = 5, xz = TRUE)
## plot(tr, axis = T)
## plot(ptr, axis = T)
## plot(ptr$X, ptr$Z, cex = 0.25, pch = 19, asp = 1)

#' 
#' ## 3-ABA
#' 

las <- readLAS("data/MixedEucaNat_normalized.laz", select = "*",  filter = "-set_withheld_flag 0")

#' 
#' #### E1.
#' 
#' Assuming that biomass is estimated using the equation `B = 0.5 * mean Z + 0.9 * 90th percentile of Z` applied on first returns only, map the biomass.
#' 

B <- pixel_metrics(las, ~0.5*mean(Z) + 0.9*quantile(Z, probs = 0.9), 10, filter = ~ReturnNumber == 1L)
plot(B, col = height.colors(50))

B <- pixel_metrics(las, .stdmetrics_z, 10)
B <- 0.5*B[["zmean"]] + 0.9*B[["zq90"]]
plot(B, col = height.colors(50))

pixel_metrics(las, ~as.list(quantile(Z), 10))

#' 
#' #### E2.
#' 
#' Map the density of ground returns at a 5 m resolution with `pixel_metrics(filter = ~Classification == LASGROUND)`.
#' 

GND <- pixel_metrics(las, ~length(Z)/25, res = 5, filter = ~Classification == LASGROUND)
plot(GND, col = heat.colors(50))

#' 
#' #### E3.
#' 
#' Map pixels that are flat (planar) using `stdshapemetrics`. These could indicate potential roads.
#' 

m <- pixel_metrics(las, .stdshapemetrics, res = 3)
plot(m[["planarity"]], col = heat.colors(50))
flat <- m[["planarity"]] > 0.85
plot(flat)

#' 
#' ## 5-DTM
#' 
#' #### E1.
#' 
#' Plot and compare these two normalized point-clouds. Why do they look different? Fix that. Hint: filter.
#' 
#' Some non ground points are below 0. It can be slightly low noise point not classified as ground by the data provider. This low points not being numerous and dark blue we hardly see them
#' 

las1 <- readLAS("data/MixedEucaNat.laz", filter = "-set_withheld_flag 0")
nlas1 <- normalize_height(las1, tin())
nlas2 <- readLAS("data/MixedEucaNat_normalized.laz", filter = "-set_withheld_flag 0")
plot(nlas1)
plot(nlas2)

nlas1 <- filter_poi(nlas1, Z > -0.1)
plot(nlas1)

#' 
#' #### E2.
#' 
#' Clip a plot somewhere in `MixedEucaNat.laz` (the non-normalized file).
#' 

circ <- clip_circle(las, 203930, 7359000, 25)
plot(circ)

#' 
#' #### E3.
#' 
#' Compute a DTM for this plot. Which method are you choosing and why?
#' 

dtm <- rasterize_terrain(circ, 0.5, kriging())
plot_dtm3d(dtm)

#' 
#' #### E4.
#' 
#' Compute a DSM (digital surface model). Hint: Look back to how you made a CHM.
#' 

dsm <- rasterize_canopy(circ, 1, p2r(0.1))
plot(dsm, col = height.colors(50))

#' 
#' #### E5.
#' 
#' Normalize the plot.
#' 

ncirc <- circ - dtm
plot(ncirc)

#' 
#' #### E6.
#' 
#' Compute a CHM.
#' 

chm <- rasterize_canopy(ncirc, 1, p2r(0.1))
plot(chm, col = height.colors(50))

#' 
#' #### E7.
#' 
#' Estimate some metrics of interest in this plot with `cloud_metrics()`.
#' 

metrics <- cloud_metrics(ncirc, .stdmetrics_z)
metrics

#' 
#' ## 6-ITS
#' 
#' Using:
#' 

las <- readLAS(files = "data/MixedEucaNat_normalized.laz", filter = "-set_withheld_flag 0")
col1 <- height.colors(50)

#' 
#' #### E1.
#' 
#' Find and count the trees.
#' 

ttops <- locate_trees(las = las, algorithm = lmf(ws = 3, hmin = 5))
x <- plot(las)
add_treetops3d(x = x, ttops = ttops)

#' 
#' #### E2.
#' 
#' Compute and map the density of trees with a 10 m resolution.
#' 

r <- terra::rast(x = ttops)
terra::res(r) <- 10
r <- terra::rasterize(x = ttops, y = r, "treeID", fun = 'count')
plot(r, col = viridis::viridis(20))

#' 
#' #### E3.
#' 
#' Segment the trees.
#' 

chm <- rasterize_canopy(las = las, res = 0.5, algorithm = p2r(subcircle = 0.15))
plot(chm, col = col1)
ttops <- locate_trees(las = chm, algorithm = lmf(ws = 2.5))
las <- segment_trees(las = las, dalponte2016(chm = chm, treetops = ttops))

plot(las, color = "treeID")

#' 
#' #### E4.
#' 
#' Assuming that a value of interest of a tree can be estimated using the crown area and the mean `Z` of the points with the formula `2.5 * area + 3 * mean Z`. Estimate the value of interest of each tree.
#' 

value_of_interest <- function(x,y,z)
{
  m <- stdtreemetrics(x,y,z)
  avgz <- mean(z)
  v <- 2.5*m$convhull_area + 3 * avgz
  return(list(V = v))
}

V <- crown_metrics(las = las, func = ~value_of_interest(X,Y,Z))
plot(x = V["V"])

#' 
#' #### E5.
#' 
#' Map the total biomass at a resolution of 10 m. The output is a mixed of ABA and ITS
#' 

Vtot <- rasterize(V, r, "V", fun = "sum")
plot(Vtot, col = viridis::viridis(20))

#' 
#' ## 7-LASCTALOG
#' 
#' This exercise is complex because it involves options not yet described. Be sure to use the [lidRbook](https://r-lidar.github.io/lidRbook/index.html) and [package documentation](https://cran.r-project.org/web/packages/lidR/lidR.pdf).
#' 
#' Using:
#' 

ctg <- readLAScatalog(folder = "data/Farm_A/")

#' 
#' #### E1.
#' 
#' Generate a raster of point density for the provided catalog. Hint: Look through the documentation for a function that will do this!
#' 

ctg <- readLAScatalog("data/Farm_A/", filter = "-drop_withheld -drop_z_below 0 -drop_z_above 40")
D1 <- rasterize_density(las = ctg, res = 4)
plot(D1, col = heat.colors(50))

#' 
#' #### E2.
#' 
#' Modify the catalog to have a point density of 10 pts/m2 using the `decimate_points()` function. If you get an error make sure to read the documentation for `decimate_points()` and try: using `opt_output_file()` to [write files to a temporary directory](https://r-lidar.github.io/lidRbook/engine.html#engine-dtm-ondisk).
#' 

## newctg <- decimate_points(las = ctg, algorithm = homogenize(density = 10, res = 5))
## #>  Error: This function requires that the LAScatalog provides an output file template.

#' 

opt_filter(ctg) <- "-drop_withheld"
opt_output_files(ctg) <- paste0(tempdir(), "/{ORIGINALFILENAME}")
newctg <- decimate_points(las = ctg, algorithm = homogenize(density = 10, res = 5))

#' 
#' #### E3.
#' 
#' Generate a raster of point density for this new decimated dataset.
#' 

opt_output_files(newctg) <- ""
D2 <- rasterize_density(las = newctg, res = 4)
plot(D2, col = heat.colors(50))

#' 
#' #### E4.
#' 
#' Read the whole decimated catalog as a single las file. The catalog isn't very big - not recommended for larger data sets!
#' 

las <- readLAS(newctg)
plot(las)

#' 
#' #### E5.
#' 
#' Read documentation for the `catalog_retile()` function and merge the dataset into larger tiles. Use `ctg` metadata to align new chunks to the lower left corner of the old ones. Hint: Visualize the chunks and use `opt_chunk_*` options.
#' 

opt_chunk_size(ctg) <- 280
opt_chunk_buffer(ctg) <- 0
opt_chunk_alignment(ctg) <- c(min(ctg$Min.X), min(ctg$Min.Y))
plot(ctg, chunk = T)

opt_output_files(ctg) <- "{tempdir()}/PRJ_A_{XLEFT}_{YBOTTOM}"
newctg <- catalog_retile(ctg = ctg)
plot(newctg)

#' 
#' ## 8-ENGINE
#' 
#' #### E1.
#' 
#' The following is a simple (and a bit naive) function to remove high noise points. - Explain what this function does - Create a user-defined function to apply using `catalog_map()` - Hint: Dont forget about buffered points... remember `lidR::filter_*` functions.
#' 

filter_noise <- function(las, sensitivity)
{
  p95 <- pixel_metrics(las, ~quantile(Z, probs = 0.95), 10)
  las <- merge_spatial(las, p95, "p95")
  las <- filter_poi(las, Z < 1+p95*sensitivity, Z > -0.5)
  las$p95 <- NULL
  return(las)
}

filter_noise_collection = function(las, sensitivity)
{
  las <- filter_noise(las, sensitivity)
  las <- filter_poi(las, buffer == 0L)
  return(las)
}

ctg = readLAScatalog("data/Farm_A/")
opt_select(ctg) <- "*"
opt_filter(ctg) <- "-drop_withheld"
opt_output_files(ctg) <- "{tempdir()}/*"
opt_chunk_buffer(ctg) <- 20
opt_chunk_size(ctg) <- 0

output <- catalog_map(ctg, filter_noise_collection, sensitivity = 1.2)

las <- readLAS(output)
plot(las)

