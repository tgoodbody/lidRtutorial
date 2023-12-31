#' ---
#' title: "LAScatalog processing engine"
#' ---
#' 
#' ## Relevant resources:
#' 
#' -   [Code](https://github.com/tgoodbody/lidRtutorial/blob/master/R/08_engine2.R)
#' -   [lidRbook section: Engine](https://r-lidar.github.io/lidRbook/engine2.html)
#' -   [lidRbook section: Thinking outside the box](https://r-lidar.github.io/lidRbook/outbox.html)
#' 
#' ## Overview
#' 
#' This code showcases the LASCATALOG PROCESSING ENGINE, which efficiently applies various functions to LiDAR catalogs in parallel. It introduces the `catalog_map()` function for processing LiDAR data in a catalog. The code includes routines to detect trees and calculate metrics on the LiDAR catalog.
#' 
#' ## Environment
#' 

# Clear environment
rm(list = ls(globalenv()))

# Load packages
library(lidR)
library(terra)
library(future)

#' 
#' ## Basic Usage
#' 
#' In this section, we will cover the basic usage of the `lidR` package, including reading LiDAR data, visualization, and inspecting metadata.
#' 
#' ### Basic Usage of `lidR` Package
#' 
#' This section introduces the basic usage of the `lidR` package for reading and visualizing LiDAR data, as well as inspecting metadata.
#' 
#' ### Reading and Visualizing LiDAR Data
#' 
#' We start by reading a LAS catalog and inspecting one of its LAS files.
#' 

# Read a LAS catalog
ctg <- readLAScatalog(folder = "data/Farm_A/")

# Inspect the first LAS file in the catalog
las_file <- ctg$filename[1]
las <- readLAS(las_file)
las

#' 
#' ### Visualizing LiDAR Data
#' 
#' We visualize the LiDAR data from the selected LAS file using a 3D plot.
#' 

# Visualize the LiDAR data in 3D
plot(las, bg = "white")

#' 
#' ### `catalog_map()` Function
#' 
#' This section demonstrates the use of the `catalog_map()` function for efficient processing of LiDAR data within a LAS catalog.
#' 
#' ### Problem Statement
#' 
#' We start by addressing a common problem - how can we apply operations to `LAS` data in a catalog?
#' 

# Read a LAS file from the catalog and filter surface points
las_file <- ctg$filename[16]
las <- readLAS(files = las_file, filter = "-drop_withheld -drop_z_below 0 -drop_z_above 40")
surflas <- filter_surfacepoints(las = las, res = 1)

#' 
#' ### Visualizing LiDAR Data
#' 
#' We visualize the selected LiDAR data, including both the original data and the surface points.
#' 

# Visualize the LiDAR data with a default color palette
plot(las, bg = "white")

# Visualize the surface points using a default color palette
plot(surflas, bg = "white")

#' 
#' ### Calculating Rumple Index
#' 
#' We calculate the rumple index using the `pixel_metrics()` function.
#' 

# Generate Area-based metrics
ri <- pixel_metrics(las = las, ~rumple_index(X,Y,Z), res = 10)
plot(ri)

#' 
#' ### Solution: `LAScatalog` Processing Engine
#' 
#' This section introduces the `LAScatalog` processing engine, a powerful tool for efficient processing of LAS data within a catalog.
#' 
#' ### Basic Usage of the `catalog_map()` Function
#' 
#' We demonstrate the basic usage of the `catalog_map()` function with a simple user-defined function.
#' 

# User-defined function for processing chunks
routine <- function(las){

  # Perform computation
  output <- pixel_metrics(las = las, func = ~max(Z), res = 20)

  return(output)
}

# Initialize parallel processing
plan(multisession)

# Specify catalog options
opt_filter(ctg) <- "-drop_withheld"

# Apply routine to catalog
out <- catalog_map(ctg = ctg, FUN = routine)

print(out)

#' 
#' ### User-Defined Functions for Processing
#' 
#' We demonstrate the use of user-defined functions to process LiDAR data within a catalog.
#' 

# User-defined function for rumple index calculation
routine_rumple <- function(las, res1 = 10, res2 = 1){
  
  # filter surface points and create rumple index
  las <- filter_surfacepoints(las = las, res = res2)
  output  <- pixel_metrics(las = las, ~rumple_index(X,Y,Z), res1)
  
  return(output)
}

# Set catalog options
opt_select(ctg) <- "xyz"
opt_filter(ctg) <- "-drop_withheld -drop_z_below 0 -drop_z_above 40"
opt_chunk_buffer(ctg) <- 0
opt_chunk_size(ctg) <- 0

# Specify options for merging
options <- list(alignment = 10)

# Apply the user-defined function to the catalog
ri <- catalog_map(ctg = ctg, FUN = routine_rumple, res1 = 10, res2 = 0.5, .options = options)

# Plot the output
plot(ri, col = height.colors(50))

#' 
#' ::: callout-note
#' ## Thinking outside the box
#' 
#' The LAScatalog engine is versatile! The functions that can be applied to LiDAR data are infinite - leverage the flexibility of `lidR` and create software that pushes the boundaries of research in forest inventory and management!
#' :::
#' 
#' ## Exercises
#' 
#' #### E1.
#' 
#' Implement Noise Filtering
#' 
#' -   Explain the purpose of the `filter_noise()` function.
#' -   Create a user-defined function to apply noise filtering using the `catalog_map()` function.
#' -   Make sure to consider buffered points when using lidR's `filter_*` functions.
#' 
#' 
#' ## Conclusion
#' 
#' This concludes the tutorial on using the `catalog_map()` function in the `lidR` package to efficiently process LAS data within a catalog.
#' 

# Instructions for cleaning up any existing .lax files
# (Note: Please replace 'path' with the appropriate path)
path <- "data/Farm_A/"
file_list <- list.files(path)
delete_lax <- file_list[grep("\\.lax$", file_list)]
file.remove(file.path(path, delete_lax))

