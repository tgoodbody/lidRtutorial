#' ---
#' title: "Read/Plot/Query/Validate"
#' ---
#' 
#' ## Relevant Resources
#' 
#' -   [Code](https://github.com/tgoodbody/lidRtutorial/blob/main/R/01_read.R)
#' 
#' -   [lidRbook section](https://r-lidar.github.io/lidRbook/index.htm)
#' 
#' ## Overview
#' 
#' Welcome to this LiDAR processing tutorial using R and the `lidR` package! In this tutorial, you will learn how to *read*, *visualize*, *query*, and *validate* LiDAR data. We'll explore basic information about a LiDAR file including the *header* and *tabular data,* as well as *visualize* point clouds using different color schemes based on attributes. We'll use the `select` argument in `readLAS()` to load specific attributes and the `filter` argument to only load points of interest or apply transformations on-the-fly. We'll validate the LiDAR data using the `las_check` function on different data files to ensure data integrity.
#' 
#' Let's get started with processing LiDAR data efficiently using `lidR` and R! Happy learning!
#' 
#' ## Environment
#' 
#' We start by loading the necessary packages, clearing our current environment, and specifying that some warnings be turned off to make our outputs clearer. We will do this for each section in the tutorial.
#' 

# Clear environment
rm(list = ls(globalenv()))

# Load packages
library(lidR)
library(sf)

#' 
#' ## Basic Usage
#' 
#' ### Load and Inspect LiDAR Data
#' 
#' Load the LiDAR point cloud data from a LAS file using the `readLAS()` function. The data is stored in the `las` object. We can inspect the header information and attributes of the `las` object.
#' 

las <- readLAS(files = "data/MixedEucaNat_normalized.laz")

# Inspect header information
las@header

# Inspect attributes of the point cloud
las@data

# Check the file size of the loaded LiDAR data
format(object.size(las), "Mb")

#' 
#' ### Visualize LiDAR Data
#' 
#' We can visualize the LiDAR data using the `plot()` function. We have several options to control the colors in the plot, such as selecting specific attributes from the data to be used as colors.
#' 

# Visualize the LiDAR data with a default color palette
plot(las, bg = "white")

# Visualize using intensity values as colors
plot(las, color = "Intensity", bg = "white")

# Visualize using the classification attribute as colors
plot(las, color = "Classification", bg = "white")

# Visualize using the scan angle rank as colors with an axis and legend
plot(las, color = "ScanAngleRank", axis = TRUE, legend = TRUE, bg = "white")

#' 
#' ## Optimized Usage
#' 
#' ### Selecting Attributes of Interest
#' 
#' The `readLAS()` function allows us to select specific attributes to be loaded into memory. This is useful to reduce memory requirements when working with large LiDAR datasets.
#' 

# Load only the xyz coordinates (X, Y, Z) and ignore other attributes
las <- readLAS(files = "data/MixedEucaNat_normalized.laz", select = "xyz")

# Inspect the loaded attributes
las@data

# Check the memory size after loading only the selected attributes
format(object.size(las), "Mb")

#' 
#' ### Filtering Points of Interest
#' 
#' We can also load only a subset of the LiDAR points based on certain criteria using the `filter` argument in `readLAS()`.
#' 

# Load only the first return points
las <- readLAS(files = "data/MixedEucaNat_normalized.laz", filter = "-keep_first")

# Inspect the loaded points
las

# Check the memory size after loading only the filtered points
format(object.size(las), "Mb")

# Visualize the filtered points
plot(las, bg = "white")

#' 
#' ### Applying Transformation on-the-fly
#' 
#' The `filter` argument in `readLAS()` can be used to apply transformations on-the-fly during loading. This can be useful for tasks such as filtering specific classifications.
#' 

# Load and visualize with an applied filter
las <- readLAS(files = "data/MixedEucaNat.laz", filter = "-keep_class 2 -keep_class 1")

plot(las, bg = "white")

#' 
#' ### Filtering Points using `filter_poi()`
#' 
#' An alternative method for filtering points is using the `filter_poi()` function, which allows filtering based on attributes of points after the pointcloud is loaded.
#' 

# Filter points with Classification == 2
class_2 <- filter_poi(las = las, Classification == 2L)

# Combine queries to filter points with Classification == 1 and ReturnNumber == 1
first_returns <- filter_poi(las = las, Classification == 1L & ReturnNumber == 1L)

plot(class_2, bg = "white")

plot(first_returns, bg = "white")

#' 
#' ### `LAS` Objects Validation
#' 
#' The `lidR` package provides a function `las_check()` to validate LAS objects for common issues.
#' 

# Load and validate LAS data
las <- readLAS(files = "data/MixedEucaNat_normalized.laz")
las_check(las)

# Visualize corrupted LAS data
las <- readLAS(files = "data/example_corrupted.laz")

plot(las, bg = "white")

#' 

# Validate corrupted LAS data
las_check(las)

#' 
#' ## Exercises and Questions
#' 
#' Using:
#' 

las <- readLAS(files = "data/MixedEucaNat_normalized.laz")

#' 
#' #### E1.
#' 
#' What are withheld points? Where are they in our point cloud?
#' 
#' #### E2.
#' 
#' Read the file dropping withheld points.
#' 
#' #### E3.
#' 
#' The withheld points seem to be legitimate points that we want to keep. Try to load the file including the withheld points but get rid of the warning (without using `suppressWarnings()`). Hint: Check available `-set_withheld` filters using `readLAS(filter = "-h")`.
#' 
#' #### E4.
#' 
#' Load only the ground points and plot the point cloud colored by the return number of the point. Do it loading the strict minimal amount of memory (4.7 Mb). Hint: use `?lidR::readLAS` and see what `select` options might help.
#' 
#' ## Conclusion
#' 
#' This concludes our tutorial on the basic usage of the `lidR` package in R for processing and analyzing LiDAR data. We covered loading LiDAR data, inspecting and visualizing the data, selecting specific attributes, filtering points, and validating LAS objects for issues.
