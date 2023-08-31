#' ---
#' title: "Regions of Interest"
#' ---
#' 
#' ## Relevant Resources
#' 
#' -   [Code](https://github.com/tgoodbody/lidRtutorial/blob/main/R/02_roi.R)
#' 
#' -   [lidRbook section](https://r-lidar.github.io/lidRbook/engine.html#engine-clip)
#' 
#' ## Overview
#' 
#' We demonstrate the selection of regions of interest (ROIs) from LiDAR data. Geometries like circles and rectangles are selected based on coordinates. Complex geometries are extracted from shapefiles to clip specific areas.
#' 
#' ## Environment
#' 

# Clear environment
rm(list = ls(globalenv()))

# Load packages
library(lidR)
library(sf)

#' 
#' ## Simple Geometries
#' 
#' ### Load LiDAR Data and Inspect
#' 
#' We start by loading some LiDAR data and inspecting its header and number of point records.
#' 

las <- readLAS(files = "data/MixedEucaNat_normalized.laz", filter = "-set_withheld_flag 0")

# Inspect the header and the number of point records
las@header
las@header$`Number of point records`

#' 
#' ### Select Circular and Rectangular Areas
#' 
#' We can select circular and rectangular areas from the LiDAR data based on specified coordinates and radii or dimensions.
#' 

# Establish coordinates
x <- 203890
y <- 7358935

# Select a circular area
circle <- clip_circle(las = las, xcenter = x, ycenter = y, radius = 30)

# Inspect the circular area and the number of point records
circle
circle@header$`Number of point records`

# Visualize the LiDAR data with a default color palette
plot(circle, bg = "white")

#' 
#' We can do the same with a rectangular area by defining corner coordinates.
#' 

# Select a rectangular area
rect <- clip_rectangle(las = las, xleft = x, ybottom = y, xright = x + 40, ytop = y + 30)

# Visualize the LiDAR data with a default color palette
plot(rect, bg = "white")

#' 
#' We can also supply multiple coordinate pairs to clip multiple ROIs.
#' 

# Select multiple random circular areas
x <- runif(2, x, x)
y <- runif(2, 7358900, 7359050)

plots <- clip_circle(las = las, xcenter = x, ycenter = y, radius = 10)

# Visualize the LiDAR data with a default color palette
plot(plots[[1]], bg = "white")

# Visualize the LiDAR data with a default color palette
plot(plots[[2]], bg = "white")

#' 
#' ## Extraction of Complex Geometries from Shapefiles
#' 
#' We demonstrate how to extract complex geometries from shapefiles using the `clip_roi()` function from the `lidR` package.
#' 
#' ::: callout-caution
#' ## Reminder about **legacy packages**
#' 
#' `maptools`, `rgdal`, and `rgeos`, underpinning the `sp` package, will retire in October 2023. Please refer to R-spatial evolution reports for details, especially https://r-spatial.org/r/2023/05/15/evolution4.html.
#' :::
#' 
#' We use the `sf` package to load an ROI and then clip to its extents.
#' 

# Load the shapefile using sf
planting <- sf::st_read(dsn = "data/shapefiles/MixedEucaNat.shp", quiet = TRUE)

# Plot the LiDAR header information without the map
plot(las@header, map = FALSE)

# Plot the planting areas on top of the LiDAR header plot
plot(planting, add = TRUE, col = "#08B5FF39")

# Extract points within the planting areas using clip_roi()
eucalyptus <- clip_roi(las = las, geometry = planting)

# Plot the extracted points within the planting areas
plot(eucalyptus, bg = "white")

#' 
#' ## Exercises and Questions
#' 
#' Now, let's read a shapefile called `MixedEucaNatPlot.shp` using `sf::st_read()` and plot it on top of the LiDAR header plot.
#' 

# Read the shapefile "MixedEucaNatPlot.shp" using st_read()
plots <- sf::st_read(dsn = "data/shapefiles/MixedEucaNatPlot.shp", quiet = TRUE)

# Plot the LiDAR header information without the map
plot(las@header, map = FALSE)

# Plot the extracted points within the planting areas
plot(plots, add = TRUE)

#' 
#' #### E1.
#' 
#' Clip the 5 plots with a radius of 11.3 m.
#' 
#' #### E2.
#' 
#' Clip a transect from A `c(203850, 7358950)` to B `c(203950, 7959000)`.
#' 
#' #### E3.
#' 
#' Clip a transect from A `c(203850, 7358950)` to B `c(203950, 7959000)` but reorient it so it is no longer on the XY diagonal. Hint = `?clip_transect`
#' 
#' ## Conclusion
#' 
#' This concludes our tutorial on selecting simple geometries and extracting complex geometries from shapefiles using the `lidR` package in R.
