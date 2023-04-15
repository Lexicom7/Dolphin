# Load necessary libraries
library(raster)
library(sf)
library(sp)
library(dismo)

# Set the working directory
getwd()
setwd("D:/0_My_Analysis/Dolphin_distribution/Dolphin")

# Define folder path
BioOracle_path <- "D:/0_My_Analysis/Dolphin_distribution/Dolphin/BIO_Cliped"
MARSPECT_path <- "D:/0_My_Analysis/Dolphin_distribution/Dolphin/MARSPECT_Cliped"

# List all .tif and .tiff files in the folder
BioOracle_files <- list.files(path = BioOracle_path, pattern = "\\.(tif|tiff)$", full.names = TRUE, ignore.case = TRUE)
print(BioOracle_files)

MARSPECT_files <- list.files(path = MARSPECT_path, pattern = "\\.(tif|tiff)$", full.names = TRUE, ignore.case = TRUE)
print(MARSPECT_files)

# Create a function to resample rasters
resample_raster <- function(raster1, raster2) {
  resampled_raster <- resample(raster1, raster2, method = "bilinear")
  return(resampled_raster)
}

# Create a new folder called "Variables"
dir.create("MARSPECT_Resample")
output_folder <- "MARSPECT_Resample"

# Read the reference BioOracle raster
reference_marspect_raster <- raster(BioOracle_files[1])

# Loop through all rasters in the MARSPECT list
for (i in seq_along(MARSPECT_files)) {
  MAR_raster <- raster(MARSPECT_files[i])
  
  # Resample MARSPECT raster to match the reference BioOracle raster resolution
  resampled_MAR_raster <- resample_raster(MAR_raster, reference_marspect_raster)
  
  # Save the resampled rasters in the "Variables" folder with the same name
  output_MAR_path <- file.path(output_folder, basename(MARSPECT_files[i]))
  writeRaster(resampled_MAR_raster, output_MAR_path, format = "GTiff")
}

# Define the new folder paths for resample MARSPECT and create a list
MARSPECTOK_path <- "D:/0_My_Analysis/Dolphin_distribution/Dolphin/MARSPECT_Resample"
MARSPECTOK_files <- list.files(path = MARSPECTOK_path, pattern = "\\.(tif|tiff)$", full.names = TRUE, ignore.case = TRUE)
print(MARSPECTOK_files)


