# Load necessary libraries
library(raster)
library(sf)
library(sp)
library(dismo)
library(mapview)

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
dir.create("Bio_Resample")
output_folder <- "Bio_Resample"

# Read the reference MARSPECT raster
reference_marspect_raster <- raster(MARSPECT_files[1])

# Loop through all rasters in the BioOracle list
for (i in seq_along(BioOracle_files)) {
  bio_raster <- raster(BioOracle_files[i])
  
  # Resample BioOracle raster to match the reference MARSPECT raster resolution
  resampled_bio_raster <- resample_raster(bio_raster, reference_marspect_raster)
  
  # Save the resampled rasters in the "Variables" folder with the same name
  output_bio_path <- file.path(output_folder, basename(BioOracle_files[i]))
  writeRaster(resampled_bio_raster, output_bio_path, format = "GTiff")
}

# Define the new folder paths for resample Bioracles and create a list
BioOracleOK_path <- "D:/0_My_Analysis/Dolphin_distribution/Dolphin/Bio_Resample"
BioOracleOK_files <- list.files(path = BioOracleOK_path, pattern = "\\.(tif|tiff)$", full.names = TRUE, ignore.case = TRUE)
print(BioOracleOK_files)


