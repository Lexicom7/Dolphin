### 1. Prepare the environment ###
# List of required packages
required_packages <- c("raster", "sf", "sp", "dismo", "maptools", "usdm", "rJava", "mapview")

# Install and load packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Load necessary libraries
library(raster)
library(sf)
library(sp)
library(dismo)
library(maptools)
library(usdm) #Optional
library(rJava)
library(mapview)
library(ggplot2)
library(dplyr)

# Set Java memory limit
options(java.parameters = "-Xmx12g") 

# Set the working directory
setwd("D:/0_My_Analysis/Dolphin_distribution/Dolphin")

### 2. Load environmental data and Region of Interest (ROI) ###
# Set paths to the environmental variables
MARSPECT_path <- "D:/0_My_Analysis/Dolphin_distribution/Dolphin/MARSPECT_Resample"
BioOracle_path <- "D:/0_My_Analysis/Dolphin_distribution/Dolphin/BIO_Cliped"

# List all .tif and .tiff files in the folders
BioOracle_files <- list.files(path = BioOracle_path, pattern = "\\.(tif|tiff)$", full.names = TRUE, ignore.case = TRUE)

MARSPECT_files <- list.files(path = MARSPECT_path, pattern = "\\.(tif|tiff)$", full.names = TRUE, ignore.case = TRUE)

# Read and stack rasters
reference_marspect_raster <- raster(MARSPECT_files[1])
bio_rasters <- stack(lapply(BioOracle_files, raster))
marspect_rasters <- stack(lapply(MARSPECT_files, raster))
Variables <- stack(bio_rasters, marspect_rasters)

# Print the image collection
names(Variables)

# Read and transform the ROI shapefile
shp_ROI <- shapefile("D:/0_My_Analysis/Dolphin_distribution/Dolphin/Vector_Layers/Study_Area.shp")
shp_ROI_transformed <- spTransform(shp_ROI, proj4string(reference_marspect_raster))



### 3. Prepare presence data ###
# Obtain species presence data from GBIF database
Stenella_clymene <- gbif(genus = "Stenella", #the genus name
                         species = "clymene", #the species name
                         ext = shp_ROI_transformed, #object to limit the geographic extent of the records
                         download=TRUE, #whether to download the records or just show the number
                         geo=TRUE, # whether to only download records that have a georeference
                         sp=FALSE) # whether to return a SpatialPointsDataFrame

# Explore the data
dim(Stenella_clymene)
colnames(Stenella_clymene)

# Visualize the presence points on a map
data(wrld_simpl)
plot(wrld_simpl, xlim=c(-100,-50), ylim=c(8,50),axes=TRUE, col="light yellow")
box() # restore the box around the map
points(Stenella_clymene$lon, Stenella_clymene$lat, col='red', cex=0.75) # add the points


# Data cleaning
# Seach empty records
lonzero = subset(Stenella_clymene, lon==0)
lonzero[, 1:13]

# Convert to SpatialPolygonsDataFrame and set CRS
coordinates(Stenella_clymene) <- ~lon+lat
crs(Stenella_clymene) <- crs(projection(reference_marspect_raster))


# Eliminate sample bias
e <- extent(Stenella_clymene)
r <- raster(e) # create a RasterLayer with the same extent
res(r) <- 0.25 # set the resolution of the cells to (for example) 1 degree
r <- extend(r, extent(r)+1)
Stenella_clymene_fix <- gridSample(Stenella_clymene, r, n=1) # sample


# Prepare data for Maxent
bg_points <- randomPoints(mask = Variables, #a Raster* object
                          n = 300, #Number of points
                          ext = extent(Stenella_clymene_fix),
                          p = Stenella_clymene_fix, #Presence points
                          excludep = TRUE) #If TRUE, presence points are excluded from background

# Convert to SpatialPointsDataFrame objects
bg_points <- SpatialPointsDataFrame(coords = bg_points, data = data.frame(id = 1:nrow(bg_points)), proj4string = CRS(proj4string(Variables)))
presence_points <- SpatialPointsDataFrame(coords = Stenella_clymene_fix, data = data.frame(id = 1:nrow(Stenella_clymene_fix)), proj4string = CRS(proj4string(Variables)))


# Split data into training and testing sets using k-fold
set.seed(0)

groupP <- kfold(presence_points, 5)
pres_train <- presence_points[groupP != 1, ]
pres_test <- presence_points[groupP == 1, ]

groupB <- kfold(bg_points, 5)
bg_train <- bg_points[groupB != 1, ]
bg_test <- bg_points[groupB == 1, ]

# Exclude variables with high correlation
pres_train_value <- extract(Variables, pres_train, 
                            method = "bilinear",  #can be set to 'simple' or 'bilinear'
                            small = FALSE, #If the function always returns a number, even when the buffer does not include the center of a single cell.
                            fun = mean, #specify a function to summarize the extracted values. max or mean
                            na.rm = TRUE) #whether to remove missing values from the data before applying the summary function (fun) to the extracted values. 

vifstep_var <- vifstep(pres_train_value)

Variables_Good <- exclude(Variables, vifstep_var)

### 4. Apply the Maxent model ###

# Function to extract variables in a larger area
my_extractor <- function(raster_layer, points) {
  extract(raster_layer, buffer(points, width=3)) %>%
    apply(2, function(x) c(mean(x), sd(x)))
  }

#Run Maxent model

maxent_model <- maxent(x=Variables_Good,
                       p=pres_train,
                       a=bg_train,
                       factors=NULL, 
                       removeDuplicates=TRUE,
                       extract=my_extractor,
                       args=c("-P", "autorun", "nowarnings"),
                       path="D:/0_My_Analysis/Dolphin_distribution/Dolphin/Maxent_model")

# Plot Maxent model results
plot(maxent_model) # Variable Contribution



# Evaluate the model with test points
evaluate_Max <- evaluate(pres_test, bg_test, maxent_model, Variables_Good)
evaluate_Max

# Make predictions using the model
Stenella_clymene_Prediction <- predict(maxent_model, Variables_Good,  ext=extent(reference_marspect_raster))



# Visualize Maxent raw values
mapview(Stenella_clymene_Prediction, 
        zcol = "value",
        map.types = "OpenStreetMap",
        legend = TRUE,
        alpha = 0.5,
        col.regions = viridis::viridis,
        main = "Maxent, raw values")

# Apply the threshold and create a presence/absence map
tr <- threshold(evaluate_Max, 'spec_sens')
presence_absence <- Stenella_clymene_Prediction > tr
mapview(presence_absence, 
        zcol = "value",
        map.types = "OpenStreetMap", 
        legend = TRUE,
        alpha = 0.3,
        col.regions = c("white", "red"),
        main = "presence/absence")

# Generate 'response plots'
.jinit() # initializes the Java Virtual Machine in R
response(maxent_model) # A response plot

### 5. Analysis of the result with Protected Areas and Exclusive Economic Zone ###

# Compare result with Exclusive Economic Zone and International Hydrographic Organization (EEZ - IHO) areas 
shp_EEZIHO <- shapefile("D:/0_My_Analysis/Dolphin_distribution/Dolphin/Vector_Layers/EEZ_IHOv42020.shp")
shp_EEZIHO_transformed <- spTransform(shp_EEZIHO, proj4string(reference_marspect_raster))

# Extract presence-absence values for each polygon
EEZIHO_ext <- extract(presence_absence, shp_EEZIHO_transformed)

# Calculate the count of presence pixels within each polygon
EEZIHO_counts <- sapply(EEZIHO_ext, function(x) sum(x == 1, na.rm = TRUE))

# Add the presence count as a new column to the shapefile
shp_EEZIHO_transformed@data$EEZIHO_counts <- EEZIHO_counts

# Calculate the sum of Presence_Count for each unique value in the IHO_SEA column
EEZIHO_sum <- aggregate(EEZIHO_counts ~ EEZ, data = shp_EEZIHO_transformed, FUN = sum)
EEZIHO_sum <- subset(EEZIHO_sum, EEZIHO_counts > 0)
EEZIHO_sum <- EEZIHO_sum[order(EEZIHO_sum$EEZIHO_counts, decreasing = TRUE), ]
print(EEZIHO_sum)


# Compare result with Protect Area (WDPA) areas
shp_WDPA <- shapefile("D:/0_My_Analysis/Dolphin_distribution/Dolphin/Vector_Layers/WDPA_Apr2023.shp")
shp_WDPA_transformed <- spTransform(shp_WDPA, proj4string(reference_marspect_raster))

# Extract presence-absence values for each polygon
WDPA_values <- extract(presence_absence, shp_WDPA_transformed)

# Calculate the count of presence pixels within each polygon
WDPA_counts <- sapply(WDPA_values, function(x) sum(x == 1, na.rm = TRUE))

# Add the presence count as a new column to the shapefile
shp_WDPA_transformed@data$WDPA_counts <- WDPA_counts


# Calculate the sum of Presence_Count for each unique value in the IHO_SEA column
WDPA_sum <- aggregate(WDPA_counts ~ ISO3, data = shp_WDPA_transformed, FUN = sum)
WDPA_sum <- subset(WDPA_sum, WDPA_counts > 0)
WDPA_sum <- WDPA_sum[order(WDPA_sum$WDPA_counts, decreasing = TRUE), ]
print(WDPA_sum)

