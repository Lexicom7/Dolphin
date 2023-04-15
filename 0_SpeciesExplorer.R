# Load required packages
library(dismo)
library(ggplot2)

# Path to the enviromental variables
MARSPECT_path <- "D:/0_My_Analysis/Dolphin_distribution/Dolphin/MARSPECT_Resample"
BioOracle_path <- "D:/0_My_Analysis/Dolphin_distribution/Dolphin/BIO_Cliped"

# List all .tif and .tiff files in the folder
BioOracle_files <- list.files(path = BioOracle_path, pattern = "\\.(tif|tiff)$", full.names = TRUE, ignore.case = TRUE)

MARSPECT_files <- list.files(path = MARSPECT_path, pattern = "\\.(tif|tiff)$", full.names = TRUE, ignore.case = TRUE)

# Read the reference MARSPECT raster
reference_marspect_raster <- raster(MARSPECT_files[1])


# Read shapefile ROI
shp_ROI <- shapefile("D:/0_My_Analysis/Dolphin_distribution/Dolphin/Study_Area.shp")
shp_ROI_transformed <- spTransform(shp_ROI, proj4string(reference_marspect_raster))



# Define the species list with genus and species names
species_list <- list(
  list(genus = "Lagenorhynchus", species = "albirostris"),
  list(genus = "Feresa", species = "attenuata"),
  list(genus = "Peponocephala", species = "electra"),
  list(genus = "Globicephala", species = "macrorhynchus"),
  list(genus = "Globicephala", species = "melas"),
  list(genus = "Pseudorca", species = "crassidens"),
  list(genus = "Grampus", species = "griseus"),
  list(genus = "Orcinus", species = "orca"),
  list(genus = "Steno", species = "bredanensis"),
  list(genus = "Stenella", species = "attenuata"),
  list(genus = "Tursiops", species = "truncatus"),
  list(genus = "Stenella", species = "coeruleoalba"),
  list(genus = "Phocoena", species = "phocoena")
)

# Initialize an empty vector to store the record counts
record_counts <- c()

# Loop through the species list and download records for each species
for (species_info in species_list) {
  species_data <- gbif(
    genus = species_info$genus,
    species = species_info$species,
    ext = shp_ROI_transformed,
    download = FALSE,
    geo = TRUE,
    removeZeros = TRUE,
    sp = FALSE
  )
  record_counts <- c(record_counts, print(species_data))
}
head(record_counts)

# Combine species names and record counts
species_records <- data.frame(Species = sapply(species_list, function(x) paste(x$genus, x$species)), RecordCount = record_counts)

# Order the species by the number of records
species_records <- species_records[order(species_records$RecordCount, decreasing = TRUE), ]

# Create a bar plot
ggplot(species_records, aes(x = reorder(Species, -RecordCount), y = RecordCount, fill = Species)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = sample(colors(), length(species_list))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Number of records downloaded for each species",
       x = "Species",
       y = "Number of records downloaded") +
  guides(fill = FALSE)