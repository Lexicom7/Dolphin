---
title: 'Caribbean and Gulf of Mexico distribution model for key dolphin species.'
author: "Luis David Almeida Famada. EAGLE - Applied Earth Observation and Geoanalysis, University of Würzburg"
date: "2023-04-15"
output: 
  html_document:
    toc: true
    toc_depth: 3
    toc_float: true
    collapsed: true
    smooth_scroll: true
    theme: journal
    highlight: kate
    df_print: paged
    code_folding: show
      
---

# Introduction

This is a final project of the course: **Introduction to Programming and Geostatistics within the Master Program Applied Earth Observation and Geoanalysis** at the University of Würzburg with Professors Dr. Martin Wegmann and Dr. Jackob Schwalb-Willmann.

The project consists of several stages and has as a long-term goal the application of ecological niche modeling tools to determine the areas of highest species richness of the _Delphinidae_ family and their interaction with the main fishing grounds and marine protected areas of the Western North Atlantic, Gulf of Mexico and the Caribbean Sea.

In this first phase, we will start with an academic exercise using the _Stenella clymene_ or Clymene dolphin species to define the basis of the algorithm to be used for all species. Clymene dolphins are also known as "short-snouted spinners" because they often spin when they jump out of the water. Although there are no known major conservation problems for this species, it is likely that there are some undocumented problems. Some dolphins are killed in directed fisheries in the Caribbean, and others incidentally in nets throughout most of their range _(Jefferson, 2018)_

![Clymene dolphin. ](SCHello.jpg){withd=80%}

The Clymene dolphin is found only in the Atlantic Ocean, in tropical to warm-temperate waters. The exact range is not well documented, especially in South Atlantic, Mid-Atlantic and West African waters.  Most sightings have been in deep, offshore waters, although Clymene dolphins are sometimes observed very close to shore where deep water approaches the coast (such as around some Caribbean islands). They are present year-round at least in the northern Gulf of Mexico and probably throughout much of their tropical range.

![Clymene dolphin distribution. Source: Jefferson (2015) ](SCDist.jpg){withd=90%}
# 1 Preparation of oceanographic variables and presence data.

After this stage, the code for the processing of the 13 most representative species of the family _Delphinidae_ in these regions will be presented. In these first stages we will work only with the **GBIF (Global Biodiversity Information Facility)** databases, then we will improve the presence databases by relying on information from Cuban institutions and other publications in the region _(Barragán-Barrera et al., 2019)_. 

The oceanographic variables used were selected from the MARSPEC _(Sbrocco et al., 2013)_ and Bio-Oracle _(Assis et al., 2018)_ datasets. These two data sets are the most widely used globally for species distribution modeling _(Melo-Merino et al., 2020)_. Please refer to the respective websites for a description of the variables and the methodology used:

**MARSPEC**: http://www.marspec.org/
**Bio-ORACLE**: https://www.bio-oracle.org/index.php

It is possible to access these layers in many ways using APIs from within the R environment. One of the most widely used is the **"sdmpredictors"** package from Samuel Bosch _(https://github.com/lifewatch/sdmpredictors)_. For the purposes of this work and to save time in the reproduction of the code, all the previous information was compiled by preprocessing each layer and preparing it for later use by the models. In this link _(https://mega.nz/folder/QNJAwBqA#rYw4nxo4wa8KJ3CEgbVy4Q)_ you can download the layers already trimmed and resampled for the study area. The resampling codes can be found in the project folder.

It was necessary to investigate the availability of presence data in the GBIF database. For this purpose, the "dismo" package of Hijmans and collaborators _(https://rspatial.org/raster/sdm)_ was mainly used.

Exploration of the availability of points of presence and species selection:
```{r message=FALSE, warning=FALSE}
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
shp_ROI <- shapefile("D:/0_My_Analysis/Dolphin_distribution/Dolphin/Vector_Layers/Study_Area.shp")
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

# Combine species names and record counts
species_records <- data.frame(Species = sapply(species_list, function(x) paste(x$genus, x$species)), RecordCount = record_counts)

# Order the species by the number of records
species_records <- species_records[order(species_records$RecordCount, decreasing = TRUE), ]

# Create a bar plot
ggplot(species_records, aes(x = reorder(Species, -RecordCount), y = RecordCount, fill = Species)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = sample(colors(), length(species_list))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Number of records available for each of the most representative dolphin species",
       x = "Species",
       y = "Number of records downloaded") +
  guides(fill = FALSE)

```