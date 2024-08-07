# Supplemental Material
Scripts for the manuscript "Depth Matters for Marine Diversity", currently in prep.

## Occurrence Search
This script created and collated a checklist of Atlantic fishes in the orders of interest, searched GBIF and OBIS for all available occurrences. Finally, the script merges all records by species according to Catalog of Fishes taxonomy while performing some basic fit-for-use filtration.
* [Code:](https://github.com/hannahlowens/3DFishRichness/blob/main/1_OccurrenceSearch.R)

## Environmental Data Processing
This script processes raw .shp data from the [World Ocean Atlas 2018](https://www.ncei.noaa.gov/products/world-ocean-atlas) and bathymetry from the [ETOPO Global Relief Model](https://www.ncei.noaa.gov/products/etopo-global-relief-model) for use in niche modeling.
* [Code:](https://github.com/hannahlowens/3DFishRichness/blob/main/2_EnvironmentalDataProcessing.R)

## ENM Workflow
This workflow generates GLM and envelope models in 2D and 3D for each species using tools from the [voluModel package](https://cran.r-project.org/package=voluModel). Steps include further occurrence cleaning, training region generation, modeling, and generating model reports for use in the Appendix II atlas.
* [Code:](https://github.com/hannahlowens/3DFishRichness/blob/main/3_ENMWorkflow.R)

## Model Statistics Comparison
This script performs several kinds of comparisons among the models generated by the ENM workflow to compare the performance of 2D versus 3D models. 
* [Code:](https://github.com/hannahlowens/3DFishRichness/blob/main/4_ModelStatComparisons.R)

## Mapping Species Richness
This script generates several 2D and 3D estimates and visualizations of species richness, including 2D maps and animated gifs to show how diversity changes with longitude across the Atlantic.
* [Code:](https://github.com/hannahlowens/3DFishRichness/blob/main/5_SpeciesRichnessMapping.R)

## Analyzing Species Richness
This script models species richness (as estimated by 2D and 3D GLMs and envelope models) as a response to temperature, latitude, and nitrate concentration. 
* [Code:](https://github.com/hannahlowens/3DFishRichness/blob/main/6_SpeciesRichnessAnalysis.R)
