library(terra)
library(voluModel)
library(ggplot2)
library(lattice)
library(gridExtra)
library(ggplotify)
library(gridGraphics)
library(gtable)
library(gifski)
library(rnaturalearth)
library(ggpubr)
library(cowplot)

setwd("~/Dropbox/MARDIGRA/data/")

# Define functions ----
# Reads in raw model outputs and puts them all in the same extent
squareUp <- function(fileName, template){
  temp <- rast(fileName)
  temp <- classify(extend(temp, 
                          y = template), 
                   matrix(c(NA,0,1,1), nrow = 2, byrow = TRUE), 
                   include.lowest = TRUE)
  return(temp)
}
  
# Finds a named layer within a raster stack
layerChooser <- function(rasterToChoose, value){
  chosenLayer <- tryCatch(expr = rasterToChoose[[which(names(rasterToChoose) %in% value)]], 
                          error = function(e) NULL)
  return(chosenLayer)
}

# Calculates 3D biodiversity volumes from input list of raster stacks with 3D species distributions
biodiversity3D <- function(inputList, outputTemplate){
  result <- outputTemplate
  values(result) <- 0
  
  for (i in names(result)){
    template <- layerChooser(outputTemplate, i)
    layerList <- lapply(inputList, FUN = function(x) layerChooser(x, i))
    layerList <- layerList[unlist(lapply(layerList, FUN = function(x) !is.null(x)))]
    divLayer <- diversityStack(layerList, template = template)
    names(divLayer) <- i
    divLayer <- mask(divLayer, mask = template)
    result[[which(names(result) %in% i)]] <- divLayer
  }
  return(result)
}

# Flattens 3D distributions into horizontal extents
flattened3Ddist <- function(dist3D){
  layerSum <- sum(dist3D)
  flatDistribution <- as.numeric(layerSum > 0)
  return(flatDistribution)
}

# Makes a vertical transect through a 3D spatial volume
verticalSample <- function(x, sampleAxis = "lon", axisValue = NA){
  samplePoints <- as.data.frame(x, xy = TRUE, na.rm=FALSE)
  
  # Takes mean axis value if none specified
  if(is.na(axisValue)){
    if(sampleAxis == "lon"){
      axisValue <-  mean(c(xmin(x), xmax(x)))
      axisValue <- samplePoints[which.min(abs(samplePoints[,"x"] - axisValue)), "x"]
    } else{
      axisValue <-  mean(c(ymin(x), ymax(x)))
      axisValue <- samplePoints[which.min(abs(samplePoints[,"y"] - axisValue)), "y"]
    }
  } else {
    if(sampleAxis == "lon"){
      axisValue <- samplePoints[which.min(abs(samplePoints[,"x"] - axisValue)), "x"]
    } else{
      axisValue <- samplePoints[which.min(abs(samplePoints[,"y"] - axisValue)), "y"]
    }
  }
  
  # Select coordinates to sample
  if(sampleAxis == "lon"){
    samplePoints <- samplePoints[samplePoints[,"x"] == axisValue,]
  } else{
    samplePoints <- samplePoints[samplePoints[,"y"] == axisValue,]
  }
  
  samplePoints <- lapply(1:nrow(samplePoints), FUN = function(x){
    cbind(rep(samplePoints[x,"x"], n = (ncol(samplePoints)-2)), 
          rep(samplePoints[x,"y"], n = (ncol(samplePoints)-2)),
          as.numeric(colnames(samplePoints[c(-1,-2)])),
          t(samplePoints[x,])[c(-1,-2),])
  })
  
  samplePoints <- as.data.frame(do.call(rbind, samplePoints))
  colnames(samplePoints) <- c("x", "y", "z", "value")
  
  if(all(is.na(samplePoints$z))){
    samplePoints <- samplePoints[, !names(samplePoints) == "z"]
  } else{  
    zVals <- sort(unique(samplePoints$z))
    zHeights <- vector(mode = "numeric", 
                       length = length(zVals))
    for(i in 1:length(zHeights)){
      if ( i == length(zHeights)){
        zHeights[[i]] <- abs(zVals[[i]] - zVals[[i-1]])
      } else{
        zHeights[[i]] <- abs(zVals[[i]] - zVals[[i+1]])
      }
    }
    
    height <- vector(mode = "numeric", length = nrow(samplePoints))
    for(i in 1:nrow(samplePoints)){
      height[[i]] <- zHeights[which(zVals %in% samplePoints$z[[i]])]
    }
    samplePoints$height <- height
}

  samplePoints <- samplePoints[complete.cases(samplePoints),]
  
  return(samplePoints)
}

# Plot transect as calculated from verticalSample()
# Adapted from https://semba-blog.netlify.app/05/10/2020/heatmaps-in-r-with-ggplot2-and-metr-packages/
transectPlot <- function(verticalTransect, 
                         scaleRange = NA,
                         verbose = FALSE,
                         plotLegend = TRUE,
                         depthLim = 7000,
                         transRange = c(-90,90),
                         transTicks = 20,
                         #land = NA, landCol = "black",
                         ...){
  #Input processing
  args <- list(...)
  
  if("option" %in% names(args)){
    option <- args$option
  } else{
    option <- "plasma"
  }
  
  if("n" %in% names(args)){
    n <- args$n
  } else{
    n <- 11
  }
  
  if("legendRound" %in% names(args)){
    legendRound <- args$legendRound
  } else{
    legendRound <- 2
  }
  
  if(!all(is.na(scaleRange))){
    if(!all(any("numeric" %in% class(scaleRange)),
            length(scaleRange) == 2)){
      warning(paste0("'scaleRange' must either be NA or\na numeric vector of length 2."))
      return(NULL)
    }
  }
  
  if(any(is.na(scaleRange))){
    begin <- min(verticalTransect$value)
    end <- max(verticalTransect$value)
  } else if(any(min(verticalTransect$value) < min(scaleRange),
                max(verticalTransect$value) > max(scaleRange))){
    begin <- min(verticalTransect$value)
    end <- max(verticalTransect$value)
    
    if(verbose){
      message(paste0("Input extremes exceed specified scale range.\n",
                     "Using input max and min instead."))
    }
  } else{
    begin <- min(scaleRange)
    end <- max(scaleRange)
  }
  
  at <- seq(from = begin, to = end, by = (end-begin)/n)
  at <- round(at, legendRound)[c(-1,-length(at))]
  
  begin <- (begin - min(verticalTransect$value)/diff(scaleRange))
  end <- max(verticalTransect$value)/end
  
  # Tiled heatmap
  transect.tile <- ggplot(verticalTransect) +
    geom_tile(aes(x = y, y = z, fill = value, height = height), na.rm = TRUE) +
    coord_cartesian(expand = FALSE) +
    labs(x = NULL, y = "Water Depth (m)") +
    scale_fill_viridis_b(option = option, begin = begin, end = end, breaks = at) +
    theme_bw() %+%
    theme(panel.background = element_rect(fill = "grey90"),
          panel.grid.major = element_line(linetype = 3, colour = "grey60"),
          axis.text = element_text(colour = 1, size = 10),
          axis.title = element_text(colour = 1, size = 12)) +
    scale_y_reverse() +
    metR::scale_x_latitude(ticks = transTicks, position = "bottom", 
                           breaks = seq(transRange[[1]], transRange[[2]], by = transTicks)) +
    guides(fill = FALSE) +
    theme(plot.margin = margin(.5,.5,.5,0, "cm")) +
    ylim(depthLim, 0) +
    xlim(transRange)
  
  if(plotLegend){
    transect.tile <- transect.tile +
      guides(fill = guide_legend(title.position = "right",direction = "vertical",
                                 title.theme = element_text(angle = 90, size = 12, colour = "black"),
                                 barheight = .5, barwidth = .95,
                                 title.hjust = 0.5, raster = FALSE,
                                 title = "Value"))
  }
  return(transect.tile) 
}

grab_grob <- function(){
  grid.echo()
  grid.grab(wrap.grobs = TRUE)
}

sideBySide <- function(map2d, map3d, transectLocation, title, 
                       land, location, count, depthLim = 7000){
  transect3D <- verticalSample(map3d, sampleAxis = "lon",
                               axisValue = transectLocation)
  
  rng <- c(minmax(map2d), transect3D$value)
  rng <- c(min(rng), max(rng))
  
  horizontalPlot <- oneRasterPlot(map2d, 
                                  title = title,
                                  land = land,
                                  scaleRange = rng, 
                                  option = "plasma",
                                  plotLegend = F,
                                  n = 11,
                                  legendRound = 0)
  abline(v = transectLocation, col = "red", lwd = 2)
  horizontalPlot <- recordPlot()
  plot.new()
  horizontalPlot <- as_ggplot(as_grob(horizontalPlot))
  
  plottedTransect <- transectPlot(transect3D, 
                                  n = 11, legendRound = 0,
                                  scaleRange = rng,
                                  option = "plasma",
                                  plotLegend = F, 
                                  depthLim = depthLim,
                                  transRange = ext(map2d)[3:4])
  
#  if(transectLocation < 0 ){layerNum <- paste0("-",layerNum)}
  png(paste0(location, title, sprintf("%03d",count), ".png"), 
      width = 800, height = 400)
  grid.arrange(horizontalPlot,
               as.grob(plottedTransect),
               nrow = 1,heights = 5)
  dev.off()
}

# Load data ----
# Atlantic Ocean shapefile
land <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf")[1]

atlanticShapefile <- vect("FAO Fishing Areas 2005/FAO_AREAS.shp")
atlanticShapefile <- aggregate(atlanticShapefile[atlanticShapefile$OCEAN=="Atlantic",], dissolve = T)
atlanticShapefile$OCEAN <- "Atlantic"

temperature <- rast("EnvironmentalData/ProcessedEnvtData/temperature.tif")

studyRegion <- crop(temperature, atlanticShapefile)
studyRegion <- crop(studyRegion, c(ext(atlanticShapefile)[1],30,-60,70))

# Loading and organizing data ----
taxonomy <- read.csv("taxaWithinAreaOfInterest_curated8Jan.csv")
taxonomy <- taxonomy[taxonomy$isAtlantic,]
gadList <- taxonomy$EschmeyerName[taxonomy$Group == "G"]
scombList <- taxonomy$EschmeyerName[taxonomy$Group == "S"]
belonList <- taxonomy$EschmeyerName[taxonomy$Group == "B"]

# Get list of species that have both 2D and 3D GLMs
glm3D <- list.files(path = "Models/GLM_3D/", pattern = ".tif", full.names = F)
glm3Dkludge <- gsub(glm3D, pattern = ".tif", replacement = "\\.asc")
glm2D <- list.files(path = "Models/GLM_2D/", pattern = "\\.asc", full.names = F)

completeList <- table(sort(c(glm3Dkludge, glm2D)))
completeList <- names(completeList[completeList == 2])
spList <- gsub(completeList, pattern = "\\.asc", replacement = "")

rm(atlanticShapefile, taxonomy, completeList, glm2D, glm3D, glm3Dkludge)

# Beloniformes ----
beloniformes <- spList[spList %in% belonList]

# Envelope
env3Dbelon <- lapply(paste0("Models/Envelope_3D/", beloniformes, ".tif"), 
                     FUN = function(x) squareUp(fileName = x, studyRegion[[1]]))

flattenedEnv3Dbelon <- lapply(env3Dbelon, FUN = function(x) flattened3Ddist(x))

env3DbelonFlat <- diversityStack(flattenedEnv3Dbelon, template = studyRegion[[1]])
terra::writeRaster(env3DbelonFlat, filename = "DiversityEstimates/BeloniformesEnvelope3D_horizontal.tif",
                   overwrite = T)
env3DbelonFlat <- rast("DiversityEstimates/BeloniformesEnvelope3D_horizontal.tif")

env3DbelonVolume <- biodiversity3D(env3Dbelon, studyRegion)
terra::writeRaster(env3DbelonVolume, filename = "DiversityEstimates/BeloniformesEnvelope3D.tif",
                   overwrite = T)
env3DbelonVolume <- rast("DiversityEstimates/BeloniformesEnvelope3D.tif")

env2Dbelon <- lapply(paste0("Models/Envelope_2D/", beloniformes, ".asc"), 
                     FUN = function(x) squareUp(fileName = x, studyRegion[[1]]))
env2Dbelon <- diversityStack(env2Dbelon, studyRegion[[1]])
terra::writeRaster(env2Dbelon, 
                   filename = "DiversityEstimates/BeloniformesEnvelope2D.tif",
                   overwrite = T)

# Test vis, Env, Beloniformes ----
env3DbelonTransect <- verticalSample(env3DbelonVolume, sampleAxis = "lon", axisValue = -30)

rng <- c(minmax(env3DbelonFlat), env3DbelonTransect$value)
rng <- c(min(rng), max(rng))

land <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf")[1]
horizontalPlot <- oneRasterPlot(env3DbelonFlat, title = "Atlantic Diversity of Beloniformes,\nEnvelope",
                                land = land, scaleRange = rng)
abline(v = -30, col = "red", lwd = 2)
transectPlot(env3DbelonTransect, n = 10, scaleRange = rng)

rm(env3Dbelon, env3DbelonFlat, env3DbelonTransect, env3DbelonVolume, 
   flattenedEnv3Dbelon, test, horizontalPlot)

# GLM, Beloniformes ----
glm3Dbelon <- lapply(paste0("Models/GLM_3D/", beloniformes, ".tif"),
                     FUN = function(x) squareUp(fileName = x, studyRegion[[1]]))

flattenedglm3Dbelon <- lapply(glm3Dbelon, FUN = function(x) flattened3Ddist(x))

glm3DbelonFlat <- diversityStack(flattenedglm3Dbelon, template = studyRegion[[1]])
terra::writeRaster(glm3DbelonFlat, filename = "DiversityEstimates/BeloniformesGLM3D_horizontal.tif",
                   overwrite = T)
glm3DbelonFlat <- rast("DiversityEstimates/BeloniformesGLM3D_horizontal.tif")

glm3DbelonVolume <- biodiversity3D(glm3Dbelon, studyRegion)
terra::writeRaster(glm3DbelonVolume, filename = "DiversityEstimates/BeloniformesGLM3D.tif",
                   overwrite = T)
glm3DbelonVolume <- rast("DiversityEstimates/BeloniformesGLM3D.tif")

glm2Dbelon <- lapply(paste0("Models/GLM_2D/", beloniformes, ".asc"), 
                     FUN = function(x) squareUp(fileName = x, studyRegion[[1]]))
glm2Dbelon <- diversityStack(glm2Dbelon, studyRegion[[1]])
terra::writeRaster(glm2Dbelon, filename = "DiversityEstimates/BeloniformesGLM2D.tif",
                   overwrite = T)

# Test vis, GLM, Beloniformes ----
glm3DbelonTransect <- verticalSample(glm3DbelonVolume, sampleAxis = "lon", axisValue = -30)

rng <- c(minmax(glm3DbelonFlat), glm3DbelonTransect$value)
rng <- c(min(rng), max(rng))

horizontalPlot <- oneRasterPlot(glm3DbelonFlat, title = "Atlantic Diversity of Beloniformes, GLM",
                                land = land, scaleRange = rng)
abline(v = -30, col = "red", lwd = 2)
transectPlot(glm3DbelonTransect, n = 10, scaleRange = rng)

rm(glm3Dbelon, glm3DbelonFlat, glm3DbelonTransect, glm3DbelonVolume, 
   flattenedglm3Dbelon, horizontalPlot)

# Scombriformes ----
scombriformes <- spList[spList %in% scombList]

# Envelope
env3Dscomb <- lapply(paste0("Models/Envelope_3D/", scombriformes, ".tif"), 
                     FUN = function(x) squareUp(fileName = x, studyRegion[[1]]))

flattenedEnv3Dscomb <- lapply(env3Dscomb, FUN = function(x) flattened3Ddist(x))

env3DscombFlat <- diversityStack(flattenedEnv3Dscomb, template = studyRegion[[1]])
terra::writeRaster(env3DscombFlat, filename = "DiversityEstimates/scombriformesEnvelope3D_horizontal.tif",
                   overwrite = T)
env3DscombFlat <- rast("DiversityEstimates/scombriformesEnvelope3D_horizontal.tif")

env3DscombVolume <- biodiversity3D(env3Dscomb, studyRegion)
terra::writeRaster(env3DscombVolume, filename = "DiversityEstimates/scombriformesEnvelope3D.tif",
                   overwrite = T)
env3DscombVolume <- rast("DiversityEstimates/scombriformesEnvelope3D.tif")

env2Dscomb <- lapply(paste0("Models/Envelope_2D/", scombriformes, ".asc"), 
                     FUN = function(x) squareUp(fileName = x, studyRegion[[1]]))
env2Dscomb <- diversityStack(env2Dscomb, studyRegion[[1]])
terra::writeRaster(env2Dscomb, 
                   filename = "DiversityEstimates/ScombriformesEnvelope2D.tif",
                   overwrite = T)

# Test vis, Env, Scombriformes ----
env3DscombTransect <- verticalSample(env3DscombVolume, sampleAxis = "lon", axisValue = -30)

rng <- c(minmax(env3DscombFlat), env3DscombTransect$value)
rng <- c(min(rng), max(rng))

horizontalPlot <- oneRasterPlot(env3DscombFlat, title = "Atlantic Diversity of Scombriformes,\nEnvelope",
                                land = land, scaleRange = rng)
abline(v = -30, col = "red", lwd = 2)
transectPlot(env3DscombTransect, n = 10, scaleRange = rng)

# GLM, Scombriformes ----
glm3Dscomb <- lapply(paste0("Models/GLM_3D/", scombriformes, ".tif"), 
                     FUN = function(x) squareUp(fileName = x, studyRegion[[1]]))

flattenedglm3Dscomb <- lapply(glm3Dscomb, FUN = function(x) flattened3Ddist(x))

glm3DscombFlat <- diversityStack(flattenedglm3Dscomb, template = studyRegion[[1]])
terra::writeRaster(glm3DscombFlat, filename = "DiversityEstimates/scombriformesGLM3D_horizontal.tif",
                   overwrite = T)
glm3DscombFlat <- rast("DiversityEstimates/scombriformesGLM3D_horizontal.tif")

glm3DscombVolume <- biodiversity3D(glm3Dscomb, studyRegion)
terra::writeRaster(glm3DscombVolume, filename = "DiversityEstimates/scombriformesGLM3D.tif",
                   overwrite = T)
glm3DscombVolume <- rast("DiversityEstimates/scombriformesGLM3D.tif")

glm2Dscomb <- lapply(paste0("Models/GLM_2D/", scombriformes, ".asc"), 
                     FUN = function(x) squareUp(fileName = x, studyRegion[[1]]))
glm2Dscomb <- diversityStack(glm2Dscomb, studyRegion[[1]])
terra::writeRaster(glm2Dscomb, 
                   filename = "DiversityEstimates/ScombriformesGLM2D.tif",
                   overwrite = T)

#How is it?
glm3DscombTransect <- verticalSample(glm3DscombVolume, sampleAxis = "lon", axisValue = -30)

oneRasterPlot(glm3DscombFlat, title = "Atlantic Diversity of Scombriformes, GLM",
                                land = land, scaleRange = rng)
abline(v = -30, col = "red", lwd = 2)
transectPlot(glm3DscombTransect, n = 10, scaleRange = rng)

rm(env3Dscomb, env3DscombFlat, env3DscombTransect, env3DscombVolume, flattenedEnv3Dscomb)

# Gadiformes ----
gadiformes <- spList[spList %in% gadList]

# Envelope
env3Dgad <- lapply(paste0("Models/Envelope_3D/", gadiformes, ".tif"), 
                   FUN = function(x) squareUp(fileName = x, studyRegion[[1]]))

flattenedEnv3Dgad <- lapply(env3Dgad, FUN = function(x) flattened3Ddist(x))

env3DgadFlat <- diversityStack(flattenedEnv3Dgad, template = studyRegion[[1]])
terra::writeRaster(env3DgadFlat, filename = "DiversityEstimates/gadiformesEnvelope3D_horizontal.tif",
                   overwrite = T)
env3DgadFlat <- rast("DiversityEstimates/gadiformesEnvelope3D_horizontal.tif")

env3DgadVolume <- biodiversity3D(env3Dgad, studyRegion)
terra::writeRaster(env3DgadVolume, filename = "DiversityEstimates/gadiformesEnvelope3D.tif",
                   overwrite = T)
env3DgadVolume <- rast("DiversityEstimates/gadiformesEnvelope3D.tif")

env2Dgad <- lapply(paste0("Models/Envelope_2D/", gadiformes, ".asc"), 
                     FUN = function(x) squareUp(fileName = x, studyRegion[[1]]))
env2Dgad <- diversityStack(env2Dgad, studyRegion[[1]])
terra::writeRaster(env2Dgad, 
                   filename = "DiversityEstimates/GadiformesEnvelope2D.tif",
                   overwrite = T)

#How is it?
env3DgadTransect <- verticalSample(env3DgadVolume, sampleAxis = "lon", axisValue = -30)

rng <- c(minmax(env3DgadFlat), env3DgadTransect$value)
rng <- c(min(rng), max(rng))

oneRasterPlot(env3DgadFlat, title = "Atlantic Diversity of gadiformes,\nEnvelope",
                                land = land, scaleRange = rng)
abline(v = -30, col = "red", lwd = 2)
transectPlot(env3DgadTransect, n = 10, scaleRange = rng)

# GLM
glm3Dgad <- lapply(paste0("Models/GLM_3D/", gadiformes, ".tif"), 
                   FUN = function(x) squareUp(fileName = x, studyRegion[[1]]))

flattenedglm3Dgad <- lapply(glm3Dgad, FUN = function(x) flattened3Ddist(x))

glm3DgadFlat <- diversityStack(flattenedglm3Dgad, template = studyRegion[[1]])
terra::writeRaster(glm3DgadFlat, filename = "DiversityEstimates/gadiformesGLM3D_horizontal.tif",
                   overwrite = T)
glm3DgadFlat <- rast("DiversityEstimates/gadiformesGLM3D_horizontal.tif")

glm3DgadVolume <- biodiversity3D(glm3Dgad, studyRegion)
terra::writeRaster(glm3DgadVolume, filename = "DiversityEstimates/gadiformesGLM3D.tif",
                   overwrite = T)
glm3DgadVolume <- rast("DiversityEstimates/gadiformesGLM3D.tif")

glm2Dgad <- lapply(paste0("Models/GLM_2D/", gadiformes, ".asc"), 
                     FUN = function(x) squareUp(fileName = x, studyRegion[[1]]))
glm2Dgad <- diversityStack(glm2Dgad, studyRegion[[1]])
terra::writeRaster(glm2Dgad, filename = "DiversityEstimates/GadiformesGLM2D.tif",
                   overwrite = T)

#How is it?
glm3DgadTransect <- verticalSample(glm3DgadVolume, sampleAxis = "lon", axisValue = -30)

land <- vect(ne_countries(scale = "medium", returnclass = "sf")[1])
land <- aggregate(land)
land <- simplifyGeom(land)
land <- crop(land, glm3DgadFlat)

rng <- c(minmax(glm3DgadFlat), glm3DgadTransect$value)
rng <- c(min(rng), max(rng))

oneRasterPlot(glm3DgadFlat, title = "Atlantic Diversity of gadiformes, GLM",
                                land = land, scaleRange = rng)
abline(v = -30, col = "red", lwd = 2)
transectPlot(glm3DgadTransect, n = 10, scaleRange = rng)

# Diversity data loading ----
bel3Dglm <- crop(rast("DiversityEstimates/BeloniformesGLM3D.tif"),
                 studyRegion)
gad3Dglm <- crop(rast("DiversityEstimates/gadiformesGLM3D.tif"),
                 studyRegion)
scomb3Dglm <- crop(rast("DiversityEstimates/scombriformesGLM3D.tif"),
                   studyRegion)

glm3DbelonFlat <- crop(rast("DiversityEstimates/BeloniformesGLM3D_horizontal.tif"), studyRegion)
glm3DgadFlat <- crop(rast("DiversityEstimates/GadiformesGLM3D_horizontal.tif"),studyRegion)
glm3DscombFlat <- crop(rast("DiversityEstimates/ScombriformesGLM3D_horizontal.tif"), studyRegion)

# Plotting ----
# Beloniformes
count <- 1
for (i in -97:20){
  sideBySide(map2d = glm3DbelonFlat, map3d = bel3Dglm,
             transectLocation = i, 
             title = "Beloniformes",
             land = land, 
             location = "../../JMIH 2023/Beloniformes/", 
             count = count,
             depthLim = 500)
  count <- count + 1
}

png_files <- list.files("../../JMIH 2023/Beloniformes/", pattern = ".*png$", full.names = TRUE)
gifski(png_files, gif_file = "../../JMIH 2023/Beloniformes/animation.gif", width = 800, height = 400, delay = .2)

# Scombriformes
count <- 1
for (i in -97:20){
  sideBySide(map2d = glm3DscombFlat, map3d = scomb3Dglm,
             transectLocation = i, title = "Scombriformes",
             land = land,
             location = "../../JMIH 2023/Scombriformes/", 
             count = count, depthLim = 6000)
  count <- count + 1
}

png_files <- list.files("../../JMIH 2023/Scombriformes/", pattern = ".*png$", full.names = TRUE)
gifski(png_files, gif_file = "../../JMIH 2023/Scombriformes/animation.gif", width = 800, height = 400, delay = .2)

# Gadiformes
count <- 1
for (i in -97:20){
  sideBySide(map2d = glm3DgadFlat, map3d = gad3Dglm,
             transectLocation = i, title = "Gadiformes",
             land = land,
             location = "../../JMIH 2023/Gadiformes/", 
             count = count, 6000)
  count <- count + 1
}

png_files <- list.files("../../JMIH 2023/Gadiformes/", pattern = ".*png$", full.names = TRUE)
gifski(png_files, gif_file = "../../JMIH 2023/Gadiformes/animation.gif", width = 800, height = 400, delay = .2)