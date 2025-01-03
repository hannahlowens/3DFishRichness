library(terra)
library(voluModel)
library(fields)
library(corrplot)
library(pscl)
library(car)
library(dplyr)
library(climateStability)
library(MASS)
library(sjPlot)
library(ggplot2)
library(ggpubr)
library(tidyr)

#setwd("YOUR_DIRECTORY")

# Modeling function
modelFunction <- function(rawData){
  result <- glm.nb(as.formula(paste0(colnames(rawData)[ncol(rawData)], "~ .")), 
                 data = rawData)
  modelReport <- NULL
  modelReport$sum <- summary(result)
  modelReport$coeff <- summary(result)$coefficients
  modelReport$performance <- performance::model_performance(result)
  modelReport$model <- result
  return(modelReport)
}

# Set up
theme_set(theme_classic())

# Data ----
land <- vect(rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")[1])
land <- aggregate(land)
land <- simplifyGeom(land)

atlanticShapefile <- vect("data/FAO Fishing Areas 2005/FAO_AREAS.shp")
atlanticShapefile <- aggregate(atlanticShapefile[atlanticShapefile$OCEAN=="Atlantic",], dissolve = T)
atlanticShapefile$OCEAN <- "Atlantic"

# Potential explanatory variables
nitrate <- crop(rast(x = "data/EnvironmentalData/ProcessedEnvtData/nitrate.tif"), 
                atlanticShapefile)[[1:100]]
temperature <- crop(rast(x = "data/EnvironmentalData/ProcessedEnvtData/temperature.tif"), 
                    atlanticShapefile)[[1:100]]

beloniformesEnv <- crop(rast(x = "data/DiversityEstimates/BeloniformesEnvelope3D.tif"), 
                        atlanticShapefile)[[1:100]]
beloniformesGLM <- crop(rast(x = "data/DiversityEstimates/BeloniformesGLM3D.tif"), 
                        atlanticShapefile)[[1:100]]
scombriformesEnv <- crop(rast(x = "data/DiversityEstimates/scombriformesEnvelope3D.tif"), 
                         atlanticShapefile)[[1:100]]
scombriformesGLM <- crop(rast(x = "data/DiversityEstimates/scombriformesGLM3D.tif"), 
                         atlanticShapefile)[[1:100]]
gadiformesEnv <- crop(rast(x = "data/DiversityEstimates/gadiformesEnvelope3D.tif"), 
                      atlanticShapefile)[[1:100]]
gadiformesGLM <- crop(rast(x = "data/DiversityEstimates/gadiformesGLM3D.tif"),
                      atlanticShapefile)[[1:100]]

beloniformesEnv2D <- crop(rast("data/DiversityEstimates/BeloniformesEnvelope2D.tif"),
                          atlanticShapefile)
beloniformesGLM2D <- crop(rast("data/DiversityEstimates/BeloniformesGLM2D.tif"),
                          atlanticShapefile)
scombriformesEnv2D <- crop(rast("data/DiversityEstimates/ScombriformesEnvelope2D.tif"),
                           atlanticShapefile)
scombriformesGLM2D <- crop(rast("data/DiversityEstimates/ScombriformesGLM2D.tif"),
                           atlanticShapefile)
gadiformesEnv2D <- crop(rast("data/DiversityEstimates/GadiformesEnvelope2D.tif"),
                        atlanticShapefile)
gadiformesGLM2D <- crop(rast("data/DiversityEstimates/GadiformesGLM2D.tif"),
                        atlanticShapefile)

bGLM3dFlat <- crop(rast(x = "data/DiversityEstimates/BeloniformesGLM3D_horizontal.tif"), 
                   atlanticShapefile)
gGLM3dFlat <- crop(rast(x = "data/DiversityEstimates/GadiformesGLM3D_horizontal.tif"), 
                   atlanticShapefile)
sGLM3dFlat <- crop(rast(x = "data/DiversityEstimates/ScombriformesGLM3D_horizontal.tif"), 
                   atlanticShapefile)

bEnv3DFlat <- crop(rast(x = "data/DiversityEstimates/BeloniformesEnvelope3D_horizontal.tif"), 
                   atlanticShapefile)
gEnv3dFlat <- crop(rast(x = "data/DiversityEstimates/GadiformesEnvelope3D_horizontal.tif"), 
                   atlanticShapefile)
sEnv3dFlat <- crop(rast(x = "data/DiversityEstimates/ScombriformesEnvelope3D_horizontal.tif"), 
                   atlanticShapefile)

# Latitudinal trends in data with depth ----
linePlotFunc <- function(x, name){
  plot(x[[1]], main = paste("Mean Latitudinal", name, "Diversity"), type = "l", col = "black", ylim = c(0,35))
  lines(x[[2]], col = "blue")
  lines(x[[3]], col = "red")
  lines(x[[4]], col = "green")
  lines(x[[5]], col = "orange")
  legend("topright", lwd=1, legend = c("Surface", "200m", "1000m", "2D", "3D"), col=c("black","blue", "red", "green", "orange"))
}

bGLMsurf <- latitudinalMean(beloniformesGLM[[1]])
colnames(bGLMsurf)[2] <- "Species Richness" 
bGLM200 <- latitudinalMean(beloniformesGLM[["200"]])
colnames(bGLM200)[2] <- "Species Richness" 
bGLM1000 <- latitudinalMean(beloniformesGLM[["1000"]])
colnames(bGLM1000)[2] <- "Species Richness"
bGLM2d <- latitudinalMean(belon2DGLM)
colnames(bGLM2d)[2] <- "Species Richness"
bGLM3dFlat <- latitudinalMean(bGLM3dFlat)
colnames(bGLM3dFlat)[2] <- "Species Richness"

pdf("~/Dropbox/MARDIGRA/Papers/beloniformesDiversityVsLatitudeGLM.pdf")
x <- list(bGLMsurf, bGLM200, bGLM1000, bGLM2d, bGLM3dFlat)
linePlotFunc(x, name = "Beloniformes")
dev.off()

sGLMsurf <- latitudinalMean(scombriformesGLM[[1]])
colnames(sGLMsurf)[2] <- "Species Richness" 
sGLM200 <- latitudinalMean(scombriformesGLM[["200"]])
colnames(sGLM200)[2] <- "Species Richness" 
sGLM1000 <- latitudinalMean(scombriformesGLM[["1000"]])
colnames(sGLM1000)[2] <- "Species Richness"
sGLM2d <- latitudinalMean(scomb2DGLM)
colnames(sGLM2d)[2] <- "Species Richness"
sGLM3dFlat <- latitudinalMean(sGLM3dFlat)
colnames(sGLM3dFlat)[2] <- "Species Richness"

pdf("~/Dropbox/MARDIGRA/Papers/scombriformesDiversityVsLatitudeGLM.pdf")
x <- list(sGLMsurf, sGLM200, sGLM1000, sGLM2d, sGLM3dFlat)
linePlotFunc(x, name = "Scombriformes")
dev.off()

gGLMsurf <- latitudinalMean(gadiformesGLM[[1]])
colnames(gGLMsurf)[2] <- "Species Richness" 
gGLM200 <- latitudinalMean(gadiformesGLM[["200"]])
colnames(gGLM200)[2] <- "Species Richness" 
gGLM1000 <- latitudinalMean(gadiformesGLM[["1000"]])
colnames(gGLM1000)[2] <- "Species Richness"
gGLM2d <- latitudinalMean(gadiformesGLM2D)
colnames(gGLM2d)[2] <- "Species Richness"
gGLM3dFlat <- latitudinalMean(gGLM3dFlat)
colnames(gGLM3dFlat)[2] <- "Species Richness"

pdf("~/Dropbox/MARDIGRA/Papers/gadiformesDiversityVsLatitudeGLM.pdf")
x <- list(gGLMsurf, gGLM200, gGLM1000, gGLM2d, gGLM3dFlat)
linePlotFunc(x, name = "Gadiformes")
dev.off()

bEnvsurf <- latitudinalMean(beloniformesEnv[[1]])
colnames(bEnvsurf)[2] <- "Species Richness" 
bEnv200 <- latitudinalMean(beloniformesEnv[["200"]])
colnames(bEnv200)[2] <- "Species Richness" 
bEnv1000 <- latitudinalMean(beloniformesEnv[["1000"]])
colnames(bEnv1000)[2] <- "Species Richness"
bEnv2d <- latitudinalMean(beloniformesEnv2D)
colnames(bEnv2d)[2] <- "Species Richness"
bEnv3dFlat <- latitudinalMean(bEnv3DFlat)
colnames(bEnv3dFlat)[2] <- "Species Richness"

pdf("~/Dropbox/MARDIGRA/Papers/beloniformesDiversityVsLatitudeEnv.pdf")
x <- list(bEnvsurf, bEnv200, bEnv1000, bEnv2d, bEnv3dFlat)
linePlotFunc(x, name = "Beloniformes")
dev.off()

sEnvsurf <- latitudinalMean(scombriformesEnv[[1]])
colnames(sEnvsurf)[2] <- "Species Richness" 
sEnv200 <- latitudinalMean(scombriformesEnv[["200"]])
colnames(sEnv200)[2] <- "Species Richness" 
sEnv1000 <- latitudinalMean(scombriformesEnv[["1000"]])
colnames(sEnv1000)[2] <- "Species Richness"
sEnv2d <- latitudinalMean(scombriformesEnv2D)
colnames(sEnv2d)[2] <- "Species Richness"
sEnv3dFlat <- latitudinalMean(sEnv3dFlat)
colnames(sEnv3dFlat)[2] <- "Species Richness"

pdf("~/Dropbox/MARDIGRA/Papers/scombriformesDiversityVsLatitudeEnv.pdf")
x <- list(sEnvsurf, sEnv200, sEnv1000, sEnv2d, sEnv3dFlat)
linePlotFunc(x, name = "Scombriformes")
dev.off()

gEnvsurf <- latitudinalMean(gadiformesEnv[[1]])
colnames(gEnvsurf)[2] <- "Species Richness" 
gEnv200 <- latitudinalMean(gadiformesEnv[["200"]])
colnames(gEnv200)[2] <- "Species Richness" 
gEnv1000 <- latitudinalMean(gadiformesEnv[["1000"]])
colnames(gEnv1000)[2] <- "Species Richness"
gEnv2d <- latitudinalMean(gadiformesEnv2D)
colnames(gEnv2d)[2] <- "Species Richness"
gEnv3dFlat <- latitudinalMean(gEnv3dFlat)
colnames(gEnv3dFlat)[2] <- "Species Richness"

pdf("~/Dropbox/MARDIGRA/Papers/gadiformesDiversityVsLatitudeEnv.pdf")
x <- list(gEnvsurf, gEnv200, gEnv1000, gEnv2d, gEnv3dFlat)
linePlotFunc(x, name = "Gadiformes")
dev.off()

# Simple 2D vs 3D model ----
# Load 2D data
beloniformesEnv2D <- crop(rast("data/DiversityEstimates/BeloniformesEnvelope2D.tif"),
                          atlanticShapefile)
beloniformesGLM2D <- crop(rast("data/DiversityEstimates/BeloniformesGLM2D.tif"),
                          atlanticShapefile)
scombriformesEnv2D <- crop(rast("data/DiversityEstimates/ScombriformesEnvelope2D.tif"),
                           atlanticShapefile)
scombriformesGLM2D <- crop(rast("data/DiversityEstimates/ScombriformesGLM2D.tif"),
                           atlanticShapefile)
gadiformesEnv2D <- crop(rast("data/DiversityEstimates/GadiformesEnvelope2D.tif"),
                        atlanticShapefile)
gadiformesGLM2D <- crop(rast("data/DiversityEstimates/GadiformesGLM2D.tif"),
                        atlanticShapefile)

# Generate samples
dummyOccs <- data.frame(0, 0, 0)
colnames(dummyOccs) <- c("longitude", "latitude", "depth")

grid2D <- mSampling2D(occs = dummyOccs, 
                      rasterTemplate = nitrate[[1]], 
                      mShp = atlanticShapefile)

grid3D <- mSampling3D(occs = dummyOccs, 
                      envBrick = nitrate, 
                      mShp = atlanticShapefile)

# Dataset assembly
data3D <- cbind(grid3D,
                xyzSample(grid3D, envBrick = nitrate, verbose = F),
                xyzSample(grid3D, envBrick = temperature, verbose = F),
                xyzSample(grid3D, envBrick = beloniformesEnv, verbose = F),
                xyzSample(grid3D, envBrick = beloniformesGLM, verbose = F),
                xyzSample(grid3D, envBrick = scombriformesEnv, verbose = F),
                xyzSample(grid3D, envBrick = scombriformesGLM, verbose = F),
                xyzSample(grid3D, envBrick = gadiformesEnv, verbose = F),
                xyzSample(grid3D, envBrick = gadiformesGLM, verbose = F))
colnames(data3D) <- c("longitude", "latitude", "depth",
                       "nitrate",
                       "temperature",
                       "belonEnv", "belonGLM", 
                       "scombEnv", "scombGLM", 
                       "gadEnv", "gadGLM")
data3D <- data3D[complete.cases(data3D),]
data3D$latitude <- abs(data3D$latitude)

data2D <- cbind(grid2D,
                terra::extract(ID = FALSE, y = grid2D, x = nitrate[[1]]),
                terra::extract(ID = FALSE, y = grid2D, x = temperature[[1]]),
                terra::extract(ID = FALSE, y = grid2D, x = beloniformesEnv2D),
                terra::extract(ID = FALSE, y = grid2D, x = beloniformesGLM2D),
                terra::extract(ID = FALSE, y = grid2D, x = scombriformesEnv2D),
                terra::extract(ID = FALSE, y = grid2D, x = scombriformesGLM2D),
                terra::extract(ID = FALSE, y = grid2D, x = gadiformesEnv2D),
                terra::extract(ID = FALSE, y = grid2D, x = gadiformesGLM2D))
colnames(data2D) <- c("longitude", "latitude",
                      "nitrate", "temperature",
                      "belonEnv", "belonGLM", 
                      "scombEnv", "scombGLM", 
                      "gadEnv", "gadGLM")
data2D <- data2D[complete.cases(data2D),]
data2D$latitude <- abs(data2D$latitude)

# Model comparisons GLM ----
belonGLMData2D <- subset(data2D, select = -c(belonEnv, scombEnv, scombGLM, gadEnv, gadGLM, 
                                             longitude))

belonGLMData3D <- subset(data3D, select = -c(belonEnv, scombEnv, scombGLM, gadEnv, gadGLM, 
                                             longitude))
# Remove irrelevant depths
belonGLMData3D <- belonGLMData3D %>% 
  filter(depth > min(belonGLMData3D$depth[belonGLMData3D$belonGLM > 0]) & 
           depth < max(belonGLMData3D$depth[belonGLMData3D$belonGLM > 0])) %>% 
  subset(select = -depth)

minSampB <- min(nrow(belonGLMData2D), nrow(belonGLMData3D), 5000)
belonGLMData2D <- belonGLMData2D[sample(1:nrow(belonGLMData2D),
                                size = minSampB,
                                replace = FALSE),]
belonGLMData3D <- belonGLMData3D[sample(1:nrow(belonGLMData3D),
                                        size = minSampB,
                                        replace = FALSE),]

bGLMmodelOutput2D <- modelFunction(belonGLMData2D)
bGLMmodelOutput3D <- modelFunction(belonGLMData3D)

scombGLMData2D <- subset(data2D, select = -c(belonEnv, belonGLM, scombEnv, gadEnv, gadGLM, 
                                             longitude))

scombGLMData3D <- subset(data3D, select = -c(belonEnv, belonGLM, scombEnv, gadEnv, gadGLM, 
                                             longitude))
scombGLMData3D <- scombGLMData3D %>% 
  filter(depth > min(scombGLMData3D$depth[scombGLMData3D$scombGLM > 0]) & 
           depth < max(scombGLMData3D$depth[scombGLMData3D$scombGLM > 0])) %>% 
  subset(select = -depth)

minSampS <- min(nrow(scombGLMData2D), nrow(scombGLMData3D), 5000)
scombGLMData2D <- scombGLMData2D[sample(1:nrow(scombGLMData2D),
                                        size = minSampS,
                                        replace = FALSE),]
scombGLMData3D <- scombGLMData3D[sample(1:nrow(scombGLMData3D),
                                        size = minSampS,
                                        replace = FALSE),]

sGLMmodelOutput2D <- modelFunction(scombGLMData2D)
sGLMmodelOutput3D <- modelFunction(scombGLMData3D)

gadGLMData2D <- subset(data2D, select = -c(belonEnv, belonGLM, scombEnv, gadEnv, scombGLM,
                                           longitude))

gadGLMData3D <- subset(data3D, select = -c(belonEnv, belonGLM, scombEnv, gadEnv, scombGLM, 
                                           longitude))
gadGLMData3D <- gadGLMData3D %>% 
  filter(depth > min(gadGLMData3D$depth[gadGLMData3D$gadGLM > 0]) & 
           depth < max(gadGLMData3D$depth[gadGLMData3D$gadGLM > 0])) %>% 
  subset(select = -depth)

minSampG <- min(nrow(gadGLMData2D), nrow(gadGLMData3D), 5000)
gadGLMData2D <- gadGLMData2D[sample(1:nrow(gadGLMData2D),
                                    size = minSampG,
                                    replace = FALSE),]
gadGLMData3D <- gadGLMData3D[sample(1:nrow(gadGLMData3D),
                                    size = minSampG,
                                    replace = FALSE),]

gGLMmodelOutput2D <- modelFunction(gadGLMData2D)
gGLMmodelOutput3D <- modelFunction(gadGLMData3D)

# Results
resultsTable <- data.frame(row.names = c("bG2D", "bG3D", 
                                         "sG2D", "sG3D", 
                                         "gG2D", "gG3D"))
resultsTable$R2 <- c(bGLMmodelOutput2D$performance$R2_Nagelkerke,
                     bGLMmodelOutput3D$performance$R2_Nagelkerke,
                     sGLMmodelOutput2D$performance$R2_Nagelkerke,
                     sGLMmodelOutput3D$performance$R2_Nagelkerke,
                     gGLMmodelOutput2D$performance$R2_Nagelkerke,
                     gGLMmodelOutput3D$performance$R2_Nagelkerke)
resultsTable$AICc <- c(bGLMmodelOutput2D$performance$AICc, 
                       bGLMmodelOutput3D$performance$AICc,
                       sGLMmodelOutput2D$performance$AICc, 
                       sGLMmodelOutput3D$performance$AICc,
                       gGLMmodelOutput2D$performance$AICc, 
                       gGLMmodelOutput3D$performance$AICc)

varsInModels <- unique(c(rownames(bGLMmodelOutput2D$coeff), 
                         rownames(bGLMmodelOutput3D$coeff), 
                         rownames(sGLMmodelOutput2D$coeff), 
                         rownames(sGLMmodelOutput3D$coeff),
                         rownames(gGLMmodelOutput2D$coeff), 
                         rownames(gGLMmodelOutput3D$coeff)))[-1]
coefTable <- matrix(nrow = nrow(resultsTable), ncol = length(varsInModels) * 2)
row.names(coefTable) <- row.names(resultsTable)
for (i in 1:(ncol(coefTable)/2)){
  variable <- varsInModels[[i]]
  coefTable[1, i] <- ifelse(0 < length(bGLMmodelOutput2D$coeff[rownames(bGLMmodelOutput2D$coeff) %in% variable,1]), 
                            bGLMmodelOutput2D$coeff[rownames(bGLMmodelOutput2D$coeff) %in% variable,1], NA)
  coefTable[1, (i+length(varsInModels))] <- ifelse(0 < length(bGLMmodelOutput2D$coeff[rownames(bGLMmodelOutput2D$coeff) %in% variable,4]), 
                                                   bGLMmodelOutput2D$coeff[rownames(bGLMmodelOutput2D$coeff) %in% variable,4], NA)
  
  coefTable[2, i] <- ifelse(0 < length(bGLMmodelOutput3D$coeff[rownames(bGLMmodelOutput3D$coeff) %in% variable,1]), 
                            bGLMmodelOutput3D$coeff[rownames(bGLMmodelOutput3D$coeff) %in% variable,1], NA)
  coefTable[2, (i+length(varsInModels))] <- ifelse(0 < length(bGLMmodelOutput3D$coeff[rownames(bGLMmodelOutput3D$coeff) %in% variable,4]), 
                                                   bGLMmodelOutput3D$coeff[rownames(bGLMmodelOutput3D$coeff) %in% variable,4], NA)
  
  coefTable[3, i] <- ifelse(0 < length(sGLMmodelOutput2D$coeff[rownames(sGLMmodelOutput2D$coeff) %in% variable,1]), 
                            sGLMmodelOutput2D$coeff[rownames(sGLMmodelOutput2D$coeff) %in% variable,1], NA)
  coefTable[3, (i+length(varsInModels))] <- ifelse(0 < length(sGLMmodelOutput2D$coeff[rownames(sGLMmodelOutput2D$coeff) %in% variable,4]), 
                                                   sGLMmodelOutput2D$coeff[rownames(sGLMmodelOutput2D$coeff) %in% variable,4], NA)
  
  coefTable[4, i] <- ifelse(0 < length(sGLMmodelOutput3D$coeff[rownames(sGLMmodelOutput3D$coeff) %in% variable,1]), 
                            sGLMmodelOutput3D$coeff[rownames(sGLMmodelOutput3D$coeff) %in% variable,1], NA)
  coefTable[4, (i+length(varsInModels))] <- ifelse(0 < length(sGLMmodelOutput3D$coeff[rownames(sGLMmodelOutput3D$coeff) %in% variable,4]), 
                                                   sGLMmodelOutput3D$coeff[rownames(sGLMmodelOutput3D$coeff) %in% variable,4], NA)
  
  coefTable[5, i] <- ifelse(0 < length(gGLMmodelOutput2D$coeff[rownames(gGLMmodelOutput2D$coeff) %in% variable,1]), 
                            gGLMmodelOutput2D$coeff[rownames(gGLMmodelOutput2D$coeff) %in% variable,1], NA)
  coefTable[5, (i+length(varsInModels))] <- ifelse(0 < length(gGLMmodelOutput2D$coeff[rownames(gGLMmodelOutput2D$coeff) %in% variable,4]), 
                                                   gGLMmodelOutput2D$coeff[rownames(gGLMmodelOutput2D$coeff) %in% variable,4], NA)
  
  coefTable[6, i] <- ifelse(0 < length(gGLMmodelOutput3D$coeff[rownames(gGLMmodelOutput3D$coeff) %in% variable,1]), 
                            gGLMmodelOutput3D$coeff[rownames(gGLMmodelOutput3D$coeff) %in% variable,1], NA)
  coefTable[6, (i+length(varsInModels))] <- ifelse(0 < length(gGLMmodelOutput3D$coeff[rownames(gGLMmodelOutput3D$coeff) %in% variable,4]), 
                                                   gGLMmodelOutput3D$coeff[rownames(gGLMmodelOutput3D$coeff) %in% variable,4], NA)
}

colnames(coefTable) <- c(varsInModels, paste0(varsInModels, "_pVal"))
finalTable <- data.frame(resultsTable, coefTable)
finalTable$SampleSize <- c(minSampB, minSampB, minSampS, minSampS, minSampG, minSampG)

write.csv(finalTable, file = "data/linearModelResults_2D3D_GLM.csv", row.names = T)

# Model comparisons Envelope ----
belonEnvData2D <- subset(data2D, select = -c(belonGLM, scombEnv, scombGLM, gadEnv, gadGLM, 
                                             longitude))

belonEnvData3D <- subset(data3D, select = -c(belonGLM, scombEnv, scombGLM, gadEnv, gadGLM, 
                                             longitude))
belonEnvData3D <- belonEnvData3D %>% 
  filter(depth > min(belonEnvData3D$depth[belonEnvData3D$belonEnv > 0]) & 
           depth < max(belonEnvData3D$depth[belonEnvData3D$belonEnv > 0])) %>% 
  subset(select = -depth)

minSampB <- min(nrow(belonEnvData2D), nrow(belonEnvData3D), 5000)
belonEnvData2D <- belonEnvData2D[sample(1:nrow(belonEnvData2D),
                                        size = minSampB,
                                        replace = FALSE),]
belonEnvData3D <- belonEnvData3D[sample(1:nrow(belonEnvData3D),
                                        size = minSampB,
                                        replace = FALSE),]

bEnvmodelOutput2D <- modelFunction(belonEnvData2D)
bEnvmodelOutput3D <- modelFunction(belonEnvData3D)

scombEnvData2D <- subset(data2D, select = -c(belonEnv, belonGLM, scombGLM, gadEnv, gadGLM, 
                                             longitude))

scombEnvData3D <- subset(data3D, select = -c(belonEnv, belonGLM, scombGLM, gadEnv, gadGLM, 
                                             longitude))
scombEnvData3D <- scombEnvData3D %>% 
  filter(depth > min(scombEnvData3D$depth[scombEnvData3D$scombEnv > 0]) & 
           depth < max(scombEnvData3D$depth[scombEnvData3D$scombEnv > 0])) %>% 
  subset(select = -depth)

minSampS <- min(nrow(scombEnvData2D), nrow(scombEnvData3D), 5000)
scombEnvData2D <- scombEnvData2D[sample(1:nrow(scombEnvData2D),
                                        size = minSampS,
                                        replace = FALSE),]
scombEnvData3D <- scombEnvData3D[sample(1:nrow(scombEnvData3D),
                                        size = minSampS,
                                        replace = FALSE),]

sEnvmodelOutput2D <- modelFunction(scombEnvData2D)
sEnvmodelOutput3D <- modelFunction(scombEnvData3D)

gadEnvData2D <- subset(data2D, select = -c(belonEnv, belonGLM, scombEnv, scombGLM, gadGLM, 
                                           longitude))

gadEnvData3D <- subset(data3D, select = -c(belonEnv, belonGLM, scombEnv, scombGLM, gadGLM, 
                                           longitude))
gadEnvData3D <- gadEnvData3D %>% 
  filter(depth > min(gadEnvData3D$depth[gadEnvData3D$gadEnv > 0]) & 
           depth < max(gadEnvData3D$depth[gadEnvData3D$gadEnv > 0])) %>% 
  subset(select = -depth)

minSampG <- min(nrow(gadEnvData2D), nrow(gadEnvData3D), 5000)
gadEnvData2D <- gadEnvData2D[sample(1:nrow(gadEnvData2D),
                                    size = minSampG,
                                    replace = FALSE),]
gadEnvData3D <- gadEnvData3D[sample(1:nrow(gadEnvData3D),
                                    size = minSampG,
                                    replace = FALSE),]

gEnvmodelOutput2D <- modelFunction(gadEnvData2D)
gEnvmodelOutput3D <- modelFunction(gadEnvData3D)

# Results
resultsTable <- data.frame(row.names = c("bG2D", "bG3D", 
                                         "sG2D", "sG3D", 
                                         "gG2D", "gG3D"))
resultsTable$R2 <- c(bEnvmodelOutput2D$performance$R2_Nagelkerke,
                     bEnvmodelOutput3D$performance$R2_Nagelkerke,
                     sEnvmodelOutput2D$performance$R2_Nagelkerke,
                     sEnvmodelOutput3D$performance$R2_Nagelkerke,
                     gEnvmodelOutput2D$performance$R2_Nagelkerke,
                     gEnvmodelOutput3D$performance$R2_Nagelkerke)
resultsTable$AICc <- c(bEnvmodelOutput2D$performance$AICc, 
                       bEnvmodelOutput3D$performance$AICc,
                       sEnvmodelOutput2D$performance$AICc, 
                       sEnvmodelOutput3D$performance$AICc,
                       gEnvmodelOutput2D$performance$AICc, 
                       gEnvmodelOutput3D$performance$AICc)

varsInModels <- unique(c(rownames(bEnvmodelOutput2D$coeff), 
                         rownames(bEnvmodelOutput3D$coeff), 
                         rownames(sEnvmodelOutput2D$coeff), 
                         rownames(sEnvmodelOutput3D$coeff),
                         rownames(gEnvmodelOutput2D$coeff), 
                         rownames(gEnvmodelOutput3D$coeff)))[-1]
coefTable <- matrix(nrow = nrow(resultsTable), ncol = length(varsInModels) * 2)
row.names(coefTable) <- row.names(resultsTable)
for (i in 1:(ncol(coefTable)/2)){
  variable <- varsInModels[[i]]
  coefTable[1, i] <- ifelse(0 < length(bEnvmodelOutput2D$coeff[rownames(bEnvmodelOutput2D$coeff) %in% variable,1]), 
                            bEnvmodelOutput2D$coeff[rownames(bEnvmodelOutput2D$coeff) %in% variable,1], NA)
  coefTable[1, (i+length(varsInModels))] <- ifelse(0 < length(bEnvmodelOutput2D$coeff[rownames(bEnvmodelOutput2D$coeff) %in% variable,4]), 
                                                   bEnvmodelOutput2D$coeff[rownames(bEnvmodelOutput2D$coeff) %in% variable,4], NA)
  
  coefTable[2, i] <- ifelse(0 < length(bEnvmodelOutput3D$coeff[rownames(bEnvmodelOutput3D$coeff) %in% variable,1]), 
                            bEnvmodelOutput3D$coeff[rownames(bEnvmodelOutput3D$coeff) %in% variable,1], NA)
  coefTable[2, (i+length(varsInModels))] <- ifelse(0 < length(bEnvmodelOutput3D$coeff[rownames(bEnvmodelOutput3D$coeff) %in% variable,4]), 
                                                   bEnvmodelOutput3D$coeff[rownames(bEnvmodelOutput3D$coeff) %in% variable,4], NA)
  
  coefTable[3, i] <- ifelse(0 < length(sEnvmodelOutput2D$coeff[rownames(sEnvmodelOutput2D$coeff) %in% variable,1]), 
                            sEnvmodelOutput2D$coeff[rownames(sEnvmodelOutput2D$coeff) %in% variable,1], NA)
  coefTable[3, (i+length(varsInModels))] <- ifelse(0 < length(sEnvmodelOutput2D$coeff[rownames(sEnvmodelOutput2D$coeff) %in% variable,4]), 
                                                   sEnvmodelOutput2D$coeff[rownames(sEnvmodelOutput2D$coeff) %in% variable,4], NA)
  
  coefTable[4, i] <- ifelse(0 < length(sEnvmodelOutput3D$coeff[rownames(sEnvmodelOutput3D$coeff) %in% variable,1]), 
                            sEnvmodelOutput3D$coeff[rownames(sEnvmodelOutput3D$coeff) %in% variable,1], NA)
  coefTable[4, (i+length(varsInModels))] <- ifelse(0 < length(sEnvmodelOutput3D$coeff[rownames(sEnvmodelOutput3D$coeff) %in% variable,4]), 
                                                   sEnvmodelOutput3D$coeff[rownames(sEnvmodelOutput3D$coeff) %in% variable,4], NA)
  
  coefTable[5, i] <- ifelse(0 < length(gEnvmodelOutput2D$coeff[rownames(gEnvmodelOutput2D$coeff) %in% variable,1]), 
                            gEnvmodelOutput2D$coeff[rownames(gEnvmodelOutput2D$coeff) %in% variable,1], NA)
  coefTable[5, (i+length(varsInModels))] <- ifelse(0 < length(gEnvmodelOutput2D$coeff[rownames(gEnvmodelOutput2D$coeff) %in% variable,4]), 
                                                   gEnvmodelOutput2D$coeff[rownames(gEnvmodelOutput2D$coeff) %in% variable,4], NA)
  
  coefTable[6, i] <- ifelse(0 < length(gEnvmodelOutput3D$coeff[rownames(gEnvmodelOutput3D$coeff) %in% variable,1]), 
                            gEnvmodelOutput3D$coeff[rownames(gEnvmodelOutput3D$coeff) %in% variable,1], NA)
  coefTable[6, (i+length(varsInModels))] <- ifelse(0 < length(gEnvmodelOutput3D$coeff[rownames(gEnvmodelOutput3D$coeff) %in% variable,4]), 
                                                   gEnvmodelOutput3D$coeff[rownames(gEnvmodelOutput3D$coeff) %in% variable,4], NA)
}

colnames(coefTable) <- c(varsInModels, paste0(varsInModels, "_pVal"))
finalTable <- data.frame(resultsTable, coefTable)
finalTable$SampleSize <- c(minSampB, minSampB, minSampS, minSampS, minSampG, minSampG)

write.csv(finalTable, file = "data/linearModelResults_2D3D_Env.csv", row.names = T)

# Plotting Coefficients ----
bPlot <- list("3D Envelope" = bEnvmodelOutput3D$model, "2D Envelope" = bEnvmodelOutput2D$model,
              "3D GLM" = bGLMmodelOutput3D$model, "2D GLM" = bGLMmodelOutput2D$model) %>% 
  dwplot(dot_args = list(aes(shape = model, 
                             color = cut(p.value, 
                                         breaks = c(0, 0.01, 0.05, 1)), 
                             size = 2)),
         whisker_args = list(aes(color = cut(p.value, 
                                             breaks = c(0, 0.01, 0.05, 1))))) %>%  
  relabel_predictors(c(latitude = "Latitude",
                       nitrate = "Nitrate",
                       temperature = "Temperature")) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + 
  xlim(-0.15,0.25) + ggtitle("Beloniformes") +
  xlab("Coefficient Estimate") + 
  scale_shape_manual(name = "Model", values = c(2,1,17,16)) +
  scale_color_manual(name = "P value",
                     labels = c("< 0.01", "0.05 - 0.01", "> 0.05"),
                     values = c("#a50f15", "#fb6a4a", "gray75")) +
  theme(plot.margin = unit(c(0,-.5,0,.5), "lines"),
        plot.background = element_blank()) +
  guides(size = "none")

sPlot <- list("3D Envelope" = sEnvmodelOutput3D$model, "2D Envelope" = sEnvmodelOutput2D$model,
              "3D GLM" = sGLMmodelOutput3D$model, "2D GLM" = sGLMmodelOutput2D$model) %>% 
  dwplot(dot_args = list(aes(shape = model, 
                             color = cut(p.value, 
                                         breaks = c(0, 0.01, 0.05, 1)), 
                             size = 2)),
         whisker_args = list(aes(color = cut(p.value, 
                                             breaks = c(0, 0.01, 0.05, 1))))) %>%  
  relabel_predictors(c(latitude = "Latitude",
                       nitrate = "Nitrate",
                       temperature = "Temperature")) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + 
  xlim(-0.15,0.25) + ggtitle("Scombriformes") +
  xlab("Coefficient Estimate") + 
  scale_shape_manual(name = "Model", values = c(2,1,17,16)) +
  scale_color_manual(name = "P value",
                     labels = c("< 0.01", "0.05 - 0.01", "> 0.05"),
                     values = c("#a50f15", "#fb6a4a", "gray75")) +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        plot.margin = unit(c(0,1,0,-1), "lines"),
        plot.background = element_blank()) +
  guides(size = "none")

gPlot <- list("3D Envelope" = gEnvmodelOutput3D$model, "2D Envelope" = gEnvmodelOutput2D$model,
              "3D GLM" = gGLMmodelOutput3D$model, "2D GLM" = gGLMmodelOutput2D$model) %>% 
  dwplot(dot_args = list(aes(shape = model, 
                             color = cut(p.value, 
                                         breaks = c(0, 0.01, 0.05, 1)), 
                             size = 2)),
         whisker_args = list(aes(color = cut(p.value, 
                                             breaks = c(0, 0.01, 0.05, 1))))) %>%  
  relabel_predictors(c(latitude = "Latitude",
                       nitrate = "Nitrate",
                       temperature = "Temperature")) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + 
  xlim(-0.15,0.25) + ggtitle("Gadiformes") +
  xlab("Coefficient Estimate") + 
  scale_shape_manual(name = "Model", values = c(2,1,17,16)) +
  scale_color_manual(name = "P value",
                     labels = c("< 0.01", "0.05 - 0.01", "> 0.05"),
                     values = c("#a50f15", "#fb6a4a", "gray75")) +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        plot.margin = unit(c(0,1,0,-1), "lines"),
        plot.background = element_blank()) +
  guides(size = "none")

allPlot <- ggarrange(bPlot, sPlot, gPlot, common.legend = F, legend = "right", 
                     ncol = 3, widths = c(2,1,1))
ggsave("~/Dropbox/MARDIGRA/Papers/Manuscript/TablesAndFigures/BioDivModels.pdf",
       allPlot, height = 6, width = 11.4, units = "cm", scale = 2)
