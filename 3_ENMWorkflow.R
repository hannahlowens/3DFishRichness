library(voluModel) 
library(terra)
library(sf)
library(fields) 
library(rnaturalearth) #
library(ENMeval) #
library(mecofun) #
library(gridExtra) #
library(ggpubr) #
library(cowplot) #
library(dplyr)
library(glm2)

setwd("~/Dropbox/MARDIGRA/data/")

# Load environmental data ----
oxygen <- scale(rast("EnvironmentalData/ProcessedEnvtData/AOU.tif"))
salinity <- scale(rast("EnvironmentalData/ProcessedEnvtData/salinity.tif"))
temperature <- scale(rast("EnvironmentalData/ProcessedEnvtData/temperature.tif"))

bathymetry <- rast("EnvironmentalData/ETOPO_2022_v1_60s_N90W180_bed.tif")
names(bathymetry) <- "bathymetry"
bathymetry <- aggregate(bathymetry, fact = 60, fun = min)
values(bathymetry)[values(bathymetry) > 0] <- NA
values(bathymetry) <- abs(values(bathymetry))

oceanOnly <- vect("FAO Fishing Areas 2005/FAO_AREAS.shp")
oceanOnly <- aggregate(oceanOnly)
oceanOnly <- project(oceanOnly, bathymetry)
bathymetry <- mask(bathymetry, mask = oceanOnly, touches = FALSE)

bottomOxygen <- scale(rast("EnvironmentalData/ProcessedEnvtData/bottomAOU.tif"))
bottomSalinity <- scale(rast("EnvironmentalData/ProcessedEnvtData/bottomSalinity.tif"))
bottomTemperature <- scale(rast("EnvironmentalData/ProcessedEnvtData/bottomTemperature.tif"))

# Occurrence processing ----
occList <- list.files(path = "FinalOccurrenceDataset/", 
                      pattern = ".csv", full.names = T)

processedOccs <- vector(mode = "list", length = length(occList))
layerNames <- as.numeric(gsub("[X]", "", names(temperature)))
j <- 1
for (sp in occList){
  occurrences <- read.csv(file = sp)
  fn <- gsub(sp, pattern = "FinalOccurrenceDataset/", replacement = "")
  occurrences <- occurrences[,c("decimalLatitude", "decimalLongitude", "depth")] # Remove metadata
  occurrences <- dplyr::distinct(occurrences)
  
  # Filter out putative depth errors
  occurrences$mapDepth <- round(terra::extract(bathymetry,
                      as.matrix(occurrences[,c("decimalLongitude","decimalLatitude")]), method = "bilinear"),3)
  
  noDataOccurrences <- occurrences[apply(is.na(occurrences[,c("depth", "mapDepth")]), MARGIN = 1, FUN = function(x) any(x)),]
  occurrences <- occurrences[!apply(is.na(occurrences[,c("depth", "mapDepth")]), MARGIN = 1, FUN = function(x) any(x)),]
  
  # Get rid of depths that are deeper or at the same depth as bathymetry
  occurrences <- occurrences[(occurrences$depth/
                                occurrences$mapDepth) < 0.95,]
  
  # Get rid of outliers
  if(nrow(occurrences) > 2){
    test <- as.numeric(scale(occurrences$depth))
    test <- between(test, -2, 2)
    test[is.na(test)] <- TRUE
    outCount <- sum(!test)
    if(outCount > 1){
      extremeOutMin <- ifelse((min(occurrences[test, "depth"], 
                                   na.rm = T) > min(occurrences[!test, "depth"], 
                                                  na.rm = T)), ((min(occurrences[test, "depth"], 
                                                                    na.rm = T)-min(occurrences[!test, "depth"], 
                                                                                   na.rm = T))/min(occurrences[test, "depth"], 
                                                                                                   na.rm = T) > 2), FALSE)
      extremeOutMax <- ifelse((max(occurrences[!test, "depth"],
                                   na.rm = T) > max(occurrences[test, "depth"],
                                                    na.rm = T)), ((max(occurrences[!test, "depth"], 
                                                                       na.rm = T) - max(occurrences[test, "depth"], 
                                                                                        na.rm = T))/max(occurrences[test, "depth"], 
                                                                                                        na.rm = T) > 2), FALSE)
      extremeOut <- any(extremeOutMin, extremeOutMax)
      occurrences <- occurrences[test,]
      while(all(extremeOut, outCount > 1)){
        test <- as.numeric(scale(occurrences$depth))
        test <- between(test, -2, 2)
        test[is.na(test)] <- TRUE
        outCount <- sum(!test)
        
        extremeOutMin <- ifelse((min(occurrences[test, "depth"], 
                                     na.rm = T) > min(occurrences[!test, "depth"], 
                                                      na.rm = T)), ((min(occurrences[test, "depth"], 
                                                                         na.rm = T)-min(occurrences[!test, "depth"], 
                                                                                        na.rm = T))/min(occurrences[test, "depth"], 
                                                                                                        na.rm = T) > 2), FALSE)
        extremeOutMax <- ifelse((max(occurrences[!test, "depth"],
                                     na.rm = T) > max(occurrences[test, "depth"],
                                                      na.rm = T)), ((max(occurrences[!test, "depth"], 
                                                                         na.rm = T) - max(occurrences[test, "depth"], 
                                                                                          na.rm = T))/max(occurrences[test, "depth"], 
                                                                                                          na.rm = T) > 2), FALSE)
        extremeOut <- any(extremeOutMin, extremeOutMax)
        occurrences <- occurrences[test,]
      }
    }
  }
  
  # Put everything back together
  occurrences <- rbind(occurrences[,1:4], noDataOccurrences)
  
  # Downsample by depth if more than 20 points
  if(nrow(occurrences) == sum(is.na(occurrences$depth))) { # Downsample to XY resolution only
    tempPoints <- occurrences
    tempPoints <- downsample(tempPoints, temperature[[1]], verbose = F)
    tempPoints$depth <- rep(NA, times = nrow(tempPoints))
    occurrences <- tempPoints
  }else if (nrow(occurrences) > 20){
    occurrences <- na.omit(occurrences)
    occurrences$index <- unlist(lapply(occurrences$depth, 
                                       FUN = function(x) which.min(abs(layerNames - x))))
    indices <- unique(occurrences$index)
    downsampledOccs <- data.frame()
    if (length(indices > 0)){
      for(i in indices){
        tempPoints <- occurrences[occurrences$index==i,]
        tempPoints <- downsample(tempPoints, temperature[[1]], verbose = F)
        tempPoints$depth <- rep(layerNames[[i]], times = nrow(tempPoints))
        downsampledOccs <- rbind(downsampledOccs, tempPoints)
      }
    }
    occurrences <- downsampledOccs
  } else { # Downsample to XY resolution only
    tempPoints <- occurrences
    tempPoints <- downsample(tempPoints, temperature[[1]], verbose = F)
    tempPoints$depth <- rep(NA, times = nrow(tempPoints))
    occurrences <- tempPoints
  }
  print(paste0(fn, "; number of occs: ", nrow(occurrences)))
  processedOccs[[j]] <- occurrences
  write.csv(processedOccs[[j]],
            file = paste0("OccurrenceDatasetCleanForModeling/", fn),
            row.names = F)
  j <- j + 1
}

rm(fn, i, indices, j, layerNames, sp, tempPoints, 
   noDataOccurrences, downsampledOccs, extremeOut, 
   extremeOutMax, extremeOutMin, test, outCount)

occList <- list.files(path = "OccurrenceDatasetCleanForModeling/", 
                      pattern = ".csv", full.names = T)
spNames <- gsub(".csv", replacement = "",
                x = gsub("OccurrenceDatasetCleanForModeling//",
                         replacement = "", x = occList))
occurrenceList <- vector(mode = "list", length = length(occList))
for(sp in occList){
  occurrences <- read.csv(sp)
  occurrenceList[[match(sp, occList)]] <- occurrences
}
names(occurrenceList) <- spNames
rm(occurrences, sp)

# Get 2D environmental data and make envelope buffers ----
# Set up vectors to fill
environmentalData2Dsurf <- c(oxygen[[1]], salinity[[1]], temperature[[1]])
environmentalData2Dbott <- c(bottomOxygen, bottomSalinity, bottomTemperature)
names(environmentalData2Dsurf) <- c("Oxygen", "Salinity", "Temperature")
names(environmentalData2Dbott) <- c("Oxygen", "Salinity", "Temperature")
occsWdata2D <- vector(mode = "list", length = length(occurrenceList))
model2DisBottom <- rep_len(x = F, length.out = length(occurrenceList))

# Which 2D are bottom?
for (i in 1:length(occurrenceList)){
  print(paste0(spNames[[i]]))
  medianDepth <- median(occurrenceList[[i]]$depth)
  if(medianDepth < 200 || is.na(medianDepth)){
    vals <- cbind(occurrenceList[[i]][,c("decimalLongitude", "decimalLatitude")],
                  raster::extract(x = environmentalData2Dsurf,
                                  y = occurrenceList[[i]][,c("decimalLongitude", "decimalLatitude")]))
  } else {
    vals <- cbind(occurrenceList[[i]][,c("decimalLongitude", "decimalLatitude")],
                  raster::extract(x = environmentalData2Dbott,
                                  y = occurrenceList[[i]][,c("decimalLongitude", "decimalLatitude")]))
    model2DisBottom[[i]] <- T
  }
}

rm(medianDepth, vals)

# Get 2D data and calculate M regions
for (i in 1:length(occurrenceList)){
  print(paste0(spNames[[i]]))
  if(model2DisBottom[[i]]){
    vals <- cbind(occurrenceList[[i]][,c("decimalLongitude", "decimalLatitude")],
                  raster::extract(x = environmentalData2Dbott,
                          y = occurrenceList[[i]][,c("decimalLongitude", "decimalLatitude")]))
  } else {
    vals <- cbind(occurrenceList[[i]][,c("decimalLongitude", "decimalLatitude")],
                  raster::extract(x = environmentalData2Dsurf,
                          y = occurrenceList[[i]][,c("decimalLongitude", "decimalLatitude")]))
  }

  vals <- vals[complete.cases(vals),]
  vals <- unique(vals)
  # Remove environmental outliers if there's enough data to do so
  if(!is.null(dim(vals))){
    if(nrow(vals) > 5){
      if(nrow(vals) > 50){
        valsMD <- as.matrix(vals[, c("Oxygen", "Temperature", "Salinity")])
        md <- mahalanobis(x = valsMD, center = colMeans(valsMD), cov = cov(valsMD))
        outs <- boxplot.stats(md)$out
        vals <- vals[!(md %in% outs),]
      }
      occsWdata2D[[i]] <- data.frame(vals)
      # Calculate M
      writeVector(x = marineBackground(occsWdata2D[[i]]), 
                  filetype = "ESRI Shapefile",
                  filename = paste0("TrainingRegions/2D/", spNames[[i]],".shp"), 
                  overwrite = TRUE) 
    } else{
      occsWdata2D[[i]] <- NA
    }
  } else {
    occsWdata2D[[i]] <- NA
  }
}
names(occsWdata2D) <- spNames

rm(md, outs, vals, valsMD, i)

# Get 3D environmental data 
occsWdata3D <- vector(mode = "list", length = length(occurrenceList))
for (i in 1:length(occurrenceList)){
  print(paste0(spNames[[i]]))
  occs <- occurrenceList[[i]]
  if(all(is.na(occs$depth))){
    occsWdata3D[[i]] <- NA
  } else {
    oxyVals <- xyzSample(occs = occs, oxygen)
    tempVals <- xyzSample(occs = occs, temperature)
    salVals <- xyzSample(occs = occs, salinity)
    vals <- cbind(occs, oxyVals, tempVals, salVals)
    colnames(vals) <- c("decimalLongitude", "decimalLatitude", "depth", "Oxygen", "Temperature", "Salinity")
    vals <- vals[complete.cases(vals),]
    row.names(vals) <- NULL
    if(!is.null(dim(vals))){
      if(nrow(vals) > 50){
        valsMD <- as.matrix(vals[, c("Oxygen", "Temperature", "Salinity")])
        md <- mahalanobis(x = valsMD, center = colMeans(valsMD), cov = cov(valsMD))
        outs <- boxplot.stats(md)$out
        vals <- vals[!(md %in% outs),]
      }
    }
    occsWdata3D[[i]] <- data.frame(vals)
    
    # Create training region
    if(nrow(occsWdata3D[[i]]) > 5){
      writeVector(x = marineBackground(occsWdata3D[[i]]), 
                  filetype = "ESRI Shapefile",
                  filename = paste0("TrainingRegions/3D/", spNames[[i]],".shp"),
                  overwrite = TRUE)
    } else{
      occsWdata3D[[i]] <- NA
    }
  }
}
names(occsWdata3D) <- spNames

# Map training regions
land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")[1]
trainingRegion2D <- list.files(path = "TrainingRegions/2D/", 
                               full.names = TRUE)
trainingRegion3D <- list.files(path = "TrainingRegions/3D/", 
                               full.names = TRUE)
pdf("trainingRegions_1Aug.pdf")
par(mfrow = c(2, 1))
for(i in 1:length(occurrenceList)){
  if(any(grepl(pattern = paste0(spNames[[i]], ".shp"), x = trainingRegion2D, ignore.case = TRUE))){
    plot(bathymetry, 
         main = paste0("2D Training Region and Occurrences, ", spNames[[i]]), 
         col = viridis(n= 11, option = "mako", direction = -1), smooth = T)
    plot(vect(trainingRegion2D[[which(grepl(pattern = paste0(spNames[[i]], ".shp"), 
                                            x = trainingRegion2D, ignore.case = TRUE))]]), 
         add = T, border = "orange", lwd = 2)
    points(occsWdata2D[[spNames[[i]]]][,c("decimalLongitude","decimalLatitude")], 
           cex = 1, pch = 20, col = "red")
  } else {
    plot(0,0, main = paste0("No 2D Training Region for ", spNames[[i]]))
  }
  if(any(grepl(pattern = paste0(spNames[[i]], ".shp"), x = trainingRegion3D, ignore.case = TRUE))){
    plot(bathymetry, 
         main = paste0("3D Training Region and Occurrences, ", spNames[[i]]), 
         col = viridis(n= 11, option = "mako", direction = -1), smooth = T)
    plot(vect(trainingRegion3D[[which(grepl(pattern = paste0(spNames[[i]], ".shp"), 
                                            x = trainingRegion3D, ignore.case = TRUE))]]), 
         add = T, border = "orange", lwd = 2)
    points(occsWdata3D[[spNames[[i]]]][,c("decimalLongitude","decimalLatitude")], 
           cex = 1, pch = 20, col = "red")
  } else {
    plot(0,0, main = paste0("No 3D Training Region for ", spNames[[i]]))
  }
}
dev.off()

# 2D envelope models ----
for(i in 1:length(occsWdata2D)){
  if(!is.na(occsWdata2D)[[spNames[[i]]]]){
    vals <- occsWdata2D[[spNames[[i]]]]
    if(any(grepl(pattern = paste0(spNames[[i]], ".shp"), x = trainingRegion2D, ignore.case = TRUE))){
      if(length(unlist(vals)) > 7){
        oxyLims <- quantile(vals$Oxygen,c(.01, .99))
        oxyLims <- matrix(c(-Inf,oxyLims[[1]],0,
                            oxyLims[[1]], oxyLims[[2]], 1,
                            oxyLims[[2]], Inf, 0), 
                          nrow = 3, byrow = TRUE)
        tempLims <- quantile(vals$Temperature,c(.01, .99))
        tempLims <- matrix(c(-Inf,tempLims[[1]],0,
                            tempLims[[1]], tempLims[[2]], 1,
                            tempLims[[2]], Inf, 0), 
                          nrow = 3, byrow = TRUE)
        salLims <- quantile(vals$Salinity,c(.01, .99))
        salLims <- matrix(c(-Inf,salLims[[1]],0,
                            salLims[[1]], salLims[[2]], 1,
                            salLims[[2]], Inf, 0), 
                          nrow = 3, byrow = TRUE)
        
        envBuffer2D <- vect(trainingRegion2D[[which(grepl(pattern = paste0(spNames[[i]], ".shp"), 
                                                          x = trainingRegion2D, ignore.case = TRUE))]])
        
        if(model2DisBottom[[i]]){
          dist_env_2d <- terra::mask(x = terra::crop(x = environmentalData2Dbott, 
                                       envBuffer2D), 
                              mask = envBuffer2D)
        } else{
          dist_env_2d <- terra::mask(x = crop(x = environmentalData2Dsurf, envBuffer2D), mask = envBuffer2D)
        }
        
        dist_env_2d$Oxygen <- classify(dist_env_2d$Oxygen, rcl = oxyLims)
        dist_env_2d$Temperature <- classify(dist_env_2d$Temperature, rcl = tempLims)
        dist_env_2d$Salinity <- classify(dist_env_2d$Salinity, rcl = salLims)
        dist_env_2d <- dist_env_2d$Oxygen*dist_env_2d$Salinity*dist_env_2d$Temperature
        dist_env_2d <- classify(dist_env_2d, matrix(c(NA, 0, 1, 1), nrow = 2, byrow = TRUE), include.lowest = T)
        terra::writeRaster(dist_env_2d, 
                           filename = paste0("Models/Envelope_2D/", spNames[[i]], ".asc"), 
                           overwrite = T)
      }
    }
  }  
}

# 3D envelope models ----
for(i in 1:length(occsWdata3D)){
  if(!is.na(occsWdata3D)[[spNames[[i]]]]){
    if(any(grepl(pattern = paste0(spNames[[i]], ".shp"), x = trainingRegion3D, ignore.case = TRUE))){
      vals <- occsWdata3D[[spNames[[i]]]]
      if(length(unlist(vals)) > 7){
        oxyLims <- quantile(vals$Oxygen,c(.01, .99))
        tempLims <- quantile(vals$Temperature,c(.01, .99))
        salLims <- quantile(vals$Salinity,c(.01, .99))
        oxygen3D <- (oxyLims[[1]] < oxygen) * (oxygen < oxyLims[[2]])
        temperature3D <- (tempLims[[1]] < temperature) * (temperature < tempLims[[2]])
        salinity3D <- (salLims[[1]] < salinity) * (salinity < salLims[[2]])
        
        envBuffer3D <- vect(trainingRegion3D[[which(grepl(pattern = paste0(spNames[[i]], ".shp"), 
                                                          x = trainingRegion3D, ignore.case = TRUE))]])
        
        # Creating a prediction stack
        depthsToProject <- sort(unique(occurrenceList[[spNames[[i]]]]$depth))
        layerNames <- names(temperature)
        index <- seq(from = match(min(depthsToProject), layerNames), 
                     to = match(max(depthsToProject), layerNames), by = 1)
        depthPred <- NULL
        for(k in index){
          depthPredLayer <- temperature3D[[k]] * salinity3D[[k]] * oxygen3D[[k]]
          crs(depthPredLayer) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
          depthPredLayer <- depthPredLayer
          depthPred[[k]] <- mask(depthPredLayer, mask = envBuffer3D)
          depthPred[[k]] <- crop(depthPred[[k]], y = envBuffer3D)
          depthPred[[k]] <- classify(depthPred[[k]], 
                                     matrix(c(NA, 0, 1, 1), 
                                            nrow = 2, byrow = TRUE), 
                                     include.lowest = T)
        }
        depthPred <- rast(depthPred)
        terra::writeRaster(depthPred, 
                           filename = paste0("Models/Envelope_3D/",
                                             spNames[[i]], ".tif"),
                    overwrite = T)
      }
    }
  }
}

# 2D GLM ----
spatialCVblocks <- function(occurrences){
  occurrences <- occurrences[order(occurrences$response, decreasing = T),]
  blockSplit <- ENMeval::get.block(occs = occurrences[occurrences$response == 1, 
                                                      c("decimalLongitude", "decimalLatitude")], 
                                   bg = occurrences[occurrences$response == 0, 
                                                    c("decimalLongitude", "decimalLatitude")])
  blockID <- c(blockSplit$occs.grp, blockSplit$bg.grp)
  return (blockID)
}

glm2D_nPres <- vector(mode="list", length(occsWdata2D))
glm2D_nAbs <- vector(mode="list", length(occsWdata2D))
glm2D_AIC <- vector(mode="list", length(occsWdata2D))
glm2D_NullDeviance <- vector(mode="list", length(occsWdata2D))
glm2d_ResidualDeviance <- vector(mode="list", length(occsWdata2D))
glm2D_AUC <- vector(mode="list", length(occsWdata2D))
glm2D_TSS <- vector(mode="list", length(occsWdata2D))
glm2d_Kappa <- vector(mode="list", length(occsWdata2D))
glm2d_threshold <- vector(mode="list", length(occsWdata2D))
glm2d_convergenceSetting <- as.list(rep(1000, times = length(occsWdata2D)))
glm2d_converged <- as.list(rep(F, times =length(occsWdata2D)))
for(i in 1:length(occsWdata2D)){
  if(all(!is.na(occsWdata2D[[spNames[[i]]]]), 
         any(grepl(pattern = paste0(spNames[[i]], ".shp"), 
                   x = trainingRegion2D, ignore.case = TRUE)))){
    vals <- occsWdata2D[[spNames[[i]]]]
    envBuffer2D <- vect(trainingRegion2D[[which(grepl(pattern = paste0(spNames[[i]], ".shp"), 
                                                      x = trainingRegion2D, ignore.case = TRUE))]])
    
    if(nrow(unique(vals[,c("decimalLongitude", 
                           "decimalLatitude")])) > 7){ # Will not model species for which there are fewer than 8 observations
      # Process presences
      response <- rep(1, times = nrow(vals))
      surfDat <- cbind(vals, response)
      suitableCentroid <- apply(surfDat[,c("Oxygen", "Salinity", "Temperature")], 
                                MARGIN = 2, FUN = mean)
      glm2D_nPres[[i]] <- nrow(surfDat)
      
      # Process absences
      if(model2DisBottom[[i]]){
        mPts <- mSampling2D(occs = surfDat, environmentalData2Dbott[[1]], mShp = envBuffer2D)
        mPts <- cbind(mPts, terra::extract(x = environmentalData2Dbott, y = mPts))
      } else{
        mPts <- mSampling2D(occs = surfDat, environmentalData2Dsurf[[1]], mShp = envBuffer2D)
        mPts <- cbind(mPts, terra::extract(x = environmentalData2Dsurf, y = mPts))
      }
      
      mPts <- mPts[complete.cases(mPts),]
      response <- rep(0, times = nrow(mPts))
      mPts <- cbind(mPts, response)
      rownames(mPts) <- NULL
      mPts$distance <- apply(mPts[,c("Oxygen", "Salinity", "Temperature")], MARGIN = 1, 
                             FUN = function(x) dist(rbind(suitableCentroid, x)))
      mPts$sampleWeight <- (mPts$distance - min(mPts$distance))/(max(mPts$distance)-min(mPts$distance))
      sampleForAbsence <- sample(x = rownames(mPts),
                                 size = round(nrow(mPts)*.75, 0), 
                                 prob = mPts$sampleWeight)
      mPts <- mPts[match(sampleForAbsence, rownames(mPts)),]
      mPts <- mPts[,colnames(surfDat)]
      glm2D_nAbs[[i]] <- nrow(mPts)
      surfDat <- rbind(surfDat, mPts)
      surfDat$spBlocks <- spatialCVblocks(surfDat)
      rownames(surfDat) <- NULL
      
      # Create model
      while(!glm2d_converged[[i]]){
        print("Not converged yet...")
        print(glm2d_convergenceSetting[[i]])
        surfModel <- glm2(formula = response ~ Temperature + I(Temperature^2) + Salinity + I(Salinity^2) + Oxygen + I(Oxygen^2), 
                         family = binomial(link = "logit"),  
                         data = surfDat, maxit = glm2d_convergenceSetting[[i]])
        glm2d_converged[[i]] <- surfModel$converged
        if(!glm2d_converged[[i]]){
          glm2d_convergenceSetting[[i]] <- glm2d_convergenceSetting[[i]] + 5000
        } 
      }
      glm2D_AIC[[i]] <- surfModel$aic
      glm2D_NullDeviance[[i]] <- surfModel$null.deviance
      glm2d_ResidualDeviance[[i]] <- surfModel$deviance
      
      # Evaluate and post-process the model
      par(mfrow=c(1,3)) 
      partial_response(surfModel, predictors = surfDat[,c("Temperature", "Salinity", "Oxygen")])
      preds_cv <- crossvalSDM(surfModel, traindat=surfDat, 
                              colname_species = 'response', kfold = surfDat[,"spBlocks"],
                              colname_pred = c("Temperature", "Salinity", "Oxygen"))
      surfEval <- evalSDM(surfDat$response, preds_cv, thresh.method = "MaxSens+Spec")
      glm2D_AUC[[i]] <- surfEval$AUC
      glm2D_TSS[[i]] <- surfEval$TSS
      glm2d_Kappa[[i]] <- surfEval$Kappa
      if(model2DisBottom[[i]]){
        surfPred <- mask(predict(environmentalData2Dbott, surfModel), envBuffer2D)
      } else {
        surfPred <- mask(predict(environmentalData2Dsurf, surfModel), envBuffer2D)
      }
      
      surfPred <- crop(surfPred, envBuffer2D)
      glm2d_threshold[[i]] <- quantile(terra::extract(surfPred, 
                                                      surfDat[surfDat$response == 1,
                                                              c("decimalLongitude", "decimalLatitude")]),
                                       .1, na.rm = T)[[1]]
      surfGLM <- surfPred > glm2d_threshold[[i]]
      surfGLM <- classify(surfGLM, matrix(c(NA, 0, FALSE, 0, TRUE, 1), 
                                          ncol = 2, byrow = TRUE))
      writeRaster(surfGLM, filename = paste0("Models/GLM_2D/", 
                                             spNames[[i]], ".asc"), overwrite = T)
    }else{
      glm2D_nPres[[i]] <- nrow(vals)
      glm2D_nAbs[[i]] <- NA
      glm2D_AUC[[i]] <- NA
      glm2D_AIC[[i]] <- NA
      glm2D_NullDeviance[[i]] <- NA
      glm2d_ResidualDeviance[[i]] <- NA
      glm2D_TSS[[i]] <- NA
      glm2d_Kappa[[i]] <- NA
      glm2d_threshold[[i]] <- NA
    }
  } else{
    glm2D_nPres[[i]] <- NA
    glm2D_nAbs[[i]] <- NA
    glm2D_AUC[[i]] <- NA
    glm2D_AIC[[i]] <- NA
    glm2D_NullDeviance[[i]] <- NA
    glm2d_ResidualDeviance[[i]] <- NA
    glm2D_TSS[[i]] <- NA
    glm2d_Kappa[[i]] <- NA
    glm2d_threshold[[i]] <- NA
  }
}
modelDiagnostics_2D <- cbind(spNames, glm2D_nPres, glm2D_nAbs, glm2d_convergenceSetting, glm2D_AUC,
                             glm2D_AIC, glm2D_NullDeviance, glm2d_ResidualDeviance,
                             glm2D_TSS, glm2d_Kappa, glm2d_threshold)
colnames(modelDiagnostics_2D) <- c("species", "nPres", "nAbs", "MaxIterations", "AUC",
                                   "AIC", "NullDeviance", "ResidualDeviance",
                                   "TSS", "Kappa", "threshold")
write.csv(x = modelDiagnostics_2D, "Models/GLM_2D/ModelDiagnosticStats.csv", row.names = F)

#3D GLM ----
glm3D_nPres <- vector(mode="list", length(occsWdata3D))
glm3D_nAbs <- vector(mode="list", length(occsWdata3D))
glm3D_AIC <- vector(mode="list", length(occsWdata3D))
glm3D_NullDeviance <- vector(mode="list", length(occsWdata3D))
glm3D_ResidualDeviance <- vector(mode="list", length(occsWdata3D))
glm3D_AUC <- vector(mode="list", length(occsWdata3D))
glm3D_TSS <- vector(mode="list", length(occsWdata3D))
glm3D_Kappa <- vector(mode="list", length(occsWdata3D))
glm3D_threshold <- vector(mode="list", length(occsWdata3D))
glm3d_convergenceSetting <- as.list(rep(1000, times = length(occsWdata3D)))
glm3d_converged <- as.list(rep(F, times =length(occsWdata3D)))
for(i in 1:length(occsWdata3D)){
  print(spNames[[i]])
  if(all(!is.na(occsWdata3D[[spNames[[i]]]]), 
         any(grepl(pattern = paste0(spNames[[i]], ".shp"), 
                   x = trainingRegion3D, ignore.case = TRUE)))){
    vals <- occsWdata3D[[spNames[[i]]]]
    envBuffer3D <- vect(trainingRegion3D[[which(grepl(pattern = paste0(spNames[[i]], ".shp"), 
                                                      x = trainingRegion3D, ignore.case = TRUE))]])
    if(nrow(unique(vals[,c("decimalLongitude", 
                           "decimalLatitude")])) > 7){ # Will not model species for which there are fewer than 8 observations
      # Process presences
      response <- rep(1, times = nrow(vals))
      depthDat <- cbind(vals, response)
      suitableCentroid <- apply(depthDat[,c("Oxygen", "Salinity", "Temperature")], 
                                MARGIN = 2, FUN = mean)
      glm3D_nPres[[i]] <- nrow(depthDat)
      
      # Process absences
      sampleDepthMin <- min(vals$depth) - diff(range(vals$depth))*.2
      sampleDepthMax <- max(vals$depth) + diff(range(vals$depth))*.2
      mPts <- mSampling3D(occs = vals, temperature, mShp = envBuffer3D, 
                          depthLimit = c(sampleDepthMin, sampleDepthMax))
      mPts$Oxygen <- xyzSample(mPts, oxygen)
      mPts$Temperature <- xyzSample(mPts, temperature)
      mPts$Salinity <- xyzSample(mPts, salinity)
      mPts <- mPts[complete.cases(mPts),]
      mPts$distance <- apply(mPts[,c("Oxygen", "Salinity", "Temperature")], MARGIN = 1, 
                             FUN = function(x) dist(rbind(suitableCentroid, x)))
      mPts$sampleWeight <- (mPts$distance - min(mPts$distance))/(max(mPts$distance)-min(mPts$distance))
      response <- rep(0, times = nrow(mPts))
      mPts <- cbind(mPts, response)
      depths <- unique(mPts$depth)
      sampleForAbsence <- NULL
      for(d in depths){
        prCount <- length(vals[vals$depth == d,])
        sampleForAbsence <- c(sampleForAbsence, 
                              sample(x = rownames(mPts[mPts$depth == d,]),
                                     size = min(round(nrow(mPts[mPts$depth == d,])*.25, 0), 
                                                sum(mPts$sampleWeight[mPts$depth == d]>0)),
                                     prob = mPts$sampleWeight[mPts$depth == d]))
      }
      mPts <- mPts[match(sampleForAbsence, rownames(mPts)),]
      glm3D_nAbs[[i]] <- nrow(mPts)
      depthDat <- rbind(depthDat, mPts[,colnames(depthDat)])
      depthDat$spBlocks <- spatialCVblocks(depthDat)
      
      # Generate model
      while(!glm3d_converged[[i]]){
        print(glm3d_convergenceSetting[[i]])
        depthModel <- glm2(formula = response ~ Temperature + I(Temperature^2) + Salinity + I(Salinity^2) + Oxygen + I(Oxygen^2), 
                          family = binomial(link = "logit"),  
                          data = depthDat, 
                          maxit = glm3d_convergenceSetting[[i]])
        glm3d_converged[[i]] <- surfModel$converged
        if(!glm3d_converged[[i]]){
          print("Not converged yet...")
          glm3d_convergenceSetting[[i]] <- glm3d_convergenceSetting[[i]] + 5000
        } 
      }
      glm3D_AIC[[i]] <- depthModel$aic
      glm3D_NullDeviance[[i]] <- depthModel$null.deviance
      glm3D_ResidualDeviance[[i]] <- depthModel$deviance
      par(mfrow=c(1,3)) 
      partial_response(depthModel, predictors = depthDat[,c("Temperature", "Salinity", "Oxygen")])
      preds_cv <- crossvalSDM(depthModel, kfold=depthDat[,"spBlocks"], 
                              traindat=depthDat, 
                              colname_species = 'response', 
                              colname_pred = c("Temperature", "Salinity", "Oxygen"))
      depthEval <- evalSDM(depthDat$response, preds_cv, thresh.method = "MaxSens+Spec")
      glm3D_AUC[[i]] <- depthEval$AUC
      glm3D_TSS[[i]] <- depthEval$TSS
      glm3D_Kappa[[i]] <- depthEval$Kappa
      
      # Project and post-process model
      depthsToProject <- sort(unique(occurrenceList[[i]]$depth))
      layerNames <- names(temperature)
      index <- seq(from = match(min(depthsToProject), layerNames), 
                   to = match(max(depthsToProject), layerNames), by = 1)
      depthPred <- NULL
      for(j in index){
        depthPreds <- c(temperature[[j]], salinity[[j]], oxygen[[j]])
        crs(depthPreds) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
        names(depthPreds) <- c("Temperature", "Salinity", "Oxygen")
        depthPred[[j]] <- mask(predict(depthPreds, depthModel), envBuffer3D)
        depthPred[[j]] <- crop(depthPred[[j]], envBuffer3D)
        names(depthPred[[j]]) <- layerNames[[j]]
      }
      depthPred <- depthPred[!unlist(lapply(depthPred, 
                                            FUN = function(x) is.null(x)))]
      depthPred <- rast(depthPred)
      glm3D_threshold[[i]] <- quantile(xyzSample(depthDat[depthDat$response == 1,c("depth", "decimalLongitude", "decimalLatitude")], 
                                                 depthPred), 
                                       .1, na.rm = T)[[1]] # MS90
      
      #Plot layers
      depthPredThresh <- depthPred > glm3D_threshold[[i]]
      depthPredRCL <- classify(depthPredThresh, matrix(c(NA, 0),
                                                       ncol = 2, byrow = TRUE))
      depthPredRCL <- depthPredRCL[[as.logical(terra::minmax(depthPredRCL)[2,])]]
      terra::writeRaster(depthPredRCL, filename = paste0("Models/GLM_3D/", spNames[[i]], ".tif"),
               overwrite = T)
    }else{
      glm3D_nPres[[i]] <- nrow(vals)
      glm3D_nAbs[[i]] <- NA
      glm3D_AUC[[i]] <- NA
      glm3D_AIC[[i]] <- NA
      glm3D_NullDeviance[[i]] <- NA
      glm3D_ResidualDeviance[[i]] <- NA
      glm3D_TSS[[i]] <- NA
      glm3D_Kappa[[i]] <- NA
      glm3D_threshold[[i]] <- NA
    }
  } else{
    glm3D_nPres[[i]] <- NA
    glm3D_nAbs[[i]] <- NA
    glm3D_AUC[[i]] <- NA
    glm3D_AIC[[i]] <- NA
    glm3D_NullDeviance[[i]] <- NA
    glm3D_ResidualDeviance[[i]] <- NA
    glm3D_TSS[[i]] <- NA
    glm3D_Kappa[[i]] <- NA
    glm3D_threshold[[i]] <- NA
  }
}
modelDiagnostics_3D <- cbind(spNames, glm3D_nPres, glm3D_nAbs, 
                             glm3d_convergenceSetting,
                             glm3D_AUC, glm3D_AIC, glm3D_NullDeviance,
                             glm3D_ResidualDeviance, glm3D_TSS, glm3D_Kappa,
                             glm3D_threshold)
colnames(modelDiagnostics_3D) <- c("species", "nPres", "nAbs",
                                   "MaxIterations",
                                   "AUC", "AIC", "NullDeviance",
                                   "ResidualDeviance", "TSS", "Kappa",
                                   "threshold")
write.csv(x = modelDiagnostics_3D, "Models/GLM_3D/ModelDiagnosticStats.csv", row.names = F)

# Visualizing individual species ----
env2D <- list.files(path = "Models/Envelope_2D/", pattern = "\\.asc", full.names = T)
env2D <- env2D[!grepl(env2D, pattern = ".aux.xml")]
env2D <- lapply(env2D, FUN = function(x){
  temporary <- rast(x)
  names(temporary) <- gsub(gsub(x, pattern = "Models/Envelope_2D//", replacement = ""), 
                           pattern = "\\.asc", replacement = "")
  return(temporary)
  })
env3D <- list.files(path = "Models/Envelope_3D/", pattern = ".tif", full.names = T)
env3D <- lapply(env3D, FUN = function(x){
  tryCatch(expr = rast(x), 
           error = function(e){
             message(e)
           },
           warning = function(w){
             message(w)
           })
})
glm2D <- list.files(path = "Models/GLM_2D/", pattern = "\\.asc", full.names = T)
glm2D <- glm2D[!grepl(glm2D, pattern = ".aux.xml")]
glm2D <- lapply(glm2D, FUN = function(x){
  temporary <- rast(x)
  names(temporary) <- gsub(gsub(x, pattern = "Models/GLM_2D//", replacement = ""), 
                           pattern = "\\.asc", replacement = "")
  return(temporary)
})
glm3D <- list.files(path = "Models/GLM_3D/", pattern = ".tif", full.names = T)
glm3D <- lapply(glm3D, FUN = function(x){
  tryCatch(expr = rast(x), 
           error = function(e){
             message(e)
           },
           warning = function(w){
             message(w)
           })
})
modelDiagnostics_2D <- read.csv("Models/GLM_2D/ModelDiagnosticStats.csv")
colnames(modelDiagnostics_2D)[c(4,7,8,11)] <- c("Max Iter.", "Null Dev.", "Res. Dev.", "Thresh.")
modelDiagnostics_3D <- read.csv("Models/GLM_3D/ModelDiagnosticStats.csv")
colnames(modelDiagnostics_3D)[c(4,7,8,11)] <- c("Max Iter.", "Null Dev.", "Res. Dev.", "Thresh.")

world <- vect(ne_countries(scale = "medium", returnclass = "sf")[1])
world <- aggregate(world)
world <- simplifyGeom(world)

for (i in 1:length(occurrenceList)){
  print(paste0(spNames[[i]], " is index ", i))
  tempTable <- as.data.frame(rbind(modelDiagnostics_2D[modelDiagnostics_2D$species == spNames[i],-1], 
                             modelDiagnostics_3D[modelDiagnostics_3D$species == spNames[i],-1]))
  if (sum(apply(tempTable[,4:10], MARGIN = 2, 
                FUN=function(x) is.na(x))) < 14){
    tempTable[, 4:10] <- round(sapply(tempTable[,4:10], 
                                      as.numeric), 3)
    label2D <- ifelse(test = model2DisBottom[[i]],
                      yes = "2D (bott.)", 
                      no = "2D (surf.)")
    stable.p <- ggtexttable(tempTable[,-6:-7], 
                            rows = c(label2D, "3D "), 
                            theme = ttheme("lBlackWhite"))
    
    indexEnv2d <- match(spNames[[i]], lapply(env2D, FUN = function(x) names(x)))
    indexEnv3d <- which(grepl(spNames[[i]], 
                              lapply(env3D, FUN = function(x) sources(x))))
    indexGLM2d <- match(spNames[[i]], lapply(glm2D, FUN = function(x) names(x)))
    indexGLM3d <- which(grepl(spNames[[i]], lapply(glm3D, 
                                                   FUN = function(x) sources(x))))
    
    env2Dtemp <- tryCatch({as.numeric(sum(env2D[[indexEnv2d]])>0)}, 
                                         error = function(e) NULL)
    env3Dtemp <- tryCatch({as.numeric(sum(env3D[[indexEnv3d]])>0)}, 
                          error = function(e) NULL)
    glm2Dtemp <- tryCatch({as.numeric(sum(glm2D[[indexGLM2d]])>0)},
                          error = function(e) NULL)
    glm3Dtemp <- tryCatch({as.numeric(sum(glm3D[[indexGLM3d]])>0)},
                          error = function(e) NULL)
    allMaps <- list(env2Dtemp, env3Dtemp, glm2Dtemp, glm3Dtemp)
    if(sum(unlist(lapply(allMaps, FUN = function(x) is.numeric(x))))>0){
      dummy <- Filter(function(x) class(x)=="SpatRaster",allMaps)[[1]]
      values(dummy) <- 0
      allMaps[is.na(allMaps)] <- dummy
      for (i in 1:4){
        if (is.numeric(allMaps[[i]])){
          allMaps[[i]] <- dummy
        } else if(is.null(allMaps[[i]])){
          allMaps[[i]] <- dummy
        }
      }
    }
    
    extentForLand <- ext(vect(lapply(Filter(function(x) class(x)=="SpatRaster",allMaps), FUN = function(y) vect(ext(y)))))
    
    land <- crop(world, extentForLand)
    map1 <- rasterComp(rast1 = allMaps[[1]], 
                       rast2 = allMaps[[2]],
                       title = paste0(spNames[[i]], " Envelope Models"), land = land,
                       rast1Name = "2D", rast2Name = "3D")
    
    map1 <- as_ggplot(as_grob(map1))
    
    map2 <- rasterComp(rast1 = allMaps[[3]], 
                       rast2 = allMaps[[4]], 
                       title = paste0(spNames[[i]], " GLM Models"), land = land,
                       rast1Name = "2D", rast2Name = "3D")
    
    map2 <- as_ggplot(as_grob(map2))
    
    space <- ggplot() + geom_blank()

    pdf(paste0("speciesReports/", spNames[[i]], ".pdf"))
    grid.arrange(map1, stable.p, space, map2, 
                 nrow = 4, heights = c(3,2, 0.3, 3), widths = 4, padding = unit(0, "line"))
    dev.off()
    rm(stable.p, map1, map2, env2Dtemp, env3Dtemp, glm2Dtemp, glm3Dtemp, allMaps, land)
  } else {
    pdf(paste0("speciesReports/", spNames[[i]], ".pdf"))
    plot.new()
    text(x=0.5, y=0.5, paste0("No models for \n", spNames[[i]], "."), font=1, cex=2)
    dev.off()
  }
}
