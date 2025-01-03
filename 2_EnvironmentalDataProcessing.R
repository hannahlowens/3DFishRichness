library(terra)
library(voluModel)
library(reshape)
library(fields)

#setwd("YOUR_DIRECTORY")

# Creating ocean volume ----
bathymetry <- rast("EnvironmentalData/ETOPO_2022_v1_60s_N90W180_bed.tif")
names(bathymetry) <- "bathymetry"
bathymetry <- aggregate(bathymetry, fact = 60, fun = min)
values(bathymetry)[values(bathymetry) > 0] <- NA
values(bathymetry) <- abs(values(bathymetry))

oceanOnly <- vect("FAO Fishing Areas 2005/FAO_AREAS.shp")
oceanOnly <- aggregate(oceanOnly)
oceanOnly <- project(oceanOnly, bathymetry)
bathymetry <- mask(bathymetry, mask = oceanOnly, touches = FALSE)

depthTemplate <- vect("EnvironmentalData/woa18_decav_t00mn01_shape/woa18_decav_t00mn01.shp")
depthTemplate <- as.numeric(gsub("[d,M_mean]", 
                                 "", names(depthTemplate)))
depthTemplate[1] <- 0
evenDeeper <- seq(from = max(depthTemplate), 
                  to = minmax(bathymetry)[[2]], 
                  by = 500)
depthTemplate <- unique(c(depthTemplate, evenDeeper))

wholeOceanTemplate <- NULL
for(i in depthTemplate){
  print(i)
  temp <- as.numeric(bathymetry > i)
  temp <- subst(temp, from = 0, to = NA)
  wholeOceanTemplate <- c(wholeOceanTemplate, temp)
}
names(wholeOceanTemplate) <- depthTemplate
wholeOceanTemplate <- rast(wholeOceanTemplate)

writeRaster(wholeOceanTemplate, filename = "EnvironmentalData/ProcessedEnvtData/wholeOceanTemplate.tif", overwrite = TRUE)

rm(bathymetry)

wholeOceanTemplate <- rast("EnvironmentalData/ProcessedEnvtData/wholeOceanTemplate.tif")

# Temperature ----
temperature <- vect("EnvironmentalData/woa18_decav_t00mn01_shape/woa18_decav_t00mn01.shp")
values(temperature)[values(temperature) == -999.999] <- NA
template <- centerPointRasterTemplate(temperature)
temperature <- rasterize(x = temperature, y = template, 
                         fun = mean, field = names(temperature))
names(temperature) <- gsub("[d,M_mean]", "", names(temperature))
names(temperature)[[1]] <- "0"

wholeOceanTemplate <- project(wholeOceanTemplate, temperature)
wholeOceanTemplate <- crop(wholeOceanTemplate, temperature)

for (i in 1:(nlyr(wholeOceanTemplate)-nlyr(temperature))) {
  temperature <- c(temperature, temperature[[nlyr(temperature)]])
}
names(temperature) <- names(wholeOceanTemplate)

# Fill in holes above and below, where possible
temperature <- approximate(temperature, 
                           method = "constant", 
                           rule = 1, 
                           z = c(0, diff(as.numeric(names(wholeOceanTemplate)))))

# Interpolate
temp <- temperature
for (i in 1:nlyr(temperature)){
  # Removes very clear outliers that are likely errors
  outliers <- boxplot(temp[[i]], range = 10)$out
  if(length(outliers) > 0){
    temp[[i]] <- subst(temp[[i]], outliers, NA)
  }
  
  # Replaces NA cells with interpolated values 
  temp[[i]] <- interpolateRaster(temp[[i]], REML = TRUE, miles = FALSE)
  temp[[i]] <- mask(crop(x = temp[[i]], 
                         y = wholeOceanTemplate[[i]]), 
                    mask = wholeOceanTemplate[[i]])
  temp[[i]] <- cover(focal(temp[[i]], w = 3, fun = "mean"), temp[[i]], values = NA)
}

names(temp) <- names(wholeOceanTemplate)
writeRaster(temp, "EnvironmentalData/ProcessedEnvtData/temperature.tif", overwrite = TRUE)
temp <- rast("EnvironmentalData/ProcessedEnvtData/temperature.tif")

bottom <- temp[[nlyr(temp)]]
for ( i in nlyr(temp):1){
  bottom <- cover(bottom, temp[[i]])
}

writeRaster(bottom, "EnvironmentalData/ProcessedEnvtData/bottomTemperature.tif", overwrite = TRUE)

template <- rast("EnvironmentalData/ProcessedEnvtData/temperature.tif")
wholeOceanTemplate <- project(wholeOceanTemplate, template)

# Salinity ----
salinity <- vect("EnvironmentalData/woa18_decav_s00mn01_shape/woa18_decav_s00mn01.shp")
values(salinity)[values(salinity) == -999.999] <- NA
salinity <- rasterize(x = salinity, y = template, 
                         fun = mean, field = names(salinity))
names(salinity) <- gsub("[d,M_mean]", "", names(salinity))
names(salinity)[[1]] <- "0"

for (i in 1:(nlyr(wholeOceanTemplate)-nlyr(salinity))) {
  salinity <- c(salinity, salinity[[nlyr(salinity)]])
}
names(salinity) <- names(wholeOceanTemplate)

# Fill in holes above and below, where possible
salinity <- approximate(salinity, 
                           method = "constant", 
                           rule = 1, 
                           z = c(0, diff(as.numeric(names(wholeOceanTemplate)))))

# Interpolate
temp <- salinity
for (i in 1:nlyr(salinity)){
  # Removes very clear outliers that are likely errors
  outliers <- boxplot(temp[[i]], range = 10)$out
  if(length(outliers) > 0){
    temp[[i]] <- subst(temp[[i]], outliers, NA)
  }
  
  # Replaces NA cells with interpolated values 
  temp[[i]] <- interpolateRaster(temp[[i]], REML = TRUE, miles = FALSE)
  temp[[i]] <- mask(crop(x = temp[[i]], 
                         y = wholeOceanTemplate[[i]]), 
                    mask = wholeOceanTemplate[[i]])
  temp[[i]] <- cover(focal(temp[[i]], w = 3, fun = "mean"), temp[[i]], values = NA)
}

names(temp) <- names(wholeOceanTemplate)
writeRaster(temp, "EnvironmentalData/ProcessedEnvtData/salinity.tif", overwrite = TRUE)
temp <- rast("EnvironmentalData/ProcessedEnvtData/salinity.tif")

bottom <- temp[[nlyr(temp)]]
for ( i in nlyr(temp):1){
  bottom <- cover(bottom, temp[[i]])
}
writeRaster(bottom, "EnvironmentalData/ProcessedEnvtData/bottomsalinity.tif", overwrite = TRUE)

# Dissolved Oxygen ----
dissOxygen <- vect("EnvironmentalData/woa18_all_o00mn01_shape/woa18_all_o00mn01.shp")
values(dissOxygen)[values(dissOxygen) == -999.999] <- NA
dissOxygen <- rasterize(x = dissOxygen, y = template, 
                      fun = mean, field = names(dissOxygen))
names(dissOxygen) <- gsub("[d,M_mean]", "", names(dissOxygen))
names(dissOxygen)[[1]] <- "0"

for (i in 1:(nlyr(wholeOceanTemplate)-nlyr(dissOxygen))) {
  dissOxygen <- c(dissOxygen, dissOxygen[[nlyr(dissOxygen)]])
}
names(dissOxygen) <- names(wholeOceanTemplate)

# Fill in holes above and below, where possible
dissOxygen <- approximate(dissOxygen, 
                        method = "constant", 
                        rule = 1, 
                        z = c(0, diff(as.numeric(names(wholeOceanTemplate)))))

# Interpolate
temp <- dissOxygen
for (i in 1:nlyr(dissOxygen)){
  # Removes very clear outliers that are likely errors
  outliers <- boxplot(temp[[i]], range = 10)$out
  if(length(outliers) > 0){
    temp[[i]] <- subst(temp[[i]], outliers, NA)
  }
  
  # Replaces NA cells with interpolated values 
  temp[[i]] <- interpolateRaster(temp[[i]], REML = TRUE, miles = FALSE)
  temp[[i]] <- mask(crop(x = temp[[i]], 
                         y = wholeOceanTemplate[[i]]), 
                    mask = wholeOceanTemplate[[i]])
  temp[[i]] <- cover(focal(temp[[i]], w = 3, fun = "mean"), temp[[i]], values = NA)
}

names(temp) <- names(wholeOceanTemplate)
writeRaster(temp, "EnvironmentalData/ProcessedEnvtData/dissOxygen.tif", overwrite = TRUE)
temp <- rast("EnvironmentalData/ProcessedEnvtData/dissOxygen.tif")

bottom <- temp[[nlyr(temp)]]
for ( i in nlyr(temp):1){
  bottom <- cover(bottom, temp[[i]])
}
writeRaster(bottom, "EnvironmentalData/ProcessedEnvtData/bottomdissOxygen.tif", overwrite = TRUE)

# Apparent Oxygen Utilization----
AOU <- vect("EnvironmentalData/woa18_all_A00mn01_shape/woa18_all_A00mn01.shp")
values(AOU)[values(AOU) == -999.999] <- NA
AOU <- rasterize(x = AOU, y = template, 
                 fun = mean, field = names(AOU))
names(AOU) <- gsub("[d,M_mean]", "", names(AOU))
names(AOU)[[1]] <- "0"

for (i in 1:(nlyr(wholeOceanTemplate)-nlyr(AOU))) {
  AOU <- c(AOU, AOU[[nlyr(AOU)]])
}
names(AOU) <- names(wholeOceanTemplate)

# Fill in holes above and below, where possible
AOU <- approximate(AOU, 
                   method = "constant", 
                   rule = 1, 
                   z = c(0, diff(as.numeric(names(wholeOceanTemplate)))))

# Interpolate side to side
temp <- AOU
for (i in 1:nlyr(AOU)){
  # Removes very clear outliers that are likely errors
  outliers <- boxplot(temp[[i]], range = 10)$out
  if(length(outliers) > 0){
    temp[[i]] <- subst(temp[[i]], outliers, NA)
  }
  
  # Replaces NA cells with interpolated values 
  temp[[i]] <- interpolateRaster(temp[[i]], REML = TRUE, miles = FALSE)
  temp[[i]] <- mask(crop(x = temp[[i]], 
                         y = wholeOceanTemplate[[i]]), 
                    mask = wholeOceanTemplate[[i]])
  temp[[i]] <- cover(focal(temp[[i]], w = 3, fun = "mean"), temp[[i]], values = NA)
}

names(temp) <- names(wholeOceanTemplate)
writeRaster(temp, "EnvironmentalData/ProcessedEnvtData/AOU.tif", overwrite = TRUE)
temp <- rast("EnvironmentalData/ProcessedEnvtData/AOU.tif")

bottom <- temp[[nlyr(temp)]]
for ( i in nlyr(temp):1){
  bottom <- cover(bottom, temp[[i]])
}
writeRaster(bottom, "EnvironmentalData/ProcessedEnvtData/bottomAOU.tif", overwrite = TRUE)

# Silicate ----
silicate <- vect("EnvironmentalData/woa18_all_i00mn01_shape/woa18_all_i00mn01.shp")
values(silicate)[values(silicate) == -999.999] <- NA
silicate <- rasterize(x = silicate, y = template, 
                 fun = mean, field = names(silicate))
names(silicate) <- gsub("[d,M_mean]", "", names(silicate))
names(silicate)[[1]] <- "0"

for (i in 1:(nlyr(wholeOceanTemplate)-nlyr(silicate))) {
  silicate <- c(silicate, silicate[[nlyr(silicate)]])
}
names(silicate) <- names(wholeOceanTemplate)

# Fill in holes above and below, where possible
silicate <- approximate(silicate, 
                   method = "constant", 
                   rule = 1, 
                   z = c(0, diff(as.numeric(names(wholeOceanTemplate)))))

# Interpolate side to side
temp <- silicate
for (i in 1:nlyr(silicate)){
  # Removes very clear outliers that are likely errors
  outliers <- boxplot(temp[[i]], range = 10)$out
  if(length(outliers) > 0){
    temp[[i]] <- subst(temp[[i]], outliers, NA)
  }
  
  # Replaces NA cells with interpolated values 
  temp[[i]] <- interpolateRaster(temp[[i]], REML = TRUE, miles = FALSE)
  temp[[i]] <- mask(crop(x = temp[[i]], 
                         y = wholeOceanTemplate[[i]]), 
                    mask = wholeOceanTemplate[[i]])
  temp[[i]] <- cover(focal(temp[[i]], w = 3, fun = "mean"), temp[[i]], values = NA)
}

names(temp) <- names(wholeOceanTemplate)
writeRaster(temp, "EnvironmentalData/ProcessedEnvtData/silicate.tif", overwrite = TRUE)
temp <- rast("EnvironmentalData/ProcessedEnvtData/silicate.tif")

bottom <- temp[[nlyr(temp)]]
for ( i in nlyr(temp):1){
  bottom <- cover(bottom, temp[[i]])
}
writeRaster(bottom, "EnvironmentalData/ProcessedEnvtData/bottomsilicate.tif", overwrite = TRUE)

# Phosphate ----
phosphate <- vect("EnvironmentalData/woa18_all_p00mn01_shape/woa18_all_p00mn01.shp")
values(phosphate)[values(phosphate) == -999.999] <- NA
phosphate <- rasterize(x = phosphate, y = template, 
                      fun = mean, field = names(phosphate))
names(phosphate) <- gsub("[d,M_mean]", "", names(phosphate))
names(phosphate)[[1]] <- "0"

for (i in 1:(nlyr(wholeOceanTemplate)-nlyr(phosphate))) {
  phosphate <- c(phosphate, phosphate[[nlyr(phosphate)]])
}
names(phosphate) <- names(wholeOceanTemplate)

# Fill in holes above and below, where possible
phosphate <- approximate(phosphate, 
                        method = "constant", 
                        rule = 1, 
                        z = c(0, diff(as.numeric(names(wholeOceanTemplate)))))

# Interpolate side to side
temp <- phosphate
for (i in 1:nlyr(phosphate)){
  # Removes very clear outliers that are likely errors
  outliers <- boxplot(temp[[i]], range = 10)$out
  if(length(outliers) > 0){
    temp[[i]] <- subst(temp[[i]], outliers, NA)
  }
  
  # Replaces NA cells with interpolated values 
  temp[[i]] <- interpolateRaster(temp[[i]], REML = TRUE, miles = FALSE)
  temp[[i]] <- mask(crop(x = temp[[i]], 
                         y = wholeOceanTemplate[[i]]), 
                    mask = wholeOceanTemplate[[i]])
  temp[[i]] <- cover(focal(temp[[i]], w = 3, fun = "mean"), temp[[i]], values = NA)
}

names(temp) <- names(wholeOceanTemplate)
writeRaster(temp, "EnvironmentalData/ProcessedEnvtData/phosphate.tif", overwrite = TRUE)
temp <- rast("EnvironmentalData/ProcessedEnvtData/phosphate.tif")

bottom <- temp[[nlyr(temp)]]
for ( i in nlyr(temp):1){
  bottom <- cover(bottom, temp[[i]])
}
writeRaster(bottom, "EnvironmentalData/ProcessedEnvtData/bottomphosphate.tif", overwrite = TRUE)

# Nitrate ----
nitrate <- vect("EnvironmentalData/woa18_all_n00mn01_shape/woa18_all_n00mn01.shp")
values(nitrate)[values(nitrate) == -999.999] <- NA
nitrate <- rasterize(x = nitrate, y = template, 
                      fun = mean, field = names(nitrate))
names(nitrate) <- gsub("[d,M_mean]", "", names(nitrate))
names(nitrate)[[1]] <- "0"

for (i in 1:(nlyr(wholeOceanTemplate)-nlyr(nitrate))) {
  nitrate <- c(nitrate, nitrate[[nlyr(nitrate)]])
}
names(nitrate) <- names(wholeOceanTemplate)

# Fill in holes above and below, where possible
nitrate <- approximate(nitrate, 
                        method = "constant", 
                        rule = 1, 
                        z = c(0, diff(as.numeric(names(wholeOceanTemplate)))))

# Interpolate side to side
temp <- nitrate
for (i in 1:nlyr(nitrate)){
  # Removes very clear outliers that are likely errors
  outliers <- boxplot(temp[[i]], range = 10)$out
  if(length(outliers) > 0){
    temp[[i]] <- subst(temp[[i]], outliers, NA)
  }
  
  # Replaces NA cells with interpolated values 
  temp[[i]] <- interpolateRaster(temp[[i]], REML = TRUE, miles = FALSE)
  temp[[i]] <- mask(crop(x = temp[[i]], 
                         y = wholeOceanTemplate[[i]]), 
                    mask = wholeOceanTemplate[[i]])
  temp[[i]] <- cover(focal(temp[[i]], w = 3, fun = "mean"), temp[[i]], values = NA)
}

names(temp) <- names(wholeOceanTemplate)
writeRaster(temp, "EnvironmentalData/ProcessedEnvtData/nitrate.tif", overwrite = TRUE)
temp <- rast("EnvironmentalData/ProcessedEnvtData/nitrate.tif")

bottom <- temp[[nlyr(temp)]]
for ( i in nlyr(temp):1){
  bottom <- cover(bottom, temp[[i]])
}
writeRaster(bottom, "EnvironmentalData/ProcessedEnvtData/bottomnitrate.tif", overwrite = TRUE)

# Mixed Layer Depth ----
mld <- vect("EnvironmentalData/woa18_A5B7_M00mn01_shape/woa18_A5B7_M00mn01.shp")
values(mld)[values(mld) == -999.999] <- NA
mld <- rasterize(x = mld, y = template, 
                     fun = mean, field = names(mld))

# Interpolate side to side
temp <- mld
# Removes very clear outliers that are likely errors
outliers <- boxplot(temp, range = 10)$out
if(length(outliers) > 0){
  temp <- subst(temp, outliers, NA)
}

# Replaces NA cells with interpolated values 
temp <- interpolateRaster(temp, REML = TRUE, miles = FALSE)
temp <- mask(crop(x = temp, y = wholeOceanTemplate[[1]]), 
             mask = wholeOceanTemplate[[1]])
temp <- cover(focal(temp, w = 3, fun = "mean"), temp, values = NA)

t2 <- wholeOceanTemplate
for (i in 1:nlyr(wholeOceanTemplate)){
  t2[[i]] <- as.numeric(temp-as.numeric(names(t2[[i]])))
}

names(t2) <- names(wholeOceanTemplate)
writeRaster(t2, "EnvironmentalData/ProcessedEnvtData/mldDistance.tif", overwrite = TRUE)
mldDist <- rast("EnvironmentalData/ProcessedEnvtData/mldDistance.tif")

rm(temp, t2)

# Density ----
density <- vect("EnvironmentalData/woa18_decav_I00mn01_shape/woa18_decav_I00mn01.shp")
values(density)[values(density) == -999.999] <- NA
density <- rasterize(x = density, y = template, 
                      fun = mean, field = names(density))
names(density) <- gsub("[d,M_mean]", "", names(density))
names(density)[[1]] <- "0"

for (i in 1:(nlyr(wholeOceanTemplate)-nlyr(density))) {
  density <- c(density, density[[nlyr(density)]])
}
names(density) <- names(wholeOceanTemplate)

# Fill in holes above and below, where possible
density <- approximate(density,
                       method = "constant",
                       rule = 1,
                       z = c(0, diff(as.numeric(names(wholeOceanTemplate)))))

# Interpolate
temp <- density
for (i in 1:nlyr(density)){
  # Removes very clear outliers that are likely errors
  outliers <- boxplot(temp[[i]], range = 10)$out
  if(length(outliers) > 0){
    temp[[i]] <- subst(temp[[i]], outliers, NA)
  }
  
  # Replaces NA cells with interpolated values 
  temp[[i]] <- interpolateRaster(temp[[i]], REML = TRUE, miles = FALSE)
  temp[[i]] <- mask(crop(x = temp[[i]], 
                         y = wholeOceanTemplate[[i]]), 
                    mask = wholeOceanTemplate[[i]])
  temp[[i]] <- cover(focal(temp[[i]], w = 3, fun = "mean"), temp[[i]], values = NA)
}

names(temp) <- names(wholeOceanTemplate)
writeRaster(temp, "EnvironmentalData/ProcessedEnvtData/density.tif", overwrite = TRUE)

# Correlation analysis
# Load environmental data ----
oxygen <- scale(rast("EnvironmentalData/ProcessedEnvtData/dissOxygen.tif"))
AOU <- scale(rast("EnvironmentalData/ProcessedEnvtData/AOU.tif"))
salinity <- scale(rast("EnvironmentalData/ProcessedEnvtData/salinity.tif"))
temperature <- scale(rast("EnvironmentalData/ProcessedEnvtData/temperature.tif"))

bathymetry <- rast("EnvironmentalData/ETOPO_2022_v1_60s_N90W180_bed.tif")
names(bathymetry) <- "bathymetry"
bathymetry <- aggregate(bathymetry, fact = 60, fun = min)
values(bathymetry)[values(bathymetry) > 0] <- NA
values(bathymetry) <- abs(values(bathymetry))

oceanOnly <- vect("FAO Fishing Areas 2005/FAO_AREAS.shp")
oceanOnly <- aggregate(oceanOnly)
oceanOnly <- project(oceanOnly, temperature)

# Testing for correlation among variables
dummyOccs <- data.frame(0, 0, 0)
colnames(dummyOccs) <- c("longitude", "latitude", "depth")

grid <- mSampling3D(occs = dummyOccs, 
                    envBrick = temperature, 
                    mShp = oceanOnly)

set.seed(seed = 200)
sampledGrid <- grid[sample(1:nrow(grid),
                           size = 50000,
                           replace = FALSE),]

allData <- cbind(xyzSample(sampledGrid, envBrick = oxygen, verbose = F),
                 xyzSample(sampledGrid, envBrick = AOU, verbose = F),
                 xyzSample(sampledGrid, envBrick = temperature, verbose = F),
                 xyzSample(sampledGrid, envBrick = salinity, verbose = F))
colnames(allData) <- c("dissolvedOxygen", "AOU", 
                       "temperature", "salinity")
allData <- allData[complete.cases(allData),]

pearson <- round(cor(allData, method = "pearson"), 2) # Remove dissolved oxygen (Correlation coefficient cutoff 0.4)

rm(allData, AOU, dummyOccs, grid, pearson, sampledGrid, oceanOnly)