# ----------------------------------------------------------------------------#
# File: compAdvPV.R                                                        ----
#                                                                           ---
# Author: Keith Williges                                                    ---
# Description: Package for numerical evaluation of electricity trade model  ---
#                                                                           ---
# Other files: Steps for calculating representative year / size of winter   ---
# hole etc. can be found in solarDataPrep.r. Python scripts                 ---
# "dataPrepScript.py" and "dataCombineScript.py" were used on the raw data  ---
# from NASA.                                                                ---
# ----------------------------------------------------------------------------#

# ----------------------------------------------------------------------------#
# R environment setup                                                    ------
# ----------------------------------------------------------------------------#

getBasemap <- function(){
  # For mapping, loads the background shapefile
  return(read_sf('./Data/basemap/Country-Land.shp'))
}

getPvMap <- function(){
  # Load the worldPVOUT raster
  worldPVraster <- raster::raster('./Data/WorldPVOUT/worldPV_resampledRasterFull')
  # mutliply the raster by 365 as we currently have data only on daily amounts
  worldPVraster <- worldPVraster * 365
  return(worldPVraster)
}

getLatMin <- function(){
  return(latMin)
}

setLatMin <- function(latMinIn){
  assign("latMin", latMinIn, envir = .GlobalEnv)
}

getLatMax <- function(){
  return(latMax)
}

setLatMax <- function(latMaxIn){
  assign("latMax", latMaxIn, envir = .GlobalEnv)
}

getVb <- function(){
  return(VbVal)
}
setVb <- function(VbSet){
  assign("VbVal", VbSet, envir = .GlobalEnv)
}

getFb <- function(){
  return(FbVal)
}
setFb <- function(FbSet){
  assign("FbVal", FbSet, envir = .GlobalEnv)
}

getCRFg <- function(){
  return(CRFgVal)
}
setCRFg <- function(CRFgSet){
  assign("CRFgVal", CRFgSet, envir = .GlobalEnv)
}

getCRFs <- function(){
  return(CRFsVal)
}
setCRFs <- function(CRFsSet){
  assign("CRFsVal", CRFsSet, envir = .GlobalEnv)
}

# calculates Cpv
CpvCalc <- function(CRFgVal, CRFsVal, latmax, latmin, interpolation,
                    worldPVraster, isAorNS, isEW){

  # Cpv calc for Autarky or N/S trade
  # Equation: (((ocf * CRFg) + (sg * CRFs)) / worldPVOUT) * 1000
  if(isAorNS == TRUE && isEW == FALSE){
    # pull the OCF and SG values from interpolated data (limit to latmin and
    # latmax as defined in function call)
    OCFval <- subset(interpolation, key == "OCF" & Lat <= latmax &
                       Lat >= latmin, select = c(Lat, value))
    SGval <- subset(interpolation, key == "SG" & Lat <= latmax &
                      Lat >= latmin, select = c(Lat, value))

    # calculate the numerator of the Cpv calc ((ocf * CRFg) + (sg*CRFs))
    Cpv_temp <- dplyr::mutate(OCFval, numerator = (OCFval$value * CRFgVal) +
                         (SGval$value * CRFsVal))

    # remove the unnecessary column from Cpv_temp
    Cpv_temp <- Cpv_temp[,-2]

    # convert the Cpv_temp vector into a matrix to use in calculations
    Cpv_raster <- replicate(360, Cpv_temp$numerator)

    # limit the worldPV raster to the correct raster::extent
    e <- raster::extent(-179, 180, latmin - 0.5, latmax + 0.5)
    worldPVraster_c <- raster::crop(worldPVraster, e)

    # calculate the resulting Cpv (output as a matrix)
    result <- (Cpv_raster / (raster::as.matrix(worldPVraster_c))) * 8760000.0

    # convert result matrix to raster with raster::extent and projection equal to
    # worldPVraster_c
    rasterResult <- raster::raster(result)
    raster::extent(rasterResult) <- e
    proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "
    raster::crs(rasterResult) <- proj

    # rename the result layer to Cpv
    names(rasterResult) <- "Cpv"
    return(rasterResult)
  }
  else if(isEW == TRUE && isAorNS == FALSE){
    # Here calculates Cpv for the E/W trade situation
    # eq: ((OCF - 1.05)/Epv) * CRFg + ((730 * SG - 1) / (730 * Epv)) * CRFs

    # pull the OCF and SG values from interpolated data (limit to latmin and
    # latmax as defined in function call)
    OCFval <- subset(interpolation, key == "OCF" & Lat <= latmax &
                       Lat >= latmin, select = c(Lat, value))
    SGval <- subset(interpolation, key == "SG" & Lat <= latmax &
                      Lat >= latmin, select = c(Lat, value))

    # calculate the first half of equation (dealing with CRFg)
    CRFg_numerator <- dplyr::mutate(OCFval, numerator = OCFval$value - 1.05)

    # remove the unnecessary column from CRFg_numerator
    CRFg_numerator <- CRFg_numerator[,-2]

    # convert CRFg_numerator to a matrix to use in calculations
    CRFg_numerator_raster <- replicate(360, CRFg_numerator$numerator)

    # define the spatial raster::extent to keep calculations consistent
    e <- raster::extent(-179, 180, latmin - 0.5, latmax + 0.5)
    worldPVraster_c <- raster::crop(worldPVraster, e)

    # calculate the first half of the equation above for Cpv
    leftSideResult <- CRFg_numerator_raster / raster::as.matrix(worldPVraster_c) *
      CRFgVal * 1000.0 * 8760.0

    # convert leftSideResult matrix to raster with extent and projection equal
    # to worldPVraster_c
    leftSideResult_raster <- raster::raster(leftSideResult)
    raster::extent(leftSideResult_raster) <- e
    proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "
    raster::crs(leftSideResult_raster) <- proj

    # Now calculate the right side of equation dealing with CRFs
    CRFs_numerator <- dplyr::mutate(SGval, numerator = 730.0 * SGval$value - 1)

    # remove the unnecessary column from CRFg_numerator
    CRFs_numerator <- CRFs_numerator[,-2]

    # convert CRFg_numerator to a matrix to use in calculations
    CRFs_numerator_raster <- replicate(360, CRFs_numerator$numerator)

    # define the spatial raster::extent to keep calculations consistent
    e <- raster::extent(-179, 180, latmin - 0.5, latmax + 0.5)
    worldPVraster_c <- raster::crop(worldPVraster, e)

    # Multiply worldPVraster_c by 730 as defined in eq. 15 for E/W trade
    CRFs_denominator <- worldPVraster_c * 730.0

    # calculate the first half of the equation above for Cpv
    rightSideResult <- (CRFs_numerator_raster / raster::as.matrix(CRFs_denominator)) *
      CRFsVal * 1000.0 * 8760.0

    # convert leftSideResult matrix to raster with extent and projection equal
    # to worldPVraster_c
    rightSideResult_raster <- raster::raster(rightSideResult)
    raster::extent(rightSideResult_raster) <- e
    proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "
    raster::crs(rightSideResult_raster) <- proj

    # Now finish calculation combining left and right sides, calculating
    # resulting Cpv (output as matrix)
    rasterResult <- leftSideResult_raster + rightSideResult_raster

    # rename the result layer to Cpv
    names(rasterResult) <- "Cpv"

    # returns result Cpv raster
    return(rasterResult)
  }
  else if(isEW == FALSE && isAorNS == FALSE) {
    # Here calculates Cpv for the GG trade situation
    # Equation: CRFg / E*pv

    CRF <- rep(CRFgVal, times = abs(latmin) + latmax + 1)
    Crf_matrix <- replicate(360, CRF)

    Cpv_raster <- raster::raster(Crf_matrix)
    e <- raster::extent(-179, 180, latmin - 0.5, latmax + 0.5)
    raster::extent(Cpv_raster) <- e
    proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "
    raster::crs(Cpv_raster) <- proj
    worldPVraster_c <- raster::crop(worldPVraster, e)
    raster::extent(worldPVraster_c) <- e

    Cpv_raster <- Cpv_raster / (median(worldPVraster_c, na.rm = TRUE)) * 1000.0 * 8760.0

    Cpv_raster_mod <- overlay(Cpv_raster, worldPVraster_c, fun = function(x, y) {
      x[is.na(y[])] <- NA
      return(x)
    })
    return(Cpv_raster_mod)
  }
}

getE <- function(gridW, basemap){
  # calculate global E for use in beta later

  # calculate e, based on formula above
  globalE <- dplyr::mutate(gridW, e = (w + 2.0) / (2.0 * (w + 1.0)))

  # create XYZ matrix to convert to raster map
  globalE <- subset(globalE, select = c(lon, lat, e))

  # create raster from XYZ above
  eMap <- raster::rasterFromXYZ(globalE)
  # assign a coordinate reference to the newly created raster
  raster::crs(eMap) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'

  return(eMap)
}


getInterpolation <- function(){
  # opens the saved interpolated SG and OCF data files; in future, new estimation
  # method can go here (e.g. new model)
  interpolation <- readRDS('./Data/newInterpolation.Rda')
  interpolation <- dplyr::mutate(interpolation, value = ifelse(value < 0, 0, value))
  return(interpolation)
}

betaCalc <- function(gridW, globalE, VbIn, FbIn, Cpv, latmin, latmax, isNSorGG){
  # calculates beta based on all the inputs, returns a raster with beta value
  # needed to calculate beta: function(w, a, e, sLCOE, k, Cb)

  # define the raster::extent of all rasters to be used in this function
  ext <- raster::extent(-179, 180, latmin - 0.5, latmax + 0.5)

  # define variables for betaCalc (in raster format), beginning with w
  w <- gridW %>% subset(lat <= latmax & lat >= latmin, select = -a) %>%
    raster::rasterFromXYZ()
  raster::extent(w) <- ext
  raster::crs(w) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "

  # define a raster
  a <- gridW %>% subset(lat <= latmax & lat >= latmin, select = -w) %>%
    raster::rasterFromXYZ()
  raster::extent(a) <- ext
  raster::crs(a) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "

  e <- globalE
  raster::extent(e) <- ext
  raster::crs(e) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "

  raster::extent(Cpv) <- ext
  raster::crs(Cpv) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "

  Vb <- rep(VbIn, times = abs(latmin) + latmax + 1)
  Vb_matrix <- replicate(360, Vb)
  Vb_raster <- raster::raster(Vb_matrix)
  raster::extent(Vb_raster) <- ext
  raster::crs(Vb_raster) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "
  names(Vb_raster) <- "Vb_raster"

  Fb <- rep(FbIn, times = abs(latmin) + latmax + 1)
  Fb_matrix <- replicate(360, Fb)
  Fb_raster <- raster::raster(Fb_matrix)
  raster::extent(Fb_raster) <- ext
  raster::crs(Fb_raster) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "
  names(Fb_raster) <- "Fb_raster"

  # calculated on the basis of equation 12

  if(isNSorGG == FALSE){
    beta <- (a * Vb_raster + a * Fb_raster * (1.0 - e) - Cpv * w * (1.0 - e)) /
      (a * Vb_raster)
    names(beta) <- "beta"

    betaStar <- overlay(w, a, e, Cpv, Vb_raster, Fb_raster, beta,
                           fun = betaFun, forcefun = TRUE)
  }
  else if(isNSorGG == TRUE){
    beta <- overlay(Fb_raster, Vb_raster, Cpv, fun = function(Fb, Vb, Cpv) {
      return(ifelse((Fb + Vb) > Cpv, 1.0, 0.0))
    })

    names(beta) <- "beta"

    betaStar <- overlay(w, a, e, Cpv, Vb_raster, Fb_raster, beta,
                           fun = betaFunNS, forcefun = TRUE)
  }
  names(betaStar) <- "betaStar"

  result_raster <- stack(w, a, e, Cpv, Vb_raster, Fb_raster, beta, betaStar)
  return(result_raster)
}

betaFun <- function(w, a, e, Cpv, Vb_raster, Fb_raster, beta){
  # actual function specified to calculate value of beta

  ifelse(is.na(Cpv) | is.na(beta), NA,
         ifelse(Cpv + a * (Fb_raster + (Vb_raster / 2.0)) > (Fb_raster + Vb_raster),
                0.0,
                ifelse(beta > 1, 1.0,
                       ifelse(beta < e, e, beta))))
}

betaFunNS <- function(w, a, e, Cpv, Vb_raster, Fb_raster, beta){
  betaStar <- beta
  names(betaStar) <- "betaStar"
  return(betaStar)
}

unitCostFun <- function(w, a, e, Cpv, Vb_raster, Fb_raster, beta, betaStar){
  ifelse(is.na(betaStar), NA,
         ifelse(betaStar <= e,
                (betaStar / e) * (Cpv + a * (Fb_raster + (Vb_raster / 2))) +
                  ((e - betaStar) / e) * (Fb_raster + Vb_raster),
                Cpv * (1 + ((betaStar - e) / (1 - e)) * w) +
                  a * (1 - ((betaStar - e) / (1 - e))) * Fb_raster +
                  a * ((1-((betaStar - e) / (1 - e)))^2) * ((Vb_raster) / 2)))
}

unitCostNSFun <- function(w, a, e, Cpv, Vb_raster, Fb_raster, beta, betaStar){

  ifelse(Cpv > (Vb_raster + Fb_raster), Vb_raster + Fb_raster, Cpv)
  #ifelse(is.na(betaStar), NA, ifelse(betaStar >= 1, Cpv, Vb_raster + Fb_raster))
}

unitCostCalc <- function(betaResult_raster, isNS, isGG){
  # 7. Unit cost calc/map ------------------------------------------------------#
  # need: beta, Cpv, Fb, Vb, e, a, w, from betaResult_rasterBrick
  #-----------------------------------------------------------------------------#
  if(isNS == FALSE && isGG == FALSE){
    unitCost_raster <- overlay(betaResult_raster, fun = unitCostFun,
                               forcefun = TRUE)
  }
  else if(isNS == TRUE || isGG == TRUE){
    unitCost_raster <- overlay(betaResult_raster, fun = unitCostNSFun,
                               forcefun = TRUE)
  }
}

tradeEW_calc <- function(){
  # get global presets for latMin and Max, w, k, Cb, CRFg, CRFs
  latmin <- getLatMin()
  latmax <- getLatMax()
  gridW <- getW(latmin, latmax)
  Vb <- getVb()
  Fb <- getFb()
  CRFg <- getCRFg()
  CRFs <- getCRFs()

  # calculate Cpv, e, beta, and unit costs in East/West trade situation

  cPVMap <- CpvCalc(CRFg, CRFs, latmax, latmin, getInterpolation(),
                    getPvMap(), isAorNS = FALSE, isEW = TRUE)

  # For B2: e calc and map
  eMap <- getE(gridW, getBasemap())

  # For B3 <- beta calculation and map
  betaMap <- betaCalc(gridW, globalE = eMap, FbIn = Fb, VbIn = Vb, cPVMap,
                      latmin, latmax, isNSorGG = FALSE)

  # For B4 <- unit cost calcs
  unitCostMap <- unitCostCalc(betaMap, isNS = FALSE, isGG = FALSE)
  resultRaster <- stack(betaMap, unitCostMap)
  names(resultRaster[[9]]) <- "unitCost_E-W_trade"
  return(resultRaster)
}

tradeNS_calc <- function(){
  # get global presets for latMin and Max, w, k, Cb, CRFg, CRFs
  latmin <- getLatMin()
  latmax <- getLatMax()
  gridW <- getW(latmin, latmax)
  Vb <- getVb()
  Fb <- getFb()
  CRFg <- getCRFg()
  CRFs <- getCRFs()

  # calculates Cpv, e, beta, and unit costs in North-South trade

  # in N/S trade, w and a become zero, so redefine gridW_global to contain only
  # 0 values, and proceed through calculation as normal
  gridW$w <- 0.0
  gridW$a <- 0.0

  cPVMap <- CpvCalc(CRFg, CRFs, latmax, latmin, getInterpolation(),
                    getPvMap(), isAorNS = TRUE, isEW = FALSE)

  # For B2: e calc and map
  eMap <- getE(gridW, getBasemap())

  # For B3 <- beta calculation and map
  betaMap <- betaCalc(gridW, globalE = eMap, FbIn = Fb, VbIn = Vb, cPVMap,
                      latmin, latmax, isNSorGG = TRUE)

  # For B4 <- unit cost calcs
  unitCostMap <- unitCostCalc(betaMap, isNS = TRUE, isGG = FALSE)
  resultRaster <- stack(betaMap, unitCostMap)
  names(resultRaster[[9]]) <- "unitCost_N-S_trade"
  return(resultRaster)
}

tradeGG_calc <- function(){
  # calculates Cpv, e, beta, and unit cost of electricity in global grid trade

  # get global presets for latMin and Max, w, k, Cb, CRFg, CRFs
  latmin <- getLatMin()
  latmax <- getLatMax()
  gridW <- getW(latmin, latmax)
  Vb <- getVb()
  Fb <- getFb()
  CRFg <- getCRFg()
  CRFs <- 0

  # in GG trade, w and a become zero, so redefine gridW_global to contain only
  # 0 values, and proceed through calculation as normal
  gridW$w <- 0
  gridW$a <- 0

  cPVMap <- CpvCalc(CRFg, CRFs, latmax, latmin, getInterpolation(),
                    getPvMap(), isAorNS = FALSE, isEW = FALSE)

  # For B2: e calc and map
  eMap <- getE(gridW, getBasemap())

  # For B3 <- beta calculation and map
  betaMap <- betaCalc(gridW, globalE = eMap, FbIn = Fb, VbIn = Vb, cPVMap,
                      latmin, latmax, isNSorGG = TRUE)

  # For B4 <- unit cost calcs
  unitCostMap <- unitCostCalc(betaMap, isNS = FALSE, isGG = TRUE)
  resultRaster <- stack(betaMap, unitCostMap)
  names(resultRaster[[9]]) <- "unitCost_G-G_trade"
  return(resultRaster)
}

autarky_calc <- function(){
  # get global presets for latMin and Max, w, k, Cb, CRFg, CRFs
  latmin <- getLatMin()
  latmax <- getLatMax()
  gridW <- getW(latmin, latmax)
  Vb <- getVb()
  Fb <- getFb()
  CRFg <- getCRFg()
  CRFs <- getCRFs()

  cPVMap <- CpvCalc(CRFg, CRFs, latmax, latmin, getInterpolation(),
                    getPvMap(), isAorNS = TRUE, isEW = FALSE)

  # For B2: e calc and map
  eMap <- getE(gridW, getBasemap())

  # For B3 <- beta calculation and map
  betaMap <- betaCalc(gridW, globalE = eMap, FbIn = Fb, VbIn = Vb, cPVMap,
                      latmin, latmax, isNSorGG = FALSE)

  # For B4 <- unit cost calcs
  unitCostMap <- unitCostCalc(betaMap, isNS = FALSE, isGG = FALSE)
  names(unitCostMap) <- "unitCost_autarky"
  return(stack(betaMap, unitCostMap))
}

getRepYear <- function(){
  # load representative year data (to calculate repYearByLat)

  ## Import cleaned data (post python scripts) - don't run this unless you want
  # to waste a lot of time and take up 20 gigabytes of space
  # insolDataLong <- read_delim(
  #"~/Documents/GitHub/PV-global-trade/Data/results/clean/solar_-70_-180.txt",
  #                          " ", escape_double = FALSE, trim_ws = TRUE)

  # 1.2 calculate representative year ------------------------------------------#
  # calculates a representative year by finding the daily average value for each
  # grid cell over the 20 year period - again, don't run unless you want to waste
  # time

  # repYear <- insolDataLong %>%
  #  group_by(MO, DY, lat, lon) %>%
  #  summarise(swv_avg =(mean(swv_dwn)))

  # save the newly created dataset
  # saveRDS(repYear, file = "./Data/repYearNew.Rda")

  # load representative year sunlight data
  repYear <- readRDS('./Data/repYearNew.Rda') %>% ungroup()

  # add a date variable to repYear
  repYear$date <- as.Date(paste(repYear$MO,"-",repYear$DY, sep = ""),
                          format = "%m-%d")

  # omit NAs, and add a day of year variable
  repYear <- na.omit(repYear)
  repYear$yday <- yday(repYear$date)
  return(repYear)
}

getRepYearByLat <- function(){
  # load repYearByLat data

  # find the average swv per latitude
  repYearByLat <- getRepYear() %>% ungroup() %>%
    group_by(date, lat) %>%
    summarise(swvByLat = (mean(swv_avg))) %>% ungroup()

  # remove missing values from final dataset before analysis
  repYearByLat <- na.omit(repYearByLat)

  # add day of year variable
  repYearByLat$yday <- yday(repYearByLat$date)
  return(repYearByLat)
}

getW <- function(lat.min, lat.max){
  # Loads global_W dataset for use in further calcs, subsets it to limit to .e.g.
  # -66:66 degrees latitude
  gridW_global <- readRDS('./Data/globalW.Rda')

  #subset if necessary
  gridW_global <- subset(gridW_global, lat >= lat.min & lat <= lat.max)
  gridW_global = gridW_global[, c(2,1,3,4)]
  return(gridW_global)
}
