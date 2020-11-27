#' Global country shapefile
#'
#' A vector shapefile with all country boundaries
#' @format: A .shp shapefile with country outline vectors
#' @source: Somewhere from the internet
"basemap/Country-Land.shp"

#' Daily direct normal insolation (DNI) per latitude
#'
#' An excel spreadsheet containing daily estimates of direct normal insolation
#' by latitude, for the time period January 1, 1985 to December 31, 2005.
#'
#' @format A .xlsx file with 7306 rows and 255 columns
#' @source the NASA Langley Research Center Atmospheric Science Data Center
#'   Surface meteorological and Solar Energy web portal supported by the NASA
#'   LaRC POWER Project
"data_insolation_per_latitude_1985-2005.xlsx"

#' Calculated size of winter hole
#'
#' A dataset containing pre-calculated winter hole coefficients for all
#' latitudes on the globe, based on DNI data from NASA. Variables are as
#' follows:
#'
#' @format A dataframe with 43,560 rows and 4 columns:
#' \describe{
#'  \item{lon}{Longitude}
#'  \item{lat}{Latitude}
#'  \item{w}{calculated size of winter hole, as described in Steininger,
#'     Grossmann, Prol, Grossmann and Williges}
#'  \item{a}{calculated vale for a, also described in (ibid)}
#'  }
"globalW"
