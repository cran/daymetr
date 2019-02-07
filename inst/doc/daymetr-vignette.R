## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

# load the library
library(daymetr)
library(ncdf4)
library(raster)
library(sp)
library(dplyr)
library(ggplot2)

# load tile outlines
tile_outlines


## ----eval = FALSE--------------------------------------------------------
#  df <- download_daymet(site = "Oak Ridge National Laboratories",
#                  lat = 36.0133,
#                  lon = -84.2625,
#                  start = 2000,
#                  end = 2010,
#                  internal = TRUE,
#                  simplify = TRUE) # return tidy data !!

## ----eval = FALSE--------------------------------------------------------
#  # code is not run
#  download_daymet_batch(file_location = 'my_sites.csv',
#                        start = 1980,
#                        end = 2010,
#                        internal = TRUE)

## ----eval = TRUE, echo = FALSE, warning=FALSE, message=FALSE, error=FALSE----
# load demo data in package
df <- read_daymet(system.file(package = "daymetr","extdata/demo_data.csv"),
                 simplify = TRUE)

## ----eval = TRUE---------------------------------------------------------
str(df)

## ----fig.width = 7, fig.height=3-----------------------------------------
# simple graph of Daymet data
df %>%
  mutate(date = as.Date(paste(year, yday, sep = "-"), "%Y-%j")) %>%
  filter(measurement == "tmax..deg.c.") %>%
  ggplot() +
  geom_line(aes(x = date, y = value)) +
  facet_wrap(~ measurement, ncol = 2)

## ----eval = FALSE--------------------------------------------------------
#  # code not run
#  df %>%
#    mutate(date = as.Date(paste(year, yday, sep = "-"), "%Y-%j")) %>%
#    group_by(date) %>%
#    filter(measurement == "tmax..deg.c." | measurement == "tmin..deg.c.") %>%
#    summarize(mean_temp = mean(value))

## ---- fig.width = 7, fig.height = 7--------------------------------------
# plot the tile outlines
# roughly painting a picture of North America
plot(tile_outlines)

## ----eval = FALSE--------------------------------------------------------
#  # code not run
#  
#  # Download tiled data for multiple years (1980 - 2012)
#  # based upon a geographic location.
#  download_daymet_tiles(location = c(36.0133,-84.2625),
#                        tiles = NULL,
#                        start = 1980,
#                        end = 2012,
#                        param = "ALL")
#  
#  # Download tiled data for multiple years (1980 - 2012)
#  # based upon a tile number (11207) restricted to the
#  # minimum temperature data only.
#  download_daymet_tiles(tiles = 11207,
#                        start = 1980,
#                        end = 2012,
#                        param = "tmin")
#  

## ---- eval = FALSE-------------------------------------------------------
#  # download monthly
#  download_daymet_ncss(location = c(34, -82, 33.75, -81.75),
#                       start = 1980,
#                       end = 1980,
#                       frequency = "monthly",
#                       param = c("tmin","tmax"),
#                       path = tempdir(),
#                       silent = TRUE)

## ---- fig.width = 7, fig.height = 7, warning=FALSE, message=FALSE--------
# read in the demo data from the package for speed
r <- raster::stack(system.file(package = "daymetr","extdata/tmin_monavg_1980_ncss.nc"))

# to set the correct projection use
raster::projection(r) <- "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +a=6378137 +b=6356752.314706705 +units=m +no_defs"

# reproject to lat lon
r <- raster::projectRaster(r, crs = "+init=epsg:4326")

# plot the monthly mean minimum temperature for 1980
plot(r)

## ---- , fig.width = 7, fig.height = 7, warning=FALSE, message=FALSE------
# plot the monthly mean minimum temperature for 1980
r_tmean <- daymet_grid_tmean(path = system.file(package = "daymetr","extdata"),
                            product = "monavg",
                            year = 1980,
                            internal = TRUE)

# plot the mean temperature raster
plot(r_tmean)

