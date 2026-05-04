library(here)
source(file=here::here("functions/fun.R"))

f_load_libraries_colors()

token_NOAA <- "ZFyQJPeWOXcEbvwIDWBvVJmfyoVPoTll"

library(worldmet)
library(lubridate)

lat <- 48.8
lon <- 2.083333

# 1) Trouver les stations NOAA/ISD les plus proches
stations <- import_isd_stations(
  lat = lat,
  lon = lon,
  n_max = 10
)

stations %>%
  dplyr::select(
    usaf, wban, code, station, ctry,
    latitude, longitude, `elev(m)`, begin, end, dist
  ) %>%
  print(n = 10)

# 2) Prendre la station la plus proche
station_code <- stations$code[1]
station_name <- stations$name[1]

station_code
station_name

# 3) Télécharger les températures horaires pour une ou plusieurs années
meteo <- import_isd_hourly(
  code = station_code,
  year = 2024
)

# 4) Garder seulement la température
temp <- meteo %>%
  dplyr::select(date, station = code, site = station, latitude, longitude, air_temp) %>%
  filter(!is.na(air_temp))

head(temp)
summary(temp$air_temp)