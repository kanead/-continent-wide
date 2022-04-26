####' Pan-African Vulture Home Range Analysis
####' Eswatini tracks - Ara Monadjem, Adam Kane, Andre Botha
####' 2 immature African white-backed vultures
####' Brownian Bridge on Regularised data for full duration and by month

####' housekeeping ----

#' remove all
rm(list = ls())
graphics.off()

#' load packages
library(tidyverse)
library(move)
library(sf)
library(amt)

####' load the data ----
trk1 <-
  read_csv("summary/reg_track_all_albers_supp.csv", col_names = TRUE)
trk1

track <- data.frame(trk1)
track <- arrange(track, id)
head(track)

#' create a move object
loc <-
  move(
    x = track$x_,
    y = track$y_,
    time = track$t_,
    proj = CRS(
      "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    ),
    data = track,
    animal = track$id
  )

#' Now create a dBBMM object
dbbmm <-
  brownian.bridge.dyn(
    object = loc,
    location.error = 20,
    window.size = 31,
    margin = 11,
    raster = 10000
  )

#' use a function to extract the areas of the brownian bridges
#' this removes the IDs which must be appended later
dbbmm_size <- function (x) {
  return(tryCatch(
    st_area(hr_isopleths(dbbmm[[x]], level = 0.95)) / 1e6,
    error = function(e)
      NULL
  ))
}
y <- c(1:length(dbbmm[1]))
areas <- lapply((y), dbbmm_size)
areas

#' add names by taking them from dbbmm object
names(dbbmm)
names(areas) <- names(dbbmm)
bb_areas <- data.frame(unlist(areas))
length(bb_areas$unlist.areas.)

#' try the function on its own for comparison
st_area(hr_isopleths(dbbmm$ID1, level = 0.95)) / 1e6 #' 13229.73
st_area(hr_isopleths(dbbmm$ID2, level = 0.95)) / 1e6 #' 51904.64

#' load in the cleaned version which is in Albers Equal Area
merged_Africa_tranform <-
  read_sf("shapefile//merged_Africa_protected_clean.shp")
st_crs(merged_Africa_tranform)

#' what's the area of overlap with the protected areas?
#' this is at 95%
#' test it on one here
intersection_1_95 <-
  st_intersection(hr_isopleths(dbbmm$ID1, level = 0.95),
                  merged_Africa_tranform$geometry)
sum(st_area(intersection_1_95)) / 1e6 #' 3031.541

#' function to calculate for overlap between BB and protected areas
dbbmm_overlap_size <- function (x) {
  return(tryCatch(
    sum(st_area(
      st_intersection(
        hr_isopleths(dbbmm[[x]], level = 0.95),
        merged_Africa_tranform$geometry
      )
    )) / 1e6 ,
    error = function(e)
      NULL
  ))
}

y <- c(1:length(dbbmm[1]))
#' run it over all of the IDs
pa_overlap <- lapply(y, dbbmm_overlap_size)
pa_overlap
names(pa_overlap) <- names(dbbmm)
pabb_areas <- data.frame(unlist(pa_overlap))
length(pabb_areas$unlist.pa_overlap.)

####' prep the data for export ----
export_data <- data.frame(cbind(bb_areas, pabb_areas))
#' rename columns
export_data <- export_data %>% rename(bb_area_95 = unlist.areas.)
export_data <-
  export_data %>% rename(bb_area_95_overlap = unlist.pa_overlap.)
export_data
export_data <- cbind(ID = rownames(export_data), export_data)
rownames(export_data) <- NULL
#' remove the X that gets added to the IDs
levels(as.factor(trk1$id))
export_data$ID <-
  gsub(pattern = "X",
       replacement = "",
       x = export_data$ID)
export_data

write.csv(export_data, file = "summary/supp_data_bb_reg_raster_supp.csv", row.names = F)

####' repeat on the monthly regularised data ----
####' load the data ----
trk1 <-
  read_csv("summary/reg_track_month_albers_supp.csv", col_names = TRUE)
trk1

track <- data.frame(trk1)
track <- arrange(track, id)
head(track)

#' create a move object
loc <-
  move(
    x = track$x_,
    y = track$y_,
    time = track$t_,
    proj = CRS(
      "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    ),
    data = track,
    animal = track$id
  )

#' Now create a dBBMM object
dbbmm <-
  brownian.bridge.dyn(
    object = loc,
    location.error = 20,
    window.size = 31,
    margin = 11,
    raster = 10000
  )

#' use a function to extract the areas of the brownian bridges
y <- c(1:length(dbbmm[1]))
areas <- lapply((y), dbbmm_size)
areas

#' add names by taking them from dbbmm object
names(dbbmm)
names(areas) <- names(dbbmm)
bb_areas <- data.frame(unlist(areas))
length(bb_areas$unlist.areas.)

#' try the function on its own for comparison
st_area(hr_isopleths(dbbmm$ID1_2015.06, level = 0.95)) / 1e6 #' 4758.755
st_area(hr_isopleths(dbbmm$ID1_2015.07, level = 0.95)) / 1e6 #' 5281.56

#' what's the area of overlap with the protected areas?
#' this is at 95%
#' test it on one here
intersection_1_95 <-
  st_intersection(hr_isopleths(dbbmm$ID1_2015.06, level = 0.95),
                  merged_Africa_tranform$geometry)
sum(st_area(intersection_1_95)) / 1e6 #' 661.5312

#' function to calculate for overlap between BB and protected areas
y <- c(1:length(dbbmm[1]))
#' run it over all of the IDs
pa_overlap <- lapply(y, dbbmm_overlap_size)
pa_overlap
names(pa_overlap) <- names(dbbmm)
pabb_areas <- data.frame(unlist(pa_overlap))
length(pabb_areas$unlist.pa_overlap.)

####' prep the monthly data for export ----
export_data <- data.frame(cbind(bb_areas, pabb_areas))
#' rename columns
export_data <- export_data %>% rename(bb_area_95 = unlist.areas.)
export_data <-
  export_data %>% rename(bb_area_95_overlap = unlist.pa_overlap.)
export_data
export_data <- cbind(ID = rownames(export_data), export_data)
rownames(export_data) <- NULL
#' replace the periods with slashes to keep the IDs consistent
levels(as.factor(trk1$id))
export_data$ID <-
  gsub(pattern = "\\.",
       replacement = "/",
       x = export_data$ID)
export_data

write.csv(export_data, file = "summary/supp_data_bb_reg_raster_month_supp.csv", row.names = F)
