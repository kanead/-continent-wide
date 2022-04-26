####' Pan-African Vulture Home Range Analysis
####' Eswatini tracks - Ara Monadjem, Adam Kane, Andre Botha
####' 2 immature African white-backed vultures
####' Create regularised data on whole track and by month
####' Calculate KDE and MCP home ranges for whole track and by month

####' housekeeping ----

#' remove all
rm(list = ls())
graphics.off()

#' Load the required packages
library(readr)
library(tidyverse)
library(amt)
library(SDLfilter)
library(ggmap)
library(sf)
library(raster)
library(adehabitatHR)
library(adehabitatLT)

####' load the data ----
data_path <- "data_eswatini"   # path to the data

files <- dir(data_path, pattern = "*.csv") # get file names
length(files)

mydata <- files %>%
  # read in all the files, appending the path before the filename
  map( ~ read_csv(file.path(data_path, .))) %>%
  reduce(rbind)

mydata$long <- as.numeric(mydata$long)
mydata$lat <- as.numeric(mydata$lat)

#' Eswatini tracks are in UTC
head(mydata)
mydata$time <-
  as.POSIXct(mydata$time, format = "%d/%m/%Y %H:%M", tz = "UTC")
head(mydata)

####' clean the data ----
#' Check for duplicated observations (ones with same lat, long, timestamp,
#'  and individual identifier).
ind2 <- mydata %>% dplyr::select(time, long, lat, id) %>%
  duplicated
sum(ind2)
#' remove them
mydata$dups <- ind2
mydata <- filter(mydata, dups == "FALSE")
mydata

#' remove duplicated timestamps 
ind3 <- mydata %>% dplyr::select(time, id) %>%
  duplicated
sum(ind3)

mydata$dups <- ind3
mydata <- filter(mydata, dups == "FALSE")
mydata


#' filter extreme data based on a speed threshold
#' based on vmax which is km/hr
#' time needs to be labelled DateTime for these functions to work
names(mydata)[names(mydata) == 'time'] <- 'DateTime'
#' needs a column with number of satellite fixes which I make up here
mydata$qi <- 6
SDLfilterData <-
  ddfilter_speed(data.frame(mydata), vmax = 100, method = 1)
length(SDLfilterData$DateTime)
head(SDLfilterData)

#' rename everything as before
mydata <- SDLfilterData
names(mydata)[names(mydata) == 'DateTime'] <- 'time'

#' select only the columns we need
mydata <- dplyr::select(mydata, time, lat, long, id, species, study)
head(mydata)

# check the minimum time and the maximum time
min_time <- mydata %>% group_by(id) %>% slice(which.min(time))
data.frame(min_time)

max_time <- mydata %>% group_by(id) %>% slice(which.max(time))
data.frame(max_time)

#' determine the length of time each bird was tracked for
duration <-
  difftime(max_time$time, min_time$time, units = "days")
duration

####' use amt package functions ----
#' convert into a track using amt
#' We can also use lat, long, which will allow us to determine
#' time of day
trk <- mk_track(
  mydata,
  .x = long,
  .y = lat,
  .t = time,
  id = id,
  crs = CRS("+init=epsg:4326")
)
trk

#' Now it is easy to calculate day/night with either movement track
trk <- trk %>% time_of_day()

#' can remove the night time points as follows:
day_trk <- filter(trk, tod_ == "day") %>% arrange(id, t_)

#' summarise the sampling rate
data_summary <-
  day_trk %>% nest(-id) %>% mutate(sr = map(data, summarize_sampling_rate)) %>%
  amt::select(id, sr) %>% unnest %>% arrange(id)
data_summary

####' data regularization ----

#' measure the time difference between points for each bird ID using dplyr
#' - Group your data by ID
#' - Compute time diffs between each timestamp in your group (the 1st time diff is NA)
#' - Create a new ID that counts no. of prior time gaps that are large
#' - Split the ID into newID by using an underscore separator at large gaps

length(levels(as.factor(trk$id)))
#' specify the time difference to break the track at
time_difference <- 30 #' the units relate to those set in difftime
#' need to add the arrange function here otherwise the order gets messed up
trk2 <- day_trk %>%
  group_by(id) %>%
  mutate(timeDiff = c(NA, difftime(tail(t_,-1), head(t_,-1), units = "mins"))) %>%
  mutate(newID = paste(id, cumsum(
    !is.na(timeDiff) &
      timeDiff > time_difference
  ), sep = "_")) %>%
  arrange(id, t_) %>%
  ungroup()
head(trk2)
tail(trk2)

#' check the number of newIDs
length(levels(as.factor(trk2$newID)))

#' create a trajectory object using adehabitatLT
trk_ltraj <-
  as.ltraj(xy = trk2[, c("x_", "y_")],
           date = trk2$t_,
           id = trk2$newID)
head(trk_ltraj)

#' rediscretization of the trajectory
#'  time step we want for the rediscretization, in seconds
tstep <- 3600 # 3600 secs = 1 hour
newtr <- redisltraj(trk_ltraj, u = tstep, type = "time")
head(newtr[1])
head(newtr[2])
class(newtr)

#' convert to class data frame
trk3 <- ld(newtr)
head(trk3)
class(trk3$date)

#' group the IDs that were split if they had big gaps back together into their original ID structure
#' this involves accessing the name of the new ID that occurs before the underscore
trk3 <- separate(trk3,
                 col = id,
                 sep = "_",
                 into = c("ID", "NA"))
head(trk3)
levels(as.factor(trk3$ID))
length(levels(as.factor(trk3$ID)))

#' remove the resultant NA column that occurs after the split
trk3 <- dplyr::select(trk3, x, y, date, ID) %>% arrange(ID, date)
head(trk3)
tail(trk3)
class(trk3$date)

#' transform back into a lat long track object
trk_lat_long <-
  mk_track(
    trk3,
    .x = x,
    .y = y,
    .t = date,
    id = ID,
    crs = CRS("+init=epsg:4326")
  )
trk_lat_long <- trk_lat_long %>% arrange(id, t_)
trk_lat_long

#' can export the regularised track with lat long coords
write.table(
  trk_lat_long,
  "summary/reg_track_all_lat_long_supp.csv",
  sep = ",",
  row.names = F
)

####' home range analyis ----
#' need data in Albers equal area

#' transform back into a track
trk_albers <-
  mk_track(
    trk3,
    .x = x,
    .y = y,
    .t = date,
    id = ID,
    crs = CRS("+init=epsg:4326")
  )  %>%
  transform_coords(
    sp::CRS(
      #' we can transform the CRS of the data to an equal area projection
      #' https://epsg.io/102022
      "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    )
  )

#' get the areas using amt
trk_albers <- arrange(trk_albers, id, t_)

#' can export the regularised track with lat long coords
write.table(
  trk_albers,
  "summary/reg_track_all_albers_supp.csv",
  sep = ",",
  row.names = F
)

#' KDE
kde <- trk_albers %>% nest(-id) %>%
  mutate(kdearea = map(data, ~ hr_kde(., levels = c(0.95)) %>% hr_area)) %>%
  dplyr::select(id, kdearea) %>% unnest()

kde$area <-  kde$area / 1000000
kde_95 <- kde %>% filter(level == 0.95) %>% arrange(id)
kde_95

#' here for mcp
mcps <- trk_albers %>% nest(-id) %>%
  mutate(mcparea = map(data, ~ hr_mcp(., levels = c(0.95)) %>% hr_area)) %>%
  dplyr::select(id, mcparea) %>% unnest()

mcps$area <- mcps$area / 1000000
mcp_95 <- mcps %>% filter(level == 0.95) %>% arrange(id)
mcp_95$area <- as.double(mcp_95$area)
mcp_95


#####' export summary stats ----
#' combine the summary stats
data_summary$kde_95_reg <- kde_95$area
data_summary$mcps_95_reg <- mcp_95$area
data_summary$species = "wb"
data_summary$study = "eswatini"
data_summary

#' #' can export this data summary
write.table(
  data_summary,
  "summary/supp_data_summary_supp.csv",
  sep = ",",
  row.names = F
)

####' create monthly groups ----
#' we need to extract monthly home ranges, some animals were tracked for over a year
#' so we must include a year-month-id grouping variable
#' first combine year and month
trk3$yr_month <- format(trk3$date, format = "%Y/%m")
trk3$yr_month <- as.factor(trk3$yr_month)

#' we also need a year-month-day variable to see what the coverage is like over the course of
#' a month for each individual
trk3$yr_month_day <- format(trk3$date, format = "%Y/%m/%d")
trk3$yr_month_day <- as.factor(trk3$yr_month_day)
head(trk3)

#' count the number of unique days when grouped by id and and month
short_months <- trk3 %>%
  group_by(ID, yr_month) %>%
  summarise(count = n_distinct(yr_month_day)) %>% filter(count < 28) %>% droplevels()
short_months #' months that had fewer than 28 days of tracking data reported
short_months$yr_month

#' we merge the two data frames and force all = TRUE so even the values that don't have a count
#' are included, this allows us to extract the tracks that have ~ a month
#' of coverage
test <- merge(short_months, trk3, all = TRUE)
length(test$ID)
length(trk3$ID)
summary(test$count)
head(test)

#' keep only the rows with the NAs which are the counts > 28 i.e. data with ~ a month
#' of coverage
day_trk_mod <- test %>% dplyr::filter(is.na(count))
head(day_trk_mod)
summary(day_trk_mod$count)

#' now create a unique identifier that has the ID, year and month
#' We will use these to build monthly home ranges
day_trk_mod$identifier <-
  paste(day_trk_mod$ID, day_trk_mod$yr_month, sep = "_")
day_trk_mod$identifier <- as.factor(day_trk_mod$identifier)
head(day_trk_mod)

trk_lat_long_month <-
  mk_track(
    day_trk_mod,
    .x = x,
    .y = y,
    .t = date,
    id = identifier,
    crs = CRS("+init=epsg:4326")
  )
trk_lat_long_month
trk_lat_long_month <- trk_lat_long_month %>% arrange(id, t_)

#' can export the regularised track by month with lat long coords
write.table(
  trk_lat_long_month,
  "summary/reg_track_month_lat_long_supp.csv",
  sep = ",",
  row.names = F
)

####' home range analyis for monthly groups ----
#' need data in Albers equal area

#' transform back into a track
trk_albers_month <-
  mk_track(
    day_trk_mod,
    .x = x,
    .y = y,
    .t = date,
    id = identifier,
    crs = CRS("+init=epsg:4326")
  )  %>%
  transform_coords(
    sp::CRS(
      #' we can transform the CRS of the data to an equal area projection
      #' https://epsg.io/102022
      "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    )
  )

#' get the areas using amt
trk_albers_month <- arrange(trk_albers_month, id, t_)
trk_albers_month

#' can export the regularised track with lat long coords
write.table(
  trk_albers_month,
  "summary/reg_track_month_albers_supp.csv",
  sep = ",",
  row.names = F
)

#' KDE
kde <- trk_albers_month %>% nest(-id) %>%
  mutate(kdearea = map(data, ~ hr_kde(., levels = c(0.95)) %>% hr_area)) %>%
  dplyr::select(id, kdearea) %>% unnest()

kde$area <-  kde$area / 1000000
kde_95 <- kde %>% filter(level == 0.95) %>% arrange(id)
kde_95


#' here for mcp
mcps <- trk_albers_month %>% nest(-id) %>%
  mutate(mcparea = map(data, ~ hr_mcp(., levels = c(0.95)) %>% hr_area)) %>%
  dplyr::select(id, mcparea) %>% unnest()

mcps$area <- mcps$area / 1000000
mcp_95 <- mcps %>% filter(level == 0.95) %>% arrange(id)
mcp_95$area <- as.double(mcp_95$area)
mcp_95

#####' export summary stats for monthly ranges ----
#' combine the summary stats
data_summary_month <- data.frame(cbind(kde_95, mcp_95$area))
data_summary_month = rename(data_summary_month,
                            kde_95_reg = area,
                            mcps_95_reg = mcp_95.area)
data_summary_month$species = "wb"
data_summary_month$study = "eswatini"
#' add the bird ID back on along with the monthly unique identifier
data_summary_month$bird <- sub("_.*", "", data_summary_month$id)
data_summary_month <- data_summary_month %>% dplyr::select(-level)
data_summary_month
#' export the monthly summary
write.table(
  data_summary_month,
  "summary/supp_data_summary_month_supp.csv",
  sep = ",",
  row.names = F
)
