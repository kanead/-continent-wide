####' Pan-African Vulture Home Range Analysis
####' Eswatini tracks - Ara Monadjem, Adam Kane, Andre Botha
####' 2 immature African white-backed vultures
####' Identification of breeding sites using recurse package
####' This is done on the raw data i.e. not the regularised tracks

####' housekeeping ----

#' remove all
rm(list = ls())
graphics.off()

#' Load the required packages
library(readr)
library(tidyverse)
library(amt)
library(SDLfilter)
library(recurse)

####' load the data ----
data_path <- "data_eswatini"   # path to the data

files <- dir(data_path, pattern = "*.csv") # get file names
length(files)

mydata <- files %>%
  # read in all the files, appending the path before the filename
  map(~ read_csv(file.path(data_path, .))) %>%
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

trk <- trk %>% arrange(id, t_)

#' Now it is easy to calculate day/night with either movement track
trk <- trk %>% time_of_day()

#' can remove the night time points as follows:
day_trk <- filter(trk, tod_ == "day") %>% arrange(id, t_)

####' create monthly groups ----
#' we need to extract monthly home ranges, some animals were tracked for over a year
#' so we must include a year-month-id grouping variable
#' first combine year and month
day_trk$yr_month <- format(day_trk$t_, format = "%Y/%m")
day_trk$yr_month <- as.factor(day_trk$yr_month)

#' we also need a year-month-day variable to see what the coverage is like over the course of
#' a month for each individual
day_trk$yr_month_day <- format(day_trk$t_, format = "%Y/%m/%d")
day_trk$yr_month_day <- as.factor(day_trk$yr_month_day)
head(day_trk)

#' count the number of unique days when grouped by id and and month
short_months <- day_trk %>%
  group_by(id, yr_month) %>%
  summarise(count = n_distinct(yr_month_day)) %>% filter(count < 28) %>% droplevels()
short_months
short_months$yr_month

#' we merge the two data frames and force all = TRUE so even the values that don't have a count
#' are included, this allows us to extract the tracks that have ~ a month
#' of coverage
test <- merge(short_months, day_trk, all = TRUE)
length(test$id)
length(day_trk$id)
summary(test$count)

#' keep only the rows with the NAs which are the counts > 28 i.e. data with ~ a month
#' of coverage
day_trk_mod <- test %>% dplyr::filter(is.na(count))
head(day_trk_mod)
summary(day_trk_mod$count)

#' now create a unique identifier that has the ID, year and month
#' We will use these to build home ranges
day_trk_mod$identifier <-
  paste(day_trk_mod$id, day_trk_mod$yr_month, sep = "_")
day_trk_mod$identifier <- as.factor(day_trk_mod$identifier)
head(day_trk_mod)

#' add the month and extract months where the bird could be breeding
day_trk_mod$month <- format(day_trk_mod$t_, format = "%m")
head(day_trk_mod)
day_trk_mod$month <- as.numeric(day_trk_mod$month)
head(day_trk_mod)

#' turn this back into a track object
day_trk <- mk_track(
  day_trk_mod,
  .x = x_,
  .y = y_,
  .t = t_,
  id = identifier,
  crs = CRS("+init=epsg:4326")
)
day_trk

##### albers equal area projection ----
#' save the track object by nesting according to id and
#' converting coords into albers equal area so units are in metres
trk1 <- day_trk %>%
  transform_coords(
    sp::CRS(
      #' we can transform the CRS of the data to an equal area projection
      #' https://epsg.io/102022
      "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    )
  )
trk1

####' Calculate revisits in 50m radius ----
memory.limit(size = 56000) # increase memory size if needs be

#' split the track into a list where each element is a unique ID-month-year identifier
#' The revisits are calculated on this basis
point_list <-
  split(trk1, trk1$id)
point_list <- discard(point_list, function(x)
  nrow(x) == 0)
point_list[1]

#' apply the recursion function to this list to give the time spent in an area in days
recursions <- lapply(point_list, function(x)
  getRecursions(data.frame(x), 50, timeunits = "days"))
summary(recursions)

#' find the max time for each of the ID-month-year identifiers
n <- c()
for (i in 1:length(recursions)) {
  n[i] <- max(recursions[[i]]$residenceTime)
}

#' have a look at the distribution, the units are in days
ggplot(data = data.frame(n), mapping = aes(x = n)) + geom_histogram()

#' combine the IDs with the max time data from the recursion analysis
recursion_time <- cbind(data.frame(names(point_list), n))
recursion_time

####' assemble the data for export ----
recursion_time <- recursion_time %>% rename(id = names.point_list.)
recursion_time <- recursion_time %>% rename(time = n)
#' pull out the ID which occurs before the first _
recursion_time$bird <-
  gsub("(.+?)(\\_.*)", "\\1", recursion_time$id)
recursion_time$study <- "eswatini"
recursion_time
write.csv(recursion_time, file = "summary/supp_data_recursion_supp.csv", row.names = FALSE)
