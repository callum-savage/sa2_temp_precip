# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#
# sa2_temp_precip.R
# 
# Author: Callum Savage <callum.savage@anu.edu.au>
# Date: October 2021
#
# This is an R script to aggregate total weekly rainfall and maximum temperature
# by Statistical Area 2 (SA2) for the years 2000 - 2018.
# 
# The general approach is as follows:
# 
#   - Downlad the SA2 shapefile from the ABS
#   
#   - Download precipitation and maximum temperature data from NCI. It is 
#     available as a 0.05 x 0.05 degree (approx 5 km^2) raster file with files 
#     for each year
#   
#   - Identify which cells overlap with which SA2 areas and calculate the 
#     proportion that each SA2 area covers the cell
#     
#   - Total precipitation calculation:
#     - For each year of data join the SA2 IDs to their grid cell precipitation
#       values
#     - Multiply precipitation values by the cell area (dependent on latitude) 
#       to calculate a daily rainfall volume
#     - Calculate a weighted sum of precipitation volume by coverage fraction,
#       and divide by SA2 area to calculate the daily rainfall in mm
#     - Sum the daily rainfall into weeks
#     
#   - Maximum temperature calculation:
#     - For each year of data join the SA2 IDs to their grid cell maximum
#       temperature values
#     - Keep the highest maximum temperature value for each SA2 area
#     - Aggregate the daily maximum temperatures into weeks by again taking the
#       maximum
#       
#   - Combine the weekly precipitation and maximum temperature data sets into 
#     a single data frame and save the output to file.
# 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Setup -------------------------------------------------------------------

# Library imports

library(here)             # Relative file paths
library(glue)             # Interpreted String Literals
library(terra)            # Spatial raster data
library(sf)               # Spatial vector data
library(exactextractr)    # Fast Extraction from Raster Datasets using Polygons
library(dplyr)            # Data manipulation
library(tidyr)            # Data cleaning
library(purrr)            # Functional programming tools
library(lubridate)        # Improved date management
library(readr)            # File read/write
library(arrow)            # Feather file format
library(assertr)          # Data quality checks
library(rsync)            # rsync interface (github INWTlab/rsync)
library(RCurl)            # For scp from NCI

# Turn off scientific notatIon

options(scipen = 999)

# Raw data import ---------------------------------------------------------

# There are two data sources for this script: the SA2 shapefile and the gridded
# temperature and precipitation data.

# The SA2 (Statistical Area 2) ESGI shapefile is avilable from the ABS website. 
# If the files are not available in the project directory they will be 
# downloaded

sa2_dir <- here("data-raw/1270055001_sa2_2016_aust_shape/")

if (!dir.exists(sa2_dir)) {
  sa2_url <- read_file(here("sa2_url.txt"))
  sa2_tempfile <- tempfile()
  download.file(sa2_url, sa2_tempfile)
  unzip(sa2_tempfile, exdir = sa2_dir)
}

# Gridded precipitation and maximum temperature data is available on NCI. To 
# access this data you will need access to NCI and also membership of the zv2 
# project.

# Define the target years we are interested in analysing. We want 2000 - 2018,
# but have to buffer by 1 due to the possibility that a week extends beyond the
# calendar year

target_years <- 1999:2019

# Check if required files are in the project directory, and raise an error if
# unavailable. If files are missing they can be copied over using scp with the
# get_raster_data.sh bash script (has to be run manually due to password input).

precip_files <- glue("agcd_v1_precip_total_r005_daily_{target_years}.nc")
tmax_files <- glue("agcd_v1_tmax_mean_r005_daily_{target_years}.nc")

get_missing_files <- function(path, file) {
 file[!file.exists(here(path, file))]
}

missing_files <- c(
  get_missing_files("data-raw/precip", precip_files),
  get_missing_files("data-raw/tmax", tmax_files)
)

if (length(missing_files) > 0) {
  stop(
    "The following files are missing and should be copied from NCI:\n\n",
    paste0(missing_files, collapse = "\n"),
    "\n\nFiles can be copied by running get_raster_data.sh"
  )
}

# Join SA2 data to raster cells -------------------------------------------

# Create a dummy raster with the same dimensions as the actual data. The 
# dimensions are very important as we will be joining by cell number.

rast_crs <- "EPSG:4283" # or maybe EPSG:9001, not 100% sure what the CRS is?

dummy_rast <- rast(
    nrows = 691,
    ncols = 886,
    xmin = 111.975,
    xmax = 156.275,
    ymin = -44.525,
    ymax = -9.975,
    crs = rast_crs
  )

# Define a function to confirm that any raster data matches the dimensions and
# extent of the dummy raster

compare_rasters <- function(rast, dummy_rast) {
  try(
    identical(dim(rast)[1:2], dim(dummy_rast)[1:2]), 
    stop("Raster dimensions don't match the dummy")
  )
  try(
    identical(ext(rast)$vector, ext(dummy_rast)$vector), 
    stop("Raster extent doesn't match the dummy")
  )
  try(
    identical(crs(rast), crs(dummy_rast)), 
    stop("Raster crs doesn't match the dummy")
  )
  return(TRUE)
}

# Read in the SA2 shapefile

sa2 <- st_read(
  here("data-raw/1270055001_sa2_2016_aust_shape/SA2_2016_AUST.shp")
)

# Remove empty geometries from sa2

sa2 <- sa2[!st_is_empty(sa2), , ]

# Reproject sa2 onto the same CRS as the raster data

sa2 <- st_transform(sa2, rast_crs)

# Join SA2 areas to overlapping raster cells and calculate the overlap fraction.
# The result, sa2_cells, is a data frame with SA2_5DIG16 as a key representing
# SA2 areas

sa2_cells <- dummy_rast %>% 
  exact_extract(sa2, include_cell = TRUE, include_cols = "SA2_5DIG16") %>% 
  bind_rows() %>% 
  select(-value) %>% 
  as_tibble()

# Join SA2 to precipitation data ------------------------------------------

# In this section we will define a function to join the SA2 IDs to a year's 
# worth of precipitation data. The output for each year will be saved locally in 
# data/sa2_precip

# Define the join function. This function is a little dodgy as it makes several
# calls outside the function environment. A potential enhancement to this
# script is to generalise this function for an arbitrary variable and 
# arbitrary summary function, and also to ensure that the function doesn't look
# outside its environment.

sa2_precip_join <- function(year,
                            sa2_cells = sa2_cells,
                            dummy_rast = dummy_rast,
                            overwrite = FALSE) {
  
  # Identify the output file path. If this file already exists, simply return
  # the path (unless overwrite = TRUE)
  
  sa2_precip_output_file <- glue("sa2_precip_{year}.feather")
  sa2_precip_output_path <- here("data", "sa2_precip", sa2_precip_output_file)
  
  if (file.exists(sa2_precip_output_path) & overwrite == FALSE) {
    return(sa2_precip_output_path)
  }
  
  # Read in the year's precipitation data
  
  precip_rast_file <- glue("agcd_v1_precip_total_r005_daily_{year}.nc")
  precip_rast_path <- here("data-raw", "precip", precip_rast_file)
  precip_rast <- rast(precip_rast_path)
  
  # Update the raster CRS definition to ensure consistency
  
  crs(precip_rast) <- rast_crs
  
  # Check that precip raster matches dummy raster
  
  compare_rasters(precip_rast, dummy_rast)
  
  # Convert the precipitation in mm to a volume by multiplying by the cell size
  # in km^2
  
  precip_rast_vol <- precip_rast * cellSize(dummy_rast, unit = "km")
  
  # Now, convert the volume raster to a dataframe using cell number as a key.
  # The day of the year is used as the column names. We only need cells which
  # have been matched to a SA2 area
  
  precip_df_wide <- precip_rast_vol %>%
    as.data.frame(cells = TRUE) %>%
    filter(cell %in% sa2_cells$cell)
  
  # Convert to tall format, so that day of year becomes a key column along with
  # cell number
  
  precip_df_tall <- pivot_longer(
    precip_df_wide,
    -cell,
    names_to = "day_of_year",
    names_pattern = "precip_(.*)",
    names_transform = list("day_of_year" = as.integer)
  )
  
  # Join sa2 cells to each day of raster data
  
  sa2_cells_precip <- left_join(sa2_cells, precip_df_tall, by = "cell")
  
  # Sum precipitation volume by SA2 by day
  
  sa2_precip <- sa2_cells_precip %>%
    mutate(cell_precip = coverage_fraction * value) %>%
    group_by(SA2_5DIG16, day_of_year) %>%
    summarise(precip = sum(cell_precip), .groups = "drop")
  
  # Save the year's precipitation data to file
  
  write_feather(sa2_precip, sa2_precip_output_path)
  
  # Return the output path
  
  return(sa2_precip_output_path)
}

# We can now run this function on each year of data

# First, check that the output directory exists

sa2_precip_dir <- here("data", "sa2_precip")

if (!dir.exists(sa2_precip_dir)) {
  dir.create(sa2_precip_dir, recursive = TRUE)
}

# Perform the join for the target years. This will take some time - around
# half an hour on my laptop. Because it is slow, the file will only be created
# if it doesn't already exist. The function returns the output file path

sa2_precip_path <- map_chr(
  target_years, 
  sa2_precip_join,
  sa2_cells = sa2_cells,
  dummy_rast = dummy_rast,
  overwrite = FALSE
)

# Aggregate precipitation data by week ------------------------------------

# We first need to import the precipitation data for each year

sa2_precip <- tibble(year = target_years, sa2_precip_path = sa2_precip_path)

sa2_precip <- sa2_precip %>%
  mutate(sa2_precip_nested = map(sa2_precip_path, read_feather)) %>% 
  select(-sa2_precip_path) %>% 
  unnest(sa2_precip_nested)
  
# Convert day of year to a date, and calculate the first day of each week

sa2_precip <- sa2_precip %>% 
  mutate(date = as_date(strptime(paste(year, day_of_year), format="%Y %j"))) %>% 
  mutate(week_starting = floor_date(date, "weeks", week_start = 1))

# Sum precipitation by week

sa2_precip_weekly <- sa2_precip %>% 
  group_by(SA2_5DIG16, week_starting) %>% 
  summarise(weekly_precip = sum(precip), .groups = "drop")

# Only keep weeks which overlap with the target years

week_filter <- function(sa2_data, target_years) {
  sa2_data %>% 
    filter(year(week_starting %m+% days(6)) > min(target_years)) %>% 
    filter(year(week_starting) < max(target_years))
}

sa2_precip_weekly <- week_filter(sa2_precip_weekly, target_years)

# Now we just need to add back in the rest of the statistical area specification

sa2_spec <- as_tibble(sa2) %>% 
  select(SA2_MAIN16, SA2_5DIG16, SA2_NAME16, AREASQKM16)

sa2_precip_weekly <- inner_join(sa2_spec, sa2_precip_weekly, by = "SA2_5DIG16")

# To calculate the precipitation in mm, we need to divide by the area of the 
# statistical area in km2

sa2_precip_weekly <- sa2_precip_weekly %>% 
  mutate(weekly_precip = weekly_precip / AREASQKM16) %>% 
  arrange(SA2_MAIN16, week_starting)

# Save the precipitation output as a safeguard

write_csv(sa2_precip_weekly, here("data/sa2_precip_weekly.csv"))

# Join SA2 to temperature data --------------------------------------------

# We now more or less repeat the process for the temperature data

sa2_tmax_join <- function(year, 
                          sa2_cells = sa2_cells, 
                          dummy_rast = dummy_rast,
                          overwrite = FALSE
                          ) {
  
  # Identify the output file path. If this file already exists, simply return
  # the path (unless overwrite = TRUE)
  
  sa2_tmax_file <- glue("sa2_tmax_{year}.feather")
  sa2_tmax_output_path <- here("data", "sa2_tmax", sa2_tmax_file)
  
  if (file.exists(sa2_tmax_output_path) & overwrite == FALSE) {
    return(sa2_tmax_output_path)
  }
  
  # Read in the year's tempearature data
  
  tmax_rast_file <- glue("agcd_v1_tmax_mean_r005_daily_{year}.nc")
  tmax_rast_path <- here("data-raw", "tmax", tmax_rast_file)
  tmax_rast <- rast(tmax_rast_path)
  
  # Update the raster CRS definition to ensure consistency
  
  crs(tmax_rast) <- rast_crs
  
  # Check that tmax raster matches dummy raster
  
  compare_rasters(tmax_rast, dummy_rast)
  
  # Now, convert the volume raster to a dataframe using cell number as a key. 
  # The day of the year is used as the column names. We only need cells which
  # have been matched to a SA2 area
  
  tmax_df_wide <- tmax_rast %>% 
    as.data.frame(cells = TRUE) %>% 
    filter(cell %in% sa2_cells$cell)
  
  # Convert to tall format, so that day of year becomes a key column along with
  # cell number
  
  tmax_df_tall <- tmax_df_wide %>% 
    pivot_longer(
      -cell,
      names_to = "day_of_year",
      names_pattern = "tmax_(.*)",
      names_transform = list("day_of_year" = as.integer),
      values_to = "tmax"
    )
  
  # Join SA2 cells to each day of raster data
  
  sa2_cells_tmax <- left_join(sa2_cells, tmax_df_tall, by = "cell")
  
  # Keep only the highest temperature identified in each SA2
  
  sa2_tmax <- sa2_cells_tmax %>% 
    group_by(SA2_5DIG16, day_of_year) %>% 
    summarise(tmax = max(tmax), .groups = "drop") %>% 
    select(SA2_5DIG16, day_of_year, tmax)
  
  # Save the year's tmax data to file
  
  write_feather(sa2_tmax, sa2_tmax_output_path)
  
  # Return the output path
  
  return(sa2_tmax_output_path)
}

# First, check that the output directory exists

sa2_tmax_dir <- here("data", "sa2_tmax")

if (!dir.exists(sa2_tmax_dir)) {
  dir.create(sa2_tmax_dir, recursive = TRUE)
}

# Perform the join for the target. This will take some time - around
# half an hour on my laptop. Because it is slow, the file will only be created
# if it doesn't already exist. The function returns the output file path

sa2_tmax_path <- map_chr(
  target_years, 
  sa2_tmax_join, 
  sa2_cells = sa2_cells, 
  dummy_rast = dummy_rast,
  overwrite = FALSE
)

# Aggregate tmax data by week ---------------------------------------------

# We first need to import the tmax data for each year

sa2_tmax <- tibble(year = target_years, sa2_tmax_path = sa2_tmax_path)

sa2_tmax <- sa2_tmax %>%
  mutate(sa2_tmax_nested = map(sa2_tmax_path, read_feather)) %>% 
  select(-sa2_tmax_path) %>% 
  unnest(sa2_tmax_nested)

# Convert day of year to a date, and calculate the first day of each week

sa2_tmax <- sa2_tmax %>% 
  mutate(date = as_date(strptime(paste(year, day_of_year), format="%Y %j"))) %>% 
  mutate(week_starting = floor_date(date, "weeks", week_start = 1))

# Find the maximum value of tmax by week

sa2_tmax_weekly <- sa2_tmax %>% 
  group_by(SA2_5DIG16, week_starting) %>% 
  summarise(weekly_tmax = max(tmax), .groups = "drop")

# Only keep weeks which overlap with the target years

sa2_tmax_weekly <- week_filter(sa2_tmax_weekly, target_years)

# Now we just need to add back in the rest of the statistical area specification

sa2_tmax_weekly <- inner_join(sa2_spec, sa2_tmax_weekly, by = "SA2_5DIG16")

# Save the weekly tmax output as a safeguard

write_csv(sa2_tmax_weekly, here("data/sa2_tmax_weekly.csv"))

# Combine precip and tmax data --------------------------------------------

# Combine the temperature and precipitation data into a single table

sa2_temp_precip <- left_join(
  sa2_precip_weekly, 
  sa2_tmax_weekly, 
  by = c(
    "SA2_MAIN16", 
    "SA2_5DIG16", 
    "SA2_NAME16", 
    "AREASQKM16", 
    "week_starting"
  )
)

# Ensure that that the table is ordered by SA2 and then week, and also round
# the values to one decimal point

sa2_temp_precip <- sa2_temp_precip %>% 
  arrange(SA2_MAIN16, week_starting) %>% 
  mutate(
    weekly_precip = round(weekly_precip, 1),
    weekly_tmax = round(weekly_tmax, 1)
  )

# Finally, perform a few cursory data quality checks

sa2_temp_precip %>% 
  assert(not_na, everything()) %>% 
  assert(within_bounds(0, 1000), weekly_precip) %>% 
  assert(within_bounds(-23, 54), weekly_tmax) %>%
  verify(wday(week_starting, week_start = 1) == 1)

# Save the final output

write_csv(sa2_temp_precip, here("data/sa2_temp_precip.csv"))
