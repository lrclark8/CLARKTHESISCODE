
# ============================================================
# WHITE SHRIMP ANALYSIS
# Maryland DNR Coastal Bays Fisheries Investigation
# Standardized September–October sampling, 1990–2018
# Any AI produced code was thoroughly checked and corrected by the author prior to implementation.
# ============================================================


# ============================================================
# 1. PACKAGES 
# ============================================================

library(tidyverse)
library(readxl)
library(lme4)
library(glmmTMB)
library(DHARMa)
library(ggeffects)
library(patchwork)
library(MuMIn)
library(patchwork)
library(janitor)
library(lubridate)
library(sf)
library(ggspatial)
library(prettymapr)
library(broom)
library(car)


# ============================================================
# 2.IMPORT RAW DATA
# ============================================================

MDNR1 <- read_excel(
  "MDNR_CBFI Data (1990-2012) DATA ORIGINAL 1.xlsx"
)

MDNR2 <- read.csv(
  "MDNR_CBFI Data (2013-2017) DATA ORIGINAL(MDNR_CBFI Data (2013-2017) DATA).csv"
)


MDNR3 <- read.csv ("CBFI_17_18 1(in).csv")


# ============================================================
# 3. PREPARE STANDARDIZED SURVEY DATA
# ============================================================

# ============================================================
# CONVERT MDNR COORDINATES
# DEGREES AND DECIMAL MINUTES TO DECIMAL DEGREES
# ============================================================

convert_mdnr_coord <- function(x, direction = "latitude") {
  
  x <- as.numeric(x)
  
  degrees <- dplyr::case_when(
    x >= 100000 & x < 1000000 ~ floor(x / 10000),
    x >= 1000000 & x < 10000000 ~ floor(x / 100000),
    TRUE ~ NA_real_
  )
  
  minutes <- dplyr::case_when(
    x >= 100000 & x < 1000000 ~ (x %% 10000) / 100,
    x >= 1000000 & x < 10000000 ~ (x %% 100000) / 1000,
    TRUE ~ NA_real_
  )
  
  result <- degrees + minutes / 60
  
  # Reject invalid minute values
  result[!is.na(minutes) & minutes >= 60] <- NA_real_
  
  if (direction == "longitude") {
    result <- -abs(result)
  }
  
  result
}

# ------------------------------------------------------------
# 3A. 1990–2012
# ------------------------------------------------------------

# Create one row per sampling event
sampling_events_1990_2012 <- MDNR1 %>%
  group_by(
    Year,
    Month,
    `Sample Date`,
    `Sample Time`,
    `Site code`,
    `Gear Code`
  ) %>%
  summarise(
    # Sampling effort
    area_covered = first(`Area Covered`) * 0.092903,
    distance_haul = first(`Distance of Haul_m`),
    duration = first(Duration),
    
    # Environmental variables
    temp_surface = first(`Water Temp Surface`),
    temp_bottom = first(`Water Temp Bottom`),
    
    salinity_surface = first(`Salinity Surface`),
    salinity_bottom = first(`Salinity Bottom`),
    
    dissolved_oxygen_surface = first(`Dissolved Oxygen Surface`),
    dissolved_oxygen_bottom = first(`Dissolved Oxygen Bottom`),
    
    secchi = first(Secchi),
    
    # Depth
    depth_start = first(`Depth Start`),
    depth_stop = first(`Depth Stop`),
    
    # Geographic coordinates
    latitude_start = first(`Latitude Start`),
    latitude_stop = first(`Latitude Stop`),
    
    longitude_start = first(`Longitude Start`),
    longitude_stop = first(`Longitude Stop`),
    
    .groups = "drop"
  ) %>%
  mutate(
    # Sampling-event identifier
    sample_id = paste(
      `Sample Date`,
      `Sample Time`,
      `Site code`,
      `Gear Code`,
      sep = "_"
    ),
    
    # Convert depth from feet to meters
    depth_start_m = depth_start * 0.3048,
    depth_stop_m = depth_stop * 0.3048,
    
    # Calculate mean sampling depth in meters
    mean_depth = rowMeans(
      cbind(depth_start_m, depth_stop_m),
      na.rm = TRUE
    ),
    
    mean_depth = ifelse(
      is.nan(mean_depth),
      NA_real_,
      mean_depth
    ),
    
    # Convert original start and stop coordinates
    latitude_start_dd = convert_mdnr_coord(latitude_start),
    latitude_stop_dd = convert_mdnr_coord(latitude_stop),
    
    longitude_start_dd = convert_mdnr_coord(
      longitude_start, "longitude"
    ),
    
    longitude_stop_dd = convert_mdnr_coord(
      longitude_stop, "longitude"
    ),
    
    # Calculate midpoint coordinates
    latitude = rowMeans(
      cbind(latitude_start_dd, latitude_stop_dd),
      na.rm = TRUE
    ),
    
    longitude = rowMeans(
      cbind(longitude_start_dd, longitude_stop_dd),
      na.rm = TRUE
    ),
    
    # Replace NaN values
    latitude = ifelse(is.nan(latitude), NA_real_, latitude),
    longitude = ifelse(is.nan(longitude), NA_real_, longitude)
    
    )
  

# Sum white shrimp abundance within each sampling event
white_shrimp_positive_1990_2012 <- MDNR1 %>%
  filter(
    `Speciment name` ==
      "White Shrimp (Litopenaeus setiferus)"
  ) %>%
  mutate(
    abundance = as.numeric(
      `Total Number of organims collected`
    )
  ) %>%
  group_by(
    Year,
    Month,
    `Sample Date`,
    `Sample Time`,
    `Site code`,
    `Gear Code`
  ) %>%
  summarise(
    abundance = sum(
      abundance,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# Join shrimp abundance to all sampling events
white_shrimp_1990_2012 <- sampling_events_1990_2012 %>%
  left_join(
    white_shrimp_positive_1990_2012,
    by = c(
      "Year",
      "Month",
      "Sample Date",
      "Sample Time",
      "Site code",
      "Gear Code"
    )
  ) %>%
  mutate(
    abundance = replace_na(abundance, 0),
    presence = as.integer(abundance > 0)
  )


# Restrict to standardized fall sampling
fall_1990_2012 <- white_shrimp_1990_2012 %>%
  filter(
    Month %in% c(9, 10),
    `Gear Code` %in% c("BTRW", "BCHS100")
  ) %>%
  transmute(
    sample_id,
    Year,
    Month,
    
    SampleDate = as.Date(`Sample Date`),
    SampleTime = `Sample Time`,
    Site = as.character(`Site code`),
    Gear = `Gear Code`,
    
    # White shrimp response variables
    presence,
    abundance,
    
    # Environmental variables
    temp_surface,
    temp_bottom,
    salinity_surface,
    salinity_bottom,
    dissolved_oxygen_surface,
    dissolved_oxygen_bottom,
    secchi,
    
    # Sampling effort
    area_covered,
    distance_haul,
    duration,
    
    # Depth
    depth_start,
    depth_stop,
    mean_depth,
    
    # Geographic coordinates
    latitude,
    longitude
  )


# ------------------------------------------------------------
# 3B. 2013–2016
# ------------------------------------------------------------

# Create sampling-event identifiers
MDNR2 <- MDNR2 %>%
  mutate(
    sample_id = paste(
      SampleDate,
      SampleTime,
      SiteNumber,
      GearCode,
      sep = "_"
    )
  )


# Create one row per sampling event
sampling_events_2013_2016 <- MDNR2 %>%
  filter(Year <= 2016) %>%
  group_by(sample_id) %>%
  summarise(
    Year = first(Year),
    Month = first(Month),
    
    SampleDate = first(SampleDate),
    SampleTime = first(SampleTime),
    
    Site = as.character(first(SiteNumber)),
    Gear = first(GearCode),
    
    # Environmental variables
    temp_surface = first(WaterTempSurface),
    temp_bottom = first(WaterTempBottom),
    
    salinity_surface = first(SalinitySurface),
    salinity_bottom = first(SalinityBottom),
    
    dissolved_oxygen_surface =
      first(DissolvedOxygenSurface),
    
    dissolved_oxygen_bottom =
      first(DissolvedOxygenBottom),
    
    secchi = first(Secchi),
    
    # Sampling effort
    area_covered = first(AreaCovered) * 0.092903,
    distance_haul = first(DistanceofHaul_m),
    duration = first(Duration),
    
    # Depth
    depth_start = first(DepthStart),
    depth_stop = first(DepthStop),
    
    # Geographic coordinates
    latitude_start = first(LatitudeStart),
    latitude_stop = first(LatitudeStop),
    
    longitude_start = first(LongitudeStart),
    longitude_stop = first(LongitudeStop),
    
    .groups = "drop"
  ) %>%
  mutate(
    # Convert depth from feet to meters
    depth_start_m = depth_start * 0.3048,
    depth_stop_m = depth_stop * 0.3048,
    
    # Calculate mean sampling depth in meters
    mean_depth = rowMeans(
      cbind(depth_start_m, depth_stop_m),
      na.rm = TRUE
    ),
    
    mean_depth = ifelse(
      is.nan(mean_depth),
      NA_real_,
      mean_depth
    ),
    
    # Convert original start and stop coordinates
    latitude_start_dd = convert_mdnr_coord(latitude_start),
    latitude_stop_dd = convert_mdnr_coord(latitude_stop),
    
    longitude_start_dd = convert_mdnr_coord(
      longitude_start, "longitude"
    ),
    
    longitude_stop_dd = convert_mdnr_coord(
      longitude_stop, "longitude"
    ),
    
    # Calculate midpoint coordinates
    latitude = rowMeans(
      cbind(latitude_start_dd, latitude_stop_dd),
      na.rm = TRUE
    ),
    
    longitude = rowMeans(
      cbind(longitude_start_dd, longitude_stop_dd),
      na.rm = TRUE
    ),
    
    # Replace NaN values
    latitude = ifelse(is.nan(latitude), NA_real_, latitude),
    longitude = ifelse(is.nan(longitude), NA_real_, longitude)
    )
  


# Sum white shrimp abundance within each sampling event
white_shrimp_positive_2013_2016 <- MDNR2 %>%
  filter(
    Year <= 2016,
    SpecimenName ==
      "White shrimp (Litopenaeus setiferus)"
  ) %>%
  group_by(sample_id) %>%
  summarise(
    abundance = sum(
      SpeciesCount,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# Join shrimp abundance to all sampling events
white_shrimp_2013_2016 <- sampling_events_2013_2016 %>%
  left_join(
    white_shrimp_positive_2013_2016,
    by = "sample_id"
  ) %>%
  mutate(
    abundance = replace_na(abundance, 0),
    presence = as.integer(abundance > 0)
  )


# Restrict to standardized fall sampling
fall_2013_2016 <- white_shrimp_2013_2016 %>%
  filter(
    Month %in% c(9, 10),
    Gear %in% c("BTRW", "BCHS100")
  ) %>%
  transmute(
    sample_id,
    Year,
    Month,
    
    SampleDate = as.Date(
      SampleDate,
      format = "%m/%d/%Y"
    ),
    
    SampleTime,
    Site,
    Gear,
    
    presence,
    abundance,
    
    temp_surface,
    temp_bottom,
    salinity_surface,
    salinity_bottom,
    dissolved_oxygen_surface,
    dissolved_oxygen_bottom,
    secchi,
    
    area_covered,
    distance_haul,
    duration,
    
    depth_start,
    depth_stop,
    mean_depth,
    
    latitude,
    longitude
  )


# ------------------------------------------------------------
# 3C. 2017–2018
# ------------------------------------------------------------

# Create sampling-event identifiers
MDNR3 <- MDNR3 %>%
  mutate(
    sample_id = paste(
      SampleDate,
      SampleTime,
      SiteNumber,
      GearCode,
      sep = "_"
    )
  )


# Create one row per sampling event
sampling_events_2017_2018 <- MDNR3 %>%
  filter(Year %in% c(2017, 2018)) %>%
  group_by(sample_id) %>%
  summarise(
    Year = first(Year),
    Month = first(Month),
    
    SampleDate = first(SampleDate),
    SampleTime = first(SampleTime),
    
    Site = as.character(first(SiteNumber)),
    Gear = first(GearCode),
    
    # Environmental variables
    temp_surface = first(WaterTempSurface),
    temp_bottom = first(WaterTempBottom),
    
    salinity_surface = first(SalinitySurface),
    salinity_bottom = first(SalinityBottom),
    
    dissolved_oxygen_surface =
      first(DissolvedOxygenSurface),
    
    dissolved_oxygen_bottom =
      first(DissolvedOxygenBottom),
    
    secchi = first(Secchi),
    
    # Sampling effort
    area_covered = first(AreaCovered) * 0.092903,
    distance_haul = first(DistanceofHaul_m),
    duration = first(Duration),
    
    # Depth
    depth_start = first(DepthStart),
    depth_stop = first(DepthStop),
    
    # Geographic coordinates
    latitude_start = first(LatitudeStart),
    latitude_stop = first(LatitudeStop),
    
    longitude_start = first(LongitudeStart),
    longitude_stop = first(LongitudeStop),
    
    .groups = "drop"
  ) %>%
  mutate(
    # Mean sampling depth
    depth_start_m = depth_start * 0.3048,
    depth_stop_m = depth_stop * 0.3048,
    
    # Calculate mean sampling depth in meters
    mean_depth = rowMeans(
      cbind(depth_start_m, depth_stop_m),
      na.rm = TRUE
    ),
    
    mean_depth = ifelse(
      is.nan(mean_depth),
      NA_real_,
      mean_depth
    ),
    
    mean_depth = ifelse(
      is.nan(mean_depth),
      NA_real_,
      mean_depth
    ),
    
    # Convert original start and stop coordinates
    latitude_start_dd = convert_mdnr_coord(latitude_start),
    latitude_stop_dd = convert_mdnr_coord(latitude_stop),
    
    longitude_start_dd = convert_mdnr_coord(
      longitude_start, "longitude"
    ),
    
    longitude_stop_dd = convert_mdnr_coord(
      longitude_stop, "longitude"
    ),
    
    # Calculate midpoint coordinates
    latitude = rowMeans(
      cbind(latitude_start_dd, latitude_stop_dd),
      na.rm = TRUE
    ),
    
    longitude = rowMeans(
      cbind(longitude_start_dd, longitude_stop_dd),
      na.rm = TRUE
    ),
    
    # Replace NaN values
    latitude = ifelse(is.nan(latitude), NA_real_, latitude),
    longitude = ifelse(is.nan(longitude), NA_real_, longitude)
    )
  


# Sum white shrimp abundance within each sampling event
white_shrimp_positive_2017_2018 <- MDNR3 %>%
  filter(
    Year %in% c(2017, 2018),
    SpecimenName ==
      "White shrimp (Litopenaeus setiferus)"
  ) %>%
  group_by(sample_id) %>%
  summarise(
    abundance = sum(
      SpeciesCount,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# Join shrimp abundance to all sampling events
white_shrimp_2017_2018 <- sampling_events_2017_2018 %>%
  left_join(
    white_shrimp_positive_2017_2018,
    by = "sample_id"
  ) %>%
  mutate(
    abundance = replace_na(abundance, 0),
    presence = as.integer(abundance > 0)
  )


# Restrict to standardized fall sampling
fall_2017_2018 <- white_shrimp_2017_2018 %>%
  filter(
    Month %in% c(9, 10),
    Gear %in% c("BTRW", "BCHS100")
  ) %>%
  transmute(
    sample_id,
    Year,
    Month,
    
    SampleDate = as.Date(
      SampleDate,
      format = "%m/%d/%Y"
    ),
    
    SampleTime,
    Site,
    Gear,
    
    presence,
    abundance,
    
    temp_surface,
    temp_bottom,
    salinity_surface,
    salinity_bottom,
    dissolved_oxygen_surface,
    dissolved_oxygen_bottom,
    secchi,
    
    area_covered,
    distance_haul,
    duration,
    
    depth_start,
    depth_stop,
    mean_depth,
    
    latitude,
    longitude
  )


# Check coordinates in the sampling-event datasets
sampling_events_1990_2012 %>%
  select(Year, latitude_start, latitude_start_dd, latitude) %>%
  head(5)

sampling_events_2013_2016 %>%
  select(Year, latitude_start, latitude_start_dd, latitude) %>%
  head(5)

sampling_events_2017_2018 %>%
  select(Year, latitude_start, latitude_start_dd, latitude) %>%
  head(5)

# Check coordinates in the standardized datasets
fall_1990_2012 %>%
  select(Year, latitude, longitude) %>%
  head(5)

fall_2013_2016 %>%
  select(Year, latitude, longitude) %>%
  head(5)

fall_2017_2018 %>%
  select(Year, latitude, longitude) %>%
  head(5)



# ------------------------------------------------------------
# 3D. COMBINE ALL THREE SURVEY PERIODS
# ------------------------------------------------------------

fall_1990_2018 <- bind_rows(
  fall_1990_2012,
  fall_2013_2016,
  fall_2017_2018
) %>%
  arrange(
    Year,
    Month,
    SampleDate,
    Site,
    Gear
  )


# ============================================================
# 3E. ADD WINTER NAO INDEX
# ============================================================

nao_data <- read.table(
  "~/UMES:MEES/thesis data/nao_index.tim.txt",
  skip = 9,
  header = FALSE,
  col.names = c("year", "month", "NAO")
)

# Calculate mean winter NAO by winter year
winter_NAO <- nao_data %>%
  mutate(
    winter_year = ifelse(month == 12, year + 1, year)
  ) %>%
  filter(
    month %in% c(12, 1, 2),
    winter_year >= 1990,
    winter_year <= 2018,
    NAO != -99.90
  ) %>%
  group_by(winter_year) %>%
  summarise(
    winter_NAO = mean(NAO, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    winter_NAO_z = as.numeric(scale(winter_NAO))
  )

fall_1990_2018 <- fall_1990_2018 %>%
  left_join(
    winter_NAO %>%
      select(winter_year, winter_NAO, winter_NAO_z),
    by = c("Year" = "winter_year")
  )


# ============================================================
# 3F. CALCULATE SAME-YEAR SPRING TEMPERATURE
# ============================================================

all_sampling_events <- bind_rows(
  sampling_events_1990_2012,
  sampling_events_2013_2016,
  sampling_events_2017_2018
)

spring_temp <- all_sampling_events %>%
  filter(
    Month %in% c(3, 4, 5),
    !is.na(temp_surface)
  ) %>%
  group_by(Year) %>%
  summarise(
    spring_temp_same_year = mean(
      temp_surface,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# Join spring temperature to fall events from the same year
fall_1990_2018 <- fall_1990_2018 %>%
  left_join(spring_temp, by = "Year")


# ============================================================
# 3G. VERIFY NAO AND SAME-YEAR SPRING TEMPERATURE
# ============================================================

# Check annual NAO values
fall_1990_2018 %>%
  distinct(Year, winter_NAO, winter_NAO_z) %>%
  arrange(Year) %>%
  print(n = Inf)

# Check annual lagged spring temperatures
fall_1990_2018 %>%
  distinct(Year, spring_temp_same_year) %>%
  arrange(Year) %>%
  print(n = Inf)

# Count missing values
fall_1990_2018 %>%
  summarise(
    missing_NAO = sum(is.na(winter_NAO_z)),
    missing_spring_temp = sum(is.na(spring_temp_same_year))
  )

# Identify years with missing lagged spring temperature
fall_1990_2018 %>%
  filter(is.na(spring_temp_same_year)) %>%
  count(Year, name = "missing_events")


# ============================================================
# 3H. VERIFY COMBINED DATASET
# ============================================================

nrow(fall_1990_2018)

sum(fall_1990_2018$presence)

summary(fall_1990_2018$latitude)
summary(fall_1990_2018$longitude)

sum(is.na(fall_1990_2018$latitude))
sum(is.na(fall_1990_2018$longitude))

# Check annual sampling effort and detections
fall_1990_2018 %>%
  group_by(Year) %>%
  summarise(
    samples = n(),
    positives = sum(presence),
    percent_positive = mean(presence) * 100,
    total_abundance = sum(abundance),
    .groups = "drop"
  ) %>%
  print(n = Inf)

# Check coordinate coverage by survey period
fall_1990_2018 %>%
  mutate(
    period = case_when(
      Year <= 2012 ~ "1990–2012",
      Year <= 2016 ~ "2013–2016",
      TRUE ~ "2017–2018"
    )
  ) %>%
  group_by(period) %>%
  summarise(
    samples = n(),
    missing_latitude = sum(is.na(latitude)),
    missing_longitude = sum(is.na(longitude)),
    .groups = "drop"
  )


# ============================================================
# 3I. COORDINATE QUALITY CONTROL
# ============================================================

coordinate_quality <- fall_1990_2018 %>%
  group_by(Year) %>%
  summarise(
    total_samples = n(),
    missing_coordinates = sum(
      is.na(latitude) | is.na(longitude)
    ),
    usable_coordinates = sum(
      !is.na(latitude) & !is.na(longitude)
    ),
    positive_detections = sum(presence),
    .groups = "drop"
  )

print(coordinate_quality, n = Inf)


# ============================================================
# 3J. VERIFY CORRECTED 1992 COORDINATES
# ============================================================

fall_1990_2018 %>%
  filter(
    Year %in% c(1992, 2006),
    Site %in% c("T001", "T002", "S005", "S007")
  ) %>%
  arrange(Site, Year, Month) %>%
  transmute(
    Site,
    Year,
    Month,
    latitude = sprintf("%.6f", latitude),
    longitude = sprintf("%.6f", longitude)
  ) %>%
  print(n = Inf, width = Inf)


# ============================================================
# 3K. WHITE SHRIMP DETECTION COORDINATE AVAILABILITY
# ============================================================

detection_coordinate_check <- fall_1990_2018 %>%
  mutate(
    coordinate_available =
      !is.na(latitude) & !is.na(longitude)
  ) %>%
  group_by(
    period = if_else(
      Year < 2008,
      "1990–2007",
      "2008–2018"
    )
  ) %>%
  summarise(
    total_events = n(),
    positive_detections = sum(presence),
    detections_with_coordinates = sum(
      presence == 1 & coordinate_available
    ),
    detections_without_coordinates = sum(
      presence == 1 & !coordinate_available
    ),
    .groups = "drop"
  )

print(detection_coordinate_check) 



# ============================================================
# 3L. ADD REFERENCE SAMPLING-SITE COORDINATES
# ============================================================

# Trawl survey sites
trawl_sites <- tibble(
  Site = sprintf("T%03d", 1:20),
  Gear = "BTRW",
  lat_min = c(
    26.243, 25.859, 23.919, 23.527, 24.425,
    21.291, 22.357, 19.418, 17.852, 14.506,
    13.024, 15.281, 10.213, 8.447, 6.370,
    4.545, 3.140, 5.257, 3.125, 1.328
  ),
  lon_min = c(
    4.747, 6.108, 5.429, 7.327, 10.514,
    8.781, 5.776, 6.018, 7.310, 9.301,
    12.396, 11.603, 13.989, 16.043, 12.454,
    17.025, 16.116, 19.494, 21.110, 20.057
  )
)

# Beach seine survey sites
seine_sites <- tibble(
  Site = sprintf("S%03d", 1:19),
  Gear = "BCHS100",
  lat_min = c(
    23.273, 24.749, 24.824, 20.388, 21.928,
    23.627, 22.557, 14.554, 19.132, 17.367,
    13.227, 12.579, 9.340, 6.432, 7.000,
    4.482, 2.162, 2.038, 18.774
  ),
  lon_min = c(
    4.380, 4.264, 6.044, 5.390, 7.017,
    6.797, 4.301, 8.581, 6.174, 7.977,
    12.054, 14.921, 16.426, 12.404, 17.578,
    17.597, 22.190, 16.619, 9.414
  )
)

# Combine and convert degrees/minutes to decimal degrees
site_coordinates <- bind_rows(
  trawl_sites,
  seine_sites
) %>%
  mutate(
    latitude = 38 + lat_min / 60,
    longitude = -(75 + lon_min / 60)
  ) %>%
  select(Site, Gear, latitude, longitude)

# Verify the reference table
print(site_coordinates, n = Inf)




# Ensure reference gear codes match the dataset
site_coordinates <- site_coordinates %>%
  mutate(
    Gear = if_else(Gear == "BCHS", "BCHS100", Gear)
  )

# Verify unique site and gear combinations
stopifnot(
  !anyDuplicated(site_coordinates[c("Gear", "Site")])
)

# Join reference coordinates without changing original data
fall_mapped <- fall_1990_2018 %>%
  left_join(
    site_coordinates %>%
      select(
        Gear, Site,
        ref_latitude = latitude,
        ref_longitude = longitude
      ),
    by = c("Gear", "Site"),
    relationship = "many-to-one"
  ) %>%
  mutate(
    # Retain complete original coordinate pairs
    original_valid =
      !is.na(latitude) &
      !is.na(longitude) &
      between(latitude, 38.0, 38.5) &
      between(longitude, -75.5, -75.0),
    
    map_latitude = if_else(
      original_valid,
      latitude,
      ref_latitude
    ),
    
    map_longitude = if_else(
      original_valid,
      longitude,
      ref_longitude
    ),
    
    coordinate_source = case_when(
      original_valid ~ "Original",
      !is.na(ref_latitude) &
        !is.na(ref_longitude) ~ "Reference",
      TRUE ~ "Missing"
    )
  )

# Check positive detections by coordinate source
fall_mapped %>%
  filter(presence == 1) %>%
  mutate(
    Period = if_else(
      Year < 2008,
      "1990–2007",
      "2008–2018"
    )
  ) %>%
  count(Period, Gear, coordinate_source) %>%
  print(n = Inf)



# ============================================================
# 4. AUDIT FINAL 1990–2018 DATASET
# ============================================================

final_data_check <- fall_1990_2018 %>%
  summarise(
    total_events = n(),
    positive_events = sum(presence),
    zero_events = sum(presence == 0),
    total_abundance = sum(abundance),
    n_years = n_distinct(Year),
    first_year = min(Year),
    last_year = max(Year),
    n_sites = n_distinct(Site)
  )

final_data_check

# Sampling effort by gear
table(fall_1990_2018$Gear)

# Sampling effort by month
table(fall_1990_2018$Month)

# Annual occurrence summary
annual_1990_2018 <- fall_1990_2018 %>%
  group_by(Year) %>%
  summarise(
    samples = n(),
    positives = sum(presence),
    percent_positive = mean(presence) * 100,
    total_abundance = sum(abundance),
    .groups = "drop"
  )

print(annual_1990_2018, n = Inf)


# ============================================================
# 5. BTRW-ONLY DATASET
# ============================================================

fall_btrw <- fall_1990_2018 %>%
  filter(
    Gear == "BTRW",
    Month %in% c(9, 10)
  ) %>%
  droplevels()

# Dataset excluding observations without lagged spring temperature
fall_btrw_same_year <- fall_btrw %>%
  filter(!is.na(spring_temp_same_year))

nrow(fall_btrw)
nrow(fall_btrw_same_year)

sum(is.na(fall_btrw_same_year$spring_temp_same_year))

range(fall_btrw_same_year$Year)


# ============================================================
# 6. BTRW SAMPLING EFFORT
# ============================================================

brtw_effort <- fall_btrw %>%
  group_by(Year, Month) %>%
  summarise(
    sampling_events = n(),
    unique_sites = n_distinct(Site),
    positive_events = sum(presence == 1, na.rm = TRUE),
    total_shrimp = sum(abundance, na.rm = TRUE),
    .groups = "drop"
  )

print(brtw_effort, n = Inf)


# ============================================================
# 7. NAO AND SPRING TEMPERATURE DATA CHECKS
# ============================================================

# Check for duplicate columns
grep(
  "NAO|spring_temp",
  names(fall_1990_2018),
  value = TRUE
)

# Check missing environmental values in BTRW samples
occ_env_check <- fall_btrw %>%
  summarise(
    total_events = n(),
    missing_temp = sum(is.na(temp_surface)),
    missing_salinity = sum(is.na(salinity_surface)),
    missing_depth = sum(is.na(mean_depth)),
    missing_NAO = sum(is.na(winter_NAO_z)),
    missing_spring_temp = sum(is.na(spring_temp_same_year)),
    missing_DO = sum(is.na(dissolved_oxygen_bottom))
  )

occ_env_check


# ============================================================
# 8. PREPARE ENVIRONMENTAL OCCURRENCE DATA
# ============================================================


occ_btrw <- fall_btrw %>%
  filter(
    !is.na(temp_surface),
    !is.na(salinity_surface),
    !is.na(mean_depth),
    !is.na(spring_temp_same_year)
  ) %>%
  mutate(
    year_decade = (Year - 2008) / 10,
    Month_f = factor(
      as.character(Month),
      levels = c("9", "10"),
      labels = c("September", "October")
    ),
    Site_f = factor(Site),
    temp_z = as.numeric(scale(temp_surface)),
    salinity_z = as.numeric(scale(salinity_surface)),
    depth_z = as.numeric(scale(mean_depth)),
    spring_temp_same_year_z = as.numeric(
      scale(spring_temp_same_year)
    )
  )

occ_btrw %>%
  summarise(
    total_events = n(),
    missing_spring = sum(is.na(spring_temp_same_year_z)),
    years = n_distinct(Year),
    sites = n_distinct(Site)
  )


# Check observations by year and month
occ_btrw %>%
  group_by(Year, Month_f) %>%
  summarise(
    sampling_events = n(),
    positive_events = sum(presence == 1),
    occurrence_percent = mean(presence == 1) * 100,
    .groups = "drop"
  ) %>%
  arrange(Year, Month_f)

# Correlations among continuous predictors
occ_correlations <- cor(
  occ_btrw %>%
    select(
      year_decade,
      temp_z,
      salinity_z,
      depth_z,
      spring_temp_same_year_z
    ),
  use = "complete.obs"
)

round(occ_correlations, 3)


# Separate dataset with complete NAO observations
occ_btrw_NAO <- occ_btrw %>%
  filter(!is.na(winter_NAO_z))

# Check NAO sample size
nrow(occ_btrw)
nrow(occ_btrw_NAO)

sum(is.na(occ_btrw$winter_NAO_z))

# Check observations by year
occ_btrw %>% count(Year)
occ_btrw_NAO %>% count(Year)

# Correlations including NAO, using NAO-complete observations
occ_correlations_NAO <- cor(
  occ_btrw_NAO %>%
    select(
      year_decade,
      temp_z,
      salinity_z,
      depth_z,
      spring_temp_same_year_z,
      winter_NAO_z
    ),
  use = "complete.obs"
)

round(occ_correlations_NAO, 3)

## Separated NAO ##
occ_btrw_NAO <- occ_btrw %>%
  filter(!is.na(winter_NAO_z))

# Main occurrence dataset
nrow(occ_btrw)

# NAO-complete dataset
nrow(occ_btrw_NAO)

# Number of observations missing NAO
sum(is.na(occ_btrw$winter_NAO_z))

# Check observations by year
occ_btrw %>% count(Year)

occ_btrw_NAO %>% count(Year)



# ------------------------------------------------------------
# Check environmental dataset
# ------------------------------------------------------------


occ_env_check <- fall_1990_2018 %>%
  filter(Gear == "BTRW") %>%
  summarise(
    total_events = n(),
    missing_temp = sum(is.na(temp_surface)),
    missing_salinity = sum(is.na(salinity_surface)),
    missing_depth = sum(is.na(mean_depth)),
    missing_NAO = sum(is.na(winter_NAO_z)),
    missing_DO = sum(is.na(dissolved_oxygen_bottom))
  )

occ_env_check


# ------------------------------------------------------------
# Check final occurrence dataset
# ------------------------------------------------------------

occ_btrw %>%
  summarise(
    total_events = n(),
    positive_events = sum(presence == 1),
    negative_events = sum(presence == 0),
    occurrence_percent = mean(presence == 1) * 100,
    years = n_distinct(Year),
    sites = n_distinct(Site)
  )


# ------------------------------------------------------------
# Annual and monthly occurrence patterns
# ------------------------------------------------------------

occ_btrw %>%
  group_by(Year, Month_f) %>%
  summarise(
    sampling_events = n(),
    positive_events = sum(presence == 1),
    occurrence_percent = mean(presence == 1) * 100,
    .groups = "drop"
  ) %>%
  arrange(Year, Month_f)


# ------------------------------------------------------------
# Predictor correlations
# ------------------------------------------------------------

##Include spring temp, exclude NAO
occ_correlations <- cor(
  occ_btrw %>%
    select(
      year_decade,
      temp_z,
      salinity_z,
      depth_z,
      spring_temp_same_year_z
    ),
  use = "complete.obs"
)

##Model with NAO
occ_correlations_NAO <- cor(
  occ_btrw %>%
    select(
      year_decade,
      temp_z,
      salinity_z,
      depth_z,
      spring_temp_same_year_z,
      winter_NAO_z
    ),
  use = "complete.obs"
)


## Check bottom and surface DO coverage due (missing DO values) ##

fall_1990_2018 %>%
  filter(Gear == "BTRW") %>%
  summarise(
    total_events = n(),
    
    bottom_DO_available =
      sum(!is.na(dissolved_oxygen_bottom)),
    
    bottom_DO_missing =
      sum(is.na(dissolved_oxygen_bottom)),
    
    bottom_DO_percent =
      mean(!is.na(dissolved_oxygen_bottom)) * 100,
    
    surface_DO_available =
      sum(!is.na(dissolved_oxygen_surface)),
    
    surface_DO_percent =
      mean(!is.na(dissolved_oxygen_surface)) * 100
  )

## Check bottom DO coverage by year ##

do_coverage_year <- fall_1990_2018 %>%
  filter(Gear == "BTRW") %>%
  group_by(Year) %>%
  summarise(
    total_events = n(),
    
    bottom_DO_available =
      sum(!is.na(dissolved_oxygen_bottom)),
    
    bottom_DO_percent =
      mean(!is.na(dissolved_oxygen_bottom)) * 100,
    
    .groups = "drop"
  ) %>%
  arrange(Year)

do_coverage_year

## Compare dataset sizes with and without bottom DO ##

do_comparison <- occ_btrw %>%
  summarise(
    events_without_DO = n(),
    positives_without_DO = sum(presence == 1),
    
    events_with_DO =
      sum(!is.na(dissolved_oxygen_bottom)),
    
    positives_with_DO =
      sum(presence == 1 &
            !is.na(dissolved_oxygen_bottom)),
    
    percent_events_retained =
      mean(!is.na(dissolved_oxygen_bottom)) * 100
  )

do_comparison

## Bottom DO coverage by month ##

occ_btrw %>%
  group_by(Month_f) %>%
  summarise(
    total_events = n(),
    DO_available = sum(!is.na(dissolved_oxygen_bottom)),
    DO_percent = mean(!is.na(dissolved_oxygen_bottom)) * 100,
    .groups = "drop"
  )

occ_btrw_NAO <- occ_btrw %>%
  filter(!is.na(winter_NAO_z))




# ============================================================
# 9. OCCURRENCE GLMM CANDIDATE MODELS
# ============================================================

# ------------------------------------------------------------
# 9A. Create dataset with post-2008 factor
# ------------------------------------------------------------

occ_post2008_btrw <- occ_btrw %>%
  mutate(
    post_2008 = factor(
      Year >= 2008,
      levels = c(FALSE, TRUE),
      labels = c("Pre-2008", "2008 and later")
    )
  )

# ------------------------------------------------------------
# 9B. GLMM control settings
# ------------------------------------------------------------

glmer_control <- glmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 2e5)
)

# ------------------------------------------------------------
# 9C. Baseline temporal models
# ------------------------------------------------------------

# Continuous year
occ_base_btrw <- glmer(
  presence ~ year_decade + (1 | Site_f),
  family = binomial,
  data = occ_btrw,
  control = glmer_control
)

# Continuous year + month
occ_base_month_btrw <- glmer(
  presence ~ year_decade + Month_f + (1 | Site_f),
  family = binomial,
  data = occ_btrw,
  control = glmer_control
)

# ------------------------------------------------------------
# 9D. Continuous-year environmental models
# ------------------------------------------------------------

# Temperature
occ_temp_btrw <- glmer(
  presence ~ year_decade + Month_f + temp_z +
    (1 | Site_f),
  family = binomial,
  data = occ_btrw,
  control = glmer_control
)

# Salinity
occ_salinity_btrw <- glmer(
  presence ~ year_decade + Month_f + salinity_z +
    (1 | Site_f),
  family = binomial,
  data = occ_btrw,
  control = glmer_control
)

# Depth
occ_depth_btrw <- glmer(
  presence ~ year_decade + Month_f + depth_z +
    (1 | Site_f),
  family = binomial,
  data = occ_btrw,
  control = glmer_control
)

# Same-year spring temperature
occ_spring_same_year_btrw <- glmer(
  presence ~ year_decade + Month_f +
    spring_temp_same_year_z + (1 | Site_f),
  family = binomial,
  data = occ_btrw,
  control = glmer_control
)

# Temperature + salinity
occ_temp_sal_btrw <- glmer(
  presence ~ year_decade + Month_f +
    temp_z + salinity_z + (1 | Site_f),
  family = binomial,
  data = occ_btrw,
  control = glmer_control
)

# Temperature + depth
occ_temp_depth_btrw <- glmer(
  presence ~ year_decade + Month_f +
    temp_z + depth_z + (1 | Site_f),
  family = binomial,
  data = occ_btrw,
  control = glmer_control
)

# Salinity + depth
occ_sal_depth_btrw <- glmer(
  presence ~ year_decade + Month_f +
    salinity_z + depth_z + (1 | Site_f),
  family = binomial,
  data = occ_btrw,
  control = glmer_control
)

# Temperature + same-year spring temperature
occ_temp_spring_same_year_btrw <- glmer(
  presence ~ year_decade + Month_f +
    temp_z + spring_temp_same_year_z +
    (1 | Site_f),
  family = binomial,
  data = occ_btrw,
  control = glmer_control
)

# Salinity + same-year spring temperature
occ_sal_spring_same_year_btrw <- glmer(
  presence ~ year_decade + Month_f +
    salinity_z + spring_temp_same_year_z +
    (1 | Site_f),
  family = binomial,
  data = occ_btrw,
  control = glmer_control
)

# Full continuous-year environmental model
occ_full_btrw <- glmer(
  presence ~ year_decade + Month_f +
    temp_z + salinity_z + depth_z +
    (1 | Site_f),
  family = binomial,
  data = occ_btrw,
  control = glmer_control
)

# ------------------------------------------------------------
# 9E. Post-2008 models
# ------------------------------------------------------------

# Post-2008 + month
occ_step_btrw <- glmer(
  presence ~ post_2008 + Month_f +
    (1 | Site_f),
  family = binomial,
  data = occ_post2008_btrw,
  control = glmer_control
)

# Post-2008 + temperature
occ_step_temp_btrw <- glmer(
  presence ~ post_2008 + Month_f + temp_z +
    (1 | Site_f),
  family = binomial,
  data = occ_post2008_btrw,
  control = glmer_control
)

# Post-2008 + salinity
occ_step_salinity_btrw <- glmer(
  presence ~ post_2008 + Month_f + salinity_z +
    (1 | Site_f),
  family = binomial,
  data = occ_post2008_btrw,
  control = glmer_control
)

# Post-2008 + depth
occ_step_depth_btrw <- glmer(
  presence ~ post_2008 + Month_f + depth_z +
    (1 | Site_f),
  family = binomial,
  data = occ_post2008_btrw,
  control = glmer_control
)

# Post-2008 + same-year spring temperature
occ_step_spring_same_year_btrw <- glmer(
  presence ~ post_2008 + Month_f +
    spring_temp_same_year_z +
    (1 | Site_f),
  family = binomial,
  data = occ_post2008_btrw,
  control = glmer_control
)

# Post-2008 + temperature + salinity
occ_step_temp_sal_btrw <- glmer(
  presence ~ post_2008 + Month_f +
    temp_z + salinity_z +
    (1 | Site_f),
  family = binomial,
  data = occ_post2008_btrw,
  control = glmer_control
)

# Post-2008 + salinity + same-year spring temperature
occ_step_sal_spring_same_year_btrw <- glmer(
  presence ~ post_2008 + Month_f +
    salinity_z + spring_temp_same_year_z +
    (1 | Site_f),
  family = binomial,
  data = occ_post2008_btrw,
  control = glmer_control
)

# Post-2008 + continuous year + month
occ_both_btrw <- glmer(
  presence ~ year_decade + post_2008 + Month_f +
    (1 | Site_f),
  family = binomial,
  data = occ_post2008_btrw,
  control = glmer_control
)

# Full post-2008 environmental model
occ_full_step_btrw <- glmer(
  presence ~ post_2008 + Month_f +
    temp_z + salinity_z + depth_z +
    (1 | Site_f),
  family = binomial,
  data = occ_post2008_btrw,
  control = glmer_control
)

# Full post-2008 model with same-year spring temperature
occ_post2008_spring_same_year <- glmer(
  presence ~ post_2008 + Month_f +
    temp_z + salinity_z + depth_z +
    spring_temp_same_year_z +
    (1 | Site_f),
  family = binomial,
  data = occ_post2008_btrw,
  control = glmer_control
)

# ------------------------------------------------------------
# 9F. Environmental model without temporal predictor
# ------------------------------------------------------------

occ_environmental_btrw <- glmer(
  presence ~ Month_f +
    temp_z + salinity_z + depth_z +
    (1 | Site_f),
  family = binomial,
  data = occ_btrw,
  control = glmer_control
)

# ------------------------------------------------------------
# 9G. Compare models with and without NAO
# ------------------------------------------------------------

# Create NAO-complete dataset
occ_btrw_NAO <- occ_btrw %>%
  filter(!is.na(winter_NAO_z))

# Environmental model without NAO
occ_full_NAO_comparison <- glmer(
  presence ~ year_decade + Month_f +
    temp_z + salinity_z + depth_z +
    (1 | Site_f),
  family = binomial,
  data = occ_btrw_NAO,
  control = glmer_control
)

# Environmental model with NAO
occ_full_NAO_btrw <- glmer(
  presence ~ year_decade + Month_f +
    temp_z + salinity_z + depth_z +
    winter_NAO_z +
    (1 | Site_f),
  family = binomial,
  data = occ_btrw_NAO,
  control = glmer_control
)

# Compare NAO models using identical observations
NAO_model_comparison <- MuMIn::AICc(
  occ_full_NAO_comparison,
  occ_full_NAO_btrw
)

print(NAO_model_comparison)

# ------------------------------------------------------------
# 9H. Compare environmental models with and without period
# ------------------------------------------------------------

occurrence_period_comparison <- MuMIn::AICc(
  occ_environmental_btrw,
  occ_full_step_btrw
)

print(occurrence_period_comparison)

# ------------------------------------------------------------
# 9I. Create candidate occurrence model list
# ------------------------------------------------------------

occurrence_models_btrw <- list(
  
  # Continuous-year models
  base = occ_base_month_btrw,
  temperature = occ_temp_btrw,
  salinity = occ_salinity_btrw,
  depth = occ_depth_btrw,
  spring_same_year = occ_spring_same_year_btrw,
  
  temp_salinity = occ_temp_sal_btrw,
  temp_depth = occ_temp_depth_btrw,
  salinity_depth = occ_sal_depth_btrw,
  
  temp_spring_same_year =
    occ_temp_spring_same_year_btrw,
  
  salinity_spring_same_year =
    occ_sal_spring_same_year_btrw,
  
  full = occ_full_btrw,
  
  # Environmental model without temporal predictor
  environmental = occ_environmental_btrw,
  
  # Post-2008 models
  post_2008 = occ_step_btrw,
  post_2008_temperature = occ_step_temp_btrw,
  post_2008_salinity = occ_step_salinity_btrw,
  post_2008_depth = occ_step_depth_btrw,
  
  post_2008_spring_same_year =
    occ_step_spring_same_year_btrw,
  
  post_2008_temp_salinity =
    occ_step_temp_sal_btrw,
  
  post_2008_salinity_spring_same_year =
    occ_step_sal_spring_same_year_btrw,
  
  year_plus_post_2008 = occ_both_btrw,
  
  full_post_2008 = occ_full_step_btrw,
  
  full_post_2008_spring_same_year =
    occ_post2008_spring_same_year
)

# ------------------------------------------------------------
# 9J. Check candidate model sample sizes
# ------------------------------------------------------------

occurrence_sample_sizes <- sapply(
  occurrence_models_btrw,
  nobs
)

print(occurrence_sample_sizes)

# Verify all models have the same number of observations
stopifnot(
  length(unique(occurrence_sample_sizes)) == 1
)

# Verify models use the same observation rows
occurrence_model_rows <- lapply(
  occurrence_models_btrw,
  function(model) rownames(model.frame(model))
)

stopifnot(
  all(vapply(
    occurrence_model_rows,
    function(x) identical(x, occurrence_model_rows[[1]]),
    logical(1)
  ))
)

# ------------------------------------------------------------
# 9K. Check convergence
# ------------------------------------------------------------

occurrence_convergence_btrw <- lapply(
  occurrence_models_btrw,
  function(model) model@optinfo$conv$lme4$messages
)

print(occurrence_convergence_btrw)

# ------------------------------------------------------------
# 9L. Check singularity
# ------------------------------------------------------------

occurrence_singularity_btrw <- sapply(
  occurrence_models_btrw,
  isSingular,
  tol = 1e-4
)

print(occurrence_singularity_btrw)

# ------------------------------------------------------------
# 9M. Compare candidate models using AICc
# ------------------------------------------------------------

occurrence_model_selection_btrw <- tibble(
  model = names(occurrence_models_btrw),
  AICc = as.numeric(sapply(
    occurrence_models_btrw,
    MuMIn::AICc
  ))
) %>%
  arrange(AICc) %>%
  mutate(
    delta_AICc = AICc - min(AICc),
    weight = exp(-0.5 * delta_AICc) /
      sum(exp(-0.5 * delta_AICc))
  )

# Display model selection results
occurrence_model_selection_btrw %>%
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 3)
    )
  ) %>%
  print(n = Inf)



# ============================================================
# 10. FINAL OCCURRENCE MODEL
# ============================================================

# ------------------------------------------------------------
# 10A. Select lowest-AICc model
# ------------------------------------------------------------

best_occurrence_model_name <-
  occurrence_model_selection_btrw$model[1]

final_occurrence_model <- occurrence_models_btrw[[
  best_occurrence_model_name
]]

# Display selected model name
print(best_occurrence_model_name)

# ------------------------------------------------------------
# 10B. View model summary
# ------------------------------------------------------------

summary(final_occurrence_model)

# Verify model formula
formula(final_occurrence_model)

# Verify AICc
MuMIn::AICc(final_occurrence_model)

# ------------------------------------------------------------
# 10C. Check observations and singularity
# ------------------------------------------------------------

nobs(final_occurrence_model)

isSingular(
  final_occurrence_model,
  tol = 1e-4
)

# Check convergence
final_occurrence_model@optinfo$conv$lme4$messages

# ------------------------------------------------------------
# 10D. Calculate coefficients, odds ratios,
#      confidence intervals, and p-values
# ------------------------------------------------------------

occurrence_results <- data.frame(
  term = names(lme4::fixef(final_occurrence_model)),
  beta = as.numeric(lme4::fixef(final_occurrence_model)),
  SE = sqrt(diag(vcov(final_occurrence_model)))
)

# Calculate 95% Wald confidence intervals
occurrence_CI <- confint(
  final_occurrence_model,
  parm = "beta_",
  method = "Wald"
)

# Calculate odds ratios and p-values
occurrence_results <- occurrence_results %>%
  mutate(
    OR = exp(beta),
    CI_low = exp(occurrence_CI[term, 1]),
    CI_high = exp(occurrence_CI[term, 2]),
    p_value = 2 * pnorm(
      abs(beta / SE),
      lower.tail = FALSE
    )
  )

# Display complete coefficient table
print(occurrence_results)

# ------------------------------------------------------------
# 10E. Create formatted coefficient table
# ------------------------------------------------------------

occurrence_coefficient_table_final <- occurrence_results %>%
  mutate(
    Predictor = dplyr::recode(
      term,
      "(Intercept)" = "Intercept",
      "year_decade" = "Year (per decade)",
      "post_20082008 and later" = "2008 and later",
      "Month_fOctober" = "October",
      "temp_z" = "Surface temperature (1 SD)",
      "salinity_z" = "Surface salinity (1 SD)",
      "depth_z" = "Mean depth (1 SD)",
      "spring_temp_same_year_z" =
        "Same-year spring temperature (1 SD)"
    ),
    Beta = round(beta, 3),
    SE = round(SE, 3),
    `Odds ratio` = round(OR, 3),
    `95% CI lower` = round(CI_low, 3),
    `95% CI upper` = round(CI_high, 3),
    `p-value` = ifelse(
      p_value < 0.001,
      "<0.001",
      sprintf("%.3f", p_value)
    )
  ) %>%
  dplyr::select(
    Predictor,
    Beta,
    SE,
    `Odds ratio`,
    `95% CI lower`,
    `95% CI upper`,
    `p-value`
  )

print(occurrence_coefficient_table_final)


# ============================================================
# 11. FINAL OCCURRENCE MODEL DIAGNOSTICS
# ============================================================

# ------------------------------------------------------------
# Simulate residuals
# ------------------------------------------------------------

set.seed(123)

sim_occurrence <- DHARMa::simulateResiduals(
  fittedModel = final_occurrence_model,
  n = 1000
)


# ------------------------------------------------------------
# Plot residual diagnostics
# ------------------------------------------------------------

plot(sim_occurrence)


# ------------------------------------------------------------
# Test residual uniformity
# ------------------------------------------------------------

DHARMa::testUniformity(sim_occurrence)


# ------------------------------------------------------------
# Test dispersion
# ------------------------------------------------------------

DHARMa::testDispersion(sim_occurrence)


# ------------------------------------------------------------
# Test for outliers
# ------------------------------------------------------------

DHARMa::testOutliers(sim_occurrence)


# ------------------------------------------------------------
# Check residual patterns
# ------------------------------------------------------------

DHARMa::plotResiduals(sim_occurrence)


# ------------------------------------------------------------
# Check collinearity
# ------------------------------------------------------------

performance::check_collinearity(
  final_occurrence_model
)


# ------------------------------------------------------------
# Check singularity
# ------------------------------------------------------------

isSingular(
  final_occurrence_model,
  tol = 1e-4
)


# ============================================================
# 12. FIGURE: PREDICTED OCCURRENCE BY PERIOD AND MONTH
# ============================================================

# Generate predicted occurrence probabilities
occ_predictions <- ggeffects::ggpredict(
  final_occurrence_model,
  terms = c("post_2008", "Month_f"),
  condition = c(
    salinity_z = 0,
    spring_temp_same_year_z = 0
  )
)

# View predictions
print(occ_predictions)

# Convert to dataframe
occ_predictions_df <- as.data.frame(occ_predictions)

# Plot predicted occurrence
fig_occurrence <- ggplot(
  occ_predictions_df,
  aes(
    x = x,
    y = predicted,
    fill = group
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_errorbar(
    aes(
      ymin = conf.low,
      ymax = conf.high
    ),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x = "Sampling period",
    y = "Predicted probability of white shrimp occurrence",
    fill = "Month"
  ) +
  theme_classic(base_size = 12)

print(fig_occurrence)


# ============================================================
# 13. FIGURE: PREDICTED OCCURRENCE ACROSS SALINITY
# ============================================================

# Generate predictions across standardized salinity
salinity_predictions <- ggeffects::ggpredict(
  final_occurrence_model,
  terms = "salinity_z [all]",
  condition = c(
    post_2008 = "2008 and later",
    Month_f = "October",
    spring_temp_same_year_z = 0
  )
)

salinity_predictions_df <- as.data.frame(
  salinity_predictions
)

# Convert standardized salinity to PSU
salinity_mean <- mean(
  occ_post2008_btrw$salinity_surface,
  na.rm = TRUE
)

salinity_sd <- sd(
  occ_post2008_btrw$salinity_surface,
  na.rm = TRUE
)

salinity_predictions_df <- salinity_predictions_df %>%
  mutate(
    salinity_PSU = x * salinity_sd + salinity_mean
  )

# Plot
fig_salinity <- ggplot(
  salinity_predictions_df,
  aes(
    x = salinity_PSU,
    y = predicted
  )
) +
  geom_ribbon(
    aes(
      ymin = conf.low,
      ymax = conf.high
    ),
    alpha = 0.2
  ) +
  geom_line(linewidth = 1) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x = "Surface salinity (PSU)",
    y = "Predicted probability of white shrimp occurrence"
  ) +
  theme_classic(base_size = 12)

print(fig_salinity) 

ggplot(
  occ_post2008_btrw,
  aes(x = salinity_surface)
) +
  geom_histogram(
    bins = 25,
    fill = "gray70",
    color = "black"
  ) +
  labs(
    x = "Surface salinity (PSU)",
    y = "Number of sampling events"
  ) +
  theme_classic()

salinity_occurrence_summary <- occ_post2008_btrw %>%
  mutate(
    salinity_group = cut(
      salinity_surface,
      breaks = c(0, 15, 20, 25, 30, 35, Inf),
      labels = c(
        "<15",
        "15–20",
        "20–25",
        "25–30",
        "30–35",
        ">35"
      ),
      right = FALSE
    )
  ) %>%
  group_by(salinity_group) %>%
  summarise(
    total_events = n(),
    positive_events = sum(presence == 1),
    occurrence_percent = round(
      100 * mean(presence == 1),
      2
    ),
    .groups = "drop"
  )

print(salinity_occurrence_summary)


# ============================================================
# 14. POSITIVE-CATCH ANALYSIS (ABUNDANCE)
# ============================================================

# ------------------------------------------------------------
# 14A. Prepare positive BTRW sampling events
# ------------------------------------------------------------

btrw_positive <- fall_1990_2018 %>%
  filter(
    Gear == "BTRW",
    abundance > 0,
    !is.na(Year),
    !is.na(Month),
    !is.na(Site),
    !is.na(area_covered),
    area_covered > 0
  )


# ------------------------------------------------------------
# Summarize positive sampling events
# ------------------------------------------------------------

positive_catch_check <- btrw_positive %>%
  summarise(
    n_positive_events = n(),
    n_years = n_distinct(Year),
    first_year = min(Year),
    last_year = max(Year),
    n_sites = n_distinct(Site),
    total_abundance = sum(abundance),
    min_area = min(area_covered),
    max_area = max(area_covered)
  )

positive_catch_check

# ------------------------------------------------------------
# Annual abundnace by sampling month
# Includes zero catches
# ------------------------------------------------------------

monthly_annual_abundance <- fall_1990_2018 %>%
  filter(
    Gear == "BTRW",
    Month %in% c(9, 10),
    !is.na(Year),
    !is.na(abundance),
    !is.na(area_covered),
    area_covered > 0
  ) %>%
  group_by(Year, Month) %>%
  summarise(
    sampling_events = n(),
    positive_events = sum(abundance > 0),
    total_shrimp = sum(abundance),
    total_area_sampled = sum(area_covered),
    catch_rate = (total_shrimp / total_area_sampled) * 1000,
    percent_positive = (positive_events / sampling_events) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    Month = factor(
      Month,
      levels = c(9, 10),
      labels = c("September", "October")
    )
  )


# View annual abundance summary
print(monthly_annual_abundance, n = Inf)


# ------------------------------------------------------------
# Plot annual abundance by sampling month
# ------------------------------------------------------------

fig_monthly_annual_abundance <- ggplot(
  monthly_annual_abundance,
  aes(
    x = Year,
    y = catch_rate,
    color = Month,
    group = Month
  )
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_color_manual(
    values = c(
      "September" = "#2878A5",
      "October" = "#D9822B"
    )
  ) +
  scale_x_continuous(
    breaks = seq(1990, 2018, by = 2)
  ) +
  labs(
    x = "Year",
    y = "White shrimp per 1,000 m² sampled",
    color = "Sampling month",
    title = "Annual White Shrimp Abundance by Month",
    subtitle = "September and October, Maryland Coastal Bays"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top"
  )

fig_monthly_annual_abundance


# ------------------------------------------------------------
# 14B. Annual abundnace by sampling month
# ------------------------------------------------------------

summary(btrw_positive$abundance)

hist(
  btrw_positive$abundance,
  main = "Distribution of positive white shrimp catches",
  xlab = "Number of white shrimp collected",
  breaks = 20
)

# Check whether abundance values are integer counts
any(
  btrw_positive$abundance !=
    floor(btrw_positive$abundance)
)

# Descriptive statistics
btrw_positive %>%
  summarise(
    mean_abundance = mean(abundance),
    median_abundance = median(abundance),
    sd_abundance = sd(abundance),
    min_abundance = min(abundance),
    max_abundance = max(abundance)
  )


# ------------------------------------------------------------
# 14C. Examine posive-catch distribution
# ------------------------------------------------------------

# Summary statistics
summary(btrw_positive$abundance)


# Histogram of positive catches
hist(
  btrw_positive$abundance,
  main = "Distribution of Positive White Shrimp Catches",
  xlab = "Number of white shrimp collected",
  ylab = "Frequency",
  breaks = 20
)


# Check whether abundance values are integer counts
any(
  btrw_positive$abundance != floor(btrw_positive$abundance)
)


# Descriptive statistics
positive_catch_summary <- btrw_positive %>%
  summarise(
    n_events = n(),
    mean_abundance = mean(abundance),
    median_abundance = median(abundance),
    sd_abundance = sd(abundance),
    min_abundance = min(abundance),
    max_abundance = max(abundance)
  )

positive_catch_summary


# ------------------------------------------------------------
# 14D. Prepare environmental monitoring dataset
# ------------------------------------------------------------

positive_env_btrw <- btrw_positive %>%
  filter(
    !is.na(temp_surface),
    !is.na(salinity_surface),
    !is.na(mean_depth)
  ) %>%
  mutate(
    Month_f = factor(
      Month,
      levels = c(9, 10),
      labels = c("September", "October")
    ),
    Site_f = factor(Site),
    
    # Time variable centered around 2008
    year_decade = (Year - 2008) / 10,
    
    # Sampling area offset
    log_area = log(area_covered),
    
    # Standardized environmental predictors
    temp_z = as.numeric(scale(temp_surface)),
    salinity_z = as.numeric(scale(salinity_surface)),
    depth_z = as.numeric(scale(mean_depth))
  )


# ------------------------------------------------------------
# Summarize environmental modeling dataset
# ------------------------------------------------------------

positive_env_check <- positive_env_btrw %>%
  summarise(
    n_events = n(),
    n_sites = n_distinct(Site),
    n_years = n_distinct(Year),
    first_year = min(Year),
    last_year = max(Year),
    total_abundance = sum(abundance)
  )

positive_env_check


# ------------------------------------------------------------
# Additional data checks
# ------------------------------------------------------------

# Abundance distribution
summary(positive_env_btrw$abundance)

# Check integer counts
any(
  positive_env_btrw$abundance !=
    floor(positive_env_btrw$abundance)
)

# Sampling area
summary(positive_env_btrw$area_covered)

# Sampling month
table(positive_env_btrw$Month_f)

# Number of events per site
table(positive_env_btrw$Site_f)


# ------------------------------------------------------------
# 14E. Fit candidate posiive-catch abundnace models
# ------------------------------------------------------------

# Baseline temporal model
pos_base <- glmmTMB::glmmTMB(
  abundance ~ year_decade + Month_f +
    offset(log_area) +
    (1 | Site_f),
  family = glmmTMB::truncated_nbinom2(),
  data = positive_env_btrw
)


# Temperature model
pos_temp <- glmmTMB::glmmTMB(
  abundance ~ year_decade + Month_f + temp_z +
    offset(log_area) +
    (1 | Site_f),
  family = glmmTMB::truncated_nbinom2(),
  data = positive_env_btrw
)


# Salinity model
pos_salinity <- glmmTMB::glmmTMB(
  abundance ~ year_decade + Month_f + salinity_z +
    offset(log_area) +
    (1 | Site_f),
  family = glmmTMB::truncated_nbinom2(),
  data = positive_env_btrw
)


# Depth model
pos_depth <- glmmTMB::glmmTMB(
  abundance ~ year_decade + Month_f + depth_z +
    offset(log_area) +
    (1 | Site_f),
  family = glmmTMB::truncated_nbinom2(),
  data = positive_env_btrw
)


# Salinity and depth model
pos_sal_depth <- glmmTMB::glmmTMB(
  abundance ~ year_decade + Month_f +
    salinity_z + depth_z +
    offset(log_area) +
    (1 | Site_f),
  family = glmmTMB::truncated_nbinom2(),
  data = positive_env_btrw
)


# Full environmental model
pos_full <- glmmTMB::glmmTMB(
  abundance ~ year_decade + Month_f +
    temp_z + salinity_z + depth_z +
    offset(log_area) +
    (1 | Site_f),
  family = glmmTMB::truncated_nbinom2(),
  data = positive_env_btrw
)



# ------------------------------------------------------------
# 14F. Check model convergence
# ------------------------------------------------------------

positive_models_btrw <- list(
  base = pos_base,
  temperature = pos_temp,
  salinity = pos_salinity,
  depth = pos_depth,
  salinity_depth = pos_sal_depth,
  full = pos_full
)


# Check convergence messages
positive_convergence_btrw <- lapply(
  positive_models_btrw,
  function(model) model$fit$convergence
)

positive_convergence_btrw


# Check model warnings
positive_warnings_btrw <- lapply(
  positive_models_btrw,
  function(model) model$fit$message
)

positive_warnings_btrw


# Check observation counts
sapply(
  positive_models_btrw,
  nobs
) 


# ------------------------------------------------------------
# 14G. Compare candidate models using AICc
# ------------------------------------------------------------

positive_model_table <- tibble(
  model = names(positive_models_btrw),
  
  AICc = sapply(
    positive_models_btrw,
    MuMIn::AICc
  )
) %>%
  arrange(AICc) %>%
  mutate(
    delta_AICc = AICc - min(AICc),
    
    weight = exp(-0.5 * delta_AICc) /
      sum(exp(-0.5 * delta_AICc))
  )


# Display model selection results
positive_model_table %>%
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 3)
    )
  )



# ------------------------------------------------------------
# 14H. Select and examine final positive-catch model
# ------------------------------------------------------------

final_positive_model <- positive_models_btrw[[
  positive_model_table$model[1]
]]


# View final model summary
summary(final_positive_model)


# Check number of observations
nobs(final_positive_model)


# ------------------------------------------------------------
# 14I. Calculate coefficients, rate ratios, and confidence intervals
# ------------------------------------------------------------

positive_coefficients <- summary(
  final_positive_model
)$coefficients$cond


positive_coefficient_table <- data.frame(
  term = rownames(positive_coefficients),
  beta = positive_coefficients[, "Estimate"],
  SE = positive_coefficients[, "Std. Error"],
  z_value = positive_coefficients[, "z value"],
  p_value = positive_coefficients[, "Pr(>|z|)"]
)


# Calculate 95% confidence intervals
positive_CI <- confint(
  final_positive_model,
  parm = "beta_",
  component = "cond",
  method = "wald"
)


# Calculate rate ratios and confidence intervals
positive_coefficient_table <- positive_coefficient_table %>%
  mutate(
    rate_ratio = exp(beta),
    CI_low = exp(positive_CI[term, 1]),
    CI_high = exp(positive_CI[term, 2])
  )


# View results
positive_coefficient_table


# ============================================================
# 14J. FINAL POSITIVE-CATCH MODEL DIAGNOSTICS
# ============================================================

# ------------------------------------------------------------
# Simulate residuals
# ------------------------------------------------------------

set.seed(123)

sim_positive <- DHARMa::simulateResiduals(
  fittedModel = final_positive_model,
  n = 1000
)


# ------------------------------------------------------------
# Plot residual diagnostics
# ------------------------------------------------------------

plot(sim_positive)


# ------------------------------------------------------------
# Test residual uniformity
# ------------------------------------------------------------

DHARMa::testUniformity(sim_positive)


# ------------------------------------------------------------
# Test dispersion
# ------------------------------------------------------------

DHARMa::testDispersion(sim_positive)


# ------------------------------------------------------------
# Test for outliers
# ------------------------------------------------------------

DHARMa::testOutliers(sim_positive)


# ------------------------------------------------------------
# Check residual patterns
# ------------------------------------------------------------

DHARMa::plotResiduals(sim_positive)


# ------------------------------------------------------------
# Check model convergence
# ------------------------------------------------------------

final_positive_model$fit$convergence

final_positive_model$fit$message



# Check the salinity model
summary(pos_salinity)

# Check the salinity + depth model
summary(pos_sal_depth)

# Check Hessian status
pos_salinity$sdr$pdHess
pos_sal_depth$sdr$pdHess



# ============================================================
# REFIT SALINITY MODEL WITH TIGHTER CONVERGENCE TOLERANCES
# ============================================================

pos_salinity_refit <- update(
  pos_salinity,
  control = glmmTMBControl(
    optimizer = nlminb,
    optCtrl = list(
      iter.max = 10000,
      eval.max = 10000,
      rel.tol = 1e-10,
      x.tol = 1e-10
    )
  )
)

# Check convergence
cat("Convergence code:\n")
print(pos_salinity_refit$fit$convergence)

cat("Convergence message:\n")
print(pos_salinity_refit$fit$message)

cat("Positive-definite Hessian:\n")
print(pos_salinity_refit$sdr$pdHess)

# Compare coefficients
coefficient_comparison <- cbind(
  Original = fixef(pos_salinity)$cond,
  Refit = fixef(pos_salinity_refit)$cond
)

print(coefficient_comparison)

# Compare log-likelihoods
cat("Original log-likelihood:\n")
print(logLik(pos_salinity))

cat("Refit log-likelihood:\n")
print(logLik(pos_salinity_refit))

# Compare AICc
cat("Original AICc:\n")
print(MuMIn::AICc(pos_salinity))

cat("Refit AICc:\n")
print(MuMIn::AICc(pos_salinity_refit))


# ============================================================
# 15. FIGURE: ANNUAL WHITE SHRIMP CATCH RATE
# ============================================================


# Keep all BTRW sampling events, including zero catches
btrw_all_events <- fall_1990_2018 %>%
  filter(
    Gear == "BTRW",
    Month %in% c(9, 10),
    !is.na(Year),
    !is.na(abundance),
    !is.na(area_covered),
    area_covered > 0
  ) %>%
  mutate(
    Month_f = factor(
      Month,
      levels = c(9, 10),
      labels = c("September", "October")
    )
  )

# Calculate annual catch rate by month
annual_catch_rate <- btrw_all_events %>%
  group_by(Year, Month_f) %>%
  summarise(
    total_shrimp = sum(abundance),
    total_area = sum(area_covered),
    n_events = n(),
    catch_rate = (total_shrimp / total_area) * 1000,
    .groups = "drop"
  )

# Review results
print(annual_catch_rate) 

fig_annual_catch_rate <- ggplot(
  annual_catch_rate,
  aes(
    x = Year,
    y = catch_rate,
    color = Month_f,
    group = Month_f
  )
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_color_manual(
    values = c(
      "September" = "#2878A5",
      "October" = "#D17A35"
    )
  ) +
  scale_x_continuous(
    breaks = seq(1990, 2018, by = 2)
  ) +
  labs(
    title = "Annual White Shrimp Catch Rate",
    x = "Year",
    y = expression("Catch rate (shrimp per 1,000 m"^2*")"),
    color = "Sampling month"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

fig_annual_catch_rate



# ============================================================
# 16. FIGURE: PREDICTED POSITIVE-CATCH ABUNDANCE ACROSS SALINITY
# ============================================================

# Median sampled area
median_area <- median(
  positive_env_btrw$area_covered,
  na.rm = TRUE
)

# Salinity mean and standard deviation
salinity_mean <- mean(
  positive_env_btrw$salinity_surface,
  na.rm = TRUE
)

salinity_sd <- sd(
  positive_env_btrw$salinity_surface,
  na.rm = TRUE
)

# Create sequence of observed salinity values
salinity_sequence <- seq(
  min(positive_env_btrw$salinity_surface, na.rm = TRUE),
  max(positive_env_btrw$salinity_surface, na.rm = TRUE),
  length.out = 100
)

# Create prediction dataset
pred_positive_direct <- positive_env_btrw[
  rep(1, length(salinity_sequence)), ,
  drop = FALSE
]

# Set predictor values
pred_positive_direct$salinity_surface <- salinity_sequence

pred_positive_direct$salinity_z <- (
  salinity_sequence - salinity_mean
) / salinity_sd

pred_positive_direct$year_decade <- 0

pred_positive_direct$Month_f <- factor(
  "October",
  levels = levels(positive_env_btrw$Month_f)
)

pred_positive_direct$log_area <- log(median_area)

# Obtain predictions and standard errors on link scale
pred_link <- predict(
  final_positive_model,
  newdata = pred_positive_direct,
  type = "link",
  se.fit = TRUE,
  re.form = NA
)

# Calculate 95% confidence limits on link scale
link_low <- pred_link$fit - 1.96 * pred_link$se.fit
link_high <- pred_link$fit + 1.96 * pred_link$se.fit

# Convert to underlying negative-binomial count parameters
mu <- exp(pred_link$fit)
mu_low <- exp(link_low)
mu_high <- exp(link_high)

# Negative-binomial dispersion parameter
theta <- sigma(final_positive_model)

# Convert underlying count parameters to
# expected abundance conditional on a positive catch
truncated_mean <- function(mu, theta) {
  p0 <- (theta / (theta + mu))^theta
  mu / (1 - p0)
}

pred_positive_direct$predicted <- truncated_mean(mu, theta)

pred_positive_direct$conf.low <- truncated_mean(
  mu_low, theta
)

pred_positive_direct$conf.high <- truncated_mean(
  mu_high, theta
)

# Check predictions
head(
  pred_positive_direct[
    , c(
      "salinity_surface",
      "predicted",
      "conf.low",
      "conf.high"
    )
  ]
) 



fig_positive_salinity <- ggplot() +
    geom_ribbon(
    data = pred_positive_direct,
    aes(
      x = salinity_surface,
      ymin = conf.low,
      ymax = conf.high
    ),
    fill = "grey70",
    alpha = 0.4
  ) +
    geom_line(
    data = pred_positive_direct,
    aes(
      x = salinity_surface,
      y = predicted
    ),
    color = "black",
    linewidth = 1
  ) +
    geom_point(
    data = positive_env_btrw,
    aes(
      x = salinity_surface,
      y = abundance
    ),
    color = "#2878A5",
    alpha = 0.55,
    size = 2
  ) +
  # Square-root transformed y-axis
  scale_y_continuous(
    trans = "sqrt",
    breaks = c(0, 1, 10, 25, 50, 100, 200, 350)
  ) +
  
  labs(
    x = "Surface salinity (PSU)",
    y = "Positive-catch abundance (shrimp)"
  ) +
  
  theme_classic(base_size = 12)

# Display figure
print(fig_positive_salinity)


# Summary of observed salinity
summary(positive_env_btrw$salinity_surface)

# Number of observations across salinity values
hist(
  positive_env_btrw$salinity_surface,
  breaks = 10,
  main = "Distribution of Surface Salinity",
  xlab = "Surface salinity (PSU)"
)


positive_env_btrw %>%
  arrange(desc(abundance)) %>%
  select(
    Year,
    Month,
    Site,
    abundance,
    salinity_surface,
    temp_surface,
    mean_depth,
    area_covered
  ) %>%
  slice_head(n = 10)

positive_env_btrw %>%
  arrange(desc(abundance)) %>%
  select(
    sample_id,
    SampleDate,
    Year,
    Month,
    Site,
    abundance,
    salinity_surface,
    temp_surface,
    mean_depth,
    area_covered
  ) %>%
  slice_head(n = 10) %>%
  print(width = Inf)


# Temporarily exclude the three largest catches
positive_sensitivity <- positive_env_btrw %>%
  arrange(desc(abundance)) %>%
  slice(-(1:3))

# Refit the same model
pos_salinity_sensitivity <- glmmTMB(
  abundance ~ year_decade + Month_f + salinity_z +
    offset(log_area) + (1 | Site_f),
  data = positive_sensitivity,
  family = truncated_nbinom2()
)

# Compare the original and sensitivity models
summary(pos_salinity_sensitivity)


# Identify the three largest observations
largest_catches <- positive_env_btrw %>%
  slice_max(abundance, n = 3, with_ties = FALSE)

# Function to refit model after removing one observation
run_sensitivity <- function(data, excluded_id) {
  
  sensitivity_data <- data %>%
    filter(sample_id != excluded_id)
  
  model <- glmmTMB(
    abundance ~ year_decade + Month_f + salinity_z +
      offset(log_area) + (1 | Site_f),
    data = sensitivity_data,
    family = truncated_nbinom2()
  )
  
  coefficients <- summary(model)$coefficients$cond
  
  tibble(
    excluded_sample = excluded_id,
    n = nobs(model),
    salinity_beta = coefficients["salinity_z", "Estimate"],
    salinity_SE = coefficients["salinity_z", "Std. Error"],
    salinity_p = coefficients["salinity_z", "Pr(>|z|)"],
    year_beta = coefficients["year_decade", "Estimate"],
    year_p = coefficients["year_decade", "Pr(>|z|)"]
  )
}

# Run one-at-a-time sensitivity analyses
individual_sensitivity <- map_dfr(
  largest_catches$sample_id,
  ~ run_sensitivity(positive_env_btrw, .x)
)

individual_sensitivity 

# ============================================================
# 17. WHITE SHRIMP LENGTH ANALYSIS
# Event-level median length, BTRW sampling
# ============================================================

# ------------------------------------------------------------
# 17A. Extract standardized individual length records
# ------------------------------------------------------------

# 1990–2012
length_MDNR1 <- MDNR1 %>%
  filter(
    grepl(
      "White Shrimp",
      `Speciment name`,
      ignore.case = TRUE
    ),
    !is.na(`Length (mm)`)
  ) %>%
  mutate(
    SampleDate_clean = as.Date(`Sample Date`),
    Year = lubridate::year(SampleDate_clean),
    Month = lubridate::month(SampleDate_clean),
    
    # Create sampling event ID
    sample_id = paste(
      `Sample Date`,
      `Sample Time`,
      `Site code`,
      `Gear Code`,
      sep = "_"
    )
  ) %>%
  filter(
    Month %in% c(9, 10),
    `Gear Code` %in% c("BTRW")
  ) %>%
  transmute(
    sample_id,
    Year,
    Month,
    SampleDate = SampleDate_clean,
    Site = as.character(`Site code`),
    Gear = `Gear Code`,
    length_mm = as.numeric(`Length (mm)`)
  )


# 2013–2016
length_MDNR2 <- MDNR2 %>%
  filter(
    Year <= 2016,
    grepl(
      "White Shrimp",
      SpecimenName,
      ignore.case = TRUE
    ),
    !is.na(SpecimenLength)
  ) %>%
  mutate(
    SampleDate_clean = as.Date(
      SampleDate,
      format = "%m/%d/%Y"
    ),
    Year = lubridate::year(SampleDate_clean),
    Month = lubridate::month(SampleDate_clean)
  ) %>%
  filter(
    Month %in% c(9, 10),
    GearCode %in% c("BTRW")
  ) %>%
  transmute(
    sample_id,
    Year,
    Month,
    SampleDate = SampleDate_clean,
    Site = as.character(SiteNumber),
    Gear = GearCode,
    length_mm = SpecimenLength
  )


# 2017–2018
length_MDNR3 <- MDNR3 %>%
  filter(
    Year %in% c(2017, 2018),
    grepl(
      "White Shrimp",
      SpecimenName,
      ignore.case = TRUE
    ),
    !is.na(SpecimenLength)
  ) %>%
  mutate(
    SampleDate_clean = as.Date(
      SampleDate,
      format = "%m/%d/%Y"
    ),
    Year = lubridate::year(SampleDate_clean),
    Month = lubridate::month(SampleDate_clean)
  ) %>%
  filter(
    Month %in% c(9, 10),
    GearCode %in% c("BTRW")
  ) %>%
  transmute(
    sample_id,
    Year,
    Month,
    SampleDate = SampleDate_clean,
    Site = as.character(SiteNumber),
    Gear = GearCode,
    length_mm = SpecimenLength
  )


# ------------------------------------------------------------
# 17B. Combine 1990–2018 length records
# ------------------------------------------------------------

length_1990_2018 <- bind_rows(
  length_MDNR1,
  length_MDNR2,
  length_MDNR3
) %>%
  arrange(
    Year,
    Month,
    SampleDate,
    Site,
    Gear
  )


# Check combined individual-level length dataset
length_check <- length_1990_2018 %>%
  summarise(
    n_lengths = n(),
    n_events = n_distinct(sample_id),
    n_sites = n_distinct(Site),
    n_years = n_distinct(Year),
    first_year = min(Year),
    last_year = max(Year),
    min_length = min(length_mm),
    median_length = median(length_mm),
    mean_length = mean(length_mm),
    max_length = max(length_mm)
  )

length_check


# ------------------------------------------------------------
# 17C. Summarize length by sampling event
# ------------------------------------------------------------

# Median length is used as the event-level response because
# the number of individual shrimp measured varies among events.

length_by_event <- length_1990_2018 %>%
  group_by(
    sample_id,
    Year,
    Month,
    Site,
    Gear
  ) %>%
  summarise(
    n_measured = n(),
    median_length = median(length_mm),
    mean_length = mean(length_mm),
    min_length = min(length_mm),
    max_length = max(length_mm),
    .groups = "drop"
  )


# Check number of shrimp measured per event
length_by_event %>%
  summarise(
    n_events = n(),
    min_measured = min(n_measured),
    median_measured = median(n_measured),
    mean_measured = mean(n_measured),
    max_measured = max(n_measured)
  )


# ------------------------------------------------------------
# FIGURE: EVENT-LEVEL MEDIAN SHRIMP LENGTH OVER TIME
# ------------------------------------------------------------

ggplot(
  length_by_event,
  aes(x = factor(Year), y = median_length)
) +
  geom_jitter(
    width = 0.15,
    size = 2.5,
    alpha = 0.7,
    color = "black"
  ) +
  labs(
    title = "Median shrimp length over time",
    x = "Year",
    y = "Median length (mm)"
  ) +
  theme_classic()


# ------------------------------------------------------------
# 17D. Join environmental variables
# ------------------------------------------------------------

length_event_env <- length_by_event %>%
  left_join(
    fall_1990_2018 %>%
      select(
        sample_id,
        temp_surface,
        salinity_surface,
        mean_depth
      ),
    by = "sample_id"
  )


# Check environmental-data availability
length_event_env %>%
  summarise(
    n_events = n(),
    temp_available = sum(!is.na(temp_surface)),
    salinity_available = sum(!is.na(salinity_surface)),
    depth_available = sum(!is.na(mean_depth))
  )


# ------------------------------------------------------------
# 17E. Prepare BTRW-only length modeling dataset
# ------------------------------------------------------------

length_btrw <- length_event_env %>%
  filter(
    Gear == "BTRW",
    !is.na(Year),
    !is.na(Month),
    !is.na(Site),
    !is.na(median_length)
  ) %>%
  mutate(
    year_decade = (Year - 2008) / 10,
    
    Month_f = factor(
      Month,
      levels = c(9, 10),
      labels = c("September", "October")
    ),
    
    Site_f = factor(Site)
  )


# ------------------------------------------------------------
# 17F. Examine BTRW length data and environmental coverage
# ------------------------------------------------------------

# Overall length summary
length_btrw %>%
  summarise(
    n_events = n(),
    n_shrimp_measured = sum(n_measured),
    n_years = n_distinct(Year),
    n_sites = n_distinct(Site),
    first_year = min(Year),
    last_year = max(Year),
    median_event_length = median(median_length),
    mean_event_length = mean(median_length),
    sd_event_length = sd(median_length),
    min_event_length = min(median_length),
    max_event_length = max(median_length)
  )


# Environmental data availability
length_btrw %>%
  summarise(
    n_events = n(),
    temp_available = sum(!is.na(temp_surface)),
    salinity_available = sum(!is.na(salinity_surface)),
    depth_available = sum(!is.na(mean_depth))
  )


# Prepare complete-case environmental modeling dataset
length_btrw_env <- length_btrw %>%
  filter(
    !is.na(median_length),
    !is.na(year_decade),
    !is.na(Month_f),
    !is.na(temp_surface),
    !is.na(salinity_surface),
    !is.na(mean_depth)
  ) %>%
  mutate(
    temp_z = as.numeric(scale(temp_surface)),
    salinity_z = as.numeric(scale(salinity_surface)),
    depth_z = as.numeric(scale(mean_depth))
  )


# Check final modeling sample size
length_btrw_env %>%
  summarise(
    n_events = n(),
    n_shrimp_measured = sum(n_measured),
    n_sites = n_distinct(Site),
    n_years = n_distinct(Year),
    first_year = min(Year),
    last_year = max(Year)
  )


# Compare September and October
length_btrw_env %>%
  group_by(Month_f) %>%
  summarise(
    n_events = n(),
    n_shrimp_measured = sum(n_measured),
    median_event_length = median(median_length),
    mean_event_length = mean(median_length),
    sd_event_length = sd(median_length),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 17G. Candidate Gaussian linear models
# ------------------------------------------------------------

# Preliminary site-random-intercept LMMs were evaluated.
# Leading models were singular, with site variance approaching
# zero, and corresponding ordinary linear models had lower AICc
# in direct comparisons. Final event-level analyses therefore
# use Gaussian linear models.



# Base model
length_lm_base <- lm(
  median_length ~
    year_decade +
    Month_f,
  data = length_btrw_env
)


# Temperature
length_lm_temp <- lm(
  median_length ~
    year_decade +
    Month_f +
    temp_z,
  data = length_btrw_env
)


# Salinity
length_lm_salinity <- lm(
  median_length ~
    year_decade +
    Month_f +
    salinity_z,
  data = length_btrw_env
)


# Depth
length_lm_depth <- lm(
  median_length ~
    year_decade +
    Month_f +
    depth_z,
  data = length_btrw_env
)


# Temperature + salinity
length_lm_temp_salinity <- lm(
  median_length ~
    year_decade +
    Month_f +
    temp_z +
    salinity_z,
  data = length_btrw_env
)


# Temperature + depth
length_lm_temp_depth <- lm(
  median_length ~
    year_decade +
    Month_f +
    temp_z +
    depth_z,
  data = length_btrw_env
)


# Salinity + depth
length_lm_salinity_depth <- lm(
  median_length ~
    year_decade +
    Month_f +
    salinity_z +
    depth_z,
  data = length_btrw_env
)


# Full environmental model
length_lm_full <- lm(
  median_length ~
    year_decade +
    Month_f +
    temp_z +
    salinity_z +
    depth_z,
  data = length_btrw_env
)

# ------------------------------------------------------------
# 17H. AICc model selection
# ------------------------------------------------------------

length_models <- list(
  base = length_lm_base,
  temperature = length_lm_temp,
  salinity = length_lm_salinity,
  depth = length_lm_depth,
  temp_salinity = length_lm_temp_salinity,
  temp_depth = length_lm_temp_depth,
  salinity_depth = length_lm_salinity_depth,
  full = length_lm_full
)


length_model_selection <- tibble(
  model = names(length_models),
  
  AICc = sapply(
    length_models,
    MuMIn::AICc
  )
) %>%
  arrange(AICc) %>%
  mutate(
    delta_AICc = AICc - min(AICc),
    
    weight = exp(-0.5 * delta_AICc) /
      sum(exp(-0.5 * delta_AICc))
  )


length_model_selection %>%
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 3)
    )
  )

# ------------------------------------------------------------
# 17I. Examine leading models
# ------------------------------------------------------------

summary(length_lm_depth)

summary(length_lm_temp_depth)

summary(length_lm_base)

summary(length_lm_salinity_depth)

confint(length_lm_depth)

confint(length_lm_temp_depth)


# ------------------------------------------------------------
# 17J. Diagnostics for lowest-AICc model
# ------------------------------------------------------------

# Diagnostic plots
par(mfrow = c(2, 2))
plot(length_lm_depth)
par(mfrow = c(1, 1))


# Normality of residuals
shapiro.test(residuals(length_lm_depth))


# Check variance inflation factors
car::vif(length_lm_depth)


# Influence diagnostics
cooks_d <- cooks.distance(length_lm_depth)


# Identify potentially influential observations
influential_obs <- which(
  cooks_d > 4 / nobs(length_lm_depth)
)


# Display flagged observations and Cook's distances
influential_obs
cooks_d[influential_obs]


# Identify five observations with the highest Cook's distances
top_influential <- head(
  order(cooks_d, decreasing = TRUE),
  5
)


# Display information for the five most influential observations
length_btrw_env %>%
  mutate(
    observation = row_number(),
    cooks_d = cooks_d
  ) %>%
  slice(top_influential) %>%
  select(
    observation,
    Year,
    Month,
    Site,
    median_length,
    n_measured,
    mean_depth,
    temp_surface,
    salinity_surface,
    cooks_d
  )


# Breusch-Pagan test for heteroscedasticity
lmtest::bptest(length_lm_depth)


# HC3 heteroscedasticity-robust standard errors
lmtest::coeftest(
  length_lm_depth,
  vcov. = sandwich::vcovHC(
    length_lm_depth,
    type = "HC3"
  )
)

# ------------------------------------------------------------
# 17K. Sensitivity analysis: candidate model selection
# Use only sampling events with at least 2 shrimp measured
# ------------------------------------------------------------

length_btrw_sensitivity <- length_btrw %>%
  filter(
    n_measured >= 2,
    !is.na(median_length),
    !is.na(year_decade),
    !is.na(Month_f),
    !is.na(temp_surface),
    !is.na(salinity_surface),
    !is.na(mean_depth)
  ) %>%
  mutate(
    temp_z = as.numeric(scale(temp_surface)),
    salinity_z = as.numeric(scale(salinity_surface)),
    depth_z = as.numeric(scale(mean_depth))
  )


# Check the sensitivity dataset
length_btrw_sensitivity %>%
  summarise(
    n_events = n(),
    n_shrimp_measured = sum(n_measured),
    n_sites = n_distinct(Site),
    n_years = n_distinct(Year),
    first_year = min(Year),
    last_year = max(Year)
  )


# Base model
length_lm_base_sensitivity <- lm(
  median_length ~ year_decade + Month_f,
  data = length_btrw_sensitivity
)


# Temperature model
length_lm_temp_sensitivity <- lm(
  median_length ~ year_decade + Month_f + temp_z,
  data = length_btrw_sensitivity
)


# Salinity model
length_lm_salinity_sensitivity <- lm(
  median_length ~ year_decade + Month_f + salinity_z,
  data = length_btrw_sensitivity
)


# Depth model
length_lm_depth_sensitivity <- lm(
  median_length ~ year_decade + Month_f + depth_z,
  data = length_btrw_sensitivity
)


# Temperature + salinity
length_lm_temp_salinity_sensitivity <- lm(
  median_length ~ year_decade + Month_f + temp_z + salinity_z,
  data = length_btrw_sensitivity
)


# Temperature + depth
length_lm_temp_depth_sensitivity <- lm(
  median_length ~ year_decade + Month_f + temp_z + depth_z,
  data = length_btrw_sensitivity
)


# Salinity + depth
length_lm_salinity_depth_sensitivity <- lm(
  median_length ~ year_decade + Month_f + salinity_z + depth_z,
  data = length_btrw_sensitivity
)


# Full model
length_lm_full_sensitivity <- lm(
  median_length ~ year_decade + Month_f +
    temp_z + salinity_z + depth_z,
  data = length_btrw_sensitivity
)


# Candidate model list
length_models_sensitivity <- list(
  base = length_lm_base_sensitivity,
  temperature = length_lm_temp_sensitivity,
  salinity = length_lm_salinity_sensitivity,
  depth = length_lm_depth_sensitivity,
  temp_salinity = length_lm_temp_salinity_sensitivity,
  temp_depth = length_lm_temp_depth_sensitivity,
  salinity_depth = length_lm_salinity_depth_sensitivity,
  full = length_lm_full_sensitivity
)


# AICc model selection
length_model_selection_sensitivity <- tibble(
  model = names(length_models_sensitivity),
  AICc = sapply(
    length_models_sensitivity,
    MuMIn::AICc
  )
) %>%
  arrange(AICc) %>%
  mutate(
    delta_AICc = AICc - min(AICc),
    weight = exp(-0.5 * delta_AICc) /
      sum(exp(-0.5 * delta_AICc))
  )


length_model_selection_sensitivity %>%
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 3)
    )
  )


# Compare original and sensitivity model selection
length_model_comparison <- length_model_selection %>%
  select(model, AICc, delta_AICc, weight) %>%
  rename(
    AICc_original = AICc,
    delta_original = delta_AICc,
    weight_original = weight
  ) %>%
  full_join(
    length_model_selection_sensitivity %>%
      select(model, AICc, delta_AICc, weight) %>%
      rename(
        AICc_sensitivity = AICc,
        delta_sensitivity = delta_AICc,
        weight_sensitivity = weight
      ),
    by = "model"
  )


length_model_comparison



# ------------------------------------------------------------
# 17L. Figure: Median length versus sampling depth
# ------------------------------------------------------------

# 1. Check that the depth model exists
stopifnot(exists("length_lm_depth"))
stopifnot(exists("length_btrw_env"))

# 2. Extract the exact observations used to fit the model
model_frame <- model.frame(length_lm_depth)

# Match the model observations to the original dataset
length_plot_data <- length_btrw_env[
  rownames(model_frame),
  ,
  drop = FALSE
]

# Check that the observations match
cat("Observations used in depth model:", nobs(length_lm_depth), "\n")
cat("Rows in plotting dataset:", nrow(length_plot_data), "\n")

if (nrow(length_plot_data) != nobs(length_lm_depth)) {
  stop("The plotting data do not match the observations used in the depth model.")
}

# 3. Check depth standardization
depth_mean <- mean(length_plot_data$mean_depth)
depth_sd <- sd(length_plot_data$mean_depth)

depth_z_check <- (length_plot_data$mean_depth - depth_mean) / depth_sd

cat("\nDepth standardization checks:\n")
cat("Mean of standardized depth:", mean(length_plot_data$depth_z), "\n")
cat("SD of standardized depth:", sd(length_plot_data$depth_z), "\n")

# Confirm that depth_z matches the raw depth values
if (!isTRUE(all.equal(
  as.numeric(length_plot_data$depth_z),
  as.numeric(depth_z_check),
  tolerance = 1e-6
))) {
  stop("depth_z does not match standardized mean_depth in the model dataset.")
}

# 4. Generate predictions from the depth model
depth_predictions <- ggpredict(
  length_lm_depth,
  terms = "depth_z [all]"
) %>%
  as.data.frame() %>%
  mutate(
    depth = x * depth_sd + depth_mean
  )

# 5. Check prediction values
cat("\nPrediction range (meters):\n")
print(range(depth_predictions$depth, na.rm = TRUE))

# 6. Create the figure
fig_length_depth <- ggplot(
  length_plot_data,
  aes(x = mean_depth, y = median_length)
) +
  geom_point(
    size = 2.5,
    alpha = 0.6
  ) +
  geom_ribbon(
    data = depth_predictions,
    aes(
      x = depth,
      ymin = conf.low,
      ymax = conf.high
    ),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  geom_line(
    data = depth_predictions,
    aes(
      x = depth,
      y = predicted
    ),
    inherit.aes = FALSE,
    linewidth = 1
  ) +
  labs(
    x = "Sampling depth (m)",
    y = "Median shrimp length (mm)"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11)
  )

# Display figure
fig_length_depth

# Display length model-selection results
length_model_selection %>%
  arrange(AICc) %>%
  mutate(
    AICc = round(AICc, 2),
    delta_AICc = round(delta_AICc, 2),
    weight = round(weight, 3)
  )


# ------------------------------------------------------------
# 17M. Diagnostic plots for the depth model
# ------------------------------------------------------------

par(mfrow = c(2, 2))

plot(length_lm_depth)

par(mfrow = c(1, 1))


# ------------------------------------------------------------
# 17N. Diagnostics for the primary depth model (91 events)
# ------------------------------------------------------------

# Model summary
summary(length_lm_depth)

# Residual diagnostic plots
par(mfrow = c(2, 2))
plot(length_lm_depth)
par(mfrow = c(1, 1))

# Residual normality
shapiro.test(residuals(length_lm_depth))

# Homogeneity of variance
lmtest::bptest(length_lm_depth)

# Multicollinearity
car::vif(length_lm_depth)


# ------------------------------------------------------------
# 17O. Influence assessment for the primary depth model
# ------------------------------------------------------------

# Recreate the exact dataset used by the model
depth_model_data <- length_btrw_env[
  rownames(model.frame(length_lm_depth)),
  ,
  drop = FALSE
]

# Calculate Cook's distances
depth_model_data$cooks_distance <- cooks.distance(length_lm_depth)
depth_model_data$model_row <- seq_len(nrow(depth_model_data))

# Identify potentially influential observations
depth_influence_threshold <- 4 / nobs(length_lm_depth)

influential_depth_events <- depth_model_data %>%
  filter(cooks_distance > depth_influence_threshold) %>%
  select(
    model_row,
    sample_id,
    Year,
    Month,
    Site,
    n_measured,
    median_length,
    mean_depth,
    cooks_distance
  ) %>%
  arrange(desc(cooks_distance))

influential_depth_events


depth_model_data_sensitivity <- depth_model_data %>%
  filter(cooks_distance <= depth_influence_threshold)

length_lm_depth_influence_check <- lm(
  median_length ~ year_decade + Month_f + depth_z,
  data = depth_model_data_sensitivity
)

summary(length_lm_depth_influence_check)

# Compare original and influence-check models
data.frame(
  Model = c("Original", "Without flagged observations"),
  Depth_estimate = c(
    coef(length_lm_depth)["depth_z"],
    coef(length_lm_depth_influence_check)["depth_z"]
  ),
  Depth_p_value = c(
    summary(length_lm_depth)$coefficients["depth_z", "Pr(>|t|)"],
    summary(length_lm_depth_influence_check)$coefficients[
      "depth_z", "Pr(>|t|)"
    ]
  ),
  R_squared = c(
    summary(length_lm_depth)$r.squared,
    summary(length_lm_depth_influence_check)$r.squared
  )
)


# ------------------------------------------------------------
# 17P. Sensitivity analysis: events with at least 2 shrimp
# ------------------------------------------------------------

# Confirm the restricted dataset
length_btrw_sensitivity %>%
  summarise(
    total_events = n(),
    total_shrimp_measured = sum(n_measured),
    minimum_measured = min(n_measured),
    maximum_measured = max(n_measured)
  )

# Model selection results
length_model_selection_sensitivity %>%
  mutate(
    across(c(AICc, delta_AICc, weight), ~ round(.x, 3))
  )

# Summary of the best-supported temperature + depth model
length_sensitivity_model <- length_models_sensitivity$temp_depth

summary(length_sensitivity_model)

# Confidence intervals
confint(length_sensitivity_model)

# HC3 robust standard errors
lmtest::coeftest(
  length_sensitivity_model,
  vcov. = sandwich::vcovHC(
    length_sensitivity_model,
    type = "HC3"
  )
)


# ------------------------------------------------------------
# 17Q. Diagnostics for the 63-event sensitivity model
# ------------------------------------------------------------

# Residual diagnostic plots
par(mfrow = c(2, 2))
plot(length_sensitivity_model)
par(mfrow = c(1, 1))

# Residual normality
shapiro.test(residuals(length_sensitivity_model))

# Homogeneity of variance
lmtest::bptest(length_sensitivity_model)

# Multicollinearity
performance::check_collinearity(length_sensitivity_model)


# ------------------------------------------------------------
# 17R. Influence assessment for the 63-event model
# ------------------------------------------------------------

# Calculate influence statistics
sensitivity_cooks_distance <- cooks.distance(length_sensitivity_model)
sensitivity_studentized_residuals <- rstudent(length_sensitivity_model)
sensitivity_leverage <- hatvalues(length_sensitivity_model)

# Screening threshold
sensitivity_influence_threshold <-
  4 / nobs(length_sensitivity_model)

# Identify potentially influential observations
influential_events <- length_btrw_sensitivity %>%
  mutate(
    observation = row_number(),
    cooks_distance = sensitivity_cooks_distance,
    studentized_residual = sensitivity_studentized_residuals,
    leverage = sensitivity_leverage
  ) %>%
  filter(cooks_distance > sensitivity_influence_threshold) %>%
  select(
    observation,
    median_length,
    n_measured,
    mean_depth,
    temp_z,
    depth_z,
    cooks_distance,
    studentized_residual,
    leverage
  ) %>%
  arrange(desc(cooks_distance))

influential_events


# ------------------------------------------------------------
# 17S. Influence sensitivity check
# ------------------------------------------------------------

# Refit the model after removing each flagged observation
# individually. Do not remove observations from the final
# model solely because they exceed the Cook's distance threshold.

flagged_observations <- influential_events$observation

influence_comparison <- lapply(
  flagged_observations,
  function(i) {
    
    model_without_event <- lm(
      formula(length_sensitivity_model),
      data = length_btrw_sensitivity[-i, ]
    )
    
    coefficients <- coef(model_without_event)
    
    data.frame(
      removed_event = i,
      temperature_estimate = coefficients["temp_z"],
      depth_estimate = coefficients["depth_z"],
      r_squared = summary(model_without_event)$r.squared
    )
  }
) %>%
  bind_rows()

influence_comparison


# ------------------------------------------------------------
# FIGURE: Predicted Median Length by Sampling Depth
# ------------------------------------------------------------

# Use the observations included in the selected model
depth_model_data <- model.frame(length_lm_depth)

depth_model_rows <- rownames(depth_model_data)

depth_original_data <- length_btrw_env[
  depth_model_rows, ,
  drop = FALSE
]

# Use the original depth standardization values
depth_mean <- mean(depth_original_data$mean_depth, na.rm = TRUE)

depth_sd <- sd(depth_original_data$mean_depth, na.rm = TRUE)

# Create prediction dataset
depth_prediction_data <- data.frame(
  depth_z = seq(
    min(depth_model_data$depth_z),
    max(depth_model_data$depth_z),
    length.out = 100
  ),
  year_decade = 0,
  Month_f = factor(
    "October",
    levels = levels(depth_model_data$Month_f)
  )
)

# Model predictions
depth_prediction_data$fit <- predict(
  length_lm_depth,
  newdata = depth_prediction_data
)

# HC3 robust covariance matrix
robust_vcov <- sandwich::vcovHC(
  length_lm_depth,
  type = "HC3"
)

# Prediction design matrix
X <- model.matrix(
  delete.response(terms(length_lm_depth)),
  data = depth_prediction_data
)

# Robust standard errors of fitted means
depth_prediction_data$SE <- sqrt(
  rowSums((X %*% robust_vcov) * X)
)

# 95% confidence intervals
critical_value <- qt(
  0.975,
  df = df.residual(length_lm_depth)
)

depth_prediction_data$lwr <-
  depth_prediction_data$fit -
  critical_value * depth_prediction_data$SE

depth_prediction_data$upr <-
  depth_prediction_data$fit +
  critical_value * depth_prediction_data$SE

# Convert standardized depth back to meters
depth_prediction_data$depth_m <-
  depth_prediction_data$depth_z * depth_sd + depth_mean

# Create figure
fig_length_depth <- ggplot() +
  geom_point(
    data = depth_original_data,
    aes(x = mean_depth, y = median_length),
    color = "#2878A5",
    alpha = 0.55,
    size = 2
  ) +
  geom_ribbon(
    data = depth_prediction_data,
    aes(x = depth_m, ymin = lwr, ymax = upr),
    fill = "grey70",
    alpha = 0.4
  ) +
  geom_line(
    data = depth_prediction_data,
    aes(x = depth_m, y = fit),
    color = "black",
    linewidth = 1
  ) +
  labs(
    x = "Sampling depth (m)",
    y = "Median shrimp length (mm)"
  ) +
  theme_classic(base_size = 12)

print(fig_length_depth)



# ============================================================
# 18. FINAL RESULTS TABLES
# ============================================================

# ------------------------------------------------------------
# 18A. Occurrence model-selection table
# ------------------------------------------------------------

occurrence_model_table <- occurrence_model_selection_btrw %>%
  mutate(
    Model = dplyr::recode(
      model,
      "base" = "Base",
      "month" = "Month",
      "temperature" = "Temperature",
      "salinity" = "Salinity",
      "depth" = "Depth",
      "temp_salinity" = "Temperature + salinity",
      "temp_depth" = "Temperature + depth",
      "salinity_depth" = "Salinity + depth",
      "full" = "Temperature + salinity + depth"
    )
  ) %>%
  select(
    Model,
    AICc,
    delta_AICc,
    weight
  ) %>%
  mutate(
    AICc = round(AICc, 2),
    delta_AICc = round(delta_AICc, 2),
    weight = round(weight, 3)
  )

print(occurrence_model_table)


# ------------------------------------------------------------
# 18B. Final occurrence model coefficients
# ------------------------------------------------------------

occurrence_coefficient_table <- occurrence_results %>%
  filter(
    term != "(Intercept)"
  ) %>%
  mutate(
    Predictor = dplyr::recode(
      term,
      "year_decade" = "Year (per decade)",
      "Month_fOctober" = "October",
      "salinity_z" = "Surface salinity (1 SD)"
    )
  ) %>%
  select(
    Predictor,
    beta,
    SE,
    OR,
    CI_low,
    CI_high,
    p_value
  ) %>%
  mutate(
    beta = round(beta, 3),
    SE = round(SE, 3),
    OR = round(OR, 3),
    CI_low = round(CI_low, 3),
    CI_high = round(CI_high, 3),
    p_value = case_when(
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    )
  ) %>%
  rename(
    Beta = beta,
    `Odds ratio` = OR,
    `95% CI lower` = CI_low,
    `95% CI upper` = CI_high,
    `p-value` = p_value
  )

print(occurrence_coefficient_table)


# ------------------------------------------------------------
# 18C. Positive-catch model-selection table
# ------------------------------------------------------------


positive_model_selection <- tibble(
  model = names(positive_models_btrw),
  AICc = sapply(
    positive_models_btrw,
    MuMIn::AICc
  )
) %>%
  arrange(AICc) %>%
  mutate(
    delta_AICc = AICc - min(AICc),
    weight = exp(-0.5 * delta_AICc) /
      sum(exp(-0.5 * delta_AICc))
  )

print(positive_model_selection)


positive_model_table <- positive_model_selection %>%
  mutate(
    Model = dplyr::recode(
      model,
      "base" = "Base",
      "temperature" = "Temperature",
      "salinity" = "Salinity",
      "depth" = "Depth",
      "salinity_depth" = "Salinity + depth",
      "full" = "Temperature + salinity + depth"
    )
  ) %>%
  select(
    Model,
    AICc,
    delta_AICc,
    weight
  ) %>%
  mutate(
    AICc = round(AICc, 2),
    delta_AICc = round(delta_AICc, 2),
    weight = round(weight, 3)
  )

print(positive_model_table)


# ------------------------------------------------------------
# 18D. Final positive-catch model coefficients
# ------------------------------------------------------------

# ------------------------------------------------------------
# Calculate final positive-catch model results
# ------------------------------------------------------------

positive_results <- data.frame(
  term = names(glmmTMB::fixef(final_positive_model)$cond),
  beta = unname(glmmTMB::fixef(final_positive_model)$cond),
  SE = sqrt(diag(vcov(final_positive_model)$cond))
) %>%
  mutate(
    rate_ratio = exp(beta),
    CI_low = exp(beta - 1.96 * SE),
    CI_high = exp(beta + 1.96 * SE),
    z_value = beta / SE,
    p_value = 2 * pnorm(abs(z_value), lower.tail = FALSE)
  )

print(positive_results)



positive_coefficient_table <- positive_results %>%
  filter(
    term != "(Intercept)"
  ) %>%
  mutate(
    Predictor = dplyr::recode(
      term,
      "year_decade" = "Year (per decade)",
      "Month_fOctober" = "October",
      "salinity_z" = "Surface salinity (1 SD)"
    )
  ) %>%
  select(
    Predictor,
    beta,
    SE,
    rate_ratio,
    CI_low,
    CI_high,
    p_value
  ) %>%
  mutate(
    beta = round(beta, 3),
    SE = round(SE, 3),
    rate_ratio = round(rate_ratio, 3),
    CI_low = round(CI_low, 3),
    CI_high = round(CI_high, 3),
    p_value = case_when(
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    )
  ) %>%
  rename(
    Beta = beta,
    `Rate ratio` = rate_ratio,
    `95% CI lower` = CI_low,
    `95% CI upper` = CI_high,
    `p-value` = p_value
  )

print(positive_coefficient_table)


# ------------------------------------------------------------
# 18E. Length model-selection table
# ------------------------------------------------------------

length_model_table <- length_model_selection %>%
  mutate(
    Model = dplyr::recode(
      model,
      "base" = "Base",
      "temperature" = "Temperature",
      "salinity" = "Salinity",
      "depth" = "Depth",
      "temp_salinity" = "Temperature + salinity",
      "temp_depth" = "Temperature + depth",
      "salinity_depth" = "Salinity + depth",
      "full" = "Temperature + salinity + depth"
    )
  ) %>%
  select(Model, AICc, delta_AICc, weight) %>%
  mutate(
    AICc = round(AICc, 2),
    delta_AICc = round(delta_AICc, 2),
    weight = round(weight, 3)
  )

print(length_model_table)


# ------------------------------------------------------------
# Print primary model tables
# ------------------------------------------------------------

print(occurrence_model_table)
print(occurrence_coefficient_table)
print(positive_model_table)
print(positive_coefficient_table)
print(length_model_table)


# ------------------------------------------------------------
# 18F. Final length model parameters (HC3 robust SE)
# ------------------------------------------------------------

# Set the selected model
length_lm_final <- length_lm_depth

# Calculate HC3 robust standard errors
length_robust_results <- lmtest::coeftest(
  length_lm_final,
  vcov. = sandwich::vcovHC(length_lm_final, type = "HC3")
)

# Convert results to a matrix
robust_matrix <- as.matrix(length_robust_results)

# Create coefficient table
length_parameter_table <- data.frame(
  term = rownames(robust_matrix),
  Beta = robust_matrix[, 1],
  SE = robust_matrix[, 2],
  p_value = robust_matrix[, 4]
) %>%
  mutate(
    critical_value = qt(
      0.975,
      df = df.residual(length_lm_final)
    ),
    CI_lower = Beta - critical_value * SE,
    CI_upper = Beta + critical_value * SE,
    Predictor = dplyr::recode(
      term,
      "(Intercept)" = "Intercept",
      "year_decade" = "Year (per decade)",
      "Month_fOctober" = "October",
      "depth_z" = "Sampling depth (1 SD)"
    )
  ) %>%
  select(
    Predictor,
    Beta,
    SE,
    CI_lower,
    CI_upper,
    p_value
  ) %>%
  mutate(
    across(
      c(Beta, SE, CI_lower, CI_upper),
      ~ round(.x, 3)
    ),
    p_value = case_when(
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    )
  )

print(length_parameter_table)


# ------------------------------------------------------------
# 18G. Length sensitivity analysis model selection
# Datset restricted to events with at least 2 shrimp
# ------------------------------------------------------------

# Confirm sensitivity dataset size and sampling restriction
length_btrw_sensitivity %>%
  summarise(
    n_events = n(),
    min_shrimp_measured = min(n_measured),
    max_shrimp_measured = max(n_measured),
    n_years = n_distinct(Year),
    first_year = min(Year),
    last_year = max(Year)
  )


names(length_model_selection_sensitivity)

head(as.data.frame(length_model_selection_sensitivity))


# ------------------------------------------------------------
# 18H. Length sensitivity analysis model selection table
# Included events with at least 2 shrimp (n = 63)
# ------------------------------------------------------------

length_sensitivity_table <- length_model_selection_sensitivity %>%
  mutate(
    Model = dplyr::recode(
      model,
      "base" = "Year + Month",
      "temperature" = "Year + Month + Temperature",
      "salinity" = "Year + Month + Salinity",
      "depth" = "Year + Month + Depth",
      "temp_salinity" = "Year + Month + Temperature + Salinity",
      "temp_depth" = "Year + Month + Temperature + Depth",
      "salinity_depth" = "Year + Month + Salinity + Depth",
      "full" = "Year + Month + Temperature + Salinity + Depth"
    )
  ) %>%
  select(
    Model,
    AICc,
    delta_AICc,
    weight
  ) %>%
  mutate(
    AICc = round(AICc, 2),
    delta_AICc = round(delta_AICc, 2),
    weight = round(weight, 3)
  )


print(length_sensitivity_table)


# ------------------------------------------------------------
# 18I. Length sensitivity analysis coefficient table
# ------------------------------------------------------------

# Extract the temperature + depth sensitivity model

length_sensitivity_best <- length_models_sensitivity$temp_depth
# Calculate HC3 robust standard errors
length_sensitivity_robust <- lmtest::coeftest(
  length_sensitivity_best,
  vcov. = sandwich::vcovHC(
    length_sensitivity_best,
    type = "HC3"
  )
)

# Convert results into a table
sensitivity_matrix <- as.matrix(length_sensitivity_robust)

length_sensitivity_table_coefficients <- data.frame(
  term = rownames(sensitivity_matrix),
  Beta = sensitivity_matrix[, 1],
  SE = sensitivity_matrix[, 2],
  p_value = sensitivity_matrix[, 4]
) %>%
  mutate(
    critical_value = qt(
      0.975,
      df = df.residual(length_sensitivity_best)
    ),
    CI_lower = Beta - critical_value * SE,
    CI_upper = Beta + critical_value * SE,
    Predictor = dplyr::recode(
      term,
      "(Intercept)" = "Intercept",
      "year_decade" = "Year (per decade)",
      "Month_fOctober" = "October",
      "temp_z" = "Temperature (1 SD)",
      "depth_z" = "Sampling depth (1 SD)"
    )
  ) %>%
  select(
    Predictor,
    Beta,
    SE,
    CI_lower,
    CI_upper,
    p_value
  ) %>%
  mutate(
    across(
      c(Beta, SE, CI_lower, CI_upper),
      ~ round(.x, 3)
    ),
    p_value = case_when(
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    )
  )

print(length_sensitivity_table_coefficients)

## Checking sensitivity diagnostics ##

# Residual diagnostics
sensitivity_residuals <- residuals(length_sensitivity_best)
sensitivity_fitted <- fitted(length_sensitivity_best)

# Normality
shapiro.test(sensitivity_residuals)

# Heteroscedasticity
lmtest::bptest(length_sensitivity_best)

# Multicollinearity
car::vif(length_sensitivity_best)

# Influence
cooks_distance_sensitivity <- cooks.distance(length_sensitivity_best)

which(
  cooks_distance_sensitivity >
    4 / nobs(length_sensitivity_best)
)

# ------------------------------------------------------------
# 18J. Influence sensitivity analysis
# ------------------------------------------------------------

# Identify observations exceeding the Cook's distance threshold
influential_rows <- which(
  cooks.distance(length_sensitivity_model) >
    4 / nobs(length_sensitivity_model)
)

print(influential_rows)

# Remove these observations temporarily
length_btrw_influence_check <- length_btrw_sensitivity %>%
  slice(-all_of(influential_rows))

# Refit the temperature + depth model
length_influence_model <- lm(
  median_length ~ year_decade + Month_f + temp_z + depth_z,
  data = length_btrw_influence_check
) 

nobs(length_influence_model) 

# Calculate HC3 robust standard errors
influence_robust <- lmtest::coeftest(
  length_influence_model,
  vcov. = sandwich::vcovHC(
    length_influence_model,
    type = "HC3"
  )
)

# Create coefficient table
influence_matrix <- as.matrix(influence_robust)

length_influence_table <- data.frame(
  term = rownames(influence_matrix),
  Beta = influence_matrix[, 1],
  SE = influence_matrix[, 2],
  p_value = influence_matrix[, 4]
) %>%
  mutate(
    critical_value = qt(
      0.975,
      df = df.residual(length_influence_model)
    ),
    CI_lower = Beta - critical_value * SE,
    CI_upper = Beta + critical_value * SE,
    Predictor = dplyr::recode(
      term,
      "(Intercept)" = "Intercept",
      "year_decade" = "Year (per decade)",
      "Month_fOctober" = "October",
      "temp_z" = "Temperature (1 SD)",
      "depth_z" = "Sampling depth (1 SD)"
    )
  ) %>%
  select(
    Predictor, Beta, SE, CI_lower, CI_upper, p_value
  ) %>%
  mutate(
    across(
      c(Beta, SE, CI_lower, CI_upper),
      ~ round(.x, 3)
    ),
    p_value = case_when(
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    )
  )

print(length_influence_table)



# ============================================================
# 19. MAPS
# ============================================================

# Include BTRW and BCHS100
spatial_data <- fall_mapped %>%
  filter(
    Gear %in% c("BTRW", "BCHS100"),
    !is.na(map_latitude),
    !is.na(map_longitude)
  ) %>%
  mutate(
    Period = case_when(
      Year < 2008 ~ "1990–2007",
      Year >= 2008 & Year <= 2012 ~ "2008–2012",
      Year >= 2013 ~ "2013–2018"
    ),
    Period = factor(
      Period,
      levels = c(
        "1990–2007",
        "2008–2012",
        "2013–2018"
      )
    ),
    Gear = factor(
      Gear,
      levels = c("BTRW", "BCHS100")
    )
  )

# Convert to spatial features
spatial_sf <- st_as_sf(
  spatial_data,
  coords = c("map_longitude", "map_latitude"),
  crs = 4326,
  remove = FALSE
)

# Check detections by period and gear
spatial_data %>%
  group_by(Period, Gear) %>%
  summarise(
    sampling_events = n(),
    positive_events = sum(presence == 1),
    .groups = "drop"
  ) %>%
  print(n = Inf)


## Basemap ##


# Define map boundaries
map_limits <- c(
  xmin = -75.25,
  xmax = -75.00,
  ymin = 37.95,
  ymax = 38.35
)

# Test basemap
test_basemap <- ggplot() +
  annotation_map_tile(
    type = "osm",
    zoomin = 0
  ) +
  geom_sf(
    data = spatial_sf,
    alpha = 0
  ) +
  coord_sf(
    crs = st_crs(4326),
    default_crs = st_crs(4326),
    xlim = c(-75.45, -75.00),
    ylim = c(37.95, 38.55),
    expand = FALSE
  ) +
  theme_bw()

test_basemap


##Map of Sampling Locations##
fig_spatial_1990_2018 <- ggplot() +
  
  # OpenStreetMap basemap
  annotation_map_tile(
    type = "osm",
    zoomin = 0
  ) +
  
  # All sampling locations
  geom_sf(
    data = spatial_sf,
    color = "gray45",
    size = 1.2,
    alpha = 0.5
  ) +
  
  # White shrimp detections colored by gear
  geom_sf(
    data = spatial_sf %>%
      filter(presence == 1),
    aes(color = Gear),
    size = 2,
    alpha = 0.9
  ) +
  
  # Gear-specific colors
  scale_color_manual(
    name = "Sampling gear",
    values = c(
      "BTRW" = "red3",
      "BCHS100" = "deepskyblue3"
    ),
    breaks = c("BTRW", "BCHS100"),
    labels = c(
      "Bottom trawl (BTRW)",
      "Beach seine (BCHS100)"
    )
  ) +
  
  # Separate maps by time period
  facet_wrap(
    ~ Period,
    ncol = 3
  ) +
  
  # Map boundaries
  coord_sf(
    crs = st_crs(4326),
    default_crs = st_crs(4326),
    xlim = c(-75.45, -75.00),
    ylim = c(37.95, 38.55),
    expand = FALSE
  ) +
  
  # Scale bar
  annotation_scale(
    location = "bl",
    width_hint = 0.25
  ) +
  
  # North arrow
  annotation_north_arrow(
    location = "tr",
    pad_x = unit(0.15, "in"),
    pad_y = unit(0.15, "in"),
    height = unit(0.40, "in"),
    width = unit(0.65, "in"),
    style = north_arrow_fancy_orienteering
  ) +
  
  # Labels
  labs(
    x = "Longitude",
    y = "Latitude",
    caption = paste(
      "Gray = sampling locations;",
      "Red = BTRW detections;",
      "Blue = BCHS100 detections"
    )
  ) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    panel.spacing = unit(0.6, "in"),
    
    axis.text.x = element_text(
      size = 9,
      angle = 45,
      hjust = 1
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    
    plot.caption = element_text(
      size = 10,
      hjust = 0.5,
      margin = margin(t = 12)
    )
  )

# Display map
fig_spatial_1990_2018

spatial_data %>%
  filter(presence == 1) %>%
  group_by(Period, Gear) %>%
  summarise(
    positive_events = n(),
    positive_sites = n_distinct(Site),
    .groups = "drop"
  )


spatial_data %>%
  group_by(Period, Gear) %>%
  summarise(
    total_sites = n_distinct(Site),
    positive_sites = n_distinct(Site[presence == 1]),
    percent_positive_sites =
      100 * positive_sites / total_sites,
    .groups = "drop"
  ) %>%
  print(n = Inf) 



fall_1990_2018 %>%
  mutate(
    Period = case_when(
      Year < 2008 ~ "1990–2007",
      Year <= 2012 ~ "2008–2012",
      Year <= 2018 ~ "2013–2018"
    )
  ) %>%
  group_by(Period, Gear) %>%
  summarise(
    total_sites = n_distinct(Site),
    positive_sites = n_distinct(Site[presence == 1]),
    percent_positive_sites =
      round(100 * positive_sites / total_sites, 1),
    .groups = "drop"
  ) %>%
  print(n = Inf)

# ============================================================
# 20. FINALIZED FIGURES
# ============================================================


print(fig_salinity)
print(fig_positive_salinity)


# ------------------------------------------------------------
# Median shrimp length over time
# ------------------------------------------------------------

length_year_summary <- length_by_event %>%
  filter(Gear == "BTRW", !is.na(median_length)) %>%
  group_by(Year) %>%
  summarise(
    annual_median = median(median_length),
    .groups = "drop"
  )

fig_length_year <- ggplot(
  length_by_event,
  aes(x = factor(Year), y = median_length)
) +
  geom_jitter(
    width = 0.15,
    height = 0,
    color = "black",
    size = 2,
    alpha = 1
  ) +
  scale_y_continuous(
    breaks = seq(20, 140, by = 10)
  ) +
  labs(
    x = "Year",
    y = "Median shrimp length (mm)"
  ) +
  theme_classic()

print(fig_length_year)


# ============================================================
# 21. FINAL TABLES
# ============================================================

white_shrimp_main_tables <- list(
  
  # Occurrence
  occurrence_model_selection = occurrence_model_table,
  occurrence_coefficients = occurrence_results,
  
  # Positive catch abundance
  positive_model_selection = positive_model_table,
  positive_coefficients = positive_results,
  
  # Length
  length_model_selection = length_model_table,
  length_coefficients = length_parameter_table
  
)


#Sensitivity tables#

white_shrimp_appendix_tables <- list(
  
  # Length sensitivity analysis
  length_sensitivity_model_selection = length_sensitivity_table,
  length_sensitivity_coefficients = length_sensitivity_table_coefficients,
  
  # Influence sensitivity analysis
  length_influence_coefficients = length_influence_table
  
) 


# Occurrence model selection
white_shrimp_main_tables$occurrence_model_selection

# Occurrence coefficients
white_shrimp_main_tables$occurrence_coefficients

# Positive catch model selection
white_shrimp_main_tables$positive_model_selection

# Positive catch coefficients
white_shrimp_main_tables$positive_coefficients

# Length model selection
white_shrimp_main_tables$length_model_selection

# Length coefficients
white_shrimp_main_tables$length_coefficients


# ============================================================
# 22. TABLE AND FIGURE DISPLAY
# ============================================================

# ------------------------------------------------------------
# OCCURRENCE
# ------------------------------------------------------------

# Updated model selection (22 candidate models)
print(occurrence_model_selection_btrw, n = Inf)

# Updated final occurrence coefficients
print(occurrence_coefficient_table_final)

# NAO model comparison
print(NAO_model_comparison)

# Predicted occurrence by period and month
print(fig_occurrence)

# Predicted occurrence across salinity
print(fig_salinity)


# ------------------------------------------------------------
# POSITIVE-CATCH ABUNDANCE
# ------------------------------------------------------------

# Positive-catch model selection
print(positive_model_table)

# Positive-catch coefficients
print(positive_coefficient_table)

# Annual catch rate by sampling month
print(fig_monthly_annual_abundance)

# Predicted positive-catch abundance across salinity
print(fig_positive_salinity)


# ------------------------------------------------------------
# SHRIMP LENGTH
# ------------------------------------------------------------

print(length_model_table)
print(length_parameter_table)

print(fig_length_year)
print(fig_length_depth)


# ------------------------------------------------------------
# SPATIAL DISTRIBUTION
# ------------------------------------------------------------

print(fig_spatial_1990_2018)


# ------------------------------------------------------------
# ADDITIONAL SENSITIVITY TABLES
# ------------------------------------------------------------

print(length_sensitivity_table)
print(length_sensitivity_table_coefficients)
print(length_influence_table)


# Exact predicted occurrence by period and month
occ_predictions_df %>%
  select(x, group, predicted, conf.low, conf.high) %>%
  print()

# Updated salinity predictions near 20 and 35 PSU
salinity_predictions_df %>%
  mutate(
    distance_20 = abs(salinity_PSU - 20),
    distance_35 = abs(salinity_PSU - 35)
  ) %>%
  filter(
    distance_20 == min(distance_20) |
      distance_35 == min(distance_35)
  ) %>%
  select(
    salinity_PSU,
    predicted,
    conf.low,
    conf.high
  ) %>%
  print()



##Sampling events, detections, and occurrence graphs/table from preliminary anaysis## 

##Calculate occurrence by period and gear  


prelim_occurrence <- fall_1990_2018 %>%
  mutate(
    Period = if_else(Year <= 2007,
                     "1990–2007",
                     "2008–2018")
  ) %>%
  group_by(Period, Gear) %>%
  summarise(
    total_events = n(),
    positive_events = sum(presence, na.rm = TRUE),
    total_abundance = sum(abundance, na.rm = TRUE),
    occurrence_percent = round(
      100 * positive_events / total_events, 2
    ),
    .groups = "drop"
  )


print(prelim_occurrence)

##Annual occurrence percentages and graph


annual_occurrence_sep <- fall_1990_2018 %>%
  filter(Month == 9) %>%
  group_by(Year, Gear) %>%
  summarise(
    total_events = n(),
    positive_events = sum(presence, na.rm = TRUE),
    occurrence_percent = 100 * positive_events / total_events,
    .groups = "drop"
  )

ggplot(
  annual_occurrence_sep,
  aes(x = Year,
      y = occurrence_percent,
      color = Gear)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  
  geom_vline(
    xintercept = 2008,
    linetype = "dashed",
    color = "gray40"
  ) +
  
  scale_color_manual(
    values = c(
      "BTRW" = "#0072B2",
      "BCHS100" = "#D55E00"
    )
  ) +
  
  scale_x_continuous(
    breaks = seq(1990, 2018, 2)
  ) +
  
  labs(
    x = "Year",
    y = "White shrimp occurrence (%)",
    color = "Sampling gear"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "top"
  )


