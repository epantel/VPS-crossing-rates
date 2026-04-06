```{r packages}
library(remotes)
install_github('ocean-tracking-network/glatos', build_vignettes = TRUE)
library(glatos)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(sf)
library(purrr)
library(ISOweek)
library(lme4)
library(arm)
library(lattice)
library(rv)
```



```{r reading in detections}
#this is the standard telemetry data
detections23<-read_glatos_detections("LENAT_detectionsWithLocs_20231128_200102.csv")

#positional data from vps
vps23detections <-read.csv("positions2023_all.csv")
vps23detections <- vps23detections %>%
  filter(FullId %in% det23$FullId)

matching_indices <- match(vps23detections$FullId, det23$FullId)

vps23detections$Species <- det23$common_name_e[matching_indices]

vps23detections <- vps23detections %>%
  filter(!(tolower(Species) %in% c("bowfin", "silver redhorse", "bigmouth buffalo")))
vps23detections$DETECTEDID <- vps23detections$FullId
vps23detections$LAT <- vps23detections$Latitude
vps23detections$LON <- vps23detections$Longitude
vps23detections$DATETIME <- as.POSIXct(vps23detections$Time, format = "%Y-%m-%d %H:%M:%S")

vps22detections <-read.csv("natives_BradysIsland_2022.csv")

vps22detections$DETECTEDID <- vps22detections$FullId
vps22detections$LAT <- vps22detections$Latitude
vps22detections$LON <- vps22detections$Longitude
vps22detections$DATETIME <- as.POSIXct(vps22detections$Time, format = "%m/%d/%Y %H:%M")



combined_vps <- rbind(vps23detections, vps22detections)

combined_vps$Month <- month(combined_vps$DATETIME)
combined_vps$Year <- year(combined_vps$DATETIME)
combined_vps$new_id <- paste(combined_vps$FullId, combined_vps$Year, sep = "-")


combined_vps <- combined_vps %>%
  filter(FullId %in% det23$FullId)

combined_vps$animal_id <- det23$animal_id[matching_indices]


write.csv(combined_vps,"combined_vps.csv", row.names=F )

```

```{r filtering false detections}
det23 <- false_detections(detections23, 3600, show_plot = F)
det23$FullId <-paste(det23$transmitter_codespace, det23$transmitter_id, sep="-")
det23<-subset(det23, passed_filter==1)
```

```{r data management}
det23$Month<-month(det23$detection_timestamp_utc)
det23$Year<-year(det23$detection_timestamp_utc)
det23$MonthYear <- as.factor(paste(det23$Month, det23$Year, sep = "-"))
```

```{r read combined vps}
combined_vps <- read.csv("combined_vps.csv")
```

```{r supplemental figures data management}
#general receier locations for supplemental figures
det23<-det23 %>% mutate(glatos_array2 =
                 case_when (glatos_array=="USR" ~ "Sandusky River - Above Deterrent",
                            station=="LSR-021" ~"Sandusky River - Above Deterrent",
                            glatos_array=="SBI" ~ "Sandusky River - Above Deterrent",
                            glatos_array=="SRVPS" ~ "Sandusky River - Below Deterrent",
                            glatos_array=="CIW" ~ "Sandusky Bay",
                            station=="CIW-001" ~ "Lake Erie and Other Tributaries",
                            station=="CIW-002" ~ "Lake Erie and Other Tributaries",
                            station=="CIW-003" ~ "Lake Erie and Other Tributaries",
                            station=="CIW-005" ~ "Lake Erie and Other Tributaries",
                            glatos_array=="SRVPS" ~ "Lake Erie and Other Tributaries",
                            glatos_array=="SDB" ~ "Sandusky Bay",
                            glatos_array=="OGN" ~ "Lake Erie and Other Tributaries",
                            glatos_array=="LSR" ~ "Sandusky River - Below Deterrent",
                            glatos_array=="LEG" ~ "Lake Erie and Other Tributaries",
                            glatos_array=="ENS" ~ "Lake Erie and Other Tributaries",
                            glatos_array=="WBA" ~ "Lake Erie and Other Tributaries",
                            glatos_array=="WRC" ~ "Lake Erie and Other Tributaries",
                            glatos_array=="TSR" ~ "Lake Erie and Other Tributaries",
                            glatos_array=="LPB" ~ "Lake Erie and Other Tributaries", 
                            glatos_array=="TRC" ~ "Lake Erie and Other Tributaries",
                            glatos_array=="MDR" ~ "Lake Erie and Other Tributaries",
                            glatos_array=="MAB" ~ "Lake Erie and Other Tributaries",
                            glatos_array=="MAU" ~ "Lake Erie and Other Tributaries",
                            glatos_array=="CUY" ~ "Lake Erie and Other Tributaries",
                            glatos_array=="HEN" ~ "Lake Erie and Other Tributaries",
                            glatos_array=="LL" ~ "Lake Erie and Other Tributaries"))
```


# White Bass Data
```{r white bass supplemental figure}
wball <- subset(det23, common_name_e=="white bass"|common_name_e =="White Bass"|common_name_e =="White bass")
values_to_remove <- c('x_071', 'x_072', 'x_076', 'x_079', 'x_080', 'x_085', 'x_104', 'x_108', 'x_109', 'x_076', 'x_111', 'x_115', 'x_121', 'x_123', 'x_125', 'x_127', 'x_161', 'x_166', 'x_169', 'x_172', 'x_183', 'x_184') #all did not live past 10 days, removed
wball <- wball[!(wball$animal_id %in% values_to_remove), ]
wball<-wball%>%
   filter(!(animal_id == "x_089" & Year==2023)) #definitely dead/expelled tag but failed the filtering
wball<-wball%>%
   filter(!(animal_id == "x_092" & Year==2023)) #definitely dead/expelled tag but failed the filtering
wb_det <- wball %>%
  group_by(glatos_array2, Month) %>%
  summarize(unique_animal_count = n_distinct(animal_id))

location_order <- c("Sandusky River - Above Deterrent", "Sandusky River - Below Deterrent", "Sandusky Bay", "Lake Erie and Other Tributaries")

wb_det$Location <- factor(wb_det$glatos_array2, levels = location_order)

ggplot(wb_det, aes(x = factor(Month, levels = 1:12, labels = substr(month.abb, 1, 3)), y = Location, size = unique_animal_count)) +
  geom_point() +
  scale_size_continuous(name = "Animal Count", breaks = seq(min(wb_det$unique_animal_count), max(wb_det$unique_animal_count), by = 5)) +
  labs(x = "Month", y = "Location", title = "Unique Detections of White Bass per Month")+
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),  
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(size = 20, hjust = 0.9),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16))

```

```{r filtering white bass monthly}
find_7_day_run <- function(dates) {
  dates <- sort(unique(as.Date(dates)))
  if (length(dates) < 7) return(FALSE)
  
  for (i in 1:(length(dates) - 6)) {
    if (all(seq(dates[i], by = "1 day", length.out = 7) %in% dates)) {
      return(TRUE)
    }
  }
  return(FALSE)
}


wbsall_monthly <- wball %>%
  mutate(
    detection_date = as.Date(detection_timestamp_utc),
    detection_month = floor_date(detection_date, "month")
  )

fish_list_wb <- wbsall_monthly %>%
  arrange(FullId, detection_date) %>%
  group_split(FullId)


fish_valid_month_wb <- map_dfr(fish_list_wb, function(fish_data) {
  fish_id <- unique(fish_data$FullId)
  
  first_date <- min(fish_data$detection_date)
  last_date <- max(fish_data$detection_date)
  all_months <- seq(floor_date(first_date, "month"), floor_date(last_date, "month"), by = "1 month")
  
  station_months <- fish_data %>%
    group_by(station, detection_month) %>%
    summarise(dates = list(unique(detection_date)), .groups = "drop")

  exclusion_month_A <- NULL
  for (i in 1:nrow(station_months)) {
    dates <- station_months$dates[[i]]
    if (length(dates) < 7) next

    month_start <- floor_date(min(dates), "month")
    first_7_days <- seq(month_start, by = "1 day", length.out = 7)
    match_dates <- intersect(dates, first_7_days)

    if (length(match_dates) >= 7 && find_7_day_run(match_dates)) {
      station_id <- station_months$station[i]
      
      # Check for future detections on other stations
      future_data <- fish_data %>% filter(detection_date > max(first_7_days))
      future_stations <- unique(future_data$station)
      
      if (length(future_stations[future_stations != station_id]) == 0) {
        exclusion_month_A <- month_start
        break
      }
    }
  }

  exclusion_month_B <- NULL
  if (is.null(exclusion_month_A)) {
    for (i in 2:length(all_months)) {
      prev_month <- all_months[i - 1]
      current_month <- all_months[i]

     if (!any(floor_date(fish_data$detection_date, "month") == prev_month)) next
      prev_month_data <- fish_data %>%
        filter(floor_date(detection_date, "month") == prev_month)

      station_groups <- prev_month_data %>%
        group_by(station) %>%
        summarise(dates = list(unique(detection_date)), .groups = "drop")

      for (j in 1:nrow(station_groups)) {
        dates <- station_groups$dates[[j]]
        if (length(dates) < 7) next

        if (find_7_day_run(dates)) {
          station_id <- station_groups$station[j]
          future_data <- fish_data %>% filter(detection_date > max(dates))
          future_stations <- unique(future_data$station)

          if (length(future_stations[future_stations != station_id]) == 0) {
            exclusion_month_B <- current_month
            break
          }
        }
      }

      if (!is.null(exclusion_month_B)) break
    }
  }

  exclusion_month <- min(c(exclusion_month_A, exclusion_month_B), na.rm = TRUE)
  if (is.infinite(exclusion_month)) exclusion_month <- NA

  valid_months <- if (is.na(exclusion_month)) {
    all_months
  } else {
    all_months[all_months < exclusion_month]
  }

  tibble(FullId = fish_id, valid_month = valid_months)
})

monthly_summary_wb <- fish_valid_month_wb %>%
  group_by(valid_month) %>%
  summarise(unique_fish = n_distinct(FullId), .groups = "drop") %>%
  arrange(valid_month)

print(monthly_summary_wb)

month_summary_wb <- fish_valid_month_wb %>%
  mutate(month_name = format(valid_month, "%B")) %>%  
  group_by(month_name) %>%
  summarise(total_fish_months = n(), .groups = "drop") %>%
  arrange(match(month_name, month.name))



fish_valid_month_wb <- fish_valid_month_wb  %>%
  mutate(year=year(valid_month))  %>%
  mutate(new_id=paste0(FullId, "-", year)) %>%
  mutate(Month=month(valid_month))

```

```{r white bass vps monthly crossing}
wbvps <- subset(combined_vps, Species=="white bass"|Species =="White Bass"|Species =="White bass")
values_to_remove <- c('x_071', 'x_072', 'x_076', 'x_079', 'x_080', 'x_085', 'x_104', 'x_108', 'x_109', 'x_076', 'x_111', 'x_115', 'x_121', 'x_123', 'x_125', 'x_127', 'x_161', 'x_166', 'x_169', 'x_172', 'x_183', 'x_184') #determined not survived
wbvps <- wbvps[!(wbvps$animal_id %in% values_to_remove), ]
wbvps<-wbvps%>%
   filter(!(animal_id == "x_089" & Year==2023)) #definitely dead/expelled tag but failed the filtering
wbvps<-wbvps%>%
   filter(!(animal_id == "x_092" & Year==2023)) #definitely dead/expelled tag but failed the filtering


#zones <- st_read("barrier buffer.shp")
wbvps_sf <- st_as_sf(wbvps, coords = c("LON", "LAT"), crs = st_crs(zones))

wb_joined <- st_join(wbvps_sf, zones)
wb_joined$Week <- week(wb_joined$DATETIME)

wbvps_joined <- wb_joined %>%
  filter(Zone != 2 & Zone != 4)

wbvps_joined<-as.data.frame(wbvps_joined)

wbvps_counts <- wbvps_joined %>%
  group_by(Month, new_id) %>%
  summarise(
    count = sum(Zone != lag(Zone), na.rm = TRUE),
    downstream_crossing = sum(Zone == 1 & lag(Zone) == 3, na.rm = TRUE),
    upstream_crossing = sum(Zone == 3 & lag(Zone) == 1, na.rm = TRUE)
  )

wb_completed <- fish_valid_month_wb %>%
  left_join(wbvps_counts, by = c("new_id", "Month")) %>%
  mutate(
    count = replace_na(count, 0),
    downstream_crossing = replace_na(downstream_crossing, 0),
    upstream_crossing = replace_na(upstream_crossing, 0)
  )
wb_completed$number_crossed <- ifelse(wb_completed$count > 0, 1, 0)
wb_completed$number_notcrossed <- ifelse(wb_completed$count == 0, 1, 0)

wbvps_crossing <- wb_completed %>% 
  group_by(Month) %>% 
  summarise(mean_crossing = mean(count, na.rm = T), SD = sd(count, na.rm = T), standard_error = sd(count, na.rm = T)/sqrt(length(count)))

wbvps_crossing <- wb_completed %>%
  group_by(Month) %>%
  summarise(
    mean_upstream = mean(upstream_crossing, na.rm = TRUE),
    mean_downstream = mean(downstream_crossing, na.rm = TRUE),
    
    var_upstream = var(upstream_crossing, na.rm = TRUE),
    var_downstream = var(downstream_crossing, na.rm = TRUE),)

```

```{r filtering whtie bass per week}

wball_clean <- wball %>%
  dplyr::select(FullId, detection_timestamp_utc) %>%
  mutate(
    week = isoweek(detection_timestamp_utc),
    year = year(detection_timestamp_utc),
    day_of_week = wday(detection_timestamp_utc, label = TRUE, abbr = TRUE)
  )


wball_clean <- wball_clean %>%
  arrange(FullId, detection_timestamp_utc)


wball_clean <- wball_clean %>%
  group_by(FullId, week, year) %>%
  mutate(
    day_diff = as.numeric(detection_timestamp_utc - lag(detection_timestamp_utc, default = first(detection_timestamp_utc))),
    consecutive_days = ifelse(day_diff == 1, 1, 0),  
    consecutive_count = cumsum(consecutive_days == 1)
  ) %>%
  ungroup()


wball_clean <- wball_clean %>%
  group_by(FullId, week, year) %>%
  mutate(
    consecutive_3_days = ifelse(consecutive_count >= 3, TRUE, FALSE) 
  ) %>%
  ungroup()


remove_dead_fish <- function(fish_data) {
  fish_data <- fish_data %>%
    arrange(detection_timestamp_utc)
  
  
  fish_data <- fish_data %>%
    mutate(
      
      dead_fish_flag = ifelse(consecutive_3_days & !any(detection_timestamp_utc > max(detection_timestamp_utc)), TRUE, FALSE)
    )
  

  fish_data <- fish_data %>%
    mutate(
      removed_from_current_week = ifelse(dead_fish_flag & wday(detection_timestamp_utc) <= 3, TRUE, FALSE)
    )
  

  fish_data <- fish_data %>%
    mutate(
      removed_from_future_weeks = ifelse(dead_fish_flag & wday(detection_timestamp_utc) >= 5, TRUE, FALSE)
    )
  
  return(fish_data)
}


wball_clean <- wball_clean %>%
  group_by(FullId, week, year) %>%
  do(remove_dead_fish(.)) %>%
  ungroup()

wball_clean_filtered <- wball_clean %>%
  filter(!(removed_from_current_week | removed_from_future_weeks))


wball_weekly_filtered <- wball_clean_filtered %>%
  distinct(FullId, year, week)


wball_weekly_filtered <- wball_weekly_filtered %>%
  mutate(week_id = year * 100 + week)


fish_week_range_filtered_wb <- wball_weekly_filtered %>%
  group_by(FullId) %>%
  summarise(
    min_week = min(week_id),
    max_week = max(week_id),
    .groups = "drop"
  )
generate_weeks <- function(min_week, max_week) {
  seq(min_week, max_week)
}

full_weeks_per_fish_filtered_wb <- fish_week_range_filtered_wb %>%
  rowwise() %>%
  mutate(weeks = list(generate_weeks(min_week, max_week))) %>%
  ungroup() %>%
  tidyr::unnest(weeks) %>%
  mutate(
    year = floor(weeks / 100),
    week = weeks %% 100
  ) %>%
  dplyr::select(FullId, year, week)

full_weeks_per_fish_filtered_wb<-full_weeks_per_fish_filtered_wb %>%
   mutate(new_id=paste0(FullId, "-", year)) %>%
  rename(Week = week)

fish_per_week_filtered_wb <- full_weeks_per_fish_filtered_wb %>%
  group_by(Week) %>%
  summarise(unique_fish = n_distinct(FullId)) %>%
  arrange(Week)

```

```{r white bass weekly crossing}
wbvps_counts_weekly <- wbvps_joined %>%
  mutate(Week = week(DATETIME)) %>% 
  mutate(Day = day(DATETIME)) %>%
  group_by(Week, new_id) %>%
  summarise(
    count = sum(Zone != lag(Zone), na.rm = TRUE),
    downstream_crossing = sum(Zone == 1 & lag(Zone) == 3, na.rm = TRUE),
    upstream_crossing = sum(Zone == 3 & lag(Zone) == 1, na.rm = TRUE))

wbvps_completed <- full_weeks_per_fish_filtered_wb %>%
  left_join(wbvps_counts_weekly, by = c("new_id", "Week")) %>%
  mutate(
    count = replace_na(count, 0),
    downstream_crossing = replace_na(downstream_crossing, 0),
    upstream_crossing = replace_na(upstream_crossing, 0)
  )

wbvps_completed$number_crossed <- ifelse(wbvps_completed$count > 0, 1, 0)
wbvps_completed$number_notcrossed <- ifelse(wbvps_completed$count == 0, 1, 0)


wbvps_crossingrate_weekly <- wbvps_completed %>%
  group_by(Week) %>%
  summarise(
    mean_upstream = mean(upstream_crossing, na.rm = TRUE),
    mean_downstream = mean(downstream_crossing, na.rm = TRUE),var_upstream = var(upstream_crossing, na.rm = TRUE),
    var_downstream = var(downstream_crossing, na.rm = TRUE),)
```


# Channel Catfish Data
```{r channel catfish supplemental figure}
channelall <- subset(det23, common_name_e=="channel catfish"|common_name_e =="Channel catfish"|common_name_e == "Channel Catfish")
channel_det <- channelall %>%
  group_by(Month,glatos_array2) %>%
  summarize(unique_animal_count = n_distinct(animal_id))

location_order <- c("Sandusky River - Above Deterrent", "Sandusky River - Below Deterrent", "Sandusky Bay", "Lake Erie and Other Tributaries")

channel_det$Location <- factor(channel_det$glatos_array2, levels = location_order)

ggplot(channel_det, aes(x = factor(Month, levels = 1:12, labels = substr(month.abb, 1, 3)), y = Location, size = unique_animal_count)) +
  geom_point() +
  scale_size_continuous(name = "Unique Animal Count", breaks = seq(min(channel_det$unique_animal_count), max(channel_det$unique_animal_count), by = 1)) +
  labs(x = "Month", y = "Location", title = "Unique Detections of Channel Catfish per Month") +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(size = 20, hjust = 0.5)
  )

library(ggplot2)
library(stringr)

ggplot(channel_det, aes(x = factor(Month, levels = 1:12, labels = substr(month.abb, 1, 3)), y = Location, size = unique_animal_count)) +
  geom_point() +
  scale_size_continuous(name = "Animal Count", breaks = seq(min(channel_det$unique_animal_count), max(channel_det$unique_animal_count), by = 1)) +
  labs(x = "Month", y = "Location", title = "Unique Detections of Channel Catfish per Month") +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),  
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(size = 20, hjust = 0.9),
    legend.title = element_text(size = 18),  
    legend.text = element_text(size = 16)    
  ) +
  scale_x_discrete(labels = function(labels) str_wrap(labels, width = 20)) + 
  scale_y_discrete(labels = function(labels) str_wrap(labels, width = 20))    
```

```{r filtering catfish monthly}
channelall <- channelall %>%
  mutate(
    detection_date = as.Date(detection_timestamp_utc),
    detection_month = floor_date(detection_date, "month")
  )

# Split data by fish
fish_list <- channelall %>%
  arrange(FullId, detection_date) %>%
  group_split(FullId)

# Process each fish
fish_valid_month_cc <- map_dfr(fish_list, function(fish_data) {
  fish_id <- unique(fish_data$FullId)
  
  # Get first and last detection months
  first_date <- min(fish_data$detection_date)
  last_date <- max(fish_data$detection_date)
  all_months <- seq(floor_date(first_date, "month"), floor_date(last_date, "month"), by = "1 month")
  
  # Pre-summarize detections per station/month
  station_months <- fish_data %>%
    group_by(station, detection_month) %>%
    summarise(dates = list(unique(detection_date)), .groups = "drop")

  ### --- Rule A: Early month 7-day run ---
  exclusion_month_A <- NULL
  for (i in 1:nrow(station_months)) {
    dates <- station_months$dates[[i]]
    if (length(dates) < 7) next

    month_start <- floor_date(min(dates), "month")
    first_7_days <- seq(month_start, by = "1 day", length.out = 7)
    match_dates <- intersect(dates, first_7_days)

    if (length(match_dates) >= 7 && find_7_day_run(match_dates)) {
      station_id <- station_months$station[i]
      
      # Check for future detections on other stations
      future_data <- fish_data %>% filter(detection_date > max(first_7_days))
      future_stations <- unique(future_data$station)
      
      if (length(future_stations[future_stations != station_id]) == 0) {
        exclusion_month_A <- month_start
        break
      }
    }
  }

  ### --- Rule B: Prior month 7-day run ---
  exclusion_month_B <- NULL
  if (is.null(exclusion_month_A)) {
    for (i in 2:length(all_months)) {
      prev_month <- all_months[i - 1]
      current_month <- all_months[i]

     if (!any(floor_date(fish_data$detection_date, "month") == prev_month)) next
      prev_month_data <- fish_data %>%
        filter(floor_date(detection_date, "month") == prev_month)

      station_groups <- prev_month_data %>%
        group_by(station) %>%
        summarise(dates = list(unique(detection_date)), .groups = "drop")

      for (j in 1:nrow(station_groups)) {
        dates <- station_groups$dates[[j]]
        if (length(dates) < 7) next

        if (find_7_day_run(dates)) {
          station_id <- station_groups$station[j]
          future_data <- fish_data %>% filter(detection_date > max(dates))
          future_stations <- unique(future_data$station)

          if (length(future_stations[future_stations != station_id]) == 0) {
            exclusion_month_B <- current_month
            break
          }
        }
      }

      if (!is.null(exclusion_month_B)) break
    }
  }

  # Decide final exclusion month
  exclusion_month <- min(c(exclusion_month_A, exclusion_month_B), na.rm = TRUE)
  if (is.infinite(exclusion_month)) exclusion_month <- NA

  # Valid months = before exclusion
  valid_months <- if (is.na(exclusion_month)) {
    all_months
  } else {
    all_months[all_months < exclusion_month]
  }

  tibble(FullId = fish_id, valid_month = valid_months)
})

monthly_summary <- fish_valid_month_cc %>%
  group_by(valid_month) %>%
  summarise(unique_fish = n_distinct(FullId), .groups = "drop") %>%
  arrange(valid_month)

print(monthly_summary)

month_summary_cc <- fish_valid_month_cc %>%
  mutate(month_name = format(valid_month, "%B")) %>%  # Extract month name
  group_by(month_name) %>%
  summarise(total_fish_months = n(), .groups = "drop") %>%
  arrange(match(month_name, month.name))  # Correct calendar order




fish_valid_month_cc <- fish_valid_month_cc  %>%
  mutate(year=year(valid_month))  %>%
  mutate(new_id=paste0(FullId, "-", year)) %>%
  mutate(Month=month(valid_month))
```

```{r catfish monthly crossings}
channelvps <- subset(combined_vps, Species=="channel catfish"|Species =="Channel catfish"|Species == "Channel Catfish")
channel_pop<-month_summary_cc$total_fish_months

channelvps_sf <- st_as_sf(channelvps, coords = c("LON", "LAT"), crs = st_crs(zones))

channel_joined <- st_join(channelvps_sf, zones)

channel_joined <- channel_joined %>%
  filter(Zone != 2 & Zone != 4)

channel_joined<-as.data.frame(channel_joined)

channelvps_counts <- channel_joined %>%
  group_by(Month, new_id) %>%
    summarise(
    count = sum(Zone != lag(Zone), na.rm = TRUE),
    downstream_crossing = sum(Zone == 1 & lag(Zone) == 3, na.rm = TRUE),
    upstream_crossing = sum(Zone == 3 & lag(Zone) == 1, na.rm = TRUE)
  )

channel_completed <- fish_valid_month_cc %>%
  left_join(channelvps_counts, by = c("new_id", "Month")) %>%
  mutate(
    count = replace_na(count, 0),
    downstream_crossing = replace_na(downstream_crossing, 0),
    upstream_crossing = replace_na(upstream_crossing, 0)
  )
channel_completed$number_crossed <- ifelse(channel_completed$count > 0, 1, 0)
channel_completed$number_notcrossed <- ifelse(channel_completed$count == 0, 1, 0)

  
channelvps_crossing <- channelvps_counts %>%
  group_by(Month) %>%
  summarise(
    mean_upstream = mean(upstream_crossing, na.rm = TRUE),
    mean_downstream = mean(downstream_crossing, na.rm = TRUE),
    var_upstream = var(upstream_crossing, na.rm = TRUE),
    var_downstream = var(downstream_crossing, na.rm = TRUE),
    )
```

# Smallmouth Buffalo Data
```{r sm buffalo supplemental figure}
smbuffall <- subset(det23, common_name_e=="smallmouth buffalo")
smbuffall<-subset(smbuffall, Year!=2021)
smbuff_det <- smbuffall %>%
  group_by(glatos_array2, Month) %>%
  summarize(unique_animal_count = n_distinct(animal_id))


smbuff_det$Location <- factor(smbuff_det$glatos_array2, levels = location_order)

ggplot(smbuff_det, aes(x = factor(Month, levels = 1:12, labels = substr(month.abb, 1, 3)), y = Location, size = unique_animal_count)) +
  geom_point() +
  scale_size_continuous(name = "Animal Count", breaks = seq(min(smbuff_det$unique_animal_count), max(smbuff_det$unique_animal_count), by = 3)) +
  labs(x = "Month", y = "Location", title = "Unique Detections of Smallmouth Buffalo per Month") +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),  
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  ) +
  scale_x_discrete(labels = function(labels) str_wrap(labels, width = 20)) +
  scale_y_discrete(labels = function(labels) str_wrap(labels, width = 20))

```

```{r filtering buffalo monthly}
smbuffall_monthly <- smbuffall %>%
  mutate(
    detection_date = as.Date(detection_timestamp_utc),
    detection_month = floor_date(detection_date, "month")
  )

fish_list_smbuff <- smbuffall_monthly %>%
  arrange(FullId, detection_date) %>%
  group_split(FullId)

fish_valid_month_smbuff <- map_dfr(fish_list_smbuff, function(fish_data) {
  fish_id <- unique(fish_data$FullId)
  
  first_date <- min(fish_data$detection_date)
  last_date <- max(fish_data$detection_date)
  all_months <- seq(floor_date(first_date, "month"), floor_date(last_date, "month"), by = "1 month")
  
  station_months <- fish_data %>%
    group_by(station, detection_month) %>%
    summarise(dates = list(unique(detection_date)), .groups = "drop")

  exclusion_month_A <- NULL
  for (i in 1:nrow(station_months)) {
    dates <- station_months$dates[[i]]
    if (length(dates) < 7) next

    month_start <- floor_date(min(dates), "month")
    first_7_days <- seq(month_start, by = "1 day", length.out = 7)
    match_dates <- intersect(dates, first_7_days)

    if (length(match_dates) >= 7 && find_7_day_run(match_dates)) {
      station_id <- station_months$station[i]
      
      future_data <- fish_data %>% filter(detection_date > max(first_7_days))
      future_stations <- unique(future_data$station)
      
      if (length(future_stations[future_stations != station_id]) == 0) {
        exclusion_month_A <- month_start
        break
      }
    }
  }

  exclusion_month_B <- NULL
  if (is.null(exclusion_month_A)) {
    for (i in 2:length(all_months)) {
      prev_month <- all_months[i - 1]
      current_month <- all_months[i]

     if (!any(floor_date(fish_data$detection_date, "month") == prev_month)) next
      prev_month_data <- fish_data %>%
        filter(floor_date(detection_date, "month") == prev_month)

      station_groups <- prev_month_data %>%
        group_by(station) %>%
        summarise(dates = list(unique(detection_date)), .groups = "drop")

      for (j in 1:nrow(station_groups)) {
        dates <- station_groups$dates[[j]]
        if (length(dates) < 7) next

        if (find_7_day_run(dates)) {
          station_id <- station_groups$station[j]
          future_data <- fish_data %>% filter(detection_date > max(dates))
          future_stations <- unique(future_data$station)

          if (length(future_stations[future_stations != station_id]) == 0) {
            exclusion_month_B <- current_month
            break
          }
        }
      }

      if (!is.null(exclusion_month_B)) break
    }
  }


  exclusion_month <- min(c(exclusion_month_A, exclusion_month_B), na.rm = TRUE)
  if (is.infinite(exclusion_month)) exclusion_month <- NA


  valid_months <- if (is.na(exclusion_month)) {
    all_months
  } else {
    all_months[all_months < exclusion_month]
  }

  tibble(FullId = fish_id, valid_month = valid_months)
})

monthly_summary_smbuff <- fish_valid_month_smbuff %>%
  group_by(valid_month) %>%
  summarise(unique_fish = n_distinct(FullId), .groups = "drop") %>%
  arrange(valid_month)


month_summary_smbuff <- fish_valid_month_smbuff %>%
  mutate(month_name = format(valid_month, "%B")) %>%  
  group_by(month_name) %>%
  summarise(total_fish_months = n(), .groups = "drop") %>%
  arrange(match(month_name, month.name))  




fish_valid_month_smbuff <- fish_valid_month_smbuff  %>%
  mutate(year=year(valid_month))  %>%
  mutate(new_id=paste0(FullId, "-", year)) %>%
  mutate(Month=month(valid_month))


```

```{r sm buffalo monthly crossing}
smbuffvps <- subset(combined_vps, Species=="smallmouth buffalo"| Species == "Smallmouth Buffalo")
smbuffvps<-subset(smbuffvps, Year!=2021)
smbuff_pop<-month_summary_smbuff$total_fish_months

smbuffvps_sf <- st_as_sf(smbuffvps, coords = c("LON", "LAT"), crs = st_crs(zones))


smbuff_joined <- st_join(smbuffvps_sf, zones)

smbuff_joined <- smbuff_joined %>%
  filter(Zone != 2 & Zone != 4)

smbuff_joined<-as.data.frame(smbuff_joined)

smbuffvps_counts <- smbuff_joined %>%
  group_by(Month, new_id) %>%
  summarise(
    count = sum(Zone != lag(Zone), na.rm = TRUE),
    downstream_crossing = sum(Zone == 1 & lag(Zone) == 3, na.rm = TRUE),
    upstream_crossing = sum(Zone == 3 & lag(Zone) == 1, na.rm = TRUE)
  )


smbuff_completed <- fish_valid_month_smbuff%>%
  left_join(smbuffvps_counts, by = c("new_id", "Month")) %>%
  mutate(
    count = replace_na(count, 0),
    downstream_crossing = replace_na(downstream_crossing, 0),
    upstream_crossing = replace_na(upstream_crossing, 0)
  )
smbuff_completed$number_crossed <- ifelse(smbuff_completed$count > 0, 1, 0)
smbuff_completed$number_notcrossed <- ifelse(smbuff_completed$count == 0, 1, 0)

  
smbuffvps_crossing <- smbuff_completed %>%
  group_by(Month) %>%
  summarise(
    mean_upstream = mean(upstream_crossing, na.rm = TRUE),
    mean_downstream = mean(downstream_crossing, na.rm = TRUE),
    var_upstream = var(upstream_crossing, na.rm = TRUE),
    var_downstream = var(downstream_crossing, na.rm = TRUE),
    )

```

# Freshwater Drum Data
```{r drum supplemental figure}
drumall <-subset(det23, common_name_e=="freshwater drum")
drumall <-subset(drumall, Year!=2021)
drum_det <- drumall %>%
  group_by(glatos_array2, Month) %>%
  summarize(unique_animal_count = n_distinct(animal_id))

drum_det$Location <- factor(drum_det$glatos_array2, levels = location_order)

ggplot(drum_det, aes(x = factor(Month, levels = 1:12, labels = substr(month.abb, 1, 3)), y = Location, size = unique_animal_count)) +
  geom_point() +
  scale_size_continuous(name = "Animal Count", breaks = seq(min(drum_det$unique_animal_count), max(drum_det$unique_animal_count), by = 2)) +
  labs(x = "Month", y = "Location", title = "Unique Detections of Freshwater Drum per Month") +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),  
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(size = 20, hjust = 0.9),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)    
  ) +
  scale_x_discrete(labels = function(labels) str_wrap(labels, width = 20)) +
  scale_y_discrete(labels = function(labels) str_wrap(labels, width = 20))  
```

```{r filtering drum monthly}
drumall_monthly <- drumall %>%
  mutate(
    detection_date = as.Date(detection_timestamp_utc),
    detection_month = floor_date(detection_date, "month")
  )

fish_list_drum <- drumall_monthly %>%
  arrange(FullId, detection_date) %>%
  group_split(FullId)


fish_valid_month_drum <- map_dfr(fish_list_drum, function(fish_data) {
  fish_id <- unique(fish_data$FullId)
  
 
  first_date <- min(fish_data$detection_date)
  last_date <- max(fish_data$detection_date)
  all_months <- seq(floor_date(first_date, "month"), floor_date(last_date, "month"), by = "1 month")
  
  
  station_months <- fish_data %>%
    group_by(station, detection_month) %>%
    summarise(dates = list(unique(detection_date)), .groups = "drop")

  
  exclusion_month_A <- NULL
  for (i in 1:nrow(station_months)) {
    dates <- station_months$dates[[i]]
    if (length(dates) < 7) next

    month_start <- floor_date(min(dates), "month")
    first_7_days <- seq(month_start, by = "1 day", length.out = 7)
    match_dates <- intersect(dates, first_7_days)

    if (length(match_dates) >= 7 && find_7_day_run(match_dates)) {
      station_id <- station_months$station[i]
      
     
      future_data <- fish_data %>% filter(detection_date > max(first_7_days))
      future_stations <- unique(future_data$station)
      
      if (length(future_stations[future_stations != station_id]) == 0) {
        exclusion_month_A <- month_start
        break
      }
    }
  }

 
  exclusion_month_B <- NULL
  if (is.null(exclusion_month_A)) {
    for (i in 2:length(all_months)) {
      prev_month <- all_months[i - 1]
      current_month <- all_months[i]

     if (!any(floor_date(fish_data$detection_date, "month") == prev_month)) next
      prev_month_data <- fish_data %>%
        filter(floor_date(detection_date, "month") == prev_month)

      station_groups <- prev_month_data %>%
        group_by(station) %>%
        summarise(dates = list(unique(detection_date)), .groups = "drop")

      for (j in 1:nrow(station_groups)) {
        dates <- station_groups$dates[[j]]
        if (length(dates) < 7) next

        if (find_7_day_run(dates)) {
          station_id <- station_groups$station[j]
          future_data <- fish_data %>% filter(detection_date > max(dates))
          future_stations <- unique(future_data$station)

          if (length(future_stations[future_stations != station_id]) == 0) {
            exclusion_month_B <- current_month
            break
          }
        }
      }

      if (!is.null(exclusion_month_B)) break
    }
  }


  exclusion_month <- min(c(exclusion_month_A, exclusion_month_B), na.rm = TRUE)
  if (is.infinite(exclusion_month)) exclusion_month <- NA

  
  valid_months <- if (is.na(exclusion_month)) {
    all_months
  } else {
    all_months[all_months < exclusion_month]
  }

  tibble(FullId = fish_id, valid_month = valid_months)
})

monthly_summary_drum <- fish_valid_month_drum %>%
  group_by(valid_month) %>%
  summarise(unique_fish = n_distinct(FullId), .groups = "drop") %>%
  arrange(valid_month)

print(monthly_summary)

month_summary_drum <- fish_valid_month_drum %>%
  mutate(month_name = format(valid_month, "%B")) %>%  
  group_by(month_name) %>%
  summarise(total_fish_months = n(), .groups = "drop") %>%
  arrange(match(month_name, month.name))  




fish_valid_month_drum <- fish_valid_month_drum  %>%
  mutate(year=year(valid_month))  %>%
  mutate(new_id=paste0(FullId, "-", year)) %>%
  mutate(Month=month(valid_month))

```

```{r drum monthly crossing}
drumvps <-subset(combined_vps, Species=="freshwater drum"|Species =="Freshwater Drum")
drumvps <- subset(drumvps, Year > 2021)
drum_pop<-month_summary_drum$total_fish_months
drumvps_sf <- st_as_sf(drumvps, coords = c("LON", "LAT"), crs = st_crs(zones))

drum_joined <- st_join(drumvps_sf, zones)

drum_joined <- drum_joined %>%
  filter(Zone != 2 & Zone != 4)

drum_joined<-as.data.frame(drum_joined)

drumvps_counts <- drum_joined %>%
  group_by(Month, new_id) %>%
    summarise(
    count = sum(Zone != lag(Zone), na.rm = TRUE),
    downstream_crossing = sum(Zone == 1 & lag(Zone) == 3, na.rm = TRUE),
    upstream_crossing = sum(Zone == 3 & lag(Zone) == 1, na.rm = TRUE)
  )


drum_completed <- fish_valid_month_drum%>%

  left_join(drumvps_counts, by = c("new_id", "Month")) %>%
  mutate(
    count = replace_na(count, 0),
    downstream_crossing = replace_na(downstream_crossing, 0),
    upstream_crossing = replace_na(upstream_crossing, 0)
  )

drum_completed$number_crossed <- ifelse(drum_completed$count > 0, 1, 0)
drum_completed$number_notcrossed <- ifelse(drum_completed$count == 0, 1, 0)

  
drumvps_crossing <- drum_completed %>%
  group_by(Month) %>%
  summarise(
    mean_upstream = mean(upstream_crossing, na.rm = TRUE),
    mean_downstream = mean(downstream_crossing, na.rm = TRUE),
    var_upstream = var(upstream_crossing, na.rm = TRUE),
    var_downstream = var(downstream_crossing, na.rm = TRUE),
    )

drumvps_crossing <- drumvps_crossing %>%
  filter(Month %in% 4:11)
print(drumvps_crossing)
```

# Shorthead Redhorse Data
```{r shorthead supplemental figure}
shortheadall <- subset(det23, common_name_e=="shorthead redhorse"|common_name_e =="Shorthead redhorse"|common_name_e =="Shorthead Redhorse")
shortheadall<-subset(shortheadall, animal_id != "x_062" &  animal_id != "x_110" & animal_id != "x_116" & animal_id != "x_128" & animal_id != "x_181") #fish that did not survive past 10 days
shortheadall<-shortheadall%>%
  filter(!(animal_id == "x_129" & Month >= 6 & Month <= 12))%>%
  filter(!(animal_id == "x_129" & Year == 2023)) #fish that died / expelled tag

shorthead_det <- shortheadall %>%
  group_by(glatos_array2, Month) %>%
  summarize(unique_animal_count = n_distinct(animal_id))

shorthead_det$Location <- factor(shorthead_det$glatos_array2, levels = location_order)

ggplot(shorthead_det, aes(x = factor(Month, levels = 1:12, labels = substr(month.abb, 1, 3)), y = Location, size = unique_animal_count)) +
  geom_point() +
  scale_size_continuous(name = "Animal Count", breaks = seq(min(shorthead_det$unique_animal_count), max(shorthead_det$unique_animal_count), by = 1)) +
  labs(x = "Month", y = "Location", title = "Unique Detections of Shorthead Redhorse per Month")+
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),  
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(size = 20, hjust = 0.9),
    legend.title = element_text(size = 18),  
    legend.text = element_text(size = 16)  
  ) 
```

```{r filtering shorthead monthly}
shortheadall_monthly <- shortheadall %>%
  mutate(
    detection_date = as.Date(detection_timestamp_utc),
    detection_month = floor_date(detection_date, "month")
  )

fish_list_shorthead <- shortheadall_monthly %>%
  arrange(FullId, detection_date) %>%
  group_split(FullId)


fish_valid_month_shorthead <- map_dfr(fish_list_shorthead, function(fish_data) {
  fish_id <- unique(fish_data$FullId)
  

  first_date <- min(fish_data$detection_date)
  last_date <- max(fish_data$detection_date)
  all_months <- seq(floor_date(first_date, "month"), floor_date(last_date, "month"), by = "1 month")
  

  station_months <- fish_data %>%
    group_by(station, detection_month) %>%
    summarise(dates = list(unique(detection_date)), .groups = "drop")

  exclusion_month_A <- NULL
  for (i in 1:nrow(station_months)) {
    dates <- station_months$dates[[i]]
    if (length(dates) < 7) next

    month_start <- floor_date(min(dates), "month")
    first_7_days <- seq(month_start, by = "1 day", length.out = 7)
    match_dates <- intersect(dates, first_7_days)

    if (length(match_dates) >= 7 && find_7_day_run(match_dates)) {
      station_id <- station_months$station[i]
      
      # Check for future detections on other stations
      future_data <- fish_data %>% filter(detection_date > max(first_7_days))
      future_stations <- unique(future_data$station)
      
      if (length(future_stations[future_stations != station_id]) == 0) {
        exclusion_month_A <- month_start
        break
      }
    }
  }


  exclusion_month_B <- NULL
  if (is.null(exclusion_month_A)) {
    for (i in 2:length(all_months)) {
      prev_month <- all_months[i - 1]
      current_month <- all_months[i]

     if (!any(floor_date(fish_data$detection_date, "month") == prev_month)) next
      prev_month_data <- fish_data %>%
        filter(floor_date(detection_date, "month") == prev_month)

      station_groups <- prev_month_data %>%
        group_by(station) %>%
        summarise(dates = list(unique(detection_date)), .groups = "drop")

      for (j in 1:nrow(station_groups)) {
        dates <- station_groups$dates[[j]]
        if (length(dates) < 7) next

        if (find_7_day_run(dates)) {
          station_id <- station_groups$station[j]
          future_data <- fish_data %>% filter(detection_date > max(dates))
          future_stations <- unique(future_data$station)

          if (length(future_stations[future_stations != station_id]) == 0) {
            exclusion_month_B <- current_month
            break
          }
        }
      }

      if (!is.null(exclusion_month_B)) break
    }
  }

  # Decide final exclusion month
  exclusion_month <- min(c(exclusion_month_A, exclusion_month_B), na.rm = TRUE)
  if (is.infinite(exclusion_month)) exclusion_month <- NA

  # Valid months = before exclusion
  valid_months <- if (is.na(exclusion_month)) {
    all_months
  } else {
    all_months[all_months < exclusion_month]
  }

  tibble(FullId = fish_id, valid_month = valid_months)
})

monthly_summary_shorthead <- fish_valid_month_shorthead %>%
  group_by(valid_month) %>%
  summarise(unique_fish = n_distinct(FullId), .groups = "drop") %>%
  arrange(valid_month)

print(monthly_summary)

month_summary_shorthead <- fish_valid_month_shorthead %>%
  mutate(month_name = format(valid_month, "%B")) %>%  # Extract month name
  group_by(month_name) %>%
  summarise(total_fish_months = n(), .groups = "drop") %>%
  arrange(match(month_name, month.name))  # Correct calendar order




fish_valid_month_shorthead <- fish_valid_month_shorthead  %>%
  mutate(year=year(valid_month))  %>%
  mutate(new_id=paste0(FullId, "-", year)) %>%
  mutate(Month=month(valid_month))

```

```{r shorthead monthly crossing}
shortheadvps <- subset(combined_vps, Species=="shorthead redhorse"|Species =="Shorthead redhorse"|Species =="Shorthead Redhorse")
shorthead<-subset(shortheadvps, animal_id != "x_062" &  animal_id != "x_110" & animal_id != "x_116" & animal_id != "x_128" & animal_id != "x_181")
shortheadvps<-shortheadvps%>%
  filter(!(animal_id == "x_129" & Month >= 6 & Month <= 12)) #these are all individuals with less than 10 days of detection or considered dead/expelled tag 

shortheadvps_sf <- st_as_sf(shortheadvps, coords = c("LON", "LAT"), crs = st_crs(zones))

# Perform spatial join
shorthead_joined <- st_join(shortheadvps_sf, zones)

shorthead_joined <- shorthead_joined %>%
  filter(Zone != 2 & Zone != 4)

shorthead_joined<-as.data.frame(shorthead_joined)


shortheadvps_counts <- shorthead_joined %>%
  group_by(Month, new_id) %>%
  summarise(
    count = sum(Zone != lag(Zone), na.rm = TRUE),
    downstream_crossing = sum(Zone == 1 & lag(Zone) == 3, na.rm = TRUE),
    upstream_crossing = sum(Zone == 3 & lag(Zone) == 1, na.rm = TRUE)
  )

print(shortheadvps_counts)

shorthead_completed <- fish_valid_month_shorthead%>%

  left_join(shortheadvps_counts, by = c("new_id", "Month")) %>%
 
  mutate(
    count = replace_na(count, 0),
    downstream_crossing = replace_na(downstream_crossing, 0),
    upstream_crossing = replace_na(upstream_crossing, 0)
  )
shorthead_completed$number_crossed <- ifelse(shorthead_completed$count > 0, 1, 0)
shorthead_completed$number_notcrossed <- ifelse(shorthead_completed$count == 0, 1, 0)

  
shortheadvps_crossing <- shortheadvps_counts %>%
  group_by(Month) %>%
  summarise(
    mean_upstream = mean(upstream_crossing, na.rm = TRUE),
    mean_downstream = mean(downstream_crossing, na.rm = TRUE),
    var_upstream = var(upstream_crossing, na.rm = TRUE),
    var_downstream = var(downstream_crossing, na.rm = TRUE),
    )

```

```{r filtering shorthead weekly}
shorthead_clean <- shortheadall %>%
  dplyr::select(FullId, detection_timestamp_utc) %>%
  mutate(
    week = isoweek(detection_timestamp_utc),
    year = year(detection_timestamp_utc),
    day_of_week = wday(detection_timestamp_utc, label = TRUE, abbr = TRUE)  
  )


shorthead_clean <- shorthead_clean %>%
  arrange(FullId, detection_timestamp_utc)


shorthead_clean <- shorthead_clean %>%
  group_by(FullId, week, year) %>%
  mutate(
    day_diff = as.numeric(detection_timestamp_utc - lag(detection_timestamp_utc, default = first(detection_timestamp_utc))),
    consecutive_days = ifelse(day_diff == 1, 1, 0), 
    consecutive_count = cumsum(consecutive_days == 1)
  ) %>%
  ungroup()


shorthead_clean <- shorthead_clean %>%
  group_by(FullId, week, year) %>%
  mutate(
    consecutive_3_days = ifelse(consecutive_count >= 3, TRUE, FALSE) 
  ) %>%
  ungroup()


remove_dead_fish <- function(fish_data) {
  fish_data <- fish_data %>%
    arrange(detection_timestamp_utc) 
  
  
  fish_data <- fish_data %>%
    mutate(
     
      dead_fish_flag = ifelse(consecutive_3_days & !any(detection_timestamp_utc > max(detection_timestamp_utc)), TRUE, FALSE)
    )
  
  
  fish_data <- fish_data %>%
    mutate(
      removed_from_current_week = ifelse(dead_fish_flag & wday(detection_timestamp_utc) <= 3, TRUE, FALSE)
    )
  
  
  fish_data <- fish_data %>%
    mutate(
      removed_from_future_weeks = ifelse(dead_fish_flag & wday(detection_timestamp_utc) >= 5, TRUE, FALSE)
    )
  
  return(fish_data)
}


shorthead_clean <- shorthead_clean %>%
  group_by(FullId, week, year) %>%
  do(remove_dead_fish(.)) %>%
  ungroup()

shorthead_clean_filtered <- shorthead_clean %>%
  filter(!(removed_from_current_week | removed_from_future_weeks))


shorthead_weekly_filtered <- shorthead_clean_filtered %>%
  distinct(FullId, year, week)


shorthead_weekly_filtered <- shorthead_weekly_filtered %>%
  mutate(week_id = year * 100 + week)


shorthead_week_range_filtered <- shorthead_weekly_filtered %>%
  group_by(FullId) %>%
  summarise(
    min_week = min(week_id),
    max_week = max(week_id),
    .groups = "drop"
  )


full_weeks_per_fish_filtered_sh <- shorthead_week_range_filtered %>%
  rowwise() %>%
 mutate(weeks = list(generate_weeks(min_week, max_week))) %>%
  ungroup() %>%
  tidyr::unnest(weeks) %>%
  mutate(
    year = floor(weeks / 100),
    week = weeks %% 100
  ) %>%
  dplyr::select(FullId, year, week)

full_weeks_per_fish_filtered_sh<-full_weeks_per_fish_filtered_sh %>%
   mutate(new_id=paste0(FullId, "-", year)) %>%
  rename(Week = week)


fish_per_week_filtered_sh <- full_weeks_per_fish_filtered_sh %>%
  group_by(Week) %>%
  summarise(unique_fish = n_distinct(FullId)) %>%
  arrange(Week)

```

```{r shorthead weekly crossing}
shorthead_counts_weekly <- shorthead_joined %>%
  mutate(Week = week(DATETIME)) %>% 
  mutate(Day = day(DATETIME)) %>%
  group_by(Week, new_id) %>%
  summarise(
    count = sum(Zone != lag(Zone), na.rm = TRUE),
    downstream_crossing = sum(Zone == 1 & lag(Zone) == 3, na.rm = TRUE),
    upstream_crossing = sum(Zone == 3 & lag(Zone) == 1, na.rm = TRUE)
  )

shvps_completed <- full_weeks_per_fish_filtered_sh %>%
  left_join(shorthead_counts_weekly, by = c("new_id", "Week")) %>%
  mutate(
    count = replace_na(count, 0),
    downstream_crossing = replace_na(downstream_crossing, 0),
    upstream_crossing = replace_na(upstream_crossing, 0)
  )

shvps_completed$number_crossed <- ifelse(shvps_completed$count > 0, 1, 0)
shvps_completed$number_notcrossed <- ifelse(shvps_completed$count == 0, 1, 0)


shorthead_crossingrate_weekly <- shvps_completed %>%
  group_by(Week) %>%
  summarise(
    mean_upstream = mean(upstream_crossing, na.rm = TRUE),
    mean_downstream = mean(downstream_crossing, na.rm = TRUE),
    var_upstream = var(upstream_crossing, na.rm = TRUE),
    var_downstream = var(downstream_crossing, na.rm = TRUE),
    )
```


#Multilevel model monthly crossing probabilities
```{r monthly crossing probabilities}
channel_completed$Species <- "Channel Catfish"
wb_completed$Species <- "White Bass"
shorthead_completed$Species <- "Shorthead Redhorse"
drum_completed$Species <- "Freshwater Drum"
smbuff_completed$Species<-"Smallmouth Buffalo"


combined_countsvps <- bind_rows(channel_completed,
  wb_completed,
  shorthead_completed,
  smbuff_completed,
  drum_completed)


combined_countsvps <- combined_countsvps %>%
  filter(Month %in% 4:11) %>%

  mutate(MonthName = case_when(
    Month == 4 ~ "April",
    Month == 5 ~ "May",
    Month == 6 ~ "June",
    Month == 7 ~ "July",
    Month == 8 ~ "August",
    Month == 9 ~ "September",
    Month == 10 ~ "October",
    Month == 11 ~ "November"
  ))
combined_countsvps <- combined_countsvps %>%
  mutate(sp_mon = paste(Species, MonthName, sep = "-"))


new.crossing.probability<- glmer(cbind(number_crossed, number_notcrossed)~1 +(1|sp_mon),
      family="binomial", data=combined_countsvps)
library(lme4)
ranef(new.crossing.probability)

dotplot(ranef(new.crossing.probability, condVar=T))
coef(new.crossing.probability)

random_effects <- ranef(new.crossing.probability)$sp_mon[, 1]
sp_mon_sd <- se.ranef(new.crossing.probability)$sp_mon[, 1]


sp_month_data <- data.frame(
  x = random_effects,
  y = factor(1:length(random_effects), labels = row.names(ranef(new.crossing.probability)$sp_mon)),
  sd = sp_mon_sd
)


sp_month_data$y <- ordered(sp_month_data$y, levels = sp_month_data$y[order(sp_month_data$x)])


dotplot(y ~ x, data = sp_month_data,
        aspect = 0.8,
        xlim = c(-4, 1.2* range((sp_month_data$x - sp_month_data$sd), (sp_month_data$x + sp_month_data$sd))[2]),
        panel = function(x, y) {
          panel.xyplot(x, y, pch = 16, col = "black")
          panel.segments((sp_month_data$x - sp_month_data$sd),
                         as.numeric(y),
                         (sp_month_data$x + sp_month_data$sd),
                         as.numeric(y), lty = 1)
        },
        xlab = "Model Random Effects",
        main = "Dot Plot with Error Bars - Sp_Mon")

sp_month_data <- data.frame(
  x = random_effects,
  y = factor(1:length(random_effects), labels = row.names(ranef(new.crossing.probability)$sp_mon)),
  sd = sp_mon_sd
)


sp_month_data <- sp_month_data %>%
  mutate(
    sp_mon = as.character(y),
    Species = sub("-(.*)", "", sp_mon),
    MonthName = sub(".*-", "", sp_mon),
    MonthNum = match(MonthName, c("April", "May", "June", "July", "August", "September", "October", "November"))
  )


sp_month_data <- sp_month_data %>%
  arrange(Species, MonthNum)

# Step 4: Reverse the y-axis order
sp_month_data$y <- factor(sp_month_data$sp_mon, levels = rev(sp_month_data$sp_mon))

sp_month_data <- sp_month_data %>%
  mutate(y_num = as.numeric(y))


sp_month_data <- sp_month_data %>%
  mutate(Species_group = as.integer(as.factor(Species)))  
plot_data <- sp_month_data

pdf("sp_mon_dotplot.pdf", width = 10, height = 9)
dotplot(y ~ x, data = sp_month_data,
        aspect = 1.5,
        xlim = c(-4, 1.2 * max(sp_month_data$x + sp_month_data$sd)),
        scales = list(
          cex = 1.2
        ),
        xlab = list("Model Random Effects", cex = 1.5),
        ylab = list("Species by Month", cex = 1.5),
        par.settings = list(
          axis.text = list(cex = 1)
        ),
        panel = function(x, y, ...) {
  unique_species <- unique(plot_data$Species)
  
  for (i in seq_along(unique_species)) {
    sp <- unique_species[i]
    rows <- which(plot_data$Species == sp)
    
    if (i %% 2 == 1) {  
      panel.rect(xleft = -100, xright = 100,
                 ybottom = min(plot_data$y_num[rows]) - 0.5,
                 ytop = max(plot_data$y_num[rows]) + 0.5,
                 col = "gray95", border = NA)
    }
  }

  
  panel.xyplot(x, y, pch = 16, col = "black", cex = 1.5)
  panel.segments(plot_data$x - plot_data$sd,
                 y,
                 plot_data$x + plot_data$sd,
                 y,
                 lty = 1, lwd = 2, col = "gray40")

        }
  )
combined_sim <-sim(new.crossing.probability, 2000) 
library(rv)

pred <- rvsims(combined_sim@fixef[,1] + combined_sim@ranef$sp_mon[,,1])

pred_sum <- summary(invlogit(pred))




new_data <- data.frame(x = pred_sum[,2],
                       y = labels(pred),
                       q2.5 = pred_sum[,5],
                       q25 = pred_sum[,6],
                       q50 = pred_sum[,7],
                       q75 = pred_sum[,8],
                       q975 = pred_sum[,9])

new_data$y <- factor(new_data$y,levels = rev(sp_month_data$sp_mon))


dotplot(y ~ x, data = new_data, xlim = c(0, 100),  
        panel = function(x, y, ...) {
          panel.points(x = 100 * new_data$q50, y = as.numeric(y), pch = 16)  
          
          panel.segments(100 * new_data$q2.5, as.numeric(y),  
                         100 * new_data$q975, as.numeric(y), ...)
          panel.segments(100 * new_data$q25, as.numeric(y),   
                         100 * new_data$q75, as.numeric(y),
                         lwd = 2.5, ...)
        },
        xlab = "Probability of Crossing (%)",
        ylab = "Species by Month")

dotplot(x ~ y, data = new_data, ylim = c(0, 100),  
        panel = function(x, y, ...) {
          panel.points(y = 100 * new_data$q50, x = as.numeric(x), pch = 16)  
          
          panel.segments(as.numeric(x), 100 * new_data$q2.5,  
                         as.numeric(x), 100 * new_data$q975, ...)
          panel.segments(as.numeric(x), 100 * new_data$q25,   
                         as.numeric(x), 100 * new_data$q75,
                         lwd = 2.5, ...)
        },
        xlab = "Species by Month",
        ylab = "Probability of Crossing (%)",
        par.settings = list(
          par.xlab.text = list(cex = 1.2),  
          par.ylab.text = list(cex = 1.2), 
          par.main.text = list(cex = 1.5)
        ),
        scales = list(
          x = list(cex = 1),  
          y = list(cex = 1)   
        ))


new_data <- new_data %>%
  mutate(
    Species = sub("-(.*)", "", as.character(y)),  
    x_pos = as.numeric(y)                        
  )

plot_data <- new_data  


dotplot(100 * q50 ~ y, data = new_data,
        # aspect = "fill",
        ylim = c(0, 100),
        xlab = "Species by Month",
        ylab = "Probability of Crossing (%)",
        par.settings = list(
          par.xlab.text = list(cex = 1.2),
          par.ylab.text = list(cex = 1.2),
          par.main.text = list(cex = 1.5)
        ),
        scales = list(
          x = list(cex = 1, rot = 90),  
          y = list(cex = 1)
        ),
        panel = function(x, y, subscripts, ...) {
          
          d <- plot_data[subscripts, ]

          
          unique_species <- unique(d$Species)
          for (i in seq_along(unique_species)) {
            sp <- unique_species[i]
            rows <- which(d$Species == sp)
            xpos <- d$x_pos[rows]

            if (i %% 2 == 1) {
              panel.rect(
                xleft = min(xpos) - 0.5,
                xright = max(xpos) + 0.5,
                ybottom = -10, ytop = 110,
                col = "gray95", border = NA
              )
            }
          }

          
          panel.points(x = d$x_pos, y = 100 * d$q50,
                       pch = 16, col = "black", cex = 1.2)

         
          panel.segments(
            x0 = d$x_pos,
            x1 = d$x_pos,
            y0 = 100 * d$q2.5,
            y1 = 100 * d$q975,
            col = "gray50", lwd = 1.5
          )

          
          panel.segments(
            x0 = d$x_pos,
            x1 = d$x_pos,
            y0 = 100 * d$q25,
            y1 = 100 * d$q75,
            col = "black", lwd = 2.5
          )
        })

plot_data <- new_data %>%
  mutate(
    Species = sub("-(.*)", "", as.character(y)),
    Month = sub(".*-", "", as.character(y))
  )

selected_species <- c("Channel Catfish", "Freshwater Drum", "Smallmouth Buffalo")

plot_subset <- plot_data %>%
  filter(Species %in% selected_species)

month_order <- c("April", "May", "June", "July", "August", "September", "October", "November")

plot_subset$Month <- factor(plot_subset$Month, levels = month_order)

# Define custom colors for each species
custom_colors <- c(
  "Channel Catfish" = "turquoise3",
  "Freshwater Drum" = "sienna1",
  "Smallmouth Buffalo" = "darkorchid4"
)



plot_subset$Month <- factor(plot_subset$Month, levels = month.name)

ggplot(plot_subset, aes(x = Month, y = 100 * q50, color = Species)) +
  
  annotate(
    "rect",
    xmin = 1.5, xmax = 5.5,
    ymin = -Inf, ymax = Inf,
    fill = "gray80",
    alpha = 0.3  
  ) +

 
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(
    aes(ymin = 100 * q2.5, ymax = 100 * q975),
    position = position_dodge(width = 0.6),
    width = 0.2,
    linewidth = 1
  ) +

  scale_color_manual(values = custom_colors) +
  scale_y_continuous(limits = c(0, 100)) +

  labs(
    x = "Month",
    y = "Probability of Crossing (%)"
  ) +

  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.75, 0.85),
    legend.background = element_rect(fill = "white", color = "gray80", size = 0.3),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


```


```{r weekly crossing probabilities}
wbvps_completed$Species <- "White Bass"
shvps_completed$Species <- "Shorthead Redhorse"


combined_weekly_countsvps <- bind_rows(wbvps_completed, shvps_completed)
combined_weekly_countsvps <- subset(combined_weekly_countsvps, Week<=26 & Week >=16)

combined_weekly_countsvps <- combined_weekly_countsvps %>%
  mutate(sp_week = paste(Species, Week, sep = "-"))


new.weekly.crossing.probability<- glmer(cbind(number_crossed, number_notcrossed)~1 +(1|sp_week),
      family="binomial", data=combined_weekly_countsvps)

ranef(new.weekly.crossing.probability)

dotplot(ranef(new.weekly.crossing.probability, condVar=T))
coef(new.weekly.crossing.probability)

random_effects1 <- ranef(new.weekly.crossing.probability)$sp_week[, 1]
sp_week_sd <- se.ranef(new.weekly.crossing.probability)$sp_week[, 1]


sp_week_data <- data.frame(
  x = random_effects1,
  y = factor(1:length(random_effects1), labels = row.names(ranef(new.weekly.crossing.probability)$sp_week)),
  sd = sp_week_sd
)


sp_week_data$y <- ordered(sp_week_data$y, levels = sp_week_data$y[order(sp_week_data$x)])


dotplot(y ~ x, data = sp_week_data,
        aspect = 0.8,
        xlim = c(-4, 1.2* range((sp_week_data$x - sp_week_data$sd), (sp_week_data$x + sp_week_data$sd))[2]),
        panel = function(x, y) {
          panel.xyplot(x, y, pch = 16, col = "black")
          panel.segments((sp_week_data$x - sp_week_data$sd),
                         as.numeric(y),
                         (sp_week_data$x + sp_week_data$sd),
                         as.numeric(y), lty = 1)
        },
        xlab = "Model Random Effects",
        main = "Dot Plot with Error Bars - Sp_Week")
sp_week_data <- sp_week_data %>%
  mutate(
    sp_week = as.character(y),
    Species = sub("-(.*)", "", sp_week),
    Week = as.numeric(sub(".*-", "", sp_week))  
  )


sp_week_data <- sp_week_data[order(sp_week_data$Species, sp_week_data$Week), ]
sp_week_data$y <- factor(sp_week_data$y, levels = rev(sp_week_data$y))  


sp_week_data$y_num <- as.numeric(sp_week_data$y)


plot_data <- sp_week_data


dotplot(y ~ x, data = sp_week_data,
        aspect = 1.5,
        xlim = c(-4, 1.2 * max(sp_week_data$x + sp_week_data$sd)),
        scales = list(cex = 1.2),
        xlab = list("Model Random Effects", cex = 1.5),
        ylab = list("Species by Week", cex = 1.5),
        par.settings = list(axis.text = list(cex = 1)),
        panel = function(x, y, ...) {
          unique_species <- unique(plot_data$Species)

          for (i in seq_along(unique_species)) {
            sp <- unique_species[i]
            rows <- which(plot_data$Species == sp)

            if (sp == "Shorthead Redhorse") {  # Only shade this species
              panel.rect(
                xleft = -100, xright = 100,
                ybottom = min(plot_data$y_num[rows]) - 0.5,
                ytop = max(plot_data$y_num[rows]) + 0.5,
                col = "gray95", border = NA
              )
            }
          }

          
          panel.segments(
            plot_data$x - plot_data$sd, y,
            plot_data$x + plot_data$sd, y,
            lty = 1, lwd = 2, col = "gray40"
          )

          
          panel.xyplot(x, y, pch = 16, col = "black", cex = 1.0)
        })

sp_week_data <- data.frame(
  x = random_effects1,
  y = factor(1:length(random_effects1), labels = row.names(ranef(new.weekly.crossing.probability)$sp_week)),
  sd = sp_week_sd
)


sp_week_data <- sp_week_data %>%
  mutate(
    sp_week = as.character(y),
    Species = sub("-(.*)", "", sp_week),
    Week = sub(".*-", "", sp_week)
  )


sp_week_data <- sp_week_data %>%
  arrange(Species, Week)


sp_week_data$y <- factor(sp_week_data$sp_week, levels = rev(sp_week_data$sp_week))

sp_week_data <- sp_week_data %>%
  mutate(y_num = as.numeric(y))


sp_week_data <- sp_week_data %>%
  mutate(Species_group = as.integer(as.factor(Species)))  

plot_data1 <- sp_week_data

pdf("sp_week_dotplot.pdf", width = 10, height = 9)
dotplot(y ~ x, data = sp_week_data,
        aspect = 1.5,
        xlim = c(-4, 1.2 * max(sp_week_data$x + sp_week_data$sd)),
        scales = list(
          cex = 1.2
        ),
        xlab = list("Model Random Effects", cex = 1.5),
        ylab = list("Species by Week", cex = 1.5),
        par.settings = list(
          axis.text = list(cex = 1)
        ),
        panel = function(x, y, ...) {
  unique_species <- unique(plot_data1$Species)
  
  for (i in seq_along(unique_species)) {
    sp <- unique_species[i]
    rows <- which(plot_data1$Species == sp)
    
    if (i %% 2 == 1) {  
      panel.rect(xleft = -100, xright = 100,
                 ybottom = min(plot_data$y_num[rows]) - 0.5,
                 ytop = max(plot_data1$y_num[rows]) + 0.5,
                 col = "gray95", border = NA)
    }
  }

  
  panel.xyplot(x, y, pch = 16, col = "black", cex = 1.5)
  panel.segments(plot_data1$x - plot_data1$sd,
                 y,
                 plot_data1$x + plot_data1$sd,
                 y,
                 lty = 1, lwd = 2, col = "gray40")

        }
  )



combined_weekly_sim <-sim(new.weekly.crossing.probability, 2000) ## 2000 random samples


pred2 <- rvsims(combined_weekly_sim@fixef[,1] + combined_weekly_sim@ranef$sp_week[,,1])

pred_sum2 <- summary(invlogit(pred2))

new_data2 <- data.frame(x = pred_sum2[,2],
                       y = labels(pred2),
                       q2.5 = pred_sum2[,5],
                       q25 = pred_sum2[,6],
                       q50 = pred_sum2[,7],
                       q75 = pred_sum2[,8],
                       q975 = pred_sum2[,9])


dotplot(y ~ x, data = new_data2, xlim = c(0, 100),  # Adjusted xlim for percentage scale
        panel = function(x, y, ...) {
          panel.points(x = 100 * new_data2$q50, y = as.numeric(y), pch = 16)  # Multiply x by 100
          
          panel.segments(100 * new_data2$q2.5, as.numeric(y),  # Multiply x by 100
                         100 * new_data2$q975, as.numeric(y), ...)
          panel.segments(100 * new_data2$q25, as.numeric(y),   # Multiply x by 100
                         100 * new_data2$q75, as.numeric(y),
                         lwd = 2.5, ...)
        },
        xlab = "Probability of Crossing (%)",
        ylab = "Species by Week")

new_data2 <- new_data2 %>%
  mutate(
    Species = sub("-(.*)", "", as.character(y)))

plot_data2 <- new_data2  


plot_data2 <- new_data2 %>%
  mutate(
    Species = sub("-(.*)", "", as.character(y)),
    Week = sub(".*-", "", as.character(y))
  )

selected_species2 <- c("Shorthead Redhorse", "White Bass")

plot_subset2 <- plot_data2 %>%
  filter(Species %in% selected_species2)


custom_colors2 <- c(
  "Shorthead Redhorse" = "darkorchid4",
  "White Bass" = "turquoise3"
)



new_labels <- c(
  "16 - 4/16–4/22", "17 - 4/23–4/29", "18 - 4/30–5/6", 
  "19 - 5/7–5/13", "20 - 5/14–5/20", "21 - 5/21–5/27", 
  "22 - 5/28–6/3", "23 - 6/4–6/10", "24 - 6/11–6/17", 
  "25 - 6/18–6/24", "26 - 6/25–7/1"
)

plot_subset2$Week <- as.numeric(as.character(plot_subset2$Week))


ggplot(plot_subset2, aes(x = Week, y = 100 * q50, color = Species)) +
  geom_errorbar(
    aes(ymin = 100 * q2.5, ymax = 100 * q975),
    position = position_dodge(width = 0.6),
    width = 0.2,
    linewidth = 1
  ) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  
  
  scale_color_manual(values = custom_colors2) +
  
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(
    breaks = 16:26,
    labels = c("4/16–4/22", "4/23–4/29", "4/30–5/6", 
               "5/7–5/13", "5/14–5/20", "5/21–5/27", 
               "5/28–6/3", "6/4–6/10", "6/11–6/17", 
               "6/18–6/24", "6/25–7/1")
  ) +
  
  labs(
    x = "Week",
    y = "Probability of Crossing (%)"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.75, 0.85),
    legend.background = element_rect(fill = "white", color = "gray80", size = 0.3),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

```
