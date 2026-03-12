# ==========================================
# PART 1: DATA LOADING & SETUP
# ==========================================


set.seed(4)
# 1. Load required libraries
library(tidyverse)
library(readr)
library(tidyr)
library(sf)     
library(writexl)

# 2. Read the canonical datasets
population <- read_csv("population.csv", show_col_types = FALSE)
cases      <- read_csv("cases.csv", show_col_types = FALSE)
movements  <- read_csv("movement.csv", show_col_types = FALSE)
activity   <- read_csv("activity.csv", show_col_types = FALSE)
prev_culls <- read_csv("prev_culls.csv", show_col_types = FALSE)

# 3. Define the simulation timeline
start_date <- as.Date("2025-09-01") 
current_date <- as.Date("2026-01-15") 
end_date <- current_date + 28         # +4 weeks
sim_days <- seq(start_date, end_date, by = "day")
num_days <- length(sim_days)
num_farms <- nrow(population)

# 4. Clean dates in historical files to prevent errors
cases <- cases %>%
  mutate(
    date_suspicious = as.Date(date_suspicious),
    date_confirmed = as.Date(date_confirmed),
    cull_start = as.Date(cull_start)
  )

prev_culls <- prev_culls %>% mutate(cull_start = as.Date(cull_start))
activity_clean <- activity %>% mutate(date_end = replace_na(date_end, end_date))


# ==========================================
# PART 2: SPATIAL RISK & DISTANCE MATRIX
# ==========================================

cat("\n--- ADDING GEOGRAPHIC & FARM RISK FACTORS ---\n")

# 1. Convert population dataframe into a spatial object (EPSG:32626)
pop_sf <- st_as_sf(population, coords = c("x", "y"), crs = 32626)

# 2. Load the High-Risk Zone (HRZ) polygon
hrz_polygon <- st_read("hrz_32626.geojson", quiet = TRUE)

# 3. Flag farms inside the HRZ
in_hrz_matrix <- st_intersects(pop_sf, hrz_polygon, sparse = FALSE)
population$in_hrz <- in_hrz_matrix[, 1]

# --- WATER HABITAT PROXIMITY ---
# 4. Load CLC data and isolate wetlands (4xx) and water bodies (5xx)
clc_data <- st_read("clc_32626.geojson", quiet = TRUE)
water_areas <- clc_data %>% 
  filter(substr(as.character(CODE_12), 1, 1) %in% c("4", "5"))

cat("Checking proximity to water habitats (this might take a few seconds)...\n")
# Check if each farm is within 1000 meters of any water polygon
# lengths() > 0 turns the output list into TRUE/FALSE
population$near_water <- lengths(st_is_within_distance(pop_sf, water_areas, dist = 1000)) > 0

cat("Number of farms within 1km of water:", sum(population$near_water), "\n")


# 5. Define Custom Background Risk (Wild Bird Spillovers)
population <- population %>%
  mutate(
    is_outdoor = if_else(production %in% c("broiler_2", "organic"), TRUE, FALSE),
    
    # We have 3 multipliers: Outdoor, HRZ, and Water Proximity.
    bg_risk = case_when(
      # The worst-case scenario: Outdoor, in HRZ, AND near water
      in_hrz == TRUE & near_water == TRUE & is_outdoor == TRUE ~ 0.00001 * 500,
      
      # Outdoor and in HRZ (but no water nearby)
      in_hrz == TRUE & near_water == FALSE & is_outdoor == TRUE ~ 0.00001 * 100,
      
      # Outdoor and near water (but not in HRZ)
      in_hrz == FALSE & near_water == TRUE & is_outdoor == TRUE ~ 0.00001 * 100,
      
      # Just an outdoor farm
      in_hrz == FALSE & near_water == FALSE & is_outdoor == TRUE ~ 0.00001 * 10,
      
      # Indoor farms (mostly safe from wild birds)
      TRUE ~ 0.00001 
    )
  )

farm_bg_risk <- population$bg_risk
names(farm_bg_risk) <- population$farm_id

cat("Farm risks calculated. Max daily background risk:", max(farm_bg_risk), "\n")

# 6. Build the Distance Matrix (in km) for spatial spread
coords <- population %>% select(x, y)
dist_matrix <- as.matrix(dist(coords)) / 1000
rownames(dist_matrix) <- population$farm_id
colnames(dist_matrix) <- population$farm_id

# ==========================================
# PART 3: INITIALIZE THE STATE MATRIX
# ==========================================
state_matrix <- matrix(0, nrow = num_farms, ncol = num_days)
rownames(state_matrix) <- population$farm_id
colnames(state_matrix) <- as.character(sim_days)

cat("\nSetup Complete! The environment is ready for the SEIR model.\n")


# to add ##
# implement how many cases are primary (from wildlife) e secondary (from farm to farm) #



# =================================================================
# PART 2: THE 4-WEEK PREDICTION (FULLY POLICY & BIOLOGY MANAGED)
# =================================================================

# --- 1. RESET THE SEIR MATRIX TO HISTORICAL REALITY ---
state_matrix[,] <- 0

for (i in 1:nrow(activity_clean)) {
  farm <- as.character(activity_clean$farm_id[i])
  start_d <- activity_clean$date_start[i]
  end_d <- activity_clean$date_end[i]
  active_days <- as.character(sim_days[sim_days >= start_d & sim_days <= end_d])
  if (length(active_days) > 0 && farm %in% rownames(state_matrix)) state_matrix[farm, active_days] <- 1
}

for (i in 1:nrow(cases)) {
  this_farm <- as.character(cases$farm_id[i])
  date_sick <- cases$date_suspicious[i]
  if (is.na(date_sick)) date_sick <- cases$date_confirmed[i]
  date_culled <- cases$cull_start[i]
  if (is.na(date_culled)) date_culled <- end_date + 1 
  if (is.na(date_sick)) next 
  
  infected_days <- as.character(sim_days[sim_days >= date_sick & sim_days < date_culled])
  if (length(infected_days) > 0 && this_farm %in% rownames(state_matrix)) state_matrix[this_farm, infected_days] <- 3 
  
  culled_days <- as.character(sim_days[sim_days >= date_culled])
  if (length(culled_days) > 0 && this_farm %in% rownames(state_matrix)) state_matrix[this_farm, culled_days] <- 4 
}

for (i in 1:nrow(prev_culls)) {
  this_farm <- as.character(prev_culls$farm_id[i])
  date_culled <- prev_culls$cull_start[i]
  if (!is.na(date_culled)) {
    culled_days <- as.character(sim_days[sim_days >= date_culled])
    if (length(culled_days) > 0 && this_farm %in% rownames(state_matrix)) state_matrix[this_farm, culled_days] <- 4 
  }
}

# --- 2. CAPACITY, VOLUME & POLICY PREP ---
median_cap <- median(population$capacity, na.rm = TRUE)
cap_mult <- population$capacity / median_cap
names(cap_mult) <- as.character(population$farm_id)
median_vol <- median(movements$volume, na.rm = TRUE)

# -----------------------------------------------------------------
# NEW: SPECIES-SPECIFIC SHEDDING AMPLIFICATION (Based on Wu et al.)
# Ducks transmit at 4.1, Chickens at 1.15. Ducks are ~3.5x more infectious.
# -----------------------------------------------------------------
species_mult <- ifelse(population$species == "duck", 3.5, 1.0)
names(species_mult) <- as.character(population$farm_id)


# -----------------------------------------------------------------
# BIOLOGICALLY DELAYS (Based on Wu et al.2026)
# -----------------------------------------------------------------
# Chickens show clinical signs rapidly (mean ~2.1 days)
chicken_delays <- 2 

# Ducks can be asymptomatic "silent spreaders" for much longer
duck_delays <- 8 

# Latency 
latency_delay <- 1  

# Track detection days to manage 10km movement bans
detection_day <- rep(NA, num_farms)
names(detection_day) <- as.character(population$farm_id)

for (i in 1:nrow(cases)) {
  farm_id <- as.character(cases$farm_id[i])
  det_date <- cases$date_suspicious[i]
  if (is.na(det_date)) det_date <- cases$date_confirmed[i]
  if (!is.na(det_date)) {
    day_idx <- which(sim_days == det_date)
    if(length(day_idx) > 0) detection_day[farm_id] <- day_idx
  }
}

# Identify Stage 1 Broilers in HRZ for Pre-Movement Testing
hrz_stage1_farms <- population %>%
  filter(in_hrz == TRUE & production == "broiler_1") %>%
  pull(farm_id) %>% as.character()

# --- 3. CULLING QUEUE PREP ---
target_cull_day <- rep(NA, num_farms)
names(target_cull_day) <- as.character(population$farm_id)
cull_type <- rep(NA, num_farms) 
names(cull_type) <- as.character(population$farm_id)

cull_capacity_per_day <- 8 # Government limit based on the narrative

# --- 4. PARAMETERS ---
beta_spatial    <- 0.007    
spatial_scale   <- 1.0        
beta_network    <- 0.04     

start_pred_idx <- which(sim_days == current_date) + 1
daily_new_cases <- rep(0, num_days)
names(daily_new_cases) <- as.character(sim_days)
predicted_sick_farms <- c() 

cat("\n--- RUNNING POLICY-MANAGED SEIR (WITH BIOLOGICAL FEEDBACK) ---\n")

# --- 5. 4-WEEK PREDICTION LOOP ---
for (t in start_pred_idx:num_days) {
  
  today_date <- sim_days[t]
  S_farms <- names(which(state_matrix[, t] == 1))
  I_farms <- names(which(state_matrix[, t] == 3)) 
  
  if (length(S_farms) == 0) next
  
  spatial_risk <- rep(0, length(S_farms))
  network_risk <- rep(0, length(S_farms))
  names(spatial_risk) <- S_farms
  names(network_risk) <- S_farms
  
  # ENFORCE 10km MOVEMENT BAN
  active_outbreaks <- names(which(detection_day <= t & detection_day >= (t - 28)))
  banned_farms <- c()
  if (length(active_outbreaks) > 0) {
    in_zone <- dist_matrix[, active_outbreaks, drop = FALSE] <= 10
    banned_farms <- rownames(in_zone)[rowSums(in_zone) > 0]
  }
  
  if (length(I_farms) > 0) {
    # Spatial Decay
    dist_I_to_S <- dist_matrix[S_farms, I_farms, drop = FALSE]
    decay_matrix <- exp(-dist_I_to_S / spatial_scale)
    
    # -------------------------------------------------------------
    # APPLY DUCK SHEDDING AMPLIFIER TO INFECTIOUS FARMS
    # -------------------------------------------------------------
    I_cap_mult <- cap_mult[I_farms] * species_mult[I_farms]
    
    scaled_decay <- sweep(decay_matrix, 2, I_cap_mult, "*")
    spatial_risk <- beta_spatial * rowSums(scaled_decay) * cap_mult[S_farms]
    
    today_moves <- movements %>% filter(date == today_date & source_farm %in% I_farms)
    
    # Apply 10km Ban
    if (length(banned_farms) > 0) {
      today_moves <- today_moves %>% filter(!(source_farm %in% banned_farms))
    }
    # Apply HRZ Pre-Movement PCR Testing
    if (nrow(today_moves) > 0) {
      today_moves <- today_moves %>% filter(!(source_farm %in% hrz_stage1_farms))
    }
    
    if (nrow(today_moves) > 0) {
      vol_by_dest <- tapply(today_moves$volume, today_moves$dest_farm, sum)
      valid_exposed <- as.character(names(vol_by_dest)[names(vol_by_dest) %in% S_farms])
      if (length(valid_exposed) > 0) {
        vol_mult <- vol_by_dest[valid_exposed] / median_vol
        network_risk[valid_exposed] <- beta_network * vol_mult
      }
    }
  }
  
  # Risk & Coin Flip
  total_lambda <- farm_bg_risk[S_farms] + spatial_risk + network_risk
  prob_infection <- 1 - exp(-total_lambda)
  new_infections <- rbinom(n = length(S_farms), size = 1, prob = prob_infection)
  
  new_I_farms <- S_farms[new_infections == 1]
  daily_new_cases[t] <- length(new_I_farms)
  predicted_sick_farms <- c(predicted_sick_farms, new_I_farms)
  
  # SEIR STATE UPDATES & QUEUE ASSIGNMENT
  if (length(new_I_farms) > 0) {
    for (farm in new_I_farms) {
      
      # BIOLOGICALLY  DELAYS (Chickens vs Ducks)
      this_species <- population$species[population$farm_id == farm]
      if (this_species == "chicken") {
        this_cull_delay <- chicken_delays
      } else {
        this_cull_delay <- duck_delays
      }
      
      exp_end_idx <- min(t + latency_delay - 1, num_days)
      sched_cull_idx <- min(t + latency_delay + this_cull_delay - 1, num_days)
      
      det_day <- max(t, sched_cull_idx - 1)
      detection_day[farm] <- det_day
      
      if (t <= exp_end_idx) state_matrix[farm, t:exp_end_idx] <- 2
      # Stuck in State 3 until the queue clears them!
      if ((exp_end_idx + 1) <= num_days) state_matrix[farm, (exp_end_idx + 1):num_days] <- 3
      
      # Add to Reactive Queue
      target_cull_day[farm] <- sched_cull_idx
      cull_type[farm] <- "reactive"
      
      # Add 1km neighbors to Preventive Queue
      neighbors_1km <- names(which(dist_matrix[farm, ] <= 1))
      neighbors_1km <- setdiff(neighbors_1km, farm)
      
      for (n_farm in neighbors_1km) {
        if (is.na(cull_type[n_farm]) || cull_type[n_farm] != "reactive") {
          target_cull_day[n_farm] <- sched_cull_idx + 1
          cull_type[n_farm] <- "preventive"
        }
      }
    }
  }
  
  # EXECUTE THE DAILY CULLING QUEUE
  eligible_reactive <- names(which(target_cull_day <= t & cull_type == "reactive" & state_matrix[, t] != 4))
  eligible_preventive <- names(which(target_cull_day <= t & cull_type == "preventive" & state_matrix[, t] != 4))
  
  budget <- cull_capacity_per_day
  culled_today <- c()
  
  # Prioritize Reactive Culls
  if (length(eligible_reactive) > 0) {
    sel <- eligible_reactive[1:min(length(eligible_reactive), budget)]
    culled_today <- c(culled_today, sel)
    budget <- budget - length(sel)
  }
  
  # Use remaining budget on Preventive Culls
  if (budget > 0 && length(eligible_preventive) > 0) {
    sel <- eligible_preventive[1:min(length(eligible_preventive), budget)]
    culled_today <- c(culled_today, sel)
  }
  
  # Apply State 4 to the culled farms
  if (length(culled_today) > 0) {
    for (c_farm in culled_today) {
      active_future_days <- which(state_matrix[c_farm, t:num_days] != 0)
      if (length(active_future_days) > 0) {
        actual_days <- (t:num_days)[active_future_days]
        state_matrix[c_farm, actual_days] <- 4
      }
    }
  }
}

# --- 6. VIEW THE FINAL RESULTS ---
cat("\n--- FINAL PREDICTION SNAPSHOT (2026-02-12) ---\n")
final_status <- table(state_matrix[, num_days])
cat("Total Farms Culled (4) by end of prediction:", ifelse(is.na(final_status["4"]), 0, final_status["4"]), "\n")

backlog_count <- sum(target_cull_day <= num_days & state_matrix[, num_days] != 4, na.rm = TRUE)
cat("WARNING: Farms stuck in the culling backlog on final day:", backlog_count, "\n")

total_predicted_cases <- sum(daily_new_cases[start_pred_idx:num_days])
cat("Total NEW cases in the next 4 weeks:", total_predicted_cases, "\n")

# --- 7. BREAKDOWN OF NEW CASES BY FARM TYPE ---
cat("\n--- BREAKDOWN OF NEW CASES DURING PREDICTION WINDOW ---\n")
if (length(predicted_sick_farms) > 0) {
  sick_farm_details <- population %>%
    filter(farm_id %in% as.numeric(predicted_sick_farms)) %>%
    group_by(species, production) %>%
    summarise(new_cases = n(), .groups = "drop") %>%
    arrange(desc(new_cases))
  
  print(sick_farm_details)
} else {
  cat("Wow! The epidemic naturally died out. No new cases predicted.\n")
}

# =================================================================
# PART 8: VISUALIZING THE PREDICTION (EPICURVE & MAP)
# =================================================================

cat("\n--- GENERATING PLOTS ---\n")

if (length(predicted_sick_farms) > 0) {
  
  # 1. Safely extract infection days including State 4 for same-day culls!
  inf_idx <- sapply(predicted_sick_farms, function(f) {
    states_in_window <- state_matrix[f, start_pred_idx:num_days]
    # We check for 2 (Exposed), 3 (Infectious), OR 4 (Culled)
    min(which(states_in_window %in% c(2, 3, 4))) + start_pred_idx - 1
  })
  
  # 2. Build the tracking dataframe
  new_cases_df <- data.frame(
    farm_id = as.integer(as.character(predicted_sick_farms)), 
    inf_date = as.Date(sim_days[inf_idx], origin = "1970-01-01") 
  ) %>%
    left_join(population, by = "farm_id") %>%
    filter(!is.na(inf_date)) # Failsafe to drop any persistent NAs
  
  # --- PLOT 1: LINE PLOT OF DAILY OUTBREAKS BY PRODUCTION TYPE ---
  start_pred_date <- as.Date(sim_days[start_pred_idx], origin = "1970-01-01")
  end_pred_date   <- as.Date(sim_days[num_days], origin = "1970-01-01")
  full_date_range <- seq.Date(from = start_pred_date, to = end_pred_date, by = "day")
  
  plot_data <- new_cases_df %>%
    group_by(inf_date, production) %>%
    summarise(cases = n(), .groups = "drop") %>%
    complete(inf_date = full_date_range, 
             production = unique(population$production), 
             fill = list(cases = 0))
  
  p1 <- ggplot(plot_data, aes(x = inf_date, y = cases, color = production, group = production)) +
    geom_line(linewidth = 1) +                  
    geom_point(size = 1.2, alpha = 0.5) +       
    scale_x_date(date_breaks = "1 week", date_labels = "%b %d") + 
    theme_minimal() +
    labs(
      title = "Epidemic Velocity: 4-Week HPAI Forecast",
      subtitle = "Daily new infections by production type",
      x = "Date", 
      y = "New Outbreaks (Daily)", 
      color = "Production Type"
    ) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  print(p1)
  
  # --- PLOT 2: SPATIAL SPREAD MAP BY WEEK ---
  
  # Mathematically assign each infection to a strict 7-day week window
  new_cases_df <- new_cases_df %>%
    mutate(
      days_since_start = as.numeric(inf_date - start_pred_date),
      week_num = ceiling((days_since_start + 1) / 7),
      infection_week = paste("Week", week_num)
    )
  
  p2 <- ggplot() +
    geom_point(data = population, aes(x = x, y = y), color = "grey85", size = 0.5, alpha = 0.5) +
    geom_point(data = new_cases_df, aes(x = x, y = y, color = infection_week), size = 2.5, alpha = 0.8) +
    scale_color_brewer(palette = "Set1") + 
    theme_void() +
    labs(
      title = "Geographic Spread of HPAI Over Time",
      subtitle = "Tracing the movement of the virus across the island by week",
      color = "Infection Wave"
    ) +
    theme(legend.position = "right")
  
  print(p2)
  
  cat("Plots generated successfully! Check your plots pane.\n")
  
} else {
  cat("No new cases to plot!\n")
}

# =================================================================
# PART 3: CALIBRATION WITH FIXED SPATIAL SCALE
# =================================================================
# --- 1. PREPARE THE CALIBRATION TARGET & SEEDS ---
seed_cases <- cases %>% filter(date_suspicious <= as.Date("2025-12-27") | date_confirmed <= as.Date("2025-12-27"))
target_total_cases <- nrow(cases) # Our target is 103

calib_start_date <- as.Date("2025-12-28")
calib_end_date <- as.Date("2026-01-15")
calib_start_idx <- which(sim_days == calib_start_date)
calib_end_idx <- which(sim_days == calib_end_date)

# -----------------------------------------------------------------
# EMPIRICAL PARAMETERS (Wu et al. & Ssematimba et al.)
# -----------------------------------------------------------------
chicken_delays <- 2  # Fixed at 2 days
duck_delays <- 8     # Fixed at 8 days
latency_delay <- 1  
spatial_scale <- 1.0 # Fixed spatial scale

species_mult <- ifelse(population$species == "duck", 3.5, 1.0)
names(species_mult) <- as.character(population$farm_id)

# --- 2. DEFINE THE PARAMETER GRID ---
param_grid <- expand.grid(
  beta_spatial = c(0.001, 0.003, 0.005, 0.007, 0.009), 
  beta_network = seq(0.01, 0.06, by = 0.01)
) 

param_grid$error <- NA 
param_grid$sim_cases <- NA

cat("\n--- STARTING FINAL MASTER CALIBRATION ---\n")
cat("Testing", nrow(param_grid), "combinations on the Queue-Managed SEIR model...\n\n")

# --- 3. CALIBRATION LOOP ---
for (row in 1:nrow(param_grid)) {
  
  test_b_spatial <- param_grid$beta_spatial[row]
  test_b_network <- param_grid$beta_network[row]
  
  # Reset the test matrix and trackers
  test_matrix <- matrix(0, nrow = num_farms, ncol = num_days)
  rownames(test_matrix) <- population$farm_id
  colnames(test_matrix) <- as.character(sim_days)
  
  test_detection_day <- rep(NA, num_farms)
  names(test_detection_day) <- as.character(population$farm_id)
  
  test_target_cull_day <- rep(NA, num_farms)
  names(test_target_cull_day) <- as.character(population$farm_id)
  
  test_cull_type <- rep(NA, num_farms)
  names(test_cull_type) <- as.character(population$farm_id)
  
  # Load Activity
  for (i in 1:nrow(activity_clean)) {
    farm <- as.character(activity_clean$farm_id[i])
    start_d <- activity_clean$date_start[i]
    end_d <- activity_clean$date_end[i]
    active_days <- as.character(sim_days[sim_days >= start_d & sim_days <= end_d])
    if (length(active_days) > 0 && farm %in% rownames(test_matrix)) test_matrix[farm, active_days] <- 1
  }
  
  # Inject Seed Cases
  for (i in 1:nrow(seed_cases)) {
    this_farm <- as.character(seed_cases$farm_id[i])
    d_sick <- seed_cases$date_suspicious[i]
    if (is.na(d_sick)) d_sick <- seed_cases$date_confirmed[i]
    if (!is.na(d_sick)) {
      inf_days <- as.character(sim_days[sim_days >= d_sick])
      if (length(inf_days) > 0 && this_farm %in% rownames(test_matrix)) test_matrix[this_farm, inf_days] <- 3
      day_idx <- which(sim_days == d_sick)
      if(length(day_idx) > 0) test_detection_day[this_farm] <- day_idx
    }
  }
  
  # Run Simulation Forward
  for (t in calib_start_idx:calib_end_idx) {
    today_date <- sim_days[t]
    S_farms <- names(which(test_matrix[, t] == 1))
    I_farms <- names(which(test_matrix[, t] == 3))
    
    if (length(S_farms) == 0) next
    
    spatial_risk <- rep(0, length(S_farms))
    network_risk <- rep(0, length(S_farms))
    names(spatial_risk) <- S_farms
    names(network_risk) <- S_farms
    
    # 10km Bans
    active_outbreaks <- names(which(test_detection_day <= t & test_detection_day >= (t - 28)))
    banned_farms <- c()
    if (length(active_outbreaks) > 0) {
      in_zone <- dist_matrix[, active_outbreaks, drop = FALSE] <= 10
      banned_farms <- rownames(in_zone)[rowSums(in_zone) > 0]
    }
    
    if (length(I_farms) > 0) {
      dist_I_to_S <- dist_matrix[S_farms, I_farms, drop = FALSE]
      decay_matrix <- exp(-dist_I_to_S / spatial_scale)
      
      # ---DUCK AMPLIFICATION IN SPATIAL RISK ---
      I_cap_mult <- cap_mult[I_farms] * species_mult[I_farms]
      
      scaled_decay <- sweep(decay_matrix, 2, I_cap_mult, "*")
      spatial_risk <- test_b_spatial * rowSums(scaled_decay) * cap_mult[S_farms]
      
      today_moves <- movements %>% filter(date == today_date & source_farm %in% I_farms)
      
      if (length(banned_farms) > 0) {
        today_moves <- today_moves %>% filter(!(source_farm %in% banned_farms))
      }
      if (nrow(today_moves) > 0) {
        today_moves <- today_moves %>% filter(!(source_farm %in% hrz_stage1_farms))
      }
      
      if (nrow(today_moves) > 0) {
        vol_by_dest <- tapply(today_moves$volume, today_moves$dest_farm, sum)
        valid_exposed <- as.character(names(vol_by_dest)[names(vol_by_dest) %in% S_farms])
        if (length(valid_exposed) > 0) {
          vol_mult <- vol_by_dest[valid_exposed] / median_vol
          network_risk[valid_exposed] <- test_b_network * vol_mult
        }
      }
    }
    
    total_lambda <- farm_bg_risk[S_farms] + spatial_risk + network_risk
    prob_infection <- 1 - exp(-total_lambda)
    new_infections <- rbinom(n = length(S_farms), size = 1, prob = prob_infection)
    new_I_farms <- S_farms[new_infections == 1]
    
    if (length(new_I_farms) > 0) {
      for (farm in new_I_farms) {
        
        # --- FIXED EMPIRICAL DELAYS BASED ON SPECIES ---
        this_species <- population$species[population$farm_id == farm]
        this_cull_delay <- if(this_species == "chicken") chicken_delays else duck_delays
        
        exp_end_idx <- min(t + latency_delay - 1, num_days)
        sched_cull_idx <- min(t + latency_delay + this_cull_delay - 1, num_days)
        
        test_detection_day[farm] <- max(t, sched_cull_idx - 1)
        
        if (t <= exp_end_idx) test_matrix[farm, t:exp_end_idx] <- 2
        if ((exp_end_idx + 1) <= num_days) test_matrix[farm, (exp_end_idx + 1):num_days] <- 3
        
        test_target_cull_day[farm] <- sched_cull_idx
        test_cull_type[farm] <- "reactive"
        
        neighbors_1km <- names(which(dist_matrix[farm, ] <= 1))
        neighbors_1km <- setdiff(neighbors_1km, farm)
        for (n_farm in neighbors_1km) {
          if (is.na(test_cull_type[n_farm]) || test_cull_type[n_farm] != "reactive") {
            test_target_cull_day[n_farm] <- sched_cull_idx + 1
            test_cull_type[n_farm] <- "preventive"
          }
        }
      }
    }
    
    # QUEUE EXECUTION
    eligible_reactive <- names(which(test_target_cull_day <= t & test_cull_type == "reactive" & test_matrix[, t] != 4))
    eligible_preventive <- names(which(test_target_cull_day <= t & test_cull_type == "preventive" & test_matrix[, t] != 4))
    
    budget <- cull_capacity_per_day
    culled_today <- c()
    
    if (length(eligible_reactive) > 0) {
      sel <- eligible_reactive[1:min(length(eligible_reactive), budget)]
      culled_today <- c(culled_today, sel)
      budget <- budget - length(sel)
    }
    
    if (budget > 0 && length(eligible_preventive) > 0) {
      sel <- eligible_preventive[1:min(length(eligible_preventive), budget)]
      culled_today <- c(culled_today, sel)
    }
    
    if (length(culled_today) > 0) {
      for (c_farm in culled_today) {
        active_future_days <- which(test_matrix[c_farm, t:num_days] != 0)
        if (length(active_future_days) > 0) {
          actual_days <- (t:num_days)[active_future_days]
          test_matrix[c_farm, actual_days] <- 4
        }
      }
    }
  }
  
  # Evaluate
  simulated_total_cases <- sum(test_matrix[, calib_end_idx] %in% c(2, 3, 4))
  error <- abs(simulated_total_cases - target_total_cases)
  
  param_grid$sim_cases[row] <- simulated_total_cases
  param_grid$error[row] <- error
  
  cat("Run", row, "of", nrow(param_grid), "complete. Sim cases:", simulated_total_cases, "(Error:", error, ")\n")
}

# --- 4. VIEW THE TOP RESULTS ---
top_5_params <- param_grid %>% 
  arrange(error) %>% 
  head(5)

print(top_5_params)

# =================================================================
# PART 4: 1000-RUN MONTE CARLO SIMULATION (UNCERTAINTY ANALYSIS & Q2)
# =================================================================

n_sims <- 1000

# Storage for macro results
mc_total_cases <- numeric(n_sims)
mc_backlog <- numeric(n_sims)
mc_type_breakdown <- list()

# TRACKERS FOR QUESTION 2 (Spatial & Temporal Prediction)
mc_daily_matrix <- matrix(0, nrow = n_sims, ncol = length(start_pred_idx:num_days))
colnames(mc_daily_matrix) <- as.character(sim_days[start_pred_idx:num_days])

farm_infection_tallies <- rep(0, num_farms)
names(farm_infection_tallies) <- as.character(population$farm_id)

# -----------------------------------------------------------------
# EMPIRICAL PARAMETERS (Wu et al. & Ssematimba et al.)
# -----------------------------------------------------------------
beta_spatial    <- 0.007    
spatial_scale   <- 1.0     
beta_network    <- 0.04     
cull_capacity_per_day <- 8

latency_delay <- 1         
chicken_delays <- 2       # Chickens show clinical signs rapidly
duck_delays <- 8         # Ducks can be asymptomatic "silent spreaders"

# Ducks transmit at 4.1, Chickens at 1.15. Ducks are ~3.5x more infectious.
species_mult <- ifelse(population$species == "duck", 3.5, 1.0)
names(species_mult) <- as.character(population$farm_id)
# -----------------------------------------------------------------

all_sim_daily_cases <- matrix(0, nrow = n_sims, ncol = num_days)
all_sim_cases_list <- list() 

# --- TAKE THE SNAPSHOT OF "TODAY" ---
state_matrix_baseline <- state_matrix
detection_day_baseline <- detection_day
target_cull_day_baseline <- target_cull_day
cull_type_baseline <- cull_type

for (sim in 1:n_sims) {
  cat("Running Simulation", sim, "of", n_sims, "...\r")
  
  # 1. Reset everything for this specific run
  mc_matrix <- state_matrix_baseline
  mc_detection_day <- detection_day_baseline
  mc_target_cull_day <- target_cull_day_baseline
  mc_cull_type <- cull_type_baseline
  
  daily_new_cases <- rep(0, num_days)
  sim_cases_list <- list() 
  sim_sick_farms <- c()
  
  # 2. The Daily Prediction Loop
  for (t in start_pred_idx:num_days) {
    today_date <- sim_days[t]
    S_farms <- names(which(mc_matrix[, t] == 1))
    I_farms <- names(which(mc_matrix[, t] == 3)) 
    
    if (length(S_farms) == 0) next
    
    spatial_risk <- rep(0, length(S_farms))
    network_risk <- rep(0, length(S_farms))
    names(spatial_risk) <- S_farms
    names(network_risk) <- S_farms
    
    active_outbreaks <- names(which(mc_detection_day <= t & mc_detection_day >= (t - 28)))
    banned_farms <- c()
    if (length(active_outbreaks) > 0) {
      in_zone <- dist_matrix[, active_outbreaks, drop = FALSE] <= 10
      banned_farms <- rownames(in_zone)[rowSums(in_zone) > 0]
    }
    
    if (length(I_farms) > 0) {
      dist_I_to_S <- dist_matrix[S_farms, I_farms, drop = FALSE]
      decay_matrix <- exp(-dist_I_to_S / spatial_scale)
      
      # --- NEW: APPLY DUCK SHEDDING AMPLIFIER ---
      I_cap_mult <- cap_mult[I_farms] * species_mult[I_farms]
      
      scaled_decay <- sweep(decay_matrix, 2, I_cap_mult, "*")
      spatial_risk <- beta_spatial * rowSums(scaled_decay) * cap_mult[S_farms]
      
      today_moves <- movements %>% filter(date == today_date & source_farm %in% I_farms)
      if (length(banned_farms) > 0) today_moves <- today_moves %>% filter(!(source_farm %in% banned_farms))
      if (nrow(today_moves) > 0) today_moves <- today_moves %>% filter(!(source_farm %in% hrz_stage1_farms))
      
      if (nrow(today_moves) > 0) {
        vol_by_dest <- tapply(today_moves$volume, today_moves$dest_farm, sum)
        valid_exposed <- as.character(names(vol_by_dest)[names(vol_by_dest) %in% S_farms])
        if (length(valid_exposed) > 0) {
          vol_mult <- vol_by_dest[valid_exposed] / median_vol
          network_risk[valid_exposed] <- beta_network * vol_mult
        }
      }
    }
    
    total_lambda <- farm_bg_risk[S_farms] + spatial_risk + network_risk
    new_infections <- rbinom(n = length(S_farms), size = 1, prob = 1 - exp(-total_lambda))
    new_I_farms <- S_farms[new_infections == 1]
    
    daily_new_cases[t] <- length(new_I_farms)
    sim_sick_farms <- c(sim_sick_farms, new_I_farms)
    
    if (length(new_I_farms) > 0) {
      sim_cases_list[[t]] <- data.frame(
        sim_id = sim,
        farm_id = as.integer(as.character(new_I_farms)),
        inf_date = as.Date(sim_days[t], origin = "1970-01-01")
      )
      
      for (farm in new_I_farms) {
        sp <- population$species[population$farm_id == farm]
        this_cull_delay <- if(sp == "chicken") chicken_delays else duck_delays        
        exp_end_idx <- min(t + latency_delay - 1, num_days)
        sched_cull_idx <- min(t + latency_delay + this_cull_delay - 1, num_days)
        
        mc_detection_day[farm] <- max(t, sched_cull_idx - 1)
        
        if (t <= exp_end_idx) mc_matrix[farm, t:exp_end_idx] <- 2
        if ((exp_end_idx + 1) <= num_days) mc_matrix[farm, (exp_end_idx + 1):num_days] <- 3
        
        mc_target_cull_day[farm] <- sched_cull_idx
        mc_cull_type[farm] <- "reactive"
        
        neighbors_1km <- setdiff(names(which(dist_matrix[farm, ] <= 1)), farm)
        for (n_farm in neighbors_1km) {
          if (is.na(mc_cull_type[n_farm]) || mc_cull_type[n_farm] != "reactive") {
            mc_target_cull_day[n_farm] <- sched_cull_idx + 1
            mc_cull_type[n_farm] <- "preventive"
          }
        }
      }
    }
    
    eligible_reactive <- names(which(mc_target_cull_day <= t & mc_cull_type == "reactive" & mc_matrix[, t] != 4))
    eligible_preventive <- names(which(mc_target_cull_day <= t & mc_cull_type == "preventive" & mc_matrix[, t] != 4))
    
    budget <- cull_capacity_per_day
    culled_today <- c()
    
    if (length(eligible_reactive) > 0) {
      sel <- eligible_reactive[1:min(length(eligible_reactive), budget)]
      culled_today <- c(culled_today, sel)
      budget <- budget - length(sel)
    }
    if (budget > 0 && length(eligible_preventive) > 0) {
      sel <- eligible_preventive[1:min(length(eligible_preventive), budget)]
      culled_today <- c(culled_today, sel)
    }
    
    if (length(culled_today) > 0) {
      for (c_farm in culled_today) {
        active_future_days <- which(mc_matrix[c_farm, t:num_days] != 0)
        if (length(active_future_days) > 0) mc_matrix[c_farm, (t:num_days)[active_future_days]] <- 4
      }
    }
  } # End of Daily Loop
  
  # --- 3. RECORD RESULTS ---
  all_sim_daily_cases[sim, ] <- daily_new_cases
  mc_total_cases[sim] <- sum(daily_new_cases[start_pred_idx:num_days])
  mc_backlog[sim] <- sum(mc_target_cull_day <= num_days & mc_matrix[, num_days] != 4, na.rm = TRUE)
  
  mc_daily_matrix[sim, ] <- daily_new_cases[start_pred_idx:num_days]
  
  if (length(sim_sick_farms) > 0) {
    farm_infection_tallies[sim_sick_farms] <- farm_infection_tallies[sim_sick_farms] + 1
    
    sim_breakdown <- population %>%
      filter(farm_id %in% as.numeric(sim_sick_farms)) %>%
      group_by(species, production) %>%
      summarise(cases = n(), .groups = "drop") %>%
      mutate(sim_id = sim)
    mc_type_breakdown[[sim]] <- sim_breakdown
  }
  
  if (length(sim_cases_list) > 0) {
    all_sim_cases_list[[sim]] <- bind_rows(sim_cases_list)
  }
}

cat("\nSimulation Complete.\n")

# --- 4. AGGREGATE AND PRINT CONFIDENCE INTERVALS ---
cat("\n======================================================\n")
cat("          MONTE CARLO PREDICTION RESULTS (4 WEEKS)      \n")
cat("======================================================\n")

ci_cases <- quantile(mc_total_cases, probs = c(0.025, 0.5, 0.975))
cat(sprintf("Total New Cases (Median): %d\n", round(ci_cases[2])))
cat(sprintf("95%% Confidence Interval: [%d , %d]\n\n", round(ci_cases[1]), round(ci_cases[3])))

ci_backlog <- quantile(mc_backlog, probs = c(0.025, 0.5, 0.975))
cat(sprintf("Farms in Culling Backlog (Median): %d\n", round(ci_backlog[2])))
cat(sprintf("95%% Confidence Interval: [%d , %d]\n\n", round(ci_backlog[1]), round(ci_backlog[3])))

if (length(mc_type_breakdown) > 0) {
  all_breakdown <- bind_rows(mc_type_breakdown) %>%
    group_by(species, production) %>%
    summarise(
      avg_cases = mean(cases),
      min_cases = min(cases),
      max_cases = max(cases),
      .groups = "drop"
    ) %>%
    arrange(desc(avg_cases))
  
  cat("--- AVERAGE BREAKDOWN BY FARM TYPE ---\n")
  print(all_breakdown)
}


all_breakdown
# --- 5. PLOT THE DISTRIBUTION ---
results_df <- data.frame(TotalCases = mc_total_cases)
p <- ggplot(results_df, aes(x = TotalCases)) +
  geom_histogram(fill = "steelblue", color = "black", bins = 15, alpha = 0.8) +
  geom_vline(aes(xintercept = ci_cases[2]), color = "red", linetype = "dashed", linewidth = 1.2) +
  geom_vline(aes(xintercept = ci_cases[1]), color = "red", linetype = "dotted", linewidth = 1) +
  geom_vline(aes(xintercept = ci_cases[3]), color = "red", linetype = "dotted", linewidth = 1) +
  labs(title = "4-Week Epidemic Prediction ",
       subtitle = "1000 Monte Carlo Simulations (Red line = Median, Dotted = 95% CI)",
       x = "Total New Cases Predicted",
       y = "Frequency (Number of Simulations)") +
  theme_minimal()

print(p)
ggsave("plot_1_4week_epidemic.png", plot = p, width = 8, height = 6, dpi = 300, bg = "white")


# =================================================================
# PART 5: VISUALIZING MONTE CARLO RESULTS (EPICURVE & MAP)
# =================================================================
if (length(all_sim_cases_list) > 0) {
  mc_cases_df <- bind_rows(all_sim_cases_list) %>% 
    left_join(population, by = "farm_id")
  
  start_pred_date <- as.Date(sim_days[start_pred_idx], origin = "1970-01-01")
  end_pred_date   <- as.Date(sim_days[num_days], origin = "1970-01-01")
  full_date_range <- seq.Date(from = start_pred_date, to = end_pred_date, by = "day")
  
  # --- PLOT 1: EXPECTED CUMULATIVE CASES (AVERAGE) ---
  plot1_data <- mc_cases_df %>%
    group_by(inf_date, production) %>%
    summarise(total_mc_cases = n(), .groups = "drop") %>%
    complete(inf_date = full_date_range, 
             production = unique(population$production), 
             fill = list(total_mc_cases = 0)) %>%
    mutate(avg_daily_cases = total_mc_cases / n_sims) %>%
    group_by(production) %>%
    arrange(inf_date) %>%
    mutate(cumulative_avg_cases = cumsum(avg_daily_cases)) %>%
    ungroup()
  
  p3 <- ggplot(plot1_data, aes(x = inf_date, y = cumulative_avg_cases, color = production, group = production)) +
    geom_line(linewidth = 1.2) +                  
    geom_point(size = 1.5, alpha = 0.6) +       
    scale_x_date(date_breaks = "1 week", date_labels = "%b %d") + 
    theme_minimal() +
    labs(
      title = "Expected Cumulative Outbreaks (1000 Monte Carlo Runs)",
      subtitle = "Average total infections over time by production type",
      x = "Date", 
      y = "Total Expected Cases (Running Average)", 
      color = "Production Type"
    ) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  
  print(p3)
  
  ggsave("plot_3_cumulative_outbreaks.png", plot = p3, width = 8, height = 6, dpi = 300, bg = "white")
  
  # --- PLOT 2: SPATIAL MAP (MEDIAN SCENARIO) ---
  median_target <- median(mc_total_cases)
  median_sim_id <- which.min(abs(mc_total_cases - median_target))
  
  cat("\nMapping Median Scenario (Simulation #", median_sim_id, "with", mc_total_cases[median_sim_id], "cases)...\n")
  
  median_cases_df <- mc_cases_df %>%
    filter(sim_id == median_sim_id) %>%
    mutate(
      days_since_start = as.numeric(inf_date - start_pred_date),
      week_num = ceiling((days_since_start + 1) / 7),
      infection_week = paste("Week", week_num)
    )
  
  p2 <- ggplot() +
    geom_point(data = population, aes(x = x, y = y), color = "grey85", size = 0.5, alpha = 0.5) +
    geom_point(data = median_cases_df, aes(x = x, y = y, color = infection_week), size = 2.5, alpha = 0.8) +
    scale_color_brewer(palette = "Set1") + 
    theme_void() +
    labs(
      title = paste("Geographic Spread (Median Scenario - Sim #", median_sim_id, ")"),
      subtitle = "Tracing the movement of the virus across the island by week",
      color = "Infection Wave"
    ) +
    theme(legend.position = "right")
  
  print(p2)
  ggsave("plot_2_map.png", plot = p2, width = 8, height = 6, dpi = 300, bg = "white")
  
  
  # =================================================================
  # QUESTION 2: PREDICTED TEMPORAL & SPATIAL EVOLUTION
  # =================================================================
  cat("\n--- GENERATING FUTURE PREDICTION VISUALS & FILES (Q2) ---\n")
  
  spatial_risk_df <- data.frame(
    farm_id = as.numeric(names(farm_infection_tallies)),
    prob_infection = as.numeric(farm_infection_tallies / n_sims)
  ) %>%
    left_join(population, by = "farm_id")
  
  write.csv(spatial_risk_df, "predicted_spatial_risk.csv", row.names = FALSE)
  cat("Exported 'predicted_spatial_risk.csv' for the Phase 1 Submission!\n")
  
  mc_daily_long <- data.frame(
    sim_id = rep(1:n_sims, each = ncol(mc_daily_matrix)),
    Date   = rep(as.Date(colnames(mc_daily_matrix)), times = n_sims),
    New_Cases = as.vector(t(mc_daily_matrix))
  )
  
  write.csv(mc_daily_long, "q2_temporal_raw_trajectories_long.csv", row.names = FALSE)
  cat("Exported q2_temporal_raw_trajectories_long.csv\n")

} else {
  cat("\nNo cases occurred in any of the 1000 simulations!\n")
}
  
# =================================================================
# PART 5: COMPARATIVE POLICY SCENARIOS (QUESTIONS 4 & 5)
# =================================================================

n_sims_scen <- 100 
scenarios <- c("Baseline", "Scenario_A_Chicken_Only", "Scenario_B_Faster_Reactive")

# Storage for the median results
results_summary <- data.frame(
  Scenario = scenarios,
  Median_Cases = numeric(3),
  Median_Backlog = numeric(3)
)
results_summary

cat(sprintf("\n--- TESTING 3 POLICY SCENARIOS (%d runs each) ---\n", n_sims_scen))

for (scen_idx in 1:3) {
  current_scenario <- scenarios[scen_idx]
  cat(sprintf("\nRunning %s...\n", current_scenario))
  
  scen_cases <- numeric(n_sims_scen)
  scen_backlog <- numeric(n_sims_scen)
  
  for (sim in 1:n_sims_scen) {
    # 1. Reset the World
    state_matrix[,] <- 0
    for (i in 1:nrow(activity_clean)) {
      farm <- as.character(activity_clean$farm_id[i])
      start_d <- activity_clean$date_start[i]
      end_d <- activity_clean$date_end[i]
      active_days <- as.character(sim_days[sim_days >= start_d & sim_days <= end_d])
      if (length(active_days) > 0 && farm %in% rownames(state_matrix)) state_matrix[farm, active_days] <- 1
    }
    
    test_detection_day <- rep(NA, num_farms)
    names(test_detection_day) <- as.character(population$farm_id)
    
    for (i in 1:nrow(cases)) {
      this_farm <- as.character(cases$farm_id[i])
      date_sick <- cases$date_suspicious[i]
      if (is.na(date_sick)) date_sick <- cases$date_confirmed[i]
      date_culled <- cases$cull_start[i]
      if (is.na(date_culled)) date_culled <- end_date + 1 
      if (is.na(date_sick)) next 
      
      infected_days <- as.character(sim_days[sim_days >= date_sick & sim_days < date_culled])
      if (length(infected_days) > 0 && this_farm %in% rownames(state_matrix)) state_matrix[this_farm, infected_days] <- 3 
      culled_days <- as.character(sim_days[sim_days >= date_culled])
      if (length(culled_days) > 0 && this_farm %in% rownames(state_matrix)) state_matrix[this_farm, culled_days] <- 4 
      
      day_idx <- which(sim_days == date_sick)
      if(length(day_idx) > 0) test_detection_day[this_farm] <- day_idx
    }
    
    test_target_cull_day <- rep(NA, num_farms)
    names(test_target_cull_day) <- as.character(population$farm_id)
    test_cull_type <- rep(NA, num_farms) 
    names(test_cull_type) <- as.character(population$farm_id)
    
    daily_new_cases <- rep(0, num_days)
    
    # 2. Run the Prediction Loop
    for (t in start_pred_idx:num_days) {
      today_date <- sim_days[t]
      S_farms <- names(which(state_matrix[, t] == 1))
      I_farms <- names(which(state_matrix[, t] == 3)) 
      
      if (length(S_farms) == 0) next
      
      spatial_risk <- rep(0, length(S_farms))
      network_risk <- rep(0, length(S_farms))
      names(spatial_risk) <- S_farms; names(network_risk) <- S_farms
      
      active_outbreaks <- names(which(test_detection_day <= t & test_detection_day >= (t - 28)))
      banned_farms <- if (length(active_outbreaks) > 0) rownames(dist_matrix[, active_outbreaks, drop = FALSE])[rowSums(dist_matrix[, active_outbreaks, drop = FALSE] <= 10) > 0] else c()
      
      if (length(I_farms) > 0) {
        dist_I_to_S <- dist_matrix[S_farms, I_farms, drop = FALSE]
        spatial_risk <- beta_spatial * rowSums(sweep(exp(-dist_I_to_S / spatial_scale), 2, cap_mult[I_farms], "*")) * cap_mult[S_farms]
        
        today_moves <- movements %>% filter(date == today_date & source_farm %in% I_farms)
        if (length(banned_farms) > 0) today_moves <- today_moves %>% filter(!(source_farm %in% banned_farms))
        if (nrow(today_moves) > 0) today_moves <- today_moves %>% filter(!(source_farm %in% hrz_stage1_farms))
        
        if (nrow(today_moves) > 0) {
          vol_by_dest <- tapply(today_moves$volume, today_moves$dest_farm, sum)
          valid_exposed <- as.character(names(vol_by_dest)[names(vol_by_dest) %in% S_farms])
          if (length(valid_exposed) > 0) network_risk[valid_exposed] <- beta_network * (vol_by_dest[valid_exposed] / median_vol)
        }
      }
      
      total_lambda <- farm_bg_risk[S_farms] + spatial_risk + network_risk
      new_infections <- rbinom(n = length(S_farms), size = 1, prob = 1 - exp(-total_lambda))
      new_I_farms <- S_farms[new_infections == 1]
      daily_new_cases[t] <- length(new_I_farms)
      
      if (length(new_I_farms) > 0) {
        for (farm in new_I_farms) {
          sp <- population$species[population$farm_id == farm]
          this_cull_delay <- if(sp == "chicken") sample(chicken_delays, 1) else sample(duck_delays, 1)
          
          # --- SCENARIO B LOGIC (Q5): SPEED UP REACTIVE CULL BY 1 DAY ---
          if (current_scenario == "Scenario_B_Faster_Reactive") {
            this_cull_delay <- max(0, this_cull_delay - 1)
          }
          
          exp_end_idx <- min(t + latency_delay - 1, num_days)
          sched_cull_idx <- min(t + latency_delay + this_cull_delay - 1, num_days)
          test_detection_day[farm] <- max(t, sched_cull_idx - 1)
          
          if (t <= exp_end_idx) state_matrix[farm, t:exp_end_idx] <- 2
          if ((exp_end_idx + 1) <= num_days) state_matrix[farm, (exp_end_idx + 1):num_days] <- 3
          
          test_target_cull_day[farm] <- sched_cull_idx
          test_cull_type[farm] <- "reactive"
          
          # --- PREVENTIVE CULLING LOGIC ---
          neighbors_1km <- setdiff(names(which(dist_matrix[farm, ] <= 1)), farm)
          
          # SCENARIO A LOGIC (Q4): ONLY CULL CHICKENS PREVENTIVELY
          if (current_scenario == "Scenario_A_Chicken_Only") {
            neighbors_1km <- neighbors_1km[population$species[population$farm_id %in% neighbors_1km] == "chicken"]
          }
          
          # SCENARIO B LOGIC (Q5): NO PREVENTIVE CULLING AT ALL
          if (current_scenario == "Scenario_B_Faster_Reactive") {
            neighbors_1km <- c() 
          }
          
          for (n_farm in neighbors_1km) {
            if (is.na(test_cull_type[n_farm]) || test_cull_type[n_farm] != "reactive") {
              test_target_cull_day[n_farm] <- sched_cull_idx + 1
              test_cull_type[n_farm] <- "preventive"
            }
          }
        }
      }
      
      # 3. Queue Execution
      eligible_reactive <- names(which(test_target_cull_day <= t & test_cull_type == "reactive" & state_matrix[, t] != 4))
      eligible_preventive <- names(which(test_target_cull_day <= t & test_cull_type == "preventive" & state_matrix[, t] != 4))
      
      budget <- cull_capacity_per_day
      culled_today <- c()
      
      if (length(eligible_reactive) > 0) {
        sel <- eligible_reactive[1:min(length(eligible_reactive), budget)]
        culled_today <- c(culled_today, sel)
        budget <- budget - length(sel)
      }
      if (budget > 0 && length(eligible_preventive) > 0) {
        sel <- eligible_preventive[1:min(length(eligible_preventive), budget)]
        culled_today <- c(culled_today, sel)
      }
      if (length(culled_today) > 0) {
        for (c_farm in culled_today) {
          active_future_days <- which(state_matrix[c_farm, t:num_days] != 0)
          if (length(active_future_days) > 0) state_matrix[c_farm, (t:num_days)[active_future_days]] <- 4
        }
      }
    }
    
    scen_cases[sim] <- sum(daily_new_cases[start_pred_idx:num_days])
    scen_backlog[sim] <- sum(test_target_cull_day <= num_days & state_matrix[, num_days] != 4, na.rm = TRUE)
  }
  
  results_summary$Median_Cases[scen_idx] <- median(scen_cases)
  results_summary$Median_Backlog[scen_idx] <- median(scen_backlog)
}

cat("\n======================================================\n")
cat("          POLICY SCENARIO RESULTS (Q4 & Q5)             \n")
cat("======================================================\n")
print(results_summary)

# Plot the comparison
# 
p_compare <- ggplot(results_summary, aes(x = Scenario, y = Median_Cases, fill = Scenario)) +
  geom_bar(stat = "identity", color = "black") +
  theme_minimal() +
  labs(
    title = "Comparison of Management Strategies",
    subtitle = "Median predicted cases over 4 weeks",
    x = "Policy Scenario",
    y = "Predicted New Cases"
  ) +
  theme(legend.position = "none")

print(p_compare)

ggsave("plot_comparison.png", plot = p_compare, width = 8, height = 6, dpi = 300, bg = "white")



# ==========================================
# PART 6: EXPORTING RESULTS TABLES
# ==========================================

# 1. Export Historical Outbreak Distribution (103 cases)
write_xlsx(summary_table, "historical_outbreak_distribution.xlsx")

# 2. Export 4-Week Prediction Breakdown by Farm Type
write_xlsx(all_breakdown, "predicted_cases_by_farm_type.xlsx")

# 3. Export Policy Scenario Comparison (Baseline vs A vs B)
write_xlsx(results_summary, "policy_scenario_comparison.xlsx")