# Clear workspace and set reproducible seed
rm(list = ls())
set.seed(2612)

# Load required packages
required_packages <- c(
  "tidyverse", "readxl", "here", "survival", "truncnorm", "forcats", "viridis",
  "tictoc", "janitor", "scales", "progressr", "gt", "jsonlite", "lubridate"
)

# Function to install missing packages
install_if_missing <- function(packages) {
  missing_packages <- packages[!packages %in% rownames(installed.packages())]
  if (length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    install.packages(missing_packages, dependencies = TRUE)
  }
}

install_if_missing(required_packages)
invisible(lapply(required_packages, library, character.only = TRUE))

# ===================================================================
# LOAD AE DATA
# ===================================================================

ae_df <- read_rds(here("data", "processed", "matched_data.rds")) %>%
  filter(drug == "bedinvetmab",
         date <= as.Date("2025-06-30"))

# ===================================================================
# 1. SETUP: Country-specific reporting rates
# ===================================================================

reporting_rates <- tribble(
  ~country_group, ~ae_per_10k_doses,
  "United States", 18.18,
  "United Kingdom", 14.61,
  "Germany", 6.06,
  "Spain", 5.41,
  "France", 3.22,
  "Italy", 3.22,
  "Canada", 18.72,
  "Australia", 13.86,
  "EU_combined", mean(c(3.22, 6.06, 5.41, 3.22)),  # Mean of FR, DE, ES, IT
  "Rest_of_World", 9.48  # Global average
)

# ===================================================================
# 2. LIFE EXPECTANCY TABLE BY SIZE AND AGE (Montoya et al. 2019 data from Supp Table 8)
# ===================================================================

life_table_by_size <- tribble(
  ~size_group, ~age_lower, ~age_upper, ~e_x_years, ~ci_lower, ~ci_upper,
  # Toy
  "Toy", 0, 1, 13.58, 13.52, 13.63,
  "Toy", 1, 2, 12.97, 12.92, 13.02,
  "Toy", 2, 3, 12.09, 12.04, 12.14,
  "Toy", 3, 4, 11.19, 11.14, 11.24,
  "Toy", 4, 5, 10.28, 10.24, 10.33,
  "Toy", 5, 6, 9.38, 9.34, 9.43,
  "Toy", 6, 7, 8.49, 8.45, 8.54,
  "Toy", 7, 8, 7.62, 7.58, 7.67,
  "Toy", 8, 9, 6.77, 6.73, 6.81,
  "Toy", 9, 10, 5.97, 5.92, 6.01,
  "Toy", 10, 11, 5.22, 5.17, 5.26,
  "Toy", 11, 12, 4.52, 4.48, 4.57,
  "Toy", 12, 13, 3.88, 3.84, 3.93,
  "Toy", 13, 14, 3.32, 3.28, 3.37,
  "Toy", 14, 15, 2.85, 2.80, 2.91,
  "Toy", 15, 16, 2.47, 2.41, 2.53,
  "Toy", 16, 17, 2.23, 2.15, 2.31,
  "Toy", 17, 99, 2.08, 1.97, 2.18,
  # Small
  "Small", 0, 1, 13.79, 13.74, 13.83,
  "Small", 1, 2, 13.02, 12.98, 13.06,
  "Small", 2, 3, 12.10, 12.07, 12.14,
  "Small", 3, 4, 11.18, 11.14, 11.21,
  "Small", 4, 5, 10.26, 10.22, 10.30,
  "Small", 5, 6, 9.35, 9.31, 9.39,
  "Small", 6, 7, 8.44, 8.41, 8.48,
  "Small", 7, 8, 7.56, 7.52, 7.59,
  "Small", 8, 9, 6.69, 6.65, 6.72,
  "Small", 9, 10, 5.88, 5.85, 5.91,
  "Small", 10, 11, 5.11, 5.08, 5.15,
  "Small", 11, 12, 4.41, 4.37, 4.44,
  "Small", 12, 13, 3.76, 3.73, 3.79,
  "Small", 13, 14, 3.20, 3.16, 3.23,
  "Small", 14, 15, 2.73, 2.69, 2.77,
  "Small", 15, 16, 2.36, 2.31, 2.40,
  "Small", 16, 17, 2.10, 2.05, 2.16,
  "Small", 17, 99, 1.97, 1.89, 2.04,
  # Medium
  "Medium", 0, 1, 12.94, 12.88, 12.99,
  "Medium", 1, 2, 12.14, 12.09, 12.20,
  "Medium", 2, 3, 11.23, 11.18, 11.28,
  "Medium", 3, 4, 10.31, 10.26, 10.36,
  "Medium", 4, 5, 9.38, 9.33, 9.44,
  "Medium", 5, 6, 8.49, 8.44, 8.54,
  "Medium", 6, 7, 7.58, 7.53, 7.63,
  "Medium", 7, 8, 6.71, 6.66, 6.76,
  "Medium", 8, 9, 5.87, 5.82, 5.91,
  "Medium", 9, 10, 5.09, 5.05, 5.14,
  "Medium", 10, 11, 4.36, 4.32, 4.41,
  "Medium", 11, 12, 3.73, 3.69, 3.77,
  "Medium", 12, 13, 3.16, 3.11, 3.20,
  "Medium", 13, 14, 2.69, 2.64, 2.73,
  "Medium", 14, 15, 2.32, 2.26, 2.38,
  "Medium", 15, 16, 2.06, 1.99, 2.13,
  "Medium", 16, 17, 1.90, 1.80, 1.99,
  "Medium", 17, 99, 1.90, 1.75, 2.05,
  # Large
  "Large", 0, 1, 11.70, 11.66, 11.73,
  "Large", 1, 2, 10.97, 10.93, 11.00,
  "Large", 2, 3, 10.08, 10.05, 10.12,
  "Large", 3, 4, 9.19, 9.15, 9.22,
  "Large", 4, 5, 8.28, 8.25, 8.31,
  "Large", 5, 6, 7.40, 7.37, 7.43,
  "Large", 6, 7, 6.54, 6.51, 6.57,
  "Large", 7, 8, 5.71, 5.68, 5.74,
  "Large", 8, 9, 4.93, 4.90, 4.96,
  "Large", 9, 10, 4.22, 4.19, 4.25,
  "Large", 10, 11, 3.62, 3.59, 3.65,
  "Large", 11, 12, 3.08, 3.05, 3.11,
  "Large", 12, 13, 2.62, 2.58, 2.65,
  "Large", 13, 14, 2.25, 2.21, 2.28,
  "Large", 14, 15, 1.98, 1.94, 2.02,
  "Large", 15, 16, 1.79, 1.73, 1.85,
  "Large", 16, 17, 1.66, 1.58, 1.75,
  "Large", 17, 99, 1.57, 1.44, 1.69,
  # Giant
  "Giant", 0, 1, 9.70, 9.54, 9.85,
  "Giant", 1, 2, 8.93, 8.78, 9.08,
  "Giant", 2, 3, 8.08, 7.93, 8.23,
  "Giant", 3, 4, 7.21, 7.06, 7.36,
  "Giant", 4, 5, 6.37, 6.22, 6.51,
  "Giant", 5, 6, 5.59, 5.44, 5.73,
  "Giant", 6, 7, 4.87, 4.72, 5.01,
  "Giant", 7, 8, 4.21, 4.07, 4.36,
  "Giant", 8, 9, 3.64, 3.50, 3.79,
  "Giant", 9, 10, 3.18, 3.02, 3.34,
  "Giant", 10, 11, 2.82, 2.64, 2.99,
  "Giant", 11, 12, 2.52, 2.31, 2.72,
  "Giant", 12, 13, 2.32, 2.06, 2.57,
  "Giant", 13, 99, 2.08, 1.75, 2.41
) %>%
  mutate(
    sd_years = (ci_upper - ci_lower) / (2 * 1.96),
    e_x_months = e_x_years * 12,
    sd_months = sd_years * 12
  )

# ===================================================================
# 3. PROCESS AE DATA BY COUNTRY
# ===================================================================

eu_countries <- c(
  "Austria", "Belgium", "Bulgaria", "Croatia", "Czech Republic",
  "Denmark", "Estonia", "Finland", "France", "Germany", "Greece",
  "Hungary", "Iceland", "Ireland", "Italy", "Lithuania", "Luxembourg",
  "Netherlands", "Norway", "Poland", "Portugal", "Slovakia", "Slovenia",
  "Spain", "Sweden", "Switzerland", "Andorra"
)

# Create analysis dataset
analysis_data <- ae_df %>%
  group_by(country) %>%
  summarise(n_aes = n(), .groups = "drop") %>%
  mutate(
    country_group = case_when(
      country == "United States" ~ "United States",
      country == "United Kingdom" ~ "United Kingdom", 
      country == "Canada" ~ "Canada",
      country == "France" ~ "France",
      country == "Germany" ~ "Germany",
      country == "Italy" ~ "Italy",
      country == "Spain" ~ "Spain",
      country == "Australia" ~ "Australia",
      country %in% eu_countries ~ "EU_combined",
      TRUE ~ "Rest_of_World"
    )
  ) %>%
  group_by(country_group) %>%
  summarise(total_aes = sum(n_aes), .groups = "drop") %>%
  left_join(reporting_rates, by = "country_group") %>%
  mutate(
    license_date = case_when(
      country_group == "United States" ~ as_date("2023-05-05"),
      country_group == "United Kingdom" ~ as_date("2020-11-10"),
      country_group == "Canada" ~ as_date("2023-03-08"),
      country_group == "Australia" ~ as_date("2022-09-01"),
      country_group %in% c("France", "Germany", "Italy", "Spain", "EU_combined") ~ as_date("2020-11-10"),
      country_group == "Rest_of_World" ~ as_date("2022-01-01")
    ),
    estimated_doses = total_aes / (ae_per_10k_doses / 10000),
    market_share = estimated_doses / sum(estimated_doses)
  ) %>%
  arrange(desc(market_share))

cat("Market composition based on AE reports + reporting rates:\n")
print(analysis_data)

# ===================================================================
# EMPIRICAL DATA (REPORTING-RATE CORRECTED)
# ===================================================================

empirical_uptake <- ae_df %>%
  mutate(
    year_month = floor_date(date, "month"),
    country_group = case_when(
      country == "United States" ~ "United States",
      country == "United Kingdom" ~ "United Kingdom",
      country == "Canada" ~ "Canada",
      country == "France" ~ "France",
      country == "Germany" ~ "Germany",
      country == "Italy" ~ "Italy",
      country == "Spain" ~ "Spain",
      country == "Australia" ~ "Australia",
      country %in% eu_countries ~ "EU_combined",
      TRUE ~ "Rest_of_World"
    )
  ) %>%
  count(country_group, year_month, name = "n_aes") %>%
  left_join(analysis_data %>% select(country_group, license_date), by = "country_group") %>%
  filter(year_month >= license_date) %>%
  left_join(reporting_rates %>% select(country_group, ae_per_10k_doses), by = "country_group") %>%
  mutate(
    months_since_launch = as.numeric(interval(license_date, year_month) / months(1)),
    estimated_doses = n_aes / (ae_per_10k_doses / 10000)
  ) %>%
  group_by(country_group) %>%
  arrange(year_month) %>%
  mutate(cumulative_doses = cumsum(estimated_doses)) %>%
  ungroup()

# ===================================================================
# FIT LOGISTIC MODEL with fixed K based on market saturation
# ===================================================================

# Calculate K externally (not fitted)
calculate_K <- function(current_cumulative, months_available) {
  # Assume markets saturate in 72 months
  t_saturation <- 72
  r_assumed <- 0.1
  t0_assumed <- t_saturation / 2
  
  # Where should we be on S-curve now?
  expected_prop <- 1 / (1 + exp(-r_assumed * (months_available - t0_assumed)))
  
  # K = current / expected_proportion
  K_estimated <- current_cumulative / expected_prop
  
  return(K_estimated)
}

# Fit only r and t0 with K fixed
fit_logistic_fixed_K <- function(data, K_fixed) {
  max_time <- max(data$months_since_launch)
  
  tryCatch({
    nls(
      cumulative_doses ~ K_fixed / (1 + exp(-r * (months_since_launch - t0))),
      data = data,
      start = list(r = 0.1, t0 = max_time * 0.5),
      control = nls.control(maxiter = 200)
    )
  }, error = function(e) NULL)
}

# Apply fixed-K approach
logistic_params <- empirical_uptake %>%
  group_by(country_group) %>%
  summarise(
    months_available = max(months_since_launch),
    current = max(cumulative_doses),
    
    # K calculated, not fitted
    K = calculate_K(current, months_available),
    
    # Fit only r and t0
    model = list(fit_logistic_fixed_K(cur_data(), K)),
    r = ifelse(!is.null(model[[1]]), coef(model[[1]])["r"], 0.1),
    t0 = ifelse(!is.null(model[[1]]), coef(model[[1]])["t0"], months_available * 0.5),
    
    pct_capacity = 100 * current / K,
    .groups = "drop"
  )

print("Fixed-K logistic parameters:")
print(logistic_params %>%
        select(country_group, months_available, K, r, t0, pct_capacity) %>%
        mutate(
          K = scales::comma(round(K)),
          r = round(r, 3),
          t0 = round(t0, 1),
          pct_capacity = sprintf("%.1f%%", pct_capacity)
        ))

# Validate fit quality
validation <- empirical_uptake %>%
  left_join(logistic_params %>% select(country_group, K, r, t0), by = "country_group") %>%
  mutate(predicted = K / (1 + exp(-r * (months_since_launch - t0)))) %>%
  group_by(country_group) %>%
  summarise(r_squared = cor(cumulative_doses, predicted, use = "complete.obs")^2, .groups = "drop")

print("\nR² with fixed K:")
print(validation %>% mutate(r_squared = round(r_squared, 3)))

fitted_predictions_k <- empirical_uptake %>%
  left_join(logistic_params %>% select(country_group, K, r, t0), by = "country_group") %>%
  filter(!is.na(K)) %>%
  mutate(
    fitted_doses = K / (1 + exp(-r * (months_since_launch - t0)))
  )

# Combine all three datasets

comparison_data_k <- fitted_predictions_k %>%
  select(country_group, year_month, months_since_launch, 
         empirical = cumulative_doses, fitted = fitted_doses) %>%
  pivot_longer(cols = c(empirical, fitted), names_to = "source", values_to = "cumulative")

# Plot the fit
logistic_comparison_plot <- comparison_data_k %>%
  filter(country_group %in% c("United States", "United Kingdom", "France", "Germany",
                              "EU_combined", "Italy", "Spain", "Rest_of_World",
                              "Australia", "Canada")) %>%
  mutate(
    country_group = recode(
      country_group,
      "EU_combined"   = "European Union",
      "Rest_of_World" = "Rest of World"
    ),
    country_group = factor(country_group, levels = c(
      "United States", "United Kingdom", "France", "Germany",
      "European Union", "Italy", "Spain", "Rest of World",
      "Australia", "Canada"
    ))
  ) %>%
  ggplot(aes(x = year_month, y = cumulative, color = source)) +
  geom_point(data = . %>% filter(source == "empirical"), size = 2, alpha = 0.8) +
  geom_line(data = . %>% filter(source == "fitted"), linewidth = 1.3, alpha = 0.9) +
  facet_wrap(~country_group, scales = "free_y", ncol = 2) +
  scale_y_continuous(labels = comma) +
  scale_color_manual(
    values = c("empirical" = "#d62728", "fitted" = "#1f77b4"),
    labels = c("Empirical (dose-corrected)", "Logistic Fit (Fixed K)")
  ) +
  labs(
    title = "Logistic Model: Empirical vs Fitted Cumulative Dose-Equivalents",
    subtitle = "Empirical data corrected for reporting rates | Four major reporting regions",
    x = "Date",
    y = "Cumulative Dose-Equivalents",
    color = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background   = element_rect(fill = "#fff1e5", color = NA),
    panel.background  = element_rect(fill = "#fff1e5", color = NA),
    strip.background  = element_rect(fill = "#f2e6dd", color = NA),
    panel.grid.major  = element_line(color = "#f0f0f0", linewidth = 0.5),
    panel.grid.minor  = element_blank(),
    text              = element_text(color = "#2d2d2d"),
    plot.title        = element_text(size = 14, face = "bold", margin = margin(b = 5)),
    plot.subtitle     = element_text(size = 11, color = "#555555", margin = margin(b = 15)),
    axis.title        = element_text(size = 11, color = "#333333"),
    axis.text         = element_text(size = 9, color = "#555555"),
    strip.text        = element_text(size = 10, face = "bold", color = "#2d2d2d"),
    legend.position   = "top",
    legend.title      = element_blank(),
    legend.text       = element_text(size = 10, color = "#2d2d2d"),
    legend.background = element_rect(fill = "#fff1e5", color = NA),
    legend.key        = element_rect(fill = "#fff1e5", color = NA),
    plot.margin       = margin(t = 15, r = 15, b = 15, l = 15),
    panel.spacing     = unit(1, "lines")
  )

logistic_comparison_plot

ggsave(here("output", "figures", "logistic_comparison_plot.png"), 
       logistic_comparison_plot, 
       width = 12, height = 8, dpi = 300, bg = "#fff1e5")

analysis_data <- analysis_data %>%
  left_join(logistic_params, by = "country_group")

# ===================================================================
# 4. VECTORIZED HELPER FUNCTIONS ⚡=
# ===================================================================

# Sample dog weight from Zoetis study distribution
# Median: 26.06 kg, IQR: 16.0-34.6 kg → SD ≈ 13.78
sample_dog_weight <- function(n = 1) {
  rtruncnorm(n, a = 1, b = 100, mean = 26.06, sd = 13.78)
}

# Assign size group based on weight
assign_size_group <- function(weights) {
  case_when(
    weights < 5.5 ~ "Toy",
    weights < 11 ~ "Small",
    weights < 26 ~ "Medium",
    weights < 45 ~ "Large",
    TRUE ~ "Giant"
  )
}

# Sample dog age (still use global distribution)
sample_dog_age <- function(n = 1) {
  rtruncnorm(n, a = 1, b = 18, mean = 11.4, sd = 2.0)
}

# VECTORIZED: Survival months based on size AND age
survival_months_from_size_and_age_vectorized <- function(ages, size_groups) {
  n <- length(ages)
  survival_months <- numeric(n)
  
  # Create lookup key for each dog
  dogs_df <- tibble(
    idx = 1:n,
    age = ages,
    size_group = size_groups
  )
  
  # For each unique combination of size and age bracket
  for (size in unique(size_groups)) {
    size_rows <- life_table_by_size %>% filter(size_group == size)
    
    for (i in 1:nrow(size_rows)) {
      mask <- (dogs_df$size_group == size) & 
        (dogs_df$age >= size_rows$age_lower[i]) & 
        (dogs_df$age < size_rows$age_upper[i])
      
      n_in_group <- sum(mask)
      
      if (n_in_group > 0) {
        survival_months[mask] <- rtruncnorm(
          n_in_group,
          a = size_rows$ci_lower[i] * 12,
          b = size_rows$ci_upper[i] * 12,
          mean = size_rows$e_x_months[i],
          sd = size_rows$sd_months[i]
        )
      }
    }
  }
  
  return(pmax(survival_months, 1))
}

# Treatment persistence 
treatment_persistence_months_vectorized <- function(survival_months) {
  n <- length(survival_months)
  max_months <- ceiling(survival_months)
  result <- numeric(n)
  u_draws <- runif(n)
  
  for (i in 1:n) {
    if (max_months[i] < 1) {
      result[i] <- 0
      next
    }
    
    months_seq <- seq_len(max_months[i])
    p_continue <- ifelse(months_seq <= 12, 0.98, 1.00)
    cum_survival <- cumprod(p_continue)
    
    discontinue <- which(cum_survival < u_draws[i])[1]
    result[i] <- ifelse(is.na(discontinue), max_months[i], discontinue)
  }
  
  return(result)
}
# ===================================================================
# 5. ITERATIVE SIMULATION TO EXACT DOSE TARGET
# ===================================================================

cutoff_date <- as_date("2025-06-30")
known_global_doses_sold <- 30e6
dose_buffer_pct <- 0.10
target_doses_consumed <- known_global_doses_sold * (1 - dose_buffer_pct)

cat(sprintf("
═══════════════════════════════════════════════════════════
ITERATIVE SIMULATION TO EXACT DOSE TARGET
═══════════════════════════════════════════════════════════

Target doses sold:                    %s
Dose buffer (10%%):                    %s  
Target doses consumed by dogs:        %s
Cutoff date:                          %s

Simulating dogs until cumulative doses = %.2fM...
═══════════════════════════════════════════════════════════\n",
            scales::comma(known_global_doses_sold),
            scales::comma(known_global_doses_sold * dose_buffer_pct),
            scales::comma(target_doses_consumed),
            as.character(cutoff_date),
            target_doses_consumed / 1e6
))

# Initialize
dogs_list <- list()
cumulative_doses <- 0
batch_num <- 0

# Start timer and progress tracking
tic()
handlers(global = TRUE)
handlers("progress")

with_progress({
  p <- progressor(steps = 99)  # 100 progress updates
  last_progress <- 0
  
  while (cumulative_doses < target_doses_consumed) {
    batch_num <- batch_num + 1
    
    # ADAPTIVE BATCH SIZE
    remaining_pct <- cumulative_doses / target_doses_consumed
    
    batch_size <- case_when(
      remaining_pct < 0.90 ~ 50000,   # Far from target: large batches
      remaining_pct < 0.96 ~ 10000,   # Getting close: medium batches
      TRUE ~ 2000                      # Very close: small batches
    )
    
    # Generate batch
    batch_dogs <- tibble(
      dog_id = (batch_num - 1) * 50000 + 1:batch_size,
      country_group = sample(
        analysis_data$country_group,
        size = batch_size,
        replace = TRUE,
        prob = analysis_data$market_share
      )
    ) %>%
      left_join(logistic_params %>% select(country_group, K, r, t0), by = "country_group") %>%
      left_join(analysis_data %>% select(country_group, license_date), 
                by = "country_group") %>%
      mutate(
        # 1️⃣ Total months market is available
        total_months = as.numeric(interval(license_date, cutoff_date) / months(1)),
        
        # 2️⃣ Cumulative maximum dogs (doses) possible by cutoff
        N_total = K / (1 + exp(-r * (total_months - t0))),
        
        # 3️⃣ Invert the logistic function to sample entry time
        u = runif(n()) * N_total,
        months_to_entry = t0 - (1 / r) * log(pmax(K / u - 1, 0.001)),
        months_to_entry = pmin(pmax(months_to_entry, 0), total_months),
        entry_date = license_date + months(round(months_to_entry)),
        
        # NEW: Weight and size
        weight_kg = sample_dog_weight(batch_size),
        age_at_start = sample_dog_age(batch_size)
      )
    
    # Assign size group
    batch_dogs$size_group <- assign_size_group(batch_dogs$weight_kg)
    
    # VECTORIZED: Survival based on size AND age
    batch_dogs$survival_months <- survival_months_from_size_and_age_vectorized(
      batch_dogs$age_at_start,
      batch_dogs$size_group
    )
    
    # Treatment persistence
    batch_dogs$treatment_months <- treatment_persistence_months_vectorized(
      batch_dogs$survival_months
    )
    
    # Finish processing
    batch_dogs <- batch_dogs %>%
      mutate(
        months_to_cutoff = as.numeric(interval(entry_date, cutoff_date) / months(1)),
        doses_by_cutoff = pmin(treatment_months, pmax(0, months_to_cutoff)) %>% ceiling()
      ) %>%
      filter(entry_date <= cutoff_date, doses_by_cutoff > 0)
    
    # Update counters
    dogs_list[[batch_num]] <- batch_dogs
    cumulative_doses <- cumulative_doses + sum(batch_dogs$doses_by_cutoff)
    
    # Update progress (cap at 99% to avoid overshoot warning)
    current_progress <- min(99, floor(100 * cumulative_doses / target_doses_consumed))
    if (current_progress > last_progress) {
      p(
        message = sprintf("Batch %d (n=%s): %.2fM doses | %.1f%% complete",
                          batch_num,
                          scales::comma(batch_size),
                          cumulative_doses / 1e6,
                          current_progress),
        amount = current_progress - last_progress
      )
      last_progress <- current_progress
    }
    
    # Safety check
    if (batch_num > 300) {
      warning("Exceeded maximum iterations")
      break
    }
  }
})

toc()

# Combine all batches
dogs_all <- bind_rows(dogs_list)

cat(sprintf("\n✓ Generated %s dogs consuming %.3fM doses\n",
            scales::comma(nrow(dogs_all)),
            cumulative_doses / 1e6))

# ===================================================================
# 5. TRIM TO EXACT TARGET
# ===================================================================

if (cumulative_doses > target_doses_consumed) {
  excess_doses <- cumulative_doses - target_doses_consumed
  cat(sprintf("⚠ Overshot by %.3fM doses (%.3f%%), trimming...\n", 
              excess_doses / 1e6, 
              100 * excess_doses / target_doses_consumed))
  
  dogs_sorted <- dogs_all %>% arrange(desc(doses_by_cutoff))
  dogs_to_remove <- 0
  running_total <- cumulative_doses
  
  for (i in 1:nrow(dogs_sorted)) {
    if (running_total - dogs_sorted$doses_by_cutoff[i] >= target_doses_consumed) {
      running_total <- running_total - dogs_sorted$doses_by_cutoff[i]
      dogs_to_remove <- dogs_to_remove + 1
    } else {
      break
    }
  }
  
  dogs_final <- dogs_all %>%
    arrange(desc(doses_by_cutoff)) %>%
    slice(-(1:dogs_to_remove))
  
  final_doses <- sum(dogs_final$doses_by_cutoff)
  
  cat(sprintf("✓ Removed %s dogs, final: %s dogs @ %.3fM doses\n\n",
              scales::comma(dogs_to_remove),
              scales::comma(nrow(dogs_final)),
              final_doses / 1e6))
} else {
  dogs_final <- dogs_all
  final_doses <- cumulative_doses
  cat("✓ Hit target exactly\n\n")
}

# ===================================================================
# 6. RESULTS WITH SIZE BREAKDOWN
# ===================================================================

# Summary by country
country_summary <- dogs_final %>%
  group_by(country_group) %>%
  summarise(
    n_dogs = n(),
    total_doses_consumed = sum(doses_by_cutoff),
    mean_age = mean(age_at_start),
    mean_weight = mean(weight_kg),
    median_weight = median(weight_kg),
    mean_duration = mean(doses_by_cutoff),
    median_duration = median(doses_by_cutoff),
    .groups = "drop"
  ) %>%
  left_join(analysis_data %>% select(country_group, total_aes), by = "country_group") %>%
  mutate(
    total_doses_with_buffer = total_doses_consumed / (1 - dose_buffer_pct),
    ae_per_1000_dogs = total_aes / n_dogs * 1000
  ) %>%
  arrange(desc(n_dogs))

# Summary by size group
size_summary <- dogs_final %>%
  group_by(size_group) %>%
  summarise(
    n_dogs = n(),
    pct_of_total = 100 * n() / nrow(dogs_final),
    mean_weight = mean(weight_kg),
    mean_age = mean(age_at_start),
    mean_survival = mean(survival_months / 12),
    mean_treatment = mean(treatment_months),
    mean_doses = mean(doses_by_cutoff),
    .groups = "drop"
  ) %>%
  arrange(factor(size_group, levels = c("Toy", "Small", "Medium", "Large", "Giant")))

cat("═══════════════════════════════════════════════════════════\n")
cat("RESULTS BY SIZE GROUP\n")
cat("═══════════════════════════════════════════════════════════\n\n")

size_summary %>%
  mutate(across(where(is.numeric), ~round(., 1))) %>%
  mutate(n_dogs = scales::comma(n_dogs)) %>%
  print(n = 20, width = Inf)

# Check why mean number of doses is the same:
dogs_final %>%
  mutate(
    time_since_entry = as.numeric(interval(entry_date, cutoff_date) / months(1))
  ) %>%
  summarise(mean_time_since_entry = mean(time_since_entry))

cat("\n═══════════════════════════════════════════════════════════\n")
cat("RESULTS BY COUNTRY\n")
cat("═══════════════════════════════════════════════════════════\n\n")

country_summary %>%
  select(country_group, n_dogs, mean_weight, mean_age, mean_duration, 
         total_aes, ae_per_1000_dogs) %>%
  mutate(
    n_dogs = scales::comma(n_dogs),
    across(c(mean_weight, mean_age, mean_duration, ae_per_1000_dogs), ~round(., 1))
  ) %>%
  print(n = 20, width = Inf)

# ===================================================================
# 7. GLOBAL SUMMARY
# ===================================================================

total_dogs <- nrow(dogs_final)
total_doses <- sum(dogs_final$doses_by_cutoff)
total_doses_with_buffer <- total_doses / (1 - dose_buffer_pct)
total_aes <- sum(country_summary$total_aes, na.rm = TRUE)

cat(sprintf("\n
═══════════════════════════════════════════════════════════
GLOBAL SUMMARY (with size/weight modeling)
═══════════════════════════════════════════════════════════

DOSES:
------
Doses consumed:                       %s
Dose buffer (10%%):                    %s
Total doses sold:                     %s ✓
Accuracy:                             %.5f%%

POPULATION:
-----------
Total unique dogs:                    %s
Mean weight:                          %.1f kg (median: %.1f kg)
Mean age:                             %.1f years
Average treatment duration:           %.1f months

SIZE DISTRIBUTION:
------------------
Toy (<5.5kg):                         %.1f%%
Small (5.5-11kg):                     %.1f%%
Medium (11-26kg):                     %.1f%%
Large (26-45kg):                      %.1f%%
Giant (≥45kg):                        %.1f%%

ADVERSE EVENTS:
---------------
Total AE reports:                     %s
Global AE rate per 1,000 dogs:        %.1f
Global AE rate per 10,000 doses:      %.1f

ASSUMPTIONS:
------------
- Weight: truncnorm(mean=26.06, sd=13.78, min=1, max=100)
- Size-specific survival curves from Montoya 2024 (2019 data)
- Dropout: 3%% monthly year 1, 0%% thereafter
- Dose buffer: 10%%
- Cutoff date: %s

═══════════════════════════════════════════════════════════
",
            scales::comma(round(total_doses)),
            scales::comma(round(total_doses_with_buffer - total_doses)),
            scales::comma(known_global_doses_sold),
            100 * abs(total_doses_with_buffer - known_global_doses_sold) / known_global_doses_sold,
            scales::comma(total_dogs),
            mean(dogs_final$weight_kg),
            median(dogs_final$weight_kg),
            mean(dogs_final$age_at_start),
            mean(dogs_final$doses_by_cutoff),
            100 * mean(dogs_final$size_group == "Toy"),
            100 * mean(dogs_final$size_group == "Small"),
            100 * mean(dogs_final$size_group == "Medium"),
            100 * mean(dogs_final$size_group == "Large"),
            100 * mean(dogs_final$size_group == "Giant"),
            scales::comma(total_aes),
            total_aes / total_dogs * 1000,
            total_aes / total_doses_with_buffer * 10000,
            as.character(cutoff_date)
))

# Optional: Save results
# write_rds(dogs_final, here("output", "simulated_dogs_with_size.rds"))
# write_csv(size_summary, here("output", "size_summary.csv"))
# write_csv(country_summary, here("output", "country_summary.csv"))

cat("\n✓ Simulation complete!\n")

exposure_summary <- dogs_final %>%
  count(country_group, name = "n_exposed")

# ===================================================================
# 2. AGGREGATE ADRs (numerator)
# ===================================================================
adr_summary <- ae_df %>%
  mutate(country_group = case_when(
    country == "United States" ~ "United States",
    country == "United Kingdom" ~ "United Kingdom",
    country == "Canada" ~ "Canada",
    country == "France" ~ "France",
    country == "Germany" ~ "Germany",
    country == "Italy" ~ "Italy",
    country == "Spain" ~ "Spain",
    country == "Australia" ~ "Australia",
    country %in% eu_countries ~ "EU_combined",
    TRUE ~ "Rest_of_World"
  )) %>%
  count(country_group, pt, name = "n_adrs")

# ===================================================================
# 3. COMBINE AND CALCULATE RATES
# ===================================================================
adr_rates <- adr_summary %>%
  left_join(exposure_summary, by = "country_group") %>%
  mutate(rate_per_10000 = 10000 * n_adrs / n_exposed)

# ===================================================================
# 4. SELECT TOP 20 ADRs GLOBALLY
# ===================================================================
top_pts <- adr_rates %>%
  filter(pt != "Lack of efficacy") %>%
  group_by(pt) %>%
  summarise(total_adrs = sum(n_adrs, na.rm = TRUE), .groups = "drop") %>%
  slice_max(total_adrs, n = 20) %>%
  pull(pt)

adr_top10 <- adr_rates %>%
  filter(pt %in% top_pts, pt != "Lack of efficacy")

# Clean names for display
adr_top10_clean <- adr_top10 %>%
  mutate(
    pt = str_to_title(pt) |> str_replace_all("_", " "),
    country_group = str_replace_all(country_group, "_", " "),
    country_group = str_to_title(country_group),
    pt = fct_reorder(pt, -n_adrs),
    country_group = fct_reorder(country_group, n_exposed),
    text_color = ifelse(rate_per_10000 > 15, "white", "black")
  )

# Plot
adr_rate_plot <- adr_top10_clean %>%
    mutate(
      pt = recode(pt,
                  "Increased Blood Urea Nitrogen" = "↑ Blood Urea Nitrogen"),
      rate_pct = rate_per_10000 / 100  # convert per 10 000 → %
    ) %>%
    ggplot(aes(x = country_group, y = pt, fill = rate_pct)) +
    geom_tile(color = "grey85", linewidth = 0.4) +
    geom_text(
      aes(label = sprintf("%.2f%%", rate_pct), color = text_color),
      size = 3,
      fontface = "bold"
    ) +
    scale_color_identity() +
    scale_fill_viridis_c(
      option = "plasma",
      direction = -1,
      name = "ADR rate (%)",
      labels = label_percent(accuracy = 0.01)
    ) +
    scale_x_discrete(position = "top", name = "Country or Region") + 
    scale_y_discrete(
      name = "Adverse Reaction (Preferred Term)",
      labels = function(x) str_wrap(x, width = 25)
    ) +
    labs(
      title = "Top 10 Adverse Drug Reaction Rates by Country or Region",
      subtitle = "Rates calculated using simulated exposure denominators (%)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(face = "bold", size = 11, vjust = 0.5),
      axis.text.y = element_text(face = "bold", size = 11),
      panel.grid = element_blank(),
      legend.position = "none",
      # legend.key.width = unit(2.5, "cm"),
      # legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      plot.margin = margin(10, 20, 10, 10)
    )
  
adr_rate_plot

ggsave(here("output", "figures", "adr_rate_plot.png"), 
       adr_rate_plot, width = 12, height = 8, dpi = 300, bg = "#fff1e5")
