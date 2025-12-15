# Pharmacovigilance Signal Detection Analysis

# Methods: PRR, MGPS (EBGM), BCPNN (IC)
# Author: Alfie Lloyd
# Date: 16/10/2025


#### Section A - Main Analysis ####

# 1. Environment Setup -------------------------------------------------------
required_packages <- c(
  "tidyverse", "readxl", "gt", "scales", "openEBGM", "rlang", "here",
  "jsonlite", "glue"
)

install_if_missing <- function(packages) {
  missing_packages <- packages[!packages %in% rownames(installed.packages())]
  if (length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    install.packages(missing_packages, dependencies = TRUE)
  }
}

install_if_missing(required_packages)

# Load packages
invisible(lapply(required_packages, library, character.only = TRUE))

options(scipen = 999, digits = 3, tibble.print_max = 20)

dir.create("output", showWarnings = FALSE, recursive = TRUE)
dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

# 2. Data Loading and Validation --------------------------------------------

# 2a. Data loading
matched <- read_rds(here("data", "processed", "matched_data.rds"))

cat("\n=== Dataset Overview ===\n")
cat("Total reactions:", format(nrow(matched), big.mark = ","), "\n")
cat("Unique dogs:", format(n_distinct(matched$dog_id), big.mark = ","), "\n")
cat("PT mapping rate:", round(mean(!is.na(matched$pt)) * 100, 1), "%\n")

# 2b. Analysis Configuration

target_drug <- "bedinvetmab"
min_reports <- 3

# 2c. Create analysis dataset

dog_pt <- matched %>%
  filter(!is.na(pt), n_study_drugs == 1) %>%
  arrange(dog_id, drug, pt) %>%
  distinct(dog_id, drug, pt, .keep_all = TRUE)

denominators <- dog_pt %>%
  mutate(is_target = drug == target_drug) %>%
  distinct(dog_id, is_target) %>%
  summarise(
    n_target = sum(is_target),
    n_comparator = sum(!is_target),
    .groups = "drop"
  )

cat("\n=== Analysis Population ===\n")
cat("Target drug:", format(denominators$n_target, big.mark = ","), "dogs\n")
cat("Comparator drugs:", format(denominators$n_comparator, big.mark = ","), "dogs\n")

# 3. Proportional Reporting Ratio (PRR) ----------------------------------------

calculate_prr_comprehensive <- function(a, n_target, b, n_comparator) {
  # Convert to numeric
  a <- as.numeric(a); b <- as.numeric(b)
  n_target <- as.numeric(n_target); n_comparator <- as.numeric(n_comparator)
  
  # Check for NAs or invalid values
  if (any(is.na(c(a, b, n_target, n_comparator))) || 
      n_target <= 0 || n_comparator <= 0 || 
      a < 0 || b < 0 || a > n_target || b > n_comparator) {
    return(list(
      cell_a = NA, cell_b = NA, cell_c = NA, cell_d = NA,
      n_target_total = n_target, n_comparator_total = n_comparator,
      prop_target = NA, prop_comparator = NA, 
      prr = NA, chi_square = NA
    ))
  }
  
  c <- n_target - a; d <- n_comparator - b
  prop_target <- a / n_target; prop_comparator <- b / n_comparator
  
  # Handle division by zero
  if (prop_comparator == 0) {
    prr <- NA
  } else {
    prr <- prop_target / prop_comparator
  }
  
  n_total <- n_target + n_comparator
  expected_a <- (a + b) * n_target / n_total
  expected_b <- (a + b) * n_comparator / n_total
  expected_c <- (c + d) * n_target / n_total
  expected_d <- (c + d) * n_comparator / n_total
  
  # Handle zero expected values in chi-square
  if (any(c(expected_a, expected_b, expected_c, expected_d) == 0)) {
    chi_square <- NA
  } else {
    chi_square <- (a - expected_a)^2 / expected_a + 
      (b - expected_b)^2 / expected_b +
      (c - expected_c)^2 / expected_c +
      (d - expected_d)^2 / expected_d
  }
  
  list(
    cell_a = a, cell_b = b, cell_c = c, cell_d = d,
    n_target_total = n_target, n_comparator_total = n_comparator,
    prop_target = prop_target, prop_comparator = prop_comparator, 
    prr = prr, chi_square = chi_square
  )
}

prr_results <- dog_pt %>%
  mutate(is_target = drug == target_drug) %>%
  distinct(dog_id, is_target, pt) %>%
  group_by(pt) %>%
  summarise(
    reports_target = sum(is_target),
    reports_comparator = sum(!is_target),
    .groups = "drop"
  ) %>%
  filter(reports_target >= min_reports, reports_comparator >= min_reports) %>%
  rowwise() %>%
  mutate(
    stats = list(calculate_prr_comprehensive(reports_target, denominators$n_target, 
                                             reports_comparator, denominators$n_comparator))
  ) %>%
  unnest_wider(stats) %>%
  ungroup() %>%
  mutate(disproportionality_signal = !is.na(prr) & !is.na(chi_square) & prr >= 2 & chi_square >= 4) %>%
  arrange(desc(prr))

cat("\n=== PRR Results ===\n")
cat("PTs analyzed:", nrow(prr_results), "\n")
cat("Signals detected:", sum(prr_results$disproportionality_signal), "\n")

# 4. Multi-item Gamma Poisson Shrinker (MGPS) ----------------------------------

# 4a. Create 2x2 matrix 
mgps_data <- dog_pt %>%
  transmute(id = dog_id, var1 = drug, var2 = pt)

mgps_processed <- processRaw(data = mgps_data, stratify = FALSE, zeroes = FALSE)

# 4b. Set initialisation parameters. Use broad priors
theta_initial <- data.frame(
  alpha1 = c(1, 1, 2),
  beta1  = c(1, 2, 3),
  alpha2 = c(2, 3, 5),
  beta2  = c(2, 3, 5),
  p      = c(0.1, 0.2, 0.3)
)

# 4c. Create hyperparameters

mgps_hyperparameters <- autoHyper(
  data = mgps_processed,
  theta_init = theta_initial,
  squashed = FALSE
)

# All positive, >0 and converge + in_bounds = TRUE with finite minimum
print(mgps_hyperparameters)

# 4d. Perform MGPS
mgps_results <- ebScores(
  processed = mgps_processed,
  hyper_estimate = mgps_hyperparameters,
  quantiles = c(5, 95)
)

# 4e. Display MGPS results. 
mgps_signals <- mgps_results$data %>%
  filter(var1 == target_drug) %>%
  transmute(
    pt = var2,
    n_observed = N,
    n_expected = E,
    relative_risk = RR,
    prr_mgps = PRR,
    ebgm = EBGM,
    eb05 = QUANT_05,
    eb95 = QUANT_95,
    mgps_signal = n_observed >= min_reports & eb05 >= 2
  ) %>%
  arrange(desc(eb05))

cat("\n=== MGPS Results ===\n")
cat("Signals detected:", sum(mgps_signals$mgps_signal), "\n")

# 5. Bayesian Confidence Propagation Neural Network (BCPNN) --------------------

calculate_bcpnn_robust <- function(n11, n1p, np1, n_total, n_simulations = 50000) {
  n10 <- n1p - n11; n01 <- np1 - n11
  n00 <- n_total - n11 - n10 - n01
  
  if (any(c(n11, n10, n01, n00) < 0)) {
    return(list(ic = NA, ic_lower = NA, ic_upper = NA))
  }
  
  tryCatch({
    alpha_vec <- c(n11, n10, n01, n00) + 0.5
    gamma_samples <- matrix(
      rgamma(n_simulations * 4, shape = rep(alpha_vec, each = n_simulations)),
      nrow = n_simulations
    )
    p_samples <- gamma_samples / rowSums(gamma_samples)
    
    p11 <- p_samples[, 1]
    p_drug <- p_samples[, 1] + p_samples[, 2]
    p_event <- p_samples[, 1] + p_samples[, 3]
    
    ic_samples <- log2(p11 / (p_drug * p_event))
    ic_samples <- ic_samples[is.finite(ic_samples)]
    
    if (length(ic_samples) == 0) {
      return(list(ic = NA, ic_lower = NA, ic_upper = NA))
    }
    
    list(
      ic = mean(ic_samples),
      ic_lower = quantile(ic_samples, 0.025),
      ic_upper = quantile(ic_samples, 0.975)
    )
  }, error = function(e) {
    list(ic = NA, ic_lower = NA, ic_upper = NA)
  })
}

n_total_dogs <- n_distinct(dog_pt$dog_id)
n_target_dogs <- denominators$n_target

pt_totals <- dog_pt %>%
  group_by(pt) %>%
  summarise(np1 = n_distinct(dog_id), .groups = "drop")

pt_target <- dog_pt %>%
  filter(drug == target_drug) %>%
  group_by(pt) %>%
  summarise(n11 = n_distinct(dog_id), .groups = "drop")

bcpnn_results <- pt_totals %>%
  left_join(pt_target, by = "pt") %>%
  mutate(n11 = replace_na(n11, 0)) %>%
  filter(n11 >= min_reports) %>%
  rowwise() %>%
  mutate(
    bcpnn_stats = list(calculate_bcpnn_robust(n11, n_target_dogs, np1, n_total_dogs))
  ) %>%
  unnest_wider(bcpnn_stats) %>%
  ungroup() %>%
  mutate(bcpnn_signal = n11 >= min_reports & ic_lower > 0) %>%
  arrange(desc(ic_lower))

cat("\n=== BCPNN Results ===\n")
cat("Signals detected:", sum(bcpnn_results$bcpnn_signal, na.rm = TRUE), "\n")

# 6. Signal Integration and Consensus ------------------------------------------

integrated_results <- prr_results %>%
  select(pt, prr_signal = disproportionality_signal, prr_value = prr, chi_square) %>%
  full_join(mgps_signals %>% select(pt, mgps_signal, ebgm, eb05, eb95), by = "pt") %>%
  full_join(bcpnn_results %>% select(pt, bcpnn_signal, ic, ic_lower, ic_upper), by = "pt") %>%
  mutate(
    methods_detecting = (prr_signal %||% FALSE) + 
      (mgps_signal %||% FALSE) + 
      (bcpnn_signal %||% FALSE),
    consensus_signal = methods_detecting >= 2,
    any_method_signal = methods_detecting >= 1,
    signal_strength = case_when(
      methods_detecting >= 2 ~ 100 + (prr_value %||% 0) + (eb05 %||% 0) + 
        (pmax(0, ic_lower, na.rm = TRUE) * 10),
      methods_detecting == 1 ~ 50 + (prr_value %||% 0) + (eb05 %||% 0) + 
        (pmax(0, ic_lower, na.rm = TRUE) * 10),
      TRUE ~ 0
    )
  ) %>%
  arrange(desc(signal_strength))

priority_signals <- integrated_results %>%
  filter(consensus_signal) %>%
  arrange(desc(signal_strength))

cat("\n=== Integrated Results ===\n")
cat("Consensus signals (≥2 methods):", sum(integrated_results$consensus_signal), "\n")
cat("Any method signals:", sum(integrated_results$any_method_signal), "\n")

# 7. Main Results Tables -------------------------------------------------------

# 7a. Top PT reports with DPA
top_pts_by_reports <- dog_pt %>%
  filter(drug == target_drug) %>%
  group_by(pt) %>%
  summarise(n_reports = n_distinct(dog_id), .groups = "drop") %>%
  arrange(desc(n_reports)) %>%
  slice_head(n = 20)

results_table_data <- top_pts_by_reports %>%
  left_join(integrated_results, by = "pt") %>%
  mutate(
    prr_display = ifelse(is.na(prr_value), "—", as.character(round(prr_value, 2))),
    chi_display = ifelse(is.na(chi_square), "—", as.character(round(chi_square, 1))),
    eb05_display = ifelse(is.na(eb05), "—", as.character(round(eb05, 2))),
    ic_lower_display = ifelse(is.na(ic_lower), "—", as.character(round(ic_lower, 2))),
    has_signal = any_method_signal %||% FALSE,
    prr_significant = !is.na(prr_value) & prr_value >= 2,
    chi_significant = !is.na(chi_square) & chi_square >= 4,
    eb05_significant = !is.na(eb05) & eb05 >= 2,
    ic_significant = !is.na(ic_lower) & ic_lower > 0,
    pt_clean = str_to_title(pt)
  )

publication_table <- results_table_data %>%
  select(`Preferred Term` = pt_clean, Reports = n_reports, PRR = prr_display, 
         `Chi-square` = chi_display, EB05 = eb05_display, `IC Floor` = ic_lower_display) %>%
  gt() %>%
  # tab_header(
  #   title = paste("Top 20 Reported Adverse Reactions:", str_to_title(target_drug)),
  #   subtitle = "Ranked by frequency with signal detection results"
  # ) %>%
  tab_source_note(
    source_note = paste("N =", format(denominators$n_target, big.mark = ","), 
                        "dogs | Significant: PRR≥2, χ²≥4, EB05≥2, IC₀₂₅>0")
  ) %>%
  fmt_number(columns = Reports, decimals = 0, use_seps = TRUE) %>%
  tab_style(
    style = list(cell_fill(color = "#fff1e5"), cell_text(color = "#2d2d2d")), 
    locations = cells_body()
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#f2e6dd"), cell_text(weight = "bold", color = "#2d2d2d")), 
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"), 
    locations = cells_body(columns = `Preferred Term`, rows = which(results_table_data$has_signal))
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
    locations = cells_body(columns = PRR, rows = which(results_table_data$prr_significant))
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
    locations = cells_body(columns = `Chi-square`, rows = which(results_table_data$chi_significant))
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
    locations = cells_body(columns = EB05, rows = which(results_table_data$eb05_significant))
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
    locations = cells_body(columns = `IC Floor`, rows = which(results_table_data$ic_significant))
  ) %>%
  cols_align(align = "center", columns = c(Reports, PRR, `Chi-square`, EB05, `IC Floor`)) %>%
  cols_align(align = "left", columns = `Preferred Term`) %>%
  tab_options(table.font.names = "Arial", table.font.size = 12)

publication_table

gtsave(publication_table, here("output", "tables", "signal_detection_summary.html"))

# 7b. Most disproportionate PTs 

top_pts_by_dpa <- dog_pt %>%
  filter(drug == target_drug) %>%
  group_by(pt) %>%
  summarise(n_reports = n_distinct(dog_id), .groups = "drop")

dpa_table_data <- top_pts_by_dpa %>%
  left_join(integrated_results, by = "pt") %>%
  arrange(desc(prr_value)) %>%
  slice_head(n = 20) %>%
  mutate(
    prr_display = ifelse(is.na(prr_value), "—", as.character(round(prr_value, 2))),
    chi_display = ifelse(is.na(chi_square), "—", as.character(round(chi_square, 1))),
    eb05_display = ifelse(is.na(eb05), "—", as.character(round(eb05, 2))),
    ic_lower_display = ifelse(is.na(ic_lower), "—", as.character(round(ic_lower, 2))),
    has_signal = any_method_signal %||% FALSE,
    prr_significant = !is.na(prr_value) & prr_value >= 2,
    chi_significant = !is.na(chi_square) & chi_square >= 4,
    eb05_significant = !is.na(eb05) & eb05 >= 2,
    ic_significant = !is.na(ic_lower) & ic_lower > 0,
    pt_clean = str_to_title(pt)
  )

dpa_publication_table <- dpa_table_data %>%
  select(`Preferred Term` = pt_clean, Reports = n_reports, PRR = prr_display, 
         `Chi-square` = chi_display, EB05 = eb05_display, `IC Floor` = ic_lower_display) %>%
  gt() %>%
  # tab_header(
  #   title = paste("Most Disproportionate Adverse Reactions:", str_to_title(target_drug)),
  #   subtitle = "Ranked by Proportional Reporting Ratio"
  # ) %>%
  tab_source_note(
    source_note = paste("N =", format(denominators$n_target, big.mark = ","), 
                        "dogs | Significant: PRR≥2, χ²≥4, EB05≥2, IC₀₂₅>0")
  ) %>%
  fmt_number(columns = Reports, decimals = 0, use_seps = TRUE) %>%
  tab_style(
    style = list(cell_fill(color = "#fff1e5"), cell_text(color = "#2d2d2d")), 
    locations = cells_body()
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#f2e6dd"), cell_text(weight = "bold", color = "#2d2d2d")), 
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"), 
    locations = cells_body(columns = `Preferred Term`, rows = which(dpa_table_data$has_signal))
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
    locations = cells_body(columns = PRR, rows = which(dpa_table_data$prr_significant))
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
    locations = cells_body(columns = `Chi-square`, rows = which(dpa_table_data$chi_significant))
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
    locations = cells_body(columns = EB05, rows = which(dpa_table_data$eb05_significant))
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
    locations = cells_body(columns = `IC Floor`, rows = which(dpa_table_data$ic_significant))
  ) %>%
  cols_align(align = "center", columns = c(Reports, PRR, `Chi-square`, EB05, `IC Floor`)) %>%
  cols_align(align = "left", columns = `Preferred Term`) %>%
  tab_options(table.font.names = "Arial", table.font.size = 12)

dpa_publication_table

gtsave(dpa_publication_table, here("output", "tables", "top_20_dpa_summary.html"))

# 8. Volcano Plot for all PTs

volcano_data <- prr_results %>%
  select(pt, prr, reports_target, disproportionality_signal, chi_square) %>%
  left_join(bcpnn_results %>% select(pt, ic_lower), by = "pt") %>%
  mutate(
    log2_prr = log2(prr),
    sqrt_chi = sqrt(chi_square),
    signal_strength = case_when(
      !is.na(ic_lower) & ic_lower > 0 & disproportionality_signal ~ 
        "Strong Signal (IC>0 & PRR≥2 & χ²≥4)",
      disproportionality_signal ~ "Moderate Signal (PRR≥2 & χ²≥4)",
      !is.na(ic_lower) & ic_lower > 0 ~ "IC Signal Only",
      TRUE ~ "No Signal"
    ),
    signal_strength = factor(signal_strength, levels = c(
      "Strong Signal (IC>0 & PRR≥2 & χ²≥4)",
      "Moderate Signal (PRR≥2 & χ²≥4)",
      "IC Signal Only", "No Signal"
    )),
    n_reports = reports_target
  ) %>%
  filter(is.finite(log2_prr), is.finite(sqrt_chi))

volcano_plot <- volcano_data %>%
  ggplot(aes(x = log2_prr, y = sqrt_chi)) +
  geom_point(aes(color = signal_strength, size = n_reports), alpha = 0.7) +
  geom_hline(yintercept = sqrt(4), linetype = "dashed", color = "#767676", linewidth = 0.5) +
  geom_vline(xintercept = log2(2), linetype = "dashed", color = "#767676", linewidth = 0.5) +
  scale_color_manual(
    values = c(
      "Strong Signal (IC>0 & PRR≥2 & χ²≥4)" = "#d62728",
      "Moderate Signal (PRR≥2 & χ²≥4)" = "#ff7f0e",
      "IC Signal Only" = "#1f77b4",
      "No Signal" = "#2ca02c"
    )
  ) +
  scale_size_continuous(
    name = "Reports",
    range = c(1, 15),
    breaks = c(3, 10, 25, 50, 100)
  ) +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "#fff1e5", color = NA),
    panel.background = element_rect(fill = "#fff1e5", color = NA),
    panel.grid.major = element_line(color = "#f0f0f0", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    text = element_text(color = "#333333"),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "bottom"
  ) +
  labs(
    title = paste("Signal Detection:", str_to_title(target_drug)),
    subtitle = "Disproportionality across all Preferred Terms",
    x = expression(log[2]("PRR")),
    y = expression(sqrt(chi^2)),
    color = "Signal Thresholds"
  )

volcano_plot

ggsave(here("output", "figures", "volcano_plot.png"), volcano_plot, 
       width = 12, height = 8, dpi = 900)

#### Section B - Secondary and Sensitivity Analyses ####

# 1. Organ System Analysis

organ_prr <- dog_pt %>%
  filter(!is.na(organ)) %>%
  mutate(is_target = drug == target_drug) %>%
  distinct(dog_id, is_target, organ) %>%
  group_by(organ) %>%
  summarise(
    reports_target = sum(is_target),
    reports_comparator = sum(!is_target),
    .groups = "drop"
  ) %>%
  filter(reports_target >= min_reports, reports_comparator >= min_reports) %>%
  rowwise() %>%
  mutate(
    stats = list(calculate_prr_comprehensive(reports_target, denominators$n_target, 
                                             reports_comparator, denominators$n_comparator))
  ) %>%
  unnest_wider(stats) %>%
  ungroup() %>%
  mutate(prr_signal = !is.na(prr) & !is.na(chi_square) & prr >= 2 & chi_square >= 4)

organ_bcpnn <- dog_pt %>%
  filter(!is.na(organ)) %>%
  mutate(is_target = drug == target_drug) %>%
  distinct(dog_id, is_target, organ) %>%
  group_by(organ) %>%
  summarise(
    n11 = sum(is_target),
    n1p = denominators$n_target,
    np1 = n(),
    .groups = "drop"
  ) %>%
  filter(n11 >= min_reports) %>%
  rowwise() %>%
  mutate(
    bcpnn_stats = list(calculate_bcpnn_robust(n11, n1p, np1, denominators$n_target + denominators$n_comparator))
  ) %>%
  unnest_wider(bcpnn_stats) %>%
  ungroup() %>%
  mutate(bcpnn_signal = n11 >= min_reports & ic_lower > 0)

organ_integrated <- organ_prr %>%
  full_join(organ_bcpnn %>% select(organ, ic, ic_lower, bcpnn_signal), by = "organ") %>%
  mutate(
    methods_detecting = (prr_signal %||% FALSE) + (bcpnn_signal %||% FALSE),
    any_signal = methods_detecting >= 1
  ) %>%
  arrange(desc(methods_detecting), desc(prr))

cat("\n=== Organ System Analysis ===\n")
cat("Systems analyzed:", nrow(organ_integrated), "\n")
cat("Systems with signals:", sum(organ_integrated$any_signal, na.rm = TRUE), "\n")

organ_table <- organ_integrated %>%
  arrange(desc(prr)) %>%
  slice_head(n = 10) %>%
  select(organ, reports_target, prr, chi_square, ic_lower) %>%
  gt() %>%
  # tab_header(
  #   title = md(paste0("**Organ System Analysis**")),
  #   subtitle = md("Disproportionality by Organ System - 10 highest PRR")
  # ) %>%
  cols_label(
    organ = "Organ System",
    reports_target = "Reports",
    prr = "PRR",
    chi_square = md("χ²"),
    ic_lower = "IC Floor"
  ) %>%
  fmt_number(columns = c(prr, chi_square, ic_lower), decimals = 2) %>%
  fmt_integer(columns = reports_target) %>%
  tab_style(
    style = cell_text(weight = "bold"), 
    locations = cells_body(rows = prr>2 & ic_lower>0)
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#fff1e5"), cell_text(color = "#2d2d2d")), 
    locations = cells_body()
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#f2e6dd"), cell_text(weight = "bold", color = "#2d2d2d")), 
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
    locations = cells_body(columns = prr, rows = !is.na(prr) & prr >= 2)
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
    locations = cells_body(columns = chi_square, rows = !is.na(chi_square) & chi_square >= 4)
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
    locations = cells_body(columns = ic_lower, rows = !is.na(ic_lower) & ic_lower > 0)
  ) %>%
  cols_align(align = "center", columns = c(reports_target, prr, chi_square, ic_lower)) %>%
  cols_align(align = "left", columns = organ) %>%
  tab_options(table.font.names = "Arial", table.font.size = 12)

organ_table

gtsave(organ_table, here("output", "tables", "organ_system_analysis.html"))

# 2. Monotherapy Sensitivity Analysis

dog_pt_mono <- dog_pt %>% filter(polypharmacy == "Monotherapy")

mono_denominators <- dog_pt_mono %>%
  mutate(is_target = drug == target_drug) %>%
  distinct(dog_id, is_target) %>%
  summarise(
    n_target = sum(is_target),
    n_comparator = sum(!is_target),
    .groups = "drop"
  )

prr_mono_results <- dog_pt_mono %>%
  mutate(is_target = drug == target_drug) %>%
  distinct(dog_id, is_target, pt) %>%
  group_by(pt) %>%
  summarise(
    reports_target = sum(is_target),
    reports_comparator = sum(!is_target),
    .groups = "drop"
  ) %>%
  filter(reports_target >= min_reports, reports_comparator >= min_reports) %>%
  rowwise() %>%
  mutate(
    stats = list(calculate_prr_comprehensive(reports_target, mono_denominators$n_target, 
                                             reports_comparator, mono_denominators$n_comparator))
  ) %>%
  unnest_wider(stats) %>%
  ungroup() %>%
  mutate(disproportionality_signal = !is.na(prr) & !is.na(chi_square) & prr >= 2 & chi_square >= 4) %>%
  arrange(desc(prr))

cat("\n=== Monotherapy Analysis ===\n")
cat("Target drug:", format(mono_denominators$n_target, big.mark = ","), "dogs\n")
cat("Signals detected:", sum(prr_mono_results$disproportionality_signal), "\n")

# Top 10 most prevalent
top_prevalent <- prr_results %>%
  arrange(desc(reports_target)) %>%
  slice_head(n = 10) %>%
  select(pt, reports_main = reports_target, prr_main = prr, chi_main = chi_square) %>%
  left_join(
    prr_mono_results %>% select(pt, reports_mono = reports_target, prr_mono = prr, chi_mono = chi_square),
    by = "pt"
  ) %>%
  filter(!is.na(prr_mono)) %>%
  mutate(section = "Most Commonly Reported PTs")

# Top 10 most disproportionate
top_disproportionate <- prr_results %>%
  filter(disproportionality_signal) %>%
  arrange(desc(prr)) %>%
  slice_head(n = 10) %>%
  select(pt, reports_main = reports_target, prr_main = prr, chi_main = chi_square) %>%
  left_join(
    prr_mono_results %>% select(pt, reports_mono = reports_target, prr_mono = prr, chi_mono = chi_square, reports_target),
    by = "pt"
  ) %>%
  filter(!is.na(prr_mono)) %>%
  mutate(section = "Most Disproportionate PT")

monotherapy_comparison <- bind_rows(top_prevalent, top_disproportionate)

mono_comparison_table <- monotherapy_comparison %>%
  gt(groupname_col = "section") %>%
  # tab_header(
  #   title = md(paste0("**Sensitivity Analysis Of Dogs Only On Bedinvetmab**")),
  #   subtitle = md("Main Analysis vs Monotherapy")
  # ) %>%
  cols_label(
    pt = "Preferred Term",
    reports_main = "Reports",
    prr_main = "PRR",
    chi_main = md("χ²"),
    reports_mono = "Reports",
    prr_mono = "PRR",
    chi_mono = md("χ²")
  ) %>%
  tab_spanner(
    label = "Main Analysis (n = {comma(denominators$n_target)})",
    columns = c(reports_main, prr_main, chi_main)
  ) %>%
  tab_spanner(
    label = "Monotherapy (n = {comma(mono_denominators$n_target)})",
    columns = c(reports_mono, prr_mono, chi_mono)
  ) %>%
  fmt_number(columns = c(prr_main, prr_mono, chi_main, chi_mono), decimals = 2) %>%
  fmt_integer(columns = c(reports_main, reports_mono)) %>%
  tab_style(
    style = list(cell_text(weight = "bold", align = "center")),
    locations = cells_row_groups()
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#fff1e5"), cell_text(color = "#2d2d2d")), 
    locations = cells_body()
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#f2e6dd"), cell_text(weight = "bold", color = "#2d2d2d")), 
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")),
    locations = cells_body(columns = prr_main, rows = !is.na(prr_main) & prr_main >= 2)
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")),
    locations = cells_body(columns = prr_mono, rows = !is.na(prr_mono) & prr_mono >= 2)
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")),
    locations = cells_body(columns = chi_main, rows = !is.na(chi_main) & chi_main >= 4)
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")),
    locations = cells_body(columns = chi_mono, rows = !is.na(chi_mono) & chi_mono >= 4)
  ) %>%
  cols_align(align = "center", columns = c(reports_main, prr_main, chi_main, reports_mono, prr_mono, chi_mono)) %>%
  cols_align(align = "left", columns = pt) %>%
  tab_options(table.font.names = "Arial", table.font.size = 12)

mono_comparison_table

gtsave(mono_comparison_table, here("output", "tables", "monotherapy_comparison.html"))

# 3. NSAID subgroup analysis

dog_pt_nsaid <- dog_pt %>% filter(drug != "grapiprant")

nsaid_denominators <- dog_pt_nsaid %>%
  mutate(is_target = drug == target_drug) %>%
  distinct(dog_id, is_target) %>%
  summarise(
    n_target = sum(is_target),
    n_comparator = sum(!is_target),
    .groups = "drop"
  )

prr_nsaid_results <- dog_pt_nsaid %>%
  mutate(is_target = drug == target_drug) %>%
  distinct(dog_id, is_target, pt) %>%
  group_by(pt) %>%
  summarise(
    reports_target = sum(is_target),
    reports_comparator = sum(!is_target),
    .groups = "drop"
  ) %>%
  filter(reports_target >= min_reports, reports_comparator >= min_reports) %>%
  rowwise() %>%
  mutate(
    stats = list(calculate_prr_comprehensive(reports_target, nsaid_denominators$n_target, 
                                             reports_comparator, nsaid_denominators$n_comparator))
  ) %>%
  unnest_wider(stats) %>%
  ungroup() %>%
  mutate(disproportionality_signal = !is.na(prr) & !is.na(chi_square) & prr >= 2 & chi_square >= 4) %>%
  arrange(desc(prr))

cat("\n=== NSAID Subgroup Analysis ===\n")
cat("Target drug:", format(nsaid_denominators$n_target, big.mark = ","), "dogs\n")
cat("Signals detected:", sum(prr_nsaid_results$disproportionality_signal), "\n")

# Top 10 most prevalent
top_prevalent <- prr_results %>%
  arrange(desc(reports_target)) %>%
  slice_head(n = 10) %>%
  select(pt, reports_main = reports_target, prr_main = prr, chi_main = chi_square) %>%
  left_join(
    prr_nsaid_results %>% select(pt, reports_nsaid = reports_target, prr_nsaid = prr, chi_nsaid = chi_square),
    by = "pt"
  ) %>%
  filter(!is.na(prr_nsaid)) %>%
  mutate(section = "Most Commonly Reported PTs")

# Top 10 most disproportionate
top_disproportionate <- prr_results %>%
  filter(disproportionality_signal) %>%
  arrange(desc(prr)) %>%
  slice_head(n = 10) %>%
  select(pt, reports_main = reports_target, prr_main = prr, chi_main = chi_square) %>%
  left_join(
    prr_nsaid_results %>% select(pt, reports_nsaid = reports_target, prr_nsaid = prr, chi_nsaid = chi_square),
    by = "pt"
  ) %>%
  filter(!is.na(prr_nsaid)) %>%
  mutate(section = "Most Disproportionate PT")

nsaid_comparison <- bind_rows(top_prevalent, top_disproportionate)

nsaid_comparison_table <- nsaid_comparison %>%
  gt(groupname_col = "section") %>%
  tab_header(
    title = md(paste0("**Sensitivity Analysis Of Dogs Only On Bedinvetmab**")),
    subtitle = md("Main Analysis vs NSAID-only")
  ) %>%
  cols_label(
    pt = "Preferred Term",
    reports_main = "Reports",
    prr_main = "PRR",
    chi_main = md("χ²"),
    reports_nsaid = "Reports",
    prr_nsaid = "PRR",
    chi_nsaid = md("χ²")
  ) %>%
  tab_spanner(
    label = "Main Analysis (n = {comma(denominators$n_target)})",
    columns = c(reports_main, prr_main, chi_main)
  ) %>%
  tab_spanner(
    label = "NSAID (n = {comma(nsaid_denominators$n_comparator)})",
    columns = c(reports_nsaid, prr_nsaid, chi_nsaid)
  ) %>%
  fmt_number(columns = c(prr_main, prr_nsaid, chi_main, chi_nsaid), decimals = 2) %>%
  fmt_integer(columns = c(reports_main, reports_nsaid)) %>%
  tab_style(
    style = list(cell_text(weight = "bold", align = "center")),
    locations = cells_row_groups()
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")),
    locations = cells_body(columns = prr_main, rows = !is.na(prr_main) & prr_main >= 2)
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")),
    locations = cells_body(columns = prr_nsaid, rows = !is.na(prr_nsaid) & prr_nsaid >= 2)
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")),
    locations = cells_body(columns = chi_main, rows = !is.na(chi_main) & chi_main >= 4)
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")),
    locations = cells_body(columns = chi_nsaid, rows = !is.na(chi_nsaid) & chi_nsaid >= 4)
  ) %>%
  cols_align(align = "center", columns = c(reports_main, prr_main, chi_main, reports_nsaid, prr_nsaid, chi_nsaid)) %>%
  cols_align(align = "left", columns = pt) %>%
  tab_options(table.font.names = "Arial", table.font.size = 12)

nsaid_comparison_table

gtsave(nsaid_comparison_table, here("output", "tables", "nsaid_comparison.html"))
# 13. Executive Summary

executive_summary <- list(
  analysis_date = as.character(Sys.Date()),
  target_drug = target_drug,
  total_dogs_analyzed = denominators$n_target + denominators$n_comparator,
  target_drug_dogs = denominators$n_target,
  comparator_drug_dogs = denominators$n_comparator,
  pts_analyzed = nrow(integrated_results),
  prr_signals = sum(integrated_results$prr_signal, na.rm = TRUE),
  mgps_signals = sum(integrated_results$mgps_signal, na.rm = TRUE),
  bcpnn_signals = sum(integrated_results$bcpnn_signal, na.rm = TRUE),
  consensus_signals = sum(integrated_results$consensus_signal, na.rm = TRUE),
  top_consensus_signals = if(nrow(priority_signals) > 0) {
    priority_signals$pt[1:min(5, nrow(priority_signals))]
  } else character(0)
)

write_json(executive_summary, here("output", "executive_summary.json"), pretty = TRUE)

cat("\n=== Analysis Complete ===\n")
cat("Consensus signals:", executive_summary$consensus_signals, "\n")
cat("Results saved to: output/\n")


# End of Analysis ==============================================================