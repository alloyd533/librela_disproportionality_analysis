# ============================================================================
# Pharmacovigilance Data Cleaning Script
# ============================================================================
# Author: Alfie Lloyd
# Date: 15/10/2025
# ============================================================================

# 1. Setup -------------------------------------------------------------------

rm(list = ls())
set.seed(2612)

required_packages <- c("tidyverse", "readxl", "here", "janitor", "scales", "gt", "jsonlite")

install_if_missing <- function(packages) {
  missing_packages <- packages[!packages %in% rownames(installed.packages())]
  if (length(missing_packages) > 0) {
    message("Installing: ", paste(missing_packages, collapse = ", "))
    install.packages(missing_packages, dependencies = TRUE)
  }
}

install_if_missing(required_packages)
invisible(lapply(required_packages, library, character.only = TRUE))

options(scipen = 999, digits = 3, tibble.print_max = 20)

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
dir.create("logs", showWarnings = FALSE, recursive = TRUE)
dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

log_file <- paste0("logs/cleaning_log_", Sys.Date(), ".txt")
cat("Cleaning started:", as.character(Sys.time()), "\n", file = log_file)

# 2. Load and Consolidate Raw Data ------------------------------------------

target_drugs <- c("bedinvetmab", "meloxicam", "carprofen", "firocoxib", 
                  "enflicoxib", "grapiprant", "robenacoxib")

raw_data_path <- here("data", "raw")
if (!dir.exists(raw_data_path)) {
  stop("Raw data directory not found: ", raw_data_path)
}

file_list <- list.files(raw_data_path, pattern = "^\\d{4}_eudra_.*\\.csv$", full.names = TRUE)
if (length(file_list) == 0) {
  stop("No yearly CSV files found in ", raw_data_path)
}

cat("Found", length(file_list), "yearly files\n")

read_and_clean_yearly <- function(file_path) {
  tryCatch({
    df <- read_csv(file_path, col_types = cols(.default = "c"), show_col_types = FALSE)
    cat("  Loaded:", basename(file_path), "(", nrow(df), "rows)\n")
    
    names(df) <- str_replace_all(names(df), " ", "_")
    
    parsed <- parse_date_time(df$Received_date, orders = c("ymd", "dmy", "ymd HMS", "dmy HMS"), quiet = TRUE)
    if (all(is.na(parsed))) {
      parsed <- suppressWarnings(as.Date(df$Received_date, format = "%Y-%m-%d"))
      if (all(is.na(parsed))) {
        parsed <- suppressWarnings(as.Date(df$Received_date, format = "%d/%m/%Y"))
      }
    }
    
    df %>%
      mutate(
        date = parsed,
        source_file = basename(file_path),
        year = str_extract(source_file, "^\\d{4}")
      ) %>%
      select(-Received_date)
    
  }, error = function(e) {
    warning("Failed to read: ", basename(file_path))
    return(NULL)
  })
}

# Process each drug
for (drug in target_drugs) {
  drug_pattern <- paste0("eudra_", drug, "\\.csv$")
  drug_files <- keep(file_list, ~ str_detect(.x, drug_pattern))
  
  if (length(drug_files) > 0) {
    cat("Processing", drug, "(", length(drug_files), "files)...\n")
    
    combined_drug_data <- map_dfr(drug_files, read_and_clean_yearly) %>%
      filter(!is.null(.))
    
    if (nrow(combined_drug_data) > 0) {
      output_path <- file.path("data", "processed", paste0(drug, "_combined.csv"))
      write_csv(combined_drug_data, output_path)
      cat("  Saved:", output_path, "\n")
    }
  }
}

# 3. Check Combined Files Exist ---------------------------------------------

cat("\nChecking for combined drug files...\n")
missing_drugs <- character()

for (drug in target_drugs) {
  drug_file <- file.path("data", "processed", paste0(drug, "_combined.csv"))
  if (!file.exists(drug_file)) {
    missing_drugs <- c(missing_drugs, drug)
  }
}

if (length(missing_drugs) > 0) {
  stop("Missing combined files for: ", paste(missing_drugs, collapse = ", "),
       "\nRun the data loading section first.")
}

cat("All drug files present\n")

# 4. Filter and Standardize -------------------------------------------------

filter_and_standardize <- function(df, drug_name) {
  initial_rows <- nrow(df)
  
  df_filtered <- df %>%
    filter(Species == "Dog") %>%
    select(-any_of(c("Case_number", "AER_form", "source_file", "Animals_treated"))) %>%
    mutate(drug = drug_name) %>%
    select(drug, year, everything())
  
  filtered_rows <- nrow(df_filtered)
  cat("  ", drug_name, ":", initial_rows, "→", filtered_rows, "rows\n")
  
  return(df_filtered)
}

all_drugs_data <- tibble()

for (drug in target_drugs) {
  drug_file <- file.path("data", "processed", paste0(drug, "_combined.csv"))
  drug_data <- read_csv(drug_file, show_col_types = FALSE) %>%
    filter_and_standardize(drug)
  all_drugs_data <- bind_rows(all_drugs_data, drug_data)
}

cat("Combined dataset:", nrow(all_drugs_data), "dog records\n")
write_csv(all_drugs_data, here("data", "processed", "complete_combined_raw.csv"))

# 5. Clean and Standardize --------------------------------------------------

all_drugs_data <- read_csv(here("data", "processed", "complete_combined_raw.csv"), show_col_types = FALSE)

complete <- all_drugs_data %>%
  rename(drug_hx = Drug) %>%
  mutate(
    death = case_when(
      is.na(Animals_died) | Animals_died == "" ~ 0,
      str_to_lower(str_trim(Animals_died)) == "no" ~ 0,
      str_to_lower(str_trim(Animals_died)) == "yes" ~ 1,
      TRUE ~ as.numeric(Animals_died)
    ),
    death = replace_na(death, 0),
    
    sex = case_when(
      str_to_upper(str_trim(Gender)) == "FEMALE" ~ "FEMALE",
      str_to_upper(str_trim(Gender)) == "MALE" ~ "MALE",
      TRUE ~ "UNKNOWN"
    ),
    
    Reaction = str_trim(Reaction),
    Country = str_to_title(str_trim(Country)),
    date = ymd(date)
  ) %>%
  select(-any_of(c("Animals_died", "Gender", "Species"))) %>%
  filter(!is.na(Reaction), Reaction != "")

cat("Data cleaning complete:", nrow(complete), "records\n")

#  Duplicate Assessment and Removal

duplicate_stats <- tibble(
  complete_duplicates = sum(duplicated(complete)),
  key_field_duplicates = sum(duplicated(complete %>% select(drug_hx, Reaction, sex, year)))
)

cat("\n=== DUPLICATE ASSESSMENT ===\n")
cat("Complete duplicates:", duplicate_stats$complete_duplicates, "\n")
cat("Key field duplicates:", duplicate_stats$key_field_duplicates, "\n")

if (duplicate_stats$complete_duplicates > 0) {
  cat("Removing", duplicate_stats$complete_duplicates, "complete duplicates\n")
  complete <- complete %>% distinct()
}

# 6. Polypharmacy Detection ------------------------------------------------

trade_names <- list(
  bedinvetmab = c("librela"),
  meloxicam = c("metacam", "loxicom", "orocam", "meloxidyl", "rheumocam", 
                "eloxioral", "melonex", "ostilox", "alloxate", "meloxivet", 
                "vivlodex", "mobic", "meloxicam", "meloxoral", "inflacam", "acticam"),
  carprofen = c("rimadyl", "novox", "vetprofen", "carprieve", "norocarp", 
                "quellin", "rovera", "carprodyl", "zinecarp", "canidryl", 
                "aventicarp", "rycarfa", "rimifin", "carpox", "tergive", 
                "levafen", "carpaquin", "carprovet", "belprofen", "movodyl", 
                "ostifen", "truprofen", "acticarp", "artriofin", "austiofen", 
                "bomazeal", "carporal", "carprocow", "carprodolor", "carprofelican", 
                "carprogesic", "carprosol", "carprotab", "carprox", "comforion", 
                "dolagis", "dolocarp", "dolox", "eurofen", "kelaprofen", 
                "librevia", "norodyl", "novocox", "prolet", "reproval", 
                "rofeniflex", "scanodyl", "xelcor", "carprofen"),
  firocoxib = c("previcox", "firox", "firocoxib", "equioxx"),
  enflicoxib = c("daxocox", "enflicoxib"),
  grapiprant = c("galliprant", "grapiprant"),
  robenacoxib = c("onsior", "robenacoxib")
)

# Pre-compile regex patterns
study_drug_patterns <- map(trade_names, ~paste0("\\b(", paste(.x, collapse = "|"), ")\\b"))

# Vectorized detection function
detect_study_drugs <- function(drug, drug_hx) {
  map2(drug, drug_hx, function(d, dh) {
    drugs_found <- character()
    for (i in seq_along(trade_names)) {
      drug_name <- names(trade_names)[i]
      if (!is.na(d) && tolower(d) == drug_name) {
        drugs_found <- c(drugs_found, drug_name)
      } else if (!is.na(dh) && str_detect(tolower(dh), regex(study_drug_patterns[[i]], ignore_case = TRUE))) {
        drugs_found <- c(drugs_found, drug_name)
      }
    }
    drugs_found
  })
}

tictoc::tic()
complete <- complete %>%
  mutate(
    study_drug_list = detect_study_drugs(drug, drug_hx),
    n_study_drugs = map_int(study_drug_list, length),
    multiple_study_drugs = as.integer(n_study_drugs > 1),
    study_drug_combination = map_chr(study_drug_list, ~if(length(.x) > 1) paste(.x, collapse = " + ") else NA_character_),
    
    all_drug_names = map(drug_hx, ~{
      if (!is.na(.x)) {
        str_trim(str_split(.x, ",")[[1]]) %>% str_extract("^[^\\s]+") %>% unique()
      } else {
        character(0)
      }
    }),
    n_total_drugs = map_int(all_drug_names, length),
    polypharmacy_binary = if_else(n_total_drugs > 1, "1", "0"),
    all_drugs_list = map_chr(all_drug_names, ~paste(.x, collapse = ", "))
  ) %>%
  mutate(
    polypharmacy = factor(polypharmacy_binary, levels = c("0", "1"), 
                          labels = c("Monotherapy", "Polypharmacy"))
  ) %>%
  select(-study_drug_list, -all_drug_names, -polypharmacy_binary)
tictoc::toc()

cat("\nPolypharmacy rates:\n")
print(complete %>%
        filter(n_study_drugs == 1) %>%
        count(polypharmacy) %>%
        mutate(percentage = round(n / sum(n) * 100, 1)))

# 7. VEDDRA Mapping ---------------------------------------------------------

veddra_path <- here("data", "combined-veddra-list-clinical-terms-reporting-suspected-adverse-events-animals-humans-veterinary-medicinal-products_en.xlsx")
if (!file.exists(veddra_path)) {
  stop("VEDDRA file not found: ", veddra_path)
}

veddra_raw <- read_excel(veddra_path)
cat("VEDDRA loaded:", nrow(veddra_raw), "terms\n")

veddra <- veddra_raw %>%
  rename_with(~ str_replace_all(.x, " ", "_")) %>%
  filter(`Current__Term_Type` != "H" | is.na(`Current__Term_Type`)) %>%
  rename_with(~ case_when(
    str_detect(.x, "Low_Level_Term.*LLT") ~ "llt",
    str_detect(.x, "Preferred_Term.*PT") ~ "pt", 
    str_detect(.x, "High_Level_Term.*HLT") ~ "hlt",
    str_detect(.x, "System_Organ_Class.*SOC") ~ "organ",
    TRUE ~ .x
  )) %>%
  select(organ, hlt, pt, llt) %>%
  mutate(llt = str_to_lower(as.character(llt))) %>%
  filter(!is.na(llt), !is.na(pt), llt != "", pt != "") %>%
  distinct()

cat("VEDDRA processed:", n_distinct(veddra$llt), "LLTs,", n_distinct(veddra$pt), "PTs\n")

# 8. Data Expansion and Reaction Mapping -----------------------------------

complete_long <- complete %>%
  mutate(dog_id = row_number()) %>%
  separate_longer_delim(Reaction, delim = ",") %>%
  mutate(
    llt = str_trim(Reaction),
    llt = str_to_lower(llt),
    llt = str_remove_all(llt, "[^a-z\\s/()\\-']"),
    llt = str_squish(llt),
    death = case_when(death == 1 ~ "Yes", death == 0 ~ "No", TRUE ~ "Unknown")
  ) %>%
  filter(!is.na(llt), llt != "", nchar(llt) >= 3) %>%
  select(-Reaction)

cat("Expanded to", nrow(complete_long), "reaction records\n")

matched <- complete_long %>%
  left_join(veddra, by = "llt", relationship = "many-to-one") %>%
  mutate(
    across(where(is.character) & !any_of("drug_hx"), as.factor),
    year = factor(year),
    death = factor(death),
    sex = case_when(sex == "FEMALE" ~ "f", sex == "MALE" ~ "m", TRUE ~ "unknown"),
    sex = factor(sex, levels = c("f", "m", "unknown"))
  ) %>%
  select(-any_of("animals_affected")) %>%
  rename_with(tolower)

cat("VEDDRA mapping rate:", round(mean(!is.na(matched$pt)) * 100, 1), "%\n")

# Save unmapped terms for review
unmapped_terms <- matched %>%
  filter(is.na(pt)) %>%
  count(llt, sort = TRUE)

if (nrow(unmapped_terms) > 0) {
  cat("Top unmapped terms:\n")
  print(head(unmapped_terms, 10))
  write_csv(unmapped_terms, here("data", "processed", "unmapped_terms_for_review.csv"))
}

# 9. Save Cleaned Data -----------------------------------------------------

write_rds(matched, here("data", "processed", "matched_data.rds"))
write_csv(matched, here("data", "processed", "matched_data.csv"))

cat("Cleaned data saved\n")

# 10. Temporal Trends -------------------------------------------------------

temporal_data <- matched %>%
  group_by(drug, year) %>%
  summarise(n_reports = n_distinct(dog_id),
            n_pts = n(),
            .groups = "drop") %>%
  mutate(
    drug_label = str_to_title(drug),
    year_numeric = as.numeric(as.character(year))
  ) %>%
  filter(!is.na(year_numeric))

# Calculate scaling factor for secondary axis
max_reports <- max(temporal_data$n_reports, na.rm = TRUE)
max_pts <- max(temporal_data$n_pts, na.rm = TRUE)
scaling_factor <- max_reports / max_pts

temporal_plot <- temporal_data %>%
  ggplot(aes(x = year_numeric, group = drug_label, color = drug_label)) +
  
  # Primary line: Number of dogs
  geom_line(aes(y = n_reports), linewidth = 1.2, alpha = 0.8) +
  geom_point(aes(y = n_reports), size = 2.5, alpha = 0.9) +
  
  # Secondary line: Number of PTs (scaled)
  geom_line(aes(y = n_pts * scaling_factor), 
            linewidth = 0.8, alpha = 0.6, linetype = "dashed") +
  geom_point(aes(y = n_pts * scaling_factor), 
             size = 1.8, alpha = 0.6, shape = 17) +
  
  scale_color_manual(
    values = c(
      "Bedinvetmab" = "#d62728",
      "Meloxicam" = "#1f77b4",
      "Carprofen" = "#ff7f0e",
      "Firocoxib" = "#2ca02c",
      "Enflicoxib" = "#9467bd",
      "Grapiprant" = "#8c564b",
      "Robenacoxib" = "#e377c2"
    )
  ) +
  
  scale_y_continuous(
    name = "Number of Reports",
    labels = comma,
    expand = expansion(mult = c(0, 0.1)),
    sec.axis = sec_axis(
      trans = ~ . / scaling_factor,
      name = "Number of PTs reported",
      labels = comma
    )
  ) +
  
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "#fff1e5", color = NA),
    panel.background = element_rect(fill = "#fff1e5", color = NA),
    panel.grid.major = element_line(color = "#f0f0f0", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    text = element_text(color = "#2d2d2d"),
    plot.title = element_text(size = 16, face = "bold", color = "#2d2d2d"),
    plot.subtitle = element_text(size = 12, color = "#555555", margin = margin(b = 20)),
    axis.title.y.left = element_text(size = 12, color = "#2d2d2d", margin = margin(r = 10)),
    axis.title.y.right = element_text(size = 12, color = "#555555", margin = margin(l = 10)),
    axis.title.x = element_text(size = 12, color = "#333333"),
    axis.text = element_text(size = 10, color = "#555555"),
    axis.text.y.right = element_text(color = "#888888"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.background = element_rect(fill = "#fff1e5", color = NA),
    plot.margin = margin(t = 20, r = 30, b = 20, l = 20)
  ) +
  
  labs(
    title = "Temporal Trends in the EudraVigilance database",
    subtitle = "Solid lines: ADR reports | Dashed lines: Preferred Terms (PTs) reported",
    x = "Year"
  ) +
  
  scale_x_continuous(breaks = seq(min(temporal_data$year_numeric, na.rm = TRUE), 
                                  max(temporal_data$year_numeric, na.rm = TRUE), by = 2))

temporal_plot

ggsave(here("output", "figures", "temporal_trends.png"), temporal_plot, 
       width = 12, height = 7, dpi = 300, bg = "#fff1e5")
# 11. Data Summary Table ----------------------------------------------------

summary_data <- matched %>%
  filter(n_study_drugs == 1) %>%
  group_by(dog_id) %>%
  slice(1) %>%
  ungroup() %>%
  select(-any_of(c("species", "drug_hx", "date", "llt", "pt", "hlt", 
                   "animals_affected", "organ", "dog_id", "drug_list", 
                   "all_drugs_list", "multiple_study_drugs", "study_drug_combination",
                   "n_total_drugs", "n_study_drugs"))) %>%
  mutate(
    sex = case_when(
      sex == "f" ~ "Female",
      sex == "m" ~ "Male",
      TRUE ~ "Unknown"
    ),
    across(c(where(is.character), where(is.factor)), 
           ~str_to_title(as.character(.x)))
  )

# Calculate PT statistics from PRE-EXPANSION data
pt_stats_raw <- complete %>%
  filter(n_study_drugs == 1) %>%
  mutate(n_reactions = str_count(Reaction, ",") + 1) %>%
  group_by(drug) %>%
  summarise(
    mean_reactions = mean(n_reactions, na.rm = TRUE),
    max_reactions = max(n_reactions, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_reactions))

data_summary_table_data <- summary_data %>%
  pivot_longer(everything(), names_to = "variable_name", values_to = "value") %>%
  group_by(variable_name) %>%
  nest() %>%
  mutate(
    summary = map2(variable_name, data, ~{
      freq <- .y %>%
        count(value, name = "n") %>%
        arrange(desc(n)) %>%
        mutate(
          percentage = n / sum(n) * 100,
          variable = paste0("  ", value)
        )
      
      n_levels <- nrow(freq)
      if (n_levels > 8) {
        freq <- freq %>%
          slice_head(n = 5) %>%
          bind_rows(tibble(
            value = NA,
            n = sum(freq$n[6:n_levels]),
            percentage = sum(freq$percentage[6:n_levels]),
            variable = paste0("  Other (", n_levels - 5, ")")
          ))
      }
      
      freq <- freq %>% select(variable, n, percentage)
      
      header <- tibble(
        variable = str_to_title(str_replace_all(.x, "_", " ")),
        n = NA_real_,
        percentage = NA_real_
      )
      
      bind_rows(header, freq)
    })
  ) %>%
  unnest(summary) %>%
  ungroup() %>%
  select(variable, n, percentage)

# Add PT statistics section with proper 3-column format
pt_header <- tibble(
  variable = "Mean number of PTs cited Per Report",
  n = NA_real_,
  percentage = NA_real_
)

pt_stats_formatted <- pt_stats_raw %>%
  mutate(
    variable = paste0("  ", str_to_title(drug)),
    n = mean_reactions,
    percentage = NA_real_
  ) %>%
  select(variable, n, percentage)

data_summary_table_data <- bind_rows(
  data_summary_table_data,
  pt_header,
  pt_stats_formatted
)

data_summary_table <- data_summary_table_data %>%
  gt() %>%
  tab_header(
    title = "EudraVigilance Database Summary",
    subtitle = paste("Total records:", format(nrow(summary_data), big.mark = ","))
  ) %>%
  cols_label(
    variable = "Variable",
    n = "Number",
    percentage = "Percentage"
  ) %>%
  fmt_number(columns = n, decimals = 0, use_seps = TRUE) %>%
  fmt_number(columns = percentage, decimals = 1) %>%
  sub_missing(columns = everything(), missing_text = "") %>%
  tab_style(
    style = cell_text(weight = "bold", color = "#1a1a1a"),
    locations = cells_body(columns = variable, rows = is.na(n))
  ) %>%
  tab_style(
    style = cell_text(color = "#333333"),
    locations = cells_body(columns = variable, rows = !is.na(n))
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#fff1e5"), cell_text(color = "#2d2d2d", size = "small")),
    locations = cells_body()
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#f2e6dd"), cell_text(weight = "bold", color = "#2d2d2d")),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#e8ddd4"), cell_text(weight = "bold", color = "#1a1a1a")),
    locations = cells_title()
  ) %>%
  cols_align(align = "left", columns = variable) %>%
  cols_align(align = "right", columns = c(n, percentage)) %>%
  tab_options(
    table.font.names = "Arial",
    table.font.size = 11,
    heading.title.font.size = 16,
    heading.subtitle.font.size = 12,
    column_labels.font.weight = "bold",
    table.border.top.style = "none",
    table.border.bottom.style = "none",
    column_labels.border.bottom.width = 2,
    column_labels.border.bottom.color = "#d62728",
    table.width = pct(100)
  ) %>%
  cols_width(
    variable ~ px(250),
    n ~ px(80),
    percentage ~ px(80)
  )

data_summary_table

gtsave(data_summary_table, here("output", "tables", "data_summary.html"))

# 12. Final Summary ---------------------------------------------------------

# Calculate on unique dogs only
dogs_summary <- matched %>%
  filter(n_study_drugs == 1) %>%
  group_by(dog_id) %>%
  slice(1) %>%
  ungroup()

final_summary <- list(
  cleaning_date = as.character(Sys.Date()),
  
  # Data coverage
  date_range = paste(min(as.numeric(as.character(matched$year)), na.rm = TRUE), "-", 
                     max(as.numeric(as.character(matched$year)), na.rm = TRUE)),
  total_dogs = n_distinct(dogs_summary$dog_id),
  total_reactions = nrow(matched %>% filter(n_study_drugs == 1)),
  
  # Drug distribution
  drugs_analyzed = length(target_drugs),
  top_drug = as.character(dogs_summary %>% count(drug, sort = TRUE) %>% slice(1) %>% pull(drug)),
  top_drug_n = dogs_summary %>% count(drug, sort = TRUE) %>% slice(1) %>% pull(n),
  
  # VEDDRA mapping
  unique_pts_identified = n_distinct(matched$pt, na.rm = TRUE),
  mapping_rate = round(sum(!is.na(matched$pt)) / nrow(matched) * 100, 1),
  
  # Clinical outcomes
  death_rate = round(mean(dogs_summary$death == "Yes", na.rm = TRUE) * 100, 1),
  
  # Polypharmacy
  polypharmacy_rate = round(mean(dogs_summary$polypharmacy == "Polypharmacy") * 100, 1),
  multiple_study_drugs_rate = round(mean(dogs_summary$multiple_study_drugs == 1) * 100, 1),
  
  # Reporting patterns
  mean_reactions_per_report = round(mean(complete %>% 
                                           filter(n_study_drugs == 1) %>% 
                                           mutate(n = str_count(Reaction, ",") + 1) %>% 
                                           pull(n), na.rm = TRUE), 1),
  max_reactions_per_report = max(complete %>% 
                                   filter(n_study_drugs == 1) %>% 
                                   mutate(n = str_count(Reaction, ",") + 1) %>% 
                                   pull(n), na.rm = TRUE)
)

write_json(final_summary, here("data", "processed", "cleaning_summary.json"), pretty = TRUE)

cat("\n=== FINAL SUMMARY ===\n")
cat("Date range:", final_summary$date_range, "\n")
cat("Total dogs (single study drug):", format(final_summary$total_dogs, big.mark = ","), "\n")
cat("Total reactions:", format(final_summary$total_reactions, big.mark = ","), "\n")
cat("Unique PTs identified:", final_summary$unique_pts_identified, "\n")
cat("VEDDRA mapping rate:", final_summary$mapping_rate, "%\n")
cat("Death rate:", final_summary$death_rate, "%\n")
cat("Polypharmacy rate:", final_summary$polypharmacy_rate, "%\n")
cat("Multiple study drugs:", final_summary$multiple_study_drugs_rate, "%\n")
cat("Mean reactions per report:", final_summary$mean_reactions_per_report, "\n")
cat("Top drug:", final_summary$top_drug, "(", format(final_summary$top_drug_n, big.mark = ","), "dogs )\n")

cat("\nCleaning completed successfully\n")
# ============================================================================
# End of Script
# ============================================================================