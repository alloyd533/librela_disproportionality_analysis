# ===============================================================================
# COMPREHENSIVE PHARMACOVIGILANCE DATA CLEANING SCRIPT
# ===============================================================================
# 
# PURPOSE:
# This script processes raw adverse drug reaction data from multiple yearly files,
# performs comprehensive data cleaning, validation, and VEDDRA terminology mapping.
# Output: Clean 'matched' dataframe ready for statistical signal detection analysis.
#
# CLINICAL CONTEXT:
# - Ensures data quality before statistical analysis
# - Standardizes adverse reaction terminology using VEDDRA
# - Handles multiple data sources and formats consistently
# - Provides comprehensive data validation and quality reporting
#
# AUTHOR: [Your Name]
# DATE: [Current Date]
# VERSION: 2.0
# ===============================================================================

# ===============================================================================
# 1. ENVIRONMENT SETUP AND DEPENDENCIES
# ===============================================================================

# Clear workspace and set reproducible seed
rm(list = ls())
set.seed(42)

# Load required packages with version tracking
required_packages <- c(
  "tidyverse",    # Data manipulation and visualization
  "readxl",       # Excel file reading
  "lubridate",    # Date handling
  "here",         # Robust file path management
  "janitor",      # Data cleaning utilities
  "scales",       # Number formatting
  "gt",           # Publication-quality tables
  "jsonlite"
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

# Load all packages
invisible(lapply(required_packages, library, character.only = TRUE))

# Set global options
options(
  scipen = 999,
  digits = 3,
  tibble.print_max = 20,
  tibble.print_min = 10
)

# Create output directories
dir.create("data", showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
dir.create("logs", showWarnings = FALSE, recursive = TRUE)

# Set up logging
log_file <- paste0("logs/cleaning_log_", Sys.Date(), ".txt")
cat("Data cleaning started at:", as.character(Sys.time()), "\n", file = log_file)

message("✓ Environment setup complete")

# ===============================================================================
# 2. RAW DATA LOADING AND CONSOLIDATION
# ===============================================================================

message("📁 Loading and consolidating raw data files...")

# Define target drugs for analysis
target_drugs <- c("bedinvetmab", "meloxicam", "carprofen", "firocoxib", 
                  "enflicoxib", "grapiprant", "robenacoxib")

# Check if raw data directory exists
raw_data_path <- here("data", "raw")
if (!dir.exists(raw_data_path)) {
  stop("❌ Raw data directory not found at: ", raw_data_path,
       "\nPlease create the directory and place yearly CSV files there.")
}

# List all yearly CSV files with improved pattern matching
file_list <- list.files(raw_data_path, 
                        pattern = "^\\d{4}_eudra_.*\\.csv$", 
                        full.names = TRUE)

if (length(file_list) == 0) {
  stop("❌ No yearly CSV files found matching pattern '^\\d{4}_eudra_.*\\.csv$' in ", raw_data_path)
}

message("✓ Found ", length(file_list), " yearly data files")

# Enhanced function to read and clean individual CSV files (simplified)
read_and_clean_yearly <- function(file_path) {
  
  tryCatch({
    # Read with all columns as character to prevent type conflicts
    df <- read_csv(file_path, col_types = cols(.default = "c"), show_col_types = FALSE)
    
    message("  ✓ Loaded: ", basename(file_path), " (", nrow(df), " rows)")
    
    # Standardize column names (replace spaces with underscores)
    names(df) <- str_replace_all(names(df), " ", "_")
    
    # Parse dates - simple approach from original script
    parsed <- parse_date_time(df$Received_date, orders = c("ymd", "dmy", "ymd HMS", "dmy HMS"), quiet = TRUE)
    if (all(is.na(parsed))) {
      parsed <- suppressWarnings(as.Date(df$Received_date, format = "%Y-%m-%d"))
      if (all(is.na(parsed))) {
        parsed <- suppressWarnings(as.Date(df$Received_date, format = "%d/%m/%Y"))
      }
    }
    
    # Add metadata and clean up
    df <- df %>%
      mutate(
        date = parsed,
        source_file = basename(file_path),
        year = str_extract(source_file, "^\\d{4}")
      ) %>%
      select(-Received_date)
    
    return(df)
    
  }, error = function(e) {
    warning("❌ Failed to read file: ", basename(file_path), " - ", e$message)
    return(NULL)
  })
}

# Process and save cleaned yearly data per drug
drug_files_processed <- 0

for (drug in target_drugs) {
  drug_pattern <- paste0("eudra_", drug, "\\.csv$")
  drug_files <- keep(file_list, ~ str_detect(.x, drug_pattern))
  
  if (length(drug_files) > 0) {
    message("Processing ", drug, " (", length(drug_files), " files)...")
    
    # Read and combine all years for this drug
    combined_drug_data <- map_dfr(drug_files, read_and_clean_yearly) %>%
      filter(!is.null(.)) # Remove any NULL results from failed reads
    
    if (nrow(combined_drug_data) > 0) {
      # Save individual drug dataset
      output_path <- file.path("data", "processed", paste0(drug, "_combined.csv"))
      write_csv(combined_drug_data, output_path)
      drug_files_processed <- drug_files_processed + 1
      
      message("  ✓ Saved: ", output_path, " (", nrow(combined_drug_data), " records)")
    }
  } else {
    warning("⚠️  No files found for drug: ", drug)
  }
}

message("✓ Processed ", drug_files_processed, " drugs successfully")

# ===============================================================================
# 3. DATA FILTERING AND STANDARDIZATION
# ===============================================================================

message("🔄 Filtering and standardizing consolidated data...")

# Simple filtering function (based on original approach)
filter_and_standardize <- function(df, drug_name) {
  
  initial_rows <- nrow(df)
  
  # Simple filtering like original script
  df_filtered <- df %>%
    filter(Species == "Dog") %>%
    select(-any_of(c("Case_number", "AER_form", "source_file", "Animals_treated"))) %>%
    mutate(drug = drug_name) %>%
    select(drug, year, everything())
  
  filtered_rows <- nrow(df_filtered)
  
  message("  ", drug_name, ": ", initial_rows, " → ", filtered_rows, " rows (", 
          initial_rows - filtered_rows, " non-dog records removed)")
  
  return(df_filtered)
}

# Load, filter, and combine all drug datasets
all_drugs_data <- tibble()

for (drug in target_drugs) {
  drug_file <- file.path("data", "processed", paste0(drug, "_combined.csv"))
  
  if (file.exists(drug_file)) {
    drug_data <- read_csv(drug_file, show_col_types = FALSE) %>%
      filter_and_standardize(drug)
    
    all_drugs_data <- bind_rows(all_drugs_data, drug_data)
  } else {
    warning("⚠️  Processed file not found: ", drug_file)
  }
}

message("✓ Combined dataset: ", nrow(all_drugs_data), " total dog records")

# Save initial combined dataset
write_csv(all_drugs_data, here("data", "processed", "complete_combined_raw.csv"))

# ===============================================================================
# 4. COMPREHENSIVE DATA QUALITY ASSESSMENT AND CLEANING
# ===============================================================================


message("🧹 Performing data cleaning and quality assessment...")
all_drugs_data <- read_csv(here("data", "processed", "complete_combined_raw.csv"))

# Clean death and sex variables AND standardize other fields
complete <- all_drugs_data %>%
  rename(drug_hx = Drug) %>%
  # Clean all essential variables and rename to final names
  mutate(
    # Handle death reporting: missing values mean no deaths occurred
    death = case_when(
      is.na(Animals_died) | Animals_died == "" ~ 0,
      str_to_lower(str_trim(Animals_died)) == "no" ~ 0,
      str_to_lower(str_trim(Animals_died)) == "yes" ~ 1,
      TRUE ~ as.numeric(Animals_died)
    ),
    death = replace_na(death, 0),
    
    # Standardize sex coding - anything not FEMALE/MALE becomes UNKNOWN  
    sex = case_when(
      str_to_upper(str_trim(Gender)) == "FEMALE" ~ "FEMALE",
      str_to_upper(str_trim(Gender)) == "MALE" ~ "MALE",
      TRUE ~ "UNKNOWN"
    ),
    
    # Clean reaction field
    Reaction = str_trim(Reaction),
    
    Country = str_to_title(str_trim(Country)),
    date = ymd(date)
  ) %>%
  # Remove original columns that were renamed
  select(-any_of(c("Animals_died", "Gender"))) %>%
  # Remove completely empty reaction records
  filter(!is.na(Reaction), Reaction != "")

message("✓ Data cleaning complete - ", nrow(complete), " records retained")

message("🔍 Performing comprehensive data quality assessment...")

# Create enhanced data quality report
generate_quality_report <- function(df) {
  
  # Basic structure metrics
  basic_stats <- tibble(
    metric = c("Total records", "Total columns", "Unique drugs", "Date range"),
    value = c(
      format(nrow(df), big.mark = ","),
      ncol(df),
      n_distinct(df$drug_hx, na.rm = TRUE),
      if("year" %in% names(df)) {
        paste(min(df$year, na.rm = TRUE), "to", max(df$year, na.rm = TRUE))
      } else "Unknown"
    )
  )
  
  # Missing data assessment
  missing_analysis <- df %>%
    summarise(across(everything(), ~sum(is.na(.x) | .x == ""))) %>%
    pivot_longer(everything(), names_to = "column", values_to = "missing_count") %>%
    mutate(
      missing_percent = round(missing_count / nrow(df) * 100, 1),
      data_quality = case_when(
        missing_percent == 0 ~ "Excellent",
        missing_percent <= 5 ~ "Good",
        missing_percent <= 15 ~ "Acceptable", 
        missing_percent <= 30 ~ "Poor",
        TRUE ~ "Critical"
      )
    ) %>%
    arrange(desc(missing_percent))
  
  # Drug distribution
  drug_distribution <- df %>%
    count(drug_hx, sort = TRUE) %>%
    mutate(percentage = round(n / sum(n) * 100, 1))
  
  # Year distribution
  year_distribution <- df %>%
    count(year, sort = TRUE) %>%
    mutate(percentage = round(n / sum(n) * 100, 1))
  
  # Gender distribution
  sex_distribution <- df %>%
    count(sex, sort = TRUE) %>%
    mutate(percentage = round(n / sum(n) * 100, 1))
  
  # Duplicate assessment
  duplicate_stats <- tibble(
    complete_duplicates = sum(duplicated(df)),
    key_field_duplicates = sum(duplicated(df %>% select(drug_hx, Reaction, sex, year)))
  )
  
  list(
    basic = basic_stats,
    missing = missing_analysis,
    drugs = drug_distribution,
    years = year_distribution,
    sex = sex_distribution,
    duplicates = duplicate_stats
  )
}

# Generate comprehensive quality report
quality_report <- generate_quality_report(complete)

# Display key quality metrics
cat("\n=== DATA QUALITY ASSESSMENT ===\n")
print(quality_report$basic)

cat("\n=== MISSING DATA ANALYSIS ===\n")
high_missing <- quality_report$missing %>% filter(missing_percent > 10)
if (nrow(high_missing) > 0) {
  print(high_missing)
} else {
  cat("✓ No columns with >10% missing data\n")
}

cat("\n=== DRUG DISTRIBUTION ===\n")
print(quality_report$drugs)

cat("\n=== DUPLICATE ASSESSMENT ===\n")
print(quality_report$duplicates)

# Remove complete duplicates if found
if (quality_report$duplicates$complete_duplicates > 0) {
  message("⚠️  Removing ", quality_report$duplicates$complete_duplicates, " complete duplicates")
  complete <- complete %>% distinct()
}

# Log quality metrics
cat("Data Quality Summary:\n", file = log_file, append = TRUE)
cat(paste("- Total records:", nrow(complete), "\n"), file = log_file, append = TRUE)
cat(paste("- Complete duplicates removed:", quality_report$duplicates$complete_duplicates, "\n"), 
    file = log_file, append = TRUE)

message("✓ Data quality assessment complete")

# ===============================================================================
# 5. ENHANCED POLYPHARMACY DETECTION
# ===============================================================================

message("💊 Implementing enhanced polypharmacy detection...")

# # Advanced polypharmacy detection function
# detect_polypharmacy_advanced <- function(drug_string) {
#   
#   if (is.na(drug_string) || str_trim(drug_string) == "") return("0")
#   
#   # Handle multiple possible delimiters and conjunctions
#   delimiters <- c(",", ";", "\\+", "&", " and ", " AND ", " with ", " WITH ")
#   
#   # Split by any delimiter
#   drug_list <- drug_string
#   for (delim in delimiters) {
#     drug_list <- unlist(strsplit(drug_list, delim, perl = TRUE))
#   }
#   
#   # Clean and extract meaningful drug names
#   drug_names <- drug_list %>%
#     str_trim() %>%                                    # Remove whitespace
#     str_to_lower() %>%                                # Standardize case
#     str_remove_all("\\(.*\\)") %>%                   # Remove parenthetical info
#     str_remove_all("\\d+\\s*(mg|ml|g|%).*$") %>%     # Remove dosing information
#     str_remove_all("\\b(tablet|capsule|injection|oral)\\b") %>%  # Remove formulation terms
#     str_extract("^[a-z][a-z0-9]*") %>%               # Extract drug name (letters + numbers)
#     .[!is.na(.) & . != "" & nchar(.) >= 3]          # Remove short/empty entries
#   
#   # Count unique meaningful drug names
#   unique_drugs <- unique(drug_names)
#   
#   polypharmacy <- if (length(unique_drugs) > 1) "1" else "0"
#   
#   # Return both indicator and cleaned list
#   return(list(
#     polypharmacy = polypharmacy,
#     drug_list = paste(unique_drugs, collapse = ", ")
#   ))
# }

# Apply polypharmacy detection with progress tracking
message("  Analyzing ", nrow(complete), " drug history records...")

# complete <- complete %>%
#   mutate(out = map(drug_hx, detect_polypharmacy_advanced)) %>%
#   unnest_wider(out, names_sep = "_") %>%                # creates out_polypharmacy, out_drug_list
#   mutate(polypharmacy = factor(out_polypharmacy,
#                                levels = c("0","1"),
#                                labels = c("Monotherapy","Polypharmacy")),
#          drug_list = out_drug_list) %>%
#   select(-out_polypharmacy, -out_drug_list)

complete <- complete |>
  mutate(
    out = drug_hx |>
      str_split(",") |>
      map(~ {
        drug_names <- str_trim(.x) |> str_extract("^[^\\s]+")
        tibble(
          polypharmacy = if (length(unique(drug_names)) > 1) "1" else "0",
          drug_list    = paste(unique(drug_names), collapse = ", ")
        )
      })
  ) |>
  unnest_wider(out) |>
  mutate(
    polypharmacy = factor(polypharmacy,
                          levels = c("0","1"),
                          labels = c("Monotherapy","Polypharmacy"))
  )

# Generate polypharmacy summary
cat("\n=== POLYPHARMACY ASSESSMENT ===\n")
print(complete %>%
  count(polypharmacy) %>%
  mutate(percentage = round(n / sum(n) * 100, 1)))

message("✓ Polypharmacy detection complete")

# ===============================================================================
# 6. VEDDRA TERMINOLOGY MAPPING WITH VALIDATION
# ===============================================================================

message("📚 Loading and validating VEDDRA terminology...")

# Load VEDDRA terminology with comprehensive error handling
veddra_path <- here("data", "combined-veddra-list-clinical-terms-reporting-suspected-adverse-events-animals-humans-veterinary-medicinal-products_en.xlsx")

if (!file.exists(veddra_path)) {
  stop("❌ VEDDRA terminology file not found at: ", veddra_path,
       "\nPlease download from EMA website and place in data/ folder.")
}

# Load and clean VEDDRA
tryCatch({
  veddra_raw <- read_excel(veddra_path)
  message("✓ VEDDRA terminology loaded: ", nrow(veddra_raw), " terms")
}, error = function(e) {
  stop("❌ Error loading VEDDRA file: ", e$message)
})

# Process VEDDRA terminology
veddra <- veddra_raw %>%
  # Standardize column names
  rename_with(~ str_replace_all(.x, " ", "_")) %>%
  # Filter for veterinary terms (exclude human-only)
  filter(`Current__Term_Type` != "H" | is.na(`Current__Term_Type`)) %>%
  # Rename key columns with meaningful names
  rename_with(~ case_when(
    str_detect(.x, "Low_Level_Term.*LLT") ~ "llt",
    str_detect(.x, "Preferred_Term.*PT") ~ "pt", 
    str_detect(.x, "High_Level_Term.*HLT") ~ "hlt",
    str_detect(.x, "System_Organ_Class.*SOC") ~ "organ",
    TRUE ~ .x
  )) %>%
  # Select only essential columns
  select(organ, hlt, pt, llt) %>%
  # Clean and standardize LLT for matching
  mutate(llt = str_to_lower(as.character(llt))) %>%
  # Remove any entries with missing critical information
  filter(!is.na(llt), !is.na(pt), llt != "", pt != "") %>%
  # Remove duplicates
  distinct()

# Validate VEDDRA structure
veddra_validation <- tibble(
  total_terms = nrow(veddra),
  unique_llt = n_distinct(veddra$llt),
  unique_pt = n_distinct(veddra$pt),  
  unique_hlt = n_distinct(veddra$hlt),
  unique_organs = n_distinct(veddra$organ),
  missing_pt_mappings = sum(is.na(veddra$pt))
)

cat("\n=== VEDDRA VALIDATION ===\n")
print(veddra_validation)

message("✓ VEDDRA terminology processed and validated")

veddra_plot <- tibble(
  level = factor(c("Organ System","High-Level Terms","Preferred Terms","Lower-Level Terms"),
                 levels = c("Organ System","High-Level Terms","Preferred Terms","Lower-Level Terms")),
  n = c(24, 231, 1090, 2777)
)

# heights proportional to counts; stack from top (Organ System) to bottom
df <- veddra_plot %>%
  mutate(
    width = n / sum(n),
    y_bot = row_number() - 1,         
    y_top = row_number(),             
    # CENTER rectangles around x = 0
    x_min = -width / 2,               # Left edge
    x_max = width / 2,                # Right edge
    x_mid = 0,                        # Center at x = 0
    y_mid = (y_bot + y_top) / 2,
    label = paste0(level, "\n", scales::comma(n))  # Add comma formatting
  )

# Financial Times theme
ft_bg   <- "#fff1e5"   
ft_grid <- "#f2e6dd"
ft_cols <- c("#0f5499", "#2e6ea6", "#5b8fc2", "#93b7d9")  

theme_ft <- function() {
  theme_minimal(base_size = 13, base_family = "Arial") +
    theme(
      plot.background  = element_rect(fill = ft_bg, colour = NA),
      panel.background = element_rect(fill = ft_bg, colour = NA),
      panel.grid.major = element_line(colour = ft_grid, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.title       = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      legend.position  = "none",
      plot.title       = element_text(face = "bold", colour = "#2d2d2d", hjust = 0, size = 16),
      plot.subtitle    = element_text(colour = "#555555", hjust = 0, size = 12),
      plot.caption     = element_text(colour = "#888888", hjust = 0, size = 10),
      plot.margin      = margin(20, 20, 20, 20)
    )
}

# Calculate symmetric limits
x_lim <- max(df$width) / 2 * 1.1

# Create the pyramid plot
veddra_pyramid <- ggplot(df) +
  geom_rect(aes(xmin = x_min, xmax = x_max, ymin = y_bot, ymax = y_top, fill = level),
            colour = "white", linewidth = 1.2, alpha = 0.9) +
  
  geom_text(aes(x = x_mid, y = y_mid, label = label), 
            fontface = "bold", 
            lineheight = 0.9, 
            colour = "#1a1a1a",
            size = 4) +
  
  scale_fill_manual(values = ft_cols) +
  scale_y_reverse(expand = c(0, 0.1), limits = c(nrow(df) + 0.1, -0.1)) +
  scale_x_continuous(expand = c(0, 0), limits = c(-x_lim, x_lim)) +
  
  labs(
    title = "VeDDRA Terminology Hierarchy",
    subtitle = "Distribution of veterinary adverse reaction terms across classification levels",
    caption = paste0("Total terms: ", scales::comma(sum(veddra_plot$n)), " • Width proportional to term count")
  ) +
  
  coord_cartesian(clip = "off") +
  theme_ft()

# Display the plot
veddra_pyramid

# ===============================================================================
# 7. DATA EXPANSION AND ADVERSE REACTION MAPPING
# ===============================================================================

message("🔄 Expanding data and mapping adverse reactions...")

# Expand reactions: one row per dog per reaction
complete_long <- complete %>%
  # Create unique dog identifier
  mutate(dog_id = row_number()) %>%
  # Split multiple reactions into separate rows
  separate_longer_delim(Reaction, delim = ",") %>%
  # Clean reaction text while preserving medical terminology
  mutate(
    llt = str_trim(Reaction),                         # Remove whitespace
    llt = str_to_lower(llt),                         # Standardize case
    # Preserve important medical syntax (/, (), -) but remove clearly problematic chars
    llt = str_remove_all(llt, "[^a-z\\s/()\\-']"),  # Keep medical syntax
    llt = str_squish(llt),                            # Remove extra whitespace
    death = case_when(
      death == 1 ~ "Yes", 
      death == 0 ~ "No",
      TRUE ~ "Unknown")
  ) %>%
  # Remove empty or very short reactions
  filter(!is.na(llt), llt != "", nchar(llt) >= 3) %>%
  # Remove original Reaction column
  select(-Reaction)

message("  Expanded to ", nrow(complete_long), " individual reaction records")

llt_counts <- complete_long %>%
  group_by(drug) %>%
  count(dog_id, name = "n_llts") %>% ungroup()

llt_counts %>%
  group_by(drug) %>%
  summarise(
    max_llts    = max(n_llts, na.rm = TRUE),
    median_llts = median(n_llts, na.rm = TRUE)
  ) %>%
  ungroup()

# Map to VEDDRA terminology
matched <- complete_long %>%
  left_join(veddra, by = "llt", relationship = "many-to-one") %>%
  # Clean and convert to appropriate data types
  mutate(
    # Convert character variables to factors (except complex drug history)
    across(where(is.character) & !any_of("drug_hx"), as.factor),
    # Ensure year is categorical for stratification
    year = factor(year),
    # Death is already cleaned and numeric, convert to factor
    death = factor(death),
    # Recode sex to short form (already cleaned as FEMALE/MALE/UNKNOWN)
    sex = case_when(
      sex == "FEMALE" ~ "f",
      sex == "MALE" ~ "m",
      TRUE ~ "unknown"
    ),
    sex = factor(sex, levels = c("f", "m", "unknown"))
  ) %>%
  # Remove unnecessary columns and standardize names
  select(-any_of("animals_affected")) %>%
  rename_with(tolower)

# Assess mapping success
cat("\n=== REACTION MAPPING ASSESSMENT ===\n")
print(tibble(
  total_reactions = nrow(matched),
  mapped_to_pt = sum(!is.na(matched$pt)),
  mapping_success_rate = round(sum(!is.na(matched$pt)) / nrow(matched) * 100, 1),
  unique_unmapped_llts = n_distinct(matched$llt[is.na(matched$pt)]),
  unique_mapped_pts = n_distinct(matched$pt[!is.na(matched$pt)]))
)

pt_counts <- matched %>%
  count(dog_id, name = "n_pts")

pt_counts %>%
  summarise(
    max_pts    = max(n_pts, na.rm = TRUE),
    median_pts = median(n_pts, na.rm = TRUE)
  )


# Show top unmapped terms for manual review
unmapped_terms <- matched %>%
  filter(is.na(pt)) %>%
  count(llt, sort = TRUE)

if (nrow(unmapped_terms) > 0) {
  cat("\nTop 10 unmapped terms (consider manual mapping):\n")
  print(unmapped_terms)
  
  # Save unmapped terms for review
  write_csv(matched %>% 
              filter(is.na(pt)) %>% 
              count(llt, sort = TRUE),
            here("data", "processed", "unmapped_terms_for_review.csv"))
}

# Save final matched dataset
write_rds(matched, here("data", "processed", "matched_data.rds"))
write_csv(matched, here("data", "processed", "matched_data.csv"))


message("Creating final display visualizations and data summary...")

# -------------------------------------------------------------------------------
# 1. TEMPORAL TRENDS PLOT
# -------------------------------------------------------------------------------

# Calculate PT reports over time by drug
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

# Create temporal trends plot with FT styling and dual metrics
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
  
  # FT color palette
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
  
  # Primary y-axis
  scale_y_continuous(
    name = "Number of Reports",
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.1)),
    # Secondary y-axis
    sec.axis = sec_axis(
      trans = ~ . / scaling_factor,
      name = "Number of PTs reported",
      labels = scales::comma
    )
  ) +
  
  # Financial Times theme
  theme_minimal() +
  theme(
    # Background
    plot.background = element_rect(fill = "#fff1e5", color = NA),
    panel.background = element_rect(fill = "#fff1e5", color = NA),
    
    # Grid
    panel.grid.major = element_line(color = "#f0f0f0", size = 0.5),
    panel.grid.minor = element_blank(),
    
    # Text
    text = element_text(family = "Arial", color = "#2d2d2d"),
    plot.title = element_text(size = 16, face = "bold", color = "#2d2d2d"),
    plot.subtitle = element_text(size = 12, color = "#555555", margin = margin(b = 20)),
    axis.title.y.left = element_text(size = 12, color = "#2d2d2d", margin = margin(r = 10)),
    axis.title.y.right = element_text(size = 12, color = "#555555", margin = margin(l = 10)),
    axis.title.x = element_text(size = 12, color = "#333333"),
    axis.text = element_text(size = 10, color = "#555555"),
    axis.text.y.right = element_text(color = "#888888"),
    
    # Legend
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.background = element_rect(fill = "#fff1e5", color = NA),
    
    # Margins
    plot.margin = margin(t = 20, r = 30, b = 20, l = 20)
  ) +
  
  # Labels
  labs(
    title = "Temporal Trends in the EudraVigalance database",
    subtitle = "Solid lines: ADR reports | Dashed lines: Preferred Terms (PTs) reported",
    x = "Year"
  ) +
  
  # X-axis scale
  scale_x_continuous(breaks = seq(min(temporal_data$year_numeric, na.rm = TRUE), 
                                  max(temporal_data$year_numeric, na.rm = TRUE), by = 2))

# Display the plot
temporal_plot

# Save the plot
ggsave(here("output", "figures", "temporal_trends.png"), temporal_plot, 
       width = 12, height = 7, dpi = 300, bg = "#fff1e5")

message("Temporal trends plot created with dual metrics")

# -------------------------------------------------------------------------------
# 2. DATA SUMMARY TABLE (SIMPLIFIED 3-COLUMN FORMAT)
# -------------------------------------------------------------------------------

# Prepare data summary excluding specified columns
summary_data <- matched %>%
  group_by(dog_id) %>%
  slice(1) %>%
  ungroup() %>%
  select(-any_of(c("species", "drug_hx", "date", "llt", "pt", "hlt", 
                   "animals_affected", "organ", "dog_id", "drug_list"))) %>%
  mutate(
    sex = case_when(
      sex == "m" ~ "Male",
      sex == "f" ~ "Female",
      TRUE ~ "Unknown"
    ),
    across(c(where(is.character), where(is.factor)), 
                ~str_to_title(as.character(.x))))

# Create 3-column summary table with tidyverse piping
data_summary_table_data <- summary_data %>%

  # Convert to long format with column names
  pivot_longer(everything(), names_to = "variable_name", values_to = "value") %>%
  # Group by variable to calculate distributions
  group_by(variable_name) %>%
  nest() %>%
  mutate(
    # Create summary rows for each variable
    summary = map2(variable_name, data, ~{
      # Get frequency counts
      freq <- .y %>%
        count(value, name = "n") %>%
        arrange(desc(n)) %>%
        mutate(
          percentage = n / sum(n) * 100,
          variable = paste0("  ", value)  # Indent level names
        )
      
      # Check if variable has >8 levels, if so only keep top 5
      n_levels <- nrow(freq)
      if (n_levels > 8) {
        freq <- freq %>%
          slice_head(n = 5) %>%
          bind_rows(tibble(
            value = NA,
            n = sum(freq$n[6:n_levels]),
            percentage = sum(freq$percentage[6:n_levels]),
            variable = paste0("  Other (", n_levels - 5, " levels)")
          ))
      }
      
      freq <- freq %>%
        select(variable, n, percentage)
      
      # Add header row with variable name
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

# Create publication-quality summary table with FT styling
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
  fmt_number(
    columns = n,
    decimals = 0,
    use_seps = TRUE
  ) %>%
  fmt_number(
    columns = percentage,
    decimals = 1,
    pattern = "{x}%"
  ) %>%
  fmt_missing(
    columns = everything(),
    missing_text = ""
  ) %>%
  # Bold the variable headers (rows where n is NA)
  tab_style(
    style = cell_text(weight = "bold", color = "#1a1a1a"),
    locations = cells_body(
      columns = variable,
      rows = is.na(n)
    )
  ) %>%
  # Style the level rows (indented)
  tab_style(
    style = cell_text(color = "#333333"),
    locations = cells_body(
      columns = variable,
      rows = !is.na(n)
    )
  ) %>%
  # Financial Times styling
  tab_style(
    style = list(
      cell_fill(color = "#fff1e5"),
      cell_text(color = "#2d2d2d", size = "small")
    ),
    locations = cells_body()
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "#f2e6dd"),
      cell_text(weight = "bold", color = "#2d2d2d")
    ),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "#e8ddd4"),
      cell_text(weight = "bold", color = "#1a1a1a")
    ),
    locations = cells_title()
  ) %>%
  tab_style(
    style = cell_text(color = "#555555", size = "small"),
    locations = cells_source_notes()
  ) %>%
  cols_align(
    align = "left",
    columns = variable
  ) %>%
  cols_align(
    align = "right",
    columns = c(n, percentage)
  ) %>%
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
    n ~ px(150),
    percentage ~ px(150)
  )

# Display the summary table
data_summary_table

# Save the summary table
gtsave(data_summary_table, here("output", "tables", "data_summary.html"))

message("Data summary table created")
# ===============================================================================
# 8. FINAL CLEANING SUMMARY AND VALIDATION
# ===============================================================================

message("📋 Generating final cleaning summary...")

# Generate comprehensive cleaning summary
final_summary <- list(
  cleaning_date = Sys.Date(),
  input_files = length(file_list),
  drugs_processed = length(target_drugs),
  initial_records = nrow(all_drugs_data),
  final_records = nrow(matched),
  total_reactions = nrow(matched),
  mapped_reactions = sum(!is.na(matched$pt)),
  mapping_success_rate = round(sum(!is.na(matched$pt)) / nrow(matched) * 100, 1),
  unique_dogs = n_distinct(matched$dog_id),
  unique_pts = n_distinct(matched$pt, na.rm = TRUE),
  data_quality_score = round(mean(quality_report$missing$missing_percent <= 5) * 100, 0),
  polypharmacy_rate = round(mean(matched$polypharmacy == "Polypharmacy") * 100, 1)
)

# Save summary as JSON
write_json(final_summary, here("data", "processed", "cleaning_summary.json"), pretty = TRUE)

# Display final summary
cat("\n=== FINAL CLEANING SUMMARY ===\n")
cat("Input files processed:", final_summary$input_files, "\n")
cat("Initial raw records:", format(final_summary$initial_records, big.mark = ","), "\n")
cat("Final cleaned records:", format(final_summary$final_records, big.mark = ","), "\n")
cat("Unique dogs:", format(final_summary$unique_dogs, big.mark = ","), "\n")
cat("VEDDRA mapping success:", final_summary$mapping_success_rate, "%\n")
cat("Polypharmacy rate:", final_summary$polypharmacy_rate, "%\n")

# Log completion
final_log <- paste0(
  "\n=== CLEANING COMPLETION ===\n",
  "Completion time: ", Sys.time(), "\n",
  "Total runtime: ", round(difftime(Sys.time(), file.info(log_file)$mtime, units = "mins"), 1), " minutes\n",
  "Status: SUCCESS\n",
  "Output: matched_data.rds/csv in data/processed/\n"
)

cat(final_log, file = log_file, append = TRUE)
cat(final_log)

message("✅ Data cleaning pipeline completed successfully!")
message("📁 Clean data saved as: data/processed/matched_data.rds")
message("📊 Ready for statistical signal detection analysis")

# ===============================================================================
# END OF CLEANING SCRIPT
# ===============================================================================