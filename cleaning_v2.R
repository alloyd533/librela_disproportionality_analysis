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

# Advanced polypharmacy detection function
detect_polypharmacy_advanced <- function(drug_string) {
  
  if (is.na(drug_string) || str_trim(drug_string) == "") return("0")
  
  # Handle multiple possible delimiters and conjunctions
  delimiters <- c(",", ";", "\\+", "&", " and ", " AND ", " with ", " WITH ")
  
  # Split by any delimiter
  drug_list <- drug_string
  for (delim in delimiters) {
    drug_list <- unlist(strsplit(drug_list, delim, perl = TRUE))
  }
  
  # Clean and extract meaningful drug names
  drug_names <- drug_list %>%
    str_trim() %>%                                    # Remove whitespace
    str_to_lower() %>%                                # Standardize case
    str_remove_all("\\(.*\\)") %>%                   # Remove parenthetical info
    str_remove_all("\\d+\\s*(mg|ml|g|%).*$") %>%     # Remove dosing information
    str_remove_all("\\b(tablet|capsule|injection|oral)\\b") %>%  # Remove formulation terms
    str_extract("^[a-z][a-z0-9]*") %>%               # Extract drug name (letters + numbers)
    .[!is.na(.) & . != "" & nchar(.) >= 3]          # Remove short/empty entries
  
  # Count unique meaningful drug names
  unique_drugs <- unique(drug_names)
  
  # Return polypharmacy indicator
  return(if (length(unique_drugs) > 1) "1" else "0")
}

# Apply polypharmacy detection with progress tracking
message("  Analyzing ", nrow(complete), " drug history records...")

complete <- complete %>%
  mutate(
    polypharmacy = map_chr(drug_hx, detect_polypharmacy_advanced),
    polypharmacy = factor(polypharmacy, 
                          levels = c("0", "1"), 
                          labels = c("Monotherapy", "Polypharmacy"))
  )

# Generate polypharmacy summary
poly_summary <- complete %>%
  count(polypharmacy) %>%
  mutate(percentage = round(n / sum(n) * 100, 1))

cat("\n=== POLYPHARMACY ASSESSMENT ===\n")
print(poly_summary)

# Cross-tabulation with drug type
poly_by_drug <- complete %>%
  count(drug_hx, polypharmacy) %>%
  group_by(drug_hx) %>%
  mutate(drug_total = sum(n),
         percentage = round(n / drug_total * 100, 1)) %>%
  filter(polypharmacy == "Polypharmacy") %>%
  arrange(desc(percentage))

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

# Alternative version with gradient colors for more visual appeal
veddra_pyramid_gradient <- ggplot(df) +
  geom_rect(aes(xmin = x_min, xmax = x_max, ymin = y_bot, ymax = y_top, fill = level),
            colour = "#ffffff", linewidth = 1.5, alpha = 0.85) +
  
  # Add subtle shadow effect
  geom_rect(aes(xmin = x_min + 0.005, xmax = x_max + 0.005, 
                ymin = y_bot + 0.05, ymax = y_top + 0.05),
            fill = "#00000015", colour = NA) +
  
  geom_rect(aes(xmin = x_min, xmax = x_max, ymin = y_bot, ymax = y_top, fill = level),
            colour = "#ffffff", linewidth = 1.5, alpha = 0.85) +
  
  geom_text(aes(x = x_mid, y = y_mid, label = level), 
            fontface = "bold", 
            colour = "white",
            size = 4.5) +
  
  geom_text(aes(x = x_mid, y = y_mid + 0.15, label = scales::comma(n)), 
            fontface = "normal", 
            colour = "white",
            size = 3.5) +
  
  scale_fill_manual(values = c("#d62728", "#ff7f0e", "#2ca02c", "#1f77b4")) +
  scale_y_reverse(expand = c(0, 0.1), limits = c(nrow(df) + 0.1, -0.1)) +
  scale_x_continuous(expand = c(0, 0), limits = c(-x_lim, x_lim)) +
  
  labs(
    title = "VeDDRA Terminology Hierarchy",
    subtitle = "Veterinary adverse reaction terms by classification level • Width ∝ number of terms",
    caption = paste0("Total: ", scales::comma(sum(veddra_plot$n)), " terms")
  ) +
  
  coord_cartesian(clip = "off") +
  theme_ft()

# Display the gradient version
veddra_pyramid_gradient
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
    llt = str_squish(llt)                            # Remove extra whitespace
  ) %>%
  # Remove empty or very short reactions
  filter(!is.na(llt), llt != "", nchar(llt) >= 3) %>%
  # Remove original Reaction column
  select(-Reaction)

message("  Expanded to ", nrow(complete_long), " individual reaction records")

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
mapping_stats <- tibble(
  total_reactions = nrow(matched),
  mapped_to_pt = sum(!is.na(matched$pt)),
  mapping_success_rate = round(sum(!is.na(matched$pt)) / nrow(matched) * 100, 1),
  unique_unmapped_llts = n_distinct(matched$llt[is.na(matched$pt)]),
  unique_mapped_pts = n_distinct(matched$pt[!is.na(matched$pt)])
)

cat("\n=== REACTION MAPPING ASSESSMENT ===\n")
print(mapping_stats)

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