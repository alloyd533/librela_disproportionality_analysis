library(tidyverse)

# Define drug names
drugs <- c("librela", "metacam", "previcox", "galliprant", "rimadyl", "onisor", "daxocox")

# List all CSV files in /data

file_list <- list.files("data", pattern = "^\\d{4}_eudra_.*\\.csv$", full.names = TRUE)

read_and_clean <- function(path) {
  df <- read_csv(path, col_types = cols(.default = "c"), show_col_types = FALSE)
  
  # Clean column names
  names(df) <- names(df) |> str_replace_all(" ", "_")
  
  # Handle date parsing
  if ("Received_date" %in% names(df)) {
    raw_date <- df$Received_date
    
    # Try parse_date_time
    parsed <- parse_date_time(raw_date, orders = c("ymd", "dmy", "ymd HMS", "dmy HMS"), quiet = TRUE)
    
    # If all NA, try as.Date fallback
    if (all(is.na(parsed))) {
      parsed <- suppressWarnings(as.Date(raw_date, format = "%Y-%m-%d"))
      if (all(is.na(parsed))) {
        parsed <- suppressWarnings(as.Date(raw_date, format = "%d/%m/%Y"))
      }
    }
    
    df <- df |> mutate(date = parsed) |> select(-Received_date)
  } else {
    df <- df |> mutate(date = NA_Date_)
  }
  
  df <- df |> mutate(
    source_file = basename(path),
    year = str_extract(source_file, "^\\d{4}")
  )
  
  return(df)
}

# Loop through drugs and process
for (drug in drugs) {
  drug_files <- keep(file_list, ~ str_detect(.x, paste0("eudra_", drug, "\\.csv$")))
  
  if (length(drug_files) > 0) {
    combined <- map_dfr(drug_files, read_and_clean)
    
    # Assign to environment with drug name
    assign(drug, combined)
    
    # Save to file in /data/
    write_csv(combined, file = file.path("data", paste0(drug, "_combined.csv")))
  }
}

### Assign the veddra 
veddra <- readxl::read_excel("data/combined-veddra-list-clinical-terms-reporting-suspected-adverse-events-animals-humans-veterinary-medicinal-products_en.xlsx") %>%
  rename_with(~ str_replace_all(.x, " ", "_"))
