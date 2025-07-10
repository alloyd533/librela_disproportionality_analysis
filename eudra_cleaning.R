library(tidyverse)

# Define drug names
drugs <- c("librela", "metacam", "previcox", "galliprant", "rimadyl", "onisor", "daxocox")

# List all CSV files in /data

file_list <- list.files("data", pattern = "^\\d{4}_eudra_.*\\.csv$", full.names = TRUE)

# Loop through drugs and process
for (drug in drugs) {
  # Find all matching files for this drug
  drug_files <- file_list |> keep(~ str_detect(.x, paste0("eudra_", drug, "\\.csv$")))
  
  # Bind all years into one dataframe
  if (length(drug_files) > 0) {
    df <- map_dfr(drug_files, ~ read_csv(.x, show_col_types = FALSE) |>
                    mutate(
                      source_file = basename(.x),
                      year = str_extract(source_file, "^\\d{4}")
                    ))
    
    # Save combined dataframe
    write_csv(df, paste0(drug, "_combined.csv"))
    
    # Remove from environment
    rm(df)
  }
}
