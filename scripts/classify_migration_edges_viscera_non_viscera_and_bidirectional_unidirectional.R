#!/usr/bin/env Rscript

# Required libraries
library(tidyverse)
library(dplyr)

# Function to classify edges based on tissue types
classify_edges <- function(data, visceral_tissues, non_visceral_tissues) {
  # Filter for edges with probability >= 0.5
  filtered_data <- data %>%
    filter(probability_across_mach2_solutions >= 0.5)
  
  # Classify each edge
  classified_data <- filtered_data %>%
    mutate(
      edge_type = case_when(
        source == "PRL" & target %in% visceral_tissues ~ "primary_to_visceral",
        source == "PRL" & target %in% non_visceral_tissues ~ "primary_to_non_visceral",
        source %in% visceral_tissues & target %in% non_visceral_tissues ~ "visceral_to_non_visceral",
        source %in% non_visceral_tissues & target %in% visceral_tissues ~ "non_visceral_to_visceral",
        source %in% visceral_tissues & target == "PRL" ~ "visceral_to_primary",
        source %in% non_visceral_tissues & target == "PRL" ~ "non_visceral_to_primary",
        source %in% visceral_tissues & target %in% visceral_tissues ~ "visceral_to_visceral",
        source %in% non_visceral_tissues & target %in% non_visceral_tissues ~ "non_visceral_to_non_visceral",
        TRUE ~ "other"
      )
    )
  
  # Print details of "other" edges to console
  other_edges <- classified_data %>%
    filter(edge_type == "other")
  
  if (nrow(other_edges) > 0) {
    cat("\nEdges classified as 'other':\n")
    print.data.frame(other_edges %>% select(mouse, cp, source, target, probability_across_mach2_solutions), row.names = FALSE)
  }
  
  # Count edges per category for each mouse and cp
  summary_data <- classified_data %>%
    group_by(mouse, cp, edge_type) %>%
    summarise(count = n(), .groups = "drop") %>%
    pivot_wider(
      names_from = edge_type,
      values_from = count,
      values_fill = 0
    )
  
  return(summary_data)
}

# Function to analyze directionality of edges
analyze_edge_directionality <- function(data) {
  # Filter for edges with probability >= 0.5
  filtered_data <- data %>%
    filter(probability_across_mach2_solutions >= 0.5)
  
  # First identify which tissue pairs are bidirectional
  bidirectional_pairs <- filtered_data %>%
    # Create a canonical tissue pair representation
    mutate(
      tissue_pair = pmap_chr(list(source, target), ~ paste(sort(c(..1, ..2)), collapse = "_"))
    ) %>%
    # Group by mouse, cp, and tissue pair to identify bidirectional pairs
    group_by(mouse, cp, tissue_pair) %>%
    summarise(
      direction_count = n_distinct(paste(source, target)),
      is_bidirectional = direction_count > 1,
      .groups = "drop"
    )
  
  # Join this information back to filtered data and count actual edges
  edge_counts <- filtered_data %>%
    # Create the same tissue pair column for joining
    mutate(
      tissue_pair = pmap_chr(list(source, target), ~ paste(sort(c(..1, ..2)), collapse = "_"))
    ) %>%
    # Join with bidirectional classification
    left_join(bidirectional_pairs %>% select(mouse, cp, tissue_pair, is_bidirectional), 
              by = c("mouse", "cp", "tissue_pair")) %>%
    # Group by mouse and cp
    group_by(mouse, cp) %>%
    # Count actual edges in each category
    summarise(
      unidirectional_edge_count = sum(!is_bidirectional),
      bidirectional_edge_count = sum(is_bidirectional),
      .groups = "drop"
    )
  
  return(edge_counts)
}

# Main function
main <- function() {

  input_file <- "/grid/siepel/home/staklins/snakemake_evotracer_machina/results/analysis_on_3_7_25_serio_pca_2p_2pr_data_3_7_25/mach2/2p_2pr_all_mach2_data.csv"
  
  # Define tissue categories
  visceral_tissues <- c("LVL", "LVL1", "LVL2", "LVM", "LVR", "LGL", "LGT", "LGR", "SMI", "DIP", "ADR", "ADL", "SPN", "SMJ", "SMD", "KDR", "KDL")
  non_visceral_tissues <- c("FMR", "FML", "TBR", "TBL", "HMR", "HML", "UNR", "UNL", "UNR1", "UNR2", "RBR", "RBL", "SPL", "SPS", "SPT", "SPC", "LNIL", "LNIR", "LNLL", "LNLR", "LNRL", "LNRR", "LNSR", "LNSL", "LNMR", "LNML", "LNAL", "LNAR", "LNBL", "LNBR", "LNCL", "LNCR", "LNYL", "LNYR", "LNDL", "LNDR")
  
  # Read the data
  data <- read_csv(input_file)
  
  # Process tissue type classifications
  result_tissue_types <- classify_edges(data, visceral_tissues, non_visceral_tissues)
  
  # Process edge directionality
  result_directionality <- analyze_edge_directionality(data)
  
  # Write the outputs
  output_file_tissue_types <- sub("\\.csv$", "_classified_viscera_non_viscera.csv", input_file)
  output_file_directionality <- sub("\\.csv$", "_edge_directionality.csv", input_file)
  
  write_csv(result_tissue_types, output_file_tissue_types)
  write_csv(result_directionality, output_file_directionality)
  
  cat("Processing complete. Results written to:\n")
  cat("1.", output_file_tissue_types, "\n")
  cat("2.", output_file_directionality, "\n")
}

# Run the main function
main() 