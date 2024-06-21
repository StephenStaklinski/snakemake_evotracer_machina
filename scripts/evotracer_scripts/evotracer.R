if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", dependencies=TRUE)
  BiocManager::install('ggtree')
} else {
  if (!requireNamespace("ggtree", quietly = TRUE)) {
    BiocManager::install('ggtree')
  }
  BiocManager::install('ggtree')
}

if (!requireNamespace("EvoTraceR", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  devtools::install_github("Nowak-Lab/EvoTraceR@v1.0.1")
}

library(EvoTraceR)
args = commandArgs(trailingOnly=TRUE)

#input_dir <- system.file("extdata", "input", package = "EvoTraceR")
input_dir <- args[1]
output_dir <- args[2]

# unzip files
# List all files in the directory
files <- list.files(input_dir)
# Filter files with ".zip" extension
zip_files <- files[grepl("\\.zip$", files)]
# Check if there are any zip files
if (length(zip_files) > 0) {
  # Loop through each zip file and unzip its contents
  for (zip_file in zip_files) {
    # Specify the full path to the zip file
    zip_file_path <- file.path(input_dir, zip_file)
    # Unzip the file
    unzip(zip_file_path, exdir = input_dir)
    cat("Unzipped:", zip_file, "\n")
  }
} else {
  cat("No zip files found in the directory.\n")
}


trimmomatic_path <- "scripts/Trimmomatic-0.39/trimmomatic-0.39.jar"
flash_path <- Sys.which("flash")

EvoTraceR_object <-
  initialize_EvoTraceR(
    input_dir = input_dir,
    output_dir = output_dir,
    trimmomatic_path = trimmomatic_path,
    flash_path = flash_path)


EvoTraceR_object <-
  asv_analysis(EvoTraceR_object = EvoTraceR_object,
               ref_name = "BC10v0",
               ref_seq = "TCTACACGCGCGTTCAACCGAGGAAAACTACACACACGTTCAACCACGGTTTTTTACACACGCATTCAACCACGGACTGCTACACACGCACTCAACCGTGGATATTTACATACTCGTTCAACCGTGGATTGTTACACCCGCGTTCAACCAGGGTCAGATACACCCACGTTCAACCGTGGTACTATACTCGGGCATTCAACCGCGGCTTTCTGCACACGCCTACAACCGCGGAACTATACACGTGCATTCACCCGTGGATC",
               ref_flank_left = "^TCTAC",
               ref_flank_right = "CCCGTGGATC$",
               ref_cut_sites = c(17, 43, 69, 95, 121, 147, 173, 199, 225, 251),
               ref_border_sites = c(1, 26, 52, 78, 104, 130, 156, 182, 208, 234),
               output_figures = TRUE,
               asv_count_cutoff = 3, # minimum number of ASVs to be counted; decided on: 03/25/22
               # pair-wise alignment parameters between un-edited barcode and edited barcode (ASV)
               pwa_type = "global", # based on AmpliCan (global = Needleman-Wunsch)
               pwa_gapOpening = -25, # based on AmpliCan: -25
               pwa_gapExtension = 0, # based on AmpliCan: 0
               pwa_match = 15, # based on AmpliCan: 15
               pwa_mismatch = -4, # based on AmpliCan: -4
               cleaning_window = c(13, 13), # cleaning window +/- from Cas9 editing size (nucleotide 17 in guide) is considered as an edit 
               batch_size = 100,
               cores = parallel::detectCores()               
               )

EvoTraceR_object <-
  analyse_mutations(EvoTraceR_object = EvoTraceR_object)

EvoTraceR_object <-
  infer_phylogeny(EvoTraceR_object = EvoTraceR_object, mutations_use = "del_ins")

EvoTraceR_object <-
  create_df_summary(EvoTraceR_object)

save.image(paste0(output_dir, "/", basename(output_dir), "_", "EvoTraceR.RData"))
