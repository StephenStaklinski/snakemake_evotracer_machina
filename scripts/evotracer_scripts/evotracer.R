# if (!requireNamespace("EvoTraceR", quietly = TRUE)) {
#   if (!requireNamespace("devtools", quietly = TRUE)) {
#     install.packages("devtools", repos = "https://cloud.r-project.org/")
#   }
#   devtools::install_github("Nowak-Lab/EvoTraceR@v1.0.1")
# }

library(EvoTraceR)
args = commandArgs(trailingOnly=TRUE)

#input_dir <- system.file("extdata", "input", package = "EvoTraceR")
input_dir <- args[1]
output_dir <- args[2]
cutoff <- as.numeric(args[3])
ref_name <- args[4]
ref_seq <- args[5]
ref_flank_left <- args[6]
ref_flank_right <- args[7]
ref_cut_sites <- as.numeric(unlist(strsplit(args[8], ",")))
ref_border_sites <- as.numeric(unlist(strsplit(args[9], ",")))

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

trimmomatic_path <- Sys.getenv("TRIMMOMATIC_PATH")
flash_path <- Sys.which("flash")

EvoTraceR_object <-
  initialize_EvoTraceR(
    input_dir = input_dir,
    output_dir = output_dir,
    trimmomatic_path = trimmomatic_path,
    flash_path = flash_path)


EvoTraceR_object <-
  asv_analysis(EvoTraceR_object = EvoTraceR_object,
               ref_name = ref_name,
               ref_seq = ref_seq,
                 ref_flank_left = paste0("^", ref_flank_left),
               ref_flank_right = paste0(ref_flank_right, "$"),
               ref_cut_sites = ref_cut_sites,
               ref_border_sites = ref_border_sites,
               output_figures = TRUE,
               asv_count_cutoff = cutoff, # minimum number of ASVs to be counted; decided on: 03/25/22
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

save.image(paste0(output_dir, "/", "evotracer.RData"))
