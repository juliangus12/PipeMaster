# Set the CRAN mirror
options(repos = c(CRAN = "https://archive.linux.duke.edu/cran/"))

# Function to install packages if they are not already installed
install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    install.packages(package)
  }
}

# List of packages required by the scripts, now including 'bda', 'Hmisc', 'crayon', and 'progress'
packages <- c("ggplot2", "plyr", "tidyr", "readr", "lubridate", "zoo", "xtable", "Rcpp", "bda", "Hmisc", "plotrix", "gplots", "evir", "dplyr", "ncdf4", "languageserver", "crayon", "progress")

# Install missing packages
invisible(sapply(packages, install_if_missing))

# Load the crayon and progress packages
library(crayon)
library(progress)

# Set the working directory to the location of the 'PipeDesign-master' directory
setwd("~/Desktop/CSO_Design_Optimization/PipeDesign-master")

# List of scripts to run in order, according to the repository structure
scripts_to_run <- c(
  "projections/CCSM4.R",
  "projections/CanESM2.R",
  "projections/GFDLESM2M.R",
  "projections/GFDLesm2m_WRF.R",
  "projections/INMCM4.R",
  "projections/MIROC5.R",
  "projections/MIROCESM.R",
  "projections/MPIESMLR_RegCM4.R",
  "projections/MPIESMLR_WRF.R",
  "projections/NAOIndex.R",
  "projections/globaltemp.R",
  "projections/localtemp.R",
  "projections/mdrtemp.R",
  "projections/tropicalstorm.R",
  "Figure/design_pipe.R",
  "Figure/precip_estimates.R",
  "Figure/precip_projections.R"
)

# Function to source scripts in order
run_scripts_in_order <- function(scripts) {
  pb <- progress_bar$new(total = length(scripts), clear = FALSE, width = 50)
  
  for (i in seq_along(scripts)) {
    script_path <- file.path(getwd(), scripts[i]) # Construct the full path using the current working directory
    
    if(file.exists(script_path)) {
      pb$tick()
      cat(green("\nRunning: "), script_path, "\n")
      
      tryCatch({
        source(script_path)
        cat(green("Success: "), script_path, "\n")
      }, error = function(e) {
        cat(red("Error in script "), script_path, ": ", e$message, "\n")
        traceback(2) # Print the stack trace of the error
      })
      
      # Add a page break after every few scripts
      
    } else {
      cat(yellow("Script not found: "), script_path, "\n")
      cat("Check the file path and ensure it is relative to the working directory: ", getwd(), "\n")
      cat("Full path attempting to read from: ", script_path, "\n")
    }
  }
}

# Run the scripts in order
run_scripts_in_order(scripts_to_run)
