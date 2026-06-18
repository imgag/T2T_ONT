#!/usr/bin/env Rscript

# Script to install SVbyEye from GitHub
message("Checking if SVbyEye is already installed...")

# Check if SVbyEye is installed
if (requireNamespace("SVbyEye", quietly = TRUE)) {
    message("SVbyEye is already installed.")
    quit(save = "no", status = 0)
}

message("SVbyEye not found. Installing from GitHub...")

# Load devtools
suppressPackageStartupMessages({
    if (!requireNamespace("devtools", quietly = TRUE)) {
        install.packages("devtools", repos = "https://cloud.r-project.org")
    }
    library(devtools)
})

# Install SVbyEye from GitHub
tryCatch({
    devtools::install_github("daewoooo/SVbyEye", branch = "master", upgrade = "never")
    message("SVbyEye installed successfully!")
}, error = function(e) {
    message("Error installing SVbyEye: ", e$message)
    quit(save = "no", status = 1)
})

# Verify installation
if (requireNamespace("SVbyEye", quietly = TRUE)) {
    message("Installation verified successfully.")
} else {
    message("Installation verification failed.")
    quit(save = "no", status = 1)
}