# Install required R packages for the stomatal correction scripts

# List of required packages
required_packages <- c("nleqslv", "ggplot2", "gridExtra")

# Check which packages are not installed
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

# Install missing packages
if(length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages, repos = "https://cloud.r-project.org/")
} else {
  cat("All required packages are already installed.\n")
}

# Load packages to verify installation
for(pkg in required_packages) {
  if(require(pkg, character.only = TRUE)) {
    cat(paste("✓", pkg, "loaded successfully\n"))
  } else {
    cat(paste("✗", pkg, "failed to load\n"))
  }
}