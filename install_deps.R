#------------------------------------------------------------------------------#
# INSTALL DEPENDENCIES FOR RUNNING DES OF
# SICK SICKER MODEL WITH RECURRENCE
# Author: Mauricio Lopez-Mendez PhD (Stanford Health Policy)
#------------------------------------------------------------------------------#
# Description ----
#* 1. Ensure package renv is installed, and restores project library version from
#* renv.lock.
#* 2. Checks for currently installed dependencies and installs any require dependency 
#* that may be missing.
#------------------------------------------------------------------------------#


# 0. Check R version 3.5 or newer 
if (getRversion() >= "3.5.0") {
  message("R version is 3.5.0 or newer.")
} else {
  stop("R version is older than 3.5.0.")
}

# 1. Ensure renv is available
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# 2.  Required packages
required_pkgs <- c(
  "data.table"  ,   # to manipulate data
  "dplyr"       ,   # to manipulate data
  "tidyr"       ,   # to manipulate data
  "reshape2"    ,   # to manipulate data
  
  "ggplot2"     ,   # to visualize data
  "ggrepel"     ,   # to visualize data
  "gridExtra"   ,   # to visualize data
  "ellipse"     ,   # to visualize data
  "ggview"      ,   # save plots
  
  "scales"      ,   # for dollar signs and commas
  "patchwork"   ,   # for combining ggplot2 figures
  "dampack"     ,   # for CEA and calculate ICERs
  
  "doParallel"  ,   # parallel processing
  "parallel"    ,   # parallel processing
  "foreach"     ,   # parallel processing
  
  "stats"       ,   # essential statistical functions
  "MethylCapSig",   # has nice multivariate lognormal random variable generator
  "survival"    ,   # core survival analysis routines
  "flexsurv"    ,   # flexible parametric survival models and multistate models
  
  "devtools"    ,   # to install packages from GitHub
  "abind"       ,   # Combine multi-dimensional arrays
  "matrixStats"     # functions operating on rows and columns of matrices
  
)

installed_pkgs   <- rownames(installed.packages())
to_install_pkgs  <- setdiff(required_pkgs, installed_pkgs)

# Install missing packages
if (length(to_install_pkgs)) {
  message("Installing extra packages: ", paste(to_install, collapse = ", "))
  failed_pkgs <- character()
  
  for (pkg in to_install_pkgs) {
    message("→ Installing ", pkg, " …")
    tryCatch(
      install.packages(pkg, dependencies = TRUE),
      error = function(e) {
        warning(sprintf("  ✗ Failed to install '%s': %s", pkg, e$message))
        failed_pkgs <<- c(failed_pkgs, pkg)
      }
    )
  }
  
  if (length(failed_pkgs)) {
    warning(
      "The following packages could not be installed:\n",
      paste0(" - ", failed_pkgs, collapse = "\n")
    )
  }
}


# 3. Restore the project library from renv.lock
message("Restoring project library via renv...")
renv::restore(prompt = T)


message("Setup complete!")



