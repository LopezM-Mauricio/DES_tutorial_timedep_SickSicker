# DES_Tutorial_beta
Beta Release of DES Tutorial. 
V.1.04.1.2025

Please cite the preprint article: https://doi.org/10.1101/2025.05.15.25327635

Structure of Repository: 

DES_Tutorial_beta/

├── DES_Tutorial.Rproj

├── install_deps.R

├── renv/ (ignore)

├── renv.lock

├── .Rhistory

├── .Rprofile

├── README.R

├── inputs/

    └── all_cause_mortality.rda
    
    └── LifeTable_USA_Mx_2015.csv    
    
├── sim_modules/
   
    └── DES_Sick_Sicker_progressive.R
    
    └── DES_Sick_Sicker_recurrence.R
    
├── analysis modules/
    
    └── Module_A_CEA.R
    
    └── Module_B_Epi_Outcomes.R
    └── Module_C_PSA.R
    
    └── Aux_Module_Cont_Time_Trace.R
    
    └── Aux_Module_Convergence.R
    
├── manuscript/
   
    └── manuscript.qmd (.docx; .pdf)
    
    └── appendix.qmd   (.docx; .pdf)
    
    └── bibliography.bib
    
    └── apa-single-spaced.csl
    
    └── figures/
    
└── R/
    
    └── DES_functions.R
    
    └── Functions.R
    
    └── Functions_cSTM_time_dep_simulation.R
    
    └── Functions_cSTM_time_dep_state_residence.R


# Workflow for collaborators

  - Clone the repo.

  - Open the project, run in RStudio Terminal:

`Rscript install_deps.R`

This will install `renv`, restore locked versions of all dependencies, and install any missing packages. Script will return a list of packages that were not successfully installed and require individual troubleshooting.

->  If all goes well, your environment should now match exactly the environment of the developer :P


# Troubleshooting Package installation

If you run into incompatibilities or mysterious errors running `install.deps.R`, verify that R version 3.5.0 or newer is installed.

You can try to install required packages manually, run the below code on a script or from the R console:

```
# Install `renv`

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Required packages

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

# Restore the project library from renv.lock
message("Restoring project library via renv...")
renv::restore(prompt = T)

```
Contact me at 

mlopezme@stanford.edu; or  mlopezme@gmail.com
