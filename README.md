# DES_Tutorial_timedep_SickSicker

## Authors:

<<<<<<< Updated upstream
<<<<<<< Updated upstream
## Citation(s):

# Quick Guide:

Follow these steps to run the simulation of the Sick-Sicker Model in continuous time and/or reproduce the CEA results and generate the exhibits from the tutorial paper.

## Set up:

1.  Clone the repository
2.  Verify R version 3.5.0 or newer is installed in your machine

## A. Lightweight Modules:

Run DES simulation, compute epidemiological and economic outcomes of the example of the tutorial.

1.  Verify that `dplyr` and `data.table` are installed.
2.  Run either of the two sim_modules:

-   `Lite_modules/DES_Sick_Sicker_progressive.R`: runs the simulation following the progressive version fo the Sick Sicker model
-   `Lite_modules/DES_Sick_Sicker_recurrence.R` : runs the simulation following the version fo the Sick Sicker model that includes recurrence to the Healthy state

3.  Compute the Epidemiological and/or Economic Outcomes based on the DES of the Sick Sicker model with recurrence to Healthy.

-   `Lite_modules/Module_A_CEA.R`: computes a cost-effectiveness analysis comparing the strategies described in the coded example of the tutorial
-   `Lite_modules/Module_B_Epi_Outcomes.R`: computes key epidemiological outcomes comparing the strategies described in the coded example of the tutorial

## B. Parallelized Modules:

We strongly suggest to use parallel engines to compute the probabilistic analysis and to do the analysis of the appropriate simulation sample size.

1.  To run distributed runs of the DES you need to have installed the standard packages for parallelization:

-   `doParallel`
-   `parallel`
-   `abind`
-   `foreach`

and a specialized package to sample from the multivariate log-normal distribution:

-   `MethylCapSig`

2.  Run distributed probabilistic analysis

-   `Parallelized_Modules/Module_C_Probabilistic_Analysis.R`: computes a probabilistic analysis including

1)  Incremental cost-effectiveness ratios (ICERs) with probabilistic output,
2)  Cost-effectiveness acceptability curves (CEACs) and frontier (CEAF),
3)  Expected Loss Curves (ELCs), and
4)  Expected value of perfect information (EVPI).

<!-- -->

3.  Run analysis of the appropriate simulation sample size

-   `Parallelized_Modules/Aux_Module_Convergence.R`: computes bootstrap distribution of prevalence estimates at each state for varying simulation samples sizes and a fixed number of replications.

## C. Code to Reproduce Plots

1.  If you are interested in reproducing the plots presented in the tutorial verify that you have installed additional packages:

-   `ggplot2`
-   `ggrepel`
-   `viridis`
-   `patchwork`
-   `ellipse`

2.  To generate the plots shown in the Tutorial paper use the following scripts:

-   `Plotting_Modules/Module_CEA_wPlots.R`
-   `Plotting_Modules/Module_Epi_Outcomes_wPlots.R`
-   `Plotting_Modules/Aux_Module_Cont_Time_Trace_wPlots.R`
-   `Plotting_Modules/Module_Probabilistic_Analysis_wPlots.R`
-   `Plotting_Modules/Aux_Module_Convergence_wPlots.R`

# Structure of Repository:
=======
=======
>>>>>>> Stashed changes
## Quick Guide for Reproducing the Coded-Example of the Tutorial: 

Follow these simple steps to reproduce part or all results and exhibits from the Tutorial Paper.

### Set up

1. Clone the repository
2. Verify R version 3.5.0 or newer is installed in your machine

### A. Run DES simulation, compute epidemiological and economic outcomes

1. Verify that `dplyr` and `data.table` are installed. 
2. Run either of the two sim_modules: 
  - `sim_modules/DES_Sick_Sicker_progressive.R`: runs the simulation following the progressive version fo the Sick Sicker model
  - `sim_modules/DES_Sick_Sicker_recurrence.R` : runs the simulation following the  version fo the Sick Sicker model that includes recurrence to the Healthy state
3. Compute the Epidemiological and/or Economic Outcomes based on the DES 
  - `analysis_modules/Module_A_CEA.R`: computes a cost-effectiveness analysis comparing the strategies described in the coded example of the tutorial
  - `analysis_modules/Module_B_Epi_Outcomes.R`: computes key epidemiological outcomes comparing the strategies described in the coded example of the tutorial

### B. Parallelized Code to run Probabilistic Analysis and Determine Simulation Sample Size 

We strongly suggest to use parallel runs to compute the probabilistic analysis and to do the analysis of the appropriate simulation sample size

1. To run distributed runs of the DES you need to verify that you have installed the standard  packages for parallelization:

  - `doParallel`
  - `parallel`  
  - `foreach`   
  
  and a specialized package to sample from the multivariate log-normal distribution: 
  
  - `MethylCapSig`
  
2. Run probabilistic analysis 
  - `analysis_modules/Module_C_PSA.R`: computes a probabilistic analysis including 
  1) Incremental cost-effectiveness ratios (ICERs) with probabilistic output,
  2)Cost-effectiveness acceptability curves (CEACs) and frontier (CEAF), 
  3) Expected Loss Curves (ELCs), and 
  4) Expected value of perfect information (EVPI).
  
3. Run analysis of the appropriate simulation sample size
  - `analysis_modules/Aux_Module_Convergence.R`: computes bootstrap distribution of prevalence estimates at each state 
  for varying simulation samples sizes and a fixed number of replications.  

### C. Code to Reproduce Publication Grade Plots

1. If you are interested in reproducing the publication-grade plots presented in the tutorial verify that you have installed additional packages:

  - `ggplot2`  
  - `viridisLite`
  - `patchwork`
  - `ggrepel`  
  - `gridExtra`
  - `ellipse`  
  - `ggview`   
  - `dampack`
  - `scales`

## Structure of Repository: 
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes

DES_Tutorial_beta/

├── DES_Tutorial.Rproj

├── install_deps.R

├── renv/ (ignore)

├── renv.lock

├── .Rhistory

├── .Rprofile

├── README.R

├── Inputs/

    └── all_cause_mortality.rda
    └── LifeTable_USA_Mx_2015.csv    
<<<<<<< Updated upstream

├── Lite_Modules/

    └── DES_Sick_Sicker_progressive.R
    └── DES_Sick_Sicker_recurrence.R
    └── Module_A_CEA.R
    └── Module_B_Epi_Outcomes.R
    └── Aux_Module_Cont_Time_Trace.R

├── Parallelized_Modules/

    └── Module_C_Probabilistic_Analysis.R
    └── Aux_Module_Convergence.R

├── Plotting_Modules/

    └── Module_CEA_wPlots.R
    └── Module_Epi_Outcomes_wPlots.R
    └── Module_Probabilistic_Analysis_wPlots.R
    └── Aux_Module_Cont_Time_Trace_wPlots.R
    └── Aux_Module_Convergence_wPlots.R

├── Manuscript/

    └── manuscript.qmd (.docx; .pdf)
    └── appendix.qmd   (.docx; .pdf)
    └── bibliography.bib
    └── apa-single-spaced.csl
    └── figures/

└── R/

    └── DES_functions.R
    └── Functions.R

## Dependencies by Module

![](table_deps.png)

# Workflow for collaborators
=======
    
├── sim_modules/
   
    └── DES_Sick_Sicker_progressive.R
    
    └── DES_Sick_Sicker_recurrence.R
    
├── analysis_modules/
    
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

## Workflow for collaborators
<<<<<<< Updated upstream
>>>>>>> Stashed changes

If you wish to contribute to this repository please follow these instructions to guarantee that you are using the same version of packages used by the developer.
=======
>>>>>>> Stashed changes

1. Run in RStudio Terminal:

`Rscript install_deps.R`

This script will install the package `renv`, which restores the 'locked' versions of all dependencies used by the developer. 
To do so, `renv` will install any missing packages and re-install versions of packages to match the version used by the developer. 
The script returns a list of packages that were not successfully installed and will require individual troubleshooting.

If all goes well, your environment should now match exactly the environment of the developer.

## Troubleshooting Package installation

<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
2. If you run into incompatibilities or mysterious errors running `install.deps.R`, verify that R version 3.5.0 or newer is installed.
=======
=======
=======
>>>>>>> Stashed changes
## Troubleshooting Package installation

>>>>>>> Stashed changes
If you run into incompatibilities or mysterious errors running `install.deps.R`, verify that R version 3.5.0 or newer is installed.
>>>>>>> Stashed changes

You can try to install required packages manually, run the below code on a script or from the R console:

    # Install `renv`

    if (!requireNamespace("renv", quietly = TRUE)) {
      install.packages("renv")
    }

    # Required packages
    required_pkgs <- c(
      "data.table"  ,   # to manipulate data
      "dplyr"       ,   # to manipulate data

<<<<<<< Updated upstream
      "ggplot2"     ,   # to visualize data
      "patchwork"   ,   # for combining ggplot2 figures
      "viridis"     ,   # color palettes
      "ellipse"     ,   # drawing ellipses and ellipse-like confidence regions
      "ggrepel"     ,   # extra geoms for ggplot
=======
required_pkgs <- c(
  "data.table"  ,   # to manipulate data
  "dplyr"       ,   # to manipulate data
  "reshape2"    ,   # to manipulate data
  
  "ggplot2"     ,   # to visualize data
  "ggrepel"     ,   # to visualize data
  "gridExtra"   ,   # to visualize data
  "ellipse"     ,   # to visualize data
  "ggview"      ,   # save plots
  
  "scales"      ,   # for dollar signs and commas
  "patchwork"   ,   # for combining ggplot2 figures
  "dampack"     ,   # for plots of CEA outcomes 
  
  "doParallel"  ,   # parallel processing
  "parallel"    ,   # parallel processing
  "foreach"     ,   # parallel processing
  
  "stats"       ,   # essential statistical functions
  "MethylCapSig",   # has nice multivariate lognormal random variable generator
  "survival"    ,   # core survival analysis routines
  "flexsurv"    ,   # flexible parametric survival models and multistate models
  
  "devtools"    ,   # to install packages from GitHub
  
  "abind"       ,   # combine multi-dimensional arrays
  "matrixStats"     # functions operating on rows and columns of matrices
)
>>>>>>> Stashed changes

      "doParallel"  ,   # parallel processing
      "parallel"    ,   # parallel processing
      "foreach"     ,   # parallel processing
      "abind"       ,   # multidimensional arrays
      
      "MethylCapSig",   # sample from multivariate lognormal
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

## Contact me at

<<<<<<< Updated upstream
mlopezme@stanford.edu and mlopezme@gmail.com
=======
[mlopezme\@stanford.edu](mailto:mlopezme@stanford.edu){.email} and [mlopezme\@gmail.com](mailto:mlopezme@gmail.com){.email}
>>>>>>> Stashed changes
