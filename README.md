# DES_Tutorial_beta
Beta Release of DES Tutorial. 
V.1.04.1.2025

Please cite the preprint article: https://doi.org/10.1101/2025.05.15.25327635

@article {Lopez-Mendez2025.05.15.25327635,
	author = {Lopez-Mendez, Mauricio and Goldhaber-Fiebert, Jeremy D. and Alarid-Escudero, Fernando},
	title = {A Tutorial on Discrete Event Simulation Models in R Using a Cost-Effectiveness Analysis Example},
	elocation-id = {2025.05.15.25327635},
	year = {2025},
	doi = {10.1101/2025.05.15.25327635},
	publisher = {Cold Spring Harbor Laboratory Press},
	abstract = {Discrete Event Simulation (DES) is a flexible and computationally efficient approach for modeling diverse processes; however, DES remains underutilized in healthcare and medical decision-making due to a lack of reliable and reproducible implementations. We developed an open-source DES framework in R to simulate individual-level state-transition models (iSTMs) in continuous time accounting for treatment effects, time dependence on state residence, and age-dependent mortality.Our DES implementation employs a modular and easily adaptable structure, with each module corresponding to a unique transition between health states. To simulate the evolution of the process (i.e., individual state transitions), we adapted the next-reaction algorithm from the stochastic chemical reactions literature. Simulation-time dependence (age-dependent mortality) and state residence time dependence (transition from Sick to Sicker) are seamlessly incorporated into the DES framework via validated non-parametric and parametric sampling routines (e.g., inversion method) of event times. Treatment effects are integrated as scaling factors of the hazard functions (proportional hazards). We illustrate the framework{\textquoteright}s benefits by implementing the Sick-Sicker Model and conduct a cost-effectiveness analysis and probabilistic sensitivity analysis. We also obtain epidemiological outcomes of interest from the DES output, such as disease prevalence, survival probabilities, and distributions of state-specific dwell times. Our DES framework offers a reliable and accessible alternative that enables the simulation of more realistic dynamics of state-transition processes at potentially lower implementation and computational costs than discrete time iSTMs.Competing Interest StatementThe authors have declared no competing interest.Funding StatementDr. Goldhaber-Fiebert and Dr. Lopez-Mendez have no relevant funding to declare. Dr. Alarid-Escudero was supported by grants U01-CA253913, and U01-CA265750 from the National Cancer Institute (NCI) as part of the Cancer Intervention and Surveillance Modeling Network (CISNET).Author DeclarationsI confirm all relevant ethical guidelines have been followed, and any necessary IRB and/or ethics committee approvals have been obtained.YesI confirm that all necessary patient/participant consent has been obtained and the appropriate institutional forms have been archived, and that any patient/participant/sample identifiers included were not known to anyone (e.g., hospital staff, patients or participants themselves) outside the research group so cannot be used to identify individuals.YesI understand that all clinical trials and any other prospective interventional studies must be registered with an ICMJE-approved registry, such as ClinicalTrials.gov. I confirm that any such study reported in the manuscript has been registered and the trial registration ID is provided (note: if posting a prospective study registered retrospectively, please provide a statement in the trial ID field explaining why the study was not registered in advance).YesI have followed all appropriate research reporting guidelines, such as any relevant EQUATOR Network research reporting checklist(s) and other pertinent material, if applicable.YesAll data and code produced are available online at https://github.com/LopezM-Mauricio/DES_Tutorial https://github.com/LopezM-Mauricio/DES_Tutorial},
	URL = {https://www.medrxiv.org/content/early/2025/05/16/2025.05.15.25327635},
	eprint = {https://www.medrxiv.org/content/early/2025/05/16/2025.05.15.25327635.full.pdf},
	journal = {medRxiv}
}

# Structure of Repository: 

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
