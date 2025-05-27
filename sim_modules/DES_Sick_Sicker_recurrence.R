#------------------------------------------------------------------------------#
#* MODULE 3 : DES of Sick Sicker Model with recurrence 
#* from Sick to Sicker S1 -> H
#------------------------------------------------------------------------------#
#* Authors: 
#* - Mauricio Lopez-Mendez <mlopezme@stanford.edu> (Stanford Health Policy)
#* https://orcid.org/0000-0002-3473-5457
#* - Jeremy Goldhaber-Fiebert <jeremygf@stanford.edu> (Stanford Health Policy)
#* https://orcid.org/0000-0002-3473-5457
#* - Fernando Alarid-Escudero <falarid@stanford.edu> (Stanford Health Policy)
#* https://orcid.org/0000-0001-5076-1172
#------------------------------------------------------------------------------#
#* Please cite the article(s) when using this code
#* MedRxiv open article
#* doi: https://doi.org/10.1101/2025.05.15.25327635 
#* bitex: 
# @article{lopez-mendez_tutorial_2025,
# title = {A {Tutorial} on {Discrete} {Event} {Simulation} {Models} in {R} {Using} a {Cost}-{Effectiveness} {Analysis} {Example}},
# url = {https://www.medrxiv.org/content/early/2025/05/16/2025.05.15.25327635},
# doi = {10.1101/2025.05.15.25327635},
# abstract = {Discrete Event Simulation (DES) is a flexible and computationally efficient approach for modeling diverse processes; however, DES remains underutilized in healthcare and medical decision-making due to a lack of reliable and reproducible implementations. We developed an open-source DES framework in R to simulate individual-level state-transition models (iSTMs) in continuous time accounting for treatment effects, time dependence on state residence, and age-dependent mortality.Our DES implementation employs a modular and easily adaptable structure, with each module corresponding to a unique transition between health states. To simulate the evolution of the process (i.e., individual state transitions), we adapted the next-reaction algorithm from the stochastic chemical reactions literature. Simulation-time dependence (age-dependent mortality) and state residence time dependence (transition from Sick to Sicker) are seamlessly incorporated into the DES framework via validated non-parametric and parametric sampling routines (e.g., inversion method) of event times. Treatment effects are integrated as scaling factors of the hazard functions (proportional hazards). We illustrate the framework’s benefits by implementing the Sick-Sicker Model and conduct a cost-effectiveness analysis and probabilistic sensitivity analysis. We also obtain epidemiological outcomes of interest from the DES output, such as disease prevalence, survival probabilities, and distributions of state-specific dwell times. Our DES framework offers a reliable and accessible alternative that enables the simulation of more realistic dynamics of state-transition processes at potentially lower implementation and computational costs than discrete time iSTMs.Competing Interest StatementThe authors have declared no competing interest.Funding StatementDr. Goldhaber-Fiebert and Dr. Lopez-Mendez have no relevant funding to declare. Dr. Alarid-Escudero was supported by grants U01-CA253913, and U01-CA265750 from the National Cancer Institute (NCI) as part of the Cancer Intervention and Surveillance Modeling Network (CISNET).Author DeclarationsI confirm all relevant ethical guidelines have been followed, and any necessary IRB and/or ethics committee approvals have been obtained.YesI confirm that all necessary patient/participant consent has been obtained and the appropriate institutional forms have been archived, and that any patient/participant/sample identifiers included were not known to anyone (e.g., hospital staff, patients or participants themselves) outside the research group so cannot be used to identify individuals.YesI understand that all clinical trials and any other prospective interventional studies must be registered with an ICMJE-approved registry, such as ClinicalTrials.gov. I confirm that any such study reported in the manuscript has been registered and the trial registration ID is provided (note: if posting a prospective study registered retrospectively, please provide a statement in the trial ID field explaining why the study was not registered in advance).YesI have followed all appropriate research reporting guidelines, such as any relevant EQUATOR Network research reporting checklist(s) and other pertinent material, if applicable.YesAll data and code produced are available online at https://github.com/LopezM-Mauricio/DES\_Tutorial https://github.com/LopezM-Mauricio/DES\_Tutorial},
# journal = {medRxiv},
# author = {Lopez-Mendez, Mauricio and Goldhaber-Fiebert, Jeremy D. and Alarid-Escudero, Fernando},
# year = {2025},
# note = {Publisher: Cold Spring Harbor Laboratory Press
#   \_eprint: https://www.medrxiv.org/content/early/2025/05/16/2025.05.15.25327635.full.pdf},
# }
#------------------------------------------------------------------------------#
#* To program this tutorial we used:
#* R version 4.2.1 (2022-06-23)
#* Platform: aarch64-apple-darwin20 (64-bit)
#* Running under: Mac OS 12.5
#* RStudio: Version 2022.07.1+554 
#------------------------------------------------------------------------------#
# Description ----

#* This code implements a Discrete Event Simulation (DES) of a Sick-Sicker model
#* with recurrence.
#* Structure: 
#* Initial parameters stored in a list and  
#* Modules 1 and 2 are wrapped in the function `sim_next_event()`
#* which updates the `dt_next_event` and `dt_event_history`
#* A while loop executes sim_next_event until everyone in the simulation has died
#* Complexity included in the model: 
#* Time Dependence: the DES incorporate two types of type dependency
#* A. Simulation time-dependence (Age specific mortality)
#* B. State residency time-dependence (Transition S1 -> S2)
#*Treatment effect acting on transition rate S1->S2
#------------------------------------------------------------------------------#
# Initial setup ---- 
rm(list = ls())    # remove any variables in R's memory 
gc()               # clean working memory
#------------------------------------------------------------------------------#
## Load packages ----
library("data.table"  )
library("dplyr"       )
library("tidyr"       )
library("reshape2"    )
library("ggplot2"     )
library("ggrepel"     )
library("gridExtra"   )
library("ellipse"     )
library("ggview"      )
library("scales"      )
library("patchwork"   )
library("dampack"     )
library("doParallel"  )
library("parallel"    )
library("foreach"     )
library("stats"       )
library("MethylCapSig")
library("survival"    )
library("flexsurv"    )
library("devtools"    )

## Load supplementary functions ----
source("R/DES_functions.R") # DES specific 
source("R/Functions.R")     # General utility functions

#------------------------------------------------------------------------------#
# Initialize Parameters ----
#* We want to pass parameters seamlessly through our function sim_next_event(), 
#* for this we can declare initial parameters in a list. Here we pass them through
#* a function `init_param` that will return such list `l_params` and also include the baseline data table
#* `dt_baseline` as an object of that list, now named `dt_next_event`. 

l_params <- init_params(
  # Sim. params
  # Time horizon
  time_init        = 25    ,
  time_end         = 100   ,
  # Sim. sample size
  n_sim            = 1e5   ,                                             
  
  # dt_baseline information        
  ## initial age
  age_init         = 25    ,                                             
  ## sex group
  Sex_grp          = "Total", 
  ## initial state
  state_init       = "H"   ,  
  
  # State Space      
  v_states         = c("H",                                              # Healthy (H)
                       "S1",                                             # Sick (S1)
                       "S2",                                             # Sicker (S2)
                       "D"),                                             # Dead (D)
  
  # Indexed Adjacency Matrix
  #* Zero entries indicate no transition between pairs of states, 
  #* Nonzero positive entries index possible transitions in arbitrary order
  m_Tr_r           = matrix(c(0,1,0,2,
                              3,0,4,5,
                              0,0,0,6,
                              0,0,0,0), 
                            nrow = 4, ncol = 4, byrow = TRUE,
                            dimnames = list(c("H", "S1","S2","D"),
                                            c("H", "S1","S2","D"))),
  # Transition Rates 
  ## Mortality 
  #* we specify subfolder in project to load Life Tables
  lifetable        = "Inputs/all_cause_mortality.Rda",
  ## Hazard Ratio while being S1
  hr_S1            =    3,
  ## Hazard Ratio while being S2
  hr_S2            =    10, 
  
  ## Other transitions
  ## H  -> S1, constant rate ~Exp(r_HS1)
  r_HS1            =    0.15,
  ## S1 -> H,  constant rate ~Exp(r_S1H)
  r_S1H            =    0.5,
  ## S1 -> S2, state-residency time dependent ~ Weibull(r_S1S2_scale_ph, r_S1S2_shape)
  ##* PH parameterization
  r_S1S2_scale_ph  =    0.08,
  r_S1S2_shape     =    1.1,
  
  # Strategy 
  #* Can be any of "SoC", "A", "B", "AB". 
  #* Only B and AB, affect transitions rates,
  strategy         =    "SoC",                                            
  
  # Tx Effects that affect state-transitions
  #* Effectiveness of strategy B hazard ratio 
  #* of becoming Sicker when Sick under treatment B
  hr_S1S2_trtB = 0.6  
)

#------------------------------------------------------------------------------#
# Simulation with Recurrence  ----
# ----------------------------------- #
# Initialize dt_next_event
dt_next_event      <- data.table()
# Initialize long-format dt_event_history 
dt_event_history   <- data.table()
# Set seed for reproducibility 
n_seed             <- 2
set.seed(n_seed)
# Max. number of events before triggering warning
# e.g., if more than 1,000 have occurred, trigger a warning message
max_num_events <- 1e3 

#------------------------------------------------------------------------------#
# Run while loop until everyone in the simulation has died
ptime             <- system.time({
  while(
    # stopping criteria: each row in dt_next_event represents
    # the number of people still alive in the simulation
    nrow(l_params$dt_next_event) >0
  ){  
    
    # Function updates next event for everyone who is still alive
    l_params$dt_next_event <- sim_next_event(l_params) 
    
    # Warning:
    if(unique(length(l_params$sim_next_event$Event_num)) > max_num_events){
      warning(
        paste0("Within a simulation time-horizon of"
               ,l_params$time_horizon[2]-l_params$time_horizon[1], "years, 
                the max number of events (i.e., transitions between health
                states) for at least one person in the simulation has 
                exceeded", max_num_events,".Reconsider if the simulation 
                should continue since these dynamics might not reflect 
                what you intended."
        )) 
    }
  }
})[3]
ptime/60
ptime

#------------------------------------------------------------------------------#
# Explore one trajectory

# Output in long format
head(dt_event_history)
# Subset to events that occurred
head(dt_event_history[status == TRUE])

# Single trajectory
dt_event_history[status == TRUE & ID == 3]

# With the current setup
# Sim person with ID:3, experienced the following path:
# H->S1->H->S1->S2->D (read from cols "from" and "to")
# Events occurred at the following times:
# 25.97768 -> 26.88503 -> 28.52374 -> 29.20695 -> 76.18308 (read from col "T_stop")
# Dwell times at each alive states:
# H(0.9776845), S1(0.9073458), H(1.6387102), S1(0.6832138), S2(46.9761233)  (read from col "tau")
#------------------------------------------------------------------------------#


