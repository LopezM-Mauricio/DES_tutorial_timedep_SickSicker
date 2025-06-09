#------------------------------------------------------------------------------#
# Post-processing Module A: Compute a Cost-effectiveness Analysis

#* Authors: 
#* - XXXXX
#* - XXXXX
#* - XXXXX
#------------------------------------------------------------------------------#
#* Please cite the article(s) when using this code
#* XXXXX open article
#* doi: XXXXX
#* bibtex: XXXXX
#------------------------------------------------------------------------------#
#* To program this tutorial we used:
#* R version 4.2.1 (2022-06-23)
#* Platform: aarch64-apple-darwin20 (64-bit)
#* Running under: Mac OS 12.5
#* RStudio: Version 2022.07.1+554 
#------------------------------------------------------------------------------#
# Description -----
#* This Module executes a Discrete Event Simulation (DES) of a Sick-Sicker model 
#* under four different strategies and computes a cost-effectiveness analysis.
#* Strategies are:
#* - Standard of Care (SoC): best available care for the patients with the 
#*   disease. This scenario reflects the natural history of the disease 
#*   progression.
#* - Strategy A: treatment A is given to patients in the Sick and Sicker states, 
#*   but only improves the quality of life of those in the Sick state.
#* - Strategy B: treatment B is given to all sick patients and reduces disease 
#*   progression from the Sick to Sicker state.
#* - Strategy AB: This strategy combines treatment A and treatment B. The disease 
#*   progression is reduced, and individuals in the Sick state have an improved 
#*   quality of life.
#* CEA Outcomes:
#*  ICERs
#*  CEA Frontier  
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
  n_sim            = 1e4   ,                                             
  
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
str(l_params)
#  Reference l_params
l_params_ref          <- l_params
#------------------------------------------------------------------------------#
## Run Simulation for all Strategies
#------------------------------------------------------------------------------#
v_strategies_names <- c( "SoC", "A", "B", "AB" )
#* since the effect on transition dynamics from SoC == A, and B == AB, 
#* we can just execute the simulation for "SOC" and "B"
v_strategies       <- c("SoC","B")        
# Set seed for reproducibility 
n_seed             <- 2
#Initialize list of event histories
l_event_history_strategies   <- list() 

# For-loop running simulation under each strategy
for(i in seq_along(v_strategies)){        
  # Initialize dt_event_history for each strategy
  dt_event_history    <- data.table()
  # Reset initial params for each iteration 
  l_params            <- l_params_ref
  # Update strategy
  l_params$strategy   <- v_strategies[i]
  # Common random numbers 
  set.seed(n_seed)
  # Run simulation under i^th strategy
  while(nrow(l_params$dt_next_event) > 0){                                         
    l_params$dt_next_event <- sim_next_event(l_params)
  }
  # Store dt_event_history for i^th strategy
  l_event_history_strategies[[i]] <-as.data.table(dt_event_history)     
}
names(l_event_history_strategies) <- v_strategies

# Duplicating event histories for equivalent strategies
# "SoC" and "A"
l_event_history_strategies$A  <- l_event_history_strategies$SoC
# "B" and "AB"
l_event_history_strategies$AB <- l_event_history_strategies $B
View(l_event_history_strategies)  

#------------------------------------------------------------------------------#
# Cost-Effectiveness Analysis (CEA) ----
#------------------------------------------------------------------------------#
# CEA Parameters ----
# -------------------------- #

#* We store the parameters necessary for the cost-effectiveness analysis
#* on a list `l_cea_params`.
#* Note we store the list of event histories in `l_cea_params`.


l_cea_params<- list(
  # list of event histories for each strategy
  l_event_history_strategies = l_event_history_strategies,
  
  ## Discounting factors ----
  d_c      = 0.03, # annual discount rate for costs
  d_e      = 0.03, # annual discount rate for QALYs
  
  ## State rewards ----
  ### Costs ----
  c_H      = 2000  , # annual cost of being Healthy
  c_S1     = 4000  , # annual cost of being Sick
  c_S2     = 15000 , # annual cost of being Sicker
  c_D      = 0     , # annual cost of being dead
  dc_trtA  = 12000 , # annual cost of receiving treatment A
  dc_trtB  = 13000 , # annual cost of receiving treatment B
  ### Utilities ----
  u_H      = 1     , # annual utility of being Healthy
  u_S1     = 0.75  , # annual utility of being Sick
  u_S2     = 0.5   , # annual utility of being Sicker
  u_D      = 0     , # annual utility of being dead
  u_trtA   = 0.95  , # annual utility when receiving treatment A
  du_trtA  = 0.95-0.75,
  ### Transition rewards ----
  du_HS1   = -0.01 , # disutility when transitioning from Healthy to Sick
  dc_HS1   = 1000  , # increase in cost when transitioning from Healthy to Sick
  dc_D     = 2000    # increase in cost when dying
  
)

# ----------------------------------- #
# Compute CEA ----
# ----------------------------------- #
## ICERs ----
#* `cea_fn` computes the cea analysis given a list of event histories, 
#* that correspond to a strategy. The details of the function are explained in the
#* manuscript subsection: Post-processing Module A: Compute a Cost-effectiveness Analysis
#* The output is the result of a competing choice analysis

table_cea<- cea_fn(l_cea_params = l_cea_params)
table_cea
## Plot CEA frontier -----
#* Function included in "R/Functions.R"; depends on the `ggplot2`  and `ggrepel` packages.
#* The latest version can be found in `dampack` package

plot.icers(table_cea, label = "all", txtsize = 16) +
  expand_limits(x = max(table_cea$Effect) + 0.1) +
  theme(legend.position = c(0.8, 0.2))


# ----------------------------------- #
# ----------------------------------- #



