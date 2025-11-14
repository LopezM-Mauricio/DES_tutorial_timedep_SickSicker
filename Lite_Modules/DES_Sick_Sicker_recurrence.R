#------------------------------------------------------------------------------#
#* MODULE 3 : DES of Sick Sicker Model with recurrence 
#* from Sick to Sicker S1 -> H
#------------------------------------------------------------------------------#
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
# Description ----

#* This code implements a Discrete Event Simulation (DES) of a Sick-Sicker model
#* with recurrence.
#* Structure: 
#* Initial parameters stored in a list , `l_params`.  
#* Modules are wrapped in the function `sim_next_event()`
#* which updates the `dt_next_event` and `dt_event_history`
#* The function `sim_next_event()` follows a modular structure : 

#------------------------------------------------------------------------------#
# General Structure and flow of the purely progressive Sick-Sicker Model
#------------------------------------------------------------------------------#
# while(nrow(dt_next_event)>0) 
# Execute

#   |-> Module  1: Currently Occupied Health States
#   |-> Module  2: Transitions between Health States
#      |-> Submodule 1: Transitions from Healthy
#       |-> Task 1. Possible Transitions from Healthy 
#       |-> Task 2. Sample latent arrival times starting from Healthy 
#           |-> 2.1: Transitions Healthy -> Sick 
#           |-> 2.2: Transitions Healthy -> Dead 
#       |-> Task 3: Predict the next state 
#       |-> Task 4: Update important variables
#      |-> Submodule 2: Transitions from Sick
#       |-> Task 1. Possible Transitions from Sick 
#       |-> Task 2. Sample latent arrival times starting from Sick 
#           |->  2.1: Transitions Sick  ->  Healthy
#           |->  2.2: Transitions Sick  ->  Sicker
#           |->  2.3: Transitions Sick  ->  Dead 
#       |-> Task 3: Predict the next state 
#       |-> Task 4: Update important variables
#      |-> Submodule 3: Transitions from Sicker
#       |-> Task 1. Possible Transitions from Sicker 
#       |-> Task 2. Sample latent arrival times starting from Sicker 
#           |-> 2.1: Transitions Sicker  ->  Dead 
#       |-> Task 3: Predict the next state 
#       |-> Task 4: Update important variables

# A while loop executes  `sim_next_event()`, which contains Module 1 and 2,  until everyone in the simulation has died

# Complexity included in the model: 
# Time Dependence: the DES incorporate two types of type dependency
# A. Simulation time-dependence (Age specific mortality)
# B. State residency time-dependence (Transition S1 -> S2)
# Treatment effect acting on transition rate S1->S2
#------------------------------------------------------------------------------#
# Initial setup ---- 
rm(list = ls())    # remove any variables in R's memory 
gc()               # clean working memory
#------------------------------------------------------------------------------#
## Load packages ----
library("data.table"  )
library("dplyr"       )

## Load supplementary functions ----
source("R/DES_functions.R") # DES specific 
source("R/Functions.R")     # General utility functions

#------------------------------------------------------------------------------#
# Initialize Parameters ----
# We want to pass parameters seamlessly through our function sim_next_event(), 
# for this we can declare initial parameters in a list. Here we pass them through
# a function `init_param` that will return such list `l_params` and also include the baseline data table
# `dt_baseline` as an object of that list, now named `dt_next_event`. 

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
  # Zero entries indicate no transition between pairs of states, 
  # Nonzero positive entries index possible transitions in arbitrary order
  m_Tr_r           = matrix(c(0,1,0,2,
                              3,0,4,5,
                              0,0,0,6,
                              0,0,0,0), 
                            nrow = 4, ncol = 4, byrow = TRUE,
                            dimnames = list(c("H", "S1","S2","D"),
                                            c("H", "S1","S2","D"))),
  # Transition Rates 
  ## Mortality 
  # we specify subfolder in project to load Life Tables
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
  # Can be any of "SoC", "A", "B", "AB". 
  # Only B and AB, affect transitions rates,
  strategy         =    "SoC",                                            
  
  # Tx Effects that affect state-transitions
  # Effectiveness of strategy B hazard ratio 
  # of becoming Sicker when Sick under treatment B
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
# Explore a single trajectory ----
# Output in long format
head(dt_event_history)
# Subset to events that occurred
head(dt_event_history[status == TRUE])

# Single trajectory
dt_event_history[status == TRUE & ID == 3]

# With the current setup, nsim = 1e5
# Sim person with ID:3, experienced the following path:
# H->S1->H->S1->S2->D (read from cols "from" and "to")
# Events occurred at the following times:
# 25.97768 -> 26.88503 -> 28.52374 -> 29.20695 -> 76.18308 (read from col "T_stop")
# Dwell times at each alive states:
# H(0.9776845), S1(0.9073458), H(1.6387102), S1(0.6832138), S2(46.9761233)  (read from col "tau")
#------------------------------------------------------------------------------#


