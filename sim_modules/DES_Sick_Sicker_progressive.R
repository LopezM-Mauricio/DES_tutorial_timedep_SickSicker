#------------------------------------------------------------------------------#
# DES of Sick Sicker Model, 
# purely progressive H->S1->S2->D, H->D, S1->D.
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

#* This code implements a Discrete Event Simulation (DES) of a purely progressive 
#* Sick-Sicker model. 
#* Structure: Code is structured in modules as described in Lopez-Mendez et al 2025.
#* Complexity included in the model: 
#* Time Dependence: the DES incorporate two types of type dependency
#* A. Simulation time-dependence (Age specific mortality)
#* B. State residency time-dependence (Transition S1 -> S2)
#*Treatment effect on transition  S1->S2
#------------------------------------------------------------------------------#
# Initial setup ---- 
rm(list = ls())    # remove any variables in R's memory 
gc()               # clean working memory
#------------------------------------------------------------------------------#
## Load packages ----
library("data.table"  )
library("dplyr"       )

## Load DES-specific and utility functions ----
source("R/DES_functions.R") # DES specific 
source("R/Functions.R")     # General utility functions

#------------------------------------------------------------------------------#
# Initialize Parameters ----

# Simulation Parameters
## Simulation Time horizon 
time_init        <-   25    
time_end         <-   100  
time_horizon     <- c(time_init,time_end)
# Simulation sample size
n_sim            <-   1e3                              # Simulation sample size

## Model State Space and Transitions    
v_states         <- c("H",                             # Healthy (H)
                      "S1",                            # Sick (S1)
                      "S2",                            # Sicker (S2)
                      "D")                             # Dead (D)
## Indexed Adjacency Matrix
# Zero entries indicate no transition between pairs of states, 
# Nonzero positive entries index possible transitions in arbitrary order
m_Tr_r           <- matrix(c(0,1,0,2,
                             3,0,4,5,
                             0,0,0,6,
                             0,0,0,0), 
                           nrow = 4, ncol = 4, byrow = TRUE,
                           dimnames = list(c("H", "S1","S2","D"),
                                           c("H", "S1","S2","D")))
## Indexed transitions, tabular
dt_trans_keys    <- data.table()
dt_Tr_r          <-  as.data.table(m_Tr_r)
dt_Tr_r[, from    := v_states]
dt_trans_keys    <- as.data.table(melt(dt_Tr_r, id.vars = "from",
                                       measure.vars = c("H" ,"S1" ,"S2" ,"D"),
                                       variable.name = "to",
                                       value.name = "transition"))
dt_trans_keys    <- dt_trans_keys[transition !=0,][order(transition)]

## Transition Rates 

### Mortality, specify subfolder in project to load Life Tables
lifetable             <- "Inputs/all_cause_mortality.Rda"
load(lifetable)
dt_bckgrd_mortality   <- as.data.table(all_cause_mortality)
dt_bckgrd_mortality   <- dt_bckgrd_mortality[,death_rate :=  Rate][Year== 2015]
### Mortality Hazard Ratio while being S1
hr_S1                 <-    3
### Mortality Hazard Ratio while being S2
hr_S2                 <-    10

### Other transitions
### H  -> S1, constant rate ~Exp(r_HS1)
r_HS1                 <-    0.15
### S1 -> H,  constant rate ~Exp(r_S1H)
r_S1H                 <-    0.5
### S1 -> S2, state-residency time dependent transition 
# tau ~ Weibull(r_S1S2_scale_ph, r_S1S2_shape)
# PH parameterization
r_S1S2_scale_ph       <-    0.08
r_S1S2_shape          <-    1.1
# Strategy 
# Can be any of "SoC", "A", "B", "AB". Only B and AB, affect transitions rates
strategy              <-    "SoC"              
## Effectiveness of strategy B and strategy AB
# hazard ratio of becoming Sicker when Sick under treatment B or AB
hr_S1S2_trtB          <- 0.6  

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Useful Data tables
# Initialize `dt_baseline`        
age_init         <- 25                          # Initial age of the cohort
## Sex group mortality rate
Sex_grp          <- "Total"  
## Initial state
state_init       <- "H"                                                 
# Init dt
dt_baseline               <- data.table(
  # Covariate info
  ID                  = 1:n_sim,                # Numerical identifiers
  Age                 = age_init,               # Age at start
  Sex                 = Sex_grp,                # Sex, only if specified
  # Event tracking          
  from                = state_init,             # Origin state, first row everyone starts in H, update accordingly at each event. 
  to                  = NA,                     # Destination state, first row everyone is NA, update accordingly at each event. 
  T_n                 = time_horizon[1],        # Simulation time, records time of events, first row all entries are zero, but will be updated to the time when transition/event  occurs 
  T_start             = time_horizon[1],        # State-residence time start, first row all entries are zero, but will be updated to the transition that occurs next
  T_stop              = 0,                      # State-residence time stop, first row all entries are zero, but will be updated to the transition that occurs next
  tau                 = 0,                      # Time-increment, first row all entries are zero, updated after each event
  status              = NA,                     # Occurrence/censoring status, transition-specific indicator, updated after each event. 
  Event_num           = 1                       # Event number, at least 1 event from H-> D over time-horizon, updated after each event.
)

#------------------------------------------------------------------------------#
# Simulation (purely progressive) ----
# Initialize dt_next_event
dt_next_event      <- data.table()
# Initialize long-format dt_event_history 
dt_event_history   <- data.table()
# Set seed for reproducibility 
n_seed             <- 2
set.seed(n_seed)

# ---- #
#  Generate state-specific life-tables 
## Inputs:   dt_bckgrd_mortality , and  hr_S1, hr_S2
## 2.1 Create l copies of the original life table, l corresponds to the number of states with a transition to dead 
## 2.2 Apply state-specific hazard ratios to the original life table to obtain state-specific life tables.
## Mortality life tables
### Mortality from Sick
dt_mortality_Sick        <- copy(dt_bckgrd_mortality)
dt_mortality_Sick[, death_rate := death_rate*hr_S1]
### Mortality from Sicker
dt_mortality_Sicker      <- copy(dt_bckgrd_mortality)
dt_mortality_Sicker[, death_rate := death_rate*hr_S2]

#------------------------------------------------------------------------------#
# Transitions between Health States 
#------------------------------------------------------------------------------#
# General Structure of the purely progressive Sick-Sicker Model:
#   |-> Module  1: Transitions from Healthy
#       |-> Module Set up
#       |-> Task 1. Possible Transitions from Healthy 
#       |-> Task 2. Sample latent arrival times starting from Healthy 
#           |-> Submodule  1.1: Transitions Healthy -> Sick 
#           |-> Submodule  1.2: Transitions Healthy -> Dead 
#       |-> Task 3: Predict the next state 
#       |-> Task 4: Update important variables
#   |-> Module  2: Transitions from Sick 
#       |-> Module Set up
#       |-> Task 1. Possible Transitions from Sick 
#       |-> Task 2. Sample latent arrival times starting from Sick 
#           |-> Submodule  2.1: Transitions Sick  ->  Healthy
#           |-> Submodule  2.2: Transitions Sick  ->  Sicker
#           |-> Submodule  2.3: Transitions Sick  ->  Dead 
#       |-> Task 3: Predict the next state 
#       |-> Task 4: Update important variables
#   |-> Module  3: Transitions from Sicker 
#       |-> Module Set up
#       |-> Task 1. Possible Transitions from Sicker 
#       |-> Task 2. Sample latent arrival times starting from Sicker 
#           |-> Submodule  3.1: Transitions Sicker  ->  Dead 
#       |-> Task 3: Predict the next state 
#       |-> Task 4: Update important variables
#------------------------------------------------------------------------------#
#  Module  1: Transitions from Healthy ----
#------------------------------------------------------------------------------#
## Module Set up ----
# Objectives:
# 1. Generate  sub_tables to be used as inputs for each state-specific Module
## Inputs: dt_baseline, or dt_main;
## 1.1 At any given time determine how many people are in each non-absorbing state using dt_main
## 1.2 Split the dt_main into k different sub-tables, one sub_dt for each of k non-absorbing states
## 1.3 Store a list of state-specific sub-tables in a list, l_subsets_main_dt.
#------------------------------------------------------------------------------#

## Read dt
dt_main  <- copy(dt_baseline) # For first event

## Determine set of current states
v_current_state       <- unique(dt_main$from) 
v_current_state_names <- paste("from", v_current_state, sep = "_")

## Subset dt_main based on origin state(s)
# For the purely progressive model and everyone starting at H, 
# this step is trivial
l_subsets_main_dt <- list()
for( i in seq_along(v_current_state)){
  l_subsets_main_dt[[i]] <-  as.data.table(dt_main[from == v_current_state[i], ])
}
names(l_subsets_main_dt) <- v_current_state_names

# ------------------------------------------ #
## Task 1. Possible Transitions from H ----
# Objective: Identify possible transitions out of the current state (H) 
# Inputs: l_subsets_main_dt$from_H; m_Tr_r; dt_trans_keys,
# ------------------------------------------ #
# Temporary copy of dt_subset from H
dt_subset        <- l_subsets_main_dt$from_H
current_state    <- unique(dt_subset$from)     # current state
# Number of possible transitions out of Healthy
transition       <- as.numeric(m_Tr_r[current_state, m_Tr_r[current_state,] !=0])   
# Expand temp to have one row for each possible transition for each ID
ID               <- unique(dt_subset$ID)# unique IDs
setkey(dt_subset, ID)      
## CJ is like expand.grid  for data.tables
dt_subset_long   <- merge(dt_subset,CJ(ID, transition), all.x = F)        
dt_subset_long[, c("from","to") := NULL]  

## Assign "from" and "to" state using dt_trans_keys
# Indexed join in data.table using "on" sets the Key automatically
dt_subset_long   <- dt_subset_long[dt_trans_keys, on = "transition", nomatch = 0L]  
dt_subset_long   <- dt_subset_long[order(ID, transition)]

# ------------------------------------------- #
## Task 2: Sample latent Event-times ----
# ------------------------------------------- #
### Sub-module 1-1: Transition #1  (H -> S1) ----
# ------------------------------------------- #

# Objective: Sample time of occurrence for transition H-> S1 
# for every person currently in H.
# Inputs: dt_subset_long; r_HS1
# ------------------------------------------- #
# Constant annual rate of becoming Sick when Healthy
# Sample from an exponential distribution
if(r_HS1 != 0){   # Case when the transition is allowed
  dt_subset_long[transition == 1, T_stop := rexp(.N, rate = r_HS1) + T_n]
} else {        # Case when transition is not allowed, assign infinite time
  dt_subset_long[transition == 1, T_stop := Inf]
}
# ------------------------------------------ #
### Sub-module 1-2: Transition #2  (H -> D) ----
# ------------------------------------------- #

# Objective: Sample time of occurrence for transition H-> D 
# for every person currently in H.
# Inputs: obtain_probs_des; dt_time2death_probs;dt_time2death_probs; 
# dt_bckgrd_mortality; dt_subset_long
# ------------------------------------------ #
# Obtain vector of probabilities of dying from life table
dt_time2death_probs <- as.data.table(obtain_probs_des(dt_bckgrd_mortality)) 
# Only keep  rows for which the transition is possible,i.e., trans #2     
dt_trans2           <- dt_subset_long[transition == 2]               

# Set key in dt_trans2
setkey(dt_trans2, Age, Sex)                                                
# Set key in dt_time2death_probs
setkey(dt_time2death_probs, Age, Sex)                                       
# Indexed merge (all.x == T performs left-join)
dt_2sample_time2death <- merge(dt_trans2[, Age := ceiling(Age)], dt_time2death_probs, all.x = TRUE)

# Generate a matrix of probabilities from rows that matched and merged.
m_probs <- as.matrix(dt_2sample_time2death[, .SD, .SDcols = patterns("^DP")])

# Sample time to death using nonparametric sampling, using `nps_nhppp` with a uniform correction
dt_trans2[, T_stop := nps_nhppp(m_probs = m_probs, correction = "uniform")]

# Write into dt_subset_long
## guarantee that the order is correct when copying
setkey(dt_subset_long, ID, T_n)
setkey(dt_trans2, ID, T_n)
dt_subset_long[transition == 2, T_stop := dt_trans2$T_stop]

# --------------------------------------- #
## Task 3: Predict the Next Event from H ----
# Objectives: Determine which transition out of H occurs first
# Inputs: dt_subset_long
# --------------------------------------- #
# Create an indicator for minimum value in T_stop by ID 
dt_subset_long[, status := as.numeric(T_stop == min(T_stop)), by = ID]
# Recover the inter-event time, tau
dt_subset_long[,tau := T_stop - T_start]
#------------------------------------------------------------------------------#
# Since everyone starts in the Healthy state and this is a progressive model with 
# no recurrence,we can skip the evaluation of the transitions out of 
# the Sick and Sicker states in the first event.
#------------------------------------------------------------------------------#
# -------------------------------------- #
## Task 4: Update Important Variables ----
# Objectives: update important variables in the simulation
# including simulation time, age, event histories.
#  Updates differ  for the global event history, dt_event_history,
#  and for the data table used for the next event, dt_next_event
# Inputs: dt_subset_long, dt_event_history, dt_next_event
# -------------------------------------- #
# Update simulation time
dt_subset_long[status == T,T_n := T_stop]
# Update age
dt_subset_long[status == T,Age := round(T_n)]
# Update dt_event_history 
dt_event_history  <- rbind(dt_event_history, dt_subset_long) 
# Guarantee consistent ordering
setkey(dt_event_history, ID, T_n, Event_num, transition)

# Update dt_next_event
dt_next_event <- copy(dt_event_history)
setkey(dt_next_event,ID, T_n, Event_num, transition)
# Keep only transitions that occurred
dt_next_event <- dt_next_event[status == TRUE] 
# only work with the last event, ignore previous events
dt_next_event <- dt_next_event[dt_next_event[, .I[Event_num == max(Event_num, na.rm = TRUE)], by = ID]$V1]

# Handling Events beyond the time horizon
# Those who are alive beyond 100 years or older, assign death status 
# (should me a negligible number with the baseline parameters)
dt_next_event[to != "D"  & Age >= time_horizon[2]-1, to :="D"]

# discard in the local data table those who died, 
# for next event we only require those who remain alive
dt_next_event <- dt_next_event[to != "D"]

# Reset and update variables of interest
dt_next_event[,':=' (T_start = T_n,
                     tau = 0,
                     status  = NA,
                     Event_num = Event_num + 1,
                     from = to)]
dt_next_event[,':=' (to = NA,
                     T_stop = 0,
                     transition = NULL)]

# ---------------------------------------------------------------------------------------- #
# 2nd Event  
# For a purely progressive model starting everyone in Healthy: 
# For the next event after the first the only relevant people to
# keep track of in the simulation are those who transition to Sick
# ---------------------------------------------------------------------------------------- #
# Module 2: Transitions from Sick ----
# ---------------------------------------------------------------------------------------- #
## Module Set up ----
# Objectives of the Module:
# 1. Generate data subtables to be used as inputs for  each state-specific Module
## Inputs: dt_baseline, or dt_main;
## 1.1 At any given time determine how many people are in each non-absorbing state using dt_main
## 1.2 Split the dt_main into k different sub-tables, one sub_dt for each of k non-absorbing states
## 1.3 Store a list of state-specific sub-tables in a list, l_subsets_main_dt.
# ---------------------------------------------------------------------------------------- #
## Read dt_next_event
dt_main  <- dt_next_event 
## Determine set of current states
v_current_state       <- unique(dt_main$from) 
v_current_state_names <- paste("from", v_current_state, sep = "_")

## Subset dt_main based on origin state(s)
# For the purely progressive model and everyone starting at H, 
# this step is trivial since everyone is in the Sick state or died after the first event. 
l_subsets_main_dt <- list()
for( i in seq_along(v_current_state)){
  l_subsets_main_dt[[i]] <-  as.data.table(dt_main[from == v_current_state[i], ])
}
names(l_subsets_main_dt) <- v_current_state_names
# ----------------------------------- #
## Task 1: Possible Transitions from Sick ----
# Objective: Identify possible transitions out of the current state (Sick) 
# Inputs: l_subsets_main_dt$from_S1; m_Tr_r; dt_trans_keys,
# ----------------------------------- #
# Temporary copy of dt_subset from S1
dt_subset        <- l_subsets_main_dt$from_S1
current_state    <- unique(dt_subset$from)     # current state
# Number of possible transitions out of S1
transition       <- as.numeric(m_Tr_r[current_state, m_Tr_r[current_state,] !=0])   
# Expand temp to have one row for each possible transition for each ID
ID               <- unique(dt_subset$ID)# unique IDs
setkey(dt_subset, ID)      
## CJ is like expand.grid  for data.tables
dt_subset_long   <- merge(dt_subset,CJ(ID, transition), all.x = F)        
dt_subset_long[, c("from","to") := NULL]  

## Assign "from" and "to" state using dt_trans_keys
# Indexed join in data.table using "on" sets the Key automatically
dt_subset_long   <- dt_subset_long[dt_trans_keys, on = "transition", nomatch = 0L]  
dt_subset_long   <- dt_subset_long[order(ID, transition)]

# ----------------------------------- #
## Task 2: Sample latent Event-times ----
# ----------------------------------- #
# ----------------------------------- #
### Sub-module 2-1: Transition #3 (S1 -> H)  ----
# Objective: Sample time of occurrence for transition S1-> H 
# for every person currently in S1
# Inputs: dt_subset_long; r_S1H
# ----------------------------------- #
# Enforce S1->H is not allowed (purely progressive)
# Sample from exponential distribution
r_S1H <- 0
if(r_S1H != 0){ # case when the transition is  allowed
  dt_subset_long[transition == 3, T_stop := rexp(.N, rate = r_S1H) + T_n]
}else {         # case when the transition is not  allowed
  dt_subset_long[transition == 3, T_stop := Inf]  
}
# ----------------------------------- #
### Sub-module 2-2: Transition #4 (S1 -> S2) ----
# Objective: Sample time of occurrence for transition S1-> H 
# for every person currently in S1
# Inputs: dt_subset_long; r_S1S2_scale_ph,r_S1S2_shape, hr_S1S2_trtB
# ----------------------------------- #
# State-residence time-dependent hazard of transition from Sick to Sicker
# Apply the treatment effect for strategy B or strategy AB
if (strategy == "B" || strategy == "AB"){
  r_S1S2_scale_ph <- r_S1S2_scale_ph * hr_S1S2_trtB    
}
# Reparameterize to match accelerated failure time specification  
r_S1S2_scale_aft <- r_S1S2_scale_ph ^(-1 / r_S1S2_shape) 

# Sample inter-event time from Weibull distribution
if(r_S1S2_shape*r_S1S2_scale_aft != 0){  # case when the transition is  allowed
  dt_subset_long[transition == 4, T_stop := rweibull(.N, shape = r_S1S2_shape, scale = r_S1S2_scale_aft) + T_n]
} else {                                # case when the transition is not  allowed
  dt_subset_long[transition == 4, T_stop := Inf]
}
# ----------------------------------- #
### Sub-module 2-2: Transition #5  (S1 -> D) ----
# Objective: Sample time of occurrence for transition S1-> D 
# for every person currently in S1.
# Inputs: obtain_probs_des; dt_time2death_probs;dt_time2death_probs; 
# dt_mortality_Sick; dt_subset_long
# ----------------------------------- #
# Obtain vector of probabilities of dying from life table
dt_time2death_probs_Sick              <- as.data.table(obtain_probs_des(dt_mortality_Sick))

# Only keep  rows for which the transition is possible,i.e., trans #5     
dt_trans5                             <- dt_subset_long[transition == 5]

# Set key in dt_trans5
setkey(dt_trans5, Age,Sex)  
# Set key in dt_time2death_probs_Sick
setkey(dt_time2death_probs_Sick, Age,Sex)                                                                                
# Indexed merge (all.x == T performs left-join)
dt_2sample_time2death_Sick          <- merge(dt_trans5[,Age := ceiling(Age)], 
                                             dt_time2death_probs_Sick, all.x = TRUE)    

# Generate a matrix of probabilities from rows that matched and merged.
m_probs                             <- as.matrix(dt_2sample_time2death_Sick[, .SD, .SDcols = patterns("^DP")])           

# Sample time to death using nonparametric sampling, using `nps_nhppp` with a uniform correction
dt_trans5[, T_stop := nps_nhppp(m_probs = m_probs, correction = "uniform")]
# Consistent ordering 
setkey(dt_subset_long,ID,T_n)
setkey(dt_trans5,ID,T_n)
# Write into dt_subset_long
dt_subset_long[transition == 5, T_stop := dt_trans5$T_stop]
dt_subset_long[transition == 5]
# ----------------------------------- #
## Task 3: Predict the Next Event from Sick ----
# Objectives: Determine which transition out of S1 occurs first
# Inputs: dt_subset_long
# ----------------------------------- #
# Create an indicator for minimum value in T_stop by ID 
dt_subset_long[, status := as.numeric(T_stop == min(T_stop)), by = ID]
# Recover the inter-event time, tau
dt_subset_long[,tau := T_stop - T_start]

# ----------------------------------- #
## Task 4: Update Important Variables ----
# Objectives: update important variables in the simulation
# including simulation time, age, event histories.
#  Updates differ  for the global event history, dt_event_history,
#  and for the data table used for the next event, dt_next_event
# Inputs: dt_subset_long, dt_event_history, dt_next_event
# ----------------------------------- #
# Update simulation time
dt_subset_long[status == T,T_n := T_stop]
# Update age
dt_subset_long[status == T,Age := round(T_n)]
# Update dt_event_history 
dt_event_history  <- rbind(dt_event_history, dt_subset_long) 
# Guarantee consistent ordering
setkey(dt_event_history, ID, T_n, Event_num, transition)

# Update dt_next_event
dt_next_event <- copy(dt_event_history)
setkey(dt_next_event,ID, T_n, Event_num, transition)
# Keep only transitions that occurred
dt_next_event <- dt_next_event[status == TRUE] 
# only work with the last event, ignore previous events
dt_next_event <- dt_next_event[dt_next_event[, .I[Event_num == max(Event_num, na.rm = TRUE)], by = ID]$V1]

# Handling Events beyond the time horizon
# Those who are alive beyond 100 years or older, assign death status
dt_next_event[to != "D"  & Age >= time_horizon[2]-1, to :="D"]
# discard those who died for next event
dt_next_event <- dt_next_event[to != "D"]

# #Reset and update variables of interest
dt_next_event[,':=' (T_start = T_n,
                     tau = 0,
                     status  = NA,
                     Event_num = Event_num + 1,
                     from = to)]
dt_next_event[,':=' (to = NA,
                     T_stop = 0,
                     transition = NULL)]

# ---------------------------------------------------------------------------------------- #
# 3rd Event  
#------------------------------------------------------------------------------#
# For a purely progressive model starting everyone in Healthy: 
# we can skip Submodule 1 (transitions from H)
# and Submodule 3 (transitions from Sicker) in the second event.
# In a purely progressive model with everyone starting in Healthy ,
# after the first event, nobody in the simulation stays Healthy.
# In a purely progressive model with everyone starting in Healthy ,
# after the second event, nobody in the simulation stays in Healthy nor Sick.
#------------------------------------------------------------------------------#
# Module 3: Transitions from Sicker ----
#------------------------------------------------------------------------------#
## Module Set up ----
# Objectives of the Module:
# 1. Generate data subtables to be used as inputs for  each state-specific Module
## Inputs: dt_baseline, or dt_main;
## 1.1 At any given time determine how many people are in each non-absorbing state using dt_main
## 1.2 Split the dt_main into k different sub-tables, one sub_dt for each of k non-absorbing states
## 1.3 Store a list of state-specific sub-tables in a list, l_subsets_main_dt.
# ---------------------------------------------------------------------------------------- #
## Read dt_next_event
dt_main  <- dt_next_event 
## Determine set of current states
v_current_state       <- unique(dt_main$from) 
v_current_state_names <- paste("from", v_current_state, sep = "_")

## Subset dt_main based on origin state(s)
# For the purely progressive model and everyone starting at H, 
# this step is trivial since everyone is in the Sicker state
l_subsets_main_dt <- list()
for( i in seq_along(v_current_state)){
  l_subsets_main_dt[[i]] <-  as.data.table(dt_main[from == v_current_state[i], ])
}
names(l_subsets_main_dt) <- v_current_state_names
# ----------------------------------- #
## Task 1: Possible Transitions from Sicker ----
# ----------------------------------- #
# Temporary copy of dt_subset from S2
dt_subset        <- l_subsets_main_dt$from_S2
current_state    <- unique(dt_subset$from)     # current state
# Number of possible transitions out of S2
transition       <- as.numeric(m_Tr_r[current_state, m_Tr_r[current_state,] !=0])   
# Expand temp to have one row for each possible transition for each ID
ID               <- unique(dt_subset$ID)# unique IDs
setkey(dt_subset, ID)      
## CJ is like expand.grid  for data.tables
dt_subset_long   <- merge(dt_subset,CJ(ID, transition), all.x = F)        
dt_subset_long[, c("from","to") := NULL]  

## Assign "from" and "to" state using dt_trans_keys
# Indexed join in data.table using "on" sets the Key automatically
dt_subset_long   <- dt_subset_long[dt_trans_keys, on = "transition", nomatch = 0L]  
dt_subset_long   <- dt_subset_long[order(ID, transition)]
# ----------------------------------- #
## Task 2: Sample latent Event-times -----
# ----------------------------------- #
### Submodule 3-1: Transition #6  (S2 -> D) ----
# Objective: Sample time of occurrence for transition S2-> D 
# for every person currently in S2.
# Inputs: obtain_probs_des; dt_time2death_probs;dt_time2death_probs; 
# dt_mortality_Sicker; dt_subset_long
# ----------------------------------- #
# Obtain vector of probabilities of dying from life table
dt_time2death_probs_Sicker              <- as.data.table(obtain_probs_des(dt_mortality_Sicker))

# Only keep  rows for which the transition is possible,i.e., trans #5     
dt_trans6                             <- dt_subset_long[transition == 6]

# Set key in dt_trans5
setkey(dt_trans6, Age,Sex)  
# Set key in dt_time2death_probs_Sick
setkey(dt_time2death_probs_Sicker, Age,Sex)                                                                                
# Indexed merge (all.x == T performs left-join)
dt_2sample_time2death_Sicker          <- merge(dt_trans6[,Age := ceiling(Age)], 
                                               dt_time2death_probs_Sicker, all.x = TRUE)    

# Generate a matrix of probabilities from rows that matched and merged.
m_probs                             <- as.matrix(dt_2sample_time2death_Sicker[, .SD, .SDcols = patterns("^DP")])           

# Sample time to death using nonparametric sampling, using `nps_nhppp` with a uniform correction
dt_trans6[, T_stop := nps_nhppp(m_probs = m_probs, correction = "uniform")]
# Consistent ordering 
setkey(dt_subset_long,ID,T_n)
setkey(dt_trans6,ID,T_n)
# Write into dt_subset_long
dt_subset_long[transition == 6, T_stop := dt_trans6$T_stop]
# ----------------------------------- #
## Task 3: Predict the Next Event from Sicker ----
# Objectives: Determine which transition out of S2 occurs first
# Inputs: dt_subset_long
# ----------------------------------- #
# Create an indicator for minimum value in T_stop by ID 
dt_subset_long[, status := as.numeric(T_stop == min(T_stop)), by = ID]
# Recover the inter-event time, tau
dt_subset_long[,tau := T_stop - T_start]

# ----------------------------------- #
## Task 4: Update Important Variables ----
# Objectives: update important variables in the simulation
# including simulation time, age, event histories.
# Updates differ  for the global event history, dt_event_history,
# and for the data table used for the next event, dt_next_event
# Inputs: dt_subset_long, dt_event_history, dt_next_event
# ----------------------------------- #
# Update simulation time
dt_subset_long[status == T,T_n := T_stop]
# Update age
dt_subset_long[status == T,Age := round(T_n)]
# Update dt_event_history 
dt_event_history  <- rbind(dt_event_history, dt_subset_long) 
# Guarantee consistent ordering
setkey(dt_event_history, ID, T_n, Event_num, transition)

# ----#
# Quick Validation ----
# ---- #
# Verify everyone follows progression H->S1->S2->D, or dies in between
dt_event_history[, from_to := paste(from, to, sep = "_")]
table(dt_event_history[status==TRUE]$from_to)
# sum of all who die must equal n_sim
(table(dt_event_history[status==TRUE]$from_to)[1]+
  table(dt_event_history[status==TRUE]$from_to)[3]+
  table(dt_event_history[status==TRUE]$from_to)[5]) == n_sim
# Note some might have reached the time-horizon before dying, 
# in this example it should be a negligible number
