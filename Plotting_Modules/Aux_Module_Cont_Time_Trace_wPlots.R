#------------------------------------------------------------------------------#
# AUXILIARY MODULE A: CONTINUOUS TIME TRACE MATRIX AND TRACE PLOT 
# SICK SICKER MODEL WITH RECURRENCE

#* Authors: 
#* - 
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

#* In this Module we illustrate how to use functions 
#* `cont_time_trace_sick_sicker`, 
#* and  `discr_time_trace_sick_sicker`, 
#* to transform the  output of the Discrete Event Simulation (DES) of 
#* the Sick-Sicker model with or without recurrence into a continuous time 
#* trace matrix or a discrete time approximation, 
#* and how to produce a trace plot from the trace matrices.

#------------------------------------------------------------------------------#
# Initial setup ---- 
rm(list = ls())    # remove any variables in R's memory 
gc()               # clean working memory
#------------------------------------------------------------------------------#
## Load packages ----
library("data.table"  )
library("dplyr"       )
library("ggplot2"     )
library("viridis")
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
  n_sim            = 1e3   ,                                             
  
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

#View(dt_event_history)


# -------------------------------------------------------------------- #
# Trace Matrix and Trace Plots ----
# -------------------------------------------------------------------- #


## A. Continuous time trace matrix ----
#* We want to recover the average state occupancy of the simulated cohort from the individual 
#* state transition trajectories. This mapping is encoded in the function `cont_time_trace_sick_sicker`
#* that returns a trace matrix. Note that in contrast to a discrete time trace matrix, where rows
#* index cycles, in the continuous time trace matrix, rows index events. 


dt_CT_trace <- cont_time_trace_sick_sicker(dt_event_history = dt_event_history,
                                           v_names_states   = l_params$v_states,  
                                           time_init        = l_params$time_horizon[1],
                                           proportions      = TRUE, 
                                           rounding_factor  = 5,
                                           long            = F)


dt_CT_trace_long <- cont_time_trace_sick_sicker(dt_event_history = dt_event_history,
                                                v_names_states   = l_params$v_states,  
                                                time_init        = l_params$time_horizon[1],
                                                proportions      = TRUE, 
                                                rounding_factor  = 5,
                                                long            = T)

## B. Trace plot ----
#* We want to plot the average state occupancy over time as a traditional trace plot
#* We can do so using `dt_CT_trace` or using an approximation of `dt_CT_trace`. 

### B.1 Continuous time Plot ----

ct_trace_plot <- ggplot(dt_CT_trace_long, aes(x = T_n, y = Proportion, color = Health_State)) +
  geom_line(linewidth = 1.3) +

  # Titles and Labels
  xlab("Time (years)") +
  ylab("State Occupancy (Proportion)") +
  labs(title = "State Occupancy over Time", color = "Health States") +

  # Grid
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    minor_breaks = seq(0, 1, by = 0.05)
  ) +

  # Colors
  scale_color_viridis_d(option = "H", begin = 0, end = 1, direction = 1) +

  # Theme
  theme_minimal(base_size = 28) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 16, angle = 0, color = "black"),
    axis.text.y = element_text(size = 16, angle = 0, color = "black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(face = "italic")
  )

ct_trace_plot


### B.2 Discrete time approximation  ----
#* Many factors like sim sample size, large number of states, long time-horizons, 
#* can contribute to extremely large continuous time trace matrices,  `dt_CT_trace`.
#* In those cases it might be more practical to use an approximation for tasks like plotting the trace. 
#* The function `discr_time_trace_sick_sicker` provides such approximation by finding the event times
#* closest to a grid of fixed times points supplied by the user. For example, if we want to compare our DES implementation of the Sick Sicker model
#* with the discrete time, time dependent cSTM implementation we can supply the same time grid used for the cSTM and find 
#* the closest evenT time to each time point on the grid. 

# Discrete Time grid for a cSTM
cycle_length <- 1/12 # cycle length equal to one year (use 1/12 for monthly)
n_age_init   <- 25  # age at baseline
n_age_max    <- 100 # maximum age of follow up
n_cycles     <- (n_age_max - n_age_init)/cycle_length # time horizon, number of cycles

v_DT_grid        <- seq(l_params$time_horizon[1],l_params$time_horizon[2], length.out = n_cycles)


dt_DT_trace      <- discr_time_trace_sick_sicker(dt_CT_trace     = dt_CT_trace ,
                                                 v_DT_grid   = v_DT_grid   , 
                                                 long        = F    )


dt_DT_trace_long <- discr_time_trace_sick_sicker(dt_CT_trace = dt_CT_trace ,
                                                 v_DT_grid   = v_DT_grid , 
                                                 long        = T     )

# Discrete Time Trace Plot
dt_trace_plot <- ggplot(dt_DT_trace_long, aes(x = T_n, y = Proportion, color = Health_State)) +
  geom_line(linewidth = 1.3) +

  # Titles and Labels
  xlab("Time (years)") +
  ylab("State Occupancy (Proportion)") +
  labs(title = "State Occupancy over Time", color = "Health States") +

  # Grid
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    minor_breaks = seq(0, 1, by = 0.05)
  ) +

  # Colors
  scale_color_viridis_d(option = "H", begin = 0, end = 1, direction = 1) +

  # Theme
  theme_minimal(base_size = 28) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 16, angle = 0, color = "black"),
    axis.text.y = element_text(size = 16, angle = 0, color = "black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(face = "italic")
  )

dt_trace_plot

