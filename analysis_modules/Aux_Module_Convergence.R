#------------------------------------------------------------------------------#
# AUXILIARY MODULE B:  CONVERGENCE ANALYSIS
# SICK SICKER MODEL WITH RECURRENCE

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

#* We illustrate how to assess the convergence of the state occupancy estimator
#* as we increase the simulation sample size. State occupancy is recovered from 
#* the output of the DES by computing the proportion of simulated people in a given state
#* at a given event time (Monte Carlo integration). As we increase the simulation sample size,
#* the Monte Carlo standard error shrinks  around the estimator of 
#* state occupancy. 

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
  n_sim            = 1e2   ,                                             
  
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
# -------------------------------- #
# Distributed Implementation ------
#  ------------------------------- #

## Cluster Setup ----

# Num. of Cores
no_cores_all  <- parallel::detectCores()-1      # recommended when not sharing computing resources
cores_90pct   <- as.integer(.90 * no_cores_all) # aggressive  when shared resources
cores_85pct   <- as.integer(.85 * no_cores_all) # aggressive  when shared resources
cores_75pct   <- as.integer(.75 * no_cores_all) # aggressive  when shared resources
cores_50pct   <- as.integer(.5  * no_cores_all) # recommended when shared resources
cores_25pct   <- as.integer(.25 * no_cores_all) # recommended when shared resources


os_info <- Sys.info()
os <- os_info[1]
no_cores <- cores_85pct
if (os == "Darwin") { # Apple Kernel
  # Initialize cluster object
  cl <- parallel::makeForkCluster(no_cores)
  # Register clusters
  doParallel::registerDoParallel(cl)
}
if (os == "linux") {
  # Initialize cluster object
  cl <- parallel::makeCluster(no_cores)
  # Register clusters
  doParallel::registerDoMC(cl)
}
if (os == "windows") {  
  # Initialize cluster object
  cl <- parallel::makeCluster(no_cores)
  # Register clusters
  doParallel::registerDoParallel(cl)
}

## Run distributed implementation ----

# Set seed for reproducibility 
n_seed            <- 2
set.seed(n_seed)

# batch and sim sample size 
v_n_sim     <- c(1e2,1e3,1e4,1e5)   # iterate over
v_n_batch   <- c(5e1,1e2,1e3)       # fix one
k <- v_n_batch[1] 

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# WARNING: DO NOT RUN THE NEXT CHUNK UNLESS YOU HAVE ENOUGH COMPUTING RESOURCES
# Should be fine with M1 Apple or similar 
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#* Function  `run_parallel_model` implements a distributed version of the DES
#* k times for a given sim sample size,  calls function `obtain_model_grid_trace` to
#* compute the trace matrix approximation in discrete time k times. 
#* Trace Matrices  are stored in an array `a_grid_trace_byState`.

#* Function `batch_summary` summarizes the k trace matrices and returns different quantiles, mean and SE of 
#* state occupancy for a give simulation sample size. 

# dir.create("analysis", showWarnings = FALSE) # Create directory if it doesn't exist
# # ------------------ #
# ## run for n_sim 1e2
l_inputs                     <- init_params(n_sim = v_n_sim[1])
a_grid_k_5e1_b_size_1e2      <- run_parallel_model(model_fn = obtain_model_grid_trace, l_inputs = l_inputs, k=k)
#Summary
batch_summary_k_5e1_b_1e2    <- batch_summary(array_grid = a_grid_k_5e1_b_size_1e2, k = k, n_sim = v_n_sim[1])
rm(a_grid_k_5e1_b_size_1e2)
# # ------------------ #
# ## run for  n_sim 1e3
l_inputs                     <- init_params(n_sim = v_n_sim[2])
a_grid_k_5e1_b_size_1e3      <- run_parallel_model(model_fn = obtain_model_grid_trace, l_inputs = l_inputs, k=k)
# Summary
batch_summary_k_5e1_b_1e3    <- batch_summary(array_grid = a_grid_k_5e1_b_size_1e3, k = k, n_sim = v_n_sim[2])
rm(a_grid_k_5e1_b_size_1e3)
# # ------------------ #
# ## run for n_sim 1e5
l_inputs <- init_params(n_sim = v_n_sim[4])
a_grid_k_5e1_b_size_1e5 <- run_parallel_model(model_fn = obtain_model_grid_trace, l_inputs = l_inputs, k=k)
# Summary
batch_summary_k_5e1_b_1e5    <- batch_summary(array_grid = a_grid_k_5e1_b_size_1e5, k = k, n_sim = v_n_sim[4])
rm(a_grid_k_5e1_b_size_1e5)

# STOP CLUSTER
stopCluster(cl)
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# Summaries ----

batch_summary_all <- as.data.frame(rbind(batch_summary_k_5e1_b_1e2,
                                         batch_summary_k_5e1_b_1e3,
                                         batch_summary_k_5e1_b_1e5)) %>% 
  mutate(
    Sim_n = ifelse(n_sim == 100, "Sim. sample size = 1e2", 
                   ifelse(n_sim == 1000, "Sim. sample size = 1e3",
                          "Sim. sample size = 1e5"))
    
  ) %>% mutate( Sim_n = factor(Sim_n, levels = c("Sim. sample size = 1e2","Sim. sample size = 1e3","Sim. sample size = 1e5") ),
                State = factor(State, levels = c("Healthy", "Sick","Sicker", "Dead")))

table(batch_summary_all$n_sim)
# Plot ----

plot_convergence1<- ggplot(batch_summary_all,aes( x= T_n, y = Mean, color = State, fill =State))+
  geom_line()+
  geom_ribbon( aes(ymin = Q025, ymax = Q975), alpha = .2, linewidth = 0.2 )+
  # Titles and Labs
  xlab("Time (years)") +
  ylab("State Occupancy (Proportion)")+
  labs(title    = "Median State Occupancy over Time by Simulation Sample Size",
       subtitle = "Confidence Bounds(2.5%,97.5%), Batch Size = 5e1")+
  
  #Grid 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), 
                     minor_breaks = seq(0, 1, by = 0.05))+
  
  # Colors
  scale_color_viridis_d(option = "H",begin = 0, end = 1,direction = 1, guide = "none") +
  scale_fill_viridis_d(option =  "H",begin = 0, end = 1,direction = 1, name = "Health States") +
  
  # Theme 
  theme_minimal(base_size = 28)+
  theme(
    strip.background = element_rect(fill = "lightgray", color = "NA", linewidth = 0.5),
    strip.text = element_text(face = "bold"),
    panel.background = element_rect(fill = "white", color = "gray90"),
    panel.grid = element_line(color = "gray90"),
    
    axis.text.x = element_text(size = 16, angle = 0, color = "black"),
    axis.text.y = element_text(size = 16, angle = 0, color = "black"),
    
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    
    plot.title    = element_text( face = "bold"),
    plot.subtitle = element_text(face = "italic"), 
    
    legend.position = "bottom"
    
  )+
  facet_wrap(~Sim_n)
plot_convergence1

require(ggview)
main_plot<- plot_convergence1 + canvas( width  = 28,
                                height = 20,
                                units  = c("in"),
                                dpi    = 200)
main_plot
