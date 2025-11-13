#------------------------------------------------------------------------------#
# Post-processing Module C: Probabilistic Analysis 

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

#* This code implements a Probabilistic  Analysis of the
#* CEA of the Sick-Sicker model under four strategies. 
#* The simulation is implemented as a Discrete Event Simulation (DES)
#* The four strategies are the following  
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
#* Probabilistic Analysis-based Economic Outcomes
#* - Cost-effectiveness scatter plot
#* - Cost-effectiveness acceptability curves (CEACs) and frontier (CEAF)
#* - Expected Loss Curves (ELCs)
#* - Expected value of perfect information (EVPI)
#------------------------------------------------------------------------------#
# Initial setup ---- 
rm(list = ls())    # remove any variables in R's memory 
gc()               # clean working memory
#------------------------------------------------------------------------------#
## Load packages ----
library("data.table"  )
library("dplyr"       )
library("ggplot2"     )
library("ellipse"     )
library("patchwork"   )

library("doParallel"  )
library("parallel"    )
library("foreach"     )

library("MethylCapSig")

## Load supplementary functions ----
source("R/DES_functions.R") # DES specific 
source("R/Functions.R")     # General utility functions
# PSA functions from previous tutorials 
# (with permission from Dr. Alarid-Escudero)
# source("R/Functions_cSTM_time_dep_state_residence.R") 
#------------------------------------------------------------------------------#

# Initialize Sim. Model Parameters ----
#* We want to pass parameters seamlessly through our function sim_next_event(), 
#* for this we can declare initial parameters in a list. Here we pass them through
#* a function `init_param` that will return such list `l_params` and also include 
#* the baseline data table `dt_baseline` as an object of that list, now named `dt_next_event`. 

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
  v_states         = c("H",                                     # Healthy (H)
                       "S1",                                    # Sick (S1)
                       "S2",                                    # Sicker (S2)
                       "D"),                                    # Dead (D)
  
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
#  Reference l_params
l_params_ref          <- l_params
# Store the parameter names into a vector
v_names_params <- names(l_params_ref)

# ------------------------------------------------------------------------ #
# Initialize Economic Parameters ----
l_cea_params<- list(
  # # list of event histories for each strategy
  l_event_history_strategies = NA,
  
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

#  Reference l_params
l_cea_params_ref          <- l_cea_params

# Store the parameter names into a vector
v_names_cea_params <- names(l_cea_params_ref)

#------------------------------------------------------------------------------#
#  Probabilistic Analysis input Dataset ----
n_sim_psa   <- 5#1e3

v_names_str <- c("Standard of care",      # store the strategy names
                 "Strategy A", 
                 "Strategy B",
                 "Strategy AB") 
n_str       <- length(v_names_str)        # number of strategies

# Generate PSA Input
dt_psa_input <- as.data.table(generate_psa_params_DES(n_sim = n_sim_psa))
# First six observations
# head(dt_psa_input)

# ### Histogram of parameters ----
# ggplot(melt(dt_psa_input, variable.name = "Parameter"), aes(x = value)) +
#   facet_wrap(~Parameter, scales = "free") +
#   geom_histogram(aes(y = after_stat(density))) +
#   ylab("") +
#   theme_bw(base_size = 16) +
#   theme(axis.text = element_text(size = 6),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.y  = element_blank(),
#         axis.ticks.y = element_blank())
#------------------------------------------------------------------------------#
#------------------- #
# Run Probabilistic Analysis  -----
#------------------- #
## Initialize data.frames with Probabilistic Analysis  output  ----
### data.frame of costs 
df_c <- as.data.frame(
  matrix(0, 
         nrow = n_sim_psa,
         ncol = n_str))
colnames(df_c) <- v_names_str

### data.frame of effectiveness 
df_e <- as.data.frame(
  matrix(0,
         nrow = n_sim_psa,
         ncol = n_str))
colnames(df_e) <- v_names_str

###  data.frame of NMB  
df_nmb <- as.data.frame(
  matrix(0,
         nrow = n_sim_psa,
         ncol = n_str))
colnames(df_nmb) <- v_names_str

v_strategies_names <- c( "SoC", "A", "B", "AB" )
#* since the effect on transition dynamics from SoC == A, and B == AB, 
#* we can just execute the simulation for "SOC" and "B"
v_strategies       <- c("SoC","B")        
# Set seed for reproducibility 
n_seed             <- 2
# ------------------------------------------------------------------------ #
# ------------------------------------------------------------------------ #

# WARNING: DO NOT RUN THE NEXT CHUNK UNLESS YOU HAVE ENOUGH COMPUTING RESOURCES

# ------------------------------------------------------------------------ #
# ------------------------------------------------------------------------ #
strategies <- v_strategies
## Non-distributed Implementation ----
# for (i in 1:n_sim_psa) {
#   print(paste("Init_PSA_iteration",i, sep = "_"))
#   l_params                  <- update_param_list(l_params_ref, dt_psa_input[i,])
#   l_cea_params              <- update_param_list(l_cea_params_ref, dt_psa_input[i,])
#   # A. Full Histories
#   l_event_history_strategies        <- list()
#   # Loop over strategies
#   for(j in seq_along(strategies)){
#     print(paste("Strategy",strategies[j],"iteration",i, sep = "_"))
#     l_params$strategy           <- strategies[j]
#     dt_event_history            <- data.table()
#     l_params$dt_next_event      <- l_params_ref$dt_next_event
#     # common random numbers for strategies
#     set.seed(n_seed)
#     # run simulation under i^th strategy
#     while(nrow(l_params$dt_next_event) >0){
#       l_params$dt_next_event          <- sim_next_event(l_params)
#     }
#     l_event_history_strategies[[j]]   <- as.data.table(dt_event_history)
#   }
#   # Store event-history for all strategies
#   l_event_history_strategies$A        <- l_event_history_strategies[[1]]
#   l_event_history_strategies$AB       <- l_event_history_strategies[[2]]
#   names(l_event_history_strategies)   <- c(strategies, "A","AB")
# 
#   # Economics Outcomes
#   l_cea_params$l_event_history_strategies <- l_event_history_strategies
#   l_cea_PSA_iter                  <- cea_psa_fn(l_cea_params = l_cea_params)
#   # PSA based outcomes
#   df_c[i, ]                       <- l_cea_PSA_iter$Cost
#   df_e[i, ]                       <- l_cea_PSA_iter$Effect
#   df_nmb[i,]                      <- l_cea_PSA_iter$NMB
#   print(paste("Fin_PSA_iteration",i, sep = "_"))
#   # # Display simulation progress
#   if (i/(n_sim_psa/100) == round(i/(n_sim_psa/100), 0)) { # display progress every 5%
#     cat('\r', paste(i/n_sim_psa * 100, "% done", sep = " "))
#   }
# }

# # ------------------------------------------------------------------------ #
# # ------------------------------------------------------------------------ #
# ## Distributed Implementation  ----
# Num. of Cores
no_cores_all  <- parallel::detectCores()-1      # recommended when not sharing computing resources
# cores_90pct   <- as.integer(.90 * no_cores_all) # aggressive  when shared resources
cores_85pct   <- as.integer(.85 * no_cores_all) # aggressive  when shared resources
# cores_75pct   <- as.integer(.75 * no_cores_all) # aggressive  when shared resources
# cores_50pct   <- as.integer(.5  * no_cores_all) # recommended when shared resources
# cores_25pct   <- as.integer(.25 * no_cores_all) # recommended when shared resources


#### Cluster  ----
os_info  <- Sys.info()
os       <- os_info[1]
no_cores <- cores_85pct

if (os == "Darwin") { # Apple Kernel
  # Initialize cluster object
  cl <- parallel::makeForkCluster(no_cores)
  # Register clusters
  doParallel::registerDoParallel(cl)
}
if (os == "Linux") {
  # Initialize cluster object
  cl <- parallel::makeCluster(no_cores)
  # Register clusters
  doParallel::registerDoMC(cl)
}
if (os == "Windows") {  
  # Initialize cluster object
  cl <- parallel::makeCluster(no_cores)
  # Register clusters
  doParallel::registerDoParallel(cl)
}

#### Run distributed Probabilistic Analysis  ----
df_ce <- foreach::foreach(i = 1:n_sim_psa, .combine = rbind,
                          .packages = c("data.table","dplyr")
) %dopar% {
  #param lists
  l_params                  <- update_param_list(l_params_ref, dt_psa_input[i,])
  l_cea_params              <- update_param_list(l_cea_params_ref, dt_psa_input[i,])
  # A. Simulate Full Event Histories
  l_event_history_strategies              <- list()
  ## Loop over strategies
  for(j in seq_along(strategies)){
    #print(paste("Strategy",strategies[j],"iteration",j, sep = "_"))
    l_params$strategy             <- strategies[j]
    dt_event_history              <- data.table()
    l_params$dt_next_event        <- l_params_ref$dt_next_event
    ### common random numbers for strategies
    set.seed(n_seed)
    ### run simulation under i^th strategy
    while(nrow(l_params$dt_next_event) >0){
      l_params$dt_next_event              <- sim_next_event(l_params)
    }
    l_event_history_strategies[[j]]       <- as.data.table(dt_event_history)
  }
  ## Store event-history for all strategies
  l_event_history_strategies$A            <- l_event_history_strategies[[1]]
  l_event_history_strategies$AB           <- l_event_history_strategies[[2]]
  names(l_event_history_strategies)       <- c(strategies, "A","AB")
  
  # B. Compute Economic Outcomes
  l_cea_params$l_event_history_strategies <- l_event_history_strategies
  l_cea_PSA_iter                          <- cea_psa_fn(l_cea_params = l_cea_params)
  df_ce                                   <- c(l_cea_PSA_iter$Cost, l_cea_PSA_iter$Effect)
}
# Extract costs and effects from the PSA dataset
df_c <- df_ce[, 1:n_str]                # columns with costs
df_e <- df_ce[, (n_str + 1):(2*n_str)]  # columns with effects
# Stop clusters
parallel::stopCluster(cl)
# ------------------------------------------------------------------------ #
# Probabilistic Analysis  Data Object ----
# #* Function included in "R/Functions.R" The latest version can be found in `dampack` package
l_psa <- make_psa_obj(cost          = df_c,
                      effectiveness = df_e,
                      parameters    = dt_psa_input,
                      strategies    = v_names_str)
l_psa$strategies              <- v_names_str
colnames(l_psa$effectiveness) <- v_names_str
colnames(l_psa$cost)          <- v_names_str
#View(l_psa)

# ------------------------------------------------------------------------ #
# Plots  -----
# ------------------------------------------------------------------------ #
# Vector with willingness-to-pay (WTP) thresholds.
v_wtp <- seq(0, 200000, by = 1000)
txtsize <- 13
#------------------- #
## Cost-Effectiveness Scatter plot ----
#------------------- #

#* Function `plot.psa` included in "R/Functions.R"; depends on `tidyr` and `ellipse` packages.
#* The latest version can be found in `dampack` package

gg_scattter <- plot.psa(l_psa, txtsize = txtsize) +
  scale_color_viridis_d( option = "magma",begin = 0.1,end = 0.8, direction = 1)+
  scale_fill_viridis_d( option = "magma",begin = 0.1,end = 0.8, direction = 1)+
  theme_bw(base_size = 18)+
  scale_y_continuous("Cost (Thousand $)",
                     breaks = number_ticks(10),
                     labels = function(x) x/1000) +
  xlab("Effectiveness (QALYs)") +
  guides(col = guide_legend(nrow = 1)) +
  theme(legend.position = "none") 
gg_scattter

#------------------- #
## Incremental cost-effectiveness ratios (ICERs) with probabilistic output ----
#------------------- #

#* Compute expected costs and effects for each strategy from the PSA
#* Function included in "R/Functions.R". The latest version can be found in `dampack` package
df_out_ce_psa <- summary(l_psa)

#* Function `calculate_icers` included in "R/Functions.R"; depends on the `dplyr` package
#* The latest version can be found in `dampack` package
df_cea_psa <-  calculate_icers(cost       = df_out_ce_psa$meanCost,
                               effect     = df_out_ce_psa$meanEffect,
                               strategies = df_out_ce_psa$Strategy)

#------------------- #
## CEA Frontier with probabilistic output ----
#------------------- #

#* Function `plot.icers` included in "R/Functions.R"; depends on the `ggplot2`  and `ggrepel` packages.
#* The latest version can be found in `dampack` package

p_icers<- plot.icers(df_cea_psa, label = "all", txtsize = txtsize) +
  #expand_limits(x = max(table_cea$QALYs) + 0.1) +
  theme(legend.position = "none")
p_icers
#------------------- #
## Cost-effectiveness acceptability curves (CEACs) and frontier (CEAF) ----
#------------------- #

#* Functions `ceac` and  `plot.ceac` included in "R/Functions.R". The latest versions can be found in `dampack` package
ceac_obj <- ceac(wtp = v_wtp, psa = l_psa)
#* Regions of highest probability of cost-effectiveness for each strategy
summary(ceac_obj)
#* CEAC & CEAF plot
gg_ceac <- plot.ceac(ceac_obj, txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 14) +
  scale_color_viridis_d( option = "magma",begin = 0.1,end = 0.8, direction = 1)+
  scale_fill_viridis_d( option = "magma",begin = 0.1,end = 0.8, direction = 1) +
  theme(legend.position = "right") + 
  theme_bw(base_size = 18)
gg_ceac

#------------------- #
## Expected Loss Curves (ELCs) ----
#------------------- #

#* Functions `calc_exp_loss` and `plot.exp_loss` included in "R/Functions.R".The latest version can be found in `dampack` package
elc_obj <- calc_exp_loss(wtp = v_wtp, psa = l_psa)
#* ELC plot
gg_elc <- plot.exp_loss(elc_obj, log_y = FALSE,
                        txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 14,
                        col = "full") +
  scale_color_viridis_d( option = "magma",begin = 0.1,end = 0.8, direction = 1)+
  scale_fill_viridis_d( option = "magma",begin = 0.1,end = 0.8, direction = 1)+
  # geom_point(aes(shape = as.name("Strategy"))) +
  scale_y_continuous("Expected Loss (Thousand $)",
                     breaks = number_ticks(10),
                     labels = function(x) x/1000) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none") 

gg_elc

#------------------- #
## Expected value of perfect information (EVPI) ----
#------------------- #

#* Functions `calc_evpi` and `plot.evpi` included in "R/Functions.R". The latest version can be found in `dampack` package
evpi <- calc_evpi(wtp = v_wtp, psa = l_psa)
#* EVPI plot
gg_evpi <- plot.evpi(evpi, effect_units = "QALY",
                     txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 14) +
  scale_y_continuous("EVPI (Thousand $)",
                     breaks = number_ticks(10),
                     labels = function(x) x/1000) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none") 
gg_evpi

#------------------- #
## Combined Plot Probabilistic Analysis  ----
#------------------- #

PA_Economic_Outcomes_plot <- (gg_scattter + gg_ceac) / (gg_elc + gg_evpi)  + 
  plot_layout(guides = "collect")+
  plot_annotation(tag_levels = 'A',
                  title = "Probabilistic Analysis-based Economic Outcomes",
                  subtitle = "")& theme(
                    plot.title = element_text(size = 22, face = "bold", hjust = 0.5)
                  ) 
PA_Economic_Outcomes_plot
# ----------------------------------- #
# ----------------------------------- #

