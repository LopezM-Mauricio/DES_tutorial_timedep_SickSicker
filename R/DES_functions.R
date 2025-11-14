# Functions to implement Discrete Event Simulations 
# of individual level state transition models. 
#* Authors: 
#* - XXXXX
#* - XXXXX
#* - XXXXX
#------------------------------------------------------------------------------#
#* Please cite the corresponding tutorials when using these functions: 
#* XXXXX open article
#* doi: XXXXX
#* bibtex: XXXXX

#* -
#* -
#* -
#* -
#* -
#----------------------------------------------------------------------------#
#----------------------------------------------------------------------------#
#----------------------------------------------------------------------------#
# Negation of % in % operator
`%nin%` = Negate(`%in%`)
#----------------------------------------------------------------------------#
#----------------------------------------------------------------------------#
#                   Nonparametric Sampling                                   #
#        non-homogeneous poisson point process                               #
#                          (NPS-NHPPP) Functions                          ----
#----------------------------------------------------------------------------#
#* When using these set of functions please cite: 
#* XXXXXX
#----------------------------------------------------------------------------#
#---------------------------------------------------- #
##  Age, Sex specific Probabilities of Event Times ----     
#---------------------------------------------------- #

#' Obtain probabilities of events by Age, and Sex from lifetable rates
#'
#' \code{obtain_rates_des} obtain probabilities of events by Age, and Sex from lifetable rates
#' 
#' The elements of each row of the data.frame should sum up to 1. 
#'
#' @param df_mortality_rate_long data.frame with mortality rates by Age and Sex
#' @return A data.frame with probabilities of events by Age and Sex
#' @export
obtain_probs_des <- function(df_mortality_rate_long){
  max_age <- max(df_mortality_rate_long$Age)
  
  df_mortrate_expanded <- data.frame()
  
  for (sex in unique(df_mortality_rate_long$Sex)) {
    # Load and filter base dataset 
    df_mortrate <- df_mortality_rate_long %>%
      # Sex has to be changed in each iteration
      filter(Sex == sex) %>%
      group_by(Year) %>% 
      mutate(death_prob_cum = 1 - exp(-cumsum(death_rate)),
             death_prob = death_prob_cum - lag(death_prob_cum)
             # death_prob = c(death_prob_cum[1], diff(death_prob_cum))
             # death_prob = diff(c(0, death_prob_cum))
      ) %>%
      ungroup() %>%
      mutate(death_prob = if_else(condition = Age == 0,
                                  true = death_prob_cum,
                                  false = death_prob))
    
    # Create a dataframe with the transposed death probabilities
    df_DP <- matrix(data = df_mortrate$death_prob,
                    ncol = length(unique(df_mortrate$Age)),
                    byrow = TRUE) %>%
      as.data.frame() %>%
      # Repeat eachrow by the total number of ages (110)
      slice(rep(1:n(),
                each = length(unique(df_mortrate$Age))))
    
    # Rename columns
    colnames(df_DP) <- paste0("DP", 0:max_age)
    
    # Create a upper triangular matrix, then df, for the 110 ages
    df_tri <- (upper.tri(x = matrix(data = 1,
                                    nrow = length(unique(df_mortrate$Age)),
                                    ncol = length(unique(df_mortrate$Age))
    ),
    diag = TRUE)*1) %>%
      as.data.frame()
    
    colnames(df_tri) <- paste0("V", 0:max_age)
    
    # stack arrays vertically
    df_tri_exp <- do.call("rbind",
                          replicate(n = (nrow(df_DP)/length(unique(df_mortrate$Age))),
                                    expr = df_tri,
                                    simplify = FALSE))
    
    # Check that all dataframes have the correct dimensions
    dim(df_mortrate)
    dim(df_DP)
    dim(df_tri_exp)
    
    name_last_call_DP  <- colnames(df_DP)[max_age + 1]
    name_last_call_tri <- colnames(df_tri)[max_age + 1]
    # Data Manipulation 
    df_mortrate_onesex <- df_mortrate %>%
      bind_cols(df_DP, df_tri_exp) %>%
      # Make multiplication to obtain
      mutate(across(DP0:eval(name_last_call_DP)) * across(V0:eval(name_last_call_tri))) %>%
      # Remove the upper triangle matrix
      select(-c(V0:eval(name_last_call_tri))) %>%
      # Obtain proportion of probability in relation to the rowsum
      mutate(across(.cols = DP0:eval(name_last_call_DP)) / rowSums(across(.cols = DP0:eval(name_last_call_DP)))) %>%
      # Replace the NA's by 0
      replace(is.na(.), 0)
    df_mortrate_expanded <- bind_rows(df_mortrate_expanded, df_mortrate_onesex)
  }
  
  return(df_mortrate_expanded)
}

#---------------------------------------------------- #
##  Non-parametric Sampling of Event Times        -----
#---------------------------------------------------- #
#' Nonparametric sampling of time to events from a discrete nonhomogeneous 
#' Poisson point process
#'
#' \code{nps_nhppp} samples states for multiple individuals simultaneously.
#' 
#' The elements of each row of the matrix `m_probs` should sum up to 1. In case
#' this does not happens, all the rows where this condition is not met will be
#' normalized.
#'
#' @param m_probs matrix with probabilities for each category. Each row 
#' represents one individual and each column an specific category. 
#' @param v_categories An optional argument. It is a vector containing the name of 
#' the categories to sample from. The length of this vector must be equal to 
#' the number of columns of `m_probs`.
#' @return A vector filled with sampled categories.
#' @export
nps_nhppp <- function(m_probs, 
                      v_categories = NULL, 
                      correction = c("none", "uniform")) {
  
  if (!is.numeric(m_probs)) {
    stop("`m_probs` must be a matrix filled with numeric elements")
  }
  
  if (!isTRUE(all.equal(target = rep(1, nrow(m_probs)), 
                        current = as.numeric(rowSums(m_probs))))) {
    
    #* Find rows where the sum of their elements is not equal to 1.
    #* We use a very small value to check that they are equal to 1 since
    #* == and != sometimes do not work as desired. The value 1.5e-8
    #* was taken as the same level of tolerance as the `all.equal` function.
    v_rows_to_norm <- which(abs(1 - rowSums(m_probs)) > 1.5e-8)
    
    warning(
      "The rows: ", 
      paste(as.character(v_rows_to_norm), collapse = ", "), 
      " in `m_probs` do not sum up to 1. The values within these rows will be normalized and then used in the sampling process. In case this behaviour is not desired modify the content of `m_probs`.")
    
    # Normalize rows
    if (length(v_rows_to_norm) == 1) {
      m_probs[v_rows_to_norm, ] <- m_probs[v_rows_to_norm, ]/sum(m_probs[v_rows_to_norm, ])
    }
    if (length(v_rows_to_norm) > 1) {
      m_probs[v_rows_to_norm, ] <- m_probs[v_rows_to_norm, ]/rowSums(m_probs[v_rows_to_norm, ])
    }
  }
  
  # Number of categories to sample from
  n_cat <- ncol(m_probs)
  
  if (is.null(v_categories)) {
    #' Generate the numeric categories based on the number of columns of 
    #' `m_probs`
    v_categories <- seq(0, (n_cat - 1))
  }
  
  correction <- match.arg(correction)
  
  # Number of time to events to draw
  n_samp <- nrow(m_probs)
  
  v_unif  <- runif(n_samp, min = 0, max = 1)
  # matrixStats is written in C and optimized for large Matrix. 
  v_sum_p <- matrixStats::rowCumsums(m_probs)
  # slower alternative using base R
  # v_sum_p <- t(apply(m_probs))
  
  v_time_to_event <- v_categories[max.col(v_sum_p >= v_unif, 
                                          ties.method = "first")]
  
  if (correction == "uniform") {
    v_time_to_event <- v_time_to_event + v_unif
  }
  return(v_time_to_event)
}


#----------------------------------------------------------------------------#
#----------------------------------------------------------------------------#
# DES Functions -----
#----------------------------------------------------------------------------#
#* Functions developed specifically for the DES tutorial
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

#---------------------------------------------------- #
##  Initial Parameters    ----   
#---------------------------------------------------- #

init_params <- function(
    # Sim
  time_init        = 25    ,
  time_end         = 100   ,
  n_sim            = 1e3   ,
  # Init DT        
  age_init         = 25    ,
  state_init       = "H"   , 
  # Sex-group Mortality
  Sex_grp          = "Total",  # Other alternatives: "Male", "Female", "Combined"
  # Model      
  v_states         = c("H",                                              # Healthy (H)
                       "S1",                                             # Sick (S1)
                       "S2",                                             # Sicker (S2)
                       "D"),                                             # Dead (D)
  
  m_Tr_r           = matrix(c(0,1,0,2,
                              3,0,4,5,
                              0,0,0,6,
                              0,0,0,0), 
                            nrow = 4, ncol = 4, byrow = TRUE,
                            dimnames = list(c("H", "S1","S2","D"),
                                            c("H", "S1","S2","D"))),
  # Transition Rates 
  ## Mortality
  lifetable        = "Inputs/all_cause_mortality.Rda",
  ### Hazard Ratio while being S1
  hr_S1            =    3,
  ### Hazard Ratio while being S2
  hr_S2            =    10, 
  
  ## H  -> S1, constant rate ~Exp(r_HS1)
  r_HS1            =    0.15,
  ## S1 -> H,  constant rate ~Exp(r_S1H)
  r_S1H            =    0.5,
  ## S1 -> S2, state-residency time dependent ~ Weibull(r_S1S2_scale_ph, r_S1S2_shape)
  ### PH parameterization
  r_S1S2_scale_ph  =    0.08,
  r_S1S2_shape     =    1.1,
  
  # Strategy 
  strategy         = "SoC",   # Can be any of "SoC", "A", "B", "AB". Only B and AB, affect transitions rates, 
  #Tx Effects that affect state-transitions
  ## Effectiveness of strategy B 
  hr_S1S2_trtB = 0.6  # hazard ratio of becoming Sicker when Sick under treatment B
){
  
  # Time horizon 
  time_horizon       <- c(time_init,time_end)
  
  # Indexed transitions
  dt_trans_keys <- data.table()
  dt_Tr_r       <-  as.data.table(m_Tr_r)
  dt_Tr_r[, from := v_states]
  dt_trans_keys <- as.data.table(melt(dt_Tr_r, id.vars = "from",
                                      measure.vars = c("H" ,"S1" ,"S2" ,"D"),
                                      variable.name = "to",
                                      value.name = "transition"))
  dt_trans_keys <- dt_trans_keys[transition !=0,][order(transition)]
  
  # load lifetable, all cause mortality
  load(lifetable)
  dt_bckgrd_mortality   <- as.data.table(all_cause_mortality) %>% mutate(death_rate  = Rate) %>% dplyr::filter(Year== 2015)
  
  # Init Data Table 
  dt_next_event               <- data.table(
    # Covariate info
    ID                  = 1:n_sim,#sample(0:n_samp*5, n_sim, replace = FALSE),                     # n_samp unique random numerical identifiers
    Age                 = age_init,#sample(25:99, n_sim, replace = TRUE),                          # for comparison we do homogenous age 25 years # Uniform sampling values, updated mechanistically
    
    
    Sex                 = Sex_grp,
    Sex_comb            = sample(c("Female","Male"), n_sim, prob = c(0.55,0.45), replace = TRUE),  # Uniform sampling values male/female,  constant covariate
    Race                = sample(c("White", "Black", "Other"), n_sim, replace = TRUE),             # Uniform sampling from "White", "Black", "Other",  constant covariate
    Smoke_6months_t0    = sample(0:6, n_sim, replace = TRUE),                                      # Uniform sampling from 0 to 6,  constant covariate
    # Transition Dynamics          
    from             = state_init,                                                                 # Origin State, first row, everyone starts in H, update accordingly at each event. 
    to                  = NA,                                                                      # Destination State, first row, everyone is NA, update accordingly at each event. 
    T_n                 = time_horizon[1],                                                                       # Simulation Time, records time of events, first chunk, all entries are zero, but will be updated to the time when transition/event  occurs 
    T_start             = 0,                                                                        # State-residence Time start, first chunk, all entries are zero, but will be updated to the transition that occurs next
    T_stop              = 0,                                                                        # State-residence Time stop, first chunk, all entries are zero, but will be updated to the transition that occurs next
    tau                 = 0,                                                                        # Time-increment, first chunk, all entries are zero, updated after each event
    status              = NA,                                                                       # Occurrence/Censoring status, transition-specific indicator, update row by row. 
    Event_num           = 1                                                                         # Event number, at least 1 event from H-> D over time-horizon, updated after each event
  )
  dt_next_event[,T_start := T_n]
  
  # Sex-specific Mortality 
  if(Sex_grp == "Combined"){
    dt_next_event[, Sex := Sex_comb][,Sex_comb := NULL]
  } else{dt_next_event[,Sex_comb := NULL]}
  if(Sex_grp %nin% c("Total","Male","Female","Combined")){
    warning("Error: Note you have specified a group that is not considered in the Life Table. 
            Allowable groups are: Total, Male, Female, Combined")
  }
  
  # Wrap list for sim
  l_params    <- list(
    time_horizon        = time_horizon,
    n_sim               = n_sim,
    #Sex_grp             = Sex_grp,
    
    age_init            = age_init, 
    state_init          = state_init,
    v_states            = v_states,
    m_Tr_r              = m_Tr_r, 
    dt_trans_keys       = dt_trans_keys,
    
    dt_next_event       = dt_next_event,
    
    lifetable           = lifetable, 
    dt_bckgrd_mortality = dt_bckgrd_mortality,
    hr_S1               = hr_S1,
    hr_S2               = hr_S2, 
    
    r_HS1               = r_HS1,
    r_S1H               = r_S1H,
    
    r_S1S2_scale_ph     = r_S1S2_scale_ph,
    r_S1S2_shape        = r_S1S2_shape,
    
    strategy            = strategy,
    
    hr_S1S2_trtB        = hr_S1S2_trtB)
  
  return(l_params)
}
#---------------------------------------------------- #
##  Update Next Event    ----   
#---------------------------------------------------- #

sim_next_event <- function(l_params){
  with(as.list(l_params), {
    # Current Event Setup ------------------------------ #
    # Determine set of current states
    v_current_state          <- unique(dt_next_event$from) 
    v_current_state_names    <- paste("from", v_current_state, sep = "_")
    
    # Subset dt_next_event based on origin state 
    l_subsets_dt_next_event <- list()
    for( i in seq_along(v_current_state)){
      l_subsets_dt_next_event[[i]] <-  as.data.table(dt_next_event[from == v_current_state[i],])
    }
    # Add names to subsets
    names(l_subsets_dt_next_event) <- v_current_state_names
    
    # Mortality lifetables
    # mortality from S1
    dt_mortality_Sick        <- copy(dt_bckgrd_mortality)
    #death_rate_S1            <- dt_bckgrd_mortality$death_rate*hr_S1
    dt_mortality_Sick[,death_rate := death_rate*hr_S1]
    # mortality from S2
    dt_mortality_Sicker      <- copy(dt_bckgrd_mortality)
    #death_rate_S2            <- dt_bckgrd_mortality$death_rate*hr_S2
    dt_mortality_Sicker[,death_rate := death_rate*hr_S2]
    
    # ----------------------------------------------------------------------------- #
    # Module 1: transitions from Healthy             ------------------------------ #
    # ----------------------------------------------------------------------------- #
    # Module Set up
    # ----------------------------------- #
    # Task:1 Possible Transitions from Healthy
    # ----------------------------------- #
    
    temp   <- l_subsets_dt_next_event$from_H
    
    # Expand each subset to have one row for each possible transition for each ID
    if(length(temp)>0){ #necessary because at some point there may not be healthy people
      
      current_state          <- unique(temp$from)
      transition             <- as.numeric(m_Tr_r[current_state,m_Tr_r[current_state,] !=0])            # number of possible transitions
      ID                     <- unique(temp$ID)                                                         # unique IDs
      setkey(temp, ID)      
      temp_long              <- merge(temp,CJ(ID,transition), all.x = F)                                # CJ is like expand.grid or expand_grid but for data.tables, can be used in more flexible ways, also more efficient
      temp_long[,c("from","to") := NULL] # clean join afterwards
      
      # Assign "from" and to" state based on possible transitions, indexed join in data.table
      temp_long              <- temp_long[dt_trans_keys, on = "transition", nomatch = 0L]               # left join trick on data.table, "on" argument sets the Key automatically
      temp_long              <- temp_long[order(ID, transition)]
      
    # ----------------------------------- #
    # Task 2. Sample latent arrival times starting from Healthy 
    # ----------------------------------- #
      # ----------------------------------- #
      ## Submodule  1.1: Transitions Healthy -> Sick   
      # ----------------------------------- #
      # constant annual rate of becoming Sick when Healthy
      # Sample from exp. distr.
      if(r_HS1 != 0){
        temp_long[transition == 1 , T_stop := rexp(.N, rate = r_HS1) + T_n]
      } else {
        temp_long[transition == 1 , T_stop := Inf]}
      # ----------------------------------- #
      ## Submodule  1.2: Transitions Healthy -> Dead
      # ----------------------------------- #
      # age-dependent, Sex stratified background mortality
      dt_time2death_probs            <- as.data.table(obtain_probs_des(dt_bckgrd_mortality))                         # obtain vector of probabilities
      
      # Subset those rows for which the transition is possible                
      dt_trans2                      <- temp_long[transition == 2]               
      
      # Set key on the identifier/strata columns (Age,Sex) for all data tables               
      setkey(dt_trans2, Age,Sex)                                                                                     # baseline reference dt
      setkey(dt_time2death_probs, Age,Sex)                                                                           # mortality probs from Health state
      # indexed merge, adding columns with interval probabilities
      dt_2sample_time2death          <- merge(dt_trans2[,Age := ceiling(Age)], dt_time2death_probs, all.x = TRUE)    # all = TRUE performs a full outer join, all.x performs left-join, all.y right-join
      
      # matrix of probabilities
      m_probs                        <- as.matrix(dt_2sample_time2death[, .SD, .SDcols = patterns("^DP")])           # choose only the probs cols
      # Sample time to death
      dt_trans2[,T_stop_aux := nps_nhppp(m_probs = m_probs, correction = "uniform")]
      # ---------------- #
      # Additional correction addressing a minor approximation error with the nps_nhppp method. (ask MLM about it)
      dt_trans2[T_n > Age, tau_aux := T_stop_aux - Age]
      dt_trans2[T_n > Age, T_stop  := T_n + tau_aux]
      dt_trans2[T_n <= Age, T_stop := T_stop_aux]
      #dt_trans2[T_start >T_stop,]
      # ---------------- #
      setkey(temp_long,ID,T_n)
      setkey(dt_trans2,ID,T_n)
      temp_long[transition == 2, T_stop := dt_trans2$T_stop]
      
      # ----------------------------------- #
      # Task 3: Predict the next state 
      # ----------------------------------- #
      # Create an indicator for minimum value in T_stop by ID and calculate tau
      temp_long[, status := as.numeric(T_stop == min(T_stop)), by = ID]
      temp_long[,tau := T_stop - T_start]
      # ----------------------------------- #
      # Task 4: Update important variables 
      # ----------------------------------- #
      
      ### Update Event Time
      temp_long[status == T,T_n := T_stop]
      ### Update Age
      temp_long[status == T,Age := round(T_n)]
      # ----------------------------------- #
      # Update dt_next_event_long 
      dt_event_history <<- rbind(dt_event_history,temp_long)
      setkey(dt_event_history,ID, T_n, Event_num, transition)
    }
    # ----------------------------------------------------------------------------- #
    # Module 1: transitions from Sick                ------------------------------ #
    # ----------------------------------------------------------------------------- #
    # Module Set up
    # ----------------------------------- #
    # Task:1 Possible Transitions from Sick
    # ----------------------------------- #
    
    temp <- l_subsets_dt_next_event$from_S1
    
    if(length(temp)>0){ #necessary because at some point there may not be S1 people
      current_state          <- unique(temp$from)
      transition             <- as.numeric(m_Tr_r[current_state,m_Tr_r[current_state,] !=0])         # possible transitions
      ID                     <- unique(temp$ID)                                                      # unique IDs
      setkey(temp, ID)      
      temp_long              <- merge(temp,CJ(ID,transition), all.x = F)                             # CJ is like expand.grid or expand_grid but for data.tables, can be used in more flexible ways, also more efficient
      temp_long[,c("from","to") := NULL] #clean join afterwards
      # Assign "from" and To" state based on possible transitions, indexed join in data.table
      temp_long              <- temp_long[dt_trans_keys, on = "transition", nomatch = 0L]            # left join trick on data.table, "on" argument sets the Key automatically
      temp_long              <- temp_long[order(ID, transition)]
    # ----------------------------------- #
    # Task 2. Sample latent arrival times starting from Sick 
    # ----------------------------------- #
      # ----------------------------------- #
      ## Submodule 2.1: Transitions Sick -> Healthy   
      # ----------------------------------- #
      # constant annual rate of becoming Health when Sick
      # Sample from exp. distr.
      if(r_S1H != 0){
        temp_long[transition == 3, T_stop := rexp(.N, rate = r_S1H) + T_n]
      }else {
        temp_long[transition == 3, T_stop := Inf]
      }
      # ----------------------------------- #
      ## Submodule 2.2: Transitions Sick -> Sicker   
      # ----------------------------------- #
      
      # State-residence time-dependent hazard of transition from Sick to Sicker
      # this transition doesn't depend on simulation time, so we can sample it at any point in time.
      
      if (strategy == "B"){
        r_S1S2_scale_ph <- r_S1S2_scale_ph* hr_S1S2_trtB # For Weibull in PH, we just multiply the scale parameter by the hazard ratio, very convenient, for other distributions we need either the non-parmaetric sampling nps approach or we can use the parametric approach via inversion using package nhppp
      }
      r_S1S2_scale_aft <- r_S1S2_scale_ph ^(-1/r_S1S2_shape) # scale, to match the parameterization used in FAEs tutorial
      
      
      if(r_S1S2_shape*r_S1S2_scale_aft != 0){ # case where transition is  allowed
        temp_long[transition == 4, T_stop := rweibull(.N, shape = r_S1S2_shape, scale = r_S1S2_scale_aft) + T_n]
      } else {                                # case where transition is not  allowed
        temp_long[transition == 4, T_stop := Inf]
      }
      
      # ----------------------------------- #
      ## Submodule 2.3: Transitions Sick -> Dead 
      # ----------------------------------- #
      
      # Stratified age-dependent Mortality when Sick              
      dt_time2death_probs_Sick              <- as.data.table(obtain_probs_des(dt_mortality_Sick))
      
      # Subset those rows for which the transition is valid 
      dt_trans5                             <- temp_long[transition == 5]
      
      # Set key on the identifier/strata columns (Age,Sex)                                                            
      setkey(dt_trans5, Age,Sex)                                                                                               # reference dt
      setkey(dt_time2death_probs_Sick, Age,Sex)                                                                                # mortality probs from Sick
      
      # indexed merge, adding columns with interval probabilities
      dt_2sample_time2death_Sick          <- merge(dt_trans5[,Age := ceiling(Age)], dt_time2death_probs_Sick, all.x = TRUE)    # all = TRUE performs a full outer join, all.x performs left-join, all.y right-join
      
      # matrix of probabilities
      m_probs                             <- as.matrix(dt_2sample_time2death_Sick[, .SD, .SDcols = patterns("^DP")])           # choose only the probs cols
      
      # Sample time to death
      dt_trans5[, T_stop_aux := nps_nhppp(m_probs = m_probs, correction = "uniform")]
      
      # ---------------- #
      # Additional correction addressing a minor approximation error with the nps_nhppp method. (ask MLM about it)
      dt_trans5[T_n > Age, tau_aux := T_stop_aux - Age]
      dt_trans5[T_n > Age, T_stop  := T_n + tau_aux]
      dt_trans5[T_n <= Age, T_stop := T_stop_aux]
      
      # ---------------- #
      setkey(temp_long,ID,T_n)
      setkey(dt_trans5,ID,T_n)
      temp_long[transition == 5, T_stop := dt_trans5$T_stop]
      
      # ----------------------------------- #
      # Task 3: Predict the next state 
      # ----------------------------------- #
      # Create an indicator for minimum value in T_stop by ID and calculate tau
      temp_long[, status := as.numeric(T_stop == min(T_stop)), by = ID]
      temp_long[, tau    := T_stop - T_start]
      # ----------------------------------- #
      # Task 4: Update Important Variables
      # ----------------------------------- #
      ### Update Event Time
      temp_long[status == T,T_n := T_stop]
      ### Update Age
      temp_long[status == T,Age := round(T_n)]
      # ----------------------------------- #
      ### Update dt_next_event_long 
      dt_event_history <<- rbind(dt_event_history,temp_long)
      setkey(dt_event_history,ID, T_n, Event_num, transition)
    }
    # ----------------------------------------------------------------------------- #
    # Module 3: transitions from Sicker               ------------------------------ #
    # ----------------------------------------------------------------------------- #
    # Module Set up
    # ----------------------------------- #
    # Task:1 Possible Transitions from Sicker
    # ----------------------------------- #
    temp <- l_subsets_dt_next_event$from_S2
    
    if(length(temp)>0){ #necessary because at some point there may not be S2 people
      
      current_state          <- unique(temp$from)
      transition             <- as.numeric(m_Tr_r[current_state,m_Tr_r[current_state,] !=0])      # possible transitions
      ID                     <- unique(temp$ID)                                                   # unique IDs
      setkey(temp, ID)      
      if(length(transition)>1){
        temp_long            <- merge(temp,CJ(ID,transition), all.x = F) 
      } else if (length(transition)==1){
        temp_long            <- merge(temp,CJ(ID,transition), all.x = T)}                       # changing to all.x = T, so that it can handle only one possible transition
      
      temp_long[,c("from","to") := NULL] #clean join afterwards
      # Assign "from" and To" state based on possible transitions, indexed join in data.table
      temp_long              <- temp_long[dt_trans_keys, on = "transition", nomatch = 0L]         # left join trick on data.table, "on" argument sets the Key automatically
      temp_long              <- temp_long[order(ID, transition)]
      
      # ----------------------------------- #
      ## Submodule  3.1: Transitions Sicker  ->  Dead
      # ----------------------------------- #
      # Stratified age-dependent Mortality when Sicker              
      dt_time2death_probs_Sicker        <- as.data.table(obtain_probs_des(dt_mortality_Sicker))
      
      # Subset those rows for which the transition is valid 
      dt_trans6                         <- temp_long[transition == 6]
      
      # Set key on the identifier/strata columns (Age,Sex)
      setkey(dt_trans6, Age,Sex)                                     # reference dt
      setkey(dt_time2death_probs_Sicker, Age,Sex)                    # mortality probs from Sick
      
      # indexed merge, adding columns with interval probabilities
      dt_2sample_time2death_Sicker      <- merge(dt_trans6[,Age := ceiling(Age)], dt_time2death_probs_Sicker, all.x = TRUE)    # all = TRUE performs a full outer join, all.x performs left-join, all.y right-join
      
      # matrix of probabilities
      m_probs                           <- as.matrix(dt_2sample_time2death_Sicker[, .SD, .SDcols = patterns("^DP")])           # choose only the probs cols
      
      # Sample time to death
      dt_trans6[,T_stop_aux := nps_nhppp(m_probs = m_probs, correction = "uniform")]
      # ---------------- #
      # Additional correction addressing a minor approximation error with the nps_nhppp method. (ask MLM about it)
      dt_trans6[T_n > Age, tau_aux := T_stop_aux - Age]
      dt_trans6[T_n > Age, T_stop  := T_n + tau_aux]
      dt_trans6[T_n <= Age, T_stop := T_stop_aux]
      
      # ---------------- #
      setkey(temp_long,ID,T_n)
      setkey(dt_trans6,ID,T_n)
      temp_long[transition == 6, T_stop := dt_trans6$T_stop]
      
      # ----------------------------------- #
      # Task 3: Predict the next state 
      # ----------------------------------- #
      # Create an indicator for minimum value in T_stop by ID and calculate tau
      temp_long[, status := as.numeric(T_stop == min(T_stop)), by = ID]
      temp_long[,tau := T_stop - T_start]
      # ----------------------------------- #
      # Task 4: Update Important Variables 
      # ----------------------------------- #
      
      ### Update Event Time
      temp_long[status == T,T_n := T_stop]
      ### Update Age
      temp_long[status == T,Age := round(T_n)]
      # ----------------------------------- #
      ### Update dt_next_event_long (global)
      dt_event_history <<- rbind(dt_event_history,temp_long)
      setkey(dt_event_history,ID,T_n,Event_num, transition)
    }
    
    #---------------------------------------------------------------------------- #
    # Next Event Set up
    #---------------------------------------------------------------------------- #
    # # Copy for next iteration
    dt_next_event <- copy(dt_event_history)
    setkey(dt_next_event,ID, T_n, Event_num, transition)
    
    # Sub-setting to work with a shorter dt in next iteration
    # Keep only transitions that occurred
    dt_next_event <- dt_next_event[status == TRUE] 
    # only work with the last event, ignore previous events,
    dt_next_event <- dt_next_event[dt_next_event[, .I[Event_num == max(Event_num, na.rm = TRUE)], by = ID]$V1] 
    #.I extracts the indices of values that meet the condition, like which.max, 
    # but in data.table it also allows for grouping
    
    # slower version, but more understandable--- #
    # # Step 1: Calculate max Event_num for each ID
    # max_event <- dt_next_event[, .(Max_Event = max(Event_num, na.rm = TRUE)), by = ID]
    # # Step 2: Filter the rows in dt_next_event by joining on ID and Event_num
    # dt_next_event <- dt_next_event[max_event, on = .(ID, Event_num = Max_Event)]
    # ------ #
    
    # Handling Events beyond the time horizon (inner)
    # # Those who are alive beyond 100 years or older, assign death status
    dt_next_event[to != "D"  & Age >= time_horizon[2]-1, to :="D"]
    # discard those who died for next event (inner)
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
    return(as.data.table(dt_next_event))
  })
}

#---------------------------------------------------- #
##  Trace matrix form DES output    ----   
#---------------------------------------------------- #
# Function to generate a continuous time trace matrix from the output of the DES, dt_Tn_trace_prop
# Instead of a discrete grid of time we have the sorted time points of all events that occurred in the simulation (event-times, T_n )
# and the proportion of people at each state at each of those event-times.
# Input : 
# dt_fullh: output of DES sim, data table in long format
# v_names_states: vector of state names of the state-transition model we are simulating 
# Output: 
# dt_Tn_trace_prop, Trace matrix, each row represents a distinct event-time. 
# Each column represent a proportion of the cohort that is at a specific state at each specific event-time
# NOTE: event-times are not equidistant to each other and do not have to be integers.

cont_time_trace_sick_sicker <- function(dt_event_history = dt_event_history,
                                        v_names_states  = v_names_states,  
                                        time_init       = time_horizon[1],
                                        proportions     = TRUE, 
                                        rounding_factor = 5, 
                                        long            = F
){
  ## Only Transitions that occurred 
  dt_event_history     <-  dt_event_history[status == TRUE]
  ## Simulation Event-Time, Tn 
  dt_event_history     <-  dt_event_history[,T_n_aux := cumsum(tau)+time_init, by = ID]
  ## Add first row for each ID 
  # Identify the first row for each ID
  first_rows   <- dt_event_history[, .SD[1], by = ID]
  # Modify the 'Value' column in the new rows,e.g., setting it to 0
  first_rows[, ':='(
    transition = 0,
    Event_num  = 0,
    T_n        = time_init, 
    T_start    = time_init,
    tau        = 0,
    status     = T, 
    from       = "H", # this function assumes everyone starts at "H", can be easily modified from any starting State.
    to         = "H", 
    T_stop     = time_init)]
  # Add the new rows for each ID at the bottom of the table
  dt_event_history     <- rbind(dt_event_history, first_rows)
  dt_event_history     <- dt_event_history[order(ID,Event_num)]
  
  ## Set of Indicators 
  # Indicators Into 
  dt_event_history[,':=' (
    into_H             = fifelse(to == "H" ,1,0),
    into_S1            = fifelse(to == "S1",1,0) ,
    into_S2            = fifelse(to == "S2",1,0) ,
    into_D             = fifelse(to == "D" ,1,0)
  ), by = ID]
  dt_event_history[,':=' (
    from_H             = fifelse(from == "H" ,1,0),
    from_S1            = fifelse(from == "S1",1,0),
    from_S2            = fifelse(from == "S2",1,0),
    from_D             = 0
  ), by = ID]
  
  # Making sure we dont count any transition at baseline
  dt_event_history[T_n == time_init,':=' (
    from_H             = 0,
    from_S1            = 0,
    from_S2            = 0,
    from_D             = 0
  ), by = ID]
  
  # Cumulative number of people into and out of each State at each time
  dt_Tn_trace  <- dt_event_history[,
                                   .(
                                     # Into
                                     into_H          = sum(into_H),
                                     into_S1         = sum(into_S1),
                                     into_S2         = sum(into_S2),
                                     into_D          = sum(into_D),
                                     # From 
                                     from_H          = sum(from_H ),
                                     from_S1         = sum(from_S1),
                                     from_S2         = sum(from_S2)
                                   ), by = T_n]
  # Sort by Sim Time 
  setkey(dt_Tn_trace, T_n)
  
  ## Counts per State, in and out 
  dt_Tn_trace[,':=' (
    H_count       = cumsum(into_H  ) - cumsum(from_H ),
    S1_count      = cumsum(into_S1 ) - cumsum(from_S1),
    S2_count      = cumsum(into_S2 ) - cumsum(from_S2),
    D_count       = cumsum(into_D  ) 
  )]
  #Total
  dt_Tn_trace[,Total := as.integer(
    H_count   +
      S1_count  +
      S2_count  +
      D_count    )]
  dt_Tn_trace <- dt_Tn_trace[,.(T_n,H_count,S1_count,S2_count,D_count,Total)]
  # Proportions (rounded)
  dt_Tn_trace_prop <- dt_Tn_trace[,':='(
    H    = round(H_count/Total,rounding_factor),
    S1   = round(S1_count/Total,rounding_factor),
    S2   = round(S2_count/Total,rounding_factor),
    D    = round(D_count/Total,rounding_factor)
  )]
  dt_Tn_trace_prop<- dt_Tn_trace_prop[,.(T_n,H,S1,S2,D,Total)]
  if(proportions == T & long ==F ){
    return(dt_Tn_trace_prop) } else if(proportions == F & long ==F) {
      return(dt_Tn_trace)}else if(proportions == T & long ==T) {
        return(as.data.table(melt(dt_Tn_trace_prop, id.vars = "T_n",
                                  measure.vars = c("H" ,"S1" ,"S2" ,"D"),
                                  variable.name = "Health_State",
                                  value.name = "Proportion")))}else{
                                    return(as.data.table(melt(dt_Tn_trace, id.vars = "T_n",
                                                              measure.vars = c("H_count" ,"S1_count" ,"S2_count" ,"D_count"),
                                                              variable.name = "Health_State",
                                                              value.name = "Count")))
                                  }
  
}
#---------------------------------------------------- #
# Function to generate a discrete time trace matrix from the continuous time matrix, dt_trace_prop
# We need a reference discrete time grid, to find the event-time closest to the discrete time points. 
# Input : 
# cont_time_trace_sick_sicker: output of fun cont_time_trace_sick_sicker, continuous time matrix or data.table in long format
# v_DT_grid : vector of arbitrary length, each element corresponds to a time point in the grid, assumed sorted from start to end.
# Output: 
# dt_Tn_trace_prop, Trace matrix, each row represents a distinct event-time. 
# Each column represent a proportion of the cohort that is at a specific state each specific event-time

discr_time_trace_sick_sicker <- function(dt_CT_trace, v_DT_grid, long = F ){
  v_index <- c()
  for(i in 1:length(v_DT_grid)){
    v_index[i] <- which.min(abs(dt_CT_trace$T_n-v_DT_grid[i]))
  }
  dt_trace_prop_red <- dt_CT_trace[v_index,]
  if(long == F){
    return(dt_trace_prop_red)} else{
      return(melt(dt_trace_prop_red, id.vars = "T_n", 
                  measure.vars = c("H" ,"S1" ,"S2" ,"D"),
                  variable.name = "Health_State",
                  value.name = "Proportion"))
    }
}
#----------------------------------------------------------------------------#
#----------------------------------------------------------------------------#
# CEA Functions ----
#----------------------------------------------------------------------------#
## CEA Analysis 
cea_fn<- function(l_cea_params){
  with(as.list(l_cea_params), {
    # Indicator of strategy and subset to events that occurred 
    dt_SoC <- l_event_history_strategies$SoC[status == T][,Trt := "SoC"]
    dt_A   <- l_event_history_strategies$A[status == T][,Trt := "A"]
    dt_B   <- l_event_history_strategies$B[status == T][,Trt := "B"]
    dt_AB  <- l_event_history_strategies$AB[status == T][,Trt := "AB"]
    # Append event histories of all strategies
    dt_all <- rbind(dt_SoC,
                    dt_A  ,
                    dt_B  ,
                    dt_AB  ) 
    # Free memory
    rm(dt_SoC,
       dt_A  ,
       dt_B  ,
       dt_AB  )
    gc()
    # --- #
    # Assigning costs
    ## for episodes
    dt_all[,Annual_Cost_ep:= fifelse(from == "H",c_H,
                                     fifelse(from == "S1", c_S1,
                                             fifelse(from == "S2", c_S2, NA)))]
    
    ## for treatments
    dt_all[,Annual_Cost_trt:= fifelse(Trt == "A" & (from == "S1" | from  == "S2"), dc_trtA,
                                      fifelse(Trt == "B" & (from == "S1" | from  == "S2"), dc_trtB,
                                              fifelse(Trt == "AB" & (from == "S1" | from  == "S2"), dc_trtA + dc_trtB, 0)))]
    ## for transitions
    dt_all[,Cost_trans:= fifelse(from == "H" & to == "S1",dc_HS1,
                                 fifelse(to == "D", dc_D,0))]
    # --- #
    #  Assigning utilities
    ## for episodes
    dt_all[,Annual_Utility_ep := fifelse(from == "H",u_H,
                                         fifelse(from == "S1", u_S1,
                                                 fifelse(from == "S2", u_S2, NA)))]
    
    ## for treatments
    dt_all[,Annual_Utility_trt:= fifelse(Trt == "A" & from =="S1",du_trtA,
                                         fifelse(Trt == "AB" & from =="S1",du_trtA,0))]
    ## for transitions
    dt_all[,Utility_trans := fifelse(from == "H" & to =="S1", du_HS1,0)]
    
    # --- #
    # Total costs and utilities per event
    dt_all[ ,Total_cost_ep     := Annual_Cost_ep    + Annual_Cost_trt     + Cost_trans]
    dt_all[ ,Total_utility_ep  := Annual_Utility_ep + Annual_Utility_trt + Utility_trans]
    
    # ----- # 
    # Discounted Cost/Utilities
    # Following equation (7) in the manuscript:
    ## Discounting Costs
    ### Exponential discounting, lower and upper bounds per event
    dt_all[,v_dwc_t1  := exp(-(d_c) * (T_start - 25)), by = .(ID,Trt)]
    dt_all[,v_dwc_t2  := exp(-(d_c) * (T_stop - 25)), by = .(ID,Trt)]
    ### Exponential discounting over the event  interval
    dt_all[,v_dwc_t12  := (v_dwc_t1 - v_dwc_t2)/d_c, by = .(ID,Trt)]
    ## Discounting effects
    ### Exponential discounting, lower and upper bounds per event
    dt_all[,v_dwe_t1  := exp(-(d_e) * (T_start - 25)), by = .(ID,Trt)]
    dt_all[,v_dwe_t2  := exp(-(d_e) * (T_stop - 25)), by = .(ID,Trt)]
    ### Exponential discounting over the event interval
    dt_all[,v_dwe_t12  := (v_dwe_t1 - v_dwe_t2)/d_e, by = .(ID,Trt)]
    ## Applying discounting 
    dt_all[ , Disc_total_cost    := Total_cost_ep    * v_dwc_t12 ]
    dt_all[ , Disc_total_utility := Total_utility_ep * v_dwe_t12 ]
    # Lifetime Discounted Costs and Utilities
    # ----- # 
    # Cumulative costs and utilities at each event for each ID
    dt_all[ , Cum_total_cost   := cumsum(Disc_total_cost   ), by = .(ID,Trt)]
    dt_all[ , Cum_total_utility:= cumsum(Disc_total_utility), by = .(ID,Trt)]
    #-------# 
    # Lifetime costs and utilities per ID
    dt_all[ , Lifetime_total_cost   := sum(Disc_total_cost   ), by = .(ID,Trt)]
    dt_all[ , Lifetime_total_utility:= sum(Disc_total_utility), by = .(ID,Trt)]
    # Collapse to a dt with Lifetime costs/utilities per row
    dt_all_lifetime      <- dt_all[, .(Lifetime_total_cost = sum(Disc_total_cost),
                                       Lifetime_total_utility = sum(Disc_total_utility)), by = .(ID, Trt)]
    #-------# 
    # Lifetime expected costs/utilities per strategy
    dt_all_mean_lifetime <- dt_all_lifetime[,.(Mean_total_cost    = mean(Lifetime_total_cost),
                                               Mean_total_utility = mean(Lifetime_total_utility)), by = Trt]
    #-------# 
    ## Incremental Cost-effectiveness Ratios (ICERs) 
    df_cea <- calculate_icers(cost       = dt_all_mean_lifetime$Mean_total_cost,
                              effect     = dt_all_mean_lifetime$Mean_total_utility,
                              strategies = dt_all_mean_lifetime$Trt)
    # Formatted ICER table
    df_cea_format <- format_table_cea(df_cea)
    #return(df_cea_format)
    return(df_cea)
  })
}

#----------------------------------------------------------------------------#

#---------------------------------------------------- #
#####  Discounting       
#---------------------------------------------------- #

#' \code{calc_discount} computes time to metastatic disease.
#'
#' @param x Quantity to be discounted.
#' @param d_factor Annual discount rate.
#' @param elapsed_time Elapsed time.
#' @return 
#' Discounted quantity
#' @export
calc_discount <- function(x, disc_factor, time_init, time_fin){
  if (disc_factor > 0) {
    # x * (1-exp(-disc_factor*(time_fin - time_init)))/disc_factor
    x * (exp(-disc_factor*time_init) - exp(-disc_factor*time_fin) )/disc_factor
  } else{
    x * (time_fin - time_init)
  }
}

# ------------------------------------------------------------------------------------ #
# Epi Outcomes Functions ----
# ------------------------------------------------------------------------------------ #
# Compute Epi Outcomes
epi_fn <- function(l_event_history_strategies = l_event_history_strategies){
  # Adding column of strategy to each dt
  for (i in seq_along(l_event_history_strategies)){
    l_event_history_strategies[[i]]$Strategy <- names(l_event_history_strategies[i])
  }
  # A.1 Dwell-times at each non-death State
  # We use full event histories
  dt_dwell <- bind_rows(l_event_history_strategies)
  dt_dwell <- dt_dwell[ status== TRUE]
  
  dt_dwell <- dt_dwell[, .(
    dwell_mean            =  mean(tau)), by = c("ID","from","Strategy")]
  # A.2 Distribution of the average dwell at each alive state under each strategy.
  dt_dwell_pop        <- dt_dwell[, .(
    mean              =  mean(dwell_mean),
    quantile_2_5pct   =  quantile(dwell_mean, probs = c(.025)),
    quantile_5pct     =  quantile(dwell_mean, probs = c(.05)),
    quantile_25pct    =  quantile(dwell_mean, probs = c(.25)),
    quantile_50pct    =  quantile(dwell_mean, probs = c(.50)),
    quantile_75pct    =  quantile(dwell_mean, probs = c(.75)),
    quantile_95pct    =  quantile(dwell_mean, probs = c(.95)),
    quantile_97_5pct  =  quantile(dwell_mean, probs = c(.975))), 
    by = c("from", "Strategy")]
  
  # B. Cohort survival probabilities and prevalence
  # Generate a continuous time trace matrix for each strategy 
  l_dt_trace_strategies <- list()
  for (i in seq_along(l_event_history_strategies)){
    l_dt_trace_strategies[[i]]  <- cont_time_trace_sick_sicker(
      dt_event_history = l_event_history_strategies[[i]]                    ,
      v_names_states   = l_params$v_names_states,
      time_init        = l_params$time_horizon[1]   ,
      proportions      = T                           ,
      rounding_factor  = 5                               ,
      long             = F
    )
    l_dt_trace_strategies[[i]] <- l_dt_trace_strategies[[i]][T_n <=100] 
    # Compute survival and prevalence from the trace
    l_dt_trace_strategies[[i]][,':='(
      # Check everything sums to 1
      Total      = H+S1+S2+D,          
      # Sum of the proportions of all alive states
      Survival   = H+S1+S2,            
      # (Sum proportions of sick-states) / (sum proportions of alive-states)
      Prevalence = (S1+S2)/(H+S1+S2)   
    )]
  }
  # Assign names of strategies
  names(l_dt_trace_strategies) <- c("SoC","B","A","AB")
  return(list(dt_dwell = dt_dwell, 
              dt_dwell_pop = dt_dwell_pop,
              l_dt_trace_strategies = l_dt_trace_strategies))
}

# Plot Epi Outcomes

plot_epi_outcomes_DES <- function(l_epi_results = epi_fn(l_event_history_strategies = l_event_history_strategies)){
  require(viridisLite)
  require(patchwork)
  # Adding column of strategy to each dt
  l_dt_trace_strategies <- l_epi_results$l_dt_trace_strategies
  for (i in seq_along(l_dt_trace_strategies)){
    l_dt_trace_strategies[[i]]$Strategy <- names(l_dt_trace_strategies[i])
  }
  dt_trace_strategies <- bind_rows(l_dt_trace_strategies)
  # A.  Plot Survival and Prevalence 
  dt_trace_strategies$Strategy <- factor(dt_trace_strategies$Strategy, levels = v_names_strategies )
  dt_trace_strategies$Survival <- round(dt_trace_strategies$Survival, 2)
  
  dt_trace_strategies$Strategy <- factor(dt_trace_strategies$Strategy, levels = v_names_strategies)
  
  # Survival
  p_surv <- ggplot(dt_trace_strategies, 
                   aes(x = T_n, y = Survival, group = Strategy)) +
    geom_line(aes(linetype = Strategy, col = Strategy), linewidth = 1.2) +
    scale_color_viridis_d( option = "magma",begin = 0.1,end = 0.8, direction = 1)+
    xlab("T_n")                      +
    ylab("Proportion")                 +
    ggtitle("Survival probabilities")  +
    theme_bw(base_size = 18)           +
    theme(legend.position = "bottom")  +
    coord_cartesian(xlim = c(25,100), ylim = c(0,1))
  # Prevalence
  p_prev <- ggplot(dt_trace_strategies, 
                   aes(x = T_n, y = Prevalence, group = Strategy)) +
    geom_line(aes(linetype = Strategy, col = Strategy), linewidth = 1.2) +
    scale_color_viridis_d( option = "magma", begin = 0.1,end = 0.8, direction = 1)+
    xlab("T_n")                      +
    ylab("Proportion") + 
    ggtitle(paste("Prevalence", "of", paste(v_names_sick_states, collapse = "&"))) + 
    theme_bw(base_size = 18)           +
    theme(legend.position = "bottom")  +
    coord_cartesian(xlim = c(25,100), ylim = c(0,1))
  
  p_surv_prev <- p_surv + p_prev + plot_layout(guides = "collect", axis = "collect") & theme(legend.position = 'bottom') 
  
  # B. IQR of Dwell times by State and Strategy 
  dt_dwell <- l_epi_results$dt_dwell
  #order Strategies
  
  dt_dwell$Strategy <- factor(dt_dwell$Strategy, levels = v_names_strategies)
  
  # Plot
  
  p_dwell <- ggplot(dt_dwell, aes(x = dwell_mean, y = Strategy, fill = Strategy)) +
    geom_violin(alpha = 0.2) + 
    geom_boxplot(alpha = 0.3,outlier.shape = NA,coef = 1.5) +
    stat_boxplot(geom ='errorbar', width = 0.5)+
    labs(title = "Interquartile Range (IQR) of Mean Dwell Time by State and Strategy",
         x = "Mean Dwell Time (yrs)",
         y = "Strategy") +
    guides(fill = "none")+
    scale_x_continuous(breaks = number_ticks(15)) + 
    theme_bw(base_size = 18) +
    scale_fill_viridis_d( option = "magma",begin = 0.1,end = 0.8, direction = 1)+
    scale_color_viridis_d( option = "magma",begin = 0.1,end = 0.8, direction = 1)+
    theme(legend.position = "") +
    facet_wrap(~from)
  
  # Combine all plots
  p_all_epi <- (p_surv_prev)/(p_dwell) + 
    plot_annotation(title = 'Epidemiological Outcomes',
                    tag_levels = 'A')& theme(
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5)
    ) 
  return(
    list(
      p_surv_prev=p_surv_prev,
      p_all_epi=p_all_epi )
  )
  
}

#------------------------------------------------------------------------------#
#             Generate a PSA input parameter dataset                        ####
#------------------------------------------------------------------------------#
#' Generate parameter sets for the probabilistic sensitivity analysis (PSA)
#'
#' \code{generate_psa_params} generates a PSA dataset of the parameters of the 
#' cost-effectiveness analysis.
#' @param n_sim Number of parameter sets for the PSA dataset
#' @param seed Seed for the random number generation
#' @return A data.frame with a PSA dataset of he parameters of the 
#' cost-effectiveness analysis
#' @export
generate_psa_params_DES <- function(n_sim = 10, seed = 071818){
  set.seed(seed) # set a seed to be able to reproduce the same results
  
  # Weibull parameters
  ## Parameters for multivariate lognormal distribution
  v_means  <- c(0.08, 1.1)    # mean vector
  v_var    <- c(0.02, 0.05)^2 # vector containing the diagonal of covariances
  coef_cor <- 0.5             # correlation coefficient
  m_cor    <- toeplitz(coef_cor^(0:1)) # correlation matrix
  m_r_S1S2_scale_shape <- MethylCapSig::mvlognormal(n = n_sim, 
                                                    Mu = v_means, 
                                                    Sigma = v_var,
                                                    R = m_cor)
  colnames(m_r_S1S2_scale_shape) <- c("r_S1S2_scale", "r_S1S2_shape")
  ## Create data.frame with PSA dataset
  df_psa <- data.frame(
    # Transition probabilities (per cycle)
    r_HS1 = rgamma(n_sim, shape = 30, rate = 170 + 30), # constant rate of becoming Sick when Healthy conditional on surviving
    r_S1H = rgamma(n_sim, shape = 60, rate = 60 + 60),  # constant rate of becoming Healthy when Sick conditional on surviving
    hr_S1 = rlnorm(n_sim, log(3), 0.01),         # rate ratio of death in S1 vs healthy 
    hr_S2 = rlnorm(n_sim, log(10), 0.02),        # rate ratio of death in S2 vs healthy 
    
    # Weibull parameters for state-residence-dependent transition probability of 
    # becoming Sicker when Sick conditional on surviving
    r_S1S2_scale = m_r_S1S2_scale_shape[, "r_S1S2_scale"],  # transition from S1 to S2 - Weibull scale parameter
    r_S1S2_shape = m_r_S1S2_scale_shape[, "r_S1S2_shape"],  # transition from S1 to S2 - Weibull shape parameter
    
    # Effectiveness of treatment B 
    hr_S1S2_trtB = rlnorm(n_sim, meanlog = log(0.6), sdlog = 0.02),    # hazard ratio of becoming Sicker when Sick under B
    
    ## State rewards
    # Costs
    c_H    = rgamma(n_sim, shape = 100,   scale = 20   ),   # cost of remaining one cycle in state H
    c_S1   = rgamma(n_sim, shape = 177.8, scale = 22.5 ),   # cost of remaining one cycle in state S1
    c_S2   = rgamma(n_sim, shape = 225,   scale = 66.7 ),   # cost of remaining one cycle in state S2
    c_trtA = rgamma(n_sim, shape = 73.5, scale = 163.3 ),   # cost of treatment A (per cycle) 
    c_trtB = rgamma(n_sim, shape = 86.2, scale = 150.8 ),   # cost of treatment B (per cycle)
    c_D    = 0,                                             # cost of being in the death state
    
    # Utilities
    u_H    = rbeta(n_sim, shape1 = 200, shape2 = 3     ),   # utility when healthy
    u_S1   = rbeta(n_sim, shape1 = 130, shape2 = 45    ),   # utility when sick
    u_S2   = rbeta(n_sim, shape1 = 230, shape2 = 230   ),   # utility when sicker
    u_D    = 0,                                             # utility when dead
    u_trtA = rbeta(n_sim, shape1 = 300, shape2 = 15    ),   # utility when being treated
    
    # Transition rewards 
    du_HS1 = rbeta(n_sim, shape1 = 11,  shape2 = 1088  ),   # disutility when transitioning from Healthy to Sick
    ic_HS1 = rgamma(n_sim, shape = 25,  scale = 40     ),   # increase in cost when transitioning from Healthy to Sick
    ic_D   = rgamma(n_sim, shape = 100, scale = 20     )    # increase in cost when dying
  )
  return(df_psa)
}
#generate_psa_params_DES(n_sim = 10, seed = 597)


# ---------- #
# PSA Economic Outcomes ----
# ---------- #


cea_psa_fn<- function(l_cea_params, n_wtp = 100000){
  with(as.list(l_cea_params), {
    # Indicator of strategy and subset to events that occurred 
    dt_SoC <- l_event_history_strategies$SoC[status == T][,Trt := "SoC"]
    dt_A   <- l_event_history_strategies$A[status == T][,Trt := "A"]
    dt_B   <- l_event_history_strategies$B[status == T][,Trt := "B"]
    dt_AB  <- l_event_history_strategies$AB[status == T][,Trt := "AB"]
    # Append event histories of all strategies
    dt_all <- rbind(dt_SoC,
                    dt_A  ,
                    dt_B  ,
                    dt_AB  ) 
    # Free memory
    rm(dt_SoC,
       dt_A  ,
       dt_B  ,
       dt_AB  )
    gc()
    # --- #
    # Assigning costs
    ## for episodes
    dt_all[,Annual_Cost_ep:= fifelse(from == "H",c_H,
                                     fifelse(from == "S1", c_S1,
                                             fifelse(from == "S2", c_S2, NA)))]
    
    ## for treatments
    dt_all[,Annual_Cost_trt:= fifelse(Trt == "A" & (from == "S1" | from  == "S2"), dc_trtA,
                                      fifelse(Trt == "B" & (from == "S1" | from  == "S2"), dc_trtB,
                                              fifelse(Trt == "AB" & (from == "S1" | from  == "S2"), dc_trtA + dc_trtB, 0)))]
    ## for transitions
    dt_all[,Cost_trans:= fifelse(from == "H" & to == "S1",dc_HS1,
                                 fifelse(to == "D", dc_D,0))]
    # --- #
    #  Assigning utilities
    ## for episodes
    dt_all[,Annual_Utility_ep := fifelse(from == "H",u_H,
                                         fifelse(from == "S1", u_S1,
                                                 fifelse(from == "S2", u_S2, NA)))]
    
    ## for treatments
    dt_all[,Annual_Utility_trt:= fifelse(Trt == "A" & from =="S1",du_trtA,
                                         fifelse(Trt == "AB" & from =="S1",du_trtA,0))]
    ## for transitions
    dt_all[,Utility_trans := fifelse(from == "H" & to =="S1", du_HS1,0)]
    
    
    # --- #
    # Totaling costs and utilities per event
    dt_all[ ,Total_cost_ep     := Annual_Cost_ep    + Annual_Cost_trt     + Cost_trans]
    dt_all[ ,Total_utility_ep  := Annual_Utility_ep + Annual_Utility_trt + Utility_trans]
    
    # ----- # 
    # Discounted Cost/Utilities
    ## Discounting Costs
    ### Exponential discounting, lower and upper bounds per event
    dt_all[,v_dwc_t1  := exp(-(d_c) * (T_start - 25)), by = .(ID,Trt)]
    dt_all[,v_dwc_t2  := exp(-(d_c) * (T_stop - 25)), by = .(ID,Trt)]
    ### Exponential discounting over the event  interval
    dt_all[,v_dwc_t12  := (v_dwc_t1 - v_dwc_t2)/d_c, by = .(ID,Trt)]
    ## Discounting effects
    ### Exponential discounting, lower and upper bounds per event
    dt_all[,v_dwe_t1  := exp(-(d_e) * (T_start - 25)), by = .(ID,Trt)]
    dt_all[,v_dwe_t2  := exp(-(d_e) * (T_stop - 25)), by = .(ID,Trt)]
    ### Exponential discounting over the event interval
    dt_all[,v_dwe_t12  := (v_dwe_t1 - v_dwe_t2)/d_e, by = .(ID,Trt)]
    ## Applying discounting 
    dt_all[ , Disc_total_cost    := Total_cost_ep    * v_dwc_t12 ]
    dt_all[ , Disc_total_utility := Total_utility_ep * v_dwe_t12 ]
    # Lifetime Discounted Costs and Utilities
    # ----- # 
    # Cumulative costs and utilities at each event for each ID
    dt_all[ , Cum_total_cost   := cumsum(Disc_total_cost   ), by = .(ID,Trt)]
    dt_all[ , Cum_total_utility:= cumsum(Disc_total_utility), by = .(ID,Trt)]
    #-------# 
    #Lifetime costs and utilities per ID
    dt_all[ , Lifetime_total_cost   := sum(Disc_total_cost   ), by = .(ID,Trt)]
    dt_all[ , Lifetime_total_utility:= sum(Disc_total_utility), by = .(ID,Trt)]
    # Collapse to a dt with Lifetime costs/utilities per row
    dt_all_lifetime      <- dt_all[, .(Lifetime_total_cost = sum(Disc_total_cost),
                                       Lifetime_total_utility = sum(Disc_total_utility)), by = .(ID, Trt)]
    #-------# 
    # Lifetime expected costs/utilities per strategy
    dt_all_mean_lifetime <- dt_all_lifetime[,.(Mean_total_cost    = mean(Lifetime_total_cost),
                                               Mean_total_utility = mean(Lifetime_total_utility)), by = Trt]
    #-------# 
    # Vector Net Monetary Benefit 
    v_nmb <- dt_all_mean_lifetime$Mean_total_utility * n_wtp - dt_all_mean_lifetime$Mean_total_cost
    #-------# 
    
    ## Incremental cost-effectiveness ratios (ICERs) 
    df_ce <- data.frame(Cost       = dt_all_mean_lifetime$Mean_total_cost,
                        Effect     = dt_all_mean_lifetime$Mean_total_utility,
                        Strategy   = dt_all_mean_lifetime$Trt, 
                        NMB        = v_nmb)
    return(df_ce)
  })
}


# Distributed Implementation ----

#---------------------------------------------------- #
#####  Run Sick-Sicker Model    ----   
#---------------------------------------------------- #

sick_sicker_model   <- function(l_params){
  while(nrow(l_params$dt_next_event) >0){                                                     
    l_params$dt_next_event <- sim_next_event(l_params)                           
  }
  return(dt_event_history)
} 

# ------------------------------------------------------------------------------------ #
# Obtain trace directly from model
# ------------------------------------------------------------------------------------ #
obtain_model_grid_trace <- function(l_params, grid_cycle_length = 6/12){
  #Generate continuous time trace matrix from model output
  dt_CT_trace <- cont_time_trace_sick_sicker(dt_event_history = as.data.table(sick_sicker_model(l_params))                    ,
                                             v_names_states   = l_params$v_names_states,
                                             time_init        = l_params$time_horizon[1]   ,
                                             proportions      = T                           ,
                                             rounding_factor  = 5                               ,
                                             long             = F )
  # Define the grid for discrete time points
  cycle_length <- grid_cycle_length # cycle length equal to one year (use 1/12 for monthly)
  n_age_init   <- l_params$time_horizon[1]  # age at baseline
  n_age_max    <- l_params$time_horizon[2] # maximum age of follow up
  n_cycles     <- (n_age_max - n_age_init)/cycle_length # time horizon, number of cycles
  grid         <- seq(from =n_age_init, to= n_age_max, length.out = n_cycles)
  
  # generate discrete time trace matrix
  dt_DT_trace <- discr_time_trace_sick_sicker(dt_CT_trace = dt_CT_trace,
                                              v_DT_grid   = grid, 
                                              long =F)
  return(dt_DT_trace)   
}

# ------------------------------------------------------------------------------------ #
# Function to run the model in parallel ----
# ------------------------------------------------------------------------------------ #

run_parallel_model <- function(model_fn, l_inputs, k) {
  start_time <- Sys.time()
  message(sprintf("Starting all runs at %s...", start_time))
  
  l_grid_trace <- foreach(i = 1:k, .combine = list,.multicombine = TRUE,
                          
                          .export = c("sick_sicker_model", "sim_next_event", 
                                      "cont_time_trace_sick_sicker", "discr_time_trace_sick_sicker",
                                      "obtain_probs_des","nps_nhppp"),
                          .packages = c("data.table","dplyr")) %dopar% {
                            dt_event_history   <- data.table()
                            l_grid_trace       <- try(as.matrix(model_fn(l_inputs) %>% 
                                                                  dplyr::select("T_n","H","S1","S2","D")))
                          }
  a_grid_trace         <- abind::abind(l_grid_trace, along = 3)
  a_grid_trace_byState <- aperm(a_grid_trace, c(1, 3, 2)) # array order by state on the third dimension, easier to compute conf int around mean estimator
  
  end_time <- Sys.time()
  message(sprintf("All runs completed at %s. Total time elapsed: %s", 
                  end_time, difftime(end_time, start_time, units = "secs")))
  #return(a_grid_trace)
  return(a_grid_trace_byState)
}

# ------------------------------------------------------------------------------------ #
# Post-processing Distributed DES
# ------------------------------------------------------------------------------------ #

batch_summary <- function(array_grid, k, n_sim){
  #Time grid
  time_grid     <- rowSums(array_grid[,,1])/k
  # Proportion at H
  dt_trace_H     <-  as.data.table(array_grid[,,2] )
  dt_trace_H[,':='(
    Mean       = rowSums(.SD)/k,
    SE         = apply(.SD, 1, function(x) sd(x) / sqrt(length(x))),
    Q025       = apply(.SD, 1, function(x) quantile(x, 0.025)),
    Q25        = apply(.SD, 1, function(x) quantile(x, 0.25)),
    Q50        = apply(.SD, 1, function(x) quantile(x, 0.5)) ,
    Q75        = apply(.SD, 1, function(x) quantile(x, 0.75)),
    Q975        = apply(.SD, 1, function(x) quantile(x, 0.975))
    
  )][,':='(
    n_sim      = n_sim,
    State      = "Healthy",
    T_n        = time_grid
  )] 
  Prop_H <- dt_trace_H[, .SD, .SDcols = c("n_sim","T_n","State","Mean","SE","Q025","Q25","Q50","Q75","Q975")]
  # Proportion at S1
  dt_trace_S1    <-  as.data.table(array_grid[,,3] )
  dt_trace_S1[,':='(
    Mean       = rowSums(.SD)/k,
    SE         = apply(.SD, 1, function(x) sd(x) / sqrt(length(x))),
    Q025       = apply(.SD, 1, function(x) quantile(x, 0.025)),
    Q25        = apply(.SD, 1, function(x) quantile(x, 0.25)),
    Q50        = apply(.SD, 1, function(x) quantile(x, 0.5)) ,
    Q75        = apply(.SD, 1, function(x) quantile(x, 0.75)),
    Q975       = apply(.SD, 1, function(x) quantile(x, 0.975))
    
  )][,':='(
    n_sim      = n_sim,
    State      = "Sick",
    T_n        = time_grid
  )] 
  Prop_S1 <- dt_trace_S1[, .SD, .SDcols = c("n_sim","T_n","State","Mean","SE","Q025","Q25","Q50","Q75","Q975")]
  # Proportion at S2
  dt_trace_S2    <-  as.data.table(array_grid[,,4] )
  dt_trace_S2[,':='(
    Mean       = rowSums(.SD)/k,
    SE         = apply(.SD, 1, function(x) sd(x) / sqrt(length(x))),
    Q025       = apply(.SD, 1, function(x) quantile(x, 0.025)),
    Q25        = apply(.SD, 1, function(x) quantile(x, 0.25)),
    Q50        = apply(.SD, 1, function(x) quantile(x, 0.5)) ,
    Q75        = apply(.SD, 1, function(x) quantile(x, 0.75)),
    Q975        = apply(.SD, 1, function(x) quantile(x, 0.975))
    
  )][,':='(
    n_sim      = n_sim,
    State      = "Sicker",
    T_n        = time_grid
  )] 
  Prop_S2 <- dt_trace_S2[, .SD, .SDcols = c("n_sim","T_n","State","Mean","SE","Q025","Q25","Q50","Q75","Q975")]
  
  # Proportion at D
  dt_trace_D     <-  as.data.table(array_grid[,,5] )
  dt_trace_D[,':='(
    Mean       = rowSums(.SD)/k,
    SE         = apply(.SD, 1, function(x) sd(x) / sqrt(length(x))),
    Q025       = apply(.SD, 1, function(x) quantile(x, 0.025)),
    Q25        = apply(.SD, 1, function(x) quantile(x, 0.25)),
    Q50        = apply(.SD, 1, function(x) quantile(x, 0.5)) ,
    Q75        = apply(.SD, 1, function(x) quantile(x, 0.75)),
    Q975        = apply(.SD, 1, function(x) quantile(x, 0.975))
    
  )][,':='(
    n_sim      = n_sim,
    State      = "Dead",
    T_n        = time_grid
  )] 
  Prop_D <- dt_trace_D[, .SD, .SDcols = c("n_sim","T_n","State","Mean","SE","Q025","Q25","Q50","Q75","Q975")]
  dt_bystate <-  as.data.table(rbind(Prop_H,Prop_S1,Prop_S2,Prop_D))
  return(dt_bystate)
}

