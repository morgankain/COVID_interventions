fixed_params        <- fixed_params[-c(1, 5, 11, 12)]

## parameters that will vary
variable_params <- sobolDesign(
  lower = c(
    E0          = 3
  , sim_start   = 14
  , int_start1  = 60
  , int_length2 = 60
  , soc_dist_level_wfh       = 0.7
  , soc_dist_level_sip = ifelse(!fit_to_sip, 0.10, NA)
  , Ca          = params[params$Parameter == "Ca", ]$lwr
  , alpha       = params[params$Parameter == "alpha", ]$lwr
  , delta       = params[params$Parameter == "delta", ]$lwr
  , mu          = params[params$Parameter == "mu", ]$lwr
  , fixed_params["rho_d"] - fixed_params["rho_d"]/3
  , fixed_params["rho_r"] - fixed_params["rho_r"]/3
  )
, upper = c(
    E0          = 6
  , sim_start   = 28
  , int_start1  = 69
  , int_length2 = 120
  , soc_dist_level_wfh       = 0.9
  , soc_dist_level_sip = ifelse(!fit_to_sip, 0.40, NA)
  , Ca          = params[params$Parameter == "Ca", ]$upr
  , alpha       = params[params$Parameter == "alpha", ]$upr
  , delta       = params[params$Parameter == "delta", ]$upr
  , mu          = params[params$Parameter == "mu", ]$upr
  , fixed_params["rho_d"] + fixed_params["rho_d"]/3
  , fixed_params["rho_r"] + fixed_params["rho_r"]/3
)
, nseq  = nparams
) %>% mutate(
  sim_start   = round(sim_start)
, int_start1  = round(int_start1)
, int_length2 = round(int_length2)
, E0 = round(E0))

if (fit_to_sip) {
 variable_params <- variable_params %>% dplyr::select(-soc_dist_level_sip) 
}

variable_params <- variable_params %>% mutate(
    sim_start  = as.Date(sim_start, origin = as.Date("2020-01-01"))
  , int_start1 = as.Date(round(int_start1), origin = as.Date("2020-01-01"))
  , int_start2 = as.Date("2020-03-17")
 ) %>% mutate(
    int_length1  = round(int_start2 - int_start1)
  ## Column for estimated beta0 (!! see lower down for desire to store likelyhood profile to define pmcmc runs)
  , beta0est     = 0
  , paramset     = seq(1, nparams)
) 

variable_params <- variable_params %>% 
  mutate(
    ## Previously had infected isolation starting after 60 days, changing this to when the SD gets lifted a bit
    iso_start = int_start2 + int_length2
  , log_lik   = 0
  , R0        = 0
  , Reff      = 0
  ) 

fixed_params <- fixed_params[-c(9, 10)]
