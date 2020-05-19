## Mot updated in a while...
fixed_params        <- fixed_params[-c(1, 5, 7, 8, 9, 11, 12)]

variable_params <- sobolDesign(
  lower = c(
    E0          = ifelse(!fit.E0
    , location_params[location_params$Parameter == "E0", ]$lwr
    , NA)
  , sim_start   = location_params[location_params$Parameter == "sim_start", ]$lwr
  , int_start1  = location_params[location_params$Parameter == "int_start1", ]$lwr
### For Santa Clara county. !!!! Model will need to be updated when we exit this reality
  , int_length2 = 76
  , soc_dist_level_wfh = location_params[location_params$Parameter == "soc_dist_level_wfh", ]$lwr
  , soc_dist_level_sip = ifelse(!fit_to_sip
    , location_params[location_params$Parameter == "soc_dist_level_sip", ]$lwr
    , NA)
  , Ca          = params[params$Parameter == "Ca", ]$lwr
  , alpha       = params[params$Parameter == "alpha", ]$lwr
  , delta       = params[params$Parameter == "delta", ]$lwr
  , mu          = params[params$Parameter == "mu", ]$lwr
  , lambda_a    = params[params$Parameter == "lambda_a", ]$lwr
  , lambda_s    = params[params$Parameter == "lambda_s", ]$lwr
  , lambda_m    = params[params$Parameter == "lambda_m", ]$lwr
  )
, upper = c(
    E0          = ifelse(!fit.E0
    , location_params[location_params$Parameter == "E0", ]$upr
    , NA)
  , sim_start   = location_params[location_params$Parameter == "sim_start", ]$upr
  , int_start1  = location_params[location_params$Parameter == "int_start1", ]$upr
### For Santa Clara county. !!!! Model will need to be updated when we exit this reality
  , int_length2 = 76
  , soc_dist_level_wfh = location_params[location_params$Parameter == "soc_dist_level_wfh", ]$upr
  , soc_dist_level_sip = ifelse(!fit_to_sip
    , location_params[location_params$Parameter == "soc_dist_level_sip", ]$upr
    , NA)
  , Ca          = params[params$Parameter == "Ca", ]$upr
  , alpha       = params[params$Parameter == "alpha", ]$upr
  , delta       = params[params$Parameter == "delta", ]$upr
  , mu          = params[params$Parameter == "mu", ]$upr
  , lambda_a    = params[params$Parameter == "lambda_a", ]$upr
  , lambda_s    = params[params$Parameter == "lambda_s", ]$upr
  , lambda_m    = params[params$Parameter == "lambda_m", ]$upr
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
if (fit.E0) {
variable_params <- variable_params %>% dplyr::select(-E0)
}

variable_params <- variable_params %>% mutate(
    sim_start  = as.Date(sim_start, origin = as.Date("2020-01-01"))
  , int_start1 = as.Date(round(int_start1), origin = as.Date("2020-01-01"))
  , int_start2 = as.Date(round(location_params[location_params$Parameter == "int_start2", ]$est)
    , origin = as.Date("2020-01-01")) 
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

if (fit_to_sip) {
variable_params <- variable_params %>% mutate(E_init = 0)
}
if (fit.E0) {
variable_params <- variable_params %>% mutate(soc_dist_level_sip = 0)
}
