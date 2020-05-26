fixed_params        <- fixed_params[-c(1, 5, 11, 12)]

## parameters that will vary
variable_params <- sobolDesign(
  lower = c(
    E0          = ifelse(!fit.E0
    , location_params[location_params$Parameter == "E0", ]$lwr
    , NA)
  , sim_start   = location_params[location_params$Parameter == "sim_start", ]$lwr
  , Ca          = params[params$Parameter == "Ca", ]$lwr
  , alpha       = params[params$Parameter == "alpha", ]$lwr
  , delta       = params[params$Parameter == "delta", ]$lwr
  , mu          = params[params$Parameter == "mu", ]$lwr
  )
, upper = c(
    E0          = ifelse(!fit.E0
    , location_params[location_params$Parameter == "E0", ]$upr
    , NA)
  , sim_start   = location_params[location_params$Parameter == "sim_start", ]$upr
  , Ca          = params[params$Parameter == "Ca", ]$upr
  , alpha       = params[params$Parameter == "alpha", ]$upr
  , delta       = params[params$Parameter == "delta", ]$upr
  , mu          = params[params$Parameter == "mu", ]$upr
)
, nseq  = nparams
) %>% mutate(
  sim_start   = round(sim_start)
, E0 = round(E0))

if (fit.E0) {
variable_params <- variable_params %>% dplyr::select(-E0)
}

variable_params <- variable_params %>% mutate(
    sim_start  = as.Date(sim_start, origin = as.Date("2020-01-01"))
 ) %>% mutate(paramset = seq(1, nparams)) 

variable_params <- variable_params %>% 
  mutate(
    log_lik        = 0
  , log_lik_se     = 0
  , beta0est       = 0
  , E_init         = 0
  , detect_k       = 0
  , detect_mid     = 0
  , theta          = 0
  , theta2         = 0
  , beta_min       = 0
  ) 
