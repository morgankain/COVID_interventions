fixed_params        <- fixed_params[!(names(fixed_params) %in% c("Ca", "alpha", "delta", "mu", "rho_d"))]

set.seed(1022)

## parameters that will vary
variable_params <- sobolDesign(
  lower = c(
    sim_start   = location_params[location_params$Parameter == "sim_start", ]$lwr
  , Ca          = params[params$Parameter == "Ca", ]$lwr
  , alpha       = params[params$Parameter == "alpha", ]$lwr
  , delta       = params[params$Parameter == "delta", ]$lwr
  , mu          = params[params$Parameter == "mu", ]$lwr
  , rho_d       = params[params$Parameter == "rho_d", ]$lwr
  )
, upper = c(
    sim_start   = location_params[location_params$Parameter == "sim_start", ]$upr
  , Ca          = params[params$Parameter == "Ca", ]$upr
  , alpha       = params[params$Parameter == "alpha", ]$upr
  , delta       = params[params$Parameter == "delta", ]$upr
  , mu          = params[params$Parameter == "mu", ]$upr
  , rho_d       = params[params$Parameter == "rho_d", ]$upr
)
, nseq  = nparams
) %>% mutate(
  sim_start   = round(sim_start),
  sim_start  = as.Date(sim_start, origin = as.Date("2020-01-01"))
 ) %>% mutate(paramset = seq(1, nparams)) %>% 
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

nparams           <- nrow(variable_params)

