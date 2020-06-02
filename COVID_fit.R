##################################
## Fit COVID epidemic with pomp ##
##################################

set.seed(10001)
fitting           <- TRUE     ## Small change in pomp objects if fitting or simulating
fit.minus         <- 0        ## Use data until X days prior to the present
more.params.uncer <- FALSE    ## Fit with more (FALSE) or fewer (TRUE) point estimates for a number of parameters
fit.E0            <- TRUE     ## Also fit initial # that starts the epidemic?
## more.params.uncer = FALSE is more supported, uses parameter ranges with more research and reacts to choice of focal.county if possible
## !!!!! For FALSE update parameters in location_params.csv
## more.params.uncer = TRUE  is less suppored, raw parameter values that can be adjusted manually
usable.cores      <- 2        ## Number of cores to use to fit
fit.with          <- "D"      ## Fit with D (deaths) or H (hospitalizations) 
fit_to_sip        <- TRUE     ## Fit beta0 and shelter in place simultaneously?
import_cases      <- FALSE    ## Use importation of cases?
n.mif_runs        <- 2        ## mif2 fitting parameters
n.mif_length      <- 100
n.mif_particles   <- 600
n.mif_rw.sd       <- 0.02
focal.county      <- "King"  ## County to fit to
## !!! Curently parameters exist for Santa Clara, Miami-Dade, New York City, King, Los Angeles
## !!! But only Santa Clara explored
# county.N        <- 1.938e6         ## County population size
## !!! Now contained within location_params.csv
nparams           <- 5               ## number of parameter sobol samples (more = longer)
nsim              <- 200             ## number of simulations for each fitted beta0 for dynamics

needed_packages <- c(
    "pomp"
  , "plyr"
  , "dplyr"
  , "ggplot2"
  , "magrittr"
  , "scales"
  , "lubridate"
  , "tidyr"
  , "foreach"
  , "doParallel"
  , "data.table")

lapply(needed_packages, require, character.only = TRUE)

## Be very careful here, adjust according to your machine
registerDoParallel(cores = usable.cores)

## Bring in pomp objects
source("COVID_pomp.R")
 
if (fit.with == "D") {
#deaths <- fread("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv")
deaths  <- read.csv("us-counties.txt")
deaths  <- deaths %>% mutate(date = as.Date(date)) %>% filter(county == focal.county)
deaths  <- deaths %>% dplyr::filter(date < max(date) - fit.minus)
} else if (fit.with == "H") {
## !! Not supported right now for SCC, but placing here for completeness
   ## !! Right now only usable for CCC
hospit     <- read.csv("contra_costa/ccc_data.csv")
hospit     <- hospit %>% 
  mutate(date = as.Date(REPORT_DATE)) %>% 
  filter(CURRENT_HOSPITALIZED != "NULL") %>% 
  mutate(ch = as.numeric(as.character(CURRENT_HOSPITALIZED))) %>% 
  dplyr::select(date, ch)
hospit    <- hospit %>% dplyr::filter(date < max(date) - fit.minus)  
}

if (!more.params.uncer) {
params <- read.csv("params.csv", stringsAsFactors = FALSE)
} else {
params <- read.csv("params_wide.csv", stringsAsFactors = FALSE) 
}
params <- params %>% mutate(Value = sapply(est, function(x) eval(parse(text = x))))

fixed_params        <- params$Value
names(fixed_params) <- params$Parameter
if (!import_cases) {fixed_params["import_rate"] <- 0}

location_params     <- read.csv("location_params.csv", stringsAsFactors = FALSE)
location_params     <- location_params %>% filter(location == focal.county)

## debug
#location_params[location_params$Parameter == "sim_start", ]$lwr <- 38
#location_params[location_params$Parameter == "sim_start", ]$upr <- 53

fixed_params        <- c(fixed_params
  , N = location_params[location_params$Parameter == "N", ]$est)

if (more.params.uncer) {
source("variable_params_more.R")
} else {
source("variable_params_less.R")
}

## Run parameters
sim_start  <- variable_params$sim_start
sim_length <- 500
sim_end    <- sim_start + sim_length

## containers for results
SEIR.sim.ss.t.ci <- data.frame(
  name     = character(0)
, lwr      = numeric(0)
, est      = numeric(0)
, upr      = numeric(0)
, paramset = numeric(0))

if (fit_to_sip) {
## beta0, soc_dist_level_sip, loglik
if (fit.E0) {
param_array <- array(
  data = 0
, dim  = c(nparams, n.mif_runs, 4))
dimnames(param_array)[[3]] <- c("beta0", "soc_dist_level_sip", "loglik", "E_init")  
} else {
param_array <- array(
  data = 0
, dim  = c(nparams, n.mif_runs, 3))
dimnames(param_array)[[3]] <- c("beta0", "soc_dist_level_sip", "loglik")  
}
} else {
## beta0, loglik
if (fit.E0) {
param_array <- array(
  data = 0
, dim  = c(nparams, n.mif_runs, 3))  
dimnames(param_array)[[3]] <- c("beta0", "loglik", "E_init")  
} else {
param_array <- array(
  data = 0
, dim  = c(nparams, n.mif_runs, 2))  
dimnames(param_array)[[3]] <- c("beta0", "loglik")  
}
}

for (i in 1:nrow(variable_params)) {

if (fit.with == "D") {
  
  ## Adjust the data for the current start date
county.data <- deaths %>% 
  mutate(day = as.numeric(date - variable_params[i, ]$sim_start)) %>% 
  filter(day > 0) %>%
  dplyr::select(day, date, deaths) %>% 
  mutate(deaths_cum = deaths) %>% 
  mutate(deaths = deaths_cum - lag(deaths_cum)) %>% 
  replace_na(list(deaths = 0)) %>%
  dplyr::select(-deaths_cum)
  
## Add days from the start of the sim to the first recorded day in the dataset
#county.data <- rbind(
#  data.frame(
#    day    = seq(1:(min(county.data$day) - 1))
#  , date   = as.Date(seq(1:(min(county.data$day) - 1)), origin = variable_params[i, ]$sim_start)
#  , deaths = 0
#  )
#, county.data
#  )

} else if (fit.with == "H") {
  
## Adjust the data for the current start date
county.data <- hospit %>% mutate(day = as.numeric(date - variable_params[i, ]$sim_start)) 
names(county.data)[2] <- "hosp"  

}

## Create intervention covariate table for the full forecast
if (fit_to_sip) {
  intervention.forecast <- with(variable_params[i, ], {

 covariate_table(
  day              = 1:sim_length
  
, intervention     = c(
      # No intervention until intervention start time
    rep(0, int_start1 - sim_start)                   
      # Intervention style 1
  , rep(1, int_length1)
      # Intervention style 2
  , rep(2, sim_length - (int_start2 - sim_start))
      # Post intervention close
)
  
, isolation        = rep(NA, sim_length)
, iso_severe_level = rep(NA, sim_length)      # % of contats that severe cases maintain
, iso_mild_level   = rep(NA, sim_length)   # % of contats that mild cases maintain

, soc_dist_level_wfh = rep(soc_dist_level_wfh, sim_length) 

, thresh_H_start   = rep(NA, sim_length)
, thresh_H_end     = rep(NA, sim_length)

, thresh_int_level = rep(NA, sim_length)
, back_int_level   = rep(NA, sim_length)
  
, order            = "constant"
, times            = "day"
  )

})
} else {
  intervention.forecast <- with(variable_params[i, ], {

 covariate_table(
  day              = 1:sim_length
  
, intervention     = c(
      # No intervention until intervention start time
    rep(0, int_start1 - sim_start)                   
      # Intervention style 1
  , rep(1, int_length1)
      # Intervention style 2
  , rep(2, sim_length - (int_start2 - sim_start))
      # Post intervention close
)
  
, isolation        = rep(NA, sim_length)
, iso_severe_level = rep(NA, sim_length)      # % of contats that severe cases maintain
, iso_mild_level   = rep(NA, sim_length)   # % of contats that mild cases maintain

, soc_dist_level_wfh = rep(soc_dist_level_wfh, sim_length) 
, soc_dist_level_sip = rep(soc_dist_level_sip, sim_length)

, thresh_H_start   = rep(NA, sim_length)
, thresh_H_end     = rep(NA, sim_length)

, thresh_int_level = rep(NA, sim_length)
, back_int_level   = rep(NA, sim_length)
  
, order            = "constant"
, times            = "day"
  )

})  
}

covid.fitting <- county.data %>%
  pomp(
    time       = "day"
  , t0         = 1
  , covar      = intervention.forecast
  , rprocess   = euler(sir_step, delta.t = 1/6)
  , rmeasure   = { 
   if (fit.with == "D") { rmeas_deaths } else if (fit.with == "H") { rmeas_hosp }
  }
  , dmeasure   = {
   if (fit.with == "D") { dmeas_deaths } else if (fit.with == "H") { dmeas_hosp }
  }
  , rinit      = sir_init
  , partrans   = par_trans
  , accumvars  = accum_names
  , paramnames = param_names
  , statenames = state_names) 

if (variable_params[i, ]$beta0est == 0) {

if (!more.params.uncer) {

checktime <- system.time({
mifs_local <- foreach(j = 1:n.mif_runs, .combine = c) %dopar%  {
    
library(pomp)
library(dplyr)

  covid.fitting %>%
  mif2(
    t0      = 1
  , params  = c(
    c(fixed_params
    , Ca    = variable_params[i, ]$Ca
    , alpha = variable_params[i, ]$alpha
    , delta = variable_params[i, ]$delta
    , mu    = variable_params[i, ]$mu
      )
  , {
    if (fit_to_sip) {
      if (fit.E0) {
    c(beta0              = rlnorm(1, log(0.7), 0.17)
    , soc_dist_level_sip = rlnorm(1, log(0.2), 0.2)
    , E_init             = 2)
      } else {
    c(beta0              = rlnorm(1, log(0.7), 0.17)
    , soc_dist_level_sip = rlnorm(1, log(0.2), 0.2)
    , E0                 = variable_params[i, ]$E0)        
      }
    } else {
      if (fit.E0) {
    c(beta0              = rlnorm(1, log(0.7), 0.17)
    , E_init             = 2)
      } else {
    c(beta0              = rlnorm(1, log(0.7), 0.17)
    , E0                 = variable_params[i, ]$E0)        
      }
    }
  }
  )
  , Np     = n.mif_particles
  , Nmif   = n.mif_length
  , cooling.fraction.50 = 0.5
  , rw.sd  = {
    if (fit_to_sip) {
      if (fit.E0) {
    rw.sd(beta0 = n.mif_rw.sd, soc_dist_level_sip = n.mif_rw.sd, E_init = n.mif_rw.sd)        
      } else {
    rw.sd(beta0 = n.mif_rw.sd, soc_dist_level_sip = n.mif_rw.sd)        
      }
    } else {
      if (fit.E0) {
    rw.sd(beta0 = n.mif_rw.sd, E_init = n.mif_rw.sd)        
      } else {
    rw.sd(beta0 = n.mif_rw.sd)        
      }
    }
  }
        )

}

})

if (fit_to_sip) {
 variable_params[i, "soc_dist_level_sip"] <- mean(coef(mifs_local)[which(dimnames(coef(mifs_local))$parameter == "soc_dist_level_sip"), ])
 param_array[i,,"soc_dist_level_sip"]     <- coef(mifs_local)[which(dimnames(coef(mifs_local))$parameter == "soc_dist_level_sip"), ]
}
if (fit.E0) {
variable_params[i, "E_init"] <- mean(coef(mifs_local)[which(dimnames(coef(mifs_local))$parameter == "E_init"), ])
param_array[i,,"E_init"]     <- coef(mifs_local)[which(dimnames(coef(mifs_local))$parameter == "E_init"), ]  
}

} else {
  
checktime  <- system.time({
mifs_local <- foreach(j = 1:n.mif_runs, .combine = c) %dopar%  {
    
library(pomp)
library(dplyr)

  covid.fitting %>%
  mif2(
    t0     = 1
  , params = c(
    c(fixed_params
      , Ca       = variable_params[i, ]$Ca
      , alpha    = variable_params[i, ]$alpha
      , delta    = variable_params[i, ]$delta
      , mu       = variable_params[i, ]$mu
      , lambda_a = variable_params[i, ]$lambda_a
      , lambda_s = variable_params[i, ]$lambda_s
      , lambda_m = variable_params[i, ]$lambda_m  
      )
  , {
    if (fit_to_sip) {
      if (fit.E0) {
    c(beta0              = rlnorm(1, log(0.7), 0.17)
    , soc_dist_level_sip = rlnorm(1, log(0.2), 0.2)
    , E_init             = rpois(1, 2))
      } else {
    c(beta0              = rlnorm(1, log(0.7), 0.17)
    , soc_dist_level_sip = rlnorm(1, log(0.2), 0.2)
    , E0                 = variable_params[i, ]$E0)        
      }
    } else {
      if (fit.E0) {
    c(beta0              = rlnorm(1, log(0.7), 0.17)
    , E_init             = rpois(1, 2))
      } else {
    c(beta0              = rlnorm(1, log(0.7), 0.17)
    , E0                 = variable_params[i, ]$E0)        
      }
    }
  }
  )
  , Np     = n.mif_particles
  , Nmif   = n.mif_length
  , cooling.fraction.50 = 0.5
  , rw.sd  = {
    if (fit_to_sip) {
      if (fit.E0) {
    rw.sd(beta0 = n.mif_rw.sd, soc_dist_level_sip = n.mif_rw.sd, E_init = n.mif_rw.sd)        
      } else {
    rw.sd(beta0 = n.mif_rw.sd, soc_dist_level_sip = n.mif_rw.sd)        
      }
    } else {
      if (fit.E0) {
    rw.sd(beta0 = n.mif_rw.sd, E_init = n.mif_rw.sd)        
      } else {
    rw.sd(beta0 = n.mif_rw.sd)        
      }
    }
  }
        )

}

})  
  
}
  
gg.fit <- mifs_local %>%
  traces() %>%
  melt() %>%
  filter(
    variable == "loglik" | 
    variable == "beta0" | 
    variable == "soc_dist_level_sip" |
    variable == "E_init") %>% 
  ggplot(aes(x = iteration, y = value, group = L1, colour = factor(L1)))+
  geom_line() +
  guides(color = FALSE) +
  facet_wrap(~variable, scales = "free_y") +
  theme_bw()

## !! Still not using uncertainty in beta which should be corrected soon
variable_params[i, "beta0est"] <- mean(coef(mifs_local)[which(dimnames(coef(mifs_local))$parameter == "beta0"), ])
param_array[i,,"beta0"]        <- coef(mifs_local)[which(dimnames(coef(mifs_local))$parameter == "beta0"), ] 

loglik.out    <- numeric(length(mifs_local))
for (k in seq_along(loglik.out)) {
loglik.out[k] <- mifs_local[[k]]@loglik
}

variable_params[i, "log_lik"]  <- mean(loglik.out)
param_array[i,,"loglik"]       <- loglik.out

}

SEIR.sim <- do.call(
  pomp::simulate
  , list(
    object         = covid.fitting
    , times        = intervention.forecast@times
    , params = {
    if (!more.params.uncer) {
      c(fixed_params, c(
      beta0 = variable_params[i, "beta0est"]
    , soc_dist_level_sip = variable_params[i, "soc_dist_level_sip"]
    , Ca    = variable_params[i, ]$Ca
    , alpha = variable_params[i, ]$alpha
    , delta = variable_params[i, ]$delta
    , mu    = variable_params[i, ]$mu
      )
      , if (fit.E0) { 
          c(E_init = variable_params[i, ]$E_init)
        } else {
          c(E0     = variable_params[i, ]$E0)     
        }
        )
      } else {
      c(fixed_params, c(
      beta0    = variable_params[i, "beta0est"]
    , soc_dist_level_sip = variable_params[i, "soc_dist_level_sip"]
      , Ca       = variable_params[i, ]$Ca
      , alpha    = variable_params[i, ]$alpha
      , lambda_a = variable_params[i, ]$lambda_a
      , lambda_s = variable_params[i, ]$lambda_s
      , lambda_m = variable_params[i, ]$lambda_m 
      )
      , if (fit.E0) { 
          c(E_init = variable_params[i, ]$E_init)
        } else {
          c(E0     = variable_params[i, ]$E0)     
        }
        ) 
      }}
    , nsim         = nsim
    , format       = "d"
    , include.data = F
    , seed         = 1001)) %>% {
      rbind(.,
         group_by(., day) %>%
           dplyr::select(-.id) %>%
           summarise_all(median) %>%
                    mutate(.id = "median"))
    }

## summarize epidemic
{
SEIR.sim.s  <- SEIR.sim  %>% 
  dplyr::group_by(.id) %>% 
  dplyr::summarize(total_H = sum(H))

SEIR.sim    <- left_join(SEIR.sim, SEIR.sim.s, by = ".id")

SEIR.sim    <- SEIR.sim %>% 
  dplyr::filter(
    total_H > 10
  ) %>% droplevels()

## Calc the summary statistics
SEIR.sim.ss <- SEIR.sim %>% 
  dplyr::filter(.id != "median") %>%
  dplyr::mutate(week   = day %/% 7) %>%
  dplyr::group_by(week, .id) %>%
  dplyr::summarize(
    ## Mean in hospitalizations by week
    sum_H = sum(H_new)
    ) %>%
  dplyr::group_by(.id) %>%  
  dplyr::mutate(
    ## Difference in hospitalizations at the level of the week
    diff_H = c(0, diff(sum_H))
    ) %>%
  dplyr::group_by(.id)

SEIR.sim.ss.1 <- SEIR.sim.ss %>% 
  dplyr::summarize(
    ## Maximum hospitalizations reached, summarized at the level of the week
    when_max_H = week[which.max(sum_H)]
    ## How many hospitalizations are reached in that week
  , max_H      = max(sum_H)
    ## First week we see a reduction in the number of hospitalizations from a runs _global_ rate peak
  , when_red_H = week[min(which(diff_H == max(diff_H)))]
    )

## Number S to calc R_eff
SEIR.sim.ss.2 <- SEIR.sim %>% 
  dplyr::filter(day == max(county.data$day)) %>%
  dplyr::select(.id, S, D)

names(SEIR.sim.ss.2)[c(2, 3)] <- c("S_now", "D_now")

SEIR.sim.ss.f <- left_join(SEIR.sim.ss.1, SEIR.sim.ss.2, by = ".id")

SEIR.sim.ss.3 <- SEIR.sim %>%
  filter(.id != "median") %>%
  dplyr::group_by(.id) %>% 
  dplyr::summarize(
    total_D = max(D)
  , total_R = max(R)
      )

SEIR.sim.ss.t   <- left_join(SEIR.sim.ss.f, SEIR.sim.ss.3, by = ".id")

## Get their CI
SEIR.sim.ss.t.s <- SEIR.sim.ss.t %>%
  pivot_longer(cols = -.id) %>%
  dplyr::group_by(name) %>%
  dplyr::summarize(
    lwr = quantile(value, c(0.025))
  , est = quantile(value, c(0.500))
  , upr = quantile(value, c(0.975))
  ) %>% 
  dplyr::mutate(paramset = variable_params[i, ]$paramset)

SEIR.sim.ss.t.s <- rbind(SEIR.sim.ss.t.s
  , data.frame(
      name = "prop_ext"
    , lwr  = nsim - nrow(SEIR.sim.ss.t)
    , est  = nsim - nrow(SEIR.sim.ss.t)
    , upr  = nsim - nrow(SEIR.sim.ss.t)
    , paramset = variable_params[i, ]$paramset))

SEIR.sim.ss.t.ci <- rbind(SEIR.sim.ss.t.ci, SEIR.sim.ss.t.s)
}

## Update variable params with R0 estimate
variable_params[i, ]$Reff <- with(variable_params[i, ], covid_R0(
  beta0est = beta0est, fixed_params = c(fixed_params, unlist(variable_params[i, ]))
  , sd_strength = soc_dist_level_sip
, prop_S = unlist(SEIR.sim.ss.t.s[SEIR.sim.ss.t.s$name == "S_now", 3]) / 
    (location_params[location_params$Parameter == "N", ]$est - 
        SEIR.sim.ss.t.s[SEIR.sim.ss.t.s$name == "total_D", 3])))

variable_params[i, ]$R0 <- with(variable_params[i, ], covid_R0(
  beta0est = beta0est, fixed_params = c(fixed_params, unlist(variable_params[i, ]))
  , sd_strength = 1, prop_S = 1))

if (((i / 20) %% 1) == 0) {
 saveRDS(
   list(
    variable_params  = variable_params
  , fixed_params     = fixed_params
  , dynamics_summary = SEIR.sim.ss.t.ci
  , param_array      = param_array
   ), paste(
     paste("output/"
       , paste(focal.county, fit_to_sip, more.params.uncer, fit.minus, Sys.Date(), "temp", sep = "_")
         , sep = "")
     , "Rds", sep = "."))
}

}

saveRDS(
   list(
    variable_params  = variable_params
  , fixed_params     = fixed_params
  , dynamics_summary = SEIR.sim.ss.t.ci
  , param_array      = param_array
   ), paste(
     paste("output/"
       , paste(focal.county, fit_to_sip, more.params.uncer, fit.minus, Sys.Date(), "final", sep = "_")
         , sep = "")
     , "Rds", sep = "."))
 