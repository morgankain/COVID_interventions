##################################################################
## Use pomp fits to simulate dynamics for various interventions ##
##################################################################
## Using continuous human movement data, fitting to cases and deaths and fitting case detection 

## Note: This script is set up to be able to be run from its own from a saved .Rds.
 ## Thus, packages are set to be loaded etc. (these lines do not need to be run if COVID_fit_cont.R was just run)

## For detailed conceivable intervention scenarios and details on how tosimulate them see:
## "potential_intervention_details_and_plots.R"

options(warning.length = 4000L)

set.seed(10001)
# source("needed_packages.R")

## Be very careful here, adjust according to your machine's capabilities
registerDoParallel(cores = usable.cores)
  
## Bring in pomp objects.
source("COVID_pomp_gammabeta.R", local = T)

  ## One special case
  if (focal.county == "Fulton") {
    
deaths <- read.csv("us-counties.txt") %>%
  mutate(date = as.Date(date)) %>%
  dplyr::filter((county == focal.county | county == "DeKalb") & state == focal.state) %>%
  group_by(date) %>% summarize(cases = sum(cases), deaths = sum(deaths)) %>%
  mutate(county = focal.county, state = focal.state)

  } else {
    
deaths  <- read.csv("us-counties.txt") %>% 
  mutate(date = as.Date(date)) %>% 
  dplyr::filter(county == focal.county) %>% 
  dplyr::filter(date < max(date) - fit.minus) %>% 
  filter(state == focal.state)
    
  }

## Load the previously saved fits
 ## If COVID_fit_cont.R was just run, use parameers already stored in the global env 
if (use.rds) {
print(rds.name)
prev.fit         <- readRDS(rds.name)
variable_params  <- prev.fit[["variable_params"]]
fixed_params     <- prev.fit[["fixed_params"]]
}
  
## drop the rows that have 0s for likelihood (in case exited prematurely) 
 ## and keep only the best fits as defined by loglik
if (loglik.max) {
  
variable_params <- variable_params %>% 
  filter(log_lik != 0) %>% 
  filter(log_lik == max(log_lik))

print(variable_params$paramset)
print(variable_params$log_lik)

} else {
  
  if (is.na(loglik.num)) {
    
variable_params <- variable_params %>% 
filter(log_lik != 0) %>% 
filter(log_lik > (max(log_lik) - loglik.thresh)) 

print(variable_params$paramset)
print(variable_params$log_lik) 

  } else {
    
variable_params <- variable_params %>% 
filter(log_lik != 0) %>% 
arrange(desc(log_lik)) %>%
slice(1:loglik.num)

print(variable_params$paramset)
print(variable_params$log_lik)

}
}

if (!params.all) {
  variable_params <- variable_params[sample(1:nrow(variable_params), nparams), ]
}

## Second bit of the intervention will affect these parameters:
fixed_params["beta0_k"]   <- int.beta0_k      ## heterogeneity value
fixed_params["catch_eff"] <- int.catch_eff   ## beta0 values caught by intervention

if (focal.county == "Fulton") {
  mobility <- readRDS(mobility.file) %>% 
    dplyr::filter((county_name == focal.county | county_name == "DeKalb") & (state_abbr == focal.state_abbr)) %>%
    dplyr::group_by(datestr) %>%
    dplyr::summarize(sip_prop = mean(sip_prop), .groups = "drop") 
} else {
  mobility <- readRDS(mobility.file) %>% 
    dplyr::filter(county_name == focal.county & state_abbr == focal.state_abbr)
}

# clean up mobility data by taking 3 day moving median and adding days since date_origin
mobility %<>% 
  dplyr::select(datestr, sip_prop) %>% 
  filter(sip_prop != 0) %>% 
  {full_join(., # join with full list of dates so NAs will occur for missing days
             data.frame(datestr = seq(pull(., datestr) %>% min, 
                                      pull(., datestr) %>% max, 
                                      by = "day")))} %>% 
  mutate(day = as.numeric(datestr - date_origin)) %>% 
  arrange(day) %>% 
  # three day moving median
  mutate(sip_prop = mutate(., 
                           sip_prop_lag = lag(sip_prop),
                           sip_prop_lead = lead(sip_prop)) %>% 
           select(contains("sip_prop")) %>%
           purrr::pmap_dbl(~median(c(...), na.rm = T))) %>% 
  dplyr::filter(!is.na(sip_prop)) 

# backfill mobility to date_origin
if (min(mobility$day) > 1) {
  mobility <- rbind(
    data.frame(
      datestr  = seq.Date(date_origin + 1, 
                          min(mobility$datestr) - 1,
                          "day")
      , sip_prop = mean(mobility$sip_prop[1:10])
    ) %>% mutate(day = as.numeric(datestr - date_origin)) 
    , mobility
  )  
}

## Adjust the cases/deaths data to anchor dates from date_origin
county.data <- deaths %>% 
  mutate(day = as.numeric(date - date_origin)) %>% 
  dplyr::select(day, date, deaths, cases) %>% 
  arrange(date) %>%
  # convert to daily deaths and cases
  mutate(deaths = deaths - lag(deaths),
         cases = cases - lag(cases)) %>%
  # convert negatve cases and deaths to zeros
  dplyr::rowwise() %>%
  mutate(cases = max(0, cases),
         deaths = max(0, deaths)) %>%
  ungroup() %>%
  dplyr::filter(!is.na(deaths), !is.na(cases)) %>%
  # fill with NAs from the latest considered sim start date to the data, 
  # NAs are to prevent any oddities with the accumulator variables
  # if the last sim_start date considered occurs when or after the data start, problems will still occur
  {if(max(variable_params[, "sim_start"]) < min(pull(., "date"))) {rbind(
    data.frame(
      date   = seq.Date(variable_params[, "sim_start"] %>% max, 
                        min(pull(., "date")) - 1,
                        "day")
      , deaths = NA
      , cases  = NA
    ) %>% 
      mutate(day = as.numeric(date - date_origin))
    , .
  )} else{ # if the data start before the last sim start date, trim the data to the last sim start and add one row of NAs to deal with accumulator variable
    filter(., date >=  variable_params[, "sim_start"] %>% max) %>% 
      {rbind(data.frame(date = min(pull(., "date")) - 1,
                       deaths = NA, 
                       cases = NA) %>% 
               mutate(day = as.numeric(date - date_origin)),
             .)}
  ## Remove dates after one week after movement data ends
  }} %>% dplyr::filter(date < (max(mobility$datestr) + 7))
  
####
## Simulate: loop over rows of variable_params
####

Reff   <- data.frame(paramset = 0, date = 0, Reff = 0)  ; Reff <- Reff[-1, ]
betat  <- data.frame(paramset = 0, date = 0, betat = 0) ; betat <- betat[-1, ]
detect <- data.frame(paramset = 0, date = 0, detect = 0); detect <- detect[-1, ]

for (i in 1:nrow(variable_params)) {
  
fixed_params["beta_catch"]  <- {                ## beta0 values caught by intervention
  if(int.beta_catch_type == "pct") {
    max(qgamma(p = (1 - int.beta_catch), 
               shape = int.beta0_k, 
               scale = variable_params[i, "beta0est"]/int.beta0_k),
        0.001)
  } else{
    max(int.beta_catch, 0.001) # crude way to ensure int.beta_catch > 0, note that even this nonzero low beta_catch will be problematic but will eventually finish
  }
}
  
####
## Set up counterfactual scenarios
####
if (counter.factual) {
  
if (cf.type == "no_int") {
  
mob.covtab <- covariate_table(
   sip_prop         = rep(mean(head(mobility$sip_prop, 7)), sim_length)
 , order            = "constant"
 , times            = seq(1, sim_length, by = 1)
 , detect_t0        = min(county.data[!is.na(county.data$cases), ]$day) - 1
 , iso_mild_level   = NA
 , iso_severe_level = NA
 , intervention     = rep(0, sim_length))  

} else if (cf.type == "delay") {
  
mob.covtab <- covariate_table(
   sip_prop         = c(rep(mean(head(mobility$sip_prop, 7)), delay.time), mobility$sip_prop
     , rep(mean(tail(mobility$sip_prop, 3)), sim_length - (length(mobility$sip_prop) + delay.time))) 
 , order            = "constant"
 , times            = seq(1, sim_length, by = 1)
 , detect_t0        = min(county.data[!is.na(county.data$cases), ]$day) - 1
 , iso_mild_level   = NA
 , iso_severe_level = NA
 , intervention     = rep(0, sim_length))  

} else if (cf.type == "may1") { 
  
int.begin    <- as.numeric((as.Date("2020-05-01") - variable_params[i, ]$sim_start))
int.phase1   <- nrow(mobility) + (int.begin - nrow(mobility))
  
mob.covtab <- covariate_table(
   sip_prop         = c(
      mobility$sip_prop[1:int.begin]
    , rep(mean(head(mobility$sip_prop, 7)), sim_length - int.phase1)
      )
 , order            = "constant"
 , times            = seq(1, sim_length, by = 1)
 , detect_t0        = min(county.data[!is.na(county.data$cases), ]$day) - 1
 , iso_mild_level   = NA
 , iso_severe_level = NA
 , intervention     = rep(0, sim_length))  
}
  
} else {
  
####
## Set up intervention scenario
####

int.begin    <- as.numeric((as.Date(int.init) - date_origin))
int.duration <- min(as.Date(int.end), variable_params[i, ]$sim_start + sim_length) - as.Date(int.init)

int.phase1 <- nrow(mobility) + (int.begin - nrow(mobility))
int.phase2 <- as.numeric(int.duration)
int.phase3 <- sim_length - int.phase1 - int.phase2

## Covariate table with all of the intervention scenarios
{
mob.covtab <- covariate_table(
  ## First, can add whatever sip_prop into the future
  
   sip_prop      = {
     if (int.movement == "pre") {
     c(
      mobility$sip_prop
    , rep(mean(tail(mobility$sip_prop, 3)), int.phase1 - length(mobility$sip_prop))
    , rep(mean(head(mobility$sip_prop, 7)), sim_length - int.phase1)
      )       
     } else if (int.movement == "post") {
     c(
      mobility$sip_prop
    , rep(mean(tail(mobility$sip_prop, 3)), int.phase1 - length(mobility$sip_prop))
    , rep(mean(tail(mobility$sip_prop, 3)), sim_length - int.phase1)
      )
     } else if (int.movement == "mid") {
     c(
      mobility$sip_prop
    , rep(mean(tail(mobility$sip_prop, 3)), int.phase1 - length(mobility$sip_prop))
    , rep((mean(head(mobility$sip_prop, 7)) + mean(tail(mobility$sip_prop, 3))) / 2, sim_length - int.phase1)
      )         
     }
    
   }
  
 , order         = "constant"
 , times         = seq(1, sim_length, by = 1)

 , iso_mild_level   = {
   if (int.type == "inf_iso") {
   iso_mild_level 
   } else {
     NA
   }
 }
 , iso_severe_level = {
   if (int.type == "inf_iso") {
   iso_severe_level
   } else {
     NA
   }
 }
  
  ## Next can add an intervention that stops large gatherings
  
 , intervention  = {
   if (int.type == "none") {
   c(
    rep(0, int.phase1)
  , rep(0, int.phase2)
  , rep(0, int.phase3)
  )     
   } else if (int.type == "tail") {
     
     if (int.continue) {
   c(
    rep(0, int.phase1)
  , rep(1, int.phase2)
  , rep(1, int.phase3)
   )   
     } else {
   c(
    rep(0, int.phase1)
  , rep(1, int.phase2)
  , rep(0, int.phase3)
   )   
     }
     
   } else if (int.type == "inf_iso") {
     
     if (int.continue) {
   c(
    rep(0, int.phase1)
  , rep(2, int.phase2)
  , rep(2, int.phase3)
   )   
     } else {
   c(
    rep(0, int.phase1)
  , rep(2, int.phase2)
  , rep(0, int.phase3)
   )   
     }
     
  }
   }
    )
}

}

covid_mobility <- pomp(
   data       = county.data %>% select(day, cases, deaths)
 , times      = "day"
 , t0         = 1
 , covar      = mob.covtab
 , rprocess   = euler(sir_step_mobility, delta.t = dt)
 , rmeasure   = {if(con_theta){rmeas_multi_logis_con}else{rmeas_multi_logis_ind}}
 , dmeasure   = {if(con_theta){dmeas_multi_logis_con}else{dmeas_multi_logis_ind}}
 , rinit      = sir_init
 , partrans   = {if(con_theta){par_trans_con}else{par_trans_ind}}
 , accumvars  = accum_names
 , paramnames = param_names
 , statenames = state_names
 , globals    = globs
)

print(int.type)

checktime <- system.time({

## Simulate from fits
SEIR.sim <- do.call(
  pomp::simulate
  , list(
    object   = covid_mobility
    , t0     = as.numeric(sim_parms$sim_start - date_origin)
  # , times  = 1:sim_length
    , times  = c(county.data$day, max(county.data$day):(max(county.data$day) + sim_length - length(county.data$day)))
    , params = c(fixed_params %>% t() %>% 
                   as.data.frame() %>% 
                   mutate(d  = variable_params[i, "alpha"] * lambda_a + 
                       (1 - variable_params[i, "alpha"]) * variable_params[i, "mu"] * (lambda_p + lambda_m) + 
                       (1 - variable_params[i, "alpha"])*(1 - variable_params[i, "mu"])*(lambda_p + lambda_s)
                   ) %>%  
                   # convert periods to rates
                   mutate(gamma    = -1/(nE*dt)*log(1-nE*dt/gamma),
                          lambda_a = -1/(nIa*dt)*log(1-nIa*dt/lambda_a),
                          lambda_p = -1/(nIp*dt)*log(1-nIp*dt/lambda_p), 
                          lambda_m = -1/(nIm*dt)*log(1-nIm*dt/lambda_m),
                          lambda_s = -1/(nIs*dt)*log(1-nIs*dt/lambda_s),
                          rho_r    = -1/dt*log(1-dt/rho_r)) %>%
                   unlist
      , c(
        beta0      = variable_params[i, "beta0est"]
      , E_init     = variable_params[i, "E_init"]
      , detect_max = variable_params[i, "detect_max"]
      , detect_mid = variable_params[i, "detect_mid"]
      , detect_k   = variable_params[i, "detect_k"]
      , theta      = variable_params[i, "theta"]
      , theta2     = variable_params[i, "theta2"]
      , beta_min   = variable_params[i, "beta_min"]
      , Ca         = variable_params[i, "Ca"]
      , alpha      = variable_params[i, "alpha"]
      , delta      = variable_params[i, "delta"]
      , mu         = variable_params[i, "mu"]
      , E_init     = variable_params[i, "E_init"]
      , rho_d      = {-1/dt*log(1-dt/variable_params[i, "rho_d"])}
      , detect_t0  = min(county.data[!is.na(county.data$cases), ]$day) - 1
        ))
    , nsim         = nsim
    , format       = "d"
    , include.data = F
    , seed         = 1001))
      
})[3]

print(checktime)

SEIR.sim <- SEIR.sim %>%
  mutate(
      date     = round(as.Date(day, origin = variable_params[i, ]$sim_start))
    , paramset = variable_params[i, ]$paramset)

if (ci.epidemic) {
  epi_ids <- SEIR.sim %>% 
    group_by(.id) %>% 
    summarise(total_infect = max(D + R)) %>% 
    filter(total_infect > ci.epidemic_cut*ceiling(variable_params[i, "E_init"])) %>% 
    pull(.id)
  print(paste0("limiting to epidemics, including ", length(epi_ids), " simulations"))
  SEIR.sim %<>% filter(.id %in% epi_ids)
}

SEIR.sim %<>% {
  rbind(.,
        group_by(., day) %>%
          dplyr::select(-.id) %>%
          summarise_all(median) %>%
          mutate(.id = "median"))
}

SEIR.sim.s <- SEIR.sim %>% 
  dplyr::filter(.id != "median") %>%
  group_by(day) %>%
  summarize(S = mean(S), D = mean(D))

betat.t      <- variable_params[i, "beta0est"] * exp(log(variable_params[i, "beta_min"])*mob.covtab@table[which(dimnames(mob.covtab@table)[[1]] == "sip_prop"), ])

Reff.t       <- with(variable_params[i, ], covid_R0(
   beta0est      = betat.t
 , fixed_params  = c(fixed_params, unlist(variable_params[i, ]))
 , sd_strength   = 1
 , prop_S        = SEIR.sim.s$S / (variable_params[i, "N"] - SEIR.sim.s$D)
  )
  )

detect.t <- c(rep(0, min(county.data[!is.na(county.data$cases), ]$day) - 1)
              , (variable_params[i, "detect_max"] / 
                   (1 + exp(-variable_params[i, "detect_k"] * (mob.covtab@times - variable_params[i, "detect_mid"])))
              )[-seq(1, min(county.data[!is.na(county.data$cases), ]$day) - 1)]
)

betat.t <- data.frame(
  paramset = variable_params[i, ]$paramset
  , date     = seq(variable_params[i, ]$sim_start
                   , variable_params[i, ]$sim_start + (length(mob.covtab@table[which(dimnames(mob.covtab@table)[[1]] == "sip_prop"), ]) - 1), by = 1)
  , betat    = betat.t
)

Reff.t <- data.frame(
  paramset = variable_params[i, ]$paramset
  , date     = seq(variable_params[i, ]$sim_start
                   , variable_params[i, ]$sim_start + (length(mob.covtab@table[which(dimnames(mob.covtab@table)[[1]] == "sip_prop"), ]) - 1), by = 1)
  , Reff     = Reff.t
)

detect.t <- data.frame(
  paramset = variable_params[i, ]$paramset
  , date     = seq(variable_params[i, ]$sim_start
                   , variable_params[i, ]$sim_start + (length(mob.covtab@table[which(dimnames(mob.covtab@table)[[1]] == "sip_prop"), ]) - 1), by = 1)
  , detect   = detect.t
)

betat  <- rbind(betat, betat.t)
Reff   <- rbind(Reff, Reff.t)
detect <- rbind(detect, detect.t)

## Stich together output
if (i == 1) {
  SEIR.sim.f <- SEIR.sim  
} else {
  SEIR.sim.f <- rbind(SEIR.sim.f, SEIR.sim) 
}  

## Keep track of progress
if (((i / 20) %% 1) == 0) {
  print(paste(round(i / nrow(variable_params), 2)*100, "% Complete", sep = ""))
}

}

## Because of the staggered start dates there is some oddity at the end of the simulation, just remove those few days
SEIR.sim.f <- SEIR.sim.f %>% filter(date < min(variable_params$sim_start + sim_length))

####
## Summary for plotting
####

SEIR.sim.f.t <- SEIR.sim.f %>% 
  dplyr::select(date, day, .id, paramset, any_of(plot_vars))  %>%
  pivot_longer(any_of(plot_vars))

SEIR.sim.f.ci <- SEIR.sim.f.t %>% 
  group_by(date, paramset, name) %>% 
  summarize(
      lwr = quantile(value, 0 + ci.stoc)
    , mid = {
      if (plot.median) {
      quantile(value, c(0.500))  
      } else {
      mean(value) 
      }
    }
    , upr = quantile(value, 1 - ci.stoc)
  ) 
