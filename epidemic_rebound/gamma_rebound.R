### Simulate trajectories under a gamma distributed beta distribution and compare it to the mean

source("epidemic_rebound/gamma_rebound_setup.R", local = T)

####
## Simulate: loop over rows of variable_params
####

Reff   <- data.frame(paramset = 0, date = 0, Reff = 0)  ; Reff <- Reff[-1, ]
betat  <- data.frame(paramset = 0, date = 0, betat = 0) ; betat <- betat[-1, ]
detect <- data.frame(paramset = 0, date = 0, detect = 0); detect <- detect[-1, ]

for (i in 1:nrow(variable_params)) {
  
county.data <- deaths %>% 
  mutate(day = as.numeric(date - variable_params[i, ]$sim_start)) %>% 
  dplyr::filter(day > 0) %>%
  dplyr::select(day, date, deaths, cases) %>% 
  mutate(deaths = deaths - lag(deaths),
         cases = cases - lag(cases)) %>% 
  dplyr::filter(!is.na(deaths), !is.na(cases))

county.data <- rbind(
  data.frame(
    day    = seq(1:(min(county.data$day) - 1))
  , date   = as.Date(seq(1:(min(county.data$day) - 1)), origin = variable_params[i, "sim_start"])
  , deaths = NA
  , cases  = NA
  )
, county.data 
  )

mobility <- readRDS("unfolded_Jun05.rds") %>% 
  dplyr::filter(county_name == focal.county & state_abbr == focal.state_abbr)  %>%
  dplyr::select(datestr, sip_prop) %>% 
  dplyr::filter(!is.na(sip_prop)) %>%
  mutate(day = as.numeric(datestr - as.Date("2019-12-31"))) %>% 
  arrange(day) %>% 
  dplyr::filter(datestr >= variable_params[i, ]$sim_start) %>%
  mutate(day = day - as.numeric((variable_params[i, ]$sim_start - as.Date("2019-12-31")))) %>%
  dplyr::filter(!is.na(sip_prop)) %>% 
  dplyr::filter(day > 0)
  
if (min(mobility$day) > 1) {
mobility <- rbind(
  data.frame(
    datestr  = as.Date(seq(1:(min(mobility$day) - 1)), origin = variable_params[i, ]$sim_start)
  , sip_prop = mean(mobility$sip_prop[1:10])
  , day      = seq(1:(min(mobility$day) - 1))
  )
, mobility
  )  
}
 
mobility <- mobility %>% filter(sip_prop != 0)

print(rds.name)
print(variable_params$paramset)
print(variable_params$log_lik)

## Remove dates after one week after movement data ends
county.data <- county.data %>% dplyr::filter(date < (max(mobility$datestr) + 7))

if (counter.factual) {
  
if (cf.type == "no_int") {
mob.covtab <- covariate_table(
   sip_prop         = rep(mean(head(mobility$sip_prop, 7)), sim_length)
 , order            = "constant"
 , times            = seq(1, sim_length, by = 1)
 , detect_t0        = min(county.data[!is.na(county.data$cases), ]$day) - 1
 , iso_mild_level   = iso_mild_level
 , iso_severe_level = iso_severe_level
 , intervention     = rep(0, sim_length))  
} else if (cf.type == "delay") {
mob.covtab <- covariate_table(
   sip_prop         = c(rep(mean(head(mobility$sip_prop, 7)), 7), mobility$sip_prop
     , rep(mean(tail(mobility$sip_prop, 3)), sim_length - (length(mobility$sip_prop) + 7))) 
 , order            = "constant"
 , times            = seq(1, sim_length, by = 1)
 , detect_t0        = min(county.data[!is.na(county.data$cases), ]$day) - 1
 , iso_mild_level   = iso_mild_level
 , iso_severe_level = iso_severe_level
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
 , iso_mild_level   = iso_mild_level
 , iso_severe_level = iso_severe_level
 , intervention     = rep(0, sim_length))  
}
  
} else {
  
####
## Set up intervention scenario
####

int.begin    <- as.numeric((as.Date(int.init) - variable_params[i, ]$sim_start))
int.duration <- as.Date(int.end) - as.Date(int.init)

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
 , detect_t0     = min(county.data[!is.na(county.data$cases), ]$day) - 1
  
 , thresh_inf       = thresh_inf.val
 , check_int        = c(rep(1, int.phase1), rep(0, int.phase2 + int.phase3))
 , iso_mild_level   = iso_mild_level
 , iso_severe_level = iso_severe_level
 , iso_mild_level_post   = iso_mild_level_post
 , iso_severe_level_post = iso_severe_level_post
  
  ## Next can add an intervention that stops large gatherings
  
 , intervention  = {
   if (int.type == "none") {
   c(
    rep(0, int.phase1)
  , rep(0, int.phase2)
  , rep(0, int.phase3)
  )     
   } else if (int.type == "tail") {
   c(
    rep(0, int.phase1)
  , rep(1, int.phase2)
  , rep(1, int.phase3)
   ) 
   } else if (int.type == "inf_iso") {
   c(
    rep(0, int.phase1)
  , rep(2, int.phase2)
  , rep(2, int.phase3)
   )   
  }
   }
    )
}

}

covid_mobility <- pomp(
   data       = county.data
 , times      = "day"
 , t0         = 1
 , covar      = mob.covtab
 , rprocess   = euler(sir_step_mobility, delta.t = 1/6)
 , rmeasure   = rmeas_multi_logis
 , dmeasure   = dmeas_multi_logis
 , rinit      = sir_init
 , partrans   = par_trans
 , accumvars  = accum_names
 , paramnames = param_names
 , statenames = state_names
 , globals    = trunc_gamma
)

## Second bit of the intervention will affect these parameters:
fixed_params["beta0_sigma"]      <- int.beta0_sigma  ## heterogeneity value
fixed_params["beta_catch"]       <- {                ## beta0 values caught by intervention
  if(int.beta_catch_type == "pct") {
    max(qgamma(p = (1 - int.beta_catch), 
               shape = int.beta0_sigma, 
               scale = variable_params[i, "beta0est"]/int.beta0_sigma),
        0.001)
  } else{
    max(int.beta_catch, 0.001) # crude way to ensure int.beta_catch > 0, note that even this nonzero low beta_catch will be problematic but will eventually finish
  }
}
fixed_params["catch_eff"]        <- int.catch_eff    ## beta0 values caught by intervention
fixed_params["sip_prop_scaling"] <- int.sip_prop_scaling

fixed_params["beta_catch_post"]  <- {                ## beta0 values caught by intervention
  if(int.beta_catch_type == "pct") {
    max(qgamma(p = (1 - int.beta_catch_post), 
               shape = int.beta0_sigma, 
               scale = variable_params[i, "beta0est"]/int.beta0_sigma),
        0.001)
  } else{
    max(int.beta_catch_post, 0.001) # crude way to ensure int.beta_catch > 0, note that even this nonzero low beta_catch will be problematic but will eventually finish
  }
}
fixed_params["catch_eff_post"]   <- int.catch_eff_post
fixed_params["beta0_sigma_post"] <- int.beta0_sigma_post  ## heterogeneity value

print(int.type)

checktime <- system.time({

SEIR.sim <- do.call(
  pomp::simulate
  , list(
    object         = covid_mobility
    , times        = mob.covtab@times
    , params = c(
      fixed_params
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
      , rho_d      = variable_params[i, "rho_d"]        
        ))
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

})[3]

print(checktime)
      
SEIR.sim <- SEIR.sim %>%
  mutate(
      date     = round(as.Date(day, origin = variable_params[i, ]$sim_start))
    , paramset = variable_params[i, ]$paramset)


## Need to think about the most principled way of simulating forward from the present. For now just take the median up until the present as the
 ## starting values, and project from there
if (sir_init.mid) {
  
new.startvals <- SEIR.sim %>% filter(.id == "median" & date == sir_init.mid.t)  
new.days      <- seq(as.numeric((as.Date(sir_init.mid.t) - variable_params[i, ]$sim_start)), sim_length)

mob.covtab.c <- mob.covtab
mob.covtab.c@times <- mob.covtab@times[1:(length(new.days))]
mob.covtab.c@table <- mob.covtab@table[, min(new.days):sim_length]

covid_mobility <- pomp(
   data       = county.data
 , times      = "day"
 , t0         = 1
 , covar      = mob.covtab.c
 , rprocess   = euler(sir_step_mobility, delta.t = 1/6)
 , rmeasure   = rmeas_multi_logis
 , dmeasure   = dmeas_multi_logis
 , rinit      = sir_init_mid
 , partrans   = par_trans
 , accumvars  = accum_names
 , paramnames = c(param_names, mid_init_param_names)
 , statenames = state_names
 , globals    = trunc_gamma
)

SEIR.sim.c <- do.call(
  pomp::simulate
  , list(
    object         = covid_mobility
    , times        = mob.covtab.c@times
    , params = c(prev.fit[["fixed_params"]]
      , c(
        beta0      = variable_params[i, "beta0est"]
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
        )
      , unlist(c(
        E_init     = round(new.startvals$E)
      , S0         = round(new.startvals$S)
      , Ia0        = round(new.startvals$Ia)
      , Ip0        = round(new.startvals$Ip)
      , Is0        = round(new.startvals$Is)
      , Im0        = round(new.startvals$Im)
      , Hr0        = round(new.startvals$Hr)
      , Hd0        = round(new.startvals$Hd)
      , R0         = round(new.startvals$R)
      , D0         = round(new.startvals$D)
      )))
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

SEIR.sim.c <- SEIR.sim.c %>%
  mutate(
      date     = round(as.Date(day, origin = sir_init.mid.t))
    , paramset = variable_params[i, ]$paramset)

# SEIR.sim.c %>% filter(date == as.Date(sir_init.mid.t) + 3, .id == "median")
# SEIR.sim %>% filter(date == as.Date(sir_init.mid.t) + 2, .id == "median")

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
 , prop_S        = SEIR.sim.s$S / (
   prev.fit[["fixed_params"]]["N"] -
     SEIR.sim.s$D
   )
 , inf_iso           = ifelse(int.type == "inf_iso", F, F)
 , iso_course_mild   = NA
 , iso_course_severe = NA
  )
 )

detect.t <- c(rep(0, mob.covtab@table[which(dimnames(mob.covtab@table)[[1]] == "detect_t0"), 1])
  , (variable_params[i, "detect_max"] / 
    (1 + exp(-variable_params[i, "detect_k"] * (mob.covtab@times - variable_params[i, "detect_mid"])))
    )[-seq(1, mob.covtab@table[which(dimnames(mob.covtab@table)[[1]] == "detect_t0"), 1])]
)

betat.t <- data.frame(
  paramset = variable_params[i, ]$paramset
, date     = seq(variable_params[i, ]$sim_start
, variable_params[i, ]$sim_start + (sim_length - 1), by = 1)
, betat    = betat.t
  )

Reff.t <- data.frame(
  paramset = variable_params[i, ]$paramset
, date     = seq(variable_params[i, ]$sim_start
, variable_params[i, ]$sim_start + (sim_length - 1), by = 1)
, Reff     = Reff.t
  )

detect.t <- data.frame(
  paramset = variable_params[i, ]$paramset
, date     = seq(variable_params[i, ]$sim_start
, variable_params[i, ]$sim_start + (sim_length - 1), by = 1)
, detect   = detect.t
  )

if (sir_init.mid) {
SEIR.sim.c <- left_join(SEIR.sim.c, detect.t[, c(2, 3)], by = "date")
SEIR.sim.c <- SEIR.sim.c %>% mutate(cases = I_new_sympt * detect)
}

betat  <- rbind(betat, betat.t)
Reff   <- rbind(Reff, Reff.t)
detect <- rbind(detect, detect.t)

## Stich together output
if (sir_init.mid) {
if (i == 1) {
 SEIR.sim.f <- SEIR.sim.c 
} else {
 SEIR.sim.f <- rbind(SEIR.sim.f, SEIR.sim.c) 
}
} else {
if (i == 1) {
 SEIR.sim.f <- SEIR.sim  
} else {
 SEIR.sim.f <- rbind(SEIR.sim.f, SEIR.sim) 
}  
}

## Keep track of progress
if (((i / 20) %% 1) == 0) {
  print(paste(round(i / nrow(variable_params), 2)*100, "% Complete", sep = ""))
}

}

SEIR.sim.f     <- SEIR.sim.f %>% filter(date < min(variable_params$sim_start + sim_length))

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

