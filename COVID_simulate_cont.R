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
# source("COVID_pomp_gammabeta.R", local = T)
# all we need from this is the R0 calculation, and we don't want to read in anything else, so lets just grab those lines
read_start = grep("covid_R0", readLines("COVID_pomp_gammabeta.R"))
read_end = grep("}", readLines("COVID_pomp_gammabeta.R"))
read_end = read_end[which(read_end > read_start)] %>% min
source(textConnection(readLines("COVID_pomp_gammabeta.R")[read_start:read_end]))

## Load the previously saved fits
 ## If COVID_fit_cont.R was just run, use parameters already stored in the global env 
if (use.rds) {
  print(rds.name)
  prev.fit         <- readRDS(rds.name)
  covid_mobility   <- prev.fit[["covid_mobility"]] # get pomp object
  if("mifs_local_v2" %in% names(prev.fit)){ # get second round of mifs if it exists
    variable_params  <- prev.fit[["mifs_local_v2"]] # otherwise get the first round
  } else{
    variable_params  <- prev.fit[["mifs_local"]] # otherwise get the first round
  }
  
  pomp_data        <- data.frame(covid_mobility) %>% 
    full_join(as.data.frame(cbind(day = covid_mobility@covar@times, 
                                  t(covid_mobility@covar@table)))) %>%
    arrange(day) %>% 
    distinct_at(vars(day), .keep_all = T) # sometimes the pomp data and covariate data disagree? not sure why..., for now default to takign the pomp data
  dt               <- covid_mobility@rprocess@delta.t
  date_origin      <- prev.fit[["date_origin"]][1]
}
  
## drop the rows that have 0s for likelihood (in case exited prematurely) 
 ## and keep only the best fits as defined by loglik
if (loglik.max) {
  
variable_params <- variable_params %>% 
  filter(log_lik != 0) %>% 
  filter(log_lik == max(log_lik))

# print(variable_params$paramset)
print(paste0("subsetting to top parameter set"))
print(variable_params$log_lik)

} else {
  
  if (is.na(loglik.num)) {
    
variable_params <- variable_params %>% 
filter(log_lik != 0) %>% 
filter(log_lik > (max(log_lik) - loglik.thresh)) 

print(paste0("subsetting to ", nrow(variable_params), " parameter sets by likelihood threshold"))
print(variable_params$log_lik)

  } else {
    
variable_params <- variable_params %>% 
filter(log_lik != 0) %>% 
arrange(desc(log_lik)) %>%
slice(1:loglik.num)

print(paste0("subsetting to ", nrow(variable_params), " parameter sets"))
print(variable_params$log_lik)

}
}

if (!params.all) {
  variable_params <- variable_params[sample(1:nrow(variable_params), nparams), ]
}


####
## Set up covariate table for counterfactual scenarios
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

  ## Covariate table with all of the intervention scenarios
    int.init = as.Date(int.init)
    sim_end = as.Date(sim_end)
    
    # some catches for problematic interventions
    if(min(int.init) < (max(pomp_data$day) + date_origin)){ # catch if interventions start before the data end
      stop("first intervention must start after data ends")
    } else if(is.unsorted(int.init)){
      stop("interventions must be specified in chronological order")
    } else if(length(int.init) != length(int.type)){
      stop("int.type and int.init must have same length")
    } else if(length(int.init) + 1 != length(int.movement)){
      stop("length(int.init) != length(int.movement[-1])\nfirst element of int.movement specifies movement type after data ends and before interventions begin")
    }
    
    keep_ints <- int.init < sim_end
    # warning for dropped interventions
    if(sum(keep_ints) < length(int.init)){ 
      print(paste("dropping last ", length(int.init) - sum(keep_ints), " interventions that begin after sim_end"))
    }
    int.init = int.init[keep_ints] # trim to interventions before the sim_end date
    int.type = int.type[keep_ints] # apply same trim to intervention types and movenent
    int.movement = int.movement[c(TRUE, keep_ints)] # need to add a T for the movement that occurs between data and first intervention
    
    print(int.type)
    
    int.phases <- c( 
      nrow(pomp_data), # first phase is the data
      # second phase is until first intervention, rest are between interventions, last is from last intervention to sim_end date
      c(int.init, sim_end) - c(max(pomp_data$day) + date_origin, int.init) 
    )
    
    post_mob <- mean(tail(pomp_data$sip_prop, 3))
    pre_mob <- mean(head(pomp_data$sip_prop, 7))
    mid_mob <- (pre_mob + post_mob)/2
    
    
    mob.covtab <- covariate_table(
      sip_prop           = c(pomp_data$sip_prop, 
                             rep(as.numeric(mapvalues(int.movement, 
                                                      from = c("pre", "mid", "post"),
                                                      to = c(pre_mob, mid_mob, post_mob))), 
                                 int.phases[-1]))       
      , order            = "constant"
      , times            = seq(min(pomp_data$day), as.numeric(sim_end - date_origin), by = 1)
      , iso_mild_level   = rep((c("none", "none", int.type) == "inf_iso")*iso_mild_level, int.phases)
      , iso_severe_level = rep((c("none", "none", int.type) == "inf_iso")*iso_severe_level, int.phases)
      , intervention      = rep(as.numeric(mapvalues(c("none", "none", int.type),
                                         from = c("none", "tail", "inf_iso"),
                                         to   = c(0, 1, 2))), int.phases)   
      )
}

# use it to make a new pomp object
covid_mobility_sim <- pomp(covid_mobility, covar = mob.covtab)

####
## Simulate: loop over rows of variable_params
####

Reff   <- data.frame(paramset = 0, date = 0, Reff = 0)  ; Reff <- Reff[-1, ]
betat  <- data.frame(paramset = 0, date = 0, betat = 0) ; betat <- betat[-1, ]
detect <- data.frame(paramset = 0, date = 0, detect = 0); detect <- detect[-1, ]

for (i in 1:nrow(variable_params)) {

  # calculate beta_catch value
  if(int.beta_catch_type == "pct") {
    beta_catch_val = max(qgamma(p = (1 - int.beta_catch), 
                                shape = int.beta0_k*dt/variable_params[i, "d"], 
                                scale = variable_params[i, "beta0"]/int.beta0_k),
                         0.001)
  } else{
    beta_catch_val = max(int.beta_catch, 0.001) # crude way to ensure int.beta_catch > 0, note that even this nonzero low beta_catch will be problematic but will eventually finish
  }
  
checktime <- system.time({

sim_times <- seq(as.numeric(variable_params[i,]$sim_start - date_origin), max(mob.covtab@times), 1) 

## Simulate from fits
    SEIR.sim <- do.call(
      pomp::simulate
      , list(
        object   = covid_mobility_sim
        , t0     = as.numeric(variable_params[i,]$sim_start - date_origin)
        , times  = sim_times
        , params = variable_params[i,] %>% unlist %>% 
                       inset("beta_catch", value = beta_catch_val) %>% 
                       inset("beta0_k", value = int.beta0_k) %>% 
                       inset("catch_eff", value = int.catch_eff)
        , nsim         = nsim
        , format       = "d"
        , include.data = F
        , seed         = 1001))
  
})[3]

print(checktime)

SEIR.sim <- SEIR.sim %>%
  mutate(
      date     = day + date_origin
    , paramset = i)

if (ci.epidemic) {
  epi_ids <- SEIR.sim %>% 
    group_by(.id) %>% 
    summarise(total_infect = max(D + R), .groups = "drop") %>% 
    filter(total_infect > ci.epidemic_cut*ceiling(variable_params[i, "E_init"])) %>% 
    filter(total_infect > ci.epidemic_cut) %>% #*ceiling(variable_params[i, "E_init"])) %>% 
    pull(.id)
  print(paste0("limiting to epidemics, including ", length(epi_ids), " simulations"))
  SEIR.sim %<>% filter(.id %in% epi_ids)
}

SEIR.sim %<>% {
  rbind(.,
        group_by(., day) %>%
          dplyr::select(-.id) %>%
          summarise_all(median, .groups = "drop") %>%
          mutate(.id = "median"))
}

SEIR.sim.s <- SEIR.sim %>% 
  dplyr::filter(.id != "median") %>%
  group_by(day) %>%
  summarise(S = mean(S), D = mean(D), .groups = "drop") %>%
  mutate(date = day + date_origin)

betat.t <- data.frame(
  paramset = i
  , date     = mob.covtab@times + date_origin
  , betat    = variable_params[i, "beta0"] * exp(log(variable_params[i, "beta_min"])*mob.covtab@table[which(dimnames(mob.covtab@table)[[1]] == "sip_prop"), ])
)

# not the cleanest, but join the sim summary and betat to get the dates aligned first, then calculate Reff
Reff.t <- full_join(SEIR.sim.s, betat.t) %>% arrange(day) %>% 
  filter(!is.na(S), !is.na(D), !is.na(betat)) %>% 
  {data.frame(
      paramset = i
      , date     = pull(., date)
      , Reff     = covid_R0(
        beta0est      = pull(., betat)
        , fixed_params  = unlist(variable_params[i, ])
        , sd_strength   = 1
        , prop_S        = pull(.,S) / (variable_params[i, "N"] - pull(.,D))
      )
    )
    }

# Reff.t <- data.frame(
#   paramset = i
#   , date     = mob.covtab@times + date_origin
#   , Reff     = covid_R0(
#     beta0est      = betat.t$betat
#     , fixed_params  = unlist(variable_params[i, ])
#     , sd_strength   = 1
#     , prop_S        = SEIR.sim.s$S / (variable_params[i, "N"] - SEIR.sim.s$D)
#   )
# )

# detect.t <- c(rep(0, min(county.data[!is.na(county.data$cases), ]$day) - 1)
#               , (variable_params[i, "detect_max"] /
#                    (1 + exp(-variable_params[i, "detect_k"] * (mob.covtab@times - variable_params[i, "detect_mid"])))
#               )[-seq(1, min(county.data[!is.na(county.data$cases), ]$day) - 1)]
# )
detect.t <- data.frame(
  paramset = i
  , date     = mob.covtab@times + date_origin
  , detect   = c(rep(0, min(pomp_data[!is.na(pomp_data$cases), "day"]) - 1)
                 , (variable_params[i, "detect_max"] /
                      (1 + exp(-variable_params[i, "detect_k"] * (mob.covtab@times - variable_params[i, "detect_mid"])))
                 )[-seq(1, min(pomp_data[!is.na(pomp_data$cases), ]$day) - 1)]))

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

####
## Summary for plotting
####

# trajectories
SEIR.sim.f.t <- SEIR.sim.f %>% 
  dplyr::select(date, day, .id, paramset, any_of(plot_vars))  %>%
  pivot_longer(any_of(plot_vars))

# confidence intervals
SEIR.sim.f.ci <- SEIR.sim.f.t %>% 
  group_by(date, paramset, name) %>% 
  summarise(
      lwr = quantile(value, 0 + ci.stoc)
    , mid = {
      if (plot.median) {
      quantile(value, c(0.500))  
      } else {
      mean(value) 
      }
    }
    , upr = quantile(value, 1 - ci.stoc),
    .groups = "drop"
  ) 

print("completed simulations")
