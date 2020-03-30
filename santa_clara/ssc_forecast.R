############################################################################
## Use pomp to simulate dynamics up until time t and then project forward ##
############################################################################

## Search !! for next steps

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
)

lapply(needed_packages, require, character.only = TRUE)
source("./scc_pomp_objs.R")
## Be very careful here, adjust according to your machine
registerDoParallel(
  cores = 2
  )

####
## Step 1: Establish a reasonable parameter set
####

 ## parameters that wont vary
fixed_params <- c(
    Ca               = 2/3
  , Cp               = 1
  , Cs               = 1
  , Cm               = 1
  , alpha            = 1/3
  , gamma            = 1/5.2
  , lambda_a         = 1/7
  , lambda_s         = 1/5
  , lambda_m         = 1/5
  , lambda_p         = 1/2
  , rho              = 1/16 #10.7
  , delta            = 0.2
  , mu               = 1-0.044 #19/20
  , N                = 1.938e6
)

 ## parameters that will vary
variable_params <- sobolDesign(
  lower = c(
    E0          = 1
  , sim_start   = 5
  , int_start1  = 60
  , int_length2 = 30
  , sd_m1       = 0.5
  , sd_m2       = 0.1
  )
, upper = c(
    E0          = 10
  , sim_start   = 15
  , int_start1  = 69
  , int_length2 = 180
  , sd_m1       = 0.9
  , sd_m2       = 0.5
)
, nseq  = 500
) %>% mutate(
  sim_start   = round(sim_start)
, int_start1  = round(int_start1)
, int_length2 = round(int_length2)
, E0 = round(E0))

variable_params <- variable_params %>% mutate(
    sim_start  = as.Date(sim_start, origin = as.Date("2020-01-01"))
  , int_start1 = as.Date(round(int_start1), origin = as.Date("2020-01-01"))
  , int_start2 = as.Date("2020-03-17")
 ) %>% mutate(
    int_length1  = round(int_start2 - int_start1)
  ## Column for estimated beta0 (!! see lower down for desire to store likelyhood profile to define pmcmc runs)
  , beta0est     = 0
) 

## !! not implemented yet, to come later: Also add infected isolation after X days?

####
## Step 2: Build the pomp object from data and simulation parameters
####

## Run parameters
sim_start  <- variable_params$sim_start
sim_length <- 300
sim_end    <- sim_start + sim_length

for (i in 1:nrow(variable_params)) {
  
## US deaths data, pull out Santa Clara County
deaths     <- read.csv(
  "NYT-us-counties.csv"
, stringsAsFactors = F) %>% 
  mutate(date = as.Date(date))

scc_deaths <- deaths %>% 
  filter(county == "Santa Clara")  %>% 
  mutate(day = as.numeric(date - variable_params[i, ]$sim_start)) %>% 
  select(day, date, deaths) %>% 
  rename(deaths_cum = deaths) %>% 
  mutate(deaths = deaths_cum - lag(deaths_cum)) %>% 
  replace_na(list(deaths = 0)) %>%
  select(-deaths_cum)

## Add days from the start of the sim to the first recorded day in the dataset
scc_deaths <- rbind(
  data.frame(
    day    = seq(1:(min(scc_deaths$day) - 1))
  , date   = as.Date(seq(1:(min(scc_deaths$day) - 1)), origin = variable_params[i, ]$sim_start)
  , deaths = 0
  )
, scc_deaths
)

## Create intervention covariate table for the full forecast
intervention.forecast <- with(variable_params[i, ], {

 covariate_table(
  day              = 1:sim_length
  
, intervention     = c(
      # No intervention until intervention start time
    rep(0, int_start1 - sim_start)                   
      # Intervention style 1
  , rep(1, int_length1)
      # Intervention style 2
  , rep(1, int_length2)
      # Post intervention close
  , rep(0, sim_length - (int_start2 - sim_start) - int_length2)
  ) 
  
, isolation        = 0
, iso_severe_level = 0  # % of contats that severe cases maintain
, iso_mild_level   = 0.1  # % of contats that mild cases maintain
  
, soc_dist_level   = c(                       # intensity of the social distancing interventions
  rep(sd_m1, int_start2 - sim_start)
     ## slightly odd way to do this, but should work
, rep(sd_m2, sim_length - (int_start2 - sim_start))
  )
  
, thresh_H_start   = c(                       # starting threshhold on total # in hospital
  2
     ## slightly odd way to do this, but should work
, 2
)  
, thresh_H_end     = c(                       # ending threshhold on total # in hospital
  10
     ## slightly odd way to do this, but should work
, 10 
)
     ## No reason to have a second parameter here, just use the same val that the user picks for social distancing
, thresh_int_level = c(                       # level of social distancing implemented with the threshhold intervention
  rep(sd_m1, int_start2 - sim_start)
     ## slightly odd way to do this, but should work
, rep(sd_m2, sim_length - (int_start2 - sim_start))
)
  
, order            = "constant"
, times            = "day"
  )

})

covid.fitting <- scc_deaths %>%
  pomp(
    time       = "day"
  , t0         = 1
  , covar      = intervention.forecast
  , rprocess   = euler(sir_step, delta.t = 1/6)
  , rmeasure   = rmeas 
  , dmeasure   = dmeas
  , rinit      = sir_init
  , partrans   = par_trans
  , accumvars  = accum_names
  , paramnames = param_names
  , statenames = state_names
  ) 

####
## Step 3: Particle filter to get a fitted beta0
####

timeoneparam <- system.time({
mifs_local <- covid.fitting %>%
  mif2(
    t0     = 1
  , params = c(fixed_params, c(beta0 = 2.5/7, E0 = variable_params[i, ]$E0))
  , Np     = 3000
  , Nmif   = 50
  , cooling.fraction.50 = 0.5
  , rw.sd  = rw.sd(beta0 = 0.02)
        )
})


##### !!
## Not currently used, but want to fit pmcmc with uncertainty in beta0 following the likelihood profile for beta0
##### !!

## super ugly way of retreiving likelihood plot, come back to this when adding pmcmc
ggout <- mifs_local %>%
  traces() %>%
  melt() %>%
  ggplot(aes(x=iteration,y=value))+
  geom_line()+
  guides(color=FALSE)+
  facet_wrap(~variable,scales="free_y")+
  theme_bw()

mif.l  <- ggout$data %>% filter(variable == "loglik" | variable == "beta0") %>% droplevels()
mif.l  <- pivot_wider(mif.l, names_from = variable, values_from = value)
ggout2 <- ggplot(mif.l, aes(beta0, loglik)) + geom_point()

variable_params[i, "beta0est"] <- coef(mifs_local)["beta0"]


####
## Step 4: Simulate from the beginning and project forward with this beta0
## !! See above note: would prefer pmcmc for uncertainty in beta0. Next step.
####

## Simulate for projections
SEIR.sim <- do.call(
  pomp::simulate
  , list(
    object         = covid.fitting
    , times        = intervention.forecast@times
    , params       = c(fixed_params, c(beta0 = variable_params[i, "beta0est"], E0 = variable_params[i, ]$E0))
    , nsim         = 100
    , format       = "d"
    , include.data = F
    , seed         = 1001)) %>% {
      rbind(.,
         group_by(., day) %>%
           select(-.id) %>%
           summarise_all(median) %>%
                    mutate(.id = "median"))
      } 

ggout3 <- ggplot(SEIR.sim) +
  geom_line(aes(x      = day, 
                y      = deaths,
                group  = .id, 
                color  = .id == "data")) + 
  xlab("Date") + ylab("New fatalities") +
  guides(color = FALSE)+
  scale_color_manual(values = c("#D5D5D3", "#24281A")) + theme_bw()

####
## Step 5: Calculate summary statistics from these runs
####

## A) Rate of change calculated with hospitalizations
## B) Maximum number of hospitalizations at one point in time
## C) Timing of maximum number of hospitalizations

## Maybe somewhat controversial [?] choice here to remove all sims that don't take off as we
 ## know these didn't happen, so will really skew our summary values
SEIR.sim.s  <- SEIR.sim  %>% 
  group_by(.id) %>% 
  summarize(total_H = sum(H))

SEIR.sim    <- left_join(SEIR.sim, SEIR.sim.s, by = ".id")

SEIR.sim    <- SEIR.sim %>% 
  filter(
    total_H > 10
  ) %>% droplevels()

SEIR.sim.ss <- SEIR.sim %>% 
  filter(.id != "median") %>%
  mutate(week   = day %/% 7) %>%
  group_by(week, .id) %>%
  summarize(
    ## Mean in hospitalizations by week
    mean_H = mean(H)
    ) %>%
  group_by(.id) %>%  
  mutate(
    ## Difference in hospitalizations at the level of the week
    diff_H = c(0, diff(mean_H))
    ) %>%
  group_by(.id) %>% 
  summarize(
    ## Maximum hospitalizations reached, summarized at the level of the week
    when_max_H = week[which.max(mean_H)]
    ## How many hospitalizations are reached in that week
  , max_H      = max(mean_H)
    ## First week we see a reduction in the number of hospitalizations from a runs _global_ rate peak
  , when_red_H = week[min(which(diff_H == max(diff_H)))]
    )
  
SEIR.sim.ss2 <- SEIR.sim %>%

  group_by()

}
