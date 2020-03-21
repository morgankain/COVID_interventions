## Santa Clara County data 
# to read in these data, should set working directory to main COVID-interventions directory
source("covid-19_setup.R")

## Data
SCC = read.csv("./SantaClara_CumCases_20200317.csv", stringsAsFactors = F) %>% 
  select(-X, -X.1) %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y"))

##' @param beta0        without intervention beta for all categories   
##' @param Ca           category specific contact rate
##' @param Cp           category specific contact rate
##' @param Cs           category specific contact rate
##' @param Cm           category specific contact rate
##' @param alpha        fraction of cases asymptomatic
##' @param gamma        1/time in exposed class
##' @param lambda_a     1/time for asympomatic to recover
##' @param lambda_s     1/time for severely symptomatic to go to the hospital
##' @param lambda_m     1/time for minorly sympomatic to recover
##' @param lambda_h     1/time  to leaving hospital  
##' @param delta        fraction of hospitalized cases that are fatal
##' @param mu           fraction of cases that are minor
##' @param rho          1/time in pre-symptomatic 
##' @param N            Population size
##' @param E0           initially exposed
##' @param intervention       1 is for social distancing, 2 is for threshhold based, currently thresh_H is based on
##' @param int_level          proportional reduction in contact rate
##' @param thresh_H_start     currently thresholding on total people in the hospital to start intervention
##' @param thresh_H_end       currently thresholding on total people in the hospital to end intervention
##' @param thresh_int_level   multiplier on beta when the thresshold causes the intervention to kick in

# Simulation parameters
covid_params <- c(
  beta0            = 0.5
, Ca               = 1
, Cp               = 1
, Cs               = 1
, Cm               = 1
, alpha            = 1/3
, gamma            = 1/5.2
, lambda_a         = 1/7
, lambda_s         = 1/4
, lambda_m         = 1/7
, lambda_p         = 1/0.5
, rho              = 1/10.7
, delta            = 0.2
, mu               = 19/20
, N                = 1937570 # (Santa Clara County) 59.02e6 (Wuhan)
, E0               = 20
#, intervention     = 2
#, thresh_H_start   = 15
#, thresh_H_end     = 2
#, thresh_int_level = 0.1
)

# simulate with a set of parameters
sim <- do.call(
  pomp::simulate
, list(
  object       = covid
, params       = covid_params
, nsim         = 50
, format       = "d"
, include.data = F
, seed         = 1001)
  ) %>% {rbind(.,
         group_by(., day) %>%
           select(-.id) %>%
           summarise_all(median) %>%
                    mutate(.id = "median"))} 

source("plotting.R")
