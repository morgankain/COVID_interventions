needed_packages <- c(
    "pomp"
  , "plyr"
  , "dplyr"
  , "ggplot2"
  , "magrittr"
  , "scales"
  , "lubridate"
  , "tidyr"
)

lapply(needed_packages, require, character.only = TRUE)

# define the changes for a time step forward
sir_step <- Csnippet("
                     // adjust betat for social distancing interventions
                     double betat;
                     if(intervention == 2 & thresh_crossed == 1){ // 2 is for threshhold intervention
                       betat =  beta0*thresh_int_level; 
                     }
                     else if(intervention == 1) betat = beta0*soc_dist_level; // 1 is for social distancing
                     else betat = beta0; // everything else is no intervention
                     
                     // adjust contact rates if isolation of symptomatic cases is in place
                     double iso_m = 1;
                     double iso_s = 1;
                     if(isolation == 1){
                        iso_m = iso_mild_level;
                        iso_s = iso_severe_level;
                     }
                    
                     // calculate transition numbers
                     double dSE = rbinom(S, 1-exp(-betat*(Ca*Ia/N + Cp*Ip/N + iso_m*Cm*Im/N + iso_s*Cs*Is/N)*dt)); 
                     double rateE[2];
                     double dE_all[2];
                     rateE[0] = alpha*gamma; // going to asymtomatic
                     rateE[1] = (1-alpha)*gamma; // going to presymptomatic
                     reulermultinom(2, E, rateE, dt, &dE_all);
                     double dEIa = dE_all[0];
                     double dEIp = dE_all[1];
                     double dIaR = rbinom(Ia, 1 - exp(-lambda_a*dt));
                     double rateIp[2];
                     double dIp_all[2];
                     rateIp[0] = mu*lambda_p; // going to minor symptomatic
                     rateIp[1] = (1-mu)*lambda_p; // going to sever symptomatic
                     reulermultinom(2, Ip, rateIp, dt, &dIp_all);
                     double dIpIm = dIp_all[0];
                     double dIpIs = dIp_all[1];
                     double dIsH = rbinom(Is, 1 - exp(-lambda_s*dt));
                     double dImR = rbinom(Im, 1 - exp(-lambda_m*dt));
                     double rateH[2];
                     double dH_all[2];
                     rateH[0] = delta*rho;
                     rateH[1] = (1-delta)*rho;
                     reulermultinom(2, H, rateH, dt, &dH_all);
                     double dHD = dH_all[0];
                     double dHR = dH_all[1];
                     
                     // update the compartments
                     S  -= dSE; // susceptible 
                     E  += dSE - dEIa - dEIp; // exposed
                     Ia += dEIa - dIaR; // infectious and asymptomatic
                     Ip += dEIp - dIpIs - dIpIm; // infectious and pre-symptomatic
                     Is += dIpIs - dIsH; // infectious and severe symptoms (that will be hospitalized)
                     Im += dIpIm - dImR; // infectious and minor symptoms
                     H  += dIsH - dHD - dHR; // hospitalized
                     R  += dHR + dImR + dIaR; // recovered
                     D  += dHD; // fatalities
                     D_new += dHD; // daily fatalities
                     if(intervention == 2 & H >= thresh_H_start) thresh_crossed = 1;
                     else if(intervention == 2 & thresh_crossed == 1 & H < thresh_H_end) thresh_crossed = 0;
                     ")

# define the initial set up, currently, every is susceptible except the exposed people
sir_init <- Csnippet("
                     S = N-E0;
                     E = E0;
                     Ia = 0;
                     Ip = 0;
                     Is = 0;
                     Im = 0;
                     H = 0;
                     R = 0;
                     D = 0;
                     D_new = 0;
                     thresh_crossed = 0;
                     ")

# define random simulator of measurement
rmeas <- Csnippet("double tol = 1e-16;
                   deaths = rpois(D_new + tol);
                  ")
# define evaluation of model prob density function
dmeas <- Csnippet("double tol = 1e-16;
                   lik = dpois(deaths, D_new + tol, give_log);
                  ")

# US deaths data 
deaths <- read.csv("NYT-us-counties.csv",
                  stringsAsFactors = F) %>% 
  mutate(date = as.Date(date))

scc_deaths = deaths %>% 
  filter(county == "Santa Clara")  %>% 
  mutate(day = as.numeric(date - as.Date("2020-01-01"))) %>% select(day, deaths) %>% 
  rename(deaths_cum = deaths) %>% 
  mutate(deaths = deaths_cum - lag(deaths_cum)) %>% 
  replace_na(list(deaths = 0)) %>%
  select(-deaths_cum)

sim_length = nrow(scc_deaths)

intervention = covariate_table(
                               day = scc_deaths$day,
                               intervention = rep(0, sim_length),
                               isolation =  rep(0, sim_length),# 1 is isolation of symptomatic
                               iso_severe_level = rep(NA, sim_length), # % of contats that severe cases maintain
                               iso_mild_level = rep(NA, sim_length), # % of contats that mild cases maintain
                               soc_dist_level = rep(NA, sim_length), # intensity of the social distancing intervention
                               thresh_H_start = rep(NA, sim_length), # starting threshhold on total # in hospital
                               thresh_H_end = rep(NA, sim_length), # ending threshhold on total # in hospital
                               thresh_int_level = rep(NA, sim_length), # level of social distancing implemented with the threshhold intervention
                               order = "constant",
                               times = "day")

# define the pomp object
covid <- 
  scc_deaths %>%
  pomp(
    times= "day",
    t0 = 1,
    covar = intervention,
    rprocess=euler(sir_step,delta.t=1/6),
    rmeasure = rmeas, 
    dmeasure = dmeas,
    rinit=sir_init,
    partrans=parameter_trans(log=c("beta0")),
    accumvars = c("D_new"),
    paramnames=c("beta0", # base transmission rate 
                 "Ca", "Cp", "Cs", "Cm", # infectious group specific contact rates (maybe symptomatic are less likely to contact others, we don't use this yet)
                 "alpha", # fraction asymptomatic
                 "mu", # fraction of cases that are mild
                 "delta", # fraction of hospitalized cases that are fatal
                 "gamma", # rate at which exposed become infectious
                 "lambda_a", "lambda_s","lambda_m", "lambda_p", # rate of movements out of infectious categories
                 "rho", # movement out of hospitals
                 "N", # population size
                 "E0"), # number of people initially exposed
    statenames=c("S","E","Ia", 
                 "Ip","Is","Im",
                 "R","H","D", "D_new",
                 "thresh_crossed")
  ) 


# Simulation parameters
fixed_params <- c(
    Ca               = 1
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
)


# simulate 
do.call(
  pomp::simulate
  , list(
    object       = covid
    , params       = c(fixed_params, c(beta0 =2.5/7, N = 100000, E0 = 10))
    , nsim         = 10
    , format       = "d"
    , include.data = T
    , seed         = 1001)) %>% 
  ggplot() +
  geom_line(aes(x=day, 
                y = deaths,
                group=.id, 
                color = .id == "data")) + 
  xlab("Date") + ylab("New fatalities") +
  guides(color=FALSE)+
  scale_color_manual(values=c("#D5D5D3", "#24281A")) + theme_bw()

# calculate LL
pf <- pfilter(covid, 
        params = c(fixed_params, c(beta0 = 2.5/7, N = 100000, E0 = 10)),
        Np = 5000)
plot(pf)
logLik(pf)

#####
## Alternative coarse brute-force method that may be a lot faster for a first pass
## given our timeframe
#####
source("ssc_sim.R")


# local search
library(foreach)
library(doParallel)
registerDoParallel()

bake(file="local_search.rds",{
  foreach(i=1:10,.combine=c) %dopar%  
    {
      library(pomp)
      library(tidyverse)
      
      covid %>%
        mif2(
          t0 = 15,
          params=c(fixed_params, c(beta0 = 2.5/7, N = 1.938e6, E0 = 3)),
          Np=2000,
          Nmif=50,
          cooling.fraction.50=0.5,
          rw.sd=rw.sd(beta0=0.02)
        )
    }
}) -> mifs_local

plot(mifs_local)

# need to evaluate the LL with more particles
bake(file="lik_local.rds",{
  
  foreach(mf=mifs_local,.combine=rbind) %dopar% 
    {
      library(pomp)
      library(tidyverse)
      evals <- replicate(2, logLik(pfilter(mf,Np=20000)))
      ll <- logmeanexp(evals,se=TRUE)
      mf %>% coef() %>% bind_rows() %>%
        bind_cols(loglik=ll[1],loglik.se=ll[2])
    }
  
}) -> results

pairs(~loglik+beta0,data=results,pch=16)

# global search 
runifDesign(
  lower=c(beta0=1.5/7),
  upper=c(beta0=6/7),
  nseq=50
) -> guesses

bake(file="global_search.rds",{
  foreach(guess=iter(guesses,"row"), 
          .combine=rbind,
          .export=c("mf1","fixed_params")
  ) %dopar% 
    {
      library(pomp)
      library(tidyverse)
      covid %>% mif2(t0 = 15,
                     params=c(beta0 = guess, fixed_params, c(N = 1.938e6, E0 = 3)),
                     Np=2000,
                     Nmif=50,
                     cooling.fraction.50=0.5,
                     rw.sd=rw.sd(beta0=0.02)) -> mf
      ll <- replicate(10,mf %>% pfilter(Np=100000) %>% logLik())
      ll <- logmeanexp(ll,se=TRUE)
      mf %>% coef() %>% bind_rows() %>%
        bind_cols(loglik=ll[1],loglik.se=ll[2], guess = guess)
    }
}) -> results

pairs(~loglik+beta0,data=results,pch=16)
plot(guesses)
