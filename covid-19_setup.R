# model adapted from Aaron King's boarding school flu example https://kingaa.github.io/pomp/FAQ.html
# https://kingaa.github.io/short-course/stochsim/stochsim.html#the-boarding-school-flu-outbreak
# and from John Drake and Pejman Rohani's model https://github.com/CEIDatUGA/ncov-wuhan-stochastic-model/blob/master/stochastic-model.pdf

needed_packages <- c(
  "pomp"
, "plyr"
, "dplyr"
, "ggplot2"
, "magrittr"
, "scales"
, "lubridate"
)

lapply(needed_packages, require, character.only = TRUE)

# To do: 
# control caused by crossing a threshhold --> seems dependent on time step because if people are reevaluating by fractions of a day, they can respond faster vs if the time step is 1 day?
# add case detection? 
# case importations?
# loss of immunity 
# update parameters to acutally reflect the parameter estimates 
# add some information on the age structure of our population? to estimate fatality rates?
# add hospital capacity affecting mortality?
# add contact tracing as an intervention which will require another category for exposed people who know they are exposed


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
                     sympt_new  +=  dIpIs + dIpIm;
                     H_new += dIsH;
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
                     sympt_new = 0;
                     H_new = 0;
                     thresh_crossed = 0;
                     ")

# define simulation length, set up data frame, right now it just has days but 0s for the observation
sim_start = as.Date("2020-01-15")
sim_end = as.Date("2020-12-01")
sim_length =  sim_end - sim_start
dat = data.frame(day= 1:sim_length, 
                 B = rep(0, sim_length))
 
# in Wuhan, the intervention started around January 23
# # use the intervention info to construct a covariate table for use in the pomp object
intervention = covariate_table(day= 1:(sim_length),
                               intervention = c(rep(0, as.Date("2020-03-17") - sim_start - 1), # start with no intervention until SCC intervened
                                                rep(1, 21),  # go to social distancing 
                                                rep(0, sim_length - (as.Date("2020-03-17")- sim_start - 1) - 21)), # threshhold based interventiosn after
                               isolation = c(rep(0, as.Date("2020-03-17") - sim_start - 1), # 0 is no isolation of cases,
                                             rep(1, sim_length - (as.Date("2020-03-17") - sim_start - 1))), # 1 is isolation of symptomatic
                               iso_severe_level = rep(0, sim_length), # % of contats that severe cases maintain
                               iso_mild_level = rep(0.1, sim_length), # % of contats that mild cases maintain
                               soc_dist_level = rep(0.2, sim_length), # intensity of the social distancing intervention
                               thresh_H_start = rep(20, sim_length), # starting threshhold on total # in hospital
                               thresh_H_end = rep(2, sim_length), # ending threshhold on total # in hospital
                               thresh_int_level = rep(0.2, sim_length), # level of social distancing implemented with the threshhold intervention
                               order = "constant",
                               times = "day")

# define the pomp object
covid <- dat %>%
  pomp(
    time= "day",
    t0=1,
    covar = intervention,
    rprocess=euler(sir_step,delta.t=1/6),
    rinit=sir_init,
    accumvars= c("sympt_new", "H_new"), # accumulate H until it gets measured, then zero it
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
                 "R", "H","D", 
                 "sympt_new", "H_new",
                 "thresh_crossed")
  ) 
