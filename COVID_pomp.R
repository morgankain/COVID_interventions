## The epidemiological model, set up as C code. 
 ## Discrete levels of social distancing
sir_step <- Csnippet("
                     // adjust betat for social distancing interventions
                     double betat;
                     if (intervention == 3 & thresh_crossed == 1) { // 2 is for threshhold intervention
                       betat = beta0*thresh_int_level; 
                     } 
                     else if (intervention == 3 & thresh_crossed == 0) {
                       betat = beta0*back_int_level;
                     } 
                     else if (intervention == 2) betat = beta0*soc_dist_level_sip; // sip is for shelter in place
                     else if (intervention == 1) betat = beta0*soc_dist_level_wfh; // wfh is for weak wfh
                     else betat = beta0; // everything else is no intervention
                     
                     // adjust contact rates if isolation of symptomatic cases is in place
                     double iso_m = 1;
                     double iso_s = 1;
                     if(isolation == 1){
                        iso_m = iso_mild_level;
                        iso_s = iso_severe_level;
                     }
                    
                     // if import rate is above zero, draw importations, assuming they are perfectly balanced with departures of susceptible individuals
                     double import = 0;
                     if(import_rate > 0){
                      import = fmin(rpois(import_rate*dt), S);
                     }
                     // tracking of total imported, removing them them from susceptibles
                     import_total += import;
                     S -= import;
                    
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
                     rateIp[0] = mu*lambda_p; // going to mild symptomatic
                     rateIp[1] = (1-mu)*lambda_p; // going to severe symptomatic
                     reulermultinom(2, Ip, rateIp, dt, &dIp_all);
                     double dIpIm = dIp_all[0];
                     double dIpIs = dIp_all[1];
                     double dImR = rbinom(Im, 1 - exp(-lambda_m*dt));
                     
                     double rateIs[2];
                     double dIs_all[2];
                     rateIs[0] = delta*lambda_s; // hospitalized ultimately going to death
                     rateIs[1] = (1-delta)*lambda_s; // hospitalized ultimately going to recovered
                     reulermultinom(2, Is, rateIs, dt, &dIs_all);
                     double dIsHd = dIs_all[0];
                     double dIsHr = dIs_all[1];
                     double dHdD = rbinom(Hd, 1 - exp(-rho_d*dt));
                     double dHrR = rbinom(Hr, 1 - exp(-rho_r*dt));
                     
                     // update the compartments
                     S  -= dSE; // susceptible 
                     E  += dSE - dEIa - dEIp + import; // exposed
                     Ia += dEIa - dIaR; // infectious and asymptomatic
                     Ip += dEIp - dIpIs - dIpIm; // infectious and pre-symptomatic
                     Is += dIpIs - dIsHd - dIsHr; // infectious and severe symptoms (that will be hospitalized)
                     Im += dIpIm - dImR; // infectious and minor symptoms
                     I   = Ia + Ip + Im + Is; // total number of infected
                     I_new_sympt += dIpIs + dIpIm; // total number of newly symptomatic
                     Hr += dIsHr - dHrR; // hospitalized that will recover
                     Hd += dIsHd - dHdD; // hospitalizations that will die
                     H   = Hr + Hd; // total hospitalizations
                     R  += dHrR + dImR + dIaR; // recovered
                     D  += dHdD; // fatalities
                     D_new += dHdD; // daily fatalities
                     H_new += dIsHr + dIsHd; // daily new hospitalizations
                     if(intervention == 3 & H >= thresh_H_start) thresh_crossed = 1;
                     else if(intervention == 3 & thresh_crossed == 1 & H < thresh_H_end) thresh_crossed = 0;
                     ")

## Discrete levels of social distancing
sir_step_mobility <- Csnippet("
                     // adjust betat for social distancing interventions
                     double betat;
                     betat = beta0*exp(log(beta_min)*sip_prop);
                     // (mobility*(1 - beta_min) + beta_min);
                     // if import rate is above zero, draw importations, assuming they are perfectly balanced with departures of susceptible individuals
                     double import = 0;
                     if(import_rate > 0){
                      import = fmin(rpois(import_rate*dt), S);
                     }
                     // tracking of total imported, removing them them from susceptibles
                     import_total += import;
                     S -= import;
                     double iso_m = 1;
                     double iso_s = 1;
                     
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
                     rateIp[0] = mu*lambda_p; // going to mild symptomatic
                     rateIp[1] = (1-mu)*lambda_p; // going to severe symptomatic
                     reulermultinom(2, Ip, rateIp, dt, &dIp_all);
                     double dIpIm = dIp_all[0];
                     double dIpIs = dIp_all[1];
                     double dImR = rbinom(Im, 1 - exp(-lambda_m*dt));
                     
                     double rateIs[2];
                     double dIs_all[2];
                     rateIs[0] = delta*lambda_s; // hospitalized ultimately going to death
                     rateIs[1] = (1-delta)*lambda_s; // hospitalized ultimately going to recovered
                     reulermultinom(2, Is, rateIs, dt, &dIs_all);
                     double dIsHd = dIs_all[0];
                     double dIsHr = dIs_all[1];
                     double dHdD = rbinom(Hd, 1 - exp(-rho_d*dt));
                     double dHrR = rbinom(Hr, 1 - exp(-rho_r*dt));
                     
                     // update the compartments
                     S  -= dSE; // susceptible 
                     E  += dSE - dEIa - dEIp + import; // exposed
                     Ia += dEIa - dIaR; // infectious and asymptomatic
                     Ip += dEIp - dIpIs - dIpIm; // infectious and pre-symptomatic
                     Is += dIpIs - dIsHd - dIsHr; // infectious and severe symptoms (that will be hospitalized)
                     Im += dIpIm - dImR; // infectious and minor symptoms
                     I   = Ia + Ip + Im + Is; // total number of infected
                     I_new_sympt += dIpIs + dIpIm; // total number of newly symptomatic
                     Hr += dIsHr - dHrR; // hospitalized that will recover
                     Hd += dIsHd - dHdD; // hospitalizations that will die
                     H   = Hr + Hd; // total hospitalizations
                     R  += dHrR + dImR + dIaR; // recovered
                     D  += dHdD; // fatalities
                     D_new += dHdD; // daily fatalities
                     H_new += dIsHr + dIsHd; // daily new hospitalizations
                     ")

## Define the initial set up, currently, everyone is susceptible except the exposed people who seed the epidemic
 ## Different sir_init are defined here depending on the parameters chosen in the parameter code block

if (fit.E0) {
  
sir_init <- Csnippet("
                     double E0 = rpois(E_init) + 1; 
                     S = N-E0;
                     E = E0;
                     Ia = 0;
                     Ip = 0;
                     Is = 0;
                     Im = 0;
                     I = 0;
                     I_new_sympt = 0;
                     Hr = 0;
                     Hd = 0;
                     H = Hd + Hr;
                     R = 0;
                     D = 0;
                     D_new = 0;
                     H_new = 0;
                     thresh_crossed = 0;
                     import_total = 0;
                     ")  
  
} else {
  
sir_init <- Csnippet("
                     S = N-E0;
                     E = E0;
                     Ia = 0;
                     Ip = 0;
                     Is = 0;
                     Im = 0;
                     I = 0;
                     I_new_sympt = 0;
                     Hr = 0;
                     Hd = 0;
                     H = Hd + Hr;
                     R = 0;
                     D = 0;
                     D_new = 0;
                     H_new = 0;
                     thresh_crossed = 0;
                     import_total = 0;
                     ")  

}

# define random simulator of measurement. Which is used will depend on parameter choices
# define evaluation of model prob density function. Which is used will depend on parameter choices
if (meas.nb) {
rmeas_deaths <- Csnippet("double tol = 1e-16;
                   deaths = rnbinom_mu(theta, D_new + tol);
                  ")
rmeas_hosp <- Csnippet("double tol = 1e-16;
                   hosp = rnbinom_mu(theta, H + tol);
                  ")
  
dmeas_deaths <- Csnippet("double tol = 1e-16;
                   lik = dnbinom_mu(deaths, theta, D_new + tol, give_log);
                  ")  
dmeas_hosp <- Csnippet("double tol = 1e-16;
                   lik = dnbinom_mu(hosp, theta, H + tol, give_log);
                  ") 
} else {
rmeas_deaths <- Csnippet("double tol = 1e-16;
                   deaths = rpois(D_new + tol);
                  ")
rmeas_hosp <- Csnippet("double tol = 1e-16;
                   hosp = rpois(H + tol);
                  ")

dmeas_deaths <- Csnippet("double tol = 1e-16;
                   lik = dpois(deaths, D_new + tol, give_log);
                  ")
dmeas_hosp <- Csnippet("double tol = 1e-16;
                   lik = dpois(hosp, H + tol, give_log);
                  ") 
}

## Parameters to fit and use. Again, which are used will depend on parameter choices 
if (fitting) {
  if (fit_to_sip) {
    
if (fit.E0) {
par_trans   <- parameter_trans(log = c("beta0", "import_rate", "E_init"),
                             logit = c("soc_dist_level_sip")) 
param_names <- c(
   "beta0"
  , "Ca", "Cp", "Cs", "Cm"
  , "alpha"
  , "mu"
  , "delta"
  , "gamma"
  , "lambda_a", "lambda_s", "lambda_m", "lambda_p"
  , "rho_r"
  , "rho_d"
  , "N"
  , "E_init"
  , "soc_dist_level_sip"
  , "import_rate"
)

} else {
par_trans   <- parameter_trans(log = c("beta0", "import_rate"),
                            logit = c("soc_dist_level_sip"))  
param_names <- c(
   "beta0"
  , "Ca", "Cp", "Cs", "Cm"
  , "alpha"
  , "mu"
  , "delta"
  , "gamma"
  , "lambda_a", "lambda_s", "lambda_m", "lambda_p"
  , "rho_r"
  , "rho_d"
  , "N"
  , "E0"
  , "soc_dist_level_sip"
  , "import_rate"
)

}

  } else {

if (fit.E0) {
par_trans   <- parameter_trans(log = c("beta0", "import_rate", "E_init"))
param_names <- c(
   "beta0"
  , "Ca", "Cp", "Cs", "Cm"
  , "alpha"
  , "mu"
  , "delta"
  , "gamma"
  , "lambda_a", "lambda_s", "lambda_m", "lambda_p"
  , "rho_r"
  , "rho_d"
  , "N"
  , "E_init"
  , "import_rate"
)

} else {
par_trans   <- parameter_trans(log = c("beta0", "import_rate"))  
param_names <- c(
   "beta0"
  , "Ca", "Cp", "Cs", "Cm"
  , "alpha"
  , "mu"
  , "delta"
  , "gamma"
  , "lambda_a", "lambda_s", "lambda_m", "lambda_p"
  , "rho_r"
  , "rho_d"
  , "N"
  , "E0"
  , "import_rate"
)

}

}

} else {
  
if (fit.E0) {
par_trans   <- parameter_trans(log = c("beta0", "E_init"))
param_names <- c(
   "beta0"
  , "Ca", "Cp", "Cs", "Cm"
  , "alpha"
  , "mu"
  , "delta"
  , "gamma"
  , "lambda_a", "lambda_s", "lambda_m", "lambda_p"
  , "rho_r"
  , "rho_d"
  , "N"
  , "E_init"
  , "import_rate"
)

} else {
par_trans   <- parameter_trans(log = c("beta0")) 
param_names <- c(
   "beta0"
  , "Ca", "Cp", "Cs", "Cm"
  , "alpha"
  , "mu"
  , "delta"
  , "gamma"
  , "lambda_a", "lambda_s", "lambda_m", "lambda_p"
  , "rho_r"
  , "rho_d"
  , "N"
  , "E0"
  , "import_rate"
)

}

}

## Variables that are zeroed after each time step (new per time step)
accum_names <- c("D_new", "H_new", "I_new_sympt")

## State variables
state_names <- c(
    "S" , "E" , "Ia"
  , "Ip", "Is", "Im"
  , "I" , "I_new_sympt"
  , "H" , "Hr", "Hd"
  , "R" , "D" 
  , "D_new", "H_new" 
  , "thresh_crossed"
  , "import_total"
)

## IN FLUX attempt for likelihood for multiple data streams:
if (fit.with == "D_C") {
rmeas_multi <- Csnippet("double tol = 1e-16;
                   deaths = rnbinom_mu(theta, D_new + tol);
                   cases  = rnbinom_mu(theta2, detect*I_new_sympt + tol);
                  ")
# define evaluation of model prob density function
dmeas_multi <- Csnippet("double tol = 1e-16;
                   lik = dnbinom_mu(deaths, theta, D_new + tol, give_log) + dnbinom_mu(cases, theta2, detect*I_new_sympt + tol, give_log);
                  ")
               ##  lik = (give_log) ? lik : exp(lik);

par_trans <- parameter_trans(log = c("beta0", "import_rate", "E_init", "theta", "theta2"),
                            logit = c("beta_min"))

sir_init <- Csnippet("
                     double E0 = rpois(E_init) + 1; 
                     S = N-E0;
                     E = E0;
                     Ia = 0;
                     Ip = 0;
                     Is = 0;
                     Im = 0;
                     I = 0;
                     I_new_sympt = 0;
                     Hr = 0;
                     Hd = 0;
                     H = Hd + Hr;
                     R = 0;
                     D = 0;
                     D_new = 0;
                     H_new = 0;
                     import_total = 0;
                     ")

param_names <- c(
   "beta0"
  , "Ca", "Cp", "Cs", "Cm"
  , "alpha"
  , "mu"
  , "delta"
  , "gamma"
  , "lambda_a", "lambda_s", "lambda_m", "lambda_p"
  , "rho_r"
  , "rho_d"
  , "N"
  , "E_init"
  , "import_rate"
  , "theta"
  , "theta2"
  , "beta_min"
  , "detect"
)

state_names = c(
    "S" , "E" , "Ia"
  , "Ip", "Is", "Im"
  , "I" , "I_new_sympt"
  , "H" , "Hr", "Hd"
  , "R" , "D" 
  , "D_new", "H_new" 
  , "import_total"
)

}

## R0 here just based on the simple transmission rate / recovery rate (weighted by the probability of going into different classes)
covid_R0 <- function (beta0est, fixed_params, sd_strength, prop_S) {
## transmission rate
 R <- beta0est * prop_S * sd_strength * 
    (                
## proportion * time in asymptomatic
      fixed_params["alpha"] * fixed_params["Ca"] * (1/fixed_params["lambda_a"]) +                  
## proportion * time in mildly symptomatic
      (1 - fixed_params["alpha"]) * fixed_params["mu"] * ((1/fixed_params["lambda_p"]) + (1/fixed_params["lambda_m"])) +    
## proportion * time in severely symptomatic
      (1 - fixed_params["alpha"]) * (1 - fixed_params["mu"]) * ((1/fixed_params["lambda_p"]) + (1/fixed_params["lambda_s"]))      
  )
 
 unlist(R)
}

