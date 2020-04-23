# for fitting E0 added E_init parameter that would need to be discarded
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
                     H_new += dIsH; // daily new hospitalizations
                     if(intervention == 2 & H >= thresh_H_start) thresh_crossed = 1;
                     else if(intervention == 2 & thresh_crossed == 1 & H < thresh_H_end) thresh_crossed = 0;
                     ")

# define the initial set up, currently, every is susceptible except the exposed people
sir_init <- Csnippet("
                     double E0 = rpois(E_init);    
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
                     H_new = 0;
                     thresh_crossed = 0;
                     ")

# define random simulator of measurement
rmeas_deaths <- Csnippet("double tol = 1e-16;
                   deaths = rpois(D_new + tol);
                  ")
# define evaluation of model prob density function
dmeas_deaths <- Csnippet("double tol = 1e-16;
                   lik = dpois(deaths, D_new + tol, give_log);
                  ")
# define random simulator of measurement
rmeas_hosp <- Csnippet("double tol = 1e-16;
                   hosp = rpois(H + tol);
                  ")
# define evaluation of model prob density function
dmeas_hosp <- Csnippet("double tol = 1e-16;
                   lik = dpois(hosp, H + tol, give_log);
                  ")

# parameters to transform
par_trans = parameter_trans(log = c("beta0", "E_init"),
                            logit = c("soc_dist_eff"))

# variables that should be zeroed after each obs
accum_names = c("D_new", "H_new")

# state variables
state_names = c(
    "S" , "E" , "Ia"
  , "Ip", "Is", "Im"
  , "R" , "H" ,"D" 
  , "D_new", "H_new" 
  , "thresh_crossed"
)
# parameter names
param_names = c(
  "beta0"
  , "Ca", "Cp", "Cs", "Cm"
  , "alpha"
  , "mu"
  , "delta"
  , "gamma"
  , "lambda_a", "lambda_s", "lambda_m", "lambda_p"
  , "rho"
  , "N"
  , "E_init"
  , "soc_dist_level"
)

## R0 here just based on the simple transmission rate / recovery rate (weighted by the probability of going into different classes)
covid_R0 <- function (beta0est, fixed_params, sd_strength) {
## transmission rate
    beta0est * sd_strength * 
    (                
## proportion * time in asymptomatic
      fixed_params["alpha"] * fixed_params["Ca"] * (1/fixed_params["lambda_a"]) +                  
## proportion * time in mildly symptomatic
      (1 - fixed_params["alpha"]) * fixed_params["mu"] * ((1/fixed_params["lambda_p"]) + (1/fixed_params["lambda_m"])) +    
## proportion * time in severely symptomatic
      (1 - fixed_params["alpha"]) * (1 - fixed_params["mu"]) * ((1/fixed_params["lambda_p"]) + (1/fixed_params["lambda_s"]))      
      )
}

