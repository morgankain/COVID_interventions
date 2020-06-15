## Discrete levels of social distancing
sir_step_mobility <- Csnippet("
                      // adjust betat for social distancing interventions
                      // double betat;
                     
                      int n_beta = round(fmax(I, 1));
                      
                      if (beta0_sigma > 0) {
                      if ((catch_eff > 0) & intervention == 1) {
                        betat = rtgamma_eff(n_beta, sqrt(beta0/ beta0_sigma), beta0, beta_catch, catch_eff)*exp(log(beta_min)*sip_prop);
                      } else {
                        betat = rgammawn(sqrt(beta0 / (beta0_sigma * n_beta)), beta0) * exp(log(beta_min)*sip_prop);
                      }
                      } else {
                        betat = beta0*exp(log(beta_min)*sip_prop);
                      }
                     
                     // if import rate is above zero, draw importations, assuming they are perfectly balanced with departures of susceptible individuals
                     double import = 0;
                     if (import_rate > 0){
                      import = fmin(rpois(import_rate*dt), S);
                     }
                     
                     // tracking of total imported, removing them them from susceptibles
                     import_total += import;
                     S -= import;

                     // adjust contact rates if isolation of symptomatic cases is in place
                     double iso_m = 1;
                     double iso_s = 1;
                     if (intervention == 2){
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

trunc_gamma <- Csnippet("
static double rtgamma_eff(int N, double sigma, double mu, 
                                 double trim, double eff){
  double out = 0; 
  for(int i = 0; i < N; ++i) {
    double ind = rgammawn(sigma, mu);
    if(runif(0, 1) <= eff){ // if we're effective
      while(ind > trim){ // and its above the cutoff, resample
        ind = rgammawn(sigma, mu);
      }
    }
    out += ind;
  }
  return out/N;
}")

sir_init_mid <- Csnippet("
                     S = S0;
                     E = ceil(E_init);
                     Ia = Ia0;
                     Ip = Ip0;
                     Is = Is0;
                     Im = Im0;
                     I = Ia + Ip + Is + Im;
                     I_new_sympt = 0;
                     Hr = Hr0;
                     Hd = Hd0;
                     H = Hd + Hr;
                     R = R0;
                     D = D0;
                     D_new = 0;
                     H_new = 0;
                     import_total = 0;
                     ")
mid_init_param_names = c("S0", "Ia0", "Ip0", "Is0", "Im0",
                         "Hr0", "Hd0", "R0", "D0")

if (fixed.E0) {
  
sir_init <- Csnippet("
                     S = N-5;
                     E = 5;
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
  
} else {

sir_init <- Csnippet("
                     double E0 = ceil(E_init);
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
    
}

rmeas_deaths <- Csnippet("double tol = 1e-16;
                   deaths = rnbinom_mu(theta, D_new + tol);
                  ")

dmeas_deaths <- Csnippet("double tol = 1e-16;
                   lik = dnbinom_mu(deaths, theta, D_new + tol, give_log);
                  ")  

if (con_theta) {

rmeas_multi_logis <- Csnippet("double tol = 1e-16;
                   double detect;
                   if (t < detect_t0) {
                   detect = 0;
                   } else {
                   detect = detect_max / (1 + exp(-detect_k * (t - detect_mid)));                   
                   }
                   deaths = rnbinom_mu(theta, D_new + tol);
                   cases  = rnbinom_mu(theta*theta2, detect*I_new_sympt + tol);
                  ")

# define evaluation of model prob density function
dmeas_multi_logis <- Csnippet("double tol = 1e-16;
                   double detect;
                   if (t < detect_t0) {
                   detect = 0;
                   } else {
                   detect = detect_max / (1 + exp(-detect_k * (t - detect_mid)));                   
                   }
            
                   if (ISNA(deaths) & ISNA(cases)) {
                   lik = 0 + 0;
                   } else if (ISNA(deaths) & !ISNA(cases)) {
                   lik = 0 + dnbinom_mu(cases, theta*theta2, detect*I_new_sympt + tol, 1);
                   } else if (!ISNA(deaths) & ISNA(cases)) {
                   lik = dnbinom_mu(deaths, theta, D_new + tol, 1) + 0;
                   } else {
                   lik = dnbinom_mu(deaths, theta, D_new + tol, 1) + dnbinom_mu(cases, theta*theta2, detect*I_new_sympt + tol, 1);
                   }
                   lik = (give_log) ? lik : exp(lik);
                  ")

  par_trans <- parameter_trans(log = c("beta0", "beta0_sigma", "import_rate"
                                       , "E_init"
                                       , "theta"
                                       , "detect_k"
                                       , "detect_mid"
  )
  , logit = c("beta_min", "detect_max", "theta2")) 

} else {
  
  rmeas_multi_logis <- Csnippet("double tol = 1e-16;
                   double detect;
                   if (t < detect_t0) {
                   detect = 0;
                   } else {
                   detect = detect_max / (1 + exp(-detect_k * (t - detect_mid)));                   
                   }
                   deaths = rnbinom_mu(theta, D_new + tol);
                   cases  = rnbinom_mu(theta2, detect*I_new_sympt + tol);
                  ")
  
  # define evaluation of model prob density function
  dmeas_multi_logis <- Csnippet("double tol = 1e-16;
                   double detect;
                   if (t < detect_t0) {
                   detect = 0;
                   } else {
                   detect = detect_max / (1 + exp(-detect_k * (t - detect_mid)));                   
                   }
            
                   if (ISNA(deaths) & ISNA(cases)) {
                   lik = 0 + 0;
                   } else if (ISNA(deaths) & !ISNA(cases)) {
                   lik = 0 + dnbinom_mu(cases, theta2, detect*I_new_sympt + tol, 1);
                   } else if (!ISNA(deaths) & ISNA(cases)) {
                   lik = dnbinom_mu(deaths, theta, D_new + tol, 1) + 0;
                   } else {
                   lik = dnbinom_mu(deaths, theta, D_new + tol, 1) + dnbinom_mu(cases, theta2, detect*I_new_sympt + tol, 1);
                   }
                   lik = (give_log) ? lik : exp(lik);
                  ")
  
par_trans <- parameter_trans(log = c("beta0", "beta0_sigma", "import_rate"
  , "E_init"
  , "theta"
  , "theta2"
  , "detect_k"
  , "detect_mid"
  )
  , logit = c("beta_min", "detect_max")) 
  
}

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
  , "beta0_sigma"
  , "detect_k"
  , "detect_mid"
  , "detect_max"
  , "beta_catch"
  , "catch_eff"
)
  

state_names = c(
    "S" , "E" , "Ia"
  , "Ip", "Is", "Im"
  , "I" , "I_new_sympt"
  , "H" , "Hr", "Hd"
  , "R" , "D" 
  , "D_new", "H_new" 
  , "import_total"
  , "betat"
)

accum_names <- c("D_new", "H_new", "I_new_sympt")

## R0 here just based on the simple transmission rate / recovery rate (weighted by the probability of going into different classes)
covid_R0   <- function (beta0est, fixed_params, sd_strength, prop_S) {
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
logis_func <- function (L, k, mid, day) {
  L / (1 + exp(-k * (day - mid)))
}
beta_func  <- function (beta0, beta_min, sip_prop) {
  beta0 * exp(log(beta_min) * sip_prop)
}
beta_func_lin   <- function (beta0, beta_min, sip_prop) {
  beta0 * (1 - sip_prop) + beta_min
}
beta_func_lin2  <- function (beta0, beta_min, sip_prop) {
beta0 * (1 - (1 - beta_min)*sip_prop)
}
beta_func_lin3  <- function (beta0, beta_min, sip_prop) {
beta0 * (1 - sip_prop*beta_min)
}



