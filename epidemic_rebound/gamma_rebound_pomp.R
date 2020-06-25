sir_step_mobility <- Csnippet("
                      // adjust betat for social distancing interventions
                      // double betat;
                     
                      int n_beta = round(fmax(I, 1));
                      
                      if (thresh_crossed == 0) {
                      
                      if (beta0_sigma > 0) {
                      if ((catch_eff > 0) & intervention == 1) {
                        betat = rtgamma_eff(n_beta, sqrt(beta0/ beta0_sigma), beta0, beta_catch, catch_eff)*exp(log(beta_min)*sip_prop);
                      } else {
                        betat = rgammawn(sqrt(beta0 / (beta0_sigma * n_beta)), beta0) * exp(log(beta_min)*sip_prop);
                      }
                      } else {
                        betat = beta0*exp(log(beta_min)*sip_prop);
                      }
                      
                      } else {
                      
                      if (beta0_sigma_post > 0) {
                      if ((catch_eff_post > 0) & intervention == 1) {
                        betat = rtgamma_eff(n_beta, sqrt(beta0/ beta0_sigma_post), beta0, beta_catch_post, catch_eff_post)*exp(log(beta_min)*(sip_prop*sip_prop_scaling));
                      } else {
                        betat = rgammawn(sqrt(beta0 / (beta0_sigma_post * n_beta)), beta0) * exp(log(beta_min)*(sip_prop*sip_prop_scaling));
                      }
                      } else {
                        betat = beta0*exp(log(beta_min)*(sip_prop*sip_prop_scaling));
                      }
                        
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
                     if (intervention == 2 & thresh_crossed == 0)  {
                      iso_m = iso_mild_level;
                      iso_s = iso_severe_level;
                     } else if (intervention == 2 & thresh_crossed == 1) {
                      iso_m = iso_mild_level_post;
                      iso_s = iso_severe_level_post; 
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
                     
                     if (I < thresh_inf & check_int == 0) thresh_crossed = 1;
            //       else if (thresh_crossed == 1 & I > thresh_inf) thresh_crossed = 0;
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
                     thresh_crossed = 0;
                     ")

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
  , "detect_k"
  , "detect_mid"
  , "theta2"
  )
  , logit = c("beta_min", "detect_max")) 

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
  , "catch_eff_post"
  , "beta_catch_post"
  , "beta0_sigma_post"
  , "sip_prop_scaling"
)

state_names <- c(
    "S" , "E" , "Ia"
  , "Ip", "Is", "Im"
  , "I" , "I_new_sympt"
  , "H" , "Hr", "Hd"
  , "R" , "D" 
  , "D_new", "H_new" 
  , "import_total"
  , "thresh_crossed"
  , "betat"
)

accum_names <- c("D_new", "H_new", "I_new_sympt")

covid_R0   <- function (beta0est, fixed_params, sd_strength, prop_S, inf_iso, iso_course_mild, iso_course_severe) {
  
if (!inf_iso) {
  
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
 
} else {
 
## transmission rate
 R <- beta0est * prop_S * sd_strength * 
    (                
## proportion * time in asymptomatic
      fixed_params["alpha"] * fixed_params["Ca"] * (1/fixed_params["lambda_a"]) +                  
## proportion * time in mildly symptomatic
      (1 - fixed_params["alpha"]) * fixed_params["mu"] * ((1/fixed_params["lambda_p"]) + (1/fixed_params["lambda_m"])) * iso_course_mild +    
## proportion * time in severely symptomatic
      (1 - fixed_params["alpha"]) * (1 - fixed_params["mu"]) * ((1/fixed_params["lambda_p"]) + (1/fixed_params["lambda_s"])) * iso_course_severe      
  ) 
  
}
 
 unlist(R)
}

check_R0 <- function (beta0est, beta_min, fixed_params, sd_strength, prop_S, desired_R) {
  
  d <- {(                
    ## proportion * time in asymptomatic
    fixed_params["alpha"] * fixed_params["Ca"] * (1/fixed_params["lambda_a"]) +                  
      ## proportion * time in mildly symptomatic
      (1 - fixed_params["alpha"]) * fixed_params["mu"] * ((1/fixed_params["lambda_p"]) + (1/fixed_params["lambda_m"])) +    
      ## proportion * time in severely symptomatic
      (1 - fixed_params["alpha"]) * (1 - fixed_params["mu"]) * ((1/fixed_params["lambda_p"]) + (1/fixed_params["lambda_s"]))      
  )}
  
  ## SIP to get desired R without any tail cutting
  log(desired_R / (d * beta0est)) / log(beta_min)
  
}
## SIP_prop for Santa Clara at the start -- 0.2386557

