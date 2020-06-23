# incremental step ----
sir_step_mobility <- Csnippet("
                     // adjust betat for social distancing interventions
                //   double betat;
                     double *E = &E1;  // Set up E boxes
                     double *Ia = &Ia1;  // Set up Ia boxes
                     double *Ip = &Ip1;  // Set up Ip boxes
                     double *Im = &Im1;  // Set up Im boxes
                     double *Is = &Is1;  // Set up Is boxes

                     double import = 0;
                     if (import_rate > 0){
                      import = fmin(rpois(import_rate*dt), S);
                     }
                     
                     // removing imported from susceptibles first
                     S -= import;
                     
                     // adjust contact rates if isolation of symptomatic cases is in place
                     double iso_m = 1;
                     double iso_s = 1;
                     if (intervention == 2) {
                      iso_m = iso_mild_level;
                      iso_s = iso_severe_level;
                     }

                     // Sum up individuals in Ia boxes
                     double sIa = 0;
                     for(int i = 0; i < nIa; i++) {
                       sIa += Ia[i];      
                     }
                     
                     // Sum up individuals in Ip boxes
                     double sIp = 0;
                     for(int i = 0; i < nIp; i++) {
                       sIp += Ip[i];      
                     }
                     
                     // Sum up individuals in Im boxes
                     double sIm = 0;
                     for(int i = 0; i < nIm; i++) {
                       sIm += Im[i];      
                     }
                     
                     // Sum up individuals in Is boxes
                     double sIs = 0;
                     for(int i = 0; i < nIs; i++) {
                       sIs += Is[i];      
                     }
                     
                    double I = sIa + sIp + sIm + sIs;
                    int n_beta = round(fmax(I, 1));
                    
                     if ((beta0_k > 0) & (I > 0)) {
                      if ((catch_eff > 0) & intervention == 1) {
                        betat = rtgamma_eff(n_beta, sqrt(beta0 / beta0_k), beta0*dt/d, beta_catch, catch_eff)*exp(log(beta_min)*sip_prop);
                        // betat = rgammawn(sqrt(beta0 / (beta0_k * I)), beta0*dt/d) * exp(log(beta_min)*sip_prop); 
                      } else {
                        betat = rgammawn(sqrt(beta0 / (beta0_k * I)), beta0*dt/d) * exp(log(beta_min)*sip_prop);                
                      }
                     } else {
                       betat = beta0*dt/d*exp(log(beta_min)*sip_prop);
                     } 
                     
                     // transition out of S to E
                     double dSE = rbinom(S, 1-exp(-betat*(Ca*sIa/N + Cp*sIp/N + iso_m*Cm*sIm/N + iso_s*Cs*sIs/N))); 
                     
                     // transitions between E boxes
                     double dE[nE];       
                     for(int i = 0; i < nE - 1; i++){
                         dE[i] = rbinom(E[i], 1.0 - exp(-gamma*dt*nE)); // number leaving E[i]
                     }
                     
                     // transitions from last E to Ia and Ip types 
                     // Rates for these transitions
                     double rateE[2];
                     rateE[0] = alpha*gamma*nE;     // going to asymtomatic
                     rateE[1] = (1-alpha)*gamma*nE; // going to presymptomatic
                     // random draws given the competing rates
                     double dE_all[2];
                     reulermultinom(2, E[nE - 1], rateE, dt, &dE_all);
                     double dEIa = dE_all[0]; // going to the first box of asymptomatic infecteds
                     double dEIp = dE_all[1]; // going to first box of presymptomatic infecteds
                     
                     // Movement among Ia boxes
                     double dIa[nIa];       
                     for(int i = 0; i < nIa; i++){
                         dIa[i] = rbinom(Ia[i], 1.0 - exp(-lambda_a*dt*nIa)); // number leaving Ia[i]
                     }
                     
                     // Movement between Ip boxes
                     double dIp[nIp];
                     for(int i = 0; i < nIp - 1; i++){
                         dIp[i] = rbinom(Ip[i], 1.0 - exp(-lambda_p*dt*nIp)); // number leaving Ip[i]
                     }
                     
                     // transitions from last Ip box to Im and Is
                     // Rates for these transitions
                     double rateIp[2];
                     double dIp_all[2];
                     rateIp[0] = mu*lambda_p*nIp; // going to mild symptomatic
                     rateIp[1] = (1-mu)*lambda_p*nIp; // going to severe symptomatic
                     // random draws given the competing rates
                     reulermultinom(2, Ip[nIp - 1], rateIp, dt, &dIp_all);
                     double dIpIm = dIp_all[0];  // going to mild symptomatic
                     double dIpIs = dIp_all[1];  // going to severe symptomatic
                     
                     // transitions among Im compartments
                     double dIm[nIm];
                     for(int i = 0; i < nIm; i++){
                         dIm[i] = rbinom(Im[i], 1.0 - exp(-lambda_m*dt*nIm)); // number leaving Im[i]
                     }
                     
                     // transitions among Is compartments
                     double dIs[nIs];
                     for(int i = 0; i < nIs - 1; i++){
                         dIs[i] = rbinom(Is[i], 1.0 - exp(-lambda_s*dt*nIs)); // number leaving Is[i]
                     }
                     
                     // transitions from last Is box to Hr and Hd
                     // Rates for these transitions
                     double rateIs[2];
                     double dIs_all[2];
                     rateIs[0] = delta*lambda_s*nIs; // hospitalized ultimately going to death
                     rateIs[1] = (1-delta)*lambda_s*nIs; // hospitalized ultimately going to recovered
                     // random draws given the competing rates
                     reulermultinom(2, Is[nIs - 1], rateIs, dt, &dIs_all);
                     double dIsHd = dIs_all[0]; // hospitalized ultimately going to death
                     double dIsHr = dIs_all[1]; // hospitalized ultimately going to recovered
                     
                     // transitions from hospitalized to recovered and dead
                     double dHdD = rbinom(Hd, 1.0 - exp(-rho_d*dt));
                     double dHrR = rbinom(Hr, 1.0 - exp(-rho_r*dt));
                     
                     // update the compartments
                     //Rprintf(\"NEW time = %g, S = %g, dSE = %g, betat = %g\\n\", t, S, dSE, betat);
                     S  -= dSE; // susceptible 
                     
                     E[0] += dSE - dE[0] + import; // exposed
                     for (int i = 1; i < nE - 1; i++) {
                          E[i] += dE[i-1] - dE[i];
                     }
                     E[nE - 1] += dE[nE - 2] - dEIa - dEIp;
                    
                     Ia[0] += dEIa - dIa[0]; // asymptomatic
                     for (int i=1; i < nIa; i++) {
                          Ia[i] += dIa[i-1] - dIa[i];
                     }
                     
                     Ip[0] += dEIp - dIp[0]; // presymptomatic
                     for (int i = 1; i < nIp - 1; i++) {
                          Ip[i] += dIp[i-1] - dIp[i];
                     }
                     Ip[nIp -1] += dIp[nIp - 2] - dIpIm - dIpIs;
                     
                     Im[0] += dIpIm - dIm[0]; // mild symptomatic
                     for (int i=1; i < nIm; i++) {
                          Im[i] += dIm[i-1] - dIm[i];
                     }

                     Is[0] += dIpIs - dIs[0]; // severe symptomatic
                     for (int i=1; i < nIs - 1; i++) {
                          Is[i] += dIs[i-1] - dIs[i];
                     }
                     Is[nIs -1] += dIs[nIs - 2] - dIsHr - dIsHd;
                     
                     I_new_sympt += dIpIs + dIpIm; // total number of newly symptomatic
                     Hr += dIsHr - dHrR; // hospitalized that will recover
                     Hd += dIsHd - dHdD; // hospitalizations that will die
                     R  += dHrR + dIm[nIm - 1] + dIa[nIa - 1]; // recovered
                     D  += dHdD; // fatalities
                     D_new += dHdD; // daily fatalities
                     ")

sir_init <- Csnippet("
                     double E0 = ceil(E_init);
                     S = N-E0;
                     
                     double *E = &E1;
                     for (int j = 0; j < nE; j++) E[j] = 0;
                     E[0] += E0;
                     
                     double *Ia = &Ia1;
                     for (int j = 0; j < nIa; j++) Ia[j] = 0;
                     
                     double *Ip = &Ip1;
                     for (int j = 0; j < nIp; j++) Ip[j] = 0;
                     
                     double *Im = &Im1;
                     for (int j = 0; j < nIm; j++) Im[j] = 0;
                     
                     double *Is = &Is1;
                     for (int j = 0; j < nIs; j++) Is[j] = 0;
                     
                     I_new_sympt = 0;
                     Hr = 0;
                     Hd = 0;
                     R = 0;
                     D = 0;
                     D_new = 0;
                     ")
    
# truncated gamma function ----
globs <- Csnippet(paste0("
int nE = ",  nE,  ";  // Number of E boxes
int nIa = ", nIa, ";  // Number of Ia boxes
int nIp = ", nIp, ";  // Number of Ip boxes
int nIm = ", nIm, ";  // Number of Im boxes
int nIs = ", nIs, ";  // Number of Is boxes

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
}

")
  )



# measures ----
rmeas_multi_logis_con <- Csnippet("double tol = 1e-16;
                   double detect;
                   if (t < detect_t0) {
                   detect = 0;
                   } else {
                   detect = detect_max / (1 + exp(-detect_k * (t - detect_mid)));                   
                   }
                   deaths = rnbinom_mu(theta, D_new + tol);
                   cases  = rnbinom_mu(theta*theta2, detect*I_new_sympt + tol);
                  ")

rmeas_multi_logis_ind <- Csnippet("double tol = 1e-16;
                   double detect;
                   if (t < detect_t0) {
                   detect = 0;
                   } else {
                   detect = detect_max / (1 + exp(-detect_k * (t - detect_mid)));                   
                   }
                   deaths = rnbinom_mu(theta, D_new + tol);
                   cases  = rnbinom_mu(theta2, detect*I_new_sympt + tol);
                  ")

# dmeasures ----
dmeas_multi_logis_con <- Csnippet("double tol = 1e-16;
                   double detect;
                   if (t < detect_t0) {
                   detect = 0;
                   } else {
                   detect = detect_max / (1 + exp(-detect_k * (t - detect_mid)));                   
                   }
            
                   if (ISNA(deaths) & ISNA(cases)) {
                   lik = 0 + 0;
                   } else if (ISNA(deaths) & (!ISNA(cases))) {
                   lik = 0 + dnbinom_mu(cases, theta*theta2, detect*I_new_sympt + tol, 1);
                   } else if ((!ISNA(deaths)) & ISNA(cases)) {
                   lik = dnbinom_mu(deaths, theta, D_new + tol, 1) + 0;
                   } else {
                   lik = dnbinom_mu(deaths, theta, D_new + tol, 1) + dnbinom_mu(cases, theta*theta2, detect*I_new_sympt + tol, 1);
                   }
                   lik = (give_log) ? lik : exp(lik);
                  ")

dmeas_multi_logis_ind <- Csnippet("double tol = 1e-16;
                   double detect;
                   if (t < detect_t0) {
                   detect = 0;
                   } else {
                   detect = detect_max / (1 + exp(-detect_k * (t - detect_mid)));                   
                   }
            
                   if (ISNA(deaths) & ISNA(cases)) {
                   lik = 0 + 0;
                   } else if (ISNA(deaths) & (!ISNA(cases))) {
                   lik = 0 + dnbinom_mu(cases, theta2, detect*I_new_sympt + tol, 1);
                   } else if ((!ISNA(deaths)) & ISNA(cases)) {
                   lik = dnbinom_mu(deaths, theta, D_new + tol, 1) + 0;
                   } else {
                   lik = dnbinom_mu(deaths, theta, D_new + tol, 1) + dnbinom_mu(cases, theta2, detect*I_new_sympt + tol, 1);
                   }
                   lik = (give_log) ? lik : exp(lik);
                  ")

# parameter transformation and names ----
par_trans_con <- parameter_trans(
  log = c("beta0", "beta0_k", "import_rate"
  , "E_init", "theta", "detect_k", "detect_mid")
  , logit = c("beta_min", "detect_max", "theta2")) 

par_trans_ind <- parameter_trans(
  log = c("beta0", "beta0_k", "import_rate"
          , "E_init", "theta", "theta2", "detect_k", "detect_mid")
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
  , "beta0_k"
  , "detect_k"
  , "detect_mid"
  , "detect_max"
  , "detect_t0"
  , "d"
  , "beta_catch"
  , "catch_eff"
)
  
state_names = c(
  "S" 
  , paste0('E',1:nE)
  , paste0('Ia',1:nIa)
  , paste0('Ip',1:nIp)
  , paste0('Im',1:nIm) 
  , paste0('Is',1:nIs)
  , "I_new_sympt"
  , "Hr", "Hd"
  , "R" , "D" 
  , "D_new"
  , "betat"
)

accum_names <- c("D_new", "I_new_sympt")

# other functions ----
## R0 here just based on the simple transmission rate / recovery rate (weighted by the probability of going into different classes)
covid_R0   <- function (beta0est, fixed_params, sd_strength, prop_S) {
  
## transmission rate
 R <- beta0est * prop_S * sd_strength * (
   
#    (                
## proportion * time in asymptomatic
#      fixed_params["alpha"] * fixed_params["Ca"] * fixed_params["lambda_a"] +                  
## proportion * time in mildly symptomatic
#      (1 - fixed_params["alpha"]) * fixed_params["mu"] * (fixed_params["lambda_p"] + fixed_params["lambda_m"]) +    
## proportion * time in severely symptomatic
#      (1 - fixed_params["alpha"]) * (1 - fixed_params["mu"]) * (fixed_params["lambda_p"] + fixed_params["lambda_s"])      
#  )
   
   fixed_params["alpha"] * fixed_params["Ca"] +
     (1 - fixed_params["alpha"]) * fixed_params["mu"] * 
     ((fixed_params["lambda_p"]/(fixed_params["lambda_p"] + fixed_params["lambda_m"])) * fixed_params["Cp"] +
     (fixed_params["lambda_m"]/(fixed_params["lambda_p"] + fixed_params["lambda_m"])) * fixed_params["Cm"]) +
     (1 - fixed_params["alpha"]) * (1 - fixed_params["mu"]) * 
     ((fixed_params["lambda_p"]/(fixed_params["lambda_p"] + fixed_params["lambda_s"])) * fixed_params["Cp"] + 
     (fixed_params["lambda_s"]/(fixed_params["lambda_p"] + fixed_params["lambda_s"])) * fixed_params["Cs"]) 
  
    ) 
 
 unlist(R)
 
}

logis_func <- function (L, k, mid, day) {
  L / (1 + exp(-k * (day - mid)))
}
beta_func  <- function (beta0, beta_min, sip_prop) {
  beta0 * exp(log(beta_min) * sip_prop)
}



