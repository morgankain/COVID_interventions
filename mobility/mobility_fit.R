library(pomp)
library(magrittr)
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)

source("./COVID_pomp_mobility.R")
all_data = data.table::fread("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv") 
fixed_params = c(
  Ca = 0.6,
  Cp = 1,
  Cs = 1, 
  Cm = 1,
  alpha = 0.43,
  gamma = 1/3.5,
  lambda_a = 1/7,
  lambda_s = 1/5.5,
  lambda_m = 1/5.5,
  lambda_p = 1/1.5,
  delta = 0.2,
  mu = 0.95,
  rho_d = 1/13.3,
  rho_r = 1/15
)

focal_county = "Santa Clara"
focal_state_abbr = "CA"
focal_state = "California"
start_date = as.Date("2020-01-15")
fixed_params["N"] <- 1927852		

data = all_data %>% 
  filter(county == focal_county & state == focal_state) %>% 
  mutate(deaths = deaths - lag(deaths),
         cases = cases - lag(cases)) %>% 
  select(date, cases, deaths) %>% 
  mutate(day = as.numeric(as.Date(date) - as.Date("2019-12-31"))) %>% 
  filter(as.Date(date) < as.Date("2020-05-10")) %>% 
  filter(!is.na(deaths)) %>% 
  select(date, day, deaths, cases) %>% 
  # 5 day moving median
  mutate(deaths = mutate(.,
                         deaths_lag1 = lag(deaths),
                         # deaths_lag2 = lag(deaths, n = 2),
                         deaths_lead1 = lead(deaths), 
                         # deaths_lead2 = lead(deaths, n =2)
                         ) %>%
           select(contains("deaths")) %>%
           purrr::pmap_dbl(~median(c(...))),
         cases = mutate(.,
                        cases_lag1 = lag(cases),
                        # cases_lag2 = lag(cases, n = 2),
                        cases_lead1 = lead(cases),
                        # cases_lead2 = lead(cases, n = 2)
                        ) %>%
           select(contains("cases")) %>%
           purrr::pmap_dbl(~median(c(...)))) %>%
  filter(!is.na(deaths))

mobility = readRDS("./unfolded_daily_clean.rds") %>% 
  filter(county_name == focal_county & state_abbr == focal_state_abbr)  %>%
  select(datestr, sip_prop) %>% 
  filter(!is.na(sip_prop)) %>%
  mutate(day = as.numeric(datestr - as.Date("2019-12-31"))) %>% 
         # sip_prop = 1 - sip_prop,
         # mobility = sip_prop - min(sip_prop),
         # mobility = mobility / max(mobility)) %>% 
  arrange(day)


# covid_mobility_deaths = pomp(
#   data = data, times = "day",
#   t0 = as.numeric(start_date - as.Date("2019-12-31")),
#   covar = covariate_table(sip_prop = mobility$sip_prop, 
#                           order = "constant",
#                           times = mobility$day),
#   rprocess = euler(sir_step_mobility, delta.t = 1/6),
#   rmeasure = rmeas_deaths_NB, 
#   dmeasure = dmeas_deaths_NB,
#   rinit = sir_init,
#   partrans = par_trans, 
#   accumvars  = accum_names,
#   paramnames = param_names,
#   statenames = state_names
# )
options(warning.length = 4000L)
covid_mobility_multi = pomp(
  data = data, times = "day",
  t0 = as.numeric(start_date - as.Date("2019-12-31")),
  covar = covariate_table(sip_prop = mobility$sip_prop, 
                          order = "constant",
                          times = mobility$day),
  rprocess = euler(sir_step_mobility, delta.t = 1/6),
  rmeasure = rmeas_multi,
  dmeasure = dmeas_multi,
  rinit = sir_init,
  partrans = par_trans, 
  accumvars  = accum_names,
  paramnames = param_names,
  statenames = state_names
)

set.seed(101)
    
library(doParallel)
mif_test = foreach(i = 1:4
                   ) %dopar% {
  library(pomp)
  mf <- mif2(covid_mobility_multi, 
             t0 = as.numeric(start_date - as.Date("2019-12-31")),
             params = c(fixed_params,
                        import_rate = 0.0,
                        beta0 = 1, 
                        beta_min = 0.1, 
                        theta_c = 5,
                        theta_d = 10,
                        E_init = 4,
                        detect_t0 = 10,
                        detect_t1 = 100,
                        detect_max = 0.1
             ),
             Np = 100, 
             Nmif = 600,
             cooling.fraction.50 = 0.75,
             rw.sd = rw.sd(beta0 = 0.02, 
                           beta_min = 0.02, 
                           theta_c = 0.02, 
                           theta_d = 0.02,
                           E_init = ivp(0.02),
                           detect_t0 = 0.02, 
                           detect_t1 = 0.02, 
                           detect_max = 0.02 
             ))
  mf
  
  # replicate(
  #   10, 
  #   mf %>% pfilter(Np=5000) %>% logLik()
  # ) %>%
  #   logmeanexp(se=TRUE) -> ll
  # 
}

# plot traces
ldply(1:4, function(i){
  data.frame(mif_no = i,
             iter = 0:600,
             traces(mif_test[[i]]))
}) %>% 
  select(mif_no, iter, loglik, beta0:detect_max) %>% 
  pivot_longer(loglik:detect_max) %>% 
  ggplot(aes(x = iter, 
             y = value, 
             color = as.factor(mif_no))) +
  geom_line() +
  facet_wrap(~name, scales = "free_y") + 
  theme_bw()

# do a terrible job of pick a best params
best_pars <- ldply(1:4, function(i){
  data.frame(mif_no = i,
             iter = 0:600,
             traces(mif_test[[i]]))
}) %>% 
  arrange(desc(loglik)) %>% 
  magrittr::extract(1,) %>% unlist

test_sim <- pomp::simulate(covid_mobility_multi, nsim = 50, 
               t0 = as.numeric(start_date - as.Date("2019-12-31")),
               params = best_pars,
               format  = "d", 
               include.data = T) %>% 
  mutate(betat = best_pars["beta0"]*
         exp(log(best_pars["beta_min"])*sip_prop),
         Reff = covid_R0(betat, fixed_params, 1, 1)) 

gridExtra::grid.arrange(ggplot(test_sim) + 
                          geom_line(aes(x = day, y = deaths, 
                                        group = .id, color = .id == "data")) + 
                          geom_point(aes(x = day, y = deaths),
                                     data = test_sim %>% filter(.id == "data")) + 
                          theme_bw() + 
                          theme(legend.position = "none"), 
                        ggplot(test_sim) + 
                          geom_line(aes(x = day, y = cases, 
                                        group = .id, color = .id == "data")) + 
                          geom_point(aes(x = day, y = cases),
                                     data = test_sim %>% filter(.id == "data")) + 
                          theme_bw() +
                          theme(legend.position = "none"),
                        ggplot(test_sim) + 
                          geom_line(aes(x = day, y = Reff, 
                                        group = .id, color = .id == "data")) +
                          geom_hline(yintercept = 1, linetype = "dashed") +
                          theme_bw() + 
                          theme(legend.position = "none"),
                        nrow = 3)
          

# plot Reff over the entire observed mobility data
# plot(mobility$day, 
#      covid_R0(coef(mif_test)["beta0"]*
#        exp(log(coef(mif_test)["beta_min"])*mobility$sip_prop),
#        fixed_params, 1, 1),
#      type = "l", 
#      xlab = "day", 
#      ylab = "Reff")


# testing functional form
# plot(mobility$day, mobility$sip_prop)
# plot(x = (0:200)/100,
#      y = exp(log(.01)*(0:200)/100),
#      type = "l")
# points(mobility$sip_prop, 
#        exp(log(0.1)*mobility$sip_prop),
#        col = "blue")
# abline(v = 1, col = "red", lty = "dashed")
# abline(h = 0.8, col = "red", lty = "dashed")