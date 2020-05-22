library(pomp)
library(magrittr)
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(doParallel)

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

focal_county = "Los Angeles"
focal_state_abbr = "CA"
focal_state = "California"
start_date = as.Date("2020-01-25")
fixed_params["N"] <- 10039107		

data = all_data %>% 
  filter(county == focal_county & state == focal_state) %>% 
  mutate(date = as.Date(date),
         day = as.numeric(date - as.Date("2019-12-31"))) %>% 
  select(date, day, cases, deaths) %>% 
  # add on zeros and NAs in the beginning 
  {rbind(data.frame(date = start_date:(min(pull(., date))-1),
                    day = as.numeric(start_date - as.Date("2019-12-31")):
                      (as.numeric(min(pull(., date)) - as.Date("2019-12-31")) - 1),
                    cases = 0,
                    deaths = NA),
         .)} %>%
  mutate(deaths = deaths - lag(deaths),
         cases = cases - lag(cases)) %>%
  filter(date < as.Date("2020-05-10")) %>% 
  select(date, day, deaths, cases) %>% 
  # 5 day moving median
  mutate(deaths = mutate(.,
                         deaths_lag1 = lag(deaths, n = 1),
                         deaths_lag2 = lag(deaths, n = 2),
                         deaths_lead1 = lead(deaths, n = 1),
                         deaths_lead2 = lead(deaths, n =2)
                         ) %>%
           select(contains("deaths")) %>%
           purrr::pmap_dbl(~median(c(...))),
         cases = mutate(.,
                        cases_lag1 = lag(cases, n = 1),
                        cases_lag2 = lag(cases, n = 2),
                        cases_lead1 = lead(cases, n = 1),
                        cases_lead2 = lead(cases, n = 2)
                        ) %>%
           select(contains("cases")) %>%
           purrr::pmap_dbl(~median(c(...)))) %>%
  filter(!is.na(cases))

mobility = readRDS("./unfolded_daily_clean.rds") %>% 
  filter(county_name == focal_county & state_abbr == focal_state_abbr)  %>%
  select(datestr, sip_prop) %>% 
  filter(!is.na(sip_prop)) %>%
  mutate(day = as.numeric(datestr - as.Date("2019-12-31"))) %>% 
  arrange(day)

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
mif_test = foreach(i = 1:4
                   ) %dopar% {
  library(pomp)
  mf <- mif2(covid_mobility_multi, 
             t0 = as.numeric(start_date - as.Date("2019-12-31")),
             # params = c(fixed_params,
             #            import_rate = 0.0,
             #            # beta0 = 1.1,
             #            # beta_min = 0.02,
             #            # theta_c = 2.5,
             #            # theta_d = 2.5,
             #            # E_init = 5,
             #            # detect_t0 = 60,
             #            # detect_t1 = 200,
             #            # detect_max = 0.7
             #            beta0 = 1,
             #            beta_min = 0.1,
             #            theta_c = 5,
             #            theta_d = 10,
             #            E_init = 4,
             #            detect_t0 = 10,
             #            detect_t1 = 200,
             #            detect_max = 0.8
             # ),
             params = best_pars,
             Np = 100, 
             Nmif = 1000,
             cooling.fraction.50 = 0.5,
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
}

# calculate logliks 
ldply(mif_test, function(mf) {
  replicate(
    10,
    mf %>% pfilter(Np=5000) %>% logLik()
  ) %>%
    logmeanexp(se=TRUE)
  })

# plot traces
ldply(1:length(mif_test), function(i){
  data.frame(mif_no = i,
             iter = 0:1000,
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

# do a terrible job of pick a best params... need to do this a better way
best_pars <- ldply(1:length(mif_test), function(i){
  data.frame(#mif_no = i,
             #iter = 0:600,
             traces(mif_test[[i]]))
}) %>% 
  arrange(desc(loglik)) %>% 
  magrittr::extract(1,) %>% unlist

best_pars <- coef(mif_test[[2]])

# simulate and plot
test_sim <- pomp::simulate(covid_mobility_multi, nsim = 50, 
               t0 = as.numeric(start_date - as.Date("2019-12-31")),
               params = best_pars,
               format  = "d", 
               include.data = T) %>% 
  mutate(betat = best_pars["beta0"]*
         exp(log(best_pars["beta_min"])*sip_prop),
         Reff = covid_R0(betat, fixed_params, 1, 1)) 

gridExtra::grid.arrange(
  # plot deaths
  ggplot(test_sim) + 
    geom_line(aes(x = day, y = deaths, 
                  group = .id, color = .id == "data"),
              alpha = 0.6) + 
    geom_point(aes(x = day, y = deaths),
               data = test_sim %>% filter(.id == "data")) + 
    scale_y_continuous(trans = "log10") + 
    scale_color_manual(values = wesanderson::wes_palette("Darjeeling1")[2:3]) + 
    theme_bw() + 
    theme(legend.position = "none"), 
  # plot cases
  ggplot(test_sim) + 
    geom_line(aes(x = day, y = cases, 
                  group = .id, color = .id == "data"),
              alpha = 0.6) + 
    geom_point(aes(x = day, y = cases),
               data = test_sim %>% filter(.id == "data")) + 
    scale_y_continuous(trans = "log10") +
    scale_color_manual(values = wesanderson::wes_palette("Darjeeling1")[2:3]) + 
    theme_bw() +
    theme(legend.position = "none"),
  # plot Reff
  ggplot(test_sim) + 
    geom_line(aes(x = day, y = Reff, 
                  group = .id)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme_bw() + 
    theme(legend.position = "none"),
  # calculate and plot detection rate
  test_sim %>% filter(.id == "data") %>%
    select(day) %>%
    mutate(
      detect = ifelse(day < best_pars["detect_t0"], 
                      0,
                      ifelse(day < best_pars["detect_t0"] + best_pars["detect_t1"],
                             best_pars["detect_max"]/best_pars["detect_t1"]*(day - best_pars["detect_t0"]) ,
                             best_pars["detect_max"]))) %>% 
    ggplot(aes(x = day, y = detect)) + 
    geom_line() + 
    theme_bw(),
  nrow = 4)

# time = 0:500
# plot(time, 
#      1 / (1 + exp(-0.02 * (time - 250))))
          
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