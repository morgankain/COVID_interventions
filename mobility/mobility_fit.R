source("./COVID_pomp_mobility.R")
all_data = fread("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv") 
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

focal_county = "Miami-Dade"
focal_state_abbr = "FL"
focal_state = "Florida"
start_date = as.Date("2020-02-22")
fixed_params["N"] <- 2716940		

data = all_data %>% 
  filter(county == focal_county & state == focal_state) %>% 
  mutate(deaths = deaths - lag(deaths),
         cases = cases - lag(cases)) %>% 
  select(date, cases, deaths) %>% 
  mutate(day = as.numeric(as.Date(date) - as.Date("2019-12-31"))) %>% 
  filter(as.Date(date) < as.Date("2020-05-14")) %>% 
  filter(!is.na(deaths)) %>% 
  select(date, day, deaths) %>% 
  # three day moving median
  mutate(deaths = mutate(., 
                         deaths_lag = lag(deaths),
                         deaths_lead = lead(deaths)) %>% 
           select(contains("deaths")) %>%
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


covid_mobility = pomp(
  data = data, times = "day",
  t0 = as.numeric(start_date - as.Date("2019-12-31")),
  covar = covariate_table(sip_prop = mobility$sip_prop, 
                          order = "constant",
                          times = mobility$day),
  rprocess = euler(sir_step_mobility, delta.t = 1/6),
  rmeasure = rmeas_deaths_NB, 
  dmeasure = dmeas_deaths_NB,
  rinit = sir_init,
  partrans = par_trans, 
  accumvars  = accum_names,
  paramnames = param_names,
  statenames = state_names
)

set.seed(101)
mif_test = mif2(covid_mobility, 
                t0 = as.numeric(start_date - as.Date("2019-12-31")),
                params = c(fixed_params,
                           theta = 5,
                           import_rate = 0.0,
                           beta0 = 1, 
                           beta_min = 0.01, 
                           E_init = 4),
                Np = 100, 
                Nmif = 400,
                cooling.fraction.50 = 0.5,
                rw.sd = rw.sd(beta0 = 0.02, 
                              beta_min = 0.02, 
                              theta = 0.005,
                              E_init = 0.02))
plot(mif_test)


pomp::simulate(covid_mobility, nsim = 11, 
               t0 = as.numeric(start_date - as.Date("2019-12-31")),
               params = coef(mif_test),
               format  = "d", 
               include.data = T) %>% 
  mutate(betat = coef(mif_test)["beta0"]*
         exp(log(coef(mif_test)["beta_min"])*sip_prop),
         Reff = covid_R0(betat, fixed_params, 1, 1)) %>%
  ggplot() + 
  geom_line(aes(x = day, y = deaths, #Reff,
                group = .id, color = .id == "data")) + 
  # geom_point(aes(x = day, y = cum_deaths),
             # data = data %>% mutate(cum_deaths = cumsum(deaths))) #+
  facet_wrap(~.id)


# plot Reff over the entire observed mobility data
plot(mobility$day, 
     covid_R0(coef(mif_test)["beta0"]*
       exp(log(coef(mif_test)["beta_min"])*mobility$sip_prop),
       fixed_params, 1, 1),
     type = "l", 
     xlab = "day", 
     ylab = "Reff")


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