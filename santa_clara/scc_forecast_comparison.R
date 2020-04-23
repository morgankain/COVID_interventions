needed_packages <- c(
  "pomp"
  , "plyr"
  , "dplyr"
  , "ggplot2"
  , "magrittr"
  , "scales"
  , "lubridate"
  , "tidyr"
  , "foreach"
  , "doParallel"
  , "data.table"
)

lapply(needed_packages, require, character.only = TRUE)

var_params = readRDS("./COVID_interventions/santa_clara/variable_params_Apr1.RDS") %>% 
  filter(beta0est > 0) %>%
  mutate(R0 = covid_R0(beta0est, fixed_params, sd_strength = 1))
# params <- read.csv("./COVID_interventions/santa_clara/params.csv", stringsAsFactors = FALSE) %>%
#   mutate(Value = sapply(Value, function(x) eval(parse(text = x))))
# fixed_params        <- params$Value
# names(fixed_params) <- params$Parameter
fixed_params <- c(
  Ca               = 2/3
  , Cp               = 1
  , Cs               = 1
  , Cm               = 1
  , alpha            = 1/3
  , gamma            = 1/5.2
  , lambda_a         = 1/7
  , lambda_s         = 1/5.76
  , lambda_m         = 1/7
  , lambda_p         = 1/2
  , rho              = 1/14.5 
  , delta            = 0.2
  , mu               = 1-0.044 
)
fixed_params        <- c(fixed_params, N = 1.938e6)

source("./COVID_interventions/santa_clara/scc_pomp_objs.R")

focal.county = "Santa Clara"
deaths     <- fread("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv")
deaths     <- deaths %>% mutate(date = as.Date(date)) %>% filter(county == focal.county) %>% 
  mutate(day = as.numeric(date - as.Date("2019-12-31"))) %>% 
  select(day, date, deaths) %>% 
  mutate(deaths_cum = deaths) %>% 
  mutate(deaths = deaths_cum - lag(deaths_cum)) %>% 
  replace_na(list(deaths = 0)) %>%
  select(-deaths_cum) %>% 
  mutate(forecasted = date > as.Date("2020-04-01"))
sim_end = deaths$date %>% max
# for(i in 1:nrow(var_params)){

model_forecast_fits = data.frame()

for(i in 1:nrow(var_params)){
  sim_start = var_params[i, "sim_start"]
  sim_length = as.numeric(sim_end - sim_start) + 1
  
  interventions <- with(var_params[i, ], {covariate_table(
    day              = as.numeric(sim_start - as.Date("2019-12-31")):as.numeric(sim_end - as.Date("2019-12-31"))
      
    , intervention     = c(rep(0, int_start1 - sim_start)
                           , rep(1, int_length1)
                           , rep(1, sim_end - int_start1 - int_length1 + 1))
    , isolation = rep(0, sim_length)
    , soc_dist_level = c(rep(0, int_start1 - sim_start)
                         , rep(sd_m1, int_length1)
                         , rep(sd_m2, sim_end - int_start1 - int_length1 + 1))
    , iso_severe_level = rep(NA, sim_length)
    , iso_mild_level = rep(NA, sim_length)
    , thresh_H_start   = rep(NA, sim_length)  
    , thresh_H_end     = rep(NA, sim_length)
    , thresh_int_level = rep(NA, sim_length)
    , order            = "constant"
    , times            = "day"
  )})
  
  covid <- pomp(data = deaths, 
                time       = "day"
                , t0         = as.numeric(sim_start - as.Date("2019-12-31"))
                , covar      = interventions
                , rprocess   = euler(sir_step, delta.t = 1/6)
                , rmeasure   = rmeas 
                , dmeasure   = dmeas
                , rinit      = sir_init
                , partrans   = par_trans
                , accumvars  = accum_names
                , paramnames = param_names
                , statenames = state_names)
  
  SEIR.sim <- pomp::simulate(covid
                             , nsim = 200
                             , params = c(fixed_params, 
                                          c(beta0 = var_params[i, "beta0est"], 
                                            E0 = var_params[i, "E0"]))
                             , format       = "d"
                             , include.data = T)
  # print(head(SEIR.sim))
  band = SEIR.sim %>%
    dplyr::mutate(
      date    = as.Date("2019-12-31") + day) %>%
    # filter(date > as.Date("2020-04-01")) %>%
    select(date,.id,deaths) %>%
    mutate(data=.id=="data") %>%
    plyr::ddply(~date+data,plyr::summarize,
                p=c(0.025,0.5,0.975),
                q=quantile(deaths,prob=p,names=FALSE)
    ) %>% 
    mutate(p=plyr::mapvalues(p,from=c(0.025,0.5,0.975),to=c("lo","med","hi")),
           data=plyr::mapvalues(data,from=c(TRUE,FALSE),to=c("data","simulation"))) %>%
    spread(p,q) %>% 
    pivot_wider(id_cols = date, 
                names_from = data, 
                values_from = c("hi", "lo", "med")) %>% 
    select(-hi_data, -lo_data, -med_simulation) %>% 
    mutate(in_CI = (med_data <= hi_simulation & lo_simulation <= med_data)) %>% 
    select(date, in_CI) %>% 
    pivot_wider(names_from = date, values_from = in_CI)
  max_deaths = SEIR.sim %>% filter(.id != "data") %>% 
    arrange(.id, day) %>%
    group_by(.id) %>% 
    summarise(max_day = which.max(deaths) + 30,
              max = max(deaths)) %>% 
    summarise(max_day_med = median(max_day),
              max_med = median(max))
  model_forecast_fits %<>% rbind(data.frame(var_params[i,],
                                            band,
                                            max_deaths))
}

model_forecast_fits_save = model_forecast_fits


test = model_forecast_fits
test %>% mutate(var_set_id = row_number()) %>% 
  pivot_longer(
  cols = X2020.01.31:X2020.04.09,
  names_to = "date",  
  values_to = "in_CI") %>% 
  mutate(date = as.Date(gsub("\\.", "-", gsub("X", "", date)))) %>% 
  mutate(in_CI = in_CI*1) %>% 
  pull(in_CI) %>%
  {left_join(., group_by(., var_set_id) %>%
               filter(date < as.Date("2020-04-02")) %>% 
            summarise(prior_good = all(in_CI)))} %>% 
  group_by(date) %>%
  summarise(avg_in_CI = mean(in_CI)) %>% View
  # filter(prior_good == T) %>%
  ggplot(aes(x = date, y = avg_in_CI)) + 
  geom_line(alpha = 0.2, group  = 1) +
  scale_x_date(labels = date_format("%Y-%b-%d"))


test %>% mutate(var_set_id = row_number()) %>% 
  select(var_set_id, max_med, max_day_med) %>% 
  ggplot() +
  geom_point(aes(x = max_day_med, y = max_med), col = "red") + 
  geom_line(aes(x = day, y = deaths), data = deaths) + 
  geom_vline(xintercept = as.Date("2020-04-01") - as.Date("2019-12-31"))


sim4 <- SEIR.sim %>%
  dplyr::mutate(
    date    = as.Date("2019-12-31") + day) %>%
  # filter(date > as.Date("2020-04-01")) %>%
  select(date,.id,deaths) %>%
  mutate(data=.id=="data") %>%
  plyr::ddply(~date+data,plyr::summarize,
              p=c(0.025,0.5,0.975),
              q=quantile(deaths,prob=p,names=FALSE)
  ) %>% 
  mutate(p=plyr::mapvalues(p,from=c(0.025,0.5,0.975),to=c("lo","med","hi")),
         data=plyr::mapvalues(data,from=c(TRUE,FALSE),to=c("data","simulation"))) %>%
  spread(p,q) %>% 
  mutate(data = factor(data, levels = c("simulation", "data"))) %>%
  ggplot(aes(x = date, y = med, ymin = lo, ymax = hi, group = data, color = data, fill  = data)) + 
  geom_ribbon(alpha = 0.5) + 
  geom_vline(xintercept = as.Date("2020-04-01")) +
  scale_fill_manual(values = c("grey", "blue")) + 
  scale_color_manual(values = c("grey", "blue")) + 
  ylab("Deaths") + 
  ggtitle(paste0("Santa Clara predictions for parameter set ", i, "\n", 
                 paste0(names(var_params[i,c(1, 2, 9, 12)]), " = ", clean_print(var_params[i,c(1, 2, 9, 12)]), collapse =  ", "),
                 "\n",
                 paste0(names(var_params[i,c(3, 5, 7, 6)]), " = ", clean_print(var_params[i,c(3, 5, 7, 6)]), collapse =  ", "))) + 
  theme_bw()

library(ggpubr)
ggarrange(sim1, sim2, sim3, sim4, nrow = 2, ncol = 2)

clean_print = function(df_row){
  classes = sapply(df_row, class)
  char_vec = c()
  for(i in 1:length(classes)){
    if(classes[i] == "numeric" | classes[i] == "integer" | classes[i] == "difftime") {
      char_vec <- c(char_vec, as.character(round(df_row[,i], 2)))
    } else char_vec <- c(char_vec, as.character(df_row[,i])) #(classes[i] == "Date")
  }
  return(char_vec)
}
