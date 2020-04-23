params = readRDS("./COVID_interventions/santa_clara/variable_params_Apr1.Rds")
params %<>% mutate(R0 = beta0est * (1/3) * (2/3) * (7) + beta0est * (2/3) * (0.5 + 5.76) * (1 - 0.956) + beta0est * (2/3) * (0.5 + 7) * 0.956) %>% 
  filter(beta0est > 0)
pairs(~R0 + E0 + sim_start + int_start1 + sd_m1 + sd_m2, data = params %>% select(-c(int_start2, paramset, int_length1)))

i = 1
scc_deaths <- fread("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv") %>% 
  mutate(date = as.Date(date)) %>%
  filter(county == "Santa Clara")  %>% 
  mutate(day = as.numeric(date - params[i, ]$sim_start)) %>% 
  select(day, date, deaths) %>% 
  rename(deaths_cum = deaths) %>% 
  mutate(deaths = deaths_cum - lag(deaths_cum)) %>% 
  replace_na(list(deaths = 0)) %>%
  select(-deaths_cum)
scc_deaths <- rbind(
  data.frame(
    day    = seq(1:(min(scc_deaths$day) - 1))
    , date   = as.Date(seq(1:(min(scc_deaths$day) - 1)), origin = params[i, ]$sim_start)
    , deaths = 0
  )
  , scc_deaths
)
sim_length <- 500
inf_iso = FALSE
intervention.forecast <- with(params[i, ], {
  
  covariate_table(
    day              = 1:sim_length
    
    , intervention     = c(
      # No intervention until intervention start time
      rep(0, int_start1 - sim_start)                   
      # Intervention style 1
      , rep(1, int_length1)
      # Intervention style 2
      , rep(1, int_length2)
      # Post intervention close
      , {
        if (!inf_iso) {
          rep(0, sim_length - (int_start2 - sim_start) - int_length2)
        } else {
          rep(1, sim_length - (int_start2 - sim_start) - int_length2)  
        }
      }
    ) 
    
    , isolation = { 
      
      if (!inf_iso) {
        
        rep(0, sim_length)
        
      } else {
        
        c(
          rep(0, iso_start - sim_start)  
          , rep(1, sim_length - (iso_start - sim_start))
        )
        
      }}
    
    , iso_severe_level = rep(0, sim_length)      # % of contats that severe cases maintain
    , iso_mild_level   = rep(0.05, sim_length)   # % of contats that mild cases maintain
    
    , soc_dist_level = {
      
      if (!inf_iso) {
        
        c(                       # intensity of the social distancing interventions
          rep(sd_m1, int_start2 - sim_start)
          ## slightly odd way to do this, but should work
          , rep(sd_m2, sim_length - (int_start2 - sim_start)))
        
      } else {
        
        c(                       # intensity of the social distancing interventions
          rep(sd_m1, int_start2 - sim_start)
          ## slightly odd way to do this, but should work
          , rep(sd_m2, iso_start - int_start2)
          ## if infected isolation is going on, increase the background contact rate
          , rep(sd_m1, sim_length - (iso_start - sim_start))
        )
        
      }}
    
    , thresh_H_start   = rep(2, sim_length)  
    , thresh_H_end     = rep(10, sim_length)
    
    ## No reason to have a second parameter here, just use the same val that the user picks for social distancing
    , thresh_int_level = c(                       # level of social distancing implemented with the threshhold intervention
      rep(sd_m1, int_start2 - sim_start)
      ## slightly odd way to do this, but should work
      , rep(sd_m2, sim_length - (int_start2 - sim_start))
    )
    
    , order            = "constant"
    , times            = "day"
  )
  
})

covid.fitting <- scc_deaths %>%
  pomp(
    time       = "day"
    , t0         = 1
    , covar      = intervention.forecast
    , rprocess   = euler(sir_step, delta.t = 1/6)
    , rmeasure   = rmeas 
    , dmeasure   = dmeas
    , rinit      = sir_init
    , partrans   = par_trans
    , accumvars  = accum_names
    , paramnames = param_names
    , statenames = state_names
  ) 
SEIR.sim <- do.call(
  pomp::simulate
  , list(
    object         = covid.fitting
    # , times        = 1:100 #intervention.forecast@times
    , params       = c(fixed_params, c(beta0 = params[i, "beta0est"], E0 = params[i, ]$E0))
    , nsim         = 100
    , format       = "d"
    , include.data = T
    , seed         = 1001)) %>% {
      rbind(.,
            group_by(., day) %>%
              select(-.id) %>%
              summarise_all(median) %>%
              mutate(.id = "median"))
    } 
# plotting with CIs
CI_data <- do.call(
  pomp::simulate
  , list(
    object         = covid.fitting
    # , times        = 1:100 #intervention.forecast@times
    , params       = c(fixed_params, c(beta0 = params[i, "beta0est"], E0 = params[i, ]$E0))
    , nsim         = 1000
    , format       = "d"
    , include.data = T
    , seed         = 1001)) %>%  
  dplyr::mutate(
    date    = params[i,"sim_start"] + day - 1) %>%
  select(date,.id,deaths) %>%
  mutate(data=.id=="data") %>%
  plyr::ddply(~date+data,plyr::summarize,
              p=c(0.025,0.5,0.975),
              q=quantile(deaths,prob=p,names=FALSE)
  ) %>% 
  mutate(p=plyr::mapvalues(p,from=c(0.025,0.5,0.975),to=c("lo","med","hi")),
         data=plyr::mapvalues(data,from=c(TRUE,FALSE),to=c("data","simulation"))) %>%
  spread(p,q) %>% 
  ggplot(aes(x=date,y=med,color=data,fill=data,ymin=lo,ymax=hi))+
  geom_ribbon(color = NA, alpha=0.5)+
  geom_line() +
  geom_point(aes(size = ifelse(data == "data" & date > as.Date("2020-03-06"), T, F)), 
             show.legend = FALSE) +
  xlab("Date") + 
  ylab("Daily Fatalities") +
  scale_color_manual(values = c(wes_palette("Royal1")[2], wes_palette("Royal1")[1])) + 
  scale_fill_manual(values = c(wes_palette("Royal1")[2], wes_palette("Royal1")[1])) + 
  scale_size_manual(values = c(0,2)) + 
  scale_x_date(labels = date_format("%Y-%b-%d")) +
  theme(legend.title = element_blank()) +
  guides(data=FALSE)

# pdf("./model_fit_Apr1.pdf", width = 12, height = 7.5)
print(CI_data)
# dev.off()

SEIR.sim <- SEIR.sim %>%
  dplyr::mutate(
    date    = params[i,"sim_start"] + day - 1
    , total_I = Is + Im + Ia + Ip
  )

SEIR.sim %>%
  ggplot() +
  geom_line(
    data = (SEIR.sim %>% filter(.id != "data" & .id != "median"))
    , aes(
      x      = date
      , y      = deaths
      , group  = .id
      , colour = inf_iso
    ), alpha = 0.4, colour = "grey60") + 
  geom_line(data = (SEIR.sim %>% filter(.id == "median"))
            , aes(  
              x     = date
              , y     = deaths)
            , colour = "black",
            , lwd = 1.5) +
  geom_line(data = (SEIR.sim %>% filter(.id == "data"))
            , aes(  
              x     = date
              , y     = deaths)
            , colour = "blue",
            , lwd = 1.5) +
  scale_x_date(labels = date_format("%Y-%b")) +
  scale_y_log10() + 
  xlab("Date") + 
  ylab("Daily Fatalities") +
  theme_bw() + 
  guides(color = F)



# SEIR.sim %>%
#   mutate(type = case_when(.id == "data" ~ "data",
#                           .id == "median" ~ "median", 
#                           T ~ "simulation"),
#          type = factor(type, levels = c("simulation", "median", "data"))) %>%
#   arrange(type) %>% View
#   ggplot() +
#   geom_line(
#     aes(
#       x      = date
#       , y      = deaths
#       , group  = .id
#       , colour = type
#       , size = as.factor(type == "simulation")
#     )) + 
#   scale_size_manual(values=c(1.5, 0.5))+
#   scale_color_manual(values = c("grey80", "black", "blue")) +
#   # geom_line(data = (SEIR.sim %>% filter(.id == "data"))
#   #           , aes(  
#   #             x     = date
#   #             , y     = deaths)
#   #           , colour = "blue",
#   #           , lwd = 1.5) +
#   scale_x_date(labels = date_format("%Y-%b")) +
#   scale_y_log10() + 
#   xlab("Date") + 
#   ylab("Daily Fatalities") +
#   theme_bw() 
