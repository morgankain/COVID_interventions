source("ggplot_theme.R")

# Figure 1: model fits ---

counties_list = {list(
  SC = list(
    focal.county = "Santa Clara",
    focal.state = "California",
    focal.state_abbr = "CA",
    rds.name = "./output/Santa Clara_TRUE_FALSE_0_2020-06-12_cont_temp_SC_ind_theta.rds",
    con_theta = F
  ),
  LA = list(
    focal.county = "Los Angeles",
    focal.state = "California",
    focal.state_abbr = "CA",
    rds.name = "output/Los Angeles_TRUE_FALSE_0_2020-06-12_cont_temp_LA_ind_theta.rds",
    con_theta = F
  ),
  KC = list(
    focal.county = "King",
    focal.state = "Washington",
    focal.state_abbr = "WA",
    rds.name = "output/King_TRUE_FALSE_0_2020-06-12_cont_temp_KC_ind_theta.rds",
    con_theta = F
  ),
  FC = list(
    focal.county = "Fulton",
    focal.state = "Georgia",
    focal.state_abbr = "GA",
    rds.name = "output/Fulton_TRUE_FALSE_0_2020-06-12_cont_temp_FC_ind_theta.rds",
    con_theta = F
  ),
  MD = list(
    focal.county = "Miami-Dade",
    focal.state = "Florida",
    focal.state_abbr = "FL",
    rds.name = "output/Miami-Dade_TRUE_FALSE_0_2020-06-12_cont_temp_MD_ind_theta.rds",
    con_theta = F
  ),
  CC = list(
    focal.county = "Contra Costa",
    focal.state = "California",
    focal.state_abbr = "CA",
    rds.name = "output/Contra_Costa_TRUE_FALSE_0_2020-06-12_cont_temp_ind_theta.rds",
    con_theta = F
  )
)}[-6]
int_vars <- list(
  counter.factual    = FALSE,
  int.movement       = "post",
  int.type           = "none",
  int.init           = "2020-06-08",
  int.end            = "2020-08-01", ## Doesnt matter, ignored for all int.type != "tail"
  sim_title          = "Reality"                                        
)
# set parameters
source("COVID_simulate_cont_params.R")
loglik.max         <- FALSE
loglik.thresh      <- 2
# plot_vars          <- c("cases", "deaths")
ci.stoc            <- 0.025
ci.epidemic        <- T
nsim               <- 500
plot_vars <- c("cases", "deaths")

# run the simulations for all locations
fig1_data <- plyr::adply(1:length(counties_list), 1, 
      function(i){
        with(c(counties_list[[i]], int_vars), {
        source("./COVID_simulate_cont.R", local = T)
        SEIR.sim.f.ci %<>% 
          full_join(county.data %>%
                      arrange(date) %>%
                      mutate(D = cumsum(ifelse(is.na(deaths), 0, deaths))*
                               ifelse(is.na(deaths), NA, 1)) %>%
                      select(date, any_of(plot_vars)) %>%
                      pivot_longer(any_of(plot_vars), values_to = "data")) %>% 
          rbind(Reff %>% mutate(name = "Reff", mid = Reff, lwr = NA, upr = NA,
                                data = NA) %>% 
                  select(-Reff)) %>% 
          mutate(intervention = sim_title,
                 county = focal.county,
                 state = focal.state_abbr)
        return(SEIR.sim.f.ci)
})}, .id = NULL)

# plot the fits and data
fig1_colors = c("#FF0000",
                "#00A08A",
                "#F2AD00",
                "#F98400",
                "#5BBCD6")[c(1,4,5,2,3)]
fig1_data %>% 
  filter(county != "Los Angeles") %>%
  mutate(county = paste0(county, " County,", state)) %>%
  filter(date < as.Date("2020-06-08")) %>%
  filter(date >= as.Date("2020-02-10")) %>%
  filter(name != "Reff") %>%
  group_by(county) %>% 
  mutate(nparams = 0.5/length(unique(paramset))) %>% 
  ggplot(aes(x = date, y = mid, ymin = lwr, ymax = upr, 
             fill  = county, color = county, 
             group = interaction(county,paramset))) +
  geom_ribbon(aes(alpha = I(nparams)), colour = NA) +
  geom_line() +
  geom_point(aes(x = date, y = data), 
             color = "black", size = 0.75) + 
  scale_y_continuous(trans = "sqrt") + 
  scale_fill_manual(guide = F, values = fig1_colors) +
  scale_color_manual(guide = F, values = fig1_colors) +
  facet_grid(name ~ county, scales = "free_y", switch = "y")  +
  ylab("") + 
  theme(
    strip.background = element_blank()
    , strip.placement = "outside"
    , strip.text = element_text(size = 16)) 

fig1_data %>% 
  filter(date < as.Date("2020-06-12")) %>%
  filter(date >= as.Date("2020-02-10")) %>%
  filter(name == "Reff") %>% 
  ggplot(aes(x = date, y = mid, color = county,
             group = interaction(paramset, county))) + 
  geom_line() + 
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(values = fig1_colors) + 
  ylab("Reff")
  
