source("needed_packages.R")
source("ggplot_theme.R")

counties_list = {list(
  SC = list(rds.name = "output/Santa Clara_0_2020-06-24_cont_final.Rds"),
  LA = list(rds.name = "output/Los Angeles_0_2020-06-24_cont_final.Rds"),
  KC = list(rds.name = "output/King_0_2020-06-24_cont_final.Rds")
  # FC = list(rds.name = "output/Fulton_0_2020-06-24_cont_round2.Rds"),
  # MD = list(rds.name = "output/Miami-Dade_0_2020-06-24_cont_round2.Rds"),
  # MC = list(rds.name = "output/Mercer_0_2020-06-24_cont_round2.Rds"),
  # BC = list(rds.name = "output/Bucks_0_2020-06-24_cont_round2.Rds"),
  # CC = list(rds.name = "output/Chester_0_2020-06-24_cont_round2.Rds")
)}

#####
## Figure 1: model fits ----
#####

int_vars <- list(
   counter.factual    = FALSE
 , int.movement       = c("post", "post")
 , int.type           = c("none")
 , int.init           = c("2020-07-01")
 , sim_end            = "2020-07-01"
 , sim_title          = "Reality"                                        
)

sim_vars = list(loglik.max     =  T,
                loglik.num     = 5,
                ci.stoc        = 0.025,
                ci.epidemic    = T,
                nsim           = 500,
                plot_vars      = c("cases", "deaths"),
                plot.median    = T)

# run the simulations for all locations
fig1_data <- plyr::ldply(counties_list, 
      function(i){
        source("COVID_simulate_cont_params.R", local = T)
        with(c(i, int_vars, sim_vars), {
        source("./COVID_simulate_cont.R", local = T)
          SEIR.sim.f.ci %<>% data.frame() %>%
          full_join(pomp_data %>%
                      arrange(day) %>%
                      mutate(D = cumsum(ifelse(is.na(deaths), 0, deaths))*
                               ifelse(is.na(deaths), NA, 1),
                             date = day + date_origin) %>%
                      select(date, any_of(plot_vars)) %>%
                      pivot_longer(any_of(plot_vars), values_to = "data")) %>%
            as.data.frame() %>%
            rbind(., (Reff %>% mutate(name = "Reff", lwr = NA, mid = Reff, upr = NA, data = NA) %>% select(-Reff) %>% dplyr::select(date, paramset, everything()))) %>%
            rbind(., (detect %>% mutate(name = "detect", lwr = NA, mid = detect, upr = NA, data = NA) %>% select(-detect) %>% dplyr::select(date, paramset, everything()))) %>%
            mutate(intervention = sim_title,
                 county = variable_params[1, "county"] %>% as.character,
                 state = variable_params[1, "state"] %>% as.character)
        return(SEIR.sim.f.ci)
})}, .id = NULL)

# plot the fits and data
fig1_colors = c("#FF0000",
                "#00A08A",
                "#F2AD00",
                "#F98400",
                "#5BBCD6",
                "black", "blue", "gold")#[c(1,4,5,2,3)]

# fig1_data$name <- sapply(fig1_data$name, simpleCap)

## LA with same axis
{
fig1_data %>% 
 # filter(county != "Los Angeles") %>%
  mutate(county = paste0(county, " County, ", state)) %>%
  filter(date < as.Date("2020-06-08")) %>%
  # filter(date >= as.Date("2020-02-10")) %>%
  # filter(!(name %in% c("Reff", "Detect"))) %>%
  group_by(county) %>% 
  mutate(nparams = 0.5/length(unique(paramset))) %>% 
  ggplot(aes(x = date, y = mid, ymin = lwr, ymax = upr, 
             fill  = county, color = county, 
             group = interaction(county,paramset))) +
  geom_ribbon(aes(alpha = I(nparams)), colour = NA) +
  geom_line() +
  geom_point(aes(x = date, y = data), 
             color = "black", size = 0.75) + 
  scale_y_continuous(trans = "log10") + 
  scale_fill_manual(guide = F, values = fig1_colors) +
  scale_color_manual(guide = F, values = fig1_colors) +
  facet_grid(name ~ county, scales = "free", switch = "y")  +
  ylab("") + 
  theme(
    strip.background = element_blank()
    , strip.placement = "outside"
    , strip.text = element_text(size = 16)
    , axis.text.x = element_text(size = 12)) +
  xlab("Date")
}

## LA with a different axis
# {
# fig1.1 <- {fig1_data %>% 
#   filter(county != "Los Angeles") %>%
#   mutate(county = paste0(county, " County, ", state)) %>%
#   filter(date < as.Date("2020-06-08")) %>%
#   filter(date >= as.Date("2020-02-10")) %>%
#   filter(!(name %in% c("Reff", "Detect"))) %>%
#   group_by(county) %>% 
#   mutate(nparams = 0.5/length(unique(paramset))) %>% 
#   ggplot(aes(x = date, y = mid, ymin = lwr, ymax = upr, 
#              fill  = county, color = county, 
#              group = interaction(county,paramset))) +
#   geom_ribbon(aes(alpha = I(nparams)), colour = NA) +
#   geom_line() +
#   geom_point(aes(x = date, y = data), 
#              color = "black", size = 0.75) + 
#   scale_y_continuous(trans = "sqrt") + 
#   scale_fill_manual(guide = F, values = fig1_colors) +
#   scale_color_manual(guide = F, values = fig1_colors) +
#   facet_grid(name ~ county, scales = "free_y", switch = "y")  +
#   ylab("") + 
#   theme(
#     strip.background = element_blank()
#     , strip.placement = "outside"
#     , strip.text.x = element_text(size = 11)
#     , strip.text.y = element_text(size = 16)
#     , axis.text.x = element_text(size = 12)
#     , axis.title.y = element_text(margin = margin(l = 1))
#     , plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm")
#     ) + 
#     xlab("") +
#     labs(tag = "C.") +
#     theme(plot.tag.position = "topleft"
#           , plot.tag = element_text(face = "bold" 
#                                    , margin = margin(b = -16, 
#                                                      r = -45)))}
# 
# fig1.2 <- {fig1_data %>% 
#   filter(county == "Los Angeles") %>%
#   mutate(county = paste0(county, " County, ", state)) %>%
#   filter(date < as.Date("2020-06-08")) %>%
#   filter(date >= as.Date("2020-02-10")) %>%
#   filter(!(name %in% c("Reff", "Detect"))) %>%
#   group_by(county) %>% 
#   mutate(nparams = 0.5/length(unique(paramset))) %>% 
#   ggplot(aes(x = date, y = mid, ymin = lwr, ymax = upr, 
#              fill  = county, color = county, 
#              group = interaction(county,paramset))) +
#   geom_ribbon(aes(alpha = I(nparams)), colour = NA) +
#   geom_line() +
#   geom_point(aes(x = date, y = data), 
#              color = "black", size = 0.75) + 
#   scale_y_continuous(trans = "sqrt", position = "right") + 
#   scale_fill_manual(guide = F, values = fig1_colors[5]) +
#   scale_color_manual(guide = F, values = fig1_colors[5]) +
#   facet_grid(name ~ county, scales = "free_y", switch = "y")  +
#   ylab("") + 
#   theme(
#     strip.background = element_blank()
#     , strip.placement = "outside"
#     , strip.text.y = element_blank()
#     , strip.text.x = element_text(size = 11)
#     , axis.text.x = element_text(size = 12)
#     , plot.margin = unit(c(0.25,0.25,0.25,0), "cm")) +
#   xlab("") 
#   # xlab("Date")
#   }
# 
# fig1.3 <- {fig1_data %>% 
#   filter(date < as.Date("2020-06-08")) %>%
#   filter(date >= as.Date("2020-02-10")) %>%
#   filter(name == "Reff") %>%
#   # ugly way to calculate 7 day smoothing 
#   arrange(county, paramset, date) %>%
#   group_by(county, paramset) %>%
#   mutate(
#     lag3 = lag(mid, 3),
#     lag2 = lag(mid, 2),
#     lag1 = lag(mid, 1),
#     lead1 = lead(mid, 1),
#     lead2 = lead(mid, 2),
#     lead3 = lead(mid, 3)
#   ) %>% 
#   rowwise() %>% 
#   mutate(mid = mean(c(lag1, lag2, lag3, mid, lead1, lead2, lead3), na.rm = T)) %>% 
#   select(-starts_with("lag"), -starts_with("lead")) %>% 
#   mutate(paramset = as.character(paramset)) %>%
#   group_by(county, state, date) %>% 
#   {rbind(., 
#          dplyr::summarise(., 
#                           med = mean(mid), 
#                           lwr = min(mid),
#                           upr = max(mid),
#                           .groups = "drop") %>% 
#            rename(mid = med) %>% 
#            mutate(paramset = "summary",
#                   name = "Detect",
#                   intervention = "Reality",
#                   data = NA))} %>% 
#   mutate(width = ifelse(paramset == "summary", 1.5, 0.25),
#          alpha = ifelse(paramset == "summary", 1, 0.2)) %>%
#   ggplot(aes(x = date, y = mid, ymin = lwr, ymax = upr, 
#              fill  = county, color = county, 
#              group = interaction(county,paramset))) +
#   geom_ribbon(alpha = 0.25, color = NA) +
#   geom_line(aes(size = I(width), alpha = I(alpha))) +
#   scale_color_manual(guide= F, values = fig1_colors[c(1, 2, 5, 3, 4)]) +
#   scale_fill_manual(guide= F, values = fig1_colors[c(1, 2, 5, 3, 4)]) +
#   geom_hline(yintercept  = 1, linetype = "dashed") + 
#   ylab(expression("R"[e])) + 
#   xlab("") +
#   scale_y_continuous(position = "left") + 
#   ggtitle("A. Reproduction number")+
#   theme(plot.title.position = "plot",
#         plot.title = element_text(face = "bold"))
#   }
# 
# fig1.4 <-  {fig1_data %>%  
#   filter(date < as.Date("2020-06-08")) %>%
#   filter(date >= as.Date("2020-02-10")) %>%
#   filter(name == "Detect") %>%
#   mutate(paramset = as.character(paramset)) %>%
#   group_by(county, state, date) %>% 
#   {rbind(., 
#          dplyr::summarise(., 
#                           med = mean(mid), 
#                           lwr = min(mid),
#                           upr = max(mid),
#                           .groups = "drop") %>% 
#            rename(mid = med) %>% 
#            mutate(paramset = "summary",
#                   name = "Detect",
#                   intervention = "Reality",
#                   data = NA))} %>% 
#   mutate(width = ifelse(paramset == "summary", 1.5, 0.25),
#          alpha = ifelse(paramset == "summary", 1, 0.2)) %>%
#   ggplot(aes(x = date, y = mid, ymin = lwr, ymax = upr, 
#              group = interaction(county, paramset),
#              fill  = county, color = county)) +
#   geom_ribbon(alpha = 0.25, color = NA) +
#   geom_line(aes(size = I(width), alpha = I(alpha))) +
#   scale_color_manual(guide= F, values = fig1_colors[c(1, 2, 5, 3, 4)]) +
#   scale_fill_manual(guide= F, values = fig1_colors[c(1, 2, 5, 3, 4)]) +
#   # geom_hline(yintercept  = 1, linetype = "dashed") + 
#   ylab("Detection probability") + 
#   xlab("") + 
#   scale_y_continuous(labels = scales::percent_format(accuracy = 1),
#                      position = "right") +
#   ggtitle("B. Syptomatic detection")+
#   theme(plot.title.position = "plot",
#         plot.title = element_text(face = "bold"))
#   }
#   
# # gridExtra::grid.arrange(fig1.1, fig1.2, ncol = 2, widths = c(4, 1.2))
# fig1 <- gridExtra::arrangeGrob(fig1.1, fig1.2, fig1.3,fig1.4,  
#                        layout_matrix = matrix(c(3, 4, 4, 1, 1, 2), 
#                                               byrow  = T, nrow = 2),
#                        widths = c(2.6, 1.3, 1.3), heights = c(1.5, 2.5)) 
# gridExtra::grid.arrange(fig1)
# # ggsave("figures/Manuscript2/Figure1.pdf", 
# #        fig1, 
# #        device = "pdf",
# #        width = 12,
# #        height = 8)
# }

#####
## Figure 2: Interventions ----
#####

int_vars <- {
list(
  lift_all = list(
  counter.factual    = FALSE,
  int.movement       = c("post", "pre"),
  int.type           = "none",
  int.init           = "2020-07-01",
  sim_end            = "2020-07-31",
  sim_title          = "Lift"                                    
),
  continue_SIP = list(
  counter.factual    = FALSE,
  int.movement       = c("post", "post"),
  int.type           = "none",
  int.init           = "2020-07-01",
  sim_end            = "2020-07-31",
  sim_title          = "Continue SIP"                                        
),
  inf_iso = list(
  counter.factual    = FALSE,
  int.movement       = c("post", "mid"),
  int.type           = "inf_iso",
  int.init           = "2020-07-01",
  sim_end            = "2020-07-31",
  sim_title          = "Isolate",
  iso_mild_level     = 0.1,
  iso_severe_level   = 0.1
),
  tail_chop1 = list(
  counter.factual    = FALSE,
  int.movement       = c("post", "mid"),
  int.type           = "tail",
  int.init           = "2020-07-01",
  sim_end            = "2020-07-31",
  sim_title          = "Superspreading_0.01_1",
  int.beta_catch     = 0.01,
  int.beta_catch_type= "pct",
  int.catch_eff      = 1.0
),
  tail_chop2 = list(
  counter.factual    = FALSE,
  int.movement       = c("post", "mid"),
  int.type           = "tail",
  int.init           = "2020-07-01",
  sim_end            = "2020-07-31",
  sim_title          = "Superspreading_0.02_1",
  int.beta_catch     = 0.02,
  int.beta_catch_type= "pct",
  int.catch_eff      = 1.0
),
  tail_chop3 = list(
  counter.factual    = FALSE,
  int.movement       = c("post", "mid"),
  int.type           = "tail",
  int.init           = "2020-07-01",
  sim_end            = "2020-07-31",
  sim_title          = "Superspreading_0.05_1",
  int.beta_catch     = 0.05,
  int.beta_catch_type= "pct",
  int.catch_eff      = 1.0
)
)
}


sim_vars = list(
loglik.max     = TRUE,
plot_vars      = c("cases", "deaths"),
ci.stoc        = 0.025,
ci.epidemic    = T,
nsim           = 200)

fig2_data <- ldply(int_vars, function(j) {
  ldply(counties_list[c("SC", "KC")], function(i) {
    source("COVID_simulate_cont_params.R", local = T)
    with(c(i, j, sim_vars), {
      print(int.beta_catch)
      source("./COVID_simulate_cont.R", local = T)
      SEIR.sim.f.ci %<>% 
        full_join(pomp_data %>%
                    arrange(day) %>%
                    mutate(D = cumsum(ifelse(is.na(deaths), 0, deaths))*
                             ifelse(is.na(deaths), NA, 1),
                           date = day + date_origin) %>%
                    select(date, any_of(plot_vars)) %>%
                    pivot_longer(any_of(plot_vars), values_to = "data")) %>%
        mutate(intervention = sim_title,
               county = variable_params[1, "county"] %>% as.character,
               state = variable_params[1, "state"] %>% as.character) 
      return(SEIR.sim.f.ci)
      
    }
    )
  }
  )
}
, .id = NULL)

fig2_colors <- c("#D67236", "dodgerblue4", "#0b775e", "magenta4")#, "magenta2")

fig2_data %>% 
  filter(date >= as.Date("2020-02-10")) %>%
  filter(intervention %in% c("Continue SIP",
                             "Isolate",
                             "Lift",
                             "Superspreading_0.01_1"#,
                             # "Superspreading_0.02_1"#,
                             # "Superspreading_0.05_1"
                             )) %>% 
  # mutate(intervention = mapvalues(intervention,
  #                                 from = c("Continue SIP",
  #                                          "Isolate",
  #                                          "Lift",
  #                                          "Superspreading_0.01_1"),
  #                                 to = c(
  #                                   "Continue Shelter in Place"
  #                                   , "Infected Isolation"
  #                                   , "Lift all Interventions"
  #                                   , "Superspreading Averted\n(1% Removed with 100% Effectiveness)"
  #                                 ))) %>%
  {ggplot(aes(x = date, y = mid, 
              ymin = lwr, ymax = upr, 
              fill  = intervention, color = intervention),
          data = .) +
      facet_grid(name ~ county, scales = "free_y", switch = "y")  +
      geom_ribbon(data = (filter(., date >= as.Date("2020-07-01")))
                  , alpha = 0.50, colour = NA) +
  geom_line(data = (filter(., date >= as.Date("2020-07-01")))) +
  geom_ribbon(data = (filter(., date <= as.Date("2020-07-01")
                             , intervention == "Lift"))
  , alpha = 0.50, colour = NA, fill = "black") +
  geom_line(data = (filter(., date <= as.Date("2020-07-01")
   , intervention == "Lift"
    # , intervention == "Superspreading Averted_A"
    ))
  , colour = "black") +
  geom_vline(xintercept = as.Date("2020-07-01"), linetype = "dashed", lwd = 0.5) + 
  geom_point(aes(x = date, y = data), 
             color = "black", size = 1) + 
  scale_y_continuous(trans = "log10", labels = trans_format('log10', math_format(10^.x))) +
  scale_fill_manual(values = fig2_colors, name = "Intervention") +
  scale_color_manual(values = fig2_colors, name = "Intervention") +
  ylab("") + 
  theme(
    strip.background = element_blank()
    , legend.position  = c(0.100, 0.925)
    , legend.background = element_blank()
    , strip.placement = "outside"
    , strip.text = element_text(size = 16)
    , axis.text.x = element_text(size = 14)) +
  xlab("Date")}

#####
## Figure 3: SIP and truncation combinations ----
#####  

# check calculation for truncated gamma distribution mean
## truncated gamma distribution 
# rtgamma = function(n, shape, scale, lower, upper, limits = "pct"){
#   if(limits == "pct"){
#     upper <- qgamma(upper, shape, scale = scale)
#     lower <- qgamma(lower, shape, scale = scale)
#   }
#   out <- rgamma(n, shape, scale = scale)
#   trunc <- which(out < lower | out > upper)
#   while(length(trunc) > 0){
#     out[trunc] <- rgamma(length(trunc), shape, scale = scale)
#     trunc <- which(out < lower | out > upper)
#   }
#   return(out)
# }
# 
# test <- adply(seq(1, 5, by = 0.1), 1, 
#               function(beta0){
#                 k = 0.16
#                 upper = 3
#                 sim_est = rtgamma(10000, 
#                                   shape = k,
#                                   scale = beta0/k, 
#                                   0, upper, limits = "abs") %>% mean
#                 alg_est = beta0/k*pgamma(upper*k/beta0, k+1)*gamma(k+1)/(pgamma(upper*k/beta0, k)*gamma(k))
#                 return(data.frame(beta0 = beta0, sim = sim_est, alg = alg_est))
#               })
# plot(test$sim, test$alg)
# abline(a = 0, b = 1, col = "red", lty = "dashed")

# here's function to calculate SIP prop needed to maintain R=1 with a given tail truncation
# works with vector catch_eff and beta_catch. output is beta_catch X catch_eff in dimensions
sip_trunc_combns = function(beta_catch, beta_catch_type = "pct", 
                            catch_eff, 
                            k, beta0, beta_min, d, dt,
                            desired_R){
  shp = k*dt/d
  scl = beta0/k
  if(beta_catch_type == "pct"){
    upper = qgamma(1 - beta_catch, shape = shp, scale = scl)
  } else{
    upper = beta_catch
  }
  
  # calculate expected value of truncated gamma dist when truncation with 100% efficacy
  #### NEED TO REDO THIS CALCULATION WITH NEW DISTRIBUTION
  E_trunc <- beta0/k*pgamma(upper/scl, shp+1)*gamma(shp+1)/(pgamma(upper/scl, shp)*gamma(shp))
  # print(E_trunc)
  # expected value is weighted sum of truncated and untruncated distributions depdending on efficacy
  out <- -(log(outer(E_trunc, catch_eff, "*") + 
                 matrix(rep(outer(beta0*dt/d, 1 - catch_eff,"*"), length(beta_catch)), 
                        byrow = T, nrow = length(beta_catch))) + log(d) - log(dt) - log(desired_R))/log(beta_min)
  return(out)
}

beta_catch_vals <- seq(0, 0.005, by = 0.00005)
catch_eff_vals <- seq(0.5, 1, by = 0.25)
loglik.max     <- F
loglik.num     <- 2

fig3_data <- adply(1:length(counties_list), 1, function(i){
  prev.fit = readRDS(counties_list[[i]]$rds.name)
  if("mifs_local_v2" %in% names(prev.fit)){ # get second round of mifs if it exists
    variable_params  <- prev.fit[["mifs_local_v2"]] # otherwise get the first round
  } else{
    variable_params  <- prev.fit[["mifs_local"]] # otherwise get the first round
  }
  if (loglik.max) {
    variable_params <- variable_params %>% 
      filter(log_lik != 0) %>% 
      filter(log_lik == max(log_lik))
  } else {
    if (is.na(loglik.num)) {
      variable_params <- variable_params %>% 
        filter(log_lik != 0) %>% 
        filter(log_lik > (max(log_lik) - loglik.thresh))  
    } else {
      variable_params <- variable_params %>% 
        filter(log_lik != 0) %>% 
        arrange(desc(log_lik)) %>%
        slice(1:loglik.num)
    }
  }
  
  covid_mobility   <- prev.fit$covid_mobility # get pomp object
  pomp_data        <- data.frame(covid_mobility) %>% 
    full_join(as.data.frame(cbind(day = covid_mobility@covar@times, 
                                  t(covid_mobility@covar@table)))) %>%
    arrange(day)
  dt               <- covid_mobility@rprocess@delta.t
  
  sip_start <- arrange(pomp_data, day) %>% head(10) %>% pull(sip_prop) %>% mean
  sip_max <- pomp_data %>% pull(sip_prop) %>% max
  
  # now loop over all the top fits
  adply(1:nrow(variable_params), 1, function(j){
    params_row = unlist(variable_params[j,])
    data.frame(
      beta_catch_pct = beta_catch_vals,
      sip = sip_trunc_combns(beta_catch_vals, 
                             "pct",
                             catch_eff_vals, 
                             k = params_row["beta0_k"],
                             beta0 = params_row["beta0"],
                             beta_min = params_row["beta_min"],
                             d = params_row["d"],
                             dt = dt,
                             desired_R = 1.0) %>% 
        magrittr::set_colnames(paste0("catch_eff_", catch_eff_vals)),
      county = as.character(variable_params[j,"county"]),
      paramset = j) %>% return
  }, .id = NULL) %>%
    # make catch efficiency a column and pivot data longer
    pivot_longer(starts_with("sip"), names_to = "catch_eff", 
                 names_prefix = "sip.catch_eff_", values_to = "sip") %>%
    # add two rows with data of sip observed
    rbind(data.frame(
      beta_catch_pct = 0,
      county = variable_params[1, "county"],
      paramset = "data",
      sip = c(sip_start, sip_max),
      catch_eff = 0
    )) %>% 
    return
}, .id = NULL) 

fig3_data %>% 
  group_by(county, catch_eff, beta_catch_pct) %>% 
  {rbind(.,
         filter(., paramset != "data") %>% 
           summarise(., 
                     lwr = min(sip),
                     upr = max(sip),
                     sip = mean(sip)) %>% 
           mutate(paramset = "summary"))} %>% 
  ungroup() %>%
  mutate(., catch_eff_label = factor(paste0(as.numeric(catch_eff)*100, "% efficiency"),
                                     levels = paste0(as.numeric(unique(pull(., catch_eff)))*100, "% efficiency"))) %>% 
  {ggplot(data = filter(., paramset == "summary"),
          mapping = aes(x = beta_catch_pct, y = sip, 
                        ymin = lwr, ymax = upr,
                        color = county, 
                        fill = county, 
                        group = interaction(paramset,county))) + 
      geom_ribbon(alpha = 0.4, color = NA) +
      geom_line(alpha = 0.75, size = 1.5) + 
      geom_point(aes(shape = type, size = I(size)), 
                 position = position_dodge(width = 0.0005),
                 alpha = 0.75, 
                 data = filter(., paramset == "data") %>% 
                   select(-catch_eff, -catch_eff_label) %>%
                   group_by(county) %>% 
                   mutate(type = ifelse(sip == max(sip), "max", "start")) %>%
                   {right_join(., 
                               expand.grid(county = pull(., county) %>% unique, 
                                           catch_eff = catch_eff_vals))}  %>%
                   mutate(., catch_eff_label = factor(paste0(as.numeric(catch_eff)*100, "% efficiency"),
                                                      levels = paste0(as.numeric(unique(pull(., catch_eff)))*100, "% efficiency"))) %>% 
                   ungroup %>% 
                   mutate(size = ifelse(type == "start", 2, 2))
      ) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 2)) +
      scale_x_continuous(labels = scales::percent_format(accuracy = 0.1)) +
      xlab("Percentile of superspreading averted") +
      ylab("Proportion sheltering-in-place") + 
      scale_color_manual(values = fig1_colors) +
      scale_fill_manual(values = fig1_colors) +
      scale_shape_manual(guide = F, values = c(19, 17)) + #95)) + 
      facet_wrap(~catch_eff_label)}

#####
## Figure 4: Epidemic rebound when rare ----
#####  

## Have to run figure 3 in advance, which is annoying and likely should be fixed

counties_list <- {
  list(
  SC =   list(
    focal.county = "Santa Clara",
    focal.state = "California",
    focal.state_abbr = "CA",
    rds.name = "output/Santa Clara_independent_theta_200_2020-06-15.rds",
    con_theta = F
  )#,
#  KC = list(
#    focal.county = "King",
#    focal.state = "Washington",
#    focal.state_abbr = "WA",
#    rds.name = "output/King_independent_theta_200_2020-06-16.rds",
#    con_theta = F
#  )
)
}

int_vars <- {
list(
  lift_all = list(
  counter.factual    = FALSE,
  int.movement       = "mid",
  int.type           = "tail",
  int.init           = "2020-06-08",
  int.end            = "2020-08-01", ## Doesnt matter, ignored for all int.type != "tail"
  sim_title          = "Continue Shelter in Place",
  thresh_inf.val     = 1,
  int.catch_eff      = 1.0,
  int.beta_catch       = 0.10, #0.10,
  int.beta0_sigma      = 0.16,
  int.beta0_sigma_post = 0.16, #0.16,
  int.beta_catch_post  = 0.05,
  int.catch_eff_post   = 0.0,
    ## Needed SIP scaling with "mid" to get R to the desired R
  int.sip_prop_scaling = check_R0(
    beta0est     = variable_params$beta0est
  , beta_min     = variable_params$beta_min
  , fixed_params = unlist(c(fixed_params, variable_params))
  , sd_strength  = 1
  , prop_S       = 1
  , desired_R    = 2.0) / 0.3463429
 ),
  continue_SIP = list(
  counter.factual    = FALSE,
  int.movement       = "mid",
  int.type           = "tail",
  int.init           = "2020-06-08",
  int.end            = "2020-08-01", ## Doesnt matter, ignored for all int.type != "tail"
  sim_title          = "Minimial Tail Cut",
  thresh_inf.val     = 1,
  int.catch_eff      = 1.0,
  int.beta_catch     = 0.10,
  int.beta0_sigma    = 0.16,
  int.beta0_sigma_post = 0.16,
  int.beta_catch_post  = 0.05,
  int.catch_eff_post   = 1.0,
  int.sip_prop_scaling = (fig3_data %>% filter(county == "Santa Clara", beta_catch_pct == 0.95, catch_eff == 1.0))$sip / 0.3463429                                       
),
  inf_iso = list(
  counter.factual    = FALSE,
  int.movement       = "mid",
  int.type           = "tail",
  int.init           = "2020-06-08",
  int.end            = "2020-08-01", ## Doesnt matter, ignored for all int.type != "tail"
  sim_title          = "Large Tail Cut",
  thresh_inf.val     = 1,
  int.catch_eff      = 1.0,
  int.beta_catch     = 0.10,
  int.beta0_sigma    = 0.16,
  int.beta0_sigma_post = 0.16,
  int.beta_catch_post  = 0.14,
  int.catch_eff_post   = 0.8,
  int.sip_prop_scaling = (fig3_data %>% filter(county == "Santa Clara", beta_catch_pct == 0.86, catch_eff == 0.8))$sip / 0.3463429
)
)}

source("epidemic_rebound/gamma_rebound_pomp.R")
source("ggplot_theme.R")
source("epidemic_rebound/gamma_rebound_params.R")

nsim           <- 100
thresh_inf.val <- 10
sim_length     <- 300     ## How many days to run the simulation
plot_vars      <- c("cases", "deaths", "I", "betat")

fig4_data <- adply(1:length(int_vars), 1, 
      function(j) {
        adply(1:length(counties_list), 1,
      function(i) {
        with(c(counties_list[[i]], int_vars[[j]]), {
        source("epidemic_rebound/gamma_rebound.R", local = T)
        SEIR.sim.f.t %<>% 
          full_join(county.data %>%
                      arrange(date) %>%
                      mutate(D = cumsum(ifelse(is.na(deaths), 0, deaths))*
                               ifelse(is.na(deaths), NA, 1)) %>%
                      select(date, any_of(plot_vars)) %>%
                      pivot_longer(any_of(plot_vars), values_to = "data")) %>% 
          mutate(intervention = sim_title,
                 county = focal.county,
                 state = focal.state_abbr) 
     #   return(SEIR.sim.f.ci)
         return(SEIR.sim.f.t)
        
        }
        )
      }
        )
      }
        , .id = NULL)

#fig4_data$name        <- sapply(fig4_data$name, simpleCap)
fig4_data$intervention <- factor(fig4_data$intervention, levels = unique(fig4_data$intervention))

fig4_colors <- c("#D67236", "dodgerblue4", "#0b775e", "magenta4")
fig4_data   <- fig4_data %>% filter(state == "CA")

check_date <- (fig4_data %>% filter(name == "I", .id == "median", intervention == "Continue Shelter in Place", date > "2020-05-01") %>% 
    filter(value == min(value)) %>% 
    filter(date == min(date)))$date

## For just I plot
fig4_data <- fig4_data %>% filter(date >= as.Date("2020-07-10")) %>% filter(name == "I")

fig4_data %>%
  ggplot(aes(x = date, y = value
 #  , ymin = lwr, ymax = upr, fill  = intervention
    , color = intervention
    , group = .id)) +
#  geom_ribbon(data = (fig4_data %>% filter(date >= check_date))
#   , alpha = 0.50, colour = NA) +
#  geom_line(data = (fig4_data %>% filter(date >= check_date))
#   ) + 
  geom_line(data = (fig4_data %>% filter(date >= check_date, .id != "median"))
    , aes(group = interaction(.id, intervention))
    , lwd = 0.25, alpha = 0.40
   ) +
  geom_line(data = (fig4_data %>% filter(date >= check_date, .id == "median"))
    , aes(group = interaction(.id, intervention))
   , lwd = 1.5) +
#  geom_ribbon(data = (fig4_data %>% filter(date <= check_date, intervention == "Continue Shelter in Place"))
#  , alpha = 0.50, colour = NA, fill = "black") +
#  geom_line(data = (fig4_data %>% filter(date <= check_date, intervention == "Continue Shelter in Place"))
#  , colour = "black") +
  geom_line(data = (fig4_data %>% filter(date <= check_date, intervention == "Continue Shelter in Place", .id != "median"))
    , lwd = 0.25, alpha = 0.40, colour = "black"
   ) +
  geom_line(data = (fig4_data %>% filter(date <= check_date, intervention == "Continue Shelter in Place", .id == "median"))
    , colour = "black"
    , lwd = 1.5) +
  geom_vline(xintercept = check_date, linetype = "dashed", lwd = 0.5) + 
  geom_point(aes(x = date, y = data), 
             color = "black", size = 1) + 
  scale_y_continuous(trans = "log10") +
  scale_fill_manual(values = fig4_colors, name = "Intervention") +
  scale_color_manual(values = fig4_colors, name = "Intervention") +
  facet_grid(name ~ county, scales = "free_y", switch = "y")  +
  ylab("") + 
  theme(
    strip.background = element_blank()
    , legend.position  = c(0.850, 0.885)
    , legend.background = element_blank()
    , strip.placement = "outside"
    , strip.text = element_text(size = 16)
    , axis.text.x = element_text(size = 14)) +
  xlab("Date")

droplevels(fig4_data) %>% 
  filter(state == "CA", .id != "median") %>% 
  filter(name == "cases") %>% 
  filter(day > 290) %>% 
  group_by(.id, intervention) %>%
  summarize(new_cases = sum(value)) %>% 
  filter(new_cases == 0) %>% 
  group_by(intervention) %>%
  summarize(prop_extinct = n() / nsim)

droplevels(fig4_data) %>% filter(state == "CA", .id != "median") %>% filter(name == "I") %>% filter(day == max(day)) %>% 
  group_by(.id, intervention) %>%
#  summarize(new_cases = sum(value)) %>% 
  group_by(intervention) %>%
  summarize(
      lwr = quantile(value, 0.025)
    , est = quantile(value, 0.500) 
    , upr = quantile(value, 0.975))

fig4_data %>%
  ggplot(aes(x = date
    # , y = value
    , y = mid
    , ymin = lwr, ymax = upr, fill  = intervention
    , color = intervention
   # , group = .id
    )) +
  geom_ribbon(data = (fig4_data %>% filter(date >= check_date))
   , alpha = 0.50, colour = NA) +
  geom_line(data = (fig4_data %>% filter(date >= check_date))
   ) + 
#  geom_line(data = (fig4_data %>% filter(date >= check_date, .id != "median"))
#    , aes(group = interaction(.id, intervention))
#    , lwd = 0.25, alpha = 0.40
#   ) +
#  geom_line(data = (fig4_data %>% filter(date >= check_date, .id == "median"))
#    , aes(group = interaction(.id, intervention))
#   , lwd = 1.5) +
  geom_ribbon(data = (fig4_data %>% filter(date <= check_date, intervention == "Continue Shelter in Place"))
  , alpha = 0.50, colour = NA, fill = "black") +
  geom_line(data = (fig4_data %>% filter(date <= check_date, intervention == "Continue Shelter in Place"))
  , colour = "black") +
#  geom_line(data = (fig4_data %>% filter(date <= check_date, intervention == "Continue Shelter in Place", .id != "median"))
#    , lwd = 0.25, alpha = 0.40, colour = "black"
#   ) +
#  geom_line(data = (fig4_data %>% filter(date <= check_date, intervention == "Continue Shelter in Place", .id == "median"))
#    , colour = "black"
#    , lwd = 1.5) +
  geom_vline(xintercept = check_date, linetype = "dashed", lwd = 0.5) + 
# geom_point(aes(x = date, y = data), color = "black", size = 1) + 
#  scale_y_continuous(trans = "log10") +
  scale_fill_manual(values = fig4_colors, name = "Intervention") +
  scale_color_manual(values = fig4_colors, name = "Intervention") +
 # facet_grid(name ~ county, scales = "free_y", switch = "y")  +
  ylab("") + 
  theme(
    strip.background = element_blank()
    , legend.position  = c(0.250, 0.885)
    , legend.background = element_blank()
    , strip.placement = "outside"
    , strip.text = element_text(size = 16)
    , axis.text.x = element_text(size = 14)) +
  xlab("Date")
