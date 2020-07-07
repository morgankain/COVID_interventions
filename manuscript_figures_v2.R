# functions and variables needed for all figures
source("needed_packages.R")
source("ggplot_theme.R")


# this is the list of counties used for figures 1 - 3
counties_list = {list(
   SC = list(rds.name = "output/Santa Clara_0_2020-06-25_cont_round2.Rds")
 , LA = list(rds.name = "output/Los Angeles_0_2020-06-25_cont_round2.Rds")
 , KC = list(rds.name = "output/King_0_2020-06-25_cont_round2.Rds")
 , FC = list(rds.name = "output/Fulton_0_2020-06-25_cont_round2.Rds")
 , MD = list(rds.name = "output/Miami-Dade_0_2020-06-25_cont_round2.Rds")
#,  MC = list(rds.name = "output/Mercer_0_2020-06-25_cont_round2.Rds")
#,  BC = list(rds.name = "output/Bucks_0_2020-06-25_cont_round2.Rds")
#,  CC = list(rds.name = "output/Chester_0_2020-06-25_cont_round2.Rds")
)}

## truncated gamma distribution simulator
rtgamma = function(n, shape, scale, lower, upper, limits = "pct"){
  if(limits == "pct"){
    upper <- qgamma(upper, shape, scale = scale)
    lower <- qgamma(lower, shape, scale = scale)
  }
  out <- rgamma(n, shape, scale = scale)
  trunc <- which(out < lower | out > upper)
  while(length(trunc) > 0){
    out[trunc] <- rgamma(length(trunc), shape, scale = scale)
    trunc <- which(out < lower | out > upper)
  }
  return(out)
}


# here's function to calculate SIP prop needed to maintain specified R with a given tail truncation
# works with vector catch_eff and beta_catch. output is beta_catch X catch_eff in dimensions
sip_trunc_combns = function(beta_catch, 
                            beta_catch_type = "pct", 
                            catch_eff, 
                            #k, beta0, beta_min, d, 
                            var_params,
                            dt,
                            desired_R){
  k = var_params["beta0_k"] %>% unlist
  beta0 = var_params["beta0"] %>% unlist
  beta_min = var_params["beta_min"] %>% unlist
  d = var_params["d"] %>% unlist
  inf =  var_params["alpha"] * var_params["Ca"] +
    (1 - var_params["alpha"]) * var_params["mu"] * 
    ((var_params["lambda_p"]/(var_params["lambda_p"] + var_params["lambda_m"])) * var_params["Cp"] +
       (var_params["lambda_m"]/(var_params["lambda_p"] + var_params["lambda_m"])) * var_params["Cm"]) +
    (1 - var_params["alpha"]) * (1 - var_params["mu"]) * 
    ((var_params["lambda_p"]/(var_params["lambda_p"] + var_params["lambda_s"])) * var_params["Cp"] + 
       (var_params["lambda_s"]/(var_params["lambda_p"] + var_params["lambda_s"])) * var_params["Cs"]) %>% 
    unlist
  
  shp = k*dt/d
  scl = beta0/k
  # print(shp)
  # print(scl)
  # print(beta_catch)
  if(beta_catch_type == "pct"){
    upper = qgamma(1 - beta_catch, shape = shp, scale = scl)
  } else{
    upper = beta_catch
  }
  
  # calculate expected value of truncated gamma dist when truncation with 100% efficacy
  E_trunc <- beta0/k*pgamma(upper/scl, shp+1)*gamma(shp+1)/(pgamma(upper/scl, shp)*gamma(shp))
  
  # expected value is weighted sum of truncated and untruncated distributions depdending on efficacy
  out <- (log(outer(E_trunc, catch_eff, "*") + outer(rep(shp*scl, length(E_trunc)), 1-catch_eff, "*")) - 
            log(desired_R) - log(dt) + log(d) + log(inf))/-log(beta_min)
  return(out)
}


fig1_colors = c(
    "#FF0000"
  , "#00A08A"
  , "#F2AD00"
  , "dodgerblue3"
  , "maroon3"
  , "#F98400"
  , "#5BBCD6"
  , "black"
)

####
## Figure 1: model fits ----
####

int_vars <- list(
   counter.factual    = FALSE
 , int.movement       = c("post", "post")
 , int.type           = c("none")
 , int.init           = c("2020-07-01")
 , sim_end            = "2020-07-01"
 , sim_title          = "Reality"                                        
)

sim_vars <- list(
   loglik.max      = F
 , loglik.num      = 10
 # , loglik.thresh   = 
 , ci.stoc         = 0.025
 , ci.epidemic     = T
 , ci.epidemic_cut = 500
 , nsim            = 500
 , plot_vars       = c("cases", "deaths")
 , plot.median     = T
  )

# run the simulations for all locations
fig1_data <- plyr::ldply(counties_list, 
      function(i){
        source("COVID_simulate_cont_params.R", local = T)
        with(c(i, int_vars, sim_vars), {
        source("./COVID_simulate_cont.R", local = T)
          SEIR.sim.f.ci %<>% data.frame() %>%
          full_join(pomp_data %>%
                      arrange(day) %>%
                      mutate(D = cumsum(ifelse(is.na(deaths), 0, deaths)) *
                               ifelse(is.na(deaths), NA, 1),
                             date = day + date_origin) %>%
                      select(date, any_of(plot_vars)) %>%
                      pivot_longer(any_of(plot_vars), values_to = "data")) %>%
            as.data.frame() %>%
            ungroup %>%
            rbind(., (Reff %>% 
                        transmute(date = date, paramset = paramset, name = "Reff", 
                                  lwr = NA, mid = Reff, upr = NA, data = NA))) %>%
            # rbind(., (detect %>% 
            #             transmute(date = date, paramset = paramset, name = "detect",
            #                       lwr = NA, mid = detect, upr = NA, data = NA))) %>%
            mutate(intervention = sim_title,
                 county = variable_params[1, "county"] %>% as.character,
                 state  = variable_params[1, "state"] %>% as.character)
        return(SEIR.sim.f.ci)
})}, .id = NULL)

# plot the fits and data
# fig1_data$name <- sapply(fig1_data$name, simpleCap)

# LA with same axis
# {
# fig1_data %>%
#  # filter(county != "Los Angeles") %>%
#   mutate(county = paste0(county, " County, ", state)) %>%
#   filter(date < as.Date("2020-06-08")) %>%
#   # filter(date >= as.Date("2020-02-10")) %>%
#   # filter(!(name %in% c("Reff", "Detect"))) %>%
#   group_by(county) %>%
#   mutate(nparams = 0.5/length(unique(paramset))) %>%
#   ggplot(aes(x = date, y = mid, ymin = lwr, ymax = upr,
#              fill  = county, color = county,
#              group = interaction(county,paramset))) +
#   geom_ribbon(aes(alpha = I(nparams)), colour = NA) +
#   geom_line() +
#   geom_point(aes(x = date, y = data),
#              color = "black", size = 0.75) +
#   scale_y_continuous(trans = "log10") +
#   scale_fill_manual(guide = F, values = fig1_colors) +
#   scale_color_manual(guide = F, values = fig1_colors) +
#   facet_grid(name ~ county, scales = "free", switch = "y")  +
#   ylab("") +
#   theme(
#     strip.background = element_blank()
#     , strip.placement = "outside"
#     , strip.text = element_text(size = 16)
#     , axis.text.x = element_text(size = 12)) +
#   xlab("Date")
# }

## LA with different axis
source("manuscript_figures_v2_1b.R")

####
## Figure 2: Interventions ----
####

int_vars <- {
list(
  lift_all = list(
  counter.factual    = FALSE
, int.movement       = c("post", "post")
, int.type           = "none"
, int.init           = "2020-07-01"
, sim_end            = "2020-08-31"
, sim_title          = "Continue SIP"                                    
),
  continue_SIP = list(
  counter.factual    = FALSE
, int.movement       = c("post", "pre")
, int.type           = "none"
, int.init           = "2020-07-01"
, sim_end            = "2020-08-31"
, sim_title          = "Lift"                                        
),
  inf_iso = list(
  counter.factual    = FALSE
, int.movement       = c("post", "mid")
, int.type           = "inf_iso"
, int.init           = "2020-07-01"
, sim_end            = "2020-08-31"
, sim_title          = "Isolate"
, iso_mild_level     = 0.1
, iso_severe_level   = 0.1
),
  tail_chop1 = list(
  counter.factual     = FALSE
, int.movement        = c("post", "mid")
, int.type            = "tail"
, int.init            = "2020-07-01"
, sim_end             = "2020-08-31"
, sim_title           = "Superspreading_05%Tail_75%Eff"
, int.beta_catch      = 0.005
, int.beta_catch_type = "pct"
, int.catch_eff       = 0.75
),
  tail_chop2 = list(
  counter.factual     = FALSE
, int.movement        = c("post", "mid")
, int.type            = "tail"
, int.init            = "2020-07-01"
, sim_end             = "2020-08-31"
, sim_title           = "Superspreading_1%Tail_50%Eff"
, int.beta_catch      = 0.01
, int.beta_catch_type = "pct"
, int.catch_eff       = 0.50
),
  tail_chop3 = list(
  counter.factual     = FALSE
, int.movement        = c("post", "mid")
, int.type            = "tail"
, int.init            = "2020-07-01"
, sim_end             = "2020-08-31"
, sim_title           = "Superspreading_1%Tail_100%Eff"
, int.beta_catch      = 0.003
, int.beta_catch_type = "pct"
, int.catch_eff       = 1.0
)
)
}

sim_vars <- list(
 loglik.max       = TRUE
, plot_vars       = c("cases", "deaths")
, ci.stoc         = 0.025
, ci.epidemic     = T
, ci.epidemic_cut = 100
, nsim            = 300
, plot.median     = F
  )

fig2_data <- ldply(int_vars, function(j) {
  ldply(counties_list[c("KC", "LA")], function (i) {
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

fig2_colors <- c("#D67236", "dodgerblue4", "#0b775e", "magenta4", "red3")
fig2_data$name <- sapply(fig2_data$name, simpleCap)

# fig2_data <- readRDS("fig2_data.Rds")

fig2_data %>% 
  filter(date >= as.Date("2020-02-10")) %>%
  filter(intervention %in% c(
    "Continue SIP"
  , "Isolate"
  , "Lift"
#  , "Superspreading_05%Tail_75%Eff"
 , "Superspreading_1%Tail_50%Eff"
 , "Superspreading_1%Tail_100%Eff"
                             )) %>% 
   mutate(intervention = mapvalues(intervention,
                                   from = c("Continue SIP",
                                            "Isolate",
                                            "Lift",
                                    #  , "Superspreading_05%Tail_75%Eff"
   "Superspreading_1%Tail_50%Eff"
 , "Superspreading_1%Tail_100%Eff"
                                     ),
                                    to  = c(
                                       "Continue shelter in place"
                                     , "Infected isolation"
                                     , "Lift all interventions"
                                     , "Superspreading Averted(Top 1% of\ndistribution removed with 50% effectiveness)"
                                     , "Superspreading Averted(Top .3% of\ndistribution removed with 100% effectiveness)"
                                      #  , "Superspreading Averted(Top .5% of\ndistribution removed with 75% effectiveness)"
                                   ))) %>%
   mutate(county = mapvalues(county,
                                   from = c("King",
                                            "Los Angeles"
                                     ),
                                    to  = c(
                                       "Seattle, WA"
                                     , "Los Angeles, CA"
                                   ))) %>% 
  mutate(name = mapvalues(name, 
    from = c("Cases", "Deaths")
  , to   = c("Daily Cases", "Daily Deaths")
  )) %>%
   mutate(intervention = factor(intervention, levels = c(
     "Lift all interventions"
   , "Continue shelter in place"
   , "Infected isolation"
  # , "Superspreading Averted(Top 1% of\ndistribution removed with 75% effectiveness)"
   , "Superspreading Averted(Top 1% of\ndistribution removed with 50% effectiveness)"
   , "Superspreading Averted(Top .3% of\ndistribution removed with 100% effectiveness)"
     ))) %>% {
    ggplot(aes(x = date, y = mid, 
              ymin = lwr, ymax = upr, 
              fill  = intervention, color = intervention),
          data = .) +
      facet_grid(name ~ county, scales = "free_y", switch = "y")  +
      geom_ribbon(data = (filter(., date >= as.Date("2020-07-01")))
                  , alpha = 0.50, colour = NA) +
  geom_line(data = (filter(., date >= as.Date("2020-07-01")))) +
  geom_ribbon(data = (filter(., date <= as.Date("2020-07-01")
                             , intervention == "Lift all interventions"))
  , alpha = 0.50, colour = NA, fill = "black") +
  geom_line(data = (filter(., date <= as.Date("2020-07-01")
   , intervention == "Lift all interventions"
    ))
  , colour = "black") +
  geom_vline(xintercept = as.Date("2020-07-01"), linetype = "dashed", lwd = 0.5) + 
  geom_point(aes(x = date, y = data), 
             color = "black", size = 1) + 
  scale_y_continuous(trans = "pseudo_log", breaks = c(1, 10, 100, 1000, 10000)) +
      #+ #, labels = trans_format('log10', math_format(10^.x))) +
  scale_fill_manual(values = fig2_colors[c(2,1,3,4,5)]
    , name = "Intervention") +
  scale_color_manual(values = fig2_colors[c(2,1,3,4,5)]
    , name = "Intervention") +
  ylab("") + 
  theme(
    strip.background = element_blank()
    , legend.position  = c(0.670, 0.905)
    , legend.text = element_text(size = 12)
    , legend.background = element_blank()
    , strip.placement = "outside"
    , strip.text = element_text(size = 16)
    , axis.text.x = element_text(size = 14)) +
  xlab("Date")}

####
## Figure 3: SIP and truncation combinations ----
####
beta_catch_vals <- seq(0, 0.01, by = 0.00005)
catch_eff_vals  <- seq(0.5, 1, by = 0.25)
loglik.max      <- F
loglik.num      <- 10

fig3_curve_data <- adply(1:length(counties_list), 1, function(i){
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
  sip_max   <- pomp_data %>% pull(sip_prop) %>% max
  
  # now loop over all the top fits
  adply(1:nrow(variable_params), 1, function(j){
    params_row = unlist(variable_params[j,])
    data.frame(
      beta_catch_pct = beta_catch_vals,
      sip = sip_trunc_combns(beta_catch_vals, 
                             "pct",
                             catch_eff_vals, 
                             var_params = params_row,
                             # k = params_row["beta0_k"],
                             # beta0 = params_row["beta0"],
                             # beta_min = params_row["beta_min"],
                             # d = params_row["d"],
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

fig3_curves <- fig3_curve_data %>% 
  group_by(county, catch_eff, beta_catch_pct) %>% 
  {rbind(.,
         filter(., paramset != "data") %>% 
           summarise(., 
                     lwr = min(sip),
                     upr = max(sip),
                     sip = mean(sip)) %>% 
           mutate(paramset = "summary"))} %>% 
  ungroup() %>%
  filter(!(catch_eff == 1 & beta_catch_pct >= 0.003)) %>%
  mutate(., catch_eff_label = factor(paste0(as.numeric(catch_eff)*100, "% efficiency"),
                                     levels = paste0(as.numeric(unique(pull(., catch_eff)))*100, "% efficiency"))) %>% 
  {ggplot(data = filter(., paramset == "summary"),
          mapping = aes(x = beta_catch_pct, y = sip, 
                        ymin = lwr, ymax = upr,
                        color = county, 
                        fill = county, 
                        group = interaction(paramset,county))) + 
      geom_ribbon(alpha = 0.15, color = NA) +
      geom_line(alpha = 1, size = 1.5) + 
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
      scale_color_manual(values = fig1_colors[c(1, 2, 5, 3, 4)]) +
      scale_fill_manual(values = fig1_colors[c(1, 2, 5, 3, 4)]) +
      scale_shape_manual(guide = F, values = c(19, 17)) + #95)) + 
      facet_wrap(~catch_eff_label) + 
      theme(legend.position = c(0.9, 0.8),
            panel.spacing = unit(0.5, "lines"))}

fig3_curves

hist_vars <- list(k = 0.16,
                  R0 = 2.5, 
                  dt = 1/6, 
                  d = 7,
                  nsim = 10000,
                  seed = 1001,
                  N = 1000,
                  N_small = 10,
                  beta_catch = 0.001)

## THE FOLLOWING ARE LARGE RANDOM NUMBER GENERATIONS, but takes < 5 mins to run 
fig3_hist_data <- list(
  base_draws = with(hist_vars, {
    set.seed(seed)
    times = d/dt
    return(array(rgamma(times*N*nsim, shape = k*dt/d, scale = R0/k), 
                 dim = c(times, N, nsim)))
  }),
  trunc_draws = with(hist_vars, {
    set.seed(seed)
    times = d/dt
    return(array(rtgamma(times*N*nsim, shape = k*dt/d, scale = R0/k, 0, 1 - beta_catch), 
                 dim = c(times, N, nsim)))
  }))

fig3_hist_df = with(fig3_hist_data, data.frame(
  ind_step_base = base_draws[1,1,],
  ind_step_trunc = trunc_draws[1,1,],
  ind_inf_base = base_draws[,1,] %>% colSums(),
  ind_inf_trunc = trunc_draws[,1,] %>% colSums(),
  small_step_base = base_draws[1, 1:hist_vars$N_small, ] %>% colMeans,
  small_step_trunc = trunc_draws[1, 1:hist_vars$N_small, ] %>% colMeans,
  small_inf_base = base_draws[, 1:hist_vars$N_small,] %>% colSums %>% colMeans,
  small_inf_trunc = trunc_draws[, 1:hist_vars$N_small,] %>% colSums %>% colMeans,
  large_step_base = base_draws[1, , ] %>% colMeans,
  large_step_trunc = trunc_draws[1, , ] %>% colMeans,
  large_inf_base = base_draws %>% colSums %>% colMeans,
  large_inf_trunc = trunc_draws %>% colSums %>% colMeans
) %>% pivot_longer(ind_step_base:large_inf_trunc) %>%
  separate(name, into = c("size", "time", "dist")))

fig3_hist <- fig3_hist_df %>%
  filter((time == "step" & size == "ind") | (time == "inf" & size == "large")) %>% 
  mutate(time = factor(mapvalues(time, from = c("step", "inf"), 
                                 to = c("4-hour time step", "Infection duration")),
                       levels = c("4-hour time step", "Infection duration"))) %>% 
  mutate( dist = mapvalues(dist, from = c("base", "trunc"), to = c("Base", paste0("Truncated upper ", hist_vars$beta_catch*100, "%")))) %>%
  mutate(size = factor(mapvalues(size, 
                                 from = c("ind", "small", "large"), 
                                 to = c("individual", 
                                        paste0(hist_vars$N_small, " infected" ), 
                                        paste0(hist_vars$N, " infected"))),
                       levels = c("individual",
                                  paste0(hist_vars$N_small, " infected" ), 
                                  paste0(hist_vars$N, " infected")))) %>% 
  unite("time_size", c("time", "size"), sep = ", ") %>% 
  ggplot(aes(x = value, y = ..count.. / sum(..count..), group = dist, fill = dist, color = dist)) + 
  geom_histogram(color = NA, position = "identity", alpha = 0.8) + 
  # geom_density(adjust = 100, alpha = 0.5) +
  geom_vline(data = data.frame(xint=with(hist_vars, qgamma(1 - beta_catch, shape = k*dt/d, scale = R0/k)),
                               time = "4-hour time step",
                               size = "individual") %>%
               unite("time_size", c("time", "size"), sep = ", "),
             aes(xintercept = xint), linetype = "dashed") +
  facet_wrap( ~ time_size, scales = "free", nrow = 2) + 
  scale_fill_manual(name = "Distribution", values = c("#F21A00", "#3B9AB2")) + 
  scale_color_manual(name = "Distribution", values = c("#F21A00", "#3B9AB2")) + 
  scale_y_continuous(trans = "sqrt", labels = scales::percent_format(accuracy = 1)) + 
  scale_x_continuous(trans = "sqrt") + #, breaks = c(1, 5, 10, 20, 40, 100)) +
  xlab("Transmission rate") + 
  ylab("Proportion") + 
  theme(legend.position = c(0.785, 0.9)) 
fig3_hist

fig3 <- gridExtra::arrangeGrob(fig3_hist, fig3_curves, 
                               layout_matrix = matrix(c(1,2), 
                                                      byrow  = T, nrow = 1),
                               widths = c(1.3, 2.6))
# gridExtra::grid.arrange(fig3)
ggsave("figures/Manuscript2/Figure3.pdf", fig3,
       width = 20, height = 9)

figS3 <- fig3_hist_df %>%
  mutate(time = factor(mapvalues(time, from = c("step", "inf"), 
                                 to = c("4-hour time step", "Infection duration")),
                       levels = c("4-hour time step", "Infection duration")),
         dist = mapvalues(dist, from = c("base", "trunc"), to = c("Base", paste0("Truncated upper ", hist_vars$beta_catch*100, "%"))),
         size = factor(mapvalues(size, 
                                 from = c("ind", "small", "large"), 
                                 to = c("Individual", 
                                        paste0(hist_vars$N_small, " infected" ), 
                                        paste0(hist_vars$N, " infected"))),
                       levels = c("Individual",
                                  paste0(hist_vars$N_small, " infected" ),
                                  paste0(hist_vars$N, " infected")),
                       ordered = T)) %>% 
  ggplot(aes(x = value, y = ..count.. / sum(..count..), group = dist, fill = dist)) + 
  geom_histogram(color = NA, position = "identity", alpha = 0.8) + 
  geom_vline(data = data.frame(xint=with(hist_vars,
                                         qgamma(1 - beta_catch, shape = k*dt/d, scale = R0/k)),
                               time = "4-hour time step",
                               size = factor("Individual", 
                                             levels = c("Individual",
                                                        paste0(hist_vars$N_small, " infected" ),
                                                        paste0(hist_vars$N, " infected")))),
             aes(xintercept = xint), linetype = "dashed") +
  facet_grid(size ~ time, scales = "free") + 
  scale_fill_manual(name = "Distribution", values = c("#F21A00", "#3B9AB2")) + 
  scale_color_manual(name = "Distribution", values = c("#F21A00", "#3B9AB2")) + 
  scale_y_continuous(trans = "sqrt" , labels = scales::percent_format(accuracy = 1)) + 
  scale_x_continuous(trans = "sqrt", breaks = c(1, 5, 10, 20, 40, 100)) +
  xlab("Transmission rate") + 
  ylab("Proportion") + 
  theme(strip.text.x = element_text(size = 12, face = "bold")
        , strip.text.y = element_text(size = 12, face = "bold")
        , axis.text.x = element_text(size = 12))
figS3

ggsave("figures/Manuscript2/FigureS3.pdf",
       figS3,  
       width = 8, height = 8)
       

## Figure 4: Epidemic rebound when rare ----
# uses data simulated in figure4_expand 
# fig4_scale <- adply(gsub(".Rds|_run2", "", list.files("./output/figure4_scale_data", pattern = "Rds")) %>% unique,
#                   1, 
#                   function(j){
#                     print(j)
#                     files = list.files("./output/figure4_scale_data", full.names = T, pattern = j)
#                     print(files)
#                     test <- adply(1:length(files), 
#                                   1, 
#                                   function(k){
#                                     readRDS(files[k]) %>% mutate(.id = paste0(.id, "_", k))
#                                   })
#                     test_sum <- test %>% 
#                       select(-data, -X1, -county, -state) %>% 
#                       filter(!is.na(.id)) %>% 
#                       pivot_wider() %>% 
#                       group_by(.id, beta_catch, catch_eff, threshold_I) %>% 
#                       arrange(day) %>% 
#                       mutate(any_cross = max(thresh_crossed)) %>%
#                       filter(any_cross > 0) %>%
#                       mutate(first_cross = first(which(thresh_crossed > 0))  + first(date),
#                              days_post = as.numeric(date - first_cross),
#                              max_post = max(days_post)) %>% 
#                       ungroup %>% 
#                       unite("intervention", beta_catch:threshold_I, remove = F) %>% 
#                       select(-any_cross, -first_cross) %>%
#                       mutate(value = I) %>%
#                       group_by(days_post, 
#                                intervention, 
#                                beta_catch,
#                                catch_eff, 
#                                threshold_I) %>%
#                       mutate(extinct = (value == 0)) %>%
#                       summarise(n = n(),
#                                 prop_extinct = sum(extinct)/n(),
#                                 lwr = min(value),
#                                 lwr_005 = quantile(ifelse(extinct, NA, value), 0.005, na.rm = T),
#                                 lwr_025 = quantile(ifelse(extinct, NA, value), 0.025, na.rm = T),
#                                 lwr_05 = quantile(ifelse(extinct, NA, value), 0.05, na.rm = T),
#                                 lwr_01 = quantile(ifelse(extinct, NA, value), 0.01, na.rm = T),
#                                 lwr_10 = quantile(ifelse(extinct, NA, value), 0.1, na.rm = T),
#                                 lwr_99_all = quantile(value, 0.005),
#                                 lwr_95_all = quantile(value, 0.025),
#                                 lwr_90_all = quantile(value, 0.05),
#                                 upr = max(value),
#                                 upr_995 = quantile(ifelse(extinct, NA, value), 0.995, na.rm = T),
#                                 upr_99 = quantile(ifelse(extinct, NA, value), 0.99, na.rm = T),
#                                 upr_975 = quantile(ifelse(extinct, NA, value), 0.975, na.rm = T),
#                                 upr_95 = quantile(ifelse(extinct, NA, value), 0.95, na.rm = T),
#                                 upr_90 = quantile(ifelse(extinct, NA, value), 0.90, na.rm = T),
#                                 upr_99_all = quantile(value, 0.995),
#                                 upr_95_all = quantile(value, 0.975),
#                                 upr_90_all = quantile(value, 0.95),
#                                 med_all = median(value),
#                                 mean_all = mean(value),
#                                 mean = mean(ifelse(extinct, NA, value), na.rm = T),
#                                 median = median(ifelse(extinct, NA, value), na.rm = T),
#                                 .groups = "drop")
#                     return(test_sum)})
# saveRDS( fig4_scale, "output/figure4_scale_data/figure4_scale_summary_stats.rds")
fig4_scale <- readRDS("output/figure4_scale_data/figure4_scale_summary_stats.rds")

<<<<<<< HEAD
int_vars <- {
list(
 just_sip = list(
     counter.factual      = FALSE
   , int.movement         = c("post", "mid")
   , int.type             = "tail"
   , int.init             = "2020-07-01"
   , sim_end              = "2021-01-31"
   , sim_title            = "Just Shelter in Place"
   , thresh_inf.val       = 2
   , int.beta_catch_type  = "pct"
   , int.catch_eff        = 0.75
   , int.beta_catch       = 0.005
   , int.beta0_k          = 0.16
   , int.beta0_k_post     = 0.16
   , int.beta_catch_post  = 0
   , int.catch_eff_post   = 1
 ),
  minor_tail = list(
   counter.factual      = FALSE
 , int.movement         = c("post", "mid")
 , int.type             = "tail"
 , int.init             = "2020-07-01"
 , sim_end              = "2021-01-31"
 , sim_title            = "Minor tail chop"
 , thresh_inf.val       = 2
 , int.beta_catch_type  = "pct"
 , int.catch_eff        = 0.75
 , int.beta_catch       = 0.005
 , int.beta0_k          = 0.16
 , int.beta0_k_post     = 0.16
 , int.beta_catch_post  = 0.001
 , int.catch_eff_post   = 1
),
  major_tail = list(
   counter.factual      = FALSE
 , int.movement         = c("post", "mid")
 , int.type             = "tail"
 , int.init             = "2020-07-01"
 , sim_end              = "2021-01-31"
 , sim_title            = "Major tail chop"
 , thresh_inf.val       = 5
 , int.beta_catch_type  = "pct"
 , int.catch_eff        = 0.75
 , int.beta_catch       = 0.005
 , int.beta0_k          = 0.16
 , int.beta0_k_post     = 0.16
 , int.beta_catch_post  = 0.002
 , int.catch_eff_post   = 1
))}

source("ggplot_theme.R")
source("epidemic_rebound/gamma_rebound_params.R")
source("epidemic_rebound/gamma_rebound_pomp2.R")

plot_vars      <- c("cases", "deaths", "I", "thresh_crossed")
nsim           <- 1000
ci.stoch       <- 0.025
ci.epidemics   <- T
plot.median    <- F

fig4_data <- adply(1:length(int_vars), 1, 
      function(j) {
        adply(1:length(counties_list), 1,
      function(i) {
        with(c(counties_list[[i]], int_vars[[j]]), {
        source("epidemic_rebound/gamma_rebound.R", local = T)
        SEIR.sim.f.ci %<>% 
          full_join(county.data %>%
                      arrange(date) %>%
                      mutate(D = cumsum(ifelse(is.na(deaths), 0, deaths))*
                                 ifelse(is.na(deaths), NA, 1)) %>%
                      select(date, any_of(plot_vars)) %>%
                      pivot_longer(any_of(plot_vars), values_to = "data")) %>% 
          mutate(intervention = sim_title,
                 county       = focal.county,
                 state        = focal.state_abbr) 
         return(SEIR.sim.f.ci)
        
        }
        )
      }
        )
      }
        , .id = NULL)

# fig4_data <- SEIR.sim.f.ci

fig4_data$name         <- sapply(fig4_data$name, simpleCap)
fig4_data$intervention <- factor(fig4_data$intervention, levels = unique(fig4_data$intervention))

fig4_colors <- c("dodgerblue4", "#D67236", "firebrick3", "#0b775e")
# fig4_colors <- c("#D67236", "dodgerblue4", "red", "#0b775e", "magenta4")

check_date <- (fig4_data %>% filter(
#  .id           == "median"
   intervention == "Shelter in Place check"
 , date          > "2020-05-01"
 , name         == "I"
  ) %>% dplyr::filter(mid <= 10) %>% 
    filter(date  == min(date)))$date

## For just I plot
# fig4_data <- fig4_data %>% filter(date >= as.Date("2020-07-10")) %>% filter(name != "betat")

#fig4_data %>%
#  ggplot(aes(x = date, y = value
 #  , ymin = lwr, ymax = upr, fill  = intervention
#    , color = intervention
#    , group = .id)) +
#  geom_ribbon(data = (fig4_data %>% filter(date >= check_date))
#   , alpha = 0.50, colour = NA) +
#  geom_line(data = (fig4_data %>% filter(date >= check_date))
#   ) + 
#  geom_line(data = (fig4_data %>% filter(date >= check_date, .id != "median"))
#    , aes(group = interaction(.id, intervention))
#    , lwd = 0.10, alpha = 0.15
#   ) +
#  geom_line(data = (fig4_data %>% filter(date >= check_date, .id == "median"))
#    , aes(group = interaction(.id, intervention))
#   , lwd = 1.0) +
#  geom_ribbon(data = (fig4_data %>% filter(date <= check_date, intervention == "Continue Shelter in Place"))
#  , alpha = 0.50, colour = NA, fill = "black") +
#  geom_line(data = (fig4_data %>% filter(date <= check_date, intervention == "Continue Shelter in Place"))
#  , colour = "black") +
#  geom_line(data = (fig4_data %>% filter(date <= check_date, intervention == "Continue Shelter in Place", .id != "median"))
#    , lwd = 0.10, alpha = 0.15, colour = "black"
#   ) +
#  geom_line(data = (fig4_data %>% filter(date <= check_date, intervention == "Continue Shelter in Place", .id == "median"))
#    , colour = "black"
#    , lwd = 1.0) +
#  geom_vline(xintercept = check_date, linetype = "dashed", lwd = 0.5) + 
#  geom_point(aes(x = date, y = data), 
#             color = "black", size = 1) + 
#  scale_fill_manual(values = fig4_colors, name = "Intervention") +
#  scale_color_manual(values = fig4_colors, name = "Intervention") +
#  facet_grid(name ~ county, scales = "free_y", switch = "y")  +
#  ylab("") + 
#  scale_y_log10() + 
#  theme(
#    strip.background = element_blank()
#    , legend.position  = c(0.850, 0.885)
#    , legend.background = element_blank()
#    , strip.placement = "outside"
#    , strip.text = element_text(size = 16)
#    , axis.text.x = element_text(size = 14)) +
#  xlab("Date")

#droplevels(fig4_data) %>% 
#  filter(state == "CA", .id != "median") %>% 
#  filter(name == "cases") %>% 
#  filter(day > 290) %>% 
#  group_by(.id, intervention) %>%
#  summarize(new_cases = sum(value)) %>% 
#  filter(new_cases == 0) %>% 
#  group_by(intervention) %>%
#  summarize(prop_extinct = n() / nsim)

#droplevels(fig4_data) %>% filter(state == "CA", .id != "median") %>% filter(name == "I") %>% filter(day == max(day)) %>% 
#  group_by(.id, intervention) %>%
#  summarize(new_cases = sum(value)) %>% 
#  group_by(intervention) %>%
#  summarize(
#     lwr = quantile(value, 0.025)
#    , est = quantile(value, 0.500) 
#    , upr = quantile(value, 0.975))

check_date <- as.Date("2020-10-15")

fig4_data %>% filter(name == "I") %>% 
  filter(date > "2020-07-01") %>% {
    ggplot(aes(x = date, y = mid), data = filter(., date >= check_date)) +
  geom_ribbon(
    aes(
      ymin = lwr, ymax = upr
    , fill  = intervention
    , color = intervention
    )
  , alpha = 0.50, colour = NA) +
      geom_line(aes(color = intervention), data = filter(., date >= check_date), 
                lwd = 0.5) +
      geom_ribbon(data = filter(., date <= check_date, 
                                intervention == "Just Shelter in Place check")
                  , aes(ymin = lwr, ymax = upr)
   , alpha = 0.50, colour = NA, fill = "black") +
#  geom_line(data = filter(., date <= check_date, 
#                          intervention == "Just Shelter in Place check"), colour = "black") +
#      geom_vline(xintercept = check_date, linetype = "dashed", lwd = 0.5) +
      geom_point(aes(x = date, y = data), color = "black", size = 1) +
#  scale_y_continuous(trans = "log10") +
  scale_fill_manual(values = fig4_colors, name = "Intervention") +
  scale_color_manual(values = fig4_colors, name = "Intervention") +
 # facet_grid(name ~ county, scales = "free_y", switch = "y")  +
  ylab("Cases") + 
  theme(
    strip.background = element_blank()
    , legend.position  = c(0.250, 0.885)
    , legend.background = element_blank()
    , strip.placement = "outside"
    , strip.text = element_text(size = 16)
    , axis.text.x = element_text(size = 14)) +
  xlab("Date") +
  scale_y_continuous(trans = "log10")
  }
=======
fig4_colors = c("tomato4", "darkgoldenrod2", "slateblue")
fig4_traj_colors = c("#DC863B", "#0B775E", "#046C9A")
fig4_scale_extinct <- {fig4_scale %>% 
    mutate(file = mapvalues(X1, from = as.character(1:9),
                            to = gsub(".Rds|_run2", "", list.files("./output/figure4_scale_data", pattern = "Rds")) %>% unique)) %>%
    filter(days_post == 6*7) %>%
    filter(beta_catch == 0.001) %>%
    filter(grepl("percent", file)) %>% 
    mutate(threshold_I = as.factor(threshold_I)) %>%
    ggplot(aes(x = catch_eff, y = prop_extinct, 
               group = threshold_I,
               # group = interaction(intervention, X1), 
               # linetype = X1,
               color = threshold_I)) +
    geom_line() + 
    geom_point() + 
    scale_color_manual(values = fig4_colors, 
                       guide = F,
                       name = "Remaining infected when interventions relaxed") +
    xlab("Efficiency of truncation") + 
    ylab("Proportion of epidemic\nsimulations extinct") + 
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0.58, 0.99)) + 
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + 
    theme(legend.position = c(0.58, 0.73),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.background = element_rect(colour = "transparent", fill = "transparent"),
          plot.margin = margin(t = 25.5, r = 9.5, b = 5.5, l = 5.5))}

fig4_scale_uprCI <- {fig4_scale %>% 
    mutate(file = mapvalues(X1, from = as.character(1:9),
                            to = gsub(".Rds|_run2", "", list.files("./output/figure4_scale_data", pattern = "Rds")) %>% unique)) %>%
    filter(days_post == 6*7) %>%
    filter(beta_catch == 0.001) %>%
    filter(grepl("percent", file)) %>% 
    mutate(threshold_I = as.factor(threshold_I)) %>%
    ggplot(aes(x = catch_eff, 
               y = upr_99, 
               group = threshold_I,
               color = threshold_I)) +
    geom_line() + 
    geom_point() +
    scale_color_manual(values = fig4_colors, 
                       name = "Remaining infected when interventions relaxed",
                       guide = F) +
    xlab("Efficiency of truncation") +
    ylab("99th percentile of concurrent infected") + 
    ylim(0, 2350) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + 
    theme(plot.margin = margin(l = 20, t = 25.5, b = 5.5, r = 9.5),
          legend.position = c(0.6, 0.73),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.background = element_rect(colour = "transparent", fill = "transparent"))}

fig4_scale_traj <- {fig4_scale %>% 
    mutate(file = mapvalues(X1, from = as.character(1:9),
                            to = gsub(".Rds|_run2", "", list.files("./output/figure4_scale_data", pattern = "Rds")) %>% unique)) %>%
    filter(beta_catch == 0.001) %>%
    filter(threshold_I == 1) %>%
    filter(grepl("percent", file)) %>% 
    filter(!(catch_eff %in% c(0.2, 0.4, 0.8))) %>% 
    mutate(catch_eff_label = factor(scales::percent(catch_eff, accuracy = 1),
                                    levels = c("0%", "60%", "100%"))) %>%
    filter(days_post >= -2*7, days_post <= 6*7) %>% 
    ggplot(aes(x = days_post, 
               y = mean, 
               ymax = upr_99, 
               ymin = lwr_01,
               color = catch_eff_label, 
               fill = catch_eff_label,
               group = catch_eff_label)) + 
    geom_ribbon(color = NA, alpha = 0.5) +
    geom_line() + 
    scale_fill_manual(values = fig4_traj_colors, 
                      name = "Efficiency of truncation") +
    scale_color_manual(values = fig4_traj_colors, 
                       name = "Efficiency of truncation") +
    theme(legend.position = c(0.45, 0.8)) +
    xlab("Days since intervention relaxation") +
    ylab("Concurrent infections") +
    theme(plot.margin = margin(t = 25.5, r = 5.5, b = 5.5, l = 5.5))}



# fig4_noscale <- {adply(gsub(".Rds|_run2", "", list.files("./output/figure4_noscale_data", pattern = "Rds")) %>% unique,
#                     1,
#                     function(j){
#                       print(j)
#                       files = list.files("./output/figure4_noscale_data", full.names = T, pattern = j)
#                       print(files)
#                       test <- adply(1:length(files),
#                                     1,
#                                     function(k){
#                                       readRDS(files[k]) %>% mutate(.id = paste0(.id, "_", k))
#                                     })
#                       test_sum <- test %>%
#                         select(-data, -X1, -county, -state) %>%
#                         filter(!is.na(.id)) %>%
#                         pivot_wider() %>%
#                         group_by(.id, beta_catch, catch_eff, threshold_I) %>%
#                         arrange(day) %>%
#                         mutate(any_cross = max(thresh_crossed)) %>%
#                         filter(any_cross > 0) %>%
#                         mutate(first_cross = first(which(thresh_crossed > 0))  + first(date),
#                                days_post = as.numeric(date - first_cross),
#                                max_post = max(days_post)) %>%
#                         ungroup %>%
#                         unite("intervention", beta_catch:threshold_I, remove = F) %>%
#                         select(-any_cross, -first_cross) %>%
#                         mutate(value = I) %>%
#                         group_by(days_post,
#                                  intervention,
#                                  beta_catch,
#                                  catch_eff,
#                                  threshold_I) %>%
#                         mutate(extinct = (value == 0)) %>%
#                         summarise(n = n(),
#                                   prop_extinct = sum(extinct)/n(),
#                                   lwr = min(value),
#                                   lwr_005 = quantile(ifelse(extinct, NA, value), 0.005, na.rm = T),
#                                   lwr_025 = quantile(ifelse(extinct, NA, value), 0.025, na.rm = T),
#                                   lwr_05 = quantile(ifelse(extinct, NA, value), 0.05, na.rm = T),
#                                   lwr_01 = quantile(ifelse(extinct, NA, value), 0.01, na.rm = T),
#                                   lwr_10 = quantile(ifelse(extinct, NA, value), 0.1, na.rm = T),
#                                   lwr_99_all = quantile(value, 0.005),
#                                   lwr_95_all = quantile(value, 0.025),
#                                   lwr_90_all = quantile(value, 0.05),
#                                   upr = max(value),
#                                   upr_995 = quantile(ifelse(extinct, NA, value), 0.995, na.rm = T),
#                                   upr_99 = quantile(ifelse(extinct, NA, value), 0.99, na.rm = T),
#                                   upr_975 = quantile(ifelse(extinct, NA, value), 0.975, na.rm = T),
#                                   upr_95 = quantile(ifelse(extinct, NA, value), 0.95, na.rm = T),
#                                   upr_90 = quantile(ifelse(extinct, NA, value), 0.90, na.rm = T),
#                                   upr_99_all = quantile(value, 0.995),
#                                   upr_95_all = quantile(value, 0.975),
#                                   upr_90_all = quantile(value, 0.95),
#                                   med_all = median(value),
#                                   mean_all = mean(value),
#                                   mean = mean(ifelse(extinct, NA, value), na.rm = T),
#                                   median = median(ifelse(extinct, NA, value), na.rm = T),
#                                   .groups = "drop") %>%
#                         mutate(file = j)
#                       return(test_sum)})}
# saveRDS( fig4_noscale, "output/figure4_noscale_data/figure4_noscale_summary_stats.rds")


# fig4_noscale <- readRDS("output/figure4_noscale_data/figure4_noscale_summary_stats.rds")
fig4_noscale_extinct <- {fig4_noscale %>% 
    filter(days_post == 6*7) %>%
    filter(beta_catch == 0.001) %>%
    filter(grepl("percent", file)) %>% 
    mutate(threshold_I = as.factor(threshold_I)) %>%
    ggplot(aes(x = catch_eff, y = prop_extinct, 
               group = threshold_I,
               # group = interaction(intervention, X1), 
               # linetype = X1,
               color = threshold_I)) +
    geom_line() + 
    geom_point() + 
    scale_color_manual(values = fig4_colors, 
                       guide = F,
                       name = "Remaining infected when interventions relaxed") +
    xlab("Efficiency of truncation") + 
    ylab("Proportion of epidemic\nsimulations extinct") + 
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0.58, 0.99)) + 
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + 
    theme(legend.position = c(0.6, 0.73),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.background = element_rect(colour = "transparent", fill = "transparent"),
          plot.margin = margin(t = 25.5, r = 9.5, b = 5.5, l = 5.5))}

fig4_noscale_uprCI <- {fig4_noscale %>% 
    filter(days_post == 6*7) %>%
    filter(beta_catch == 0.001) %>%
    filter(grepl("percent", file)) %>% 
    mutate(threshold_I = as.factor(threshold_I)) %>%
    ggplot(aes(x = catch_eff, 
               y = upr_99, 
               group = threshold_I,
               color = threshold_I)) +
    geom_line() + 
    geom_point() +
    scale_color_manual(values = fig4_colors, 
                       name = "Remaining infected when\ninterventions relaxed") +
    xlab("Efficiency of truncation") +
    ylab("99th percentile of concurrent infected") + 
    ylim(0, 2350) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + 
    theme(plot.margin = margin(l = 20, t = 25.5, b = 5.5, r = 9.5),
          legend.position = c(0.6, 0.77),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.background = element_rect(colour = "transparent", fill = "transparent"))}

fig4_noscale_traj <- {fig4_noscale %>% 
    filter(beta_catch == 0.001) %>%
    filter(threshold_I == 1) %>%
    filter(grepl("percent", file)) %>% 
    filter(!(catch_eff %in% c(0.2, 0.4, 0.8))) %>% 
    mutate(catch_eff_label = factor(scales::percent(catch_eff, accuracy = 1),
                                    levels = c("0%", "60%", "100%"))) %>%
    filter(days_post >= -2*7, days_post <= 6*7) %>% 
    ggplot(aes(x = days_post, 
               y = mean, 
               ymax = upr_99, 
               ymin = lwr_01,
               color = catch_eff_label, 
               fill = catch_eff_label,
               group = catch_eff_label)) + 
    geom_ribbon(color = NA, alpha = 0.5) +
    geom_line() + 
    scale_fill_manual(values = fig4_traj_colors, 
                      name = "Efficiency of truncation") +
    scale_color_manual(values = fig4_traj_colors, 
                       name = "Efficiency of truncation") +
    theme(legend.position = c(0.45, 0.8)) +
    xlab("Days since intervention relaxation") +
    ylab("Concurrent infections") +
    theme(plot.margin = margin(t = 25.5, r = 5.5, b = 5.5, l = 5.5))}

fig4 <- gridExtra::arrangeGrob(
  gridExtra::arrangeGrob(fig4_noscale_traj, fig4_noscale_extinct, fig4_noscale_uprCI, 
                         layout_matrix = matrix(c(1, 2, 3), 
                                                byrow  = T, nrow = 1),
                         widths = c(1.5, 1.5, 1.5), 
                         clip = "off",
                         padding = unit(2, "line"),
                         top = grid::textGrob("Truncation with no shelter-in-place adjustment, resulting in variable mean", 
                                              x = 0.05, y = 0.7, hjust = 0, gp = gpar(fontsize=18, fontface = "bold"))),
  gridExtra::arrangeGrob(fig4_scale_traj, fig4_scale_extinct, fig4_scale_uprCI, 
                         layout_matrix = matrix(c(1, 2, 3), 
                                                byrow  = T, nrow = 1),
                         widths = c(1.5, 1.5, 1.5), 
                         clip = "off",
                         padding = unit(2, "line"),
                         top = grid::textGrob("Truncation with shelter-in-place adjustment, resulting in fixed mean", 
                                              x = 0.05, hjust = 0, y = 0.7, gp = gpar(fontsize=18, fontface = "bold"))))

# gridExtra::grid.arrange(fig4)

ggsave("figures/Manuscript2/Figure4_v2.pdf", 
       fig4, 
       width = 13, 
       height = 10)



# Figure S4 ----
figS4_extinct <- {fig4_scale %>% 
    mutate(file = mapvalues(X1, from = as.character(1:9),
                            to = gsub(".Rds|_run2", "", list.files("./output/figure4_scale_data", pattern = "Rds")) %>% unique)) %>%
    filter(days_post == 6*7) %>%
    filter(catch_eff == 1) %>%
    filter(grepl("betacatch", file)) %>% 
    mutate(threshold_I = as.factor(threshold_I)) %>%
    ggplot(aes(x = beta_catch, y = prop_extinct, 
               group = threshold_I,
               color = threshold_I)) +
    geom_line() + 
    geom_point() + 
    scale_color_manual(values = fig4_colors, 
                       name = "Remaining infected when interventions relaxed") +
    xlab("Efficiency of truncation") + 
    ylab("Proportion of epidemic\nsimulations extinct") + 
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
    scale_x_continuous(labels = scales::percent_format(accuracy = 0.001)) + 
    theme(legend.position = c(0.6, 0.72),
          legend.text = element_text(size = 12),
          legend.background = element_rect(colour = "transparent", fill = "transparent"),
          legend.title = element_text(size = 12),
          plot.margin = margin(l = 5.5, r = 11, t = 20.5, b = 5.5))}
# figS4_extinct

figS4_uprCI <- {fig4_scale %>% 
    mutate(file = mapvalues(X1, from = as.character(1:9),
                            to = gsub(".Rds|_run2", "", list.files("./output/figure4_scale_data", pattern = "Rds")) %>% unique)) %>%
    filter(days_post == 6*7) %>%
    filter(catch_eff == 1) %>%
    filter(grepl("betacatch", file)) %>% 
    mutate(threshold_I = as.factor(threshold_I)) %>%
    ggplot(aes(x = beta_catch, 
               y = upr_99, 
               group = threshold_I,
               color = threshold_I)) +
    geom_line() + 
    geom_point() + 
    scale_color_manual(values = fig4_colors, 
                       name = "Remaining infected when\ninterventions are relaxed",
                       guide = F) +
    xlab("Upper percentile of transmission rates averted") +
    ylab("99th percentile of concurrent infected") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 0.001)) + 
    theme(plot.margin = margin(l = 20.5, r = 11, t = 20.5, b = 5.5))}

figS4 <- gridExtra::arrangeGrob(figS4_extinct, figS4_uprCI, 
                                layout_matrix = matrix(c(1, 2), 
                                                       byrow  = F, nrow = 2),
                                widths = c(1.5))
ggsave("figures/Manuscript2/FigureS4.pdf", 
       figS4, 
       width = 6, 
       height = 7)

# Figure 4 version with no SIP scaling ----

fig5_data <- adply(list.files("./output/figure4_noscale_data", 
                              full.names = T, pattern = "Rds"),
                   1,
                   function(j){
                     test <- readRDS(j)
                     
                     test_sum <- test %>% select(-data, -X1, -county, -state) %>% 
                       filter(!is.na(.id)) %>% 
                       pivot_wider() %>% 
                       group_by(.id, beta_catch, catch_eff, threshold_I) %>% 
                       arrange(day) %>% 
                       mutate(any_cross = max(thresh_crossed)) %>%
                       filter(any_cross > 0) %>%
                       mutate(first_cross = first(which(thresh_crossed > 0))  + first(date),
                              days_post = as.numeric(date - first_cross),
                              max_post = max(days_post)) %>% 
                       ungroup %>% 
                       unite("intervention", beta_catch:threshold_I, remove = F) %>% 
                       select(-any_cross, -first_cross) %>%
                       mutate(value = I) %>%
                       group_by(days_post, 
                                intervention, 
                                beta_catch,
                                catch_eff, 
                                threshold_I) %>%
                       mutate(extinct = (value == 0)) %>%
                       summarise(n = n(),
                                 prop_extinct = sum(extinct)/n(),
                                 lwr = min(value),
                                 lwr_005 = quantile(ifelse(extinct, NA, value), 0.005, na.rm = T),
                                 lwr_025 = quantile(ifelse(extinct, NA, value), 0.025, na.rm = T),
                                 lwr_05 = quantile(ifelse(extinct, NA, value), 0.05, na.rm = T),
                                 lwr_01 = quantile(ifelse(extinct, NA, value), 0.01, na.rm = T),
                                 lwr_10 = quantile(ifelse(extinct, NA, value), 0.1, na.rm = T),
                                 lwr_99_all = quantile(value, 0.005),
                                 lwr_95_all = quantile(value, 0.025),
                                 lwr_90_all = quantile(value, 0.05),
                                 upr = max(value),
                                 upr_995 = quantile(ifelse(extinct, NA, value), 0.995, na.rm = T),
                                 upr_99 = quantile(ifelse(extinct, NA, value), 0.99, na.rm = T),
                                 upr_975 = quantile(ifelse(extinct, NA, value), 0.975, na.rm = T),
                                 upr_95 = quantile(ifelse(extinct, NA, value), 0.95, na.rm = T),
                                 upr_90 = quantile(ifelse(extinct, NA, value), 0.90, na.rm = T),
                                 upr_99_all = quantile(value, 0.995),
                                 upr_95_all = quantile(value, 0.975),
                                 upr_90_all = quantile(value, 0.95),
                                 med_all = median(value),
                                 mean_all = mean(value),
                                 mean = mean(ifelse(extinct, NA, value), na.rm = T),
                                 median = median(ifelse(extinct, NA, value), na.rm = T),
                                 .groups = "drop")
                     return(test_sum)
                   })

fig5_data %>% 
  filter(days_post == 6*7) %>%
  mutate(threshold_I = as.factor(threshold_I)) %>%
  ggplot(aes(x = catch_eff, 
             y = upr_99, 
             group = threshold_I,
             color = threshold_I)) +
  geom_line() + 
  geom_point()
  

# Figure S5: supplement figure showing effect of multiple gamma distributed infection periods on variance----
pdf("figures/Manuscript2/FigureS5.pdf",
    width = 12, height = 5)
par(mfrow = c(1,3))
hist_vars = list(R0 = 2.5,
                 k = 0.16, 
                 nstage = 7,
                 dt = 1/6,
                 d_avg = 7,
                 nsim = 1e6)

set.seed(1001)
with(hist_vars, replicate(nsim, {rgamma(1, k, scale = R0/k)}) %>% 
  {hist(., breaks = 100,main = paste0("A. entire infection Gamma(",k, ", ", R0, "/", k, ")\nmean = ", mean(.), "\nvariance = ", var(.)))})

with(hist_vars, replicate(nsim, {
  d <- rgeom(1, dt/d_avg)+1
  sum(rgamma(d, k*dt/d_avg, scale = R0/k))}) %>% 
  {hist(., breaks = 100,main = paste0("B. 1 geometrically distributed infection period\nmean = ", 
                                      mean(.), "\nvariance = ", var(.)))})

with(hist_vars, replicate(nsim, {
  d <- sum(rgeom(nstage, nstage*1/d_avg*dt)+1)
  sum(rgamma(d, k*dt/d_avg, scale = R0/k))}) %>% 
  {hist(., 
        breaks = 100,
        main = paste0("C. ", nstage, 
                      " geometrically distributed infection periods\nmean = ", 
                      mean(.), 
                      "\nvariance = ", 
                      var(.)))})
dev.off()
>>>>>>> b0f1377cef9c454b268494c5df28f4ed8a341f43
