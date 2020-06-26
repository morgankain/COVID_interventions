# functions and variables needed for all figures
source("needed_packages.R")
source("ggplot_theme.R")

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
 , loglik.num      = 3
 , ci.stoc         = 0.025
 , ci.epidemic     = T
 , ci.epidemic_cut = 500
 , nsim            = 500
 , plot_vars       = c("cases", "deaths")
 , plot.median     = F
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
            rbind(., (detect %>% 
                        transmute(date = date, paramset = paramset, name = "detect",
                                  lwr = NA, mid = detect, upr = NA, data = NA) %>%
                select(-detect) %>% dplyr::select(date, paramset, everything()))) %>%
            mutate(intervention = sim_title,
                 county = variable_params[1, "county"] %>% as.character,
                 state  = variable_params[1, "state"] %>% as.character)
        return(SEIR.sim.f.ci)
})}, .id = NULL)

# plot the fits and data
fig1_data$name <- sapply(fig1_data$name, simpleCap)

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
, sim_title           = "Superspreading_1%Tail_75%Eff"
, int.beta_catch      = 0.01
, int.beta_catch_type = "pct"
, int.catch_eff       = 0.75
),
  tail_chop2 = list(
  counter.factual     = FALSE
, int.movement        = c("post", "mid")
, int.type            = "tail"
, int.init            = "2020-07-01"
, sim_end             = "2020-08-31"
, sim_title           = "Superspreading_2%Tail_50%Eff"
, int.beta_catch      = 0.02
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
, int.beta_catch      = 0.01
, int.beta_catch_type = "pct"
, int.catch_eff       = 1.0
)
)
}

sim_vars <- list(
 loglik.max       = TRUE
, plot_vars       = c("cases", "deaths")
, ci.stoc         = 0.10
, ci.epidemic     = T
, ci.epidemic_cut = 500
, nsim            = 100
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

fig2_colors <- c("#D67236", "dodgerblue4", "#0b775e", "magenta4")
fig2_data$name <- sapply(fig2_data$name, simpleCap)

fig2_data %>% 
  filter(date >= as.Date("2020-02-10")) %>%
  filter(intervention %in% c(
    "Continue SIP"
  , "Isolate"
  , "Lift"
  , "Superspreading_1%Tail_75%Eff"
# , "Superspreading_2%Tail_50%Eff"
# , "Superspreading_1%Tail_100%Eff"
                             )) %>% 
   mutate(intervention = mapvalues(intervention,
                                   from = c("Continue SIP",
                                            "Isolate",
                                            "Lift",
                                            "Superspreading_1%Tail_75%Eff"),
                                    to  = c(
                                       "Continue shelter in place"
                                     , "Infected isolation"
                                     , "Lift all interventions"
                                     , "Superspreading Averted(Top 1% of\ndistribution removed with 75% effectiveness)"
                                   ))) %>%
   mutate(county = mapvalues(county,
                                   from = c("King",
                                            "Los Angeles"
                                     ),
                                    to  = c(
                                       "Seattle, WA"
                                     , "Los Angeles, CA"
                                   ))) %>% 
   mutate(intervention = factor(intervention, levels = c(
     "Lift all interventions"
   , "Continue shelter in place"
   , "Infected isolation"
   , "Superspreading Averted(Top 1% of\ndistribution removed with 75% effectiveness)"
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
  scale_y_continuous(trans = "log10", labels = trans_format('log10', math_format(10^.x))) +
  scale_fill_manual(values = fig2_colors[c(2,1,3,4)]
    , name = "Intervention") +
  scale_color_manual(values = fig2_colors[c(2,1,3,4)]
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
loglik.num      <- 3

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

fig3_curves <- fig3_data %>% 
  group_by(county, catch_eff, beta_catch_pct) %>% 
  {rbind(.,
         filter(., paramset != "data") %>% 
           summarise(., 
                     lwr = min(sip),
                     upr = max(sip),
                     sip = mean(sip)) %>% 
           mutate(paramset = "summary"))} %>% 
  ungroup() %>%
  filter(!(catch_eff == 1 & beta_catch_pct > 0.005)) %>%
  mutate(., catch_eff_label = factor(paste0(as.numeric(catch_eff)*100, "% efficiency"),
                                     levels = paste0(as.numeric(unique(pull(., catch_eff)))*100, "% efficiency"))) %>% 
  {ggplot(data = filter(., paramset == "summary"),
          mapping = aes(x = beta_catch_pct, y = sip, 
                        ymin = lwr, ymax = upr,
                        color = county, 
                        fill = county, 
                        group = interaction(paramset,county))) + 
      geom_ribbon(alpha = 0.25, color = NA) +
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
system.time({base_draws = with(hist_vars, {
  set.seed(seed)
  times = d/dt
  return(array(rgamma(times*N*nsim, shape = k*dt/d, scale = R0/k), 
        dim = c(times, N, nsim)))
})})

system.time({trunc_draws = with(hist_vars, {
  set.seed(seed)
  times = d/dt
  return(array(rtgamma(times*N*nsim, shape = k*dt/d, scale = R0/k, 0, 1 - beta_catch), 
               dim = c(times, N, nsim)))
})})

fig3_hist_data = data.frame(
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
  separate(name, into = c("size", "time", "dist"))

fig3_hist <- fig3_hist_data %>%
  filter((time == "step" & size == "ind") | (time == "inf" & size == "large")) %>% 
  mutate(time = factor(mapvalues(time, from = c("step", "inf"), 
                                 to = c("4-hour time step", "Infection duration")),
                       levels = c("4-hour time step", "Infection duration"))) %>% 
  mutate( dist = mapvalues(dist, from = c("base", "trunc"), to = c("Base", paste0("Truncated upper ", hist_vars$beta_catch*100, "%")))) %>%
  mutate(size = factor(mapvalues(size, 
                                 from = c("ind", "small", "large"), 
                                 to = c("individual", paste0(hist_vars$N_small, " infected" ), 
                                        paste0(hist_vars$N, " infected"))),
                       levels = c("individual", paste0(hist_vars$N_small, " infected" ), 
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

# Only a few small comments, bottom left is number infected? 
# If so the title is a bit confusing. 
# Scale the y axis in the histograms? 
# Other than that just tiny axis cosmetics and stuff 
# but I think this is about all we need to show in the
# main text and the big 6 panel figure in the supp will 
# fill in the rest


fig3 <- gridExtra::arrangeGrob(fig3_hist, fig3_curves, 
                               layout_matrix = matrix(c(1,2), 
                                                      byrow  = T, nrow = 1),
                               widths = c(1.3, 2.6))
# gridExtra::grid.arrange(fig3)
ggsave("figures/Manuscript2/Figure3.pdf", fig3,
       width = 20, height = 9)

figS3 <- fig3_hist_data %>%
  mutate(time = factor(mapvalues(time, from = c("step", "inf"), to = c("4-hour time step", "infection")),
                                 levels = c("4-hour time step", "infection"))) %>% 
  mutate(size = factor(mapvalues(size, 
                                 from = c("ind", "small", "large"), 
                                 to = c("individual", paste0(hist_vars$N_small, " infected" ), 
                                        paste0(hist_vars$N, " infected"))),
                       levels = c("individual", paste0("small population, N = ", hist_vars$N_small), 
                                  paste0("large population, N = ", hist_vars$N)))) %>%
  ggplot(aes(x = value, y = ..density.., group = dist, fill = dist, color = dist)) + 
  geom_histogram(color = NA, position = "identity", alpha = 0.5) + 
  # geom_density(adjust = 100, alpha = 0.5) + 
  geom_vline(data = data.frame(xint=with(hist_vars, qgamma(1 - beta_catch, shape = k*dt/d, scale = R0/k)),
                               time = "4-hour time step", 
                               size = "individual"), 
             aes(xintercept = xint), linetype = "dashed") +
  facet_grid(size ~ time, scales = "free") + 
  scale_fill_manual(values = c("red", "blue")) + 
  scale_color_manual(values = c("red", "blue")) + 
  scale_y_continuous(trans = "sqrt") + 
  scale_x_continuous(trans = "sqrt", breaks = c(1, 5, 10, 20, 40, 100)) +
  xlab("infection rate")

figS3

## Figure 4: Epidemic rebound when rare ----
####  

## Have to run figure 3 in advance, which is annoying and likely should be fixed

counties_list <- {
  list(
  FC = list(
    focal.county     = "Fulton"
  , focal.state      = "Georgia"
  , focal.state_abbr = "GA"
  , rds.name         = "output/Fulton_0_2020-06-25_cont_round2.rds"
  , con_theta        = F
  )
)
}

int_vars <- {
list(
 just_sip_check = list(
     counter.factual      = FALSE
     , int.movement         = c("post", "mid")
     , int.type             = "tail"
     , int.init             = "2020-07-01"
     , sim_end              = "2021-01-31"
     , sim_title            = "Just Shelter in Place check"
     , thresh_inf.val       = 5
     , int.beta_catch_type  = "pct"
     , int.catch_eff        = 0.75
     , int.beta_catch       = 0.005
     , int.beta0_k          = 0.16
     , int.beta0_k_post     = 0.16
     , int.beta_catch_post  = 0
     , int.catch_eff_post   = 1
   ),
 just_sip = list(
   counter.factual      = FALSE
   , int.movement         = c("post", "mid")
   , int.type             = "tail"
   , int.init             = "2020-07-01"
   , sim_end              = "2021-01-31"
   , sim_title            = "Just Shelter in Place"
   , thresh_inf.val       = 5
   , int.beta_catch_type  = "pct"
   , int.catch_eff        = 0.75
   , int.beta_catch       = 0.005
   , int.beta0_k          = 0.16
   , int.beta0_k_post     = 0.16
   , int.beta_catch_post  = 0
   , int.catch_eff_post   = 0
 ),
  minor_tail = list(
   counter.factual      = FALSE
 , int.movement         = c("post", "mid")
 , int.type             = "tail"
 , int.init             = "2020-07-01"
 , sim_end              = "2021-01-31"
 , sim_title            = "Minor tail chop"
 , thresh_inf.val       = 5
 , int.beta_catch_type  = "pct"
 , int.catch_eff        = 0.75
 , int.beta_catch       = 0.005
 , int.beta0_k          = 0.16
 , int.beta0_k_post     = 0.16
 , int.beta_catch_post  = 0.0001
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
 , int.beta_catch_post  = 0.001
 , int.catch_eff_post   = 1
))}

source("ggplot_theme.R")
source("epidemic_rebound/gamma_rebound_params.R")
source("epidemic_rebound/gamma_rebound_pomp2.R")

nsim           <- 200
plot_vars      <- c("cases", "deaths", "I")
ci.stoch       <- 0.0005
ci.epidemics   <- F
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

# fig4_colors <- c("#D67236", "dodgerblue4", "#0b775e", "magenta4")
fig4_colors <- c("#D67236", "dodgerblue4", "red", "#0b775e", "magenta4")

check_date <- (fig4_data %>% filter(
#  .id           == "median"
   intervention == "Just Shelter in Place"
 , date          > "2020-05-01"
  ) %>% dplyr::filter(mid == 0) %>% 
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

fig4_data %>% 
  filter(name == "Cases") %>% 
  filter(date > "2020-07-01") %>%
  {ggplot(aes(x = date, y = mid), data = filter(., date >= check_date)) +
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
                                intervention == "Just Shelter in Place")
                  , aes(ymin = lwr, ymax = upr)
   , alpha = 0.50, colour = NA, fill = "black") +
  geom_line(data = filter(., date <= check_date, 
                          intervention == "Just Shelter in Place"), colour = "black") +
      geom_vline(xintercept = check_date, linetype = "dashed", lwd = 0.5) +
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
  scale_y_continuous(trans = "log10")}
