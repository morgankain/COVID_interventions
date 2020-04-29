# have to source the COVID_pomp script first because this calls covid_R0
fitting        <- F
source("COVID_pomp.R")
# read all the outputs files ----
all_fits = alply(list.files("./output", full.names = T, pattern = ".rds"),
                  1, 
                  function(file_path){
                    print(file_path)
                    fit <- readRDS(file_path)
                    full_Reff <- matrix(nrow = dim(fit$param_array)[1], 
                                        ncol = dim(fit$param_array)[2], 
                                        data = 0)
                    full_R0 <- matrix(nrow = dim(fit$param_array)[1], 
                                        ncol = dim(fit$param_array)[2], 
                                        data = 0)
 # still using double forloop...                     
                    for (i in 1:dim(full_Reff)[1]) {
                      temp_dyn <- fit$dynamics_summary %>% filter(paramset == i)
                      for (j in 1:dim(full_Reff)[2]) {
                        full_Reff[i, j] <- with(fit$variable_params[i, ], covid_R0(
                          beta0est     = fit$param_array[i, j, 1]
                          , fixed_params = c(fit$fixed_params, unlist(fit$variable_params[i, ]))
                          , sd_strength  = fit$param_array[i, j, 2]
                          , prop_S = unlist(temp_dyn[temp_dyn$name == "S_now", 3]) /
                            (fit$fixed_params["N"] - temp_dyn[temp_dyn$name == "total_D", 3]))
                        )
                        full_R0[i, j] <- with(fit$variable_params[i, ], covid_R0(
                          beta0est     = fit$param_array[i, j, 1]
                          , fixed_params = c(fit$fixed_params, unlist(fit$variable_params[i, ]))
                          , sd_strength  = 1
                          , prop_S = 1)
                        )
                      }
                    }
                    fit$variable_params$county = gsub("./output/", "", 
                                                           strsplit(file_path, split = "_")[[1]][1])
                    fit$variable_params$fit_sip_eff = strsplit(file_path, split = "_")[[1]][2]
                    fit$variable_params$days_minus = as.numeric(strsplit(file_path, split = "_")[[1]][4])
                    fit$variable_params %<>% mutate(date = as.Date("2020-04-23") - days_minus)
                    fit$full_Reff = full_Reff
                    fit$full_R0 = full_R0
                    return(fit)
                  })

names(all_fits) = gsub(".rds", "",
  sapply(strsplit(list.files("./output", full.names = T, pattern = ".rds"), split = "/"), "[", 3))

# figure 1: Reff vs R0 -----
# # last_fit = all_fits[["Santa Clara_TRUE_FALSE_0_2020-04-25_final"]]
# # 
# rbind(all_fits[["Contra Costa_TRUE_FALSE_0_2020-04-28_final_D"]]$variable_params,
#       all_fits[["Contra Costa_FALSE_FALSE_0_2020-04-28_final_D"]]$variable_params) %>% 
#   filter(log_lik > -39.5) %>%
#   ggplot(aes(x = Reff, y = ..density.., fill = fit_sip_eff)) +
#   geom_histogram(position = "identity", alpha = 0.5) + 
#   geom_density(alpha = 0.5)

# really ugly way to set this up...
Rt_data <- rbind(data.frame(name = "Reff", 
                            value = as.vector(all_fits[["Santa Clara_TRUE_FALSE_0_2020-04-25_final"]]$full_Reff),
                            county = "Santa Clara"),
                 data.frame(name = "R0", 
                            value = as.vector(all_fits[["Santa Clara_TRUE_FALSE_0_2020-04-25_final"]]$full_R0),
                            county = "Santa Clara")
                 ) 
Rt_data %>%
  ggplot(aes(x = value, fill = name)) + 
  geom_histogram(aes(y = ..density..,),
                 position = "identity",
                 colour = "black", bins = 50,
                 lwd = 0.2, alpha = 0.4) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 1, lty = "dashed", alpha = 0.8) + 
  scale_fill_manual(name = "", 
                    values = c("#5BBCD6", "#FF0000", "#F98400", "#00A08A",    "#F2AD00", "black")) +
                    #values = c("firebrick3", "steelblue3", "green", "black")) +
  # facet_wrap(~county, ncol = 1) +
  theme_bw() + 
  theme(legend.position = "none",
        axis.text=element_text(size=15),
        axis.title=element_text(size=18)) +
  annotate("text", x = 1.4, y = 3.3, 
           hjust = 0.5,
           size = 6,
           parse = TRUE, 
           # label = as.character(expression("R"[e]^{}*" under shelter in place"))
           # label = expression(paste("R"[e]*"under\nshelter in place"))
           label = expression(atop(paste("R" [e]*" under"), "shelter-in-place"))
  ) +
  annotate("text", x = 3, y = 2.1, 
           hjust = 0.5,
           size = 6,
           parse = TRUE, label = as.character(expression("R"[0]))) +
  xlab(expression("R"[0]^{}*" estimate")) + 
  ylab("Density") 