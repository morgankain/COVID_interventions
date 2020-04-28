# R0 vs Reff estimats ----
fit <- readRDS("./output/Santa Clara_TRUE_FALSE_0_2020-04-25_final.rds")
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
as.vector(full_Reff) %>% median
as.vector(full_Reff) %>% quantile(probs = c(0.025, 0.975))


