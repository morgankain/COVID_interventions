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

scc_wide <- readRDS("./COVID_interventions/santa_clara/variable_params_scc_wide_params.rds")  
  # mutate(R0_baseline = covid_R0(beta0est, 
  #                               c(Ca = Ca, Cp = 1, Cs = 1, Cm  = 1,
  #                                 alpha = alpha, gamma  = 1/4, 
  #                                 lambda_a = lambda_a, lambda_s = lambda_s,
  #                                 lambda_m = lambda_m, lambda_p = 1, 
  #                                 rho = 0.01428571, mu = 0.95600000), 
  #                               sd_strength = 1))
pairs(~E0 + sim_start + int_start1 + sd_m1 + sd_m2 + Ca + alpha + lambda_a +
        lambda_s + lambda_m + beta0est + log_lik , 
      data = scc_wide %>% 
        filter(log_lik > max(log_lik) - 2))

library(rlang)
library(ggpubr) 

scc_periods <- scc_wide %>% 
  mutate(lambda_m_inv = 1/lambda_m,
         lambda_a_inv = 1/lambda_a,
         lambda_s_inv = 1/lambda_s)

profiles <- mlply(data.frame(prof_var = c("Ca", "alpha", "lambda_a_inv", 
                                          "lambda_m_inv", "lambda_s_inv", 
                                          "sim_start", "E0", "beta0est",
                                          "sd_m2"),
                             prof_var_max = c("Ca[which.max(log_lik)]", 
                                              "alpha[which.max(log_lik)]", 
                                              "lambda_a_inv[which.max(log_lik)]", 
                                              "lambda_m_inv[which.max(log_lik)]", 
                                              "lambda_s_inv[which.max(log_lik)]",
                                              "sim_start[which.max(log_lik)]", 
                                              "E0[which.max(log_lik)]",
                                              "beta0est[which.max(log_lik)]",
                                              "sd_m2[which.max(log_lik)]"),
                             stringsAsFactors = F),
      function(prof_var, prof_var_max) {
        print(prof_var)
        ll_profile <- scc_periods %>% 
          filter(log_lik > -100) %>%
          mutate(sim_start = as.numeric(difftime(sim_start, as.Date("2019-12-31"))),
                 var_bin = if(prof_var == "sim_start" | prof_var == "E0") 
                   as.character(!!parse_quosure(prof_var)) else cut_number(!!parse_quosure(prof_var), n = 40))%>% 
          group_by(var_bin) %>% 
          summarise(max_loglik = max(log_lik),
                    var_max = !!parse_quosure(prof_var_max),
                    mean_var = mean(!!parse_quosure(prof_var)),
                    mean_loglik = mean(log_lik))
        
        return(list(plot = ggplot(ll_profile) +
                      ggtitle(prof_var) + 
                      xlab("variable values") + 
                      ylab("log likelihood") +
                      ylim(-59.5, -56.5) +
                      # geom_smooth(aes(x = var_max, y = max_loglik), col = "red") +
                      geom_line(aes(x = var_max, y = max_loglik), col = "black") +
                      geom_point(aes(x = var_max, y = max_loglik), col = "black") +
                      # geom_smooth(aes(x = mean_var, y = mean_loglik), col = "red") +
                      geom_line(aes(x = mean_var, y = mean_loglik), col = "blue") +
                      geom_point(aes(x = mean_var, y = mean_loglik), col = "blue"),
                    data = ll_profile))
      })

ggarrange(profiles[[1]][[1]],
          profiles[[2]][[1]],
          profiles[[3]][[1]],
          profiles[[4]][[1]],
          profiles[[5]][[1]], 
          profiles[[6]][[1]],
          profiles[[7]][[1]], 
          profiles[[8]][[1]], 
          profiles[[9]][[1]], 
          nrow = 3, ncol = 3)



fails <- mlply(data.frame(prof_var = c("Ca", "alpha", "lambda_a_inv", 
                                       "lambda_m_inv", "lambda_s_inv", 
                                       "sim_start", "E0", "beta0est",
                                       "sd_m2"),
                          prof_var_max = c("Ca[which.max(log_lik)]", 
                                           "alpha[which.max(log_lik)]", 
                                           "lambda_a_inv[which.max(log_lik)]", 
                                           "lambda_m_inv[which.max(log_lik)]", 
                                           "lambda_s_inv[which.max(log_lik)]",
                                           "sim_start[which.max(log_lik)]", 
                                           "E0[which.max(log_lik)]",
                                           "beta0est[which.max(log_lik)]",
                                           "sd_m2[which.max(log_lik)]"),
                          stringsAsFactors = F),
                  function(prof_var, prof_var_max) {
                    print(prof_var)
                    ll_profile <- scc_periods %>% 
                      mutate(sim_start = as.numeric(difftime(sim_start, as.Date("2019-12-31"))),
                             var_bin = if(prof_var == "sim_start" | prof_var == "E0") 
                               as.character(!!parse_quosure(prof_var)) else cut_number(!!parse_quosure(prof_var), n = 40) %>% as.character,
                             fail = ifelse(log_lik < -100, "fail", "success"))%>% 
                      group_by(var_bin, fail) %>% 
                      summarise(n = n()) %>% 
                      ungroup %>%
                      pivot_wider(names_from = fail, values_from =  n,
                                  values_fill = list(n = 0)) %>% 
                      mutate(success_frac = success/(success + fail),
                             total = success + fail) %>% 
                      separate(var_bin, into = c("bin_lower", "bin_upper"), sep = ",",
                               remove = F) %>%
                      mutate(bin_lower = gsub("\\(|\\[", "", bin_lower) %>% as.numeric,
                             bin_upper = gsub("]", "", bin_upper) %>% as.numeric,
                             bin_med = if(prof_var == "sim_start" | prof_var == "E0") as.numeric(var_bin) else (bin_lower+ bin_upper)/2)
                    
                    # head(ll_profile) %>% print
                    return(list(plot = ggplot(ll_profile) +
                             ggtitle(prof_var) + 
                             xlab("bin median") + 
                             ylab("pct successes") +
                             geom_line(aes(x = bin_med, y = success_frac), col = "black") +
                             geom_point(aes(x = bin_med, y = success_frac), col = "black"),
                             data = ll_profile))
                  })
ggarrange(fails[[1]][[1]],
          fails[[2]][[1]],
          fails[[3]][[1]],
          fails[[4]][[1]],
          fails[[5]][[1]],
          fails[[6]][[1]],
          fails[[7]][[1]],
          fails[[8]][[1]], 
          fails[[9]][[1]], 
          nrow = 3, ncol = 3)

# heatmap 

scc_periods %>% 
  # filter(log_lik > max(log_lik) - 100) %>%
  group_by(sim_start, E0) %>% 
  summarise(avg_log_lik = mean(log_lik),
            med_log_lik = median(log_lik),
            max_log_lik = max(log_lik)) %>% 
  ggplot(aes(x = sim_start, y = E0, fill = med_log_lik)) + 
  geom_tile()


par1 = "Ca"
par2 = "alpha"

scc_periods %>% 
  # filter(log_lik > max(log_lik) - 100) %>%
  mutate(par1_bin = cut(!!parse_quosure(par1), breaks = 10),
         par2_bin = cut(!!parse_quosure(par2), breaks = 10)) %>%
  group_by(par1_bin, par2_bin) %>% 
  summarise(avg_log_lik = mean(log_lik),
            med_log_lik = median(log_lik),
            max_log_lik = max(log_lik),
            num = n()) %>% 
  ungroup %>%
  mutate_at(vars(par1_bin, par2_bin), ~gsub("\\(|\\[|]", "", as.character(.))) %>% 
  separate(par1_bin, c("par1_lower", "par1_upper"), sep = ",") %>% 
  separate(par2_bin, c("par2_lower", "par2_upper"), sep = ",") %>% 
  mutate_at(vars(starts_with("par1"), starts_with("par2")), as.numeric) %>% 
  ggplot(aes(xmin = par1_lower, xmax = par1_upper,
             ymin = par2_lower, ymax = par2_upper, 
             fill = med_log_lik)) + 
  xlab(par1) + ylab(par2) +
  geom_rect()
