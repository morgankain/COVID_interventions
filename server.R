needed_packages <- c(
  "shiny"
, "ggplot2"
, "deSolve"
, "plyr"
, "magrittr"
, "scales"
, "lubridate"
, "dplyr"
, "pomp"
, "reshape2"
, "shinythemes"
, "shinydashboard"
, "gridExtra"
)

# Commented out because package installation is not possible on shinyapps.io 
# If you are running locally, uncomment and run these lines to install missing packages

## Check if the packages are installed. *If they are not install them*, then load them
# if (length(setdiff(needed_packages, rownames(installed.packages()))) > 0) {
#   install.packages(setdiff(needed_packages, rownames(installed.packages())))
# }

lapply(needed_packages, require, character.only = TRUE)

# these packages need to be manually loaded for shinyapps.io 
library("scales")
library("lubridate")
library("gridExtra")

source("ggplot_theme.R")
source("covid-19_setup.R")

SCC    <- read.csv("./SantaClara_CumCases_20200317.csv", stringsAsFactors = F) %>% 
  select(-X, -X.1) %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y"))

params <- read.csv("covid_params.csv")
params <- params %>% 
  filter(!(Symbol %in% c("N", "E0", "Multiple"))) %>% 
  mutate(Value = as.numeric(as.character(Value)))

sir_step <- Csnippet("
                     // adjust betat for social distancing interventions
                     double betat;
                     if(intervention == 2 & thresh_crossed == 1){ // 2 is for threshhold intervention
                       betat =  beta0*thresh_int_level; 
                     }
                     else if(intervention == 1) betat = beta0*soc_dist_level; // 1 is for social distancing
                     else betat = beta0; // everything else is no intervention
                     
                     // adjust contact rates if isolation of symptomatic cases is in place
                     double iso_m = 1;
                     double iso_s = 1;
                     if(isolation == 1){
                        iso_m = iso_mild_level;
                        iso_s = iso_severe_level;
                     }
                    
                     // calculate transition numbers
                     double dSE = rbinom(S, 1-exp(-betat*(Ca*Ia/N + Cp*Ip/N + iso_m*Cm*Im/N + iso_s*Cs*Is/N)*dt)); 
                     double rateE[2];
                     double dE_all[2];
                     rateE[0] = alpha*gamma; // going to asymtomatic
                     rateE[1] = (1-alpha)*gamma; // going to presymptomatic
                     reulermultinom(2, E, rateE, dt, &dE_all);
                     double dEIa = dE_all[0];
                     double dEIp = dE_all[1];
                     double dIaR = rbinom(Ia, 1 - exp(-lambda_a*dt));
                     double rateIp[2];
                     double dIp_all[2];
                     rateIp[0] = mu*lambda_p; // going to minor symptomatic
                     rateIp[1] = (1-mu)*lambda_p; // going to sever symptomatic
                     reulermultinom(2, Ip, rateIp, dt, &dIp_all);
                     double dIpIm = dIp_all[0];
                     double dIpIs = dIp_all[1];
                     double dIsH = rbinom(Is, 1 - exp(-lambda_s*dt));
                     double dImR = rbinom(Im, 1 - exp(-lambda_m*dt));
                     double rateH[2];
                     double dH_all[2];
                     rateH[0] = delta*rho;
                     rateH[1] = (1-delta)*rho;
                     reulermultinom(2, H, rateH, dt, &dH_all);
                     double dHD = dH_all[0];
                     double dHR = dH_all[1];
                     
                     // update the compartments
                     S  -= dSE; // susceptible 
                     E  += dSE - dEIa - dEIp; // exposed
                     Ia += dEIa - dIaR; // infectious and asymptomatic
                     Ip += dEIp - dIpIs - dIpIm; // infectious and pre-symptomatic
                     Is += dIpIs - dIsH; // infectious and severe symptoms (that will be hospitalized)
                     Im += dIpIm - dImR; // infectious and minor symptoms
                     H  += dIsH - dHD - dHR; // hospitalized
                     R  += dHR + dImR + dIaR; // recovered
                     D  += dHD; // fatalities
                     sympt_new  +=  dIpIs + dIpIm;
                     H_new += dIsH;
                     if(intervention == 2 & H >= thresh_H_start) thresh_crossed = 1;
                     else if(intervention == 2 & thresh_crossed == 1 & H < thresh_H_end) thresh_crossed = 0;
                     ")

# define the initial set up, currently, every is susceptible except the exposed people
sir_init <- Csnippet("
                     S = N - E0;
                     E = E0;
                     Ia = 0;
                     Ip = 0;
                     Is = 0;
                     Im = 0;
                     H = 0;
                     R = 0;
                     D = 0;
                     sympt_new = 0;
                     H_new = 0;
                     thresh_crossed = 0;
                     ")

function(input, output, session) {
  
## Simulation parameters
covid_params        <- params$Value
names(covid_params) <- params$Symbol

## Update all of the sliders based on the user's first intervention parameters so dates wont break
observe({
  
updateSliderInput(
  session
, "int_start2"
, min   = input$int_start1 + input$int_length1
, max   = input$sim_len
, value = input$int_start1 + input$int_length1
, step  = 2
)
  
})

observe({
  
updateSliderInput(
  session
, "int_length2"
, min   = 0
, max   = input$sim_len - input$int_start2
, value = input$sim_len - input$int_start2
, step  = 2
)
  
})

observe({

updateSliderInput(
  session
, "iso_start"
, min   = 0
, max   = input$sim_len - 1
, value = round((input$sim_len - input$int_start1 + input$int_length1) / 2)
, step  = 2
)
  
})

observe({
  
updateSliderInput(
  session
, "iso_length"
, min   = 1
, max   = input$sim_len - input$iso_start
, value = input$sim_len - input$iso_start
, step  = 2
)
  
})

epi.dat <- eventReactive(input$do, {
    
covid_params        <- c(
    covid_params
  , N  = input$pop_size
  , E0 = ifelse(input$pop_size > 100000, 10, 1)
  )

## Update covid_params with user input, sort of an ugly way to do this but reduces
 ## need for input$ statements and reasonable way to be able to modify with new parameters 
covid_params["beta0"]    <- input$beta0
covid_params["Ca"]       <- input$Ca
covid_params["alpha"]    <- input$alpha
covid_params["lambda_a"] <- input$lambda_a
covid_params["lambda_s"] <- input$lambda_s
covid_params["lambda_m"] <- input$lambda_m
covid_params["lambda_p"] <- input$lambda_p
covid_params["rho"]      <- input$rho
covid_params["delta"]    <- input$delta
covid_params["mu"]       <- input$mu
  
sim_start  <- as.Date("2020-01-01")
sim_length <- input$sim_len
sim_end    <- sim_start + sim_length

dat        <- data.frame(
  day = 1:sim_length
, B   = rep(0, sim_length)
  )

intervention <- covariate_table(
  day              = 1:sim_length
  
, intervention     = c(
      # No intervention until intervention start time
    rep(0, input$int_start1)                   
      # Intervention style 1
  , rep(as.numeric(input$int_type1), input$int_length1)
      # Intermediate period if desired between interventions
  , rep(0, input$int_start2 - input$int_length1 - input$int_start1)
      # Intervention style 2
  , rep(as.numeric(input$int_type2), input$int_length2)
      # Post intervention close
  , rep(0, sim_length - input$int_start2 - input$int_length2)
  ) 
  
, isolation        = 
   { if (input$iso == 2) {
    c(
      # 0 is no isolation of cases
    rep(0, input$iso_start)      
      # 1 is isolation of symptomatic
  , rep(1, input$iso_length)  
      # Post isolation close  
  , rep(0, sim_length - input$iso_length - input$iso_start)
  ) 
   } else {
    c(
    rep(0, sim_length)      
    )     
   }
   }
, iso_severe_level = rep(input$iso_sm, sim_length)  # % of contats that severe cases maintain
, iso_mild_level   = rep(input$iso_mm, sim_length)  # % of contats that mild cases maintain
  
, soc_dist_level   = c(                       # intensity of the social distancing interventions
  rep(input$sd_m1, input$int_start1 + input$int_length1)
     ## slightly odd way to do this, but should work
, rep(input$sd_m2, sim_length - input$int_start1 - input$int_length1)
  )
  
, thresh_H_start   = c(                       # starting threshhold on total # in hospital
  rep(input$t_s1, input$int_start1 + input$int_length1)
     ## slightly odd way to do this, but should work
, rep(input$t_s2, sim_length - input$int_start1 - input$int_length1)
)  
, thresh_H_end     = c(                       # ending threshhold on total # in hospital
  rep(input$t_e1, input$int_start1 + input$int_length1)
     ## slightly odd way to do this, but should work
, rep(input$t_e2, sim_length - input$int_start1 - input$int_length1)   
)
     ## No reason to have a second parameter here, just use the same val that the user picks for social distancing
, thresh_int_level = c(                       # level of social distancing implemented with the threshhold intervention
  rep(input$sd_m1, input$int_start1 + input$int_length1)
     ## slightly odd way to do this, but should work
, rep(input$sd_m2, sim_length - input$int_start1 - input$int_length1)    
)
  
, order            = "constant"
, times            = "day"
  )

covid <- dat %>%
  pomp(
    time       = "day"
  , t0         = 1
  , covar      = intervention
  , rprocess   = euler(sir_step, delta.t = 1/6)
  , rinit      = sir_init
  , accumvars  = c("sympt_new", "H_new") # accumulate H until it gets measured, then zero it
  , paramnames = c(
      "beta0"
    , "Ca", "Cp", "Cs", "Cm"
    , "alpha"
    , "mu"
    , "delta"
    , "gamma"
    , "lambda_a", "lambda_s", "lambda_m", "lambda_p"
    , "rho"
    , "N"
    , "E0"
  )
  , statenames = c(
      "S" , "E" , "Ia"
    , "Ip", "Is", "Im"
    , "R" , "H" ,"D" 
    , "sympt_new", "H_new"
    , "thresh_crossed"
      )
  ) 

## Simulate with a set of parameters
epi.out <- do.call(
  pomp::simulate
, list(
  object       = covid
, params       = covid_params
, nsim         = input$num_sims
, format       = "d"
, include.data = F
, seed         = 1001)
  ) %>% {rbind(.,
         group_by(., day) %>%
           select(-.id) %>%
           summarise_all(median) %>%
                    mutate(.id = "median"))} 

epi.out <- epi.out %>%
  dplyr::mutate(
      date    = sim_start + day - 1
    , total_I = Is + Im + Ia + Ip
    )

epi.out <- epi.out %>% 
  group_by(.id) %>%
  mutate(
      new_deaths     = c(0, diff(D))
    , new_recoveries = c(0, diff(R))    
  )


d1.1 <- data.frame(
  x1  = c(sim_start + input$int_start1)
, x2  = c(sim_start + input$int_start1 + input$int_length1)
, y1  = c(0)
, y2  = c(max(
  epi.out[input$plotval]
  ))
, Intervention = "First Intervention"
)

d1.2 <- data.frame(
  x1  = c(sim_start + input$int_start2)
, x2  = c(sim_start + input$int_start2 + input$int_length2)
, y1  = c(0)
, y2  = c(max(
  epi.out[input$plotval]
  ))
, Intervention = "Second Intervention"
)

d1.3 <- data.frame(
  x1 = c(sim_start + input$iso_start)
, x2 = c(sim_start + input$iso_start + input$iso_length)
, y1 = c(max(epi.out[input$plotval]) - max(epi.out[input$plotval]) / ifelse(input$pscale == 1, 20, 5))
, y2 = c(max(epi.out[input$plotval]))
, Intervention = "Infected Isolation"
)

dtot     <- rbind(d1.1, d1.2, d1.3)

## return
  list(
    epi.out = epi.out
  , dtot    = dtot
    )
  
  })
  
  output$graph1 <- renderPlot({ 
    
epi.dat.s <- epi.dat()[["epi.out"]] %>%
  dplyr::group_by(.id) %>%
  dplyr::mutate(
    daily_cases = c(0, diff(total_I))
    )

## update Y axis
dtot <- epi.dat()[["dtot"]]
dtot <- mutate(dtot, y2 = max(epi.dat()[["epi.out"]][input$plotval]))

## Not used, for later...
adj_thresh <- epi.dat.s %>%
  dplyr::group_by(date) %>% 
  dplyr::summarize(thresh_c = mean(thresh_crossed)) %>% 
  mutate(adj_thresh = ifelse(thresh_c < 0.5, 0, 1))

col1   <- ifelse(input$int_type1 == "1", "seagreen4", "cadetblue4")
alpha1 <- ifelse(input$int_type1 == "1", 0.8, 0.6)
col2   <- ifelse(input$int_type2 == "1", "seagreen4", "cadetblue4")
alpha2 <- ifelse(input$int_type2 == "1", 0.8, 0.6)

col3   <- "brown4"
alpha3 <- ifelse(input$iso == 2, 0.5, 0)
    
coltot   <- c(col1, col2, col3)
alphatot <- c(alpha1, alpha2, alpha3)

ynames <- c(
      "total_I" = "Total Infected"
    , "H"       = "Hospitalized"
    , "Ia"      = "Asymptomatic Infections"
    , "Is"      = "Severe Infections"
    , "Im"             = "Minor Infections"
    , "new_recoveries" = "Recoveries"
    , "new_deaths"     = "Deaths")

gg1 <- epi.dat.s %>% ggplot() +
  
  geom_rect(data = dtot, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill = Intervention), color = NA, alpha = alphatot) +
  
  scale_fill_manual(values = coltot) +
  
  geom_line(
    data = (epi.dat.s %>% filter(.id != "median"))
  , aes(
    x     = date
  , y     = get(input$plotval)
  , group = .id
    ), alpha = 0.5, colour = "grey90") + 
  
  geom_line(data = (epi.dat.s %>% filter(.id == "median"))
  , aes(  
    x     = date
  , y     = get(input$plotval))
  , colour = "black", lwd = 1.5) +
  
  scale_x_date(labels = date_format("%Y-%b")) +
  scale_y_continuous(limits = c(0, max(epi.dat()[["epi.out"]][input$plotval]))) + 
  xlab("Date") + 
  ylab(ynames[input$plotval]) +
  guides(color = FALSE)

if (input$pscale == 2) {
  
  yscale1 <- c(seq(1, 10, by = 10) %o% 10^(1:ceiling(log10(max(epi.dat.s[input$plotval])))))
  gg1     <- gg1 + scale_y_continuous(trans = "pseudo_log", breaks = yscale1, labels = comma)
  
if (input$plotval == "H") {
  
  gg1     <- gg1 + geom_hline(yintercept = input$hosp_cap, lwd = 1, linetype = "dashed", colour = "firebrick4")
  
  tdf     <- data.frame(
      x = c(sim_start + 200)
    , y = c(input$hosp_cap) / 2
    , text = c("Hospital Beds Available")
)
  
  gg1     <- gg1 + geom_text(data = tdf, aes(x = x, y = y, label = text), colour = "firebrick4")
  
  }
  
}
  
#  gg1     <- gg1 + geom_line(
#    data = adj_thresh
#  , aes(x = date, y = adj_thresh$adj_thres * c(max(epi.dat.s[input$plotval]) - max(epi.dat.s[input$plotval]) / ifelse(input$pscale == 1, 20, 5))
#  ), lwd = 1, colour = "cadetblue4")

  gg1
  
    })
  
  output$graph2 <- renderPlot({ 
    
epi.dat.s <- epi.dat()[["epi.out"]] %>%
  dplyr::mutate(
      date    = sim_start + day - 1
    , total_I = Is + Im + Ia + Ip
    ) %>%
  dplyr::group_by(.id)

maxval <- max(c(max(epi.dat.s$R), max(epi.dat.s$D)))

col1   <- ifelse(input$int_type1 == "1", "seagreen4", "cadetblue4")
alpha1 <- ifelse(input$int_type1 == "1", 0.8, 0.6)
col2   <- ifelse(input$int_type2 == "1", "seagreen4", "cadetblue4")
alpha2 <- ifelse(input$int_type2 == "1", 0.8, 0.6)

col3 <- "brown4"
alpha3 <- ifelse(input$iso == 2, 0.5, 0)
    
coltot   <- c(col1, col2, col3)
alphatot <- c(alpha1, alpha2, alpha3)

dtot    <- epi.dat()[["dtot"]]
dtot$y2 <- maxval

epi.dat.s.r <- epi.dat.s %>% 
  dplyr::select(day, .id, R, D) %>% 
  mutate(date = sim_start + day - 1) %>% 
#  filter(.id == "median") %>%
  tidyr::pivot_longer(c(R, D), names_to = "Class", values_to = "Count")

gg2 <- epi.dat.s.r %>% ggplot() +
  geom_rect(data = dtot, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill = Intervention), color = NA, alpha = alphatot) +
  scale_fill_manual(values = coltot) +
  
  geom_line(
    data = (epi.dat.s.r %>% filter(.id == "median"))
  , aes(
    x        = date
  , y        = Count
  , colour   = Class
  , linetype = Class
  , group    = interaction(.id, Class)
    ), lwd = 1.5) + 
  
  geom_line(
    data = (epi.dat.s.r %>% filter(.id != "median"))
  , aes(
    x        = date
  , y        = Count
  , colour   = Class
  , linetype = Class
  , group    = interaction(.id, Class)
    ), alpha = 0.4, lwd = 0.5) + 
  
  scale_x_date(labels = date_format("%Y-%b")) +
  xlab("Date") + 
  ylab("Count") +
  scale_color_manual(values = c("black", "deeppink3")) +
  scale_linetype_manual(values = c("solid", "dashed"))

if (input$pscale == 2) {
  
  yscale1 <- c(seq(1, 10, by = 10) %o% 10^(1:ceiling(log10(maxval))))
  gg2     <- gg2 + scale_y_continuous(trans = "pseudo_log", breaks = yscale1, labels = comma)
  
}  

  gg2
    
    })
  
 output$datadown <- downloadHandler(

    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
		  paste(paste("epi.out", paste(strsplit(as.character(Sys.Date()), split = " +")[[1]], collapse = "_"), sep = "_"), "csv", sep = ".")
	  },

    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      write.csv(epi.dat()[["epi.out"]], file, row.names = FALSE)
    }
  )
  
}

