
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
)

## Check if the packages are installed. *If they are not install them*, then load them
if (length(setdiff(needed_packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(needed_packages, rownames(installed.packages())))  
}

lapply(needed_packages, require, character.only = TRUE)

source("ggplot_theme.R")

SCC <- read.csv("./SantaClara_CumCases_20200317.csv", stringsAsFactors = F) %>% 
  select(-X, -X.1) %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y"))

sir_step <- Csnippet("
                     double betat;
                     if(intervention == 2 & thresh_crossed == 1){ // 2 is for threshhold intervention
                       betat =  beta0*thresh_int_level; 
                     }
                     else if(intervention == 1) betat = beta0*contact; // 1 is for social distancing
                     else betat = beta0; // everything else is no intervention
                     double dSE = rbinom(S, 1-exp(-betat*(Ca*Ia/N + Cp*Ip/N + Cm*Im/N + Cs*Is/N)*dt)); 
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
                     rateIp[0] = mu*rho; // going to minor symptomatic
                     rateIp[1] = (1-mu)*rho; // going to sever symptomatic
                     reulermultinom(2, Ip, rateIp, dt, &dIp_all);
                     double dIpIm = dIp_all[0];
                     double dIpIs = dIp_all[1];
                     double dIsH = rbinom(Is, 1 - exp(-lambda_s*dt));
                     double dImR = rbinom(Im, 1 - exp(-lambda_m*dt));
                     double rateH[2];
                     double dH_all[2];
                     rateH[0] = delta*lambda_h;
                     rateH[1] = (1-delta)*lambda_h;
                     reulermultinom(2, H, rateH, dt, &dH_all);
                     double dHD = dH_all[0];
                     double dHR = dH_all[1];
                     
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
                     S = N-E0;
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

function(input, output) {
  
  epi.dat <- reactive({
  
sim_length <- as.Date("2020-12-01") - as.Date("2019-12-01")
dat        <- data.frame(
  day = 0:sim_length
, B   = rep(0, sim_length + 1)
  )

## Intervention start and end date (only used if intervention == 1)

int_start  <- (sim_start + input$int_start) - sim_start
int_length <- (sim_start + input$int_len)   - sim_start
int_level  <- input$int_size

## use the intervention info to construct a covariate table for use in the pomp object
contact_rate <- covariate_table(
  day     = 0:(sim_length)
, contact = c(rep(1, int_start)
    , rep(int_level, int_length)
    , rep(1, sim_length - int_start - int_length + 1)
  )
, order   = "constant"
, times   = "day"
  )  

covid <- dat %>%
  pomp(
    time       = "day"
  , t0         = 0
  , covar      = contact_rate
  , rprocess   = euler(sir_step, delta.t = 1/6)
  , rinit      = sir_init
  , accumvars  = c("sympt_new", "H_new"), # accumulate H until it gets measured, then zero it
    paramnames = c(
      "beta0"
    , "Ca", "Cp", "Cs", "Cm"
    , "alpha"
    , "gamma"
    , "lambda_a", "lambda_s","lambda_m", "lambda_h"
    , "delta"
    , "mu"
    , "rho"
    , "N"
    , "E0"
    , "intervention"
    , "thresh_H_start"
    , "thresh_H_end"
    , "thresh_int_level"
  )
  , statenames = c(
   "S","E","Ia", "Ip","Is","Im"
  ,"R", "H","D", "sympt_new", "H_new"
  , "thresh_crossed"
    )
  ) 

pop_size <- 1937570

## Simulation parameters
covid_params <- c(
  beta0            = 0.5
, Ca               = 1
, Cp               = 1
, Cs               = 1
, Cm               = 1
, alpha            = 1/3
, gamma            = 1/5.2
, lambda_a         = 1/7
, lambda_s         = 1/4
, lambda_m         = 1/7
, lambda_h         = 1/10.7
, delta            = 0.2
, mu               = 19/20
, rho              = 1/0.5
, N                = pop_size # (Santa Clara County) 59.02e6 (Wuhan)
, E0               = 1
, intervention     = as.numeric(input$int_type)
, thresh_H_start   = input$int_start_t
, thresh_H_end     = input$int_end_t
, thresh_int_level = input$int_size_t
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

## return
  list(
    epi.out = epi.out
    )
  
  })

  output$graph1 <- renderPlot({ 
    
  epi.dat()[["epi.out"]] %>% 
   mutate(date = sim_start + day - 1) %>%
   ggplot() + geom_line(
     aes(x = date
       , y = H # y = Is + Im + Ia + Ip,
       , group = .id
       , color = .id == "median")
     ) + 
  scale_x_date(labels = date_format("%Y-%b")) +
  guides(color = FALSE) +
  scale_color_manual(values = c("#D5D5D3", "#24281A")) 
    
    })
  
  output$graph2 <- renderPlot({ 
    ggplot()    
    })
  
  output$graph3 <- renderPlot({ 
    
  epi.out.s <- epi.dat()[["epi.out"]] %>% 
  filter(.id != "median") %>%
  group_by(.id) %>% 
  dplyr::summarize(
    peak_val  = max(Ia + Ip + Is)
  , peak_time = which((Ia + Ip + Is) == max(Ia + Ip + Is))[1]
  ) %>% tidyr::pivot_longer(c(peak_val, peak_time), names_to = "Est")
    
  epi.out.s %>% ggplot(aes(Est, value)) + geom_boxplot()
     
    })
  
  output$graph4 <- renderPlot({ 
    ggplot()   
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

