## Santa Clara County data 
# to read in these data, should set working directory to main COVID-interventions directory
source("./covid-19_setup.R")

SCC = read.csv("./SantaClara_CumCases_20200317.csv", stringsAsFactors = F) %>% 
  select(-X, -X.1) %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y"))


# simulate with a set of parameters
sim = covid %>%
  simulate(params=c(beta0 = 0.5, # without intervention beta for all categories
                    Ca = 1, Cp = 1, Cs = 1, Cm = 1, # category specific contact rates
                    alpha = 1/3, # fraction of cases asymptomatic
                    gamma = 1/5.2, # 1 over time in exposed class
                    lambda_a = 1/7, # 1/time for asympomatic to recover
                    lambda_s = 1/4, # 1/time for severely symptomatic to go to the hospitl 
                    lambda_m = 1/7, # 1/time for minorly sympomatic to recover
                    lambda_h = 1/10.7, # 1/time  to leaving hospital  
                    delta = 0.2, # fraction of hospitalized cases that are fatal
                    mu = 19/20, # fraction of cases that are minor
                    rho = 1/0.5, # 1/time in pre-symptomatic 
                    # N=59.02e6, # population size in Wuhan
                    N = 1937570, # Santa Clara County population
                    E0 = 1,  # initially exposed
                    intervention = 2, # 1 is for social distancing, 2 is for threshhold based, currently threshH is based on 
                    thresh_H_start = 15, # currently thresholding on total people in the hospital to start intervention
                    thresh_H_end = 2, # currently thresholding on total people in the hospital to end intervention
                    thresh_int_level = 0.1), # multiplier on beta when the thresshold causes the intervention to kick in
           nsim=50,format="d",include.data=F) %>%
# calulate the median of the simulations
  {rbind(.,
         group_by(., day) %>%
           select(-.id) %>%
           summarise_all(median) %>%
                    mutate(.id = "median"))} 

sim %>% 
  mutate(date = sim_start + day - 1) %>%
  # filter(.id == 3) %>%
ggplot() +
  geom_line(aes(x=date, 
                # y = Is + Im + Ia + Ip,
                y = H,
                group=.id, 
                color = .id == "median")) + 
  scale_x_date(labels = date_format("%Y-%b")) +
  guides(color=FALSE)+
  scale_color_manual(values=c("#D5D5D3", "#24281A")) +
  # geom_step(data = SCC, aes(x = Date, y = Cumulative_Deaths),
  #           color = "red") +
  # geom_vline(xintercept = as.Date("2020-03-17"), col = "blue") + # shelter in place date for SCC
  # geom_vline(xintercept = as.Date("2020-01-23"), col = "red") + # lockdown date for Wuhan
  # ylim(0, 100) + 
  theme_bw()

