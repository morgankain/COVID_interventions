if (fit.with == "D_C" | fit.with == "D") {
## Scrape death data from the NYT github repo to stay up to date or load a previously saved dataset
#deaths <- fread("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv")
deaths  <- read.csv("us-counties.txt")
deaths  <- deaths %>% mutate(date = as.Date(date)) %>% dplyr::filter(county == focal.county)
deaths  <- deaths %>% dplyr::filter(date < max(date) - fit.minus)

## One special case
if (focal.county == "Fulton") {
deaths  <- read.csv("us-counties.txt")
deaths  <- deaths %>% mutate(date = as.Date(date)) %>% 
  dplyr::filter((county == focal.county | county == "DeKalb") & state == "Georgia") %>%
  group_by(date) %>% summarize(cases = sum(cases), deaths = sum(deaths))
deaths <- deaths %>% mutate(county = focal.county, state = "Georgia")
}

} else if (fit.with == "H") {
## Not supported right now for SCC, ony currently have data for Contra Costa County.
 ## To use H, supply your own data and change the path
hospit     <- read.csv("contra_costa/ccc_data.csv")
hospit     <- hospit %>% 
  mutate(date = as.Date(REPORT_DATE)) %>% 
  filter(CURRENT_HOSPITALIZED != "NULL") %>% 
  mutate(ch = as.numeric(as.character(CURRENT_HOSPITALIZED))) %>% 
  dplyr::select(date, ch)
hospit    <- hospit %>% dplyr::filter(date < max(date) - fit.minus)  
}
  
## Load the previously saved fits
 ## If COVID_fit_cont.R was just run, use parameers already stored in the global env 
if (use.rds) {
prev.fit         <- readRDS(rds.name)
variable_params  <- prev.fit[["variable_params"]]
fixed_params     <- prev.fit[["fixed_params"]]
}
  
## drop the rows that have 0s for likelihood (in case exited prematurely) 
 ## and keep only the best fits as defined by loglik
variable_params <- variable_params %>% 
  filter(log_lik != 0) %>%
  filter(log_lik == max(log_lik)) ## only used for manuscript dyanmic trajectory plots
