
needed_packages <- c(
  "shiny"
, "shinythemes"
, "shinydashboard"
)

# Commented out because package installation is not possible on shinyapps.io 
# If you are running locally, uncomment and run these lines to install missing packages

## Check if the packages are installed. *If they are not install them*, then load them
# if (length(setdiff(needed_packages, rownames(installed.packages()))) > 0) {
#    install.packages(setdiff(needed_packages, rownames(installed.packages())))  
# }

lapply(needed_packages, require, character.only = TRUE)

# these packages need to be manually loaded for shinyapps.io 
library("shinythemes")
library("shinydashboard")

fluidPage(
  
  ## theme
   theme = shinytheme("lumen")
  
  ## Custom styling of any widgets
  , tags$style(
    " .checkbox {font-size: 18px}"
    )

  ## Title
  , headerPanel("Predicting the effects of COVID intervention strategies")
  
, h5("Starting parameter values have been gathered from previous work. 
However, because parameter values are likely to be variable across space and time, do not interpret our results as quantitative predictions for any specific location. 
Our goal is to show, qualitatively, the robust result that a second peak arises under a wide variety of intervention scenarios. 
Our ability to make quantiative predictions will improve as more data becomes available; parameter values will continue to be updated.")

, sidebarPanel(
  
    numericInput("pop_size"
      , label = h5("Starting population size")
      , min   = NA
      , max   = NA
      , value = 1937570
    )
  
        , radioButtons("int_type1"
                , h5("FIRST intervention strategy")
                , choices = list(
                    "Social Distancing" = 1
                  , "Threshold Based
('Lightswitch Method')"   = 2)
                  , selected  = 1
                 )
  
        , radioButtons("int_type2"
                , h5("SECOND intervention strategy")
                , choices = list(
                    "Social Distancing" = 1
                  , "Threshold Based
('Lightswitch Method')"   = 2)
                  , selected  = 2
                 )
      
  , radioButtons("iso"
    , h5("Isolation of symptomatic infected individuals?")
    , choices = list(
        "No"  = 1
      , "Yes" = 2)
       , selected  = 1
      )
        
   , conditionalPanel(condition = "input.iso == '2'"
   
  ## these values get updated according to other parameters, so they are just placeholders        
   , sliderInput("iso_start"
    , label = h5("Start date of symptomatic infected isolation")
    , min   = 0
    , max   = 200
    , value = 100
    , step  = 2
    ) 
           
   , sliderInput("iso_length"
    , label = h5("Length of symptomatic infected isolation")
    , min   = 0
    , max   = 100
    , value = 30
    , step  = 2
    ) 
           
      )
  
  , sliderInput("num_sims"
    , label = h5("Number of simulations")
    , min   = 5
    , max   = 200
    , value = 50
    , step  = 5
    )  
  
  , selectInput(
    "plotval", "Class to plot:"
    , c(
      "Total Infected"          = "total_I"
    , "Hospitalized"            = "H"
    , "Asymptomatic Infections" = "Ia"
    , "Severe Infections"       = "Is"
    , "Minor Infections"        = "Im"
    , "Recoveries"              = "new_recoveries"
    , "Dead"                    = "new_deaths")
    )
  
, tags$head(
    tags$style(HTML('#do{background-color:orange}'))
  )
  
, h5("Wait 3 seconds after the application loads before
clicking simulate to avoid an error")
, actionButton("do", "Simulate", width = '100%')

, radioButtons("pscale"
                , p("")
                , choices = list(
                    "Linear Scale" = 1
                  , "Log Scale"    = 2)
                  , selected  = 1
                 )
  
  )
    
    , mainPanel(
    
       column(6
         
  , h3("First Intervention", style = "color:maroon")

  , sliderInput("int_start1"
    , label = h5("Start date of intervention (days since first case)")
    , min   = 0
    , max   = 100
    , value = 50
    , step  = 2
    )  

  , sliderInput("int_length1"
    , label = h5("Length of intervention (number of days)")
    , min   = 0
    , max   = 200
    , value = 50
    , step  = 2
    )

  , sliderInput("sd_m1"
    , label = h5("Proportion of baseline contact rate (0 - 1)")
    , min   = 0
    , max   = 1
    , value = 0.25
    , step  = 0.05
    ) 
      
      , conditionalPanel(condition = "input.int_type1 == '2'"

  , sliderInput("t_s1"
    , label = h5("Threshold quantity: number of daily hospitalized cases before intervention STARTS")
    , min   = 0
    , max   = 100
    , value = 15
    , step  = 1
    )  
        
  , sliderInput("t_e1"
    , label = h5("Threshold quantity: number of daily hospitalized cases before intervention ENDS")
    , min   = 0
    , max   = 100
    , value = 2
    , step  = 1
    )  
        
      )  
        
        ) , column(6
          
    , h3("Second Intervention", style = "color:maroon")
       
    ## dates get updated according to choice in int_start1: see server.R
  , sliderInput("int_start2"
    , label = h5("Start date of intervention (days since first case)")
    , min   = 100
    , max   = 550
    , value = 100
    , step  = 2
    ) 
        
    ## dates get updated according to choice in int_start1: see server.R       
  , sliderInput("int_length2"
    , label = h5("Length of intervention (number of days)")
    , min   = 0
    , max   = 450
    , value = 450
    , step  = 2
    )
          
  , sliderInput("sd_m2"
    , label = h5("Proportion of baseline contact rate (0 - 1)")
    , min   = 0
    , max   = 1
    , value = 0.25
    , step  = 0.05
    )  
    
      , conditionalPanel(condition = "input.int_type2 == '2'"
        
  , sliderInput("t_s2"
    , label = h5("Threshold quantity: number of daily hospitalized cases before intervention STARTS")
    , min   = 0
    , max   = 100
    , value = 20
    , step  = 1
    )  
        
  , sliderInput("t_e2"
    , label = h5("Threshold quantity: number of daily hospitalized cases before intervention ENDS")
    , min   = 0
    , max   = 100
    , value = 8
    , step  = 1
    )  
        
      )  
        )
#          )
      
  ## Box to for the parameters
  , box(width = 12
      
  , tabsetPanel(type = "tabs"
    
    , tabPanel("Dynamics"
      , plotOutput("graph1")
      )
    
    , tabPanel("Cumulative"
      , plotOutput("graph2")
      )
    
    , tabPanel("Remaining model parameters"
      
    ,  column(6

  , sliderInput("sim_len"
    , label = h5("Length of simulation")
    , min   = 200
    , max   = 800
    , value = 550
    , step  = 5
    )  
      
  , sliderInput("hosp_cap"
    , label = h5("Available hospital beds")
    , min   = 200
    , max   = 10000
    , value = 4400
    , step  = 50
    )  
      
  , sliderInput("iso_mm"
    , label = h5("Proportion of baseline contact rate (0 - 1) for isolation intervention on minor infection")
    , min   = 0
    , max   = 1
    , value = 0.1
    , step  = 0.05
    )
      
  , sliderInput("iso_sm"
    , label = h5("Proportion of baseline contact rate (0 - 1) for isolation intervention on severe infection")
    , min   = 0
    , max   = 1
    , value = 0
    , step  = 0.05
    )
    
  , sliderInput("alpha"
    , label = h5("Fraction of cases asymptomatic")
    , min   = 0
    , max   = 1
    , value = 0.333
    , step  = 0.05
    ) 
      
  , sliderInput("mu"
    , label = h5("Fraction of cases that are minor")
    , min   = 0
    , max   = 1
    , value = 0.956
    , step  = 0.05
    )   
      
  , sliderInput("Ca"
    , label = h5("Transmission rate modifier for asymptomatic infection")
    , min   = 0
    , max   = 1
    , value = 0.667
    , step  = 0.05
    )    
      
  , sliderInput("delta"
    , label = h5("Fraction of hospitalized cases that are fatal")
    , min   = 0
    , max   = 1
    , value = 0.2
    , step  = 0.05
    )  
      
      )
      
  , column(6
        
  , sliderInput("beta0"
    , label = h5("Daily transmission Rate")
    , min   = 0
    , max   = 1
    , value = 0.5
    , step  = 0.05
    )    
    
  , sliderInput("gamma"
    , label = h5("1/time in exposed class")
    , min   = 0
    , max   = 1
    , value = 0.1923
    , step  = 0.05
    ) 
    
  , sliderInput("lambda_a"
    , label = h5("1/time for asymptomatics to recover")
    , min   = 0
    , max   = 1
    , value = 0.143
    , step  = 0.05
    )  
    
  , sliderInput("lambda_s"
    , label = h5("1/time for severely symptomatic to be hospitalized")
    , min   = 0
    , max   = 1
    , value = 0.174
    , step  = 0.05
    )  
    
  , sliderInput("lambda_m"
    , label = h5("1/time for minorly symptomatic to recover")
    , min   = 0
    , max   = 1
    , value = 0.143
    , step  = 0.05
    )  
    
  , sliderInput("lambda_p"
    , label = h5("1/time in pre-symptomatic")
    , min   = 0
    , max   = 2
    , value = 2
    , step  = 0.1
    )     
 
  , sliderInput("rho"
    , label = h5("1/time for leaving hospital")
    , min   = 0
    , max   = 1
    , value = 0.0689
    , step  = 0.05
    )  
  
  )
  
  ### Add: Cp, Cs, Cm?

      )
    
    , tabPanel("Download Current Run"
      , downloadButton("datadown", "Download")
    )
    
      )
    )
  )
)
  
        

        
        
  