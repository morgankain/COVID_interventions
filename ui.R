
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
  
    , h5("Model development: Marissa Childs, Devin Kirk, Morgan Kain, Mallory Harris, Nicole Nova")
    , h5("Model coding: Marissa Childs")
    , h5("Shiny app:", tags$a("Morgan Kain", href = "http://www.morgankain.weebly.com"))
    , h5("Parameter estimate search: Mallory Harris, Lisa Couper")
  
  # , h4("Click on the tab 'Model Details' below for a description of the model and this shiny app")
  
    , mainPanel(
    
    ## Box to hold a series of optional changes to the model, which change what parameters users can alter
      box(width = 20
      , column(3
        
  , actionButton("do", "Simulate")
        
  , sliderInput("num_sims"
    , label = h5("Number of simulations")
    , min   = 5
    , max   = 200
    , value = 50
    , step  = 5
    )  
        
  , sliderInput("sim_len"
    , label = h5("Length of simulation")
    , min   = 200
    , max   = 600
    , value = 300
    , step  = 5
    )
      )
   
    , column(3     
        
        , radioButtons("int_type1"
                , h5("Initial intervention strategy")
                , choices = list(
                    "Social Distancing" = 1
                  , "Threshold Based"   = 2)
                  , selected  = 1
                 )
        
        , radioButtons("int_type2"
                , h5("Secondary intervention strategy")
                , choices = list(
                    "Social Distancing" = 1
                  , "Threshold Based"   = 2)
                  , selected  = 1
                 )
 
        , radioButtons("iso"
                , h5("Quarantine of symptomatic infected individuals?")
                , choices = list(
                    "No"  = 1
                  , "Yes" = 2)
                  , selected  = 1
                 )
        
         , conditionalPanel(condition = "input.iso == '2'"
           
   , sliderInput("iso_start"
    , label = h5("Start date of symptomatic infected quarantine")
    , min   = 0
    , max   = 100
    , value = 30
    , step  = 2
    ) 
           
   , sliderInput("iso_length"
    , label = h5("Length of symptomatic infected quarantine")
    , min   = 0
    , max   = 100
    , value = 30
    , step  = 2
    ) 
           
      )
      
    )
         , column(3
    
  , sliderInput("int_start1"
    , label = h5("Start date of FIRST intervention (days since first case entered exposed box)")
    , min   = 0
    , max   = 100
    , value = 30
    , step  = 2
    )  

  , sliderInput("int_length1"
    , label = h5("Length of FIRST intervention (number of days)")
    , min   = 0
    , max   = 100
    , value = 30
    , step  = 2
    )

  , sliderInput("sd_m1"
    , label = h5("Proportion of baseline contact rate (0 - 1)")
    , min   = 0
    , max   = 1
    , value = 0.3
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
        
        ) , column(3
          
    ## dates get updated according to choice in int_start1: see server.R
  , sliderInput("int_start2"
    , label = h5("Start date of SECOND intervention (days since first case entered exposed box)")
    , min   = 50
    , max   = 200
    , value = 80
    , step  = 2
    ) 
        
    ## dates get updated according to choice in int_start1: see server.R       
  , sliderInput("int_length2"
    , label = h5("Length of SECOND intervention (number of days)")
    , min   = 0
    , max   = 200
    , value = 150
    , step  = 2
    )
          
  , sliderInput("sd_m2"
    , label = h5("Proportion of baseline contact rate (0 - 1)")
    , min   = 0
    , max   = 1
    , value = 0.3
    , step  = 0.05
    )  
    
      , conditionalPanel(condition = "input.int_type2 == '2'"
        
  , sliderInput("t_s2"
    , label = h5("Threshold quantity: number of daily hospitalized cases before intervention STARTS")
    , min   = 0
    , max   = 100
    , value = 15
    , step  = 1
    )  
        
  , sliderInput("t_e2"
    , label = h5("Threshold quantity: number of daily hospitalized cases before intervention ENDS")
    , min   = 0
    , max   = 100
    , value = 2
    , step  = 1
    )  
        
      )  
        )
          )
      
  ## Box to for the parameters
  , box(width = 12
      
  , tabsetPanel(type = "tabs"
    
    , tabPanel("Dynamics"
      , plotOutput("graph1")
    #  , plotOutput("graph2")
      )
    
    , tabPanel("Summary Statistics"
      , plotOutput("graph3")
   #   , plotOutput("graph4")
      )
    
    , tabPanel("Model Details"
       , h4("Click the link below to see details about the model. Click refresh at the top of this page to return.")
       , tags$a(href='info.pdf', 'Model Details')
      )
    
    , tabPanel("Download Current Run"
      , downloadButton("datadown", "Download")
    )
    
      )
    )
  )
)
  
        

        
        
  