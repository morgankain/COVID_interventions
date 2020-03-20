
needed_packages <- c(
  "shiny"
, "shinythemes"
, "shinydashboard"
)

## Check if the packages are installed. *If they are not install them*, then load them
if (length(setdiff(needed_packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(needed_packages, rownames(installed.packages())))  
}

lapply(needed_packages, require, character.only = TRUE)


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
  
    , h4("Click on the tab 'Model Details' below for a description of the model and this shiny app")
  
    , mainPanel(
    
    ## Box to hold a series of optional changes to the model, which change what parameters users can alter
      box(width = 20
      , column(8
        
  , sliderInput("num_sims"
    , label = h5("Number of simulations")
    , min   = 5
    , max   = 200
    , value = 50
    , step  = 5
    )  
        
        , radioButtons("int_type"
                , h5("Intervention strategy")
                , choices = list(
               #    "None"              = "None"
                    "Social Distancing" = 1
                  , "Threshold Based"   = 2)
                  , selected  = 1
                 )
        
      , conditionalPanel(condition = "input.int_type == '1'"
        
  , sliderInput("int_start"
    , label = h5("Start date of intervention (days since first case entered exposed box)")
    , min   = 0
    , max   = 100
    , value = 30
    , step  = 2
    )  
        
  , sliderInput("int_size"
    , label = h5("Proportion of baseline contact rate (0 - 1)")
    , min   = 0
    , max   = 1
    , value = 0.3
    , step  = 0.05
    )  
        
  , sliderInput("int_len"
    , label = h5("Length of intervention (number of days)")
    , min   = 0
    , max   = 100
    , value = 30
    , step  = 2
    )
       
      )
         
      , conditionalPanel(condition = "input.int_type == '2'"
        
  , sliderInput("int_start_t"
    , label = h5("Threshold quantity: number of daily hospitalized cases before intervention STARTS")
    , min   = 0
    , max   = 100
    , value = 15
    , step  = 1
    )  
        
  , sliderInput("int_end_t"
    , label = h5("Threshold quantity: number of daily hospitalized cases before intervention ENDS")
    , min   = 0
    , max   = 100
    , value = 2
    , step  = 1
    )  
        
  , sliderInput("int_size_t"
    , label = h5("Proportion of baseline contact rate (0 - 1)")
    , min   = 0
    , max   = 1
    , value = 0.1
    , step  = 0.05
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
  
        

        
        
  