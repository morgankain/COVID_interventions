
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
  
    , h5("Model development and coding: Marissa Childs, Devin Kirk, Morgan Kain, Nicole Nova")
    , h5("Shiny app:", tags$a("Morgan Kain", href = "http://www.morgankain.weebly.com"))
    , h5("Planning and parameter estimate search: Mallory Harris, Lisa Couper")
  
    , h4("Click on the tab 'Model Details' below for a description of the model and this shiny app")
  
    , mainPanel(
    
    ## Box to hold a series of optional changes to the model, which change what parameters users can alter
      box(width = 20
      , column(8
        
        , radioButtons("int_type"
                , h5("Intervention strategy")
                , choices = list(
                    "None"              = "None"
                  , "Shelter in Place"  = "STP"
                  , "Social Distancing" = "SD")
                  , selected  = "None")
        
      , conditionalPanel(condition = "input.int_type != 'None'"
        
          , numericInput("int_size"
           , label = h5("Magnitude of intervention (0 - 1)")
           , value = 0)       
        
          , numericInput("int_start"
           , label = h5("Start date of intervention (day since first case enters exposed box)")
           , value = 0)          
        
          , numericInput("int_len"
           , label = h5("Length of intervention (number of days)")
           , value = 0)  
        
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
        

        
        
  