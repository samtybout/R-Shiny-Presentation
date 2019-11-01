# Plotting function ----
lineplot = function(slope, intercept, color){
  x = seq(0,10,0.1)
  y = slope*x + intercept
  plot(x, y, type = "l", col = color, lwd = 3)
}

# Shiny ----
ui <- fluidPage(
  titlePanel("Simple Plot Example"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      sliderInput(inputId = "slope",
                  label = "Slope",
                  min = -10,
                  max = 10,
                  value = 1),
      
      sliderInput(inputId = "intercept",
                  label = "Intercept",
                  min = 0,
                  max = 10,
                  value = 0),
      
      selectInput(inputId = "color",
                  label = "Color",
                  choices = c("Red", "Blue", "Yellow", "Green"))
    ),
    
    mainPanel(
      plotOutput(outputId = "projection_plot")     
    )
  )
)


# Define server logic ----



server <- function(input, output) {
  
  output$projection_plot = renderPlot(
    
    lineplot(input$slope, input$intercept, input$color)
    
  )
  
}

# Run the app ----
shinyApp(ui = ui, server = server)