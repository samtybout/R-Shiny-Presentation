library(ggplot2)
library(shiny)

# Simulation Functions ----
sim_path = function(spp, ext, D_start = 1, t_max = 100, D_max = 10^8, i_max = 10^4){
  t = 0
  D = D_start
  i = 0
  record = data.frame(t = t, D = D)
  while (t < t_max && D > 0 && D < D_max && i <= i_max){
    i = i + 1
    t_spp = rexp(1, D*spp)
    t_ext = rexp(1, D*ext)
    if (t_spp < t_ext){
      t = t + t_spp
      D = D + 1
    }
    else{
      t = t + t_ext
      D = D - 1
    }
    step = data.frame(t = t, D = D)
    record = rbind(record, step)
  }
  record
}

sim_plot_base = function(){
  data = data.frame()
  plt = ggplot(data = data)
  plt = plt + theme_minimal()
  plt = plt + xlim(0,10) + ylim(0,100)
  plt = plt + xlab("Time") + ylab("Diversity")
  #plt = plt + sim_plot_line(1,1)
  return(plt)
}

sim_plot_line = function(spp, ext, D_start = 50, t_max = 100, D_max = 10^8, i_max = 10^3){
  data = sim_path(spp, ext, D_start, t_max, D_max, i_max)
  geom_line(data = data, aes(x = t, y = D))
}

multiplot = function(spp, ext, D_start, ntrials){
  plt = sim_plot_base()
  for (trial in 1:ntrials){
    plt = plt + sim_plot_line(spp, ext, D_start)
  }
  plt
}

# Initialize plot ----
plt = sim_plot_base

# Shiny ----
ui <- fluidPage(
  titlePanel("Volatility Simulation"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      sliderInput(inputId = "spp",
                  label = "Speciation Rate",
                  min = 0.01,
                  max = 5,
                  value = 1),
      
      
      sliderInput(inputId = "ext",
                  label = "Extinction Rate",
                  min = 0.01,
                  max = 5,
                  value = 1),
      
      sliderInput(inputId = "n0",
                  label = "Starting Diversity",
                  min = 1,
                  max = 100,
                  value = 25),
      
      sliderInput(inputId = "n_trials",
                  label = "Number of trials",
                  min = 1,
                  max = 10,
                  value = 1),

    ),
    
    mainPanel(
      plotOutput(outputId = "simulation_plot")     
    )
  )
)




# Define server logic ----
server <- function(input, output) {
  
  output$simulation_plot = renderPlot(
    
    multiplot(spp = input$spp, ext = input$ext, D_start = input$n0, ntrials = input$n_trials)
    
  )
  
  
}

# Run the app ----
shinyApp(ui = ui, server = server)