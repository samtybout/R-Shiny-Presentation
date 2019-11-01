library(ggplot2)
library(shiny)

# Probability Density Functions ----

pnt = function(n0,n,t,spp,ext,alpha = NULL, beta = NULL){
  
  if(n0 < 0){
    print("Initial diversity cannot be negative")
    return(NaN)
  }
  
  if(spp < 0 | ext < 0){
    print("Speciation and extinction rates cannot be negative")
    return(NaN)
  }
  
  # If initial diversity is 0, diversity will always be 0
  if(n0 == 0){
    if(n == 0){
      return(1)
    }
    else{
      return(0)
    }
  }
  
  if(n < 0){
    return(0)
  }
  
  if(t == 0){
    return(pnt_F(n0,n))
  }
  
  if(spp == 0 & ext == 0){
    return(pnt_F(n0,n))
  }
  
  if(spp == 0 & ext > 0){
    return(pnt_D(n0,n,t,ext))
  }
  
  if(spp > 0 & ext ==0){
    return(pnt_E(n0,n,t,spp))
  }
  
  if(spp == ext){
    return(pnt_A(n0,n,t,spp))
  }
  
  return(pnt_B(n0,n,t,spp,ext))
}

pnt_F = function(n0,n){
  # Speciation = Extinction = 0
  # When both rates are 0, diversity remains constant indefinitely
  if(n == n0){
    return(1)
  }
  else{
    return(0)
  }
}

pnt_D = function(n0,n,t,ext){
  # Speciation = 0
  # Extinction > 0
  # A pure death process; each extinction can be treated as an independent random
  # variable with an exponential distribution
  if (n > n0){
    return(0)
  }
  else{
    p_ext = 1 - exp(-ext * t)
    return(choose(n0,n) * (p_ext^(n0 - n)) * ((1 - p_ext)^n))
  }
}

pnt_E = function(n0,n,t,spp){
  # Speciation > 0
  # Extinction = 0
  # Yule process; equation from Stochastic Processes by Sheldon Ross, 1983
  if (n < n0){
    return(0)
  }
  else{
    return(
      choose(n - 1, n0 - 1) * exp(-1 * spp * t * n0) * (1 - exp(-spp * t))^(n - n0)
    )
  }
}

pnt_A = function(n0,n,t,vol){
  # Speciation = Extinction
  # Equations from Raup 1985
  
  # Equation A12
  if (n == 0){
    return(
      ((vol * t) / (1 + vol * t))^n0
    )
  }
  
  # Equation A16
  else{
    j = 1:min(n0,n)
    return(
      (((vol * t) / (1 + vol * t))^(n0 + n)) *
        sum(choose(n0,j) * choose(n-1,j-1) * ((vol * t)^(-2 * j)))
    )
  }
}

pnt_B = function(n0,n,t,spp,ext){
  # Speciation =/= Extinction
  # Equations from Raup 1985
  
  # Equation A13
  alpha = (ext*(expm1((spp-ext)*t)))/(spp*exp((spp-ext)*t)-ext)
  
  # Equation A14
  if(n == 0){
    return(alpha^n0)
  }
  
  # Equation A18
  else{
    beta = alpha * spp / ext
    j = 0:min(n0,n)
    p = sum(
      choose(n0, j) * choose(n0+n-j-1, n0-1) * 
        (alpha^(n0-j)) * (beta^(n-j)) * ((1-alpha-beta)^j)
    )
    if (p < 0 | p > 1){
      return(0)
    }
    else{
      return(p)
    }
  }
}
p_dist = function(n0,max,t,spp,ext){
  p = c()
  d = 0:max
  for (n in d){
    p = c(p,pnt(n0,n,t,spp,ext))
  }
  return(data.frame(D = d, p = p))
}

mode_n = function(n0,t,spp,ext){
  expected = n0 * exp((spp - ext) * t)
  range = 0:(round(expected) * 2)
  return(mode_n_sub(n0,t,spp,ext,range))
}

mode_n_sub = function(n0,t,spp,ext,range){
  bot = min(range)
  top = max(range)
  # print(c(bot,top))
  if(length(range) < 3){
    p_low = pnt(n0,bot,t,spp,ext)
    p_high = pnt(n0,top,t,spp,ext)
    if(p_low > p_high){
      return(bot)
    }
    if(p_low < p_high){
      return(top)
    }
    else{
      return((bot+top)/2)
    }
  }  
  
  mid = (bot + top) %/% 2
  low = (bot + mid) %/% 2
  high = (mid + top) %/% 2
  
  p_low = pnt(n0,low,t,spp,ext)
  p_high = pnt(n0,high,t,spp,ext)
  
  if(p_low > p_high){
    return(mode_n_sub(n0,t,spp,ext,bot:mid))
  }
  else{
    return(mode_n_sub(n0,t,spp,ext,mid:top))
  }
}

conf_pnt = function(n0,t,spp,ext,int = 0.95){
  center = mode_n(n0,t,spp,ext)
  p = pnt(n0,center,t,spp,ext)
  lower = center - 1
  upper = center + 1
  while(p < int){
    p_upper = pnt(n0,upper,t,spp,ext)
    p_lower = pnt(n0,lower,t,spp,ext)
    if(p_upper > p_lower){
      p = p + p_upper
      upper = upper + 1
    }
    else{
      p = p + p_lower
      lower = lower - 1
    }
  }
  return(data.frame(lower = lower + 1, upper = upper - 1, p = p))
}

conf_series = function(n0,t_series,spp,ext,int = 0.95){
  data = data.frame()
  for (t in t_series){
    step = cbind(t = t, conf_pnt(n0,t,spp,ext,int))
    data = rbind(data, step)
  }
  return(data)
}

density_plot = function(n0,t_series,spp,ext,n_max,interpolate=FALSE, contours = FALSE,
                        confints = c(), breaks = c(0.05,0.1,0.25,0.5,0.9,0.95)){
  data = data.frame()
  for(time in t_series){
    for(n in 0:n_max){
      data = rbind(data, data.frame(t = time, n = n, p = pnt(n0,n,time,spp,ext)))
    }
  }
  data$p[data$p == 0] = min(data$p)
  # print(data$p)
  # print(log(abs(data)))
  plt = ggplot(data,aes(x = t, y = n))
  plt = plt + geom_raster(aes(fill = log(p)),interpolate=interpolate)
  plt = plt + scale_fill_gradient()
  if(contours){
    plt = plt + geom_contour(aes(z = p),colour = "white",breaks = breaks,
                             show.legend = TRUE)
  }
  for(interval in confints){
    conf_data = conf_series(n0,t_series,spp,ext,interval)
    plt = plt + geom_line(data=conf_data,aes(x = t, y = lower),colour = "white") + 
      geom_line(data=conf_data,aes(x = t, y = upper),colour = "white")
  }
  plt = plt + theme_minimal()
  plt
}


# Shiny ----
ui <- fluidPage(
  titlePanel("title panel"),
  
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
                  max = 50,
                  value = 25)
    ),
    
    mainPanel(
         plotOutput(outputId = "projection_plot")     
      )
  )
)


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

# Define server logic ----
server <- function(input, output) {
  
  output$projection_plot = renderPlot(
    
    density_plot(input$n0,seq(0.1,30),input$spp,input$ext,50,contours = FALSE, interpolate = FALSE)
    
  )
  
}

# Run the app ----
shinyApp(ui = ui, server = server)