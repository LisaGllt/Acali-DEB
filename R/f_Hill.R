library(here)
source(file = here::here("functions/fun.R"))
source(file = here::here("functions/fun_MCSim.R"))
f_load_libraries_colors()

EC <- 0.52
Slope <- 0.208

f_Hill <- function(x, EC, Slope){
  Y = x^Slope/(EC^Slope + x^Slope)
  return(Y)
}

df <- data.frame(
  x = seq(0,1,0.001)
) |> 
  mutate(
    y = f_Hill(1000*x, EC, Slope)
  )


ggplot(
  data = df,
  aes(
    x = x,
    y = y
  )
) +
geom_line() +
#scale_x_log10()+
theme_minimal()


EC <- 0.52
Slope <- 0.608
Slope <- 3

f_Hill <- function(x, EC, Slope){
  Y = exp(Slope*(x-EC))/(1+exp(Slope*(x-EC)))
  return(Y)
}

df <- data.frame(
  x = seq(0,1,0.001)
) |> 
  mutate(
    y = f_Hill(x, EC, Slope)
  )


ggplot(
  data = df,
  aes(
    x = x,
    y = y
  )
) +
  geom_line() +
  #scale_x_log10()+
  theme_minimal()


# Homographie tq F(0)=0 et f(1)=1


Slope <- 0.08

f_Hill <- function(x, EC, Slope){
  Y = x/(Slope * (1-x) + x)
  return(Y)
}

df <- data.frame(
  x = seq(0,1,0.001)
) |> 
  mutate(
    y = f_Hill(x, EC, Slope)
  )


ggplot(
  data = df,
  aes(
    x = x,
    y = y
  )
) +
  geom_line() +
  #scale_x_log10()+
  theme_minimal()





Slope <- 0.4

f_Hill <- function(x, EC, Slope){
  Y = x/(1 + x)
  return(Y)
}

df <- data.frame(
  x = seq(0,2,0.001)
) |> 
  mutate(
    y = f_Hill(x, EC, Slope)
  )


ggplot(
  data = df,
  aes(
    x = x,
    y = y
  )
) +
  geom_line() +
  lims(x=c(0,2), y=c(0,1))+
  #scale_x_log10()+
  theme_minimal()





z <- 0.1

f_Hill <- function(x, z){
  Y = ifelse(x <= z, z, x)
  return(Y)
}

df <- data.frame(
  x = seq(0,1,0.001)
) |> 
  mutate(
    y = f_Hill(x, z)
  )


ggplot(
  data = df,
  aes(
    x = x,
    y = y
  )
) +
  geom_line() +
  lims(x=c(0,1), y=c(0,1))+
  #scale_x_log10()+
  theme_minimal()








Slope1 <- 1
Slope2 <- 0.8

f_Hill <- function(x, Slope1, Slope2){
  Y = ifelse(x <= 0.5, x*Slope1/((Slope1-1) * x + 1), x/(Slope1 * (1-x) + x)) 
  return(Y)
}

df <- data.frame(
  x = seq(0,1,0.01)
) |> 
  mutate(
    y = f_Hill(x, Slope1, Slope2)
  )


ggplot(
  data = df,
  aes(
    x = x,
    y = y
  )
) +
  geom_line() +
  #scale_x_log10()+
  theme_minimal()




K <- 1
Alpha <- 2

f_Hill <- function(x, K, Alpha){
  Y = x^Alpha/(x^Alpha+K*(1-x)^Alpha) 
  return(Y)
}

df <- data.frame(
  x = seq(0,1,0.05)
) |> 
  mutate(
    y = f_Hill(x, K, Alpha)
  )


ggplot(
  data = df,
  aes(
    x = x,
    y = y
  )
) +
  geom_line() +
  #scale_x_log10()+
  theme_minimal()


b <- 10

f_Hill <- function(x, K, Alpha){
  Y = (1+b)*x/(1+b*x) 
  return(Y)
}

df <- data.frame(
  x = seq(0,1,0.01)
) |> 
  mutate(
    y = f_Hill(x, K, Alpha)
  )


ggplot(
  data = df,
  aes(
    x = x,
    y = y
  )
) +
  geom_line() +
  #scale_x_log10()+
  theme_minimal()

a <- 1
b <- 5

f_Hill <- function(x, K, Alpha){
  Y = x^a/(x^a+(1-x)^b) 
  return(Y)
}

df <- data.frame(
  x = seq(0,1,0.01)
) |> 
  mutate(
    y = f_Hill(x, K, Alpha)
  )


ggplot(
  data = df,
  aes(
    x = x,
    y = y
  )
) +
  geom_line() +
  #scale_x_log10()+
  theme_minimal()



{
z <- 0.7
a <- 0.3
K <- 0.378
x_s = z * K /(1+K*z-z)
b = z - a * x_s

x <- 0.2

f_Hill <- function(x, a, z, K){
  Y = x/(K*(1-x)+x) 
  x_s = z * K /(1+K*z-z)
  b = z - a * x_s
  Y_real = ifelse(Y<=z, a*x+b, Y)
  return(Y_real)
}

df <- data.frame(
  x = seq(0,1,0.01)
) |> 
  mutate(
    y = f_Hill(x, a, z, K),
    y2 = a*x+b
  )


ggplot() +
  geom_line(
    data = df,
    aes(
      x = x,
      y = y
    ),
    color = "blue"
  ) +
  #scale_x_log10()+
  theme_minimal()
}



alpha <- 2
delta <- 7
t_A <- 100
kapJ <- 0.6
kapA <- 0.8

f_T <- function(x, alpha){
  
  y <- ifelse(x < 0, 0, ifelse(x > 1, 1, x^alpha/(x^alpha +(1-x)^alpha)))

  return(y)
}

df_fT <- data.frame(
  x = seq(0, 1, 0.01)
) |> 
  mutate(y = f_T(x, alpha))

ggplot() +
  geom_line(
    data = df_fT,
    aes(
      x = x,
      y = y
    ),
    color = "blue"
  ) +
  theme_minimal()

t <- seq(0,300, 0.1)

df <- data.frame(
  t = t
) |> 
  mutate(
    x = (t -(t_A))/(delta),
    y_T = f_T(x, alpha),
    kap = f_T(x, alpha)*kapA + (1-f_T(x, alpha))*kapJ
  )

ggplot() +
  geom_line(
    data = df,
    aes(
      x = x,
      y = y_T
    ),
    color = Nord_polar[1]
  ) +
  lims(x=c(0,10))+
  theme_minimal()


ggplot() +
  geom_abline(
    intercept = kapA,
    slope = 0,
    color = Nord_aurora[1],
    linetype = "dashed"
  )+
  geom_abline(
    intercept = kapJ,
    slope = 0,
    color = Nord_aurora[4],
    linetype = "dashed"
  )+
  
  annotate(
    geom = "text",
    x = 75,
    y = 0.63,
    fontface = "bold",
    label = "kapJ",
    color = Nord_aurora[4],
    size = 4
  ) +
  
  annotate(
    geom = "text",
    x = 75,
    y = 0.77,
    fontface = "bold",
    label = "kapA",
    color = Nord_aurora[1],
    size = 4
  ) +

  geom_line(
    data = df,
    aes(
      x = t,
      y = kap
    ),
    color = Nord_polar[1],
    linewidth = 1
  ) +
  
  geom_vline(
    xintercept = t_A,
    color = Nord_frost[2],
    linetype = "dashed",
    linewidth = 1
  )+
  
  annotate(
    geom = "text",
    x = 125,
    y = 0,
    fontface = "bold",
    label = "tA tq Eh(tA)=Ehp",
    color = Nord_frost[2],
    size = 4
  ) +
  
  lims(x=c(75,150), y = c(0, 1))+
  labs(y = "Kapreal", x = "Time (days)")+
  theme_minimal()


