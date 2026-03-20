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


Slope <- 0.4

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
