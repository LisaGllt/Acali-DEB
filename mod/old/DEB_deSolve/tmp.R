
f_Tc <- function(Texp){ # Texp en °C
  
  TA <- 7.97600e+03
  TH <- 2.93200e+02
  TAH <- 2.87500e+04
  Tref <- 2.93150e+02 
  
  # Coefficient for temperature correction (Texp in °C)
  sA = exp(TA/Tref-TA/(Texp+273.15))
  srH = (1 + exp(TAH/TH-TAH/Tref))/(1+exp(TAH/TH-TAH/(Texp+273.15)))
  Tc  = sA*((Texp+273.15>=Tref)*srH + (Texp+273.15<Tref))
  
  return(Tc)
}

TA <- 7.97600e+03
TH <- 2.93200e+02
TAH <- 2.87500e+04
Tref <- 2.93150e+02 

Temp <- seq(-10, 40, 0.1)
Tc <- f_Tc(Temp)

df_Tc <- data.frame(
  Temperature = Temp, 
  Factor = Tc
) 

p <- ggplot()+
  geom_line(
    data = df_Tc,
    aes(
      x = Temperature,
      y = Factor
    ),
    color = Nord_polar[4]
  )+
  # geom_vline(
  #   xintercept = (TA - 273.15),
  #   color = col_red
  # )+
  geom_vline(
    xintercept = (TH - 273.15),
    color = Nord_aurora[1],
    linewidth = 0.5
  )+
  geom_vline(
    xintercept = (Tref - 273.15),
    color = Nord_aurora[2],
    linewidth = 0.5
  )+
  # geom_vline(
  #   xintercept = (TAH - 273.15),
  #   color = col_red
  # )+
  labs(
    x = "Temperature (°C)",
    y = "Correction factor"
  )+
  lims(x=c(17, 23))+
  theme_minimal()
p










