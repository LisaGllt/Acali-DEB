

# Predictions IMD ----


# With the chosen model, we aim to predict the impact of a chronic exposure to imidacloprid on earthworm growth and reproduction in comparison to a non-exposed earthworm.
# 
# Scenario No 1 :
#   
#   -   Versailles meadow "Les Closeaux" in 2024
# -   OM = 3.26%
# -   Density : 3
# -   T°C : From the closest meteo station

path_pred <- file.path(path_F_IMD_AG, "Simulations")

Year <- c(2024)                  # Years to get
Start_day <- as.Date("24-01-01") # Start of simulation
Sim_duration <- 100              # Days
step_sim <- 0.5
No_scenario <- 1
bol_degradation <- 0 # No degradation of the substance in the soil

# Controles, median froger_2023, median pelosi_2021, max silva_2019, max pelosi_2021 (= Former application dose)

#l_concentrations <- c(0, 3, 15, 60, 160) 

l_concentrations <- c(0, 3, 15, 60) 

col_concentrations <- c(Nord_aurora[4], Nord_aurora[3], Nord_aurora[2], Nord_aurora[1], Nord_aurora[5])

f_In_sim_predictions(path_pred, Year, Start_day, Sim_duration, l_concentrations, bol_degradation, step_sim, No_scenario)

# for j in $(seq 1 5); do
# ./mcsim.DEB-IMD DEB_IMD_prediction_1_${j}.in
# done

# for j in $(seq 1 3); do
# ./mcsim.DEB-IMD DEB_IMD_prediction_1_${j}.in
# done

#./mcsim.DEB-IMD DEB_IMD_prediction_1_4.in

df_Temp_Closeaux <- f_get_temperatures_Closeaux(Year, Start_day, Sim_duration)


f_Tc <- function(Texp){
  TAH = 28750
  TH = 293.2  
  TA = 7976   
  Tref = 293.15
  
  sA  = exp(TA/Tref-TA/(Texp+273.15));
  srH = (1 + exp(TAH/TH-TAH/Tref))/(1+exp(TAH/TH-TAH/(Texp+273.15)))
  Tc  = sA*((Texp+273.15>=Tref)*srH + (Texp+273.15<Tref))
  return(Tc)
  
}

df_Temp_Closeaux <- df_Temp_Closeaux |> 
  mutate(
    Tc = f_Tc(air_temp_mean_day)
  )

pT <- ggplot(
  data = df_Temp_Closeaux,
  aes(
    x = Time,
    y = air_temp_mean_day
  )
) +
  geom_line(color = Nord_polar[1])+
  labs(
    x = "Time (days)",
    y = "Temperature (°C)"
  )+
  theme_minimal()

pTc <- ggplot(
  data = df_Temp_Closeaux,
  aes(
    x = Time,
    y = Tc
  )
) +
  geom_line(color = Nord_polar[1])+
  labs(
    x = "Time (days)",
    y = "Tc"
  )+
  theme_minimal()


df_sim <- f_Read_predictions(path_pred, No_scenario, l_concentrations, Sim_duration, step_sim)

df_sim_W <- df_sim |>
  filter(Endpt == "Weight")

df_sim_R <- df_sim |>
  filter(Endpt == "Reproduction")

df_sim_S <- df_sim |>
  filter(Endpt == "Stress")

df_sim_OM <- df_sim |>
  filter(Endpt == "Organic_matter")


alpha_ribbon <- 0.15

pW <- ggplot(
  data = df_sim_W,
  aes(
    x = Time,
    y = predict.endpoint
  )
) +
  labs(
    x = "Time (days)",
    y = "Weight (g)"
  )+
  geom_line(
    aes(
      color = as.factor(Ce)
    )
  )+
  geom_ribbon(
    aes(
      x = Time,
      ymax = up,
      ymin = low,
      fill = as.factor(Ce)
    ),
    alpha = alpha_ribbon
  )+
  scale_color_manual(values = col_concentrations, name = "Concentration (ng/g)")+
  scale_fill_manual(values = col_concentrations, name = "Concentration (ng/g)")+
  theme_minimal()

pOM <- ggplot(
  data = df_sim_OM,
  aes(
    x = Time,
    y = predict.endpoint,
    color = as.factor(Ce)
  )
) +
  labs(
    x = "Time (days)",
    y = "Organic matter (%)"
  )+
  geom_line()+
  theme_minimal()

pS <- ggplot(
  data = df_sim_S,
  aes(
    x = Time,
    y = predict.endpoint,
    color = as.factor(Ce)
  )
) +
  labs(
    x = "Time (days)",
    y = "Stress"
  )+
  geom_line()+
  theme_minimal()

p <- pT + pW + pOM + pS + plot_layout(ncol=1)

