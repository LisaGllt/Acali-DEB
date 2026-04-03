library(here)
source(file = here::here("functions/fun.R"))
f_load_libraries_colors() 

# Data ----

df_data <- read.csv2(here::here("data/Data_mineralisation.csv")) |> 
  mutate(exposition = paste0(polymere, concentration))

mod_lm <- lm(CO2 ~ 0 + exposition:jour, data = df_data) 
mod_null <- lm(CO2 ~ 0 + jour, data = df_data)

anova(mod_null, mod_lm)
emtrends(mod_lm, pairwise ~ exposition, var = "jour")



df_data_J11 <- read.csv2(here::here("data/Data_mineralisation.csv")) |> 
  mutate(exposition = paste0(polymere, concentration)) |> 
  filter(jour <=11)

mod_lm_J11 <- lm(CO2 ~ 0 + exposition:jour, data = df_data_J11) 
mod_null_J11 <- lm(CO2 ~ 0 + jour, data = df_data_J11)

anova(mod_null_J11, mod_lm_J11)
emtrends(mod_lm_J11, pairwise ~ exposition, var = "jour")

plot(mod_lm)




