library(here)
source(file = here::here("functions/fun.R"))
source(file = here::here("functions/fun_MCSim.R"))
f_load_libraries_colors()

# 1. Data ----

l_Experiments_carac <- f_get_experiments()
df_carac_exp <- l_Experiments_carac$df_carc_exp

l_Experiments <- c(
  "Bart2019_XP1_1", "Bart2019_XP1_2",  
  #"Bart2019_XP2",
  "Bart2019_XP3",   
  "Bart2019_XP4_1", "Bart2019_XP4_2",
  "Bart2019_XP5", "Bart2020",
  #"Gollot2026_lufa"
  "Gollot2026_dens_D1", "Gollot2026_dens_D2", "Gollot2026_dens_D2b",
  "Gollot2026_dens_D5", "Gollot2026_dens_D7"
)

df_data <- f_import_data_DEB(l_Experiments, Adults_alone = FALSE) |> 
  mutate(No_experiment = No_sim)
df_data_long <- df_data |>
  dplyr::select(c(Time, Weight, Reproduction, No_sim, No_experiment)) |> 
  pivot_longer(
    cols = -c("Time", "No_sim", "No_experiment"),
    names_to = "Variable",
    values_to = "Value"
  ) 

df_No_sim <- df_data |> 
  mutate(No_experiment = No_sim) |> 
  dplyr::select(ID_experiment, No_experiment)

time_limits <- df_data %>%
  mutate(No_experiment = as.numeric(No_experiment)) |> 
  group_by(No_experiment) %>%
  summarise(time_cutoff = max(Time, na.rm = TRUE))

# text_simulations <- f_In_simulation_MC(l_Experiments, TRUE, TRUE)
# 
# writeLines(
#   text_simulations,
#   here::here("mod/DEB_MonteCarlo/Text_simulations.in")
# )

# 2. MC ----

File_path <- here::here("mod/DEB_SASLI")
seed <- 121212
Nb_Iter <- 200
RTOL <- 1e-4
ATOL <- 1e-3
Adults_alone <- FALSE
OM_diff <- TRUE



# Simulations simples ----
{
  text_param <- "
  
    pAm = 140;
    Fm = 0.642;
    r_pAm_pM = 0.1;
    kap = 0.899;
    Ehp = 198.68;
    E_coc = 63.7;
    rOM_ClxHorse = 0.00001;
    f_Slope = 1.86;
  "
  
  dv <- 2
  w = 8000
  
  Endpoints_print_sim <- Endpoints_print <- "Length_struct,
    Energy"
  
  f_In_SimpleSimulation(File_path, RTOL, ATOL, text_param, l_Experiments, Endpoints_print_sim, Adults_alone, OM_diff)
  
  system("cd mod/DEB_SASLI; ./mcsim.DEBcaliSLI DEB_SimpleSim.in")
  
  df_sim <- f_MCSim_read_sim(here::here(paste0(File_path, "/sim.out"))) |> 
    mutate(
      No_experiment = No_sim,
      Weight = dv * Length_struct^3 + Energy/w
      ) |> 
    pivot_longer(
      cols = -c(Time, No_sim, No_experiment),
      names_to = c("Variable"),
      values_to = "Value"
    ) 
  
  Experiments_keep <- c(3, 8)
  df_sim_plot <- df_sim |> 
    filter(No_experiment %in% Experiments_keep)
  df_data_long_plot <- df_data_long |> 
    filter(No_experiment %in% Experiments_keep)
  
  
  p_sim <- ggplot()+
    geom_line(
      data = df_sim_plot,
      mapping = aes(
        x = Time,
        y = Value,
        color = Variable
      )
    )+
    geom_point(
      data = df_data_long_plot, 
      mapping = aes(
        x = Time,
        y = Value
      ),
      color = Nord_polar[4],
      alpha = 0.6
    )+
    facet_wrap(
      ~ Variable + No_experiment,
      scales = "free",
      ncol = length(Experiments_keep)
    )+
    scale_color_manual(values = c(Nord_aurora[3], Nord_aurora[4], Nord_frost[4], Nord_aurora[1], Nord_frost[1]))+
    scale_x_continuous(
      breaks = seq(0, max(df_sim_plot$Time), by = 28),
      minor_breaks = seq(0, max(df_sim_plot$Time), by = 14)
    ) +
    theme_minimal()+
    theme(
      legend.position = "none"
    )
  print(p_sim)
}





