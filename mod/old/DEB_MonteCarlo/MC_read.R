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

df_data <- f_import_data_DEB(l_Experiments, Adults_alone = TRUE) |> 
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

File_path <- here::here("mod/DEB_MonteCarlo")
seed <- 121212
Nb_Iter <- 200
RTOL <- 1e-7
ATOL <- 1e-6
Adults_alone <- TRUE
OM_diff <- TRUE

Endpoints_print <- "Weight,
    Reproduction"

text_param <- "
  rOM_ClxHorse = 0.002;
  Fm = 0.3;
  v = 0.18;
  pAm = 100;
  r_pAm_pM = 0.8;
  kap = 0.7;
  Eg = 4183;
  kJ = 0.0001;      #= 0.002793
  Ehp = 80;
  E_coc = 100;
  "

text_distrib <- "
Distrib (rOM_ClxHorse, LogUniform,               0.001,      0.01);
Distrib (Fm,           Uniform,               0.2,      0.4);
Distrib (pAm,          Uniform,               50,      400);
Distrib (r_pAm_pM,     Uniform,               0.6,      1);
Distrib (kap,          Uniform,               0.4,     0.8);
Distrib (kJ,           LogUniform,            0.0001,     0.003);
Distrib (Ehp,          Uniform,              50,      120);
Distrib (E_coc,       Uniform,              30,      100);
"

Nb_param <- 7

f_In_MC_tot(File_path, seed, Nb_Iter, RTOL, ATOL, Endpoints_print, text_distrib, l_Experiments, Adults_alone, OM_diff)

# makemcsim DEBcali.model 
# for j in $(seq 1 12); do
# ./mcsim.DEBcali DEB_MC_${j}.in
# done

# for j in $(seq 1 12); do
# ./mcsim.DEBcaliS DEB_MC_${j}.in
# done

# for j in $(seq 1 12); do
# ./mcsim.DEBcaliF DEB_MC_${j}.in
# done

times_keep <- seq(0, 400, 1)

df_MC_sim <- f_read_all_MC(path_file, l_Experiments, Nb_param, times_sim, times_keep, time_limits)
saveRDS(df_MC_sim, file = here::here("mod/DEB_MonteCarlo/df_MC_sim_S200_3.rds"))
#df_MC_sim <- readRDS(here::here("mod/DEB_MonteCarlo/df_MC_sim_1000_pM.rds"))
cat("Simulations réussies :", length(unique(df_MC_sim$Iter)))

# Plot de certaines itérations
{ 
  Iter_keep <- sample(seq(1, Nb_Iter, 1), 29, replace = FALSE)
  df_data_plot <- df_MC_sim |> 
    filter(Iter %in% Iter_keep)
  
  p <- ggplot()+
    geom_line(
      data = df_data_plot,
      mapping = aes(
        x = Time,
        y = Value,
        group = Iter,
        color = Variable
      )
    )+
    geom_point(
      data = df_data_long, 
      mapping = aes(
        x = Time,
        y = Value
      ),
      color = Nord_polar[4],
      alpha = 0.6
    )+
    facet_grid(
      rows = vars(Variable), 
      cols = vars(No_experiment),
      scales = "free"
      )+
    scale_color_manual(values = c(Nord_aurora[4], Nord_frost[4]))+
    scale_x_continuous(
      breaks = seq(0, max(df_data_plot$Time), by = 28),
      minor_breaks = seq(0, max(df_data_plot$Time), by = 14)
    ) +
    theme_minimal()+
    theme(
      legend.position = "none"
    )
  
  p
  rm(df_data_plot)
}

df_errors <- f_MC_RMSE_calc(df_data, df_MC_sim)

Nb_best <-3

Best_RMSE_W_Iter <- (df_errors |> 
  arrange(RMSE_Weight))$Iter[1:Nb_best]

Best_RMSE_R_Iter <- (df_errors |> 
  arrange(RMSE_Reproduction))$Iter[1:Nb_best]

Best_RMSE_tot_Iter <- (df_errors |> 
  arrange(NRMSE_tot))$Iter[1:Nb_best]

Best_RMSE_W_Iter
Best_RMSE_R_Iter
Best_RMSE_tot_Iter


# Plot de certaines itérations
{ 
  Iter_keep <- Best_RMSE_tot_Iter
  #Iter_keep <- c(65)
  df_data_plot <- df_MC_sim |> 
    filter(Iter %in% Iter_keep)
  
  p_tot <- ggplot()+
    geom_line(
      data = df_data_plot,
      mapping = aes(
        x = Time,
        y = Value,
        group = Iter,
        color = Variable
      )
    )+
    geom_point(
      data = df_data_long, 
      mapping = aes(
        x = Time,
        y = Value
      ),
      color = Nord_polar[4],
      alpha = 0.6
    )+
    facet_grid(
      rows = vars(Variable), 
      cols = vars(No_experiment),
      scales = "free"
    )+
    scale_color_manual(values = c(Nord_aurora[4], Nord_frost[4]))+
    scale_x_continuous(
      breaks = seq(0, max(df_data_plot$Time), by = 28),
      minor_breaks = seq(0, max(df_data_plot$Time), by = 14)
    ) +
    theme_minimal()+
    theme(
      legend.position = "none"
    )
  print(p_tot)
  rm(df_data_plot)
}




df_param_best <- df_MC_sim |> 
  filter(Iter %in% Iter_keep) |> 
  filter(Time == 14 & Variable == "Weight" & No_experiment == 1)
df_param_best


# Plot de certaines itérations
{ 
  Iter_keep2 <- c(106)
  #Iter_keep <- c(65)
  df_data_plot2 <- df_MC_sim |> 
    filter(Iter %in% Iter_keep2)
  
  p_tot2 <- ggplot()+
    geom_line(
      data = df_data_plot2,
      mapping = aes(
        x = Time,
        y = Value,
        group = Iter,
        color = Variable
      )
    )+
    geom_point(
      data = df_data_long, 
      mapping = aes(
        x = Time,
        y = Value
      ),
      color = Nord_polar[4],
      alpha = 0.6
    )+
    facet_grid(
      rows = vars(Variable), 
      cols = vars(No_experiment),
      scales = "free"
    )+
    scale_color_manual(values = c(Nord_aurora[4], Nord_frost[4]))+
    scale_x_continuous(
      breaks = seq(0, max(df_data_plot2$Time), by = 28),
      minor_breaks = seq(0, max(df_data_plot2$Time), by = 14)
    ) +
    theme_minimal()+
    theme(
      legend.position = "none"
    )
  p_tot2
  rm(df_data_plot2)
}

# Simulations simples ----
{
  text_param <- "
  rOM_ClxHorse = 0.002;
  Fm = 0.3;
  v = 0.18;
  pAm = 100;
  r_pAm_pM = 0.8;
  kap = 0.7;
  Eg = 4183;
  kJ = 0.0001;      #= 0.002793
  Ehp = 80;
  E_coc = 100;
  "
  
  Endpoints_print_sim <- Endpoints_print <- "Weight,
    Reproduction,
    Maturity,"
  
  f_In_SimpleSimulation(File_path, RTOL, ATOL, text_param, l_Experiments, Endpoints_print_sim, Adults_alone, OM_diff)
  
  system("cd mod/DEB_MonteCarlo; ./mcsim.DEBcaliS DEB_SimpleSim.in")
  
  df_sim <- f_MCSim_read_sim(here::here(paste0(File_path, "/sim.out"))) |> 
    mutate(No_experiment = No_sim) |> 
    pivot_longer(
      cols = -c(Time, No_sim, No_experiment),
      names_to = c("Variable"),
      values_to = "Value"
    ) 
  
  Experiments_keep <- c(1, 3, 5, 8, 9)
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
      nrow = 3
    )+
    scale_color_manual(values = c(Nord_aurora[3], Nord_aurora[4], Nord_frost[4]))+
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





