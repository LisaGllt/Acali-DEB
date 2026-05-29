# MCSim to R ----

f_get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

f_get_experiments <- function() {
  l_Experiments_tot <- c(
    "Bart2019_XP1_1",  
    "Bart2019_XP1_2",  
    "Bart2019_XP2",
    "Bart2019_XP3",   
    "Bart2019_XP4_1",
    "Bart2019_XP4_2",
    "Bart2019_XP5", 
    "Bart2020",
    "Gollot2026_lufa",
    "Gollot2026_dens_D1",
    "Gollot2026_dens_D2",
    "Gollot2026_dens_D2b",
    "Gollot2026_dens_D5",
    "Gollot2026_dens_D7"
  )
  
  l_carac <- c(
    "3 g/ind/14 days - Closeaux soil (15°C)",
    "3 g/ind/14 days - Closeaux soil (15°C)",
    "3 g/ind/14 days - Closeaux soil (15°C)",
    "0 or 3 g/ind/14 days - Closeaux soil (15°C)",
    "1 g/ind/14 days - Closeaux soil (15°C)",
    "1 g/ind/14 days - Closeaux soil (15°C)",
    "1.5 g/ind/14 days - Closeaux soil (15°C)",
    "3 g/ind/14 days - Closeaux soil (15°C)",
    "3 g/ind/14 days - LUFA 2.4 soil (18°C)",
    "3 g/14 days - 1/cosm - Closeaux soil (18°C)",
    "3 g/14 days - 2/cosm - Closeaux soil (18°C)",
    "6 g/14 days - 2/cosm - Closeaux soil (18°C)",
    "3 g/14 days - 5/cosm - Closeaux soil (18°C)",
    "3 g/14 days - 7/cosm - Closeaux soil (18°C)"
  )
  
  df_carc_exp <- data.frame(
    ID_experiment = l_Experiments_tot,
    carac = l_carac
  )
  
  res <- list(
    l_Experiments_tot = l_Experiments_tot,
    l_carac = l_carac,
    df_carc_exp = df_carc_exp
  )
  
  return(res)
}

f_import_data_DEB <- function(l_Experiments, Adults_alone = FALSE, fit = FALSE){
  
  XP_no_R_data <- c("Bart2019_XP1_1", "Bart2019_XP4_1", "Bart2019_XP3", "Bart2020", "Gollot2026_D1")
  
  if(!Adults_alone){
    df_data <- read.csv(here::here("data/Data_DEB_mean.csv")) |> 
      filter((ID_experiment %in% c("Bart2019_XP1_2", "Bart2019_XP4_2"))|!(Status == "A" & Nb_ind == 1) & ID_experiment %in% l_Experiments)
  }else{
    df_data <- read.csv(here::here("data/Data_DEB_mean.csv")) |> 
      filter(ID_experiment %in% l_Experiments)
  }
  
  WHC_clx <- 39.14 # mL (100%)
  WHC_lufa24 <- 47.1 # mL (100%)
  
  Per_C_org_lufa24 <- 1.83 # C %
  # Assumption : There is 58% of carbon in soil organic matter
  tx_OM_clx <- 0.0326 # g/g
  tx_OM_lufa24 <- Per_C_org_lufa24 * 1/0.58/100 # g/g
  tx_OM_horse_dung <- 0.9 # g/g of OM in horse dung (!!!!!!!!!!!!!!!!!!!!!!)
  
  df_data <- df_data |> 
    mutate(
      
      Density = Nb_ind,
      
      # Corresponding humidity calculation (%)
      Moisture = case_when(
        Soil == "Clx" ~ WHC * WHC_clx,
        Soil == "LUFA2.4" ~ WHC * WHC_lufa24
      ),
      
      # g of OM from soil in the entire cosm 
      OM_soil = case_when(
        Soil == "Clx" ~ tx_OM_clx*Soil_w,
        Soil == "LUFA2.4" ~ tx_OM_lufa24*Soil_w
      ),
      
      # g soil of OM from horse dung in the entire cosm 
      OM_horse = Food*tx_OM_horse_dung*Density, 
      
      Soil_type = case_when(
        Soil == "Clx" ~ 1,
        Soil == "LUFA2.4" ~ 2
      )
    )
  
  df_DEB_mean <- df_data |> 
    mutate(Weight = Weight/1000) |> # conversion en g
    dplyr::select(ID_experiment, Time, Weight, Reproduction, Environment, OM_soil, OM_horse, Density, Texp, Soil_type, Soil_w) |> 
    mutate(No_sim = as.factor(match(ID_experiment, unique(ID_experiment))))
  
  if (fit == TRUE){
    df_DEB_mean <- df_DEB_mean %>%
      group_by(ID_experiment) %>%
      mutate(
        Reproduction = if_else(
          ID_experiment %in% XP_no_R_data &
            Time != 0,
          NA_real_,
          Reproduction
        )
      ) %>%
      ungroup()
  }
  
  return(df_DEB_mean)
}

# Monte Carlo ----


f_In_simulation_MC <- function(l_Experiments, Endpoints_print, Adults_alone, OM_diff){
  
  df_DEB_mean <- f_import_data_DEB(l_Experiments, Adults_alone) |> 
    mutate(
      OM_soil = signif(OM_soil, 3),
      OM_horse = signif(OM_horse, 3)
    )  |> 
    mutate(across(where(is.numeric), ~ replace_na(.x, -1)))
  
  char_final <- ""
  
  for(ID_experiment_i in unique(df_DEB_mean$ID_experiment)){
    
    # Data corresponding to this specific cosm
    df_DEB_i <- subset(df_DEB_mean, ID_experiment == ID_experiment_i)
    
    Time_sim_i <- max(df_DEB_i$Time)
    
    # Soil OM data
    df_OM_soil_i <- df_DEB_i |> 
      arrange(Time) |>  
      dplyr::select(Time, OM_soil, Environment) |> 
      filter(Environment == "Soil_change") 
    
    char_Replace_tOM_soil <- df_OM_soil_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_Replace_OM_soil  <- df_OM_soil_i |> pull(OM_soil) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Replace_soil <- paste(rep("Replace, ", length(df_OM_soil_i$Time)), collapse = "")
    
    df_OM_horse_Replace_i <- df_DEB_i |> 
      arrange(Time) |>  
      dplyr::select(Time, OM_horse, Environment) |> 
      filter(Environment == "Soil_change") 
    
    char_Replace_tOM_horse <- df_OM_horse_Replace_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_Replace_OM_horse  <- df_OM_horse_Replace_i |> pull(OM_horse) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Replace_horse <- paste(rep("Replace, ", length(df_OM_horse_Replace_i$Time)), collapse = "")
    
    df_OM_horse_Add_i <- df_DEB_i |> 
      arrange(Time) |>  
      dplyr::select(Time, OM_horse, Environment) |> 
      filter(Environment == "Food") 
    
    char_Add_tOM_horse <- df_OM_horse_Add_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_Add_OM_horse  <- df_OM_horse_Add_i |> pull(OM_horse) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Add_horse <- paste(rep("Add, ", length(df_OM_horse_Add_i$Time)), collapse = "")
    
    # Density data
    df_density_i <- df_DEB_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Density != lag(Density))
    
    char_tdensity <- df_density_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_density  <- df_density_i |> pull(Density) |>  paste(collapse = ", ") # Nb individual in a cosm
    
    # Soil quantity in cosm (g)
    df_Soilw_i <- df_DEB_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Soil_w != lag(Soil_w))
    
    char_tsoilw <- df_Soilw_i |> pull(Time)   |>  paste(collapse = ", ") # Time in days
    char_soilw  <- df_Soilw_i |> pull(Soil_w) |>  paste(collapse = ", ") # Nb individual in a cosm
    
    # MCSim syntaxe
    
    text_init <- paste("    Init= 1;", sep="")
    
    text_Winit <- paste("    Winit=", subset(df_DEB_i, Time==0)$Weight, ";", sep="")
    
    text_Soil <- paste("    Soil=", subset(df_DEB_i, Time==0)$Soil_type, ";", sep="")
    
    text_OM_soil_init <- paste("    OM_soil_t0=", subset(df_DEB_i, Time==0)$OM_soil, ";", sep="")
    text_OM_horse_init <- paste("    OM_horse_t0=", subset(df_DEB_i, Time==0)$OM_horse, ";", sep="")
    
    
    text_OM_soil_Replace <- paste("    OM_soil_replace=Events(OM_soil,\n                          ", 
                                  length(df_OM_soil_i$Time), ",\n                          ", 
                                  char_Replace_tOM_soil,",\n                          ", 
                                  char_Replace_soil, "\n                           ",
                                  char_Replace_OM_soil, ");", sep="")
    
    
    
    text_OM_horse_Replace <- paste("    OM_horse_replace=Events(OM_horse,\n                           ", 
                                   length(df_OM_horse_Replace_i$Time), ",\n                           ", 
                                   char_Replace_tOM_horse,",\n                           ", 
                                   char_Replace_horse, "\n                            ",
                                   char_Replace_OM_horse, ");", sep="")
    
    text_OM_horse_Add <- paste("    OM_horse_add=Events(OM_horse,\n                           ", 
                               length(df_OM_horse_Add_i$Time), ",\n                           ", 
                               char_Add_tOM_horse,",\n                           ", 
                               char_Add_horse, "\n                            ",
                               char_Add_OM_horse, ");
                               ", sep="")
    
    
    
    if(OM_diff){
      text_OM_final <- paste(
        text_OM_soil_Replace,
        text_OM_horse_Replace,
        text_OM_horse_Add,
        sep="\n"
      )
    }else {
      text_OM_final <- text_OM_tot
    }
    
    
    
    if(length(df_density_i$Density) == 1){
      text_density <- paste("    Dens=",char_density, ";", sep="")
    }else{
      text_density <- paste("    Dens=NDoses(", 
                            length(df_density_i$Time), ",\n                ", 
                            char_density,",\n                ", 
                            char_tdensity, ");", sep="")
    }
    
    if (length(df_Soilw_i$Time)==1) {
      text_WeightSoilCosm <- paste0("    WeightSoilCosm = ", df_Soilw_i$Soil_w,";")
    }else {
      text_WeightSoilCosm <- paste("    WeightSoilCosm=NDoses(", 
                                   length(df_Soilw_i$Time), ",\n                ", 
                                   char_soilw,",\n                ", 
                                   char_tsoilw, ");", sep="")
    }
    
    text_Texp <- paste("    Texp=",subset(df_DEB_i, Time==0)$Texp, ";", sep="")
    
    text_print <- paste0("    PrintStep(", Endpoints_print, ",
				0, ", Time_sim_i, ", 0.1);")
    
    
    char_i <- paste(
      paste("Simulation { #", ID_experiment_i), 
      
      text_init,
      text_WeightSoilCosm,
      text_Winit,
      text_OM_soil_init,
      text_OM_horse_init,
      text_Soil,
      text_Texp,
      text_density,
      text_OM_final,
      text_print,
      
      paste("}"),
      sep="\n"
      
    )
    
    char_final <- paste(char_final, char_i, sep="\n")
    
  }
  
  return(char_final)
}

f_create_MC_block <- function(name, seed, Nb_Iter, RTOL, ATOL) {
  
  sprintf(
    'MonteCarlo( "%s.out",  # output file
      %d,           # iterations
      %d);          # random seed

Integrate (Lsodes,  %g, %g, 0); # Integrate(Solver, RTOL, ATOL, ITOL);
    ',
    name, Nb_Iter, seed, RTOL, ATOL
  )
}

f_In_SimpleSimulation <- function(file_path, RTOL, ATOL, text_param, l_Experiments, Endpoints_print, Adults_alone, OM_diff) {
  
  text_Simulations <- f_In_simulation_MC(l_Experiments, Endpoints_print, Adults_alone, OM_diff)
  
  text_start <- sprintf(
    '#### DEB A. caliginosa
#===============================================

Integrate (Lsodes,  %g, %g, 0); # Integrate(Solver, RTOL, ATOL, ITOL);

########## Experiments ################################################

    ',
    RTOL, ATOL
  )
  
  text_end <- "
  End."
  
  text_full <- paste(
    text_start,
    text_param,
    text_Simulations,
    text_end,
    sep = "\n"
  )
  
  writeLines(
    text_full,
    here::here(file_path, "DEB_SimpleSim.in")
  )
}

f_In_MC_tot <- function(File_path, seed, Nb_Iter, RTOL, ATOL, Endpoints_print, text_distrib, l_Experiments, Adults_alone, OM_diff) {
  
  for (i in 1:length(l_Experiments)) {
    
    text_mc_block_i <- f_create_MC_block(paste0("DEB_MC_", i), seed, Nb_Iter, RTOL, ATOL)
    
    expe_i <- l_Experiments[i]
    text_Simulation_i <- f_In_simulation_MC(expe_i, Endpoints_print, Adults_alone, OM_diff)
    
    text_end <- "
End."
    
    text_full_i <- paste0(text_mc_block_i, text_distrib, text_Simulation_i, text_end, sep="\n")
    
    writeLines(
      text_full_i,
      here::here(File_path, paste0("DEB_MC_", i, ".in"))
    )
  }
  
}

f_read_MC <- function(path_file, nb_param, times_sim){
  
  df_data <- data.table::fread(
    path_file,
    header = TRUE,
    showProgress = TRUE,
    nThread = 2
  )
  
  df_data_param <- df_data[,1:(nb_param+1)] |>
    pivot_longer(
      cols = -Iter,              # garde Iter comme identifiant
      names_to = "Parameter",    # noms des paramètres
      values_to = "Value"        # valeurs associées
    ) |>
    mutate(Type = "MC")
  
  state_cols <- names(df_data)[grepl("_\\d+\\.\\d+$", names(df_data))]
  param_cols <- setdiff(names(df_data), c("Iter", state_cols))
  
  df_data_long <- df_data %>%
    pivot_longer(
      cols = all_of(state_cols),
      names_to = c("Variable", "No", "Time_idx"),
      names_pattern = "([A-Za-z_]+)_(\\d+)\\.(\\d+)",
      values_to = "Value"
    ) %>%
    mutate(
      No = as.integer(No),
      Time_idx = as.integer(Time_idx)
    ) %>%
    rowwise() %>%
    mutate(
      Time = times_sim[[No]][Time_idx]
    ) %>%
    ungroup() %>%
    # Garder les colonnes dans l'ordre voulu
    dplyr::select(Iter, No, Time_idx, Time, all_of(param_cols), everything())
  
  
  res <- list(
    df_MC_tot = df_data_long,
    df_MC = df_data_param
  )
  
  return(res)
}

f_read_all_MC <- function(path_file, l_Experiments, nb_param, times_sim, times_keep, time_limits){
  
  df_res <- data.frame()
  
  times_sim <- list(
    seq(0, time_limits$time_cutoff[1], 0.1),
    seq(0, time_limits$time_cutoff[2], 0.1),
    seq(0, time_limits$time_cutoff[3], 0.1),
    seq(0, time_limits$time_cutoff[4], 0.1),
    seq(0, time_limits$time_cutoff[5], 0.1),
    seq(0, time_limits$time_cutoff[6], 0.1),
    seq(0, time_limits$time_cutoff[7], 0.1),
    seq(0, time_limits$time_cutoff[8], 0.1),
    seq(0, time_limits$time_cutoff[9], 0.1),
    seq(0, time_limits$time_cutoff[10], 0.1),
    seq(0, time_limits$time_cutoff[11], 0.1),
    seq(0, time_limits$time_cutoff[12], 0.1)
  )
  
  for (i in 1:length(l_Experiments)) {
    
    print(l_Experiments[i])
    path_MC_i <- paste0(path_file, "/DEB_MC_", i, ".out")
    times_sim_i <- list()
    times_sim_i[[1]] <- times_sim[[i]]
    df_MC_i <- f_read_MC(path_MC_i, nb_param, times_sim_i)$df_MC_tot|> 
      filter(Time %in% times_keep) |> 
      mutate(No_experiment = i)  |> 
      left_join(time_limits, by = "No_experiment") |> 
      filter(Time <= time_cutoff) |> 
      dplyr::select(-time_cutoff)
    
    df_res <- rbind(df_res, df_MC_i)
  }
  return(df_res)
}

f_MC_RMSE_calc <- function(df_data, df_MonteCarlo_sim){
  
  df_data_simple <- df_data |> 
    dplyr::select(No_sim, Time, Reproduction, Weight) |> 
    mutate(
      No_experiment = as.numeric(No_sim),
      Time = as.numeric(Time)
    )
  
  df_MonteCarlo_sim_wide <- df_MonteCarlo_sim |> 
    distinct() |> 
    filter(Variable %in% c("Weight", "Reproduction")) |> 
    pivot_wider(
      id_cols = c(No_experiment, Time, Iter),
      names_from = Variable, 
      values_from = Value
    ) |> 
    mutate(
      No_experiment = as.numeric(No_experiment)
    )
  
  Sd_Weight <- sd(df_data_simple$Weight, na.rm=TRUE)
  Sd_Reproduction <- sd(df_data_simple$Reproduction, na.rm=TRUE)
  
  df_error <- df_MonteCarlo_sim_wide %>%
    inner_join(
      df_data_simple,
      by = c("No_experiment", "Time")
    ) %>%
    mutate(
      err_weight = (Weight.x - Weight.y)^2,
      err_repro  = (Reproduction.x - Reproduction.y)^2
    ) %>%
    summarise(
      RMSE_Weight = sqrt(mean(err_weight, na.rm = TRUE)),
      RMSE_Reproduction = sqrt(mean(err_repro, na.rm = TRUE)),
      .by = Iter
    ) |> 
    mutate(
      NRMSE_Weight = RMSE_Weight/Sd_Weight,
      NRMSE_Reproduction = RMSE_Reproduction/Sd_Reproduction,
      NRMSE_tot = sqrt((NRMSE_Weight^2+NRMSE_Reproduction^2)/2)
    )
  
  return(df_error)
}

# MCMC ----

f_In_experiments <- function(l_Experiments, Adult_alone, OM_diff){
  
  df_DEB_mean <- f_import_data_DEB(l_Experiments, Adults_alone) |> 
    mutate(
      OM_soil = signif(OM_soil, 3),
      OM_horse = signif(OM_horse, 3)
    )  |> 
    mutate(across(where(is.numeric), ~ replace_na(.x, -1)))
  
  char_final <- ""
  
  XP_no_R_data <- c("Bart2019_XP1_1", "Bart2019_XP4_1", "Bart2019_XP3", "Bart2020", "Gollot2026_D1")
  
  for(ID_experiment_i in unique(df_DEB_mean$ID_experiment)){
    
    
    bol_XP3 <- 0
    if (ID_experiment_i == "Bart2019_XP3"){
      bol_XP3 <- 1
    }
    
    # Data corresponding to this specific cosm
    df_DEB_i <- subset(df_DEB_mean, ID_experiment == ID_experiment_i)
    
    # Weight data (remove Weigth(t=0))
    df_Ww_i <- df_DEB_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Weight != lag(Weight))
    
    char_tW <- df_Ww_i |> pull(Time)   |>  paste(collapse = ", ") # Time in days
    char_Ww <- df_Ww_i |> pull(Weight) |>  paste(collapse = ", ") # Weight (g)
    
    # Reproduction data
    df_R_i <- df_DEB_i |> 
       arrange(Time) 
    
    Time_xp_i <- max(df_DEB_i$Time)
    
    if (ID_experiment_i %in% XP_no_R_data) {
      char_tR <- paste0("0,", Time_xp_i)
      char_R <- "0, -1"
    } else {
      char_tR <- df_R_i |> pull(Time)         |>  paste(collapse = ", ") # Time in days
      char_R  <- df_R_i |> pull(Reproduction) |>  paste(collapse = ", ") # Cumulated reproduction (#)
    }
    
    # Soil OM data
    df_OM_soil_i <- df_DEB_i |> 
      arrange(Time) |>  
      dplyr::select(Time, OM_soil, Environment) |> 
      filter(Environment == "Soil_change") 
    
    char_Replace_tOM_soil <- df_OM_soil_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_Replace_OM_soil  <- df_OM_soil_i |> pull(OM_soil) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Replace_soil <- paste(rep("Replace, ", length(df_OM_soil_i$Time)), collapse = "")
    
    # Horse OM data
    df_OM_horse_Replace_i <- df_DEB_i |> 
      arrange(Time) |>  
      dplyr::select(Time, OM_horse, Environment) |> 
      filter(Environment == "Soil_change") 
    
    char_Replace_tOM_horse <- df_OM_horse_Replace_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_Replace_OM_horse  <- df_OM_horse_Replace_i |> pull(OM_horse) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Replace_horse <- paste(rep("Replace, ", length(df_OM_horse_Replace_i$Time)), collapse = "")
    
    df_OM_horse_Add_i <- df_DEB_i |> 
      arrange(Time) |>  
      dplyr::select(Time, OM_horse, Environment) |> 
      filter(Environment == "Food") 
    
    char_Add_tOM_horse <- df_OM_horse_Add_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_Add_OM_horse  <- df_OM_horse_Add_i |> pull(OM_horse) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Add_horse <- paste(rep("Add, ", length(df_OM_horse_Add_i$Time)), collapse = "")
    
    # OM tot data
    df_OM_tot_i <- df_DEB_i |> 
      arrange(Time) |>                    
      dplyr::select(Time, OM_horse, OM_soil, Environment) |> 
      filter(Environment %in% c("Soil_change", "Food")) |> 
      mutate(
        OM_soil = case_when(
          Environment == "Soil_change" ~ OM_soil,
          Environment == "Food" ~ 0
        ),
        OM_tot = OM_soil+OM_horse,
        Word = case_when(
          Environment == "Soil_change" ~ "Replace",
          Environment == "Food" ~ "Add"
        )
      )
    
    char_tOM_tot<- df_OM_tot_i  |> pull(Time)     |>  paste(collapse = ", ") # Time in days
    char_OM_tot  <- df_OM_tot_i |> pull(OM_tot) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Add_tot <- df_OM_tot_i |> pull(Word)     |>  paste(collapse = ", ")
    
    # Density data
    df_density_i <- df_DEB_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Density != lag(Density))
    
    char_tdensity <- df_density_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_density  <- df_density_i |> pull(Density) |>  paste(collapse = ", ") # Nb individual in a cosm
    
    # Soil quantity in cosm (g)
    df_Soilw_i <- df_DEB_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Soil_w != lag(Soil_w))
    
    char_tsoilw <- df_Soilw_i |> pull(Time)   |>  paste(collapse = ", ") # Time in days
    char_soilw  <- df_Soilw_i |> pull(Soil_w) |>  paste(collapse = ", ") # Nb individual in a cosm
    
    # MCSim syntaxe
    
    text_init <- paste("    Init= 1;", sep="")
    
    text_Winit <- paste("    Winit=", subset(df_DEB_i, Time==0)$Weight, ";", sep="")
    
    text_OM_soil_init <- paste("    OM_soil_t0=", subset(df_DEB_i, Time==0)$OM_soil, ";", sep="")
    text_OM_horse_init <- paste("    OM_horse_t0=", subset(df_DEB_i, Time==0)$OM_horse, ";", sep="")
    
    text_Ww <- paste(
      paste("    Print(Weight,", char_tW, ");" , sep=""),
      paste("    Data(Weight,", char_Ww,  ");" , sep=""), 
      sep = "\n"
    )
    text_Soil <- paste("    Soil=", subset(df_DEB_i, Time==0)$Soil_type, ";", sep="")
    
    text_R <- paste(
      paste("    Print(Reproduction,", char_tR, ");" , sep=""),
      paste("    Data(Reproduction,", char_R,  ");" , sep=""), 
      sep = "\n")
    
    
    text_OM_soil_Replace <- paste("    OM_soil_replace=Events(OM_soil,\n                          ", 
                                  length(df_OM_soil_i$Time), ",\n                          ", 
                                  char_Replace_tOM_soil,",\n                          ", 
                                  char_Replace_soil, "\n                           ",
                                  char_Replace_OM_soil, ");", sep="")
    
    
    
    text_OM_horse_Replace <- paste("    OM_horse_replace=Events(OM_horse,\n                           ", 
                                   length(df_OM_horse_Replace_i$Time), ",\n                           ", 
                                   char_Replace_tOM_horse,",\n                           ", 
                                   char_Replace_horse, "\n                            ",
                                   char_Replace_OM_horse, ");", sep="")
    
    text_OM_horse_Add <- paste("    OM_horse_add=Events(OM_horse,\n                           ", 
                               length(df_OM_horse_Add_i$Time), ",\n                           ", 
                               char_Add_tOM_horse,",\n                           ", 
                               char_Add_horse, "\n                            ",
                               char_Add_OM_horse, ");", sep="")
    
    text_OM_tot <- paste("    event_OM_tot=Events(OM_tot,\n                           ", 
                         length(df_OM_tot_i$Time), ",\n                           ", 
                         char_tOM_tot,",\n                           ", 
                         char_Add_tot, ",\n                            ",
                         char_OM_tot, ");", sep="")
    text_bol_XP3 <- paste("    bol_XP3=", bol_XP3, ";")
    
    if(OM_diff){
      text_OM_final <- paste(
        text_OM_soil_Replace,
        text_OM_horse_Replace,
        text_OM_horse_Add,
       # text_bol_XP3,
        sep="\n"
      )
    }else {
      text_OM_final <- text_OM_tot
    }
    
    if(length(df_density_i$Density) == 1){
      text_density <- paste("    Dens=",char_density, ";", sep="")
    }else{
      text_density <- paste("    Dens=NDoses(", 
                            length(df_density_i$Time), ",\n                ", 
                            char_density,",\n                ", 
                            char_tdensity, ");", sep="")
    }
    
    if (length(df_Soilw_i$Time)==1) {
      text_WeightSoilCosm <- paste0("    WeightSoilCosm = ", df_Soilw_i$Soil_w,";")
    }else {
      text_WeightSoilCosm <- paste("    WeightSoilCosm=NDoses(", 
                                   length(df_Soilw_i$Time), ",\n                ", 
                                   char_soilw,",\n                ", 
                                   char_tsoilw, ");", sep="")
    }
    
    text_Texp <- paste("    Texp=",subset(df_DEB_i, Time==0)$Texp, ";", sep="")
    
    
    char_i <- paste(
      paste("Experiment { #", ID_experiment_i), 
      
      text_init,
      text_WeightSoilCosm,
      text_Winit,
      text_OM_soil_init,
      text_OM_horse_init,
      text_Soil,
      text_Texp,
      text_density,
      text_OM_final,
      text_Ww,
      text_R,
      
      paste("}"),
      sep="\n"
      
    )
    
    char_final <- paste(char_final, char_i, sep="\n")
    
  }
  
  return(char_final)
}


f_In_experiments_sim <- function(l_Experiments, Adult_alone, OM_diff){
  
  df_DEB_mean <- read.csv(here::here("data/DEB_sim_data.csv")) |> 
    mutate(
      OM_soil = signif(OM_soil, 3),
      OM_horse = signif(OM_horse, 3)
    )  |> 
    mutate(across(where(is.numeric), ~ replace_na(.x, -1)))
  
  char_final <- ""
  
  for(ID_experiment_i in unique(df_DEB_mean$ID_experiment)){
    
    # Data corresponding to this specific cosm
    df_DEB_i <- subset(df_DEB_mean, ID_experiment == ID_experiment_i)
    
    # Weight data (remove Weigth(t=0))
    df_Ww_i <- df_DEB_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Weight != lag(Weight))
    
    char_tW <- df_Ww_i |> pull(Time)   |>  paste(collapse = ", ") # Time in days
    char_Ww <- df_Ww_i |> pull(Weight) |>  paste(collapse = ", ") # Weight (g)
    
    # Reproduction data
    df_R_i <- df_DEB_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Reproduction != lag(Reproduction))
    
    char_tR <- df_R_i |> pull(Time)         |>  paste(collapse = ", ") # Time in days
    char_R  <- df_R_i |> pull(Reproduction) |>  paste(collapse = ", ") # Cumulated reproduction (#)
    
    # Soil OM data
    df_OM_soil_i <- df_DEB_i |> 
      arrange(Time) |>  
      dplyr::select(Time, OM_soil, Environment) |> 
      filter(Environment == "Soil_change") 
    
    char_tOM_soil <- df_OM_soil_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_OM_soil  <- df_OM_soil_i |> pull(OM_soil) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Add_soil <- paste(rep("Replace, ", length(df_OM_soil_i$Time)), collapse = "")
    
    # Horse dung OM data
    df_OM_horse_i <- df_DEB_i |> 
      arrange(Time) |>                    
      dplyr::select(Time, OM_horse, Environment) |> 
      filter(Environment %in% c("Soil_change", "Food")) |> 
      mutate(
        Word = case_when(
          Environment == "Soil_change" ~ "Replace",
          Environment == "Food" ~ "Add"
        )
      )
    
    char_tOM_horse<- df_OM_horse_i  |> pull(Time)     |>  paste(collapse = ", ") # Time in days
    char_OM_horse  <- df_OM_horse_i |> pull(OM_horse) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Add_horse <- df_OM_horse_i |> pull(Word)     |>  paste(collapse = ", ")
    
    # OM tot data
    df_OM_tot_i <- df_DEB_i |> 
      arrange(Time) |>                    
      dplyr::select(Time, OM_horse, OM_soil, Environment) |> 
      filter(Environment %in% c("Soil_change", "Food")) |> 
      mutate(
        OM_soil = case_when(
          Environment == "Soil_change" ~ OM_soil,
          Environment == "Food" ~ 0
        ),
        OM_tot = OM_soil+OM_horse,
        Word = case_when(
          Environment == "Soil_change" ~ "Replace",
          Environment == "Food" ~ "Add"
        )
      )
    
    char_tOM_tot<- df_OM_tot_i  |> pull(Time)     |>  paste(collapse = ", ") # Time in days
    char_OM_tot  <- df_OM_tot_i |> pull(OM_tot) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Add_tot <- df_OM_tot_i |> pull(Word)     |>  paste(collapse = ", ")
    
    # Density data
    df_density_i <- df_DEB_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Density != lag(Density))
    
    char_tdensity <- df_density_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_density  <- df_density_i |> pull(Density) |>  paste(collapse = ", ") # Nb individual in a cosm
    
    # Soil quantity in cosm (g)
    df_Soilw_i <- df_DEB_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Soil_w != lag(Soil_w))
    
    char_tsoilw <- df_Soilw_i |> pull(Time)   |>  paste(collapse = ", ") # Time in days
    char_soilw  <- df_Soilw_i |> pull(Soil_w) |>  paste(collapse = ", ") # Nb individual in a cosm
    
    # MCSim syntaxe
    
    text_init <- paste("    Init= 1;", sep="")
    
    text_Winit <- paste("    Winit=", subset(df_DEB_i, Time==0)$Weight, ";", sep="")
    text_Ww <- paste(
      paste("    Print(Weight,", char_tW, ");" , sep=""),
      paste("    Data(Weight,", char_Ww,  ");" , sep=""), 
      sep = "\n"
    )
    text_Soil <- paste("    Soil=", subset(df_DEB_i, Time==0)$Soil_type, ";", sep="")
    
    text_R <- paste(
      paste("    Print(Reproduction,", char_tR, ");" , sep=""),
      paste("    Data(Reproduction,", char_R,  ");" , sep=""), 
      sep = "\n")
    
    
    text_OM_soil <- paste("    event_OM_soil=Events(OM_soil,\n                          ", 
                          length(df_OM_soil_i$Time), ",\n                          ", 
                          char_tOM_soil,",\n                          ", 
                          char_Add_soil, "\n                           ",
                          char_OM_soil, ");", sep="")
    
    text_OM_horse <- paste("    event_OM_horse=Events(OM_horse,\n                           ", 
                           length(df_OM_horse_i$Time), ",\n                           ", 
                           char_tOM_horse,",\n                           ", 
                           char_Add_horse, ",\n                            ",
                           char_OM_horse, ");", sep="")
    
    text_OM_tot <- paste("    event_OM_tot=Events(OM_tot,\n                           ", 
                         length(df_OM_tot_i$Time), ",\n                           ", 
                         char_tOM_tot,",\n                           ", 
                         char_Add_tot, ",\n                            ",
                         char_OM_tot, ");", sep="")
    
    if(OM_diff){
      text_OM_final <- paste(
        text_OM_soil,
        text_OM_horse,
        sep="\n"
      )
    }else {
      text_OM_final <- text_OM_tot
    }
    
    if(length(df_density_i$Density) == 1){
      text_density <- paste("    Dens=",char_density, ";", sep="")
    }else{
      text_density <- paste("    Dens=NDoses(", 
                            length(df_density_i$Time), ",\n                ", 
                            char_density,",\n                ", 
                            char_tdensity, ");", sep="")
    }
    
    if (length(df_Soilw_i$Time)==1) {
      text_WeightSoilCosm <- paste0("    WeightSoilCosm = ", df_Soilw_i$Soil_w,";")
    }else {
      text_WeightSoilCosm <- paste("    WeightSoilCosm=NDoses(", 
                                   length(df_Soilw_i$Time), ",\n                ", 
                                   char_soilw,",\n                ", 
                                   char_tsoilw, ");", sep="")
    }
    
    text_Texp <- paste("    Texp=",subset(df_DEB_i, Time==0)$Texp, ";", sep="")
    
    
    char_i <- paste(
      paste("Experiment { #", ID_experiment_i), 
      
      text_init,
      text_WeightSoilCosm,
      text_Winit,
      text_Soil,
      text_Texp,
      text_OM_final,
      text_density,
      text_Ww,
      text_R,
      
      paste("}"),
      sep="\n"
      
    )
    
    char_final <- paste(char_final, char_i, sep="\n")
    
  }
  
  return(char_final)
}

f_create_mcmc_block <- function(name, seeds, Nb_Iter, RTOL, ATOL) {
  
  sprintf(
    'MCMC( "%s.out",  # output file
      "",                # name of restart file
      "",                # name of data file
      %d, 0,             # iterations, print prediction flag
      1, 10000,          # printing frequency, iters to print
      %d);               # random seed

Integrate (Lsodes,  %g, %g, 0); # Integrate(Solver, RTOL, ATOL, ITOL);
    ',
    name, Nb_Iter, seeds, RTOL, ATOL
  )
}

f_In_tot <- function(file_path, l_Experiments, Adults_alone = FALSE, OM_diff, text_priors, text_likelihood, NbIter, seeds, RTOL, ATOL) {
  
  text_Level_global <- 
    'Level{ # Global
    '
  
  text_Level_exp <- '
Level{

  ############## Experiments ###################
  '
  
  text_experiment <- f_In_experiments(l_Experiments, Adult_alone, OM_diff)
  
  text_end <- '
} # End
} # End global

End.'
  
  # .in generation and saving
  for (id in names(seeds)) {
    text_start <- f_create_mcmc_block(paste0("DEB_", id), seeds[[id]], NbIter, RTOL, ATOL)
    
    text_full <- paste(
      text_start,
      text_Level_global,
      text_priors,
      text_likelihood,
      text_Level_exp,
      text_experiment,
      text_end,
      sep = "\n"
    )
    
    writeLines(
      text_full,
      here::here(file_path, paste0("DEB_", id, ".in"))
    )
  }
}

f_In_tot_sim <- function(file_path, l_Experiments, Adults_alone = FALSE, OM_diff, text_priors, text_likelihood, NbIter, seeds) {
  
  text_Level_global <- 
    'Level{ # Global
    '
  
  text_Level_exp <- '
Level{

  ############## Experiments ###################
  '
  
  text_experiment <- f_In_experiments_sim(l_Experiments, Adult_alone, OM_diff)
  
  text_end <- '
} # End
} # End global

End.'
  
  # .in generation and saving
  for (id in names(seeds)) {
    text_start <- f_create_mcmc_block(paste0("DEB_", id), seeds[[id]], NbIter)
    
    text_full <- paste(
      text_start,
      text_Level_global,
      text_priors,
      text_likelihood,
      text_Level_exp,
      text_experiment,
      text_end,
      sep = "\n"
    )
    
    writeLines(
      text_full,
      here::here(file_path, paste0("DEB_", id, ".in"))
    )
  }
}

f_In_simulations <- function(file_path, l_Experiments, Adults_alone, OM_diff, text_param, Endpoints_print, step_sim){

  
  df_DEB_mean <- f_import_data_DEB(l_Experiments, Adults_alone) |> 
    mutate(
      OM_soil = signif(OM_soil, 3),
      OM_horse = signif(OM_horse, 3)
    )  |> 
    mutate(across(where(is.numeric), ~ replace_na(.x, -1)))
  
  char_experiments <- ""
  
  for(ID_experiment_i in unique(df_DEB_mean$ID_experiment)){
    
    # Data corresponding to this specific cosm
    df_DEB_i <- subset(df_DEB_mean, ID_experiment == ID_experiment_i)
    
    Time_sim_i <- max(df_DEB_i$Time)
    
    # Soil OM data
    df_OM_soil_i <- df_DEB_i |> 
      arrange(Time) |>  
      dplyr::select(Time, OM_soil, Environment) |> 
      filter(Environment == "Soil_change") 
    
    char_Replace_tOM_soil <- df_OM_soil_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_Replace_OM_soil  <- df_OM_soil_i |> pull(OM_soil) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Replace_soil <- paste(rep("Replace, ", length(df_OM_soil_i$Time)), collapse = "")
    
    df_OM_horse_Replace_i <- df_DEB_i |> 
      arrange(Time) |>  
      dplyr::select(Time, OM_horse, Environment) |> 
      filter(Environment == "Soil_change") 
    
    char_Replace_tOM_horse <- df_OM_horse_Replace_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_Replace_OM_horse  <- df_OM_horse_Replace_i |> pull(OM_horse) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Replace_horse <- paste(rep("Replace, ", length(df_OM_horse_Replace_i$Time)), collapse = "")
    
    df_OM_horse_Add_i <- df_DEB_i |> 
      arrange(Time) |>  
      dplyr::select(Time, OM_horse, Environment) |> 
      filter(Environment == "Food") 
    
    char_Add_tOM_horse <- df_OM_horse_Add_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_Add_OM_horse  <- df_OM_horse_Add_i |> pull(OM_horse) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Add_horse <- paste(rep("Add, ", length(df_OM_horse_Add_i$Time)), collapse = "")
    
    # Density data
    df_density_i <- df_DEB_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Density != lag(Density))
    
    char_tdensity <- df_density_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_density  <- df_density_i |> pull(Density) |>  paste(collapse = ", ") # Nb individual in a cosm
    
    # Soil quantity in cosm (g)
    df_Soilw_i <- df_DEB_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Soil_w != lag(Soil_w))
    
    char_tsoilw <- df_Soilw_i |> pull(Time)   |>  paste(collapse = ", ") # Time in days
    char_soilw  <- df_Soilw_i |> pull(Soil_w) |>  paste(collapse = ", ") # Nb individual in a cosm
    
    # MCSim syntaxe
    
    text_init <- paste("    Init= 1;", sep="")
    
    text_Winit <- paste("    Winit=", subset(df_DEB_i, Time==0)$Weight, ";", sep="")
    
    text_Soil <- paste("    Soil=", subset(df_DEB_i, Time==0)$Soil_type, ";", sep="")
    
    text_OM_soil_init <- paste("    OM_soil_t0=", subset(df_DEB_i, Time==0)$OM_soil, ";", sep="")
    text_OM_horse_init <- paste("    OM_horse_t0=", subset(df_DEB_i, Time==0)$OM_horse, ";", sep="")
    
    
    text_OM_soil_Replace <- paste("    OM_soil_replace=Events(OM_soil,\n                          ", 
                                  length(df_OM_soil_i$Time), ",\n                          ", 
                                  char_Replace_tOM_soil,",\n                          ", 
                                  char_Replace_soil, "\n                           ",
                                  char_Replace_OM_soil, ");", sep="")
    
    
    
    text_OM_horse_Replace <- paste("    OM_horse_replace=Events(OM_horse,\n                           ", 
                                   length(df_OM_horse_Replace_i$Time), ",\n                           ", 
                                   char_Replace_tOM_horse,",\n                           ", 
                                   char_Replace_horse, "\n                            ",
                                   char_Replace_OM_horse, ");", sep="")
    
    text_OM_horse_Add <- paste("    OM_horse_add=Events(OM_horse,\n                           ", 
                               length(df_OM_horse_Add_i$Time), ",\n                           ", 
                               char_Add_tOM_horse,",\n                           ", 
                               char_Add_horse, "\n                            ",
                               char_Add_OM_horse, ");
                               ", sep="")
    
    
    
    if(OM_diff){
      text_OM_final <- paste(
        text_OM_soil_Replace,
        text_OM_horse_Replace,
        text_OM_horse_Add,
        sep="\n"
      )
    }else {
      text_OM_final <- text_OM_tot
    }
    
    
    
    if(length(df_density_i$Density) == 1){
      text_density <- paste("    Dens=",char_density, ";", sep="")
    }else{
      text_density <- paste("    Dens=NDoses(", 
                            length(df_density_i$Time), ",\n                ", 
                            char_density,",\n                ", 
                            char_tdensity, ");", sep="")
    }
    
    if (length(df_Soilw_i$Time)==1) {
      text_WeightSoilCosm <- paste0("    WeightSoilCosm = ", df_Soilw_i$Soil_w,";")
    }else {
      text_WeightSoilCosm <- paste("    WeightSoilCosm=NDoses(", 
                                   length(df_Soilw_i$Time), ",\n                ", 
                                   char_soilw,",\n                ", 
                                   char_tsoilw, ");", sep="")
    }
    
    text_Texp <- paste("    Texp=",subset(df_DEB_i, Time==0)$Texp, ";", sep="")
    
    text_print <- paste0("    PrintStep(", Endpoints_print, ",
				0, ", Time_sim_i, ", ", step_sim,");")
    
    
    char_i <- paste(
      paste("Simulation { #", ID_experiment_i), 
      
      text_init,
      text_WeightSoilCosm,
      text_Winit,
      text_OM_soil_init,
      text_OM_horse_init,
      text_Soil,
      text_Texp,
      text_density,
      text_OM_final,
      text_print,
      
      paste("}"),
      sep="\n"
      
    )
    
    char_experiments <- paste(char_experiments, char_i, sep="\n")
    
  }
  
  text_start <- "#### DEB A. caliginosa
#===============================================

Integrate(Lsodes, 1E-5, 1E-6, 1);

########## Experiments ################################################
"
  
  text_end <- "
  End."
  
  text_full <- paste(
    text_start,
    text_param,
    char_experiments,
    text_end,
    sep = "\n"
  )
  
  writeLines(
    text_full,
    here::here(file_path, "DEB_sim_full.in")
  )
}

f_In_predobs <- function(file_path, l_Experiments, Adults_alone, OM_diff, l_param_name){
  
  XP_no_R_data <- c("Bart2019_XP1_1", "Bart2019_XP4_1", "Bart2019_XP3", "Bart2020", "Gollot2026_D1")
  
  df_DEB_mean <- f_import_data_DEB(l_Experiments, Adults_alone, fit = TRUE) |> 
    mutate(
      OM_soil = signif(OM_soil, 3),
      OM_horse = signif(OM_horse, 3)
    )  |> 
    mutate(across(where(is.numeric), ~ replace_na(.x, -1)))
  
  compteur_exp <- 0 
  
  for(ID_experiment_i in unique(df_DEB_mean$ID_experiment)){
    text_full_i <- ""
    char_i <- ""
    compteur_exp <- compteur_exp + 1
    
    text_start <- paste0('#### DEB A. caliginosa
#===============================================

SetPoints("DEB_setpoint_predobs_', compteur_exp, '.out", "tab_setpoint.out", 0,', l_param_name|>  paste(collapse = ", "),');

Integrate(Lsodes, 1E-5, 1E-6, 1);

########## Experiments ################################################
')
    
    text_end <- "
  End."
    
    # Data corresponding to this experiment
    df_DEB_i <- subset(df_DEB_mean, ID_experiment == ID_experiment_i)
    
    Time_sim_i <- max(df_DEB_i$Time)
    
    # Soil OM data
    df_OM_soil_i <- df_DEB_i |> 
      arrange(Time) |>  
      dplyr::select(Time, OM_soil, Environment) |> 
      filter(Environment == "Soil_change") 
    
    char_Replace_tOM_soil <- df_OM_soil_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_Replace_OM_soil  <- df_OM_soil_i |> pull(OM_soil) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Replace_soil <- paste(rep("Replace, ", length(df_OM_soil_i$Time)), collapse = "")
    
    df_OM_horse_Replace_i <- df_DEB_i |> 
      arrange(Time) |>  
      dplyr::select(Time, OM_horse, Environment) |> 
      filter(Environment == "Soil_change") 
    
    char_Replace_tOM_horse <- df_OM_horse_Replace_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_Replace_OM_horse  <- df_OM_horse_Replace_i |> pull(OM_horse) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Replace_horse <- paste(rep("Replace, ", length(df_OM_horse_Replace_i$Time)), collapse = "")
    
    df_OM_horse_Add_i <- df_DEB_i |> 
      arrange(Time) |>  
      dplyr::select(Time, OM_horse, Environment) |> 
      filter(Environment == "Food") 
    
    char_Add_tOM_horse <- df_OM_horse_Add_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_Add_OM_horse  <- df_OM_horse_Add_i |> pull(OM_horse) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Add_horse <- paste(rep("Add, ", length(df_OM_horse_Add_i$Time)), collapse = "")
    
    # Density data
    df_density_i <- df_DEB_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Density != lag(Density))
    
    char_tdensity <- df_density_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_density  <- df_density_i |> pull(Density) |>  paste(collapse = ", ") # Nb individual in a cosm
    
    # Soil quantity in cosm (g)
    df_Soilw_i <- df_DEB_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Soil_w != lag(Soil_w))
    
    char_tsoilw <- df_Soilw_i |> pull(Time)   |>  paste(collapse = ", ") # Time in days
    char_soilw  <- df_Soilw_i |> pull(Soil_w) |>  paste(collapse = ", ") # Nb individual in a cosm
    
    # MCSim syntaxe
    
    text_init <- paste("    Init= 1;", sep="")
    
    text_Winit <- paste("    Winit=", subset(df_DEB_i, Time==0)$Weight, ";", sep="")
    
    text_Soil <- paste("    Soil=", subset(df_DEB_i, Time==0)$Soil_type, ";", sep="")
    
    text_OM_soil_init <- paste("    OM_soil_t0=", subset(df_DEB_i, Time==0)$OM_soil, ";", sep="")
    text_OM_horse_init <- paste("    OM_horse_t0=", subset(df_DEB_i, Time==0)$OM_horse, ";", sep="")
    
    
    text_OM_soil_Replace <- paste("    OM_soil_replace=Events(OM_soil,\n                          ", 
                                  length(df_OM_soil_i$Time), ",\n                          ", 
                                  char_Replace_tOM_soil,",\n                          ", 
                                  char_Replace_soil, "\n                           ",
                                  char_Replace_OM_soil, ");", sep="")
    
    
    
    text_OM_horse_Replace <- paste("    OM_horse_replace=Events(OM_horse,\n                           ", 
                                   length(df_OM_horse_Replace_i$Time), ",\n                           ", 
                                   char_Replace_tOM_horse,",\n                           ", 
                                   char_Replace_horse, "\n                            ",
                                   char_Replace_OM_horse, ");", sep="")
    
    text_OM_horse_Add <- paste("    OM_horse_add=Events(OM_horse,\n                           ", 
                               length(df_OM_horse_Add_i$Time), ",\n                           ", 
                               char_Add_tOM_horse,",\n                           ", 
                               char_Add_horse, "\n                            ",
                               char_Add_OM_horse, ");
                               ", sep="")
    
    
    
    if(OM_diff){
      text_OM_final <- paste(
        text_OM_soil_Replace,
        text_OM_horse_Replace,
        text_OM_horse_Add,
        sep="\n"
      )
    }else {
      text_OM_final <- text_OM_tot
    }
    
    
    
    if(length(df_density_i$Density) == 1){
      text_density <- paste("    Dens=",char_density, ";", sep="")
    }else{
      text_density <- paste("    Dens=NDoses(", 
                            length(df_density_i$Time), ",\n                ", 
                            char_density,",\n                ", 
                            char_tdensity, ");", sep="")
    }
    
    if (length(df_Soilw_i$Time)==1) {
      text_WeightSoilCosm <- paste0("    WeightSoilCosm = ", df_Soilw_i$Soil_w,";")
    }else {
      text_WeightSoilCosm <- paste("    WeightSoilCosm=NDoses(", 
                                   length(df_Soilw_i$Time), ",\n                ", 
                                   char_soilw,",\n                ", 
                                   char_tsoilw, ");", sep="")
    }
    
    text_Texp <- paste("    Texp=",subset(df_DEB_i, Time==0)$Texp, ";", sep="")
    
#     text_print <- paste0("    PrintStep(", Endpoints_print, ",
# 				0, ", Time_sim_i, ", ", step_sim,");")
    
    # Weight data times
    df_Ww_i <- df_DEB_i |> 
      arrange(Time)             
    char_tW <- df_Ww_i |> pull(Time)   |>  paste(collapse = ", ") # Time in days
    
    
    # Reproduction data times
    df_R_i <- df_DEB_i |> 
      arrange(Time)
    Time_xp_i <- max(df_DEB_i$Time)
    
    if (ID_experiment_i %in% XP_no_R_data) {
      char_tR <- "0"
    } else {
      char_tR <- df_R_i |> pull(Time)         |>  paste(collapse = ", ") # Time in days
    }
    
    
    
    text_print <- paste0(
      "    Print(Weight,", char_tW, ");
    Print(Reproduction,", char_tR,");"
    )
    
    
    char_i <- paste(
      paste("Simulation { #", ID_experiment_i), 
      
      text_init,
      text_WeightSoilCosm,
      text_Winit,
      text_OM_soil_init,
      text_OM_horse_init,
      text_Soil,
      text_Texp,
      text_density,
      text_OM_final,
      text_print,
      
      paste("}"),
      sep="\n"
      
    )

    text_full_i <- paste(
      text_start,
      char_i,
      text_end,
      sep = "\n"
    )
    
    writeLines(
      text_full_i,
      here::here(file_path, paste0("DEB_setpoint_predobs_",compteur_exp,".in"))
    )
    
  } # End for loop

}

f_read_predobs <- function(file_path, Nb_experiment, df_data_obs){
  
  Sim.Res.Exp <- lapply(1:Nb_experiment, function(j){
    read_tsv(file.path(file_path, paste0("DEB_setpoint_predobs_", j, ".out"))) |> 
      mutate(No_sim = j)
  }) |> 
    bind_rows()
  
  df_pred <- NULL
  
  for (i in 1:Nb_experiment){
    
    ID_i <- unique(df_data_obs$ID_experiment)[i]
    
    Sim.Res.Exp_i  <- subset(Sim.Res.Exp, No_sim == i) |> 
      dplyr::select(where(~ !all(is.na(.x))))
    df_data_obs_i <- subset(df_data_obs, ID_experiment == ID_i)
    
    Time_pred = df_data_obs_i$Time # !!!
    Nom_Endpoints = c("Weight", "Reproduction")
    Nsortie = 2
    MPV = IC_min = IC_max = Endpoint = index = NULL
    
    for (j in 1:length(Nom_Endpoints)){
      
      index =grep(paste0("^",Nom_Endpoints[j],"_") , colnames(Sim.Res.Exp_i), fixed=FALSE)
      
      
      Data_sim = as.matrix( Sim.Res.Exp_i[ , index] )
      
      MPV    = c(  MPV,  as.numeric( Data_sim[nrow(Data_sim), ] ) )
      IC_min = c(IC_min,  apply(Data_sim, 2, quantile,  0.025, na.rm=TRUE) )
      IC_max = c(IC_max,  apply(Data_sim, 2, quantile,  0.975, na.rm=TRUE) ) 
      
      Endpoint = c(Endpoint, rep(Nom_Endpoints[j], ncol(Data_sim)))
    }
    
    df_pred_i <- data.frame(
      predict.endpoint = MPV,
      Time    = rep(Time_pred, Nsortie),
      low     = IC_min,
      up      = IC_max,
      Endpt   = Endpoint
    )
    
    df_pred_i <- df_pred_i %>%
      mutate(ID_experiment=ID_i)
    
    df_pred <- rbind(df_pred, df_pred_i)
  }
  
  return(df_pred)
}

f_In_Setpoint_full <- function(file_path, l_Experiments, Adults_alone, OM_diff, l_param_name, step_sim, Endpoints_print){
  
  
  df_DEB_mean <- f_import_data_DEB(l_Experiments, Adults_alone) |> 
    mutate(
      OM_soil = signif(OM_soil, 3),
      OM_horse = signif(OM_horse, 3)
    )  |> 
    mutate(across(where(is.numeric), ~ replace_na(.x, -1)))
  
  compteur_exp <- 0 
  
  for(ID_experiment_i in unique(df_DEB_mean$ID_experiment)){
    text_full_i <- ""
    char_i <- ""
    compteur_exp <- compteur_exp + 1
    
    text_start <- paste0('#### DEB A. caliginosa
#===============================================

SetPoints("DEB_setpoint_full_', compteur_exp, '.out", "tab_setpoint.out", 0,', l_param_name|>  paste(collapse = ", "),');

Integrate(Lsodes, 1E-5, 1E-6, 1);

########## Experiments ################################################
')
    
    text_end <- "
  End."
    
    # Data corresponding to this experiment
    df_DEB_i <- subset(df_DEB_mean, ID_experiment == ID_experiment_i)
    
    Time_sim_i <- max(df_DEB_i$Time)
    
    # Soil OM data
    df_OM_soil_i <- df_DEB_i |> 
      arrange(Time) |>  
      dplyr::select(Time, OM_soil, Environment) |> 
      filter(Environment == "Soil_change") 
    
    char_Replace_tOM_soil <- df_OM_soil_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_Replace_OM_soil  <- df_OM_soil_i |> pull(OM_soil) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Replace_soil <- paste(rep("Replace, ", length(df_OM_soil_i$Time)), collapse = "")
    
    df_OM_horse_Replace_i <- df_DEB_i |> 
      arrange(Time) |>  
      dplyr::select(Time, OM_horse, Environment) |> 
      filter(Environment == "Soil_change") 
    
    char_Replace_tOM_horse <- df_OM_horse_Replace_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_Replace_OM_horse  <- df_OM_horse_Replace_i |> pull(OM_horse) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Replace_horse <- paste(rep("Replace, ", length(df_OM_horse_Replace_i$Time)), collapse = "")
    
    df_OM_horse_Add_i <- df_DEB_i |> 
      arrange(Time) |>  
      dplyr::select(Time, OM_horse, Environment) |> 
      filter(Environment == "Food") 
    
    char_Add_tOM_horse <- df_OM_horse_Add_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_Add_OM_horse  <- df_OM_horse_Add_i |> pull(OM_horse) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Add_horse <- paste(rep("Add, ", length(df_OM_horse_Add_i$Time)), collapse = "")
    
    # Density data
    df_density_i <- df_DEB_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Density != lag(Density))
    
    char_tdensity <- df_density_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_density  <- df_density_i |> pull(Density) |>  paste(collapse = ", ") # Nb individual in a cosm
    
    # Soil quantity in cosm (g)
    df_Soilw_i <- df_DEB_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Soil_w != lag(Soil_w))
    
    char_tsoilw <- df_Soilw_i |> pull(Time)   |>  paste(collapse = ", ") # Time in days
    char_soilw  <- df_Soilw_i |> pull(Soil_w) |>  paste(collapse = ", ") # Nb individual in a cosm
    
    # MCSim syntaxe
    
    text_init <- paste("    Init= 1;", sep="")
    
    text_Winit <- paste("    Winit=", subset(df_DEB_i, Time==0)$Weight, ";", sep="")
    
    text_Soil <- paste("    Soil=", subset(df_DEB_i, Time==0)$Soil_type, ";", sep="")
    
    text_OM_soil_init <- paste("    OM_soil_t0=", subset(df_DEB_i, Time==0)$OM_soil, ";", sep="")
    text_OM_horse_init <- paste("    OM_horse_t0=", subset(df_DEB_i, Time==0)$OM_horse, ";", sep="")
    
    
    text_OM_soil_Replace <- paste("    OM_soil_replace=Events(OM_soil,\n                          ", 
                                  length(df_OM_soil_i$Time), ",\n                          ", 
                                  char_Replace_tOM_soil,",\n                          ", 
                                  char_Replace_soil, "\n                           ",
                                  char_Replace_OM_soil, ");", sep="")
    
    
    
    text_OM_horse_Replace <- paste("    OM_horse_replace=Events(OM_horse,\n                           ", 
                                   length(df_OM_horse_Replace_i$Time), ",\n                           ", 
                                   char_Replace_tOM_horse,",\n                           ", 
                                   char_Replace_horse, "\n                            ",
                                   char_Replace_OM_horse, ");", sep="")
    
    text_OM_horse_Add <- paste("    OM_horse_add=Events(OM_horse,\n                           ", 
                               length(df_OM_horse_Add_i$Time), ",\n                           ", 
                               char_Add_tOM_horse,",\n                           ", 
                               char_Add_horse, "\n                            ",
                               char_Add_OM_horse, ");
                               ", sep="")
    
    
    
    if(OM_diff){
      text_OM_final <- paste(
        text_OM_soil_Replace,
        text_OM_horse_Replace,
        text_OM_horse_Add,
        sep="\n"
      )
    }else {
      text_OM_final <- text_OM_tot
    }
    
    
    
    if(length(df_density_i$Density) == 1){
      text_density <- paste("    Dens=",char_density, ";", sep="")
    }else{
      text_density <- paste("    Dens=NDoses(", 
                            length(df_density_i$Time), ",\n                ", 
                            char_density,",\n                ", 
                            char_tdensity, ");", sep="")
    }
    
    if (length(df_Soilw_i$Time)==1) {
      text_WeightSoilCosm <- paste0("    WeightSoilCosm = ", df_Soilw_i$Soil_w,";")
    }else {
      text_WeightSoilCosm <- paste("    WeightSoilCosm=NDoses(", 
                                   length(df_Soilw_i$Time), ",\n                ", 
                                   char_soilw,",\n                ", 
                                   char_tsoilw, ");", sep="")
    }
    
    text_Texp <- paste("    Texp=",subset(df_DEB_i, Time==0)$Texp, ";", sep="")
    
    #     text_print <- paste0("    PrintStep(", Endpoints_print, ",
    # 				0, ", Time_sim_i, ", ", step_sim,");")
    
    text_print <- paste0("    PrintStep(", Endpoints_print, ",
				0, ", Time_sim_i, ", ", step_sim,");")
    
    
    char_i <- paste(
      paste("Simulation { #", ID_experiment_i), 
      
      text_init,
      text_WeightSoilCosm,
      text_Winit,
      text_OM_soil_init,
      text_OM_horse_init,
      text_Soil,
      text_Texp,
      text_density,
      text_OM_final,
      text_print,
      
      paste("}"),
      sep="\n"
      
    )
    
    text_full_i <- paste(
      text_start,
      char_i,
      text_end,
      sep = "\n"
    )
    
    writeLines(
      text_full_i,
      here::here(file_path, paste0("DEB_setpoint_full_",compteur_exp,".in"))
    )
    
  } # End for loop
  
}

f_read_Setpoint_full <- function(file_path, Nb_experiment, select_times, step_sim){
  
  df_pred_full <- NULL
  
  for (i in 1:Nb_experiment){
    
    Sim.Res.Exp_i <- readr::read_tsv(
      file.path(file_path, paste0("DEB_setpoint_full_", i, ".out"))
    ) |>
      dplyr::mutate(No_sim = i) |> 
      dplyr::select(where(~ !all(is.na(.x))))
    
    ID_i <- unique(df_data_obs$ID_experiment)[i]
    
    df_data_obs_i <- subset(df_data_obs, ID_experiment == ID_i)
    
    Time_pred = seq(0, max(df_data_obs_i$Time, na.rm=TRUE), step_sim) # !!!
    
    
    Nom_Endpoints = c("Weight", "Reproduction")
    Nsortie = 2
    MPV = IC_min = IC_max = Endpoint = index = NULL
    
    for (j in 1:length(Nom_Endpoints)){
      
      index =grep(paste0("^",Nom_Endpoints[j],"_") , colnames(Sim.Res.Exp_i), fixed=FALSE)
      
      
      Data_sim = as.matrix( Sim.Res.Exp_i[ , index] )
      
      MPV    = c(  MPV,  as.numeric( Data_sim[nrow(Data_sim), ] ) )
      IC_min = c(IC_min,  apply(Data_sim, 2, quantile,  0.025, na.rm=TRUE) )
      IC_max = c(IC_max,  apply(Data_sim, 2, quantile,  0.975, na.rm=TRUE) ) 
      
      Endpoint = c(Endpoint, rep(Nom_Endpoints[j], ncol(Data_sim)))
    }
    
    df_pred_i <- data.frame(
      predict.endpoint = MPV,
      Time    = rep(Time_pred, Nsortie),
      low     = IC_min,
      up      = IC_max,
      Endpt   = Endpoint
    )
    
    df_pred_i <- df_pred_i %>%
      mutate(ID_experiment=ID_i) |> 
      filter (Time %in% select_times) |> 
      mutate(
        predict.endpoint = case_when(
         ( predict.endpoint <= 5e-6 & Endpt == "Reproduction") ~ 0,
          .default = predict.endpoint
        )
        )
    
    df_pred_full <- rbind(df_pred_full, df_pred_i)
  }
  
  saveRDS(df_pred_full, file.path(file_path, paste0("Sim.Res.Exp.full.rds")))
  
  return(df_pred_full)
}

f_In_simulations_sim <- function(file_path, l_Experiments, Adults_alone, OM_diff, text_param, Endpoints_print){
  
  df_DEB_mean <- read.csv(here::here("data/DEB_sim_data.csv")) |> 
    mutate(
      OM_soil = signif(OM_soil, 3),
      OM_horse = signif(OM_horse, 3)
    )  |> 
    mutate(across(where(is.numeric), ~ replace_na(.x, -1)))
  
  char_experiments <- ""
  
  for(ID_experiment_i in unique(df_DEB_mean$ID_experiment)){
    
    # Data corresponding to this specific cosm
    df_DEB_i <- subset(df_DEB_mean, ID_experiment == ID_experiment_i)
    
    # Soil OM data
    df_OM_soil_i <- df_DEB_i |> 
      arrange(Time) |>  
      dplyr::select(Time, OM_soil, Environment) |> 
      filter(Environment == "Soil_change") 
    
    char_tOM_soil <- df_OM_soil_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_OM_soil  <- df_OM_soil_i |> pull(OM_soil) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Add_soil <- paste(rep("Replace, ", length(df_OM_soil_i$Time)), collapse = "")
    
    # Horse dung OM data
    df_OM_horse_i <- df_DEB_i |> 
      arrange(Time) |>                    
      dplyr::select(Time, OM_horse, Environment) |> 
      filter(Environment %in% c("Soil_change", "Food")) |> 
      mutate(
        Word = case_when(
          Environment == "Soil_change" ~ "Replace",
          Environment == "Food" ~ "Add"
        )
      )
    
    char_tOM_horse<- df_OM_horse_i  |> pull(Time)     |>  paste(collapse = ", ") # Time in days
    char_OM_horse  <- df_OM_horse_i |> pull(OM_horse) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Add_horse <- df_OM_horse_i |> pull(Word)     |>  paste(collapse = ", ")
    
    # OM tot data
    df_OM_tot_i <- df_DEB_i |> 
      arrange(Time) |>                    
      dplyr::select(Time, OM_horse, OM_soil, Environment) |> 
      filter(Environment %in% c("Soil_change", "Food")) |> 
      mutate(
        OM_soil = case_when(
          Environment == "Soil_change" ~ OM_soil,
          Environment == "Food" ~ 0
        ),
        OM_tot = OM_soil+OM_horse,
        Word = case_when(
          Environment == "Soil_change" ~ "Replace",
          Environment == "Food" ~ "Add"
        )
      )
    
    char_tOM_tot<- df_OM_tot_i  |> pull(Time)     |>  paste(collapse = ", ") # Time in days
    char_OM_tot  <- df_OM_tot_i |> pull(OM_tot) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Add_tot <- df_OM_tot_i |> pull(Word)     |>  paste(collapse = ", ")
    
    # Density data
    df_density_i <- df_DEB_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Density != lag(Density))
    
    char_tdensity <- df_density_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_density  <- df_density_i |> pull(Density) |>  paste(collapse = ", ") # Nb individual in a cosm
    
    # Soil quantity in cosm (g)
    df_Soilw_i <- df_DEB_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Soil_w != lag(Soil_w))
    
    char_tsoilw <- df_Soilw_i |> pull(Time)   |>  paste(collapse = ", ") # Time in days
    char_soilw  <- df_Soilw_i |> pull(Soil_w) |>  paste(collapse = ", ") # Nb individual in a cosm
    
    # MCSim syntaxe
    
    text_init <- paste("    Init= 1;", sep="")
    
    text_Winit <- paste("    Winit=", subset(df_DEB_i, Time==0)$Weight, ";", sep="")
    text_Soil <- paste("    Soil=", subset(df_DEB_i, Time==0)$Soil_type, ";", sep="")
    
    text_OM_soil <- paste("    event_OM_soil=Events(OM_soil,\n                          ", 
                          length(df_OM_soil_i$Time), ",\n                          ", 
                          char_tOM_soil,",\n                          ", 
                          char_Add_soil, "\n                           ",
                          char_OM_soil, ");", sep="")
    
    text_OM_horse <- paste("    event_OM_horse=Events(OM_horse,\n                           ", 
                           length(df_OM_horse_i$Time), ",\n                           ", 
                           char_tOM_horse,",\n                           ", 
                           char_Add_horse, ",\n                            ",
                           char_OM_horse, ");", sep="")
    
    text_OM_tot <- paste("    event_OM_tot=Events(OM_tot,\n                           ", 
                         length(df_OM_tot_i$Time), ",\n                           ", 
                         char_tOM_tot,",\n                           ", 
                         char_Add_tot, ",\n                            ",
                         char_OM_tot, ");", sep="")
    
    if(OM_diff){
      text_OM_final <- paste(
        text_OM_soil,
        text_OM_horse,
        sep="\n"
      )
    }else {
      text_OM_final <- text_OM_tot
    }
    
    if(length(df_density_i$Density) == 1){
      text_density <- paste("    Dens=",char_density, ";", sep="")
    }else{
      text_density <- paste("    Dens=NDoses(", 
                            length(df_density_i$Time), ",\n                ", 
                            char_density,",\n                ", 
                            char_tdensity, ");", sep="")
    }
    
    if (length(df_Soilw_i$Time)==1) {
      text_WeightSoilCosm <- paste0("    WeightSoilCosm = ", df_Soilw_i$Soil_w,";")
    }else {
      text_WeightSoilCosm <- paste("    WeightSoilCosm=NDoses(", 
                                   length(df_Soilw_i$Time), ",\n                ", 
                                   char_soilw,",\n                ", 
                                   char_tsoilw, ");", sep="")
    }
    
    text_Texp <- paste("    Texp=",subset(df_DEB_i, Time==0)$Texp, ";", sep="")
    
    char_i <- paste(
      paste("Simulation { #", ID_experiment_i), 
      
      text_init,
      text_WeightSoilCosm,
      text_Winit,
      text_Soil,
      text_Texp,
      text_density,
      text_OM_final,
      text_print,
      
      paste("}"),
      sep="\n"
      
    )
    
    char_experiments <- paste(char_experiments, char_i, sep="\n")
    
  }
  
  text_start <- "#### DEB A. caliginosa
#===============================================

Integrate(Lsodes, 1E-6, 1E-8, 1);

########## Experiments ################################################
"
  
  text_end <- "
  End."
  
  text_full <- paste(
    text_start,
    text_param,
    char_experiments,
    text_end,
    sep = "\n"
  )
  
  writeLines(
    text_full,
    here::here(file_path, "DEB_sim_full.in")
  )
}




f_MCSim <- function(path_mod) {
  
  # 1. Read files ----
  
  Files_in <- list.files(path = path_mod, pattern = "C[0-8]\\.in")
  Files_out <- list.files(path = path_mod, pattern = "C[0-8]\\.out")
  
  tmp <- grep(
    pattern = "out.kernel",
    Files_out
  )
  if (length(tmp) > 0) {
    Files_out <- Files_out[-tmp]
  }
  
  NC <- length(Files_out)
  
  # 2. Priors ----
  
  Initialize <- readLines(paste(path_mod, Files_in[1], sep = "/"))
  
  Exl.Par1 <- grep(pattern = "*\\#Distrib*", Initialize)
  Exl.Par2 <- grep(pattern = "*\\# Distrib*", Initialize)
  Parameter <- grep(pattern = "*Distrib*", Initialize)
  Parameter <- Parameter[!Parameter %in% Exl.Par1]
  Parameter <- Parameter[!Parameter %in% Exl.Par2]
  
  np <- length(Parameter) # Used just to be able to count the number of parameters with individual variation
  
  Experiment <- grep(pattern = "Experiment", Initialize, value = T)
  nexperiement <- length(Experiment)
  ParTable <- as.data.frame(matrix(NA, np, 6))
  
  colnames(ParTable) <- c("Nom", "Distribution", "P1", "P2", "P3", "P4")
  for (i in 1:np) {
    temp <- gsub("\t", "", Initialize[Parameter[i]])
    temp <- strsplit(temp, "\\#")[[1]][1]
    temp <- strsplit(temp, ".\\(|\\)|\\,| ")[[1]]
    temp <- temp[temp != ""]
    
    for (j in 2:(length(temp) - 1)) {
      ParTable[i, j - 1] <- temp[j]
    }
  }
  
  ParTable[, 3:ncol(ParTable)] <- apply(ParTable[, 3:ncol(ParTable)], 2, as.numeric)
  
  # 3. Results ----
  
  Data_T <- Data_V <- Data_L <- vector("list", length = NC) # sous forme liste
  
  for (i in 1:NC) {
    data <- read.table(eval(paste(path_mod, Files_out[i], sep = "/")), header = TRUE)
    
    if (i == 1) {
      times <- data[, 1]
    }
    tmp <- as.matrix(data[, -1])
    
    Data_L[[i]] <- tmp # Fichier complet pour les NC chaines
    Data_T[[i]] <- tmp[, 1:np] # Fichier parametres pour les NC chaines
    Data_V[[i]] <- tmp[, (np + 1):ncol(tmp)] # Fichier Vraisemblance pour les NC chaines
    names(Data_L)[i] <- paste0("Chain", i)
    names(Data_T)[i] <- paste0("Chain", i)
    names(Data_V)[i] <- paste0("Chain", i)
  }
  
  # Chains
  
  Niter <- nrow(Data_T[[1]]) # number of lines
  PI <- as.numeric(times[nrow(Data_T[[1]])] - times[(nrow(Data_T[[1]]) - 1)]) # pas de l'iteration
  
  nb <- 10000
  if (nb >= Niter) {
    nb <- Niter / 2
  } # control nb < Niter
  
  IterSelect <- ((Niter - nb / PI):Niter)
  
  Data_All <- NULL
  for (i in 1:NC) {
    Data_T[[i]] <- mcmc(Data_T[[i]], start = 1, end = Niter, thin = PI)
    Data_All <- rbind(Data_All, Data_L[[i]][IterSelect, ])
  }
  res.mcmc <- coda::as.mcmc.list(Data_T)
  df_res.mcmc <- ggs(res.mcmc)
  resmc <- summary(res.mcmc)
  
  Mode <- which.max(Data_All[, ncol(Data_All)]) # MPV computation
  Result_mode <- Data_All[Mode[1], ]
  
  Min <- apply(Data_All, 2, min)
  Max <- apply(Data_All, 2, max)
  Result_sd <- apply(Data_All, 2, sd)
  
  Result_mean <- apply(Data_All, 2, mean)
  Result_quantil <- apply(Data_All, 2, function(x) {
    quantile(x, c(0.025, 0.975))
  })
  
  Res <- rbind(Result_mode, Result_quantil, Min, Max, Result_mean, Result_sd)
  
  # 4. AIC and BIC calculations ----
  
  
  # 5. tab_setpoint.out construction ----
  
  ## 5.1. tab_setpoint pop ----
  
  Niter_chain <- 333 # Niter_chain last lines selected to setpoint.in / doit pas être trop grand sinon impossible de faire les graphs (marche au moins avec 1000)
  #Niter_chain <- 1
  
  Selected_chain <- NULL
  
  for (i in 1:NC) {
    Selected_chain <- rbind(Selected_chain, Data_L[[i]][(nrow(Data_L[[i]]) - Niter_chain + 1):nrow(Data_L[[i]]), ]) # Ajout du +1
  }
  
  tmp <- -c((ncol(Selected_chain) - 2):ncol(Selected_chain))
  Selected_chain <- Selected_chain[, tmp]
  MVP <- Result_mode[tmp] # "Meilleure" solution dans toutes les chaines
  
  Selected_chain <- rbind(Selected_chain, as.list(MVP)) # Ajout du as.list()
  df_Selected_chain <- as.data.frame(Selected_chain)
  
  tab_setpoint <- df_Selected_chain[, 1:np]
  
  # 6. Likelihood and deviance ratio ----
  
  tmp.names <- Devc <- LnData <- LnPost_start <- LnPost <- NULL
  
  for (i in 1:NC) {
    temp <- -2 * Data_V[[i]][IterSelect, ncol(Data_V[[i]])] # -2log(vraisemblance), tableau
    Devc <- cbind(Devc, temp)
    LnPost_start <- cbind(LnPost_start, Data_V[[i]][seq(100, 2 * Niter / 4), ncol(Data_V[[i]]), drop = FALSE])
    LnPost <- cbind(LnPost, Data_V[[i]][IterSelect, ncol(Data_V[[i]])])
    LnData <- cbind(LnData, Data_V[[i]][IterSelect, (ncol(Data_V[[i]])-1)])
    tmp.names <- c(tmp.names, paste("Chain", i))
  }
  
  df_LnPost_start <- data.frame(iteration = (100):(2 * Niter / 4), LnPost_start) %>%
    pivot_longer(-iteration, names_to = "Chain", values_to = "LnPosterior")
  
  df_LnPost <- data.frame(iteration = IterSelect[1:nrow(Devc)], LnPost) %>%
    pivot_longer(-iteration, names_to = "Chain", values_to = "LnPosterior")
  
  df_Devc <- data.frame(iteration = IterSelect[1:nrow(Devc)], RatioDeviance = Devc / min(Devc)) %>%
    pivot_longer(-iteration, names_to = "Chain", values_to = "RatioDeviance")
  
  df_LnData <- data.frame(iteration = IterSelect[1:nrow(Devc)], LnData) %>%
    pivot_longer(-iteration, names_to = "Chain", values_to = "LnData")
  
  # 7. Outputs ----
  tab_setpoint <- tab_setpoint |>
    mutate(across(where(is.list), ~ simplify2array(.) |> unlist()))
  
  write.table(tab_setpoint, file = paste(path_mod, "tab_setpoint.out", sep = "/"), quote = FALSE)
  
  
  MCMC_out <- list(
    NC              = NC,
    Files_in        = Files_in,
    Files_out       = Files_out,
    Nb_param        = np,
    Nb_experiement  = nexperiement,
    Priors          = ParTable,
    Nb_iter_kept    = nb,         # Last iterations used to describe the posterior distributions
    Chains          = df_res.mcmc,
    df_LnPost_start = df_LnPost_start,
    df_LnPost       = df_LnPost,
    df_LnData       = df_LnData,
    df_Devc         = df_Devc,
    Summary_res     = Res,
    N_iter_setpoint = Niter_chain # last lines selected to setpoint.in (cannot be too big)
  )
  
  return(MCMC_out)
}





f_MCSim_read_sim <- function(path_sim_out){
  library(readr)
  
  # Read file
  lines <- read_lines(path_sim_out)
  
  # Find the begining of each simulation
  indices_sim <- which(str_detect(lines, "^Results of Simulation"))
  
  if (length(indices_sim)==0) {
    indices_sim <- c(0)
  }
  
  n_simulations <- length(indices_sim)
  
  # Initialisation of the list for the dataframes of each simulation
  list_simulations <- list()
  
  idx_begin <- indices_sim + 1
  idx_end <- c((indices_sim[-1] - 1), length(lines))
  
  # For each simulation
  for (i in seq_along(indices_sim)) {
    
    bloc <- lines[idx_begin[i]:idx_end[i]]
    
    # To data.frame
    df <- read_table2(paste(bloc, collapse = "\n"))
    df$No_sim <- i  # Ad the ID of the simulation
    list_simulations[[i]] <- df
  }
  
  # Combine data.frames
  df_total <- bind_rows(list_simulations) |> 
    mutate(No_sim = as.factor(No_sim))
  
  return(df_total)
  
}

f_MCSim_read_MCMC_out <- function(path_sim_out, times_pred){
  library(readr)
  
  # Read file
  lines <- read_lines(path_sim_out)
  
  # Find the begining of each simulation
  indices_sim <- which(str_detect(lines, "^Results of Simulation"))
  
  if (length(indices_sim)==0) {
    indices_sim <- c(0)
  }
  
  n_simulations <- length(indices_sim)
  
  # Initialisation of the list for the dataframes of each simulation
  list_simulations <- list()
  
  idx_begin <- indices_sim + 1
  idx_end <- c((indices_sim[-1] - 1), length(lines))
  
  # For each simulation
  for (i in seq_along(indices_sim)) {
    
    bloc <- lines[idx_begin[i]:idx_end[i]]
    
    # To data.frame
    df <- read_table2(paste(bloc, collapse = "\n"))
    df$No_sim <- i  # Ad the ID of the simulation
    list_simulations[[i]] <- df
  }
  
  # Combine data.frames
  df_total <- bind_rows(list_simulations) |> 
    mutate(No_sim = as.factor(No_sim))
  

  Time_pred <- times_pred
  
  Nom_Endpoints <- c("freal", "kapreal", "Maturity", "Energy",
                     "Organic_matter", "Weight", "Reproduction")
  
  df_long <- lapply(Nom_Endpoints, function(nom) {
    
    index <- grep(paste0("^", nom, "_"), colnames(df_total))
    
    df_MC_light_param <- df_MC_light[,1:8] |>
      pivot_longer(
        cols = -Iter,              # garde Iter comme identifiant
        names_to = "Parameter",     # noms des paramètres
        values_to = "Value"        # valeurs associées
      ) |>
      mutate(Type = "MC")
    
    df_bornes <- data.frame(
      Parameter = unique(df_MC_light_param$Parameter),
      Mean = c(1117.2, 1.5, NA, 0.5, 100, 200, 0.2)
    )
      dplyr::select(No_sim, all_of(index)) |>
      group_by(No_sim) |>
      mutate(Sim = row_number()) |>
      ungroup() |>
      pivot_longer(
        cols = -c(No_sim, Sim),
        names_to = "Time",
        names_pattern = paste0("^", nom, "_(.*)$"),
        names_transform = list(Time = as.numeric),
        values_to = "Value"
      ) |>
      mutate(
        Param = nom,
        Time = rep(Time_pred, max(Sim))
        )
    
  }) |>
    bind_rows() |>
    relocate(No_sim, Sim, Param, Time, Value)
  
  df_final <- df_long |>
    rename(Endpt = Param) |>
    dplyr::select(Time, Endpt, Value, No_sim, Sim)
  
  
  return(df_final)
  
}

f_MCSim_ind <- function(path_mod){
  
  # Read files
  
  Files_in  = list.files(path = path_mod, pattern = "[0-3]\\.in")
  Files_out = list.files(path = path_mod, pattern = "[0-3]\\.out")
  
  tmp = grep(
    pattern="out.kernel", 
    Files_out
  )
  if(length(tmp)>0){ Files_out = Files_out [ -tmp ] }
  
  NC = length(Files_out)
  
  # Priors
  
  Initialize = readLines(paste(path_mod,Files_in[1], sep="/"))
  
  Exl.Par1   = grep(pattern="*\\#Distrib*", Initialize)  
  Exl.Par2   = grep(pattern="*\\# Distrib*", Initialize) 
  Parameter = grep(pattern="*Distrib*", Initialize)
  Parameter =  Parameter[!Parameter %in% Exl.Par1]
  Parameter =  Parameter[!Parameter %in% Exl.Par2]
  
  np_tmp = length(Parameter) # Used just to be able to count the number of parameters with individual variation
  
  Experiment = grep(pattern="Experiment", Initialize, value =T)
  nexperiement =  length(Experiment)
  ParTable =as.data.frame(matrix(NA,np_tmp,6))
  
  colnames(ParTable)=c("Nom","Distribution","P1","P2","P3","P4")
  for (i in 1:np_tmp){ 
    
    temp = gsub("\t","", Initialize[Parameter[i]] )
    temp = strsplit( temp, "\\#")[[1]][1]
    temp = strsplit( temp, ".\\(|\\)|\\,| ")[[1]]
    temp = temp[temp!=""]
    
    for(j in 2:(length(temp)-1) ) {  ParTable[i,j-1] = temp[j] }
  }
  
  Nb_param_IndVar <- length(gsub("^Vr_", "", ParTable$Nom[grepl("^Vr_", ParTable$Nom)]))
  np <- np_tmp - Nb_param_IndVar
  
  ParTable = head(ParTable,-Nb_param_IndVar)
  
  ParTable[,3:ncol(ParTable)]=apply(ParTable[,3:ncol(ParTable)],2,as.numeric)
  
  # Results 
  
  Data_T = Data_V = Data_L =  vector("list", length=NC) # sous forme liste
  
  for (i in 1:NC){
    
    data <- read.table( eval(paste(path_mod,Files_out[i], sep="/"))   , header = TRUE)
    
    if(i ==1){  times = data[,1] }
    tmp = as.matrix(data[,-1]) 
    
    Data_L[[i]] = tmp                                     # Fichier complet pour les NC chaines
    Data_T[[i]] = tmp[,1:np]                              # Fichier parametres pour les NC chaines
    Data_V[[i]] = tmp[,(np+1):ncol(tmp)]                  # Fichier Vraisemblance pour les NC chaines
    names(Data_L)[i]= paste0("Chain",i) 
    names(Data_T)[i]= paste0("Chain",i) 
    names(Data_V)[i]= paste0("Chain",i)
    
  }
  
  # Chains
  
  Niter = nrow(Data_T[[1]]) # number of lines
  PI    = as.numeric(times[nrow(Data_T[[1]])] - times[(nrow(Data_T[[1]])-1)]) # pas de l'iteration
  
  nb = 10000
  if(nb>=Niter){nb=Niter/2} #control nb < Niter
  
  IterSelect = (( Niter- nb / PI ): Niter  )
  
  Data_All = NULL
  for (i in 1:NC){
    Data_T[[i]] =  mcmc(Data_T[[i]], start = 1, end = Niter , thin = PI)
    Data_All = rbind(Data_All, Data_L[[i]][IterSelect,])  }
  res.mcmc =  coda::as.mcmc.list(Data_T)
  df_res.mcmc = ggs(res.mcmc)
  resmc  = summary(res.mcmc)
  
  Mode    = which.max(Data_All[,ncol(Data_All)]) # MPV computation
  Result_mode = Data_All[ Mode[1], ] 
  
  Min = apply(Data_All,2,min)
  Max = apply(Data_All,2,max)
  Result_sd     = apply( Data_All , 2, sd)
  Result_mean   = apply( Data_All , 2, mean)
  Result_quantil = apply( Data_All , 2, function(x){ quantile(x, c(0.025, 0.975) ) }  )
  
  Res = rbind(Result_mode, Result_quantil, Min, Max, Result_mean, Result_sd )
  
  # AIC and BIC calculations
  
  # MaxVraiss = Data_All[ Mode[1], ncol(Data_All)] 
  # 
  # AIC <- -2*MaxVraiss+2*np
  # Ndata = 13 # number of observed data used
  # BIC <- -2*MaxVraiss+log(Ndata)*np
  
  # tab_setpoint.out construction
  
  # tab_setpoint pop
  
  Niter_chain = 333 # Niter_chain last lines selected to setpoint.in / doit pas être trop grand sinon impossible de faire les graphs (marche au moins avec 1000)
  
  Selected_chain = NULL
  
  for (i in 1:NC){
    Selected_chain = rbind( Selected_chain , Data_L[[i]][ (nrow (Data_L[[i]])- Niter_chain+1) : nrow (Data_L[[i]]),] ) # Ajout du +1
  }
  
  tmp = -c( (ncol(Selected_chain)-2) : ncol(Selected_chain)  )
  Selected_chain <-Selected_chain[,tmp]
  MVP =  Result_mode[ tmp] # "Meilleure" solution dans toutes les chaines
  
  Selected_chain <-rbind( Selected_chain, as.list(MVP)) # Ajout du as.list()
  df_Selected_chain <- as.data.frame(Selected_chain)
  
  tab_setpoint <- df_Selected_chain[,1:np]
  
  
  # Setpoint individuel / per experiment
  col_ind <- df_Selected_chain[, (np+1):length(df_Selected_chain)]
  
  colnames(tab_setpoint) <- gsub("\\.1\\.$", "", colnames(tab_setpoint))
  Names_par <- gsub("\\.1\\.$", "", names(tab_setpoint))
  names_param_ind <- grep("^Vr_", Names_par, value = TRUE)
  names_param_sigma <- grep("^Sigma_", Names_par, value = TRUE)
  names_param_pop_ind <- sub("^Vr_", "", names_param_ind)
  
  col_pop <- tab_setpoint |> 
    dplyr::select(-all_of(c(names_param_sigma, names_param_ind, names_param_pop_ind)))
  
  # Version LG
  col_ind_piled <- NULL
  
  for (i in 1:nexperiement){
    
    start = Nb_param_IndVar*(i-1)+1
    end = Nb_param_IndVar*(i)
    
    col_ind_i <- col_ind[, start:end]
    names(col_ind_i) <- names_param_pop_ind
    col_ind_piled <- rbind(col_ind_piled, col_ind_i)
  }
  
  n_repeats <- nrow(col_ind_piled) / nrow(col_pop)
  
  col_pop_rep <- col_pop[rep(1:nrow(col_pop), each = n_repeats), ]
  
  tab_setpoint_ind_tmp <- cbind(col_pop_rep, col_ind_piled)
  
  # Version RB
  # tab_setpoint_ind_tmp <- cbind(col_pop,col_ind[,1:20]) 
  # colnames(tab_setpoint_ind_tmp) <- gsub("\\.1\\.(\\d+)\\.", "_ind_\\1", colnames(tab_setpoint_ind_tmp))
  
  cols_growth <- grep("^a_growth", colnames(tab_setpoint_ind_tmp), value = TRUE)
  cols_fasting <- grep("^a_fasting", colnames(tab_setpoint_ind_tmp), value = TRUE)
  cols_other <- setdiff(colnames(tab_setpoint_ind_tmp), c(cols_growth, cols_fasting))
  
  # Réorganisation
  tab_setpoint_ind <- tab_setpoint_ind_tmp[, c(cols_other, cols_growth, cols_fasting)]
  
  
  # Likelihood and deviance ratio 
  
  tmp.names <- Devc <- LnData <- LnPost_start <- LnPost <- NULL
  
  for (i in 1:NC) {
    temp <- -2 * Data_V[[i]][IterSelect, 3] # -2log(vraisemblance), tableau
    Devc <- cbind(Devc, temp)
    LnPost_start <- cbind(LnPost_start, Data_V[[i]][seq(100, 2 * Niter / 4), ncol(Data_V[[i]]), drop = FALSE])
    LnPost <- cbind(LnPost, Data_V[[i]][IterSelect, ncol(Data_V[[i]])])
    LnData <- cbind(LnData, Data_V[[i]][IterSelect, (ncol(Data_V[[i]])-1)])
    tmp.names <- c(tmp.names, paste("Chain", i))
  }
  
  df_LnPost_start <- data.frame(iteration = (100):(2 * Niter / 4), LnPost_start) %>%
    pivot_longer(-iteration, names_to = "Chain", values_to = "LnPosterior")
  
  df_LnPost <- data.frame(iteration = IterSelect[1:nrow(Devc)], LnPost) %>%
    pivot_longer(-iteration, names_to = "Chain", values_to = "LnPosterior")
  
  df_Devc <- data.frame(iteration = IterSelect[1:nrow(Devc)], RatioDeviance = Devc / min(Devc)) %>%
    pivot_longer(-iteration, names_to = "Chain", values_to = "RatioDeviance")
  
  df_LnData <- data.frame(iteration = IterSelect[1:nrow(Devc)], LnData) %>%
    pivot_longer(-iteration, names_to = "Chain", values_to = "LnData")
  
  # Outputs 
  tab_setpoint <- tab_setpoint |> 
    mutate(across(where(is.list), ~ simplify2array(.) |> unlist()))
  tab_setpoint_ind <- tab_setpoint_ind |> 
    mutate(across(where(is.list), ~ simplify2array(.) |> unlist()))
  
  write.table(tab_setpoint, file=paste(path_mod, "tab_setpoint.out", sep="/"),quote=FALSE)
  write.table(tab_setpoint_ind, file=paste(path_mod, "tab_setpoint_ind.out", sep="/"),quote=FALSE)
  
  MCMC_out <- list(
    NC = NC,
    Files_in = Files_in,
    Files_out = Files_out,
    Nb_param = np,
    Nb_experiement = nexperiement,
    Priors = ParTable, 
    Nb_iter_kept = nb, # Last iterations used to describe the posterior distributions
    Chains = df_res.mcmc,
    df_LnPost_start = df_LnPost_start,
    df_LnPost = df_LnPost,
    df_LnData = df_LnData,
    df_Devc = df_Devc,
    Summary_res = Res,
    N_iter_setpoint = Niter_chain # last lines selected to setpoint.in (cannot be too big)
  )
  
  return(MCMC_out)
  
}


f_AS_Sobol_design <- function(Sobol_type, seed, Nb_sim, CV, Processors, Nboot, Nout, l_params_median, Temperature_bol){
  
  # Nb_sim : Number of simulations
  # Nboot : Number of replicates
  # Nout : Number of outputs
  
  if (!Temperature_bol){
    l_params_median <- l_params_median[1:(length(l_params_median)-4)] # remove TAH TH TA Tref from the analysis
  }
  
  library(sensitivity)
  set.seed(seed)
  
  l_param_names <- as.character(names(l_params_median))
  Nb_param_AS <- length(l_param_names)
  
  # Design parameters
  
  l_Lower <- l_params_median - l_params_median*CV
  l_Upper <- l_params_median + l_params_median*CV
  l_Lower <- pmax(l_Lower, 0)
  
  X1 <- matrix(
    runif(
      Nb_sim * Nb_param_AS,
      min = rep(l_Lower, each = Nb_sim),
      max = rep(l_Upper, each = Nb_sim)
    ),
    nrow = Nb_sim, ncol = Nb_param_AS
  )
  
  X2 <- matrix(
    runif(
      Nb_sim * Nb_param_AS,
      min = rep(l_Lower, each = Nb_sim),
      max = rep(l_Upper, each = Nb_sim)
    ),
    nrow = Nb_sim, ncol = Nb_param_AS
  )
  
  colnames(X1) <- colnames(X2) <- l_param_names
  
  # Simulation experimental design
  if (Sobol_type == "Jansen"){
    Sobol <- soboljansen(model = NULL, X1, X2, nboot = Nboot, conf = 0.95)
  } else if (Sobol_type == "Martinez"){
    Sobol <- sobolmartinez(model = NULL, X1, X2, nboot = Nboot, conf = 0.95)
  }
  
  Sobol_param <- Sobol$X
  
  res <- apply(Sobol_param, 1, function(x) {
    as.list(x)
  })
  
  return(
    list(
      Sobol = Sobol,
      Sobol_param = res,
      X1 = X1,
      X2 = X2
    )
  )
}

f_AS_Sobol_index <- function(df_AS_sim, Sobol_info){
  
  df_AS_sim_mean <- df_AS_sim |> 
    summarise(
      across(
        .cols = is.numeric, 
        .fns = list(Mean = mean),
        na.rm = TRUE,
        .names = "{col}"
      )
    )
  
  df_AS_sim <- df_AS_sim |>
    mutate(
      Weight_Juv = as.numeric(Weight_Juv),
      across(
        c(
          Weight_Juv, Weight_Adult,
          Energy_Juv, Energy_Adult,
          Maturity_Juv, Maturity_Adult,
          Reproduction, Organic_matter
        ),
        ~ coalesce(.x, df_AS_sim_mean[[cur_column()]][1])
      )
    )
  
  # Sobol calculations
  
  Nb_out <- length(df_AS_sim)
  df_res <- data.frame()  # liste vide
  
  for (i in seq_along(df_AS_sim)) {
    
    Y <- dplyr::pull(df_AS_sim, i)
    tell(x = Sobol_info, y = Y, nboot =  Nboot, conf = 0.95)
    
    df_res_i <- data.frame(
      Variable = colnames(df_AS_sim)[i],
      FOI = Sobol_info$S[,1],
      FOI.borninf = Sobol_info$S[,4],
      FOI.bornsup = Sobol_info$S[,5],
      TI = Sobol_info$T[,1],
      TI.borninf = Sobol_info$T[,4],
      TI.bornsup = Sobol_info$T[,5]
    )
    
    df_res <- rbind(df_res, df_res_i)
  }
  
  df_res <- df_res %>%
    mutate(Parameter = rep(rownames(Sobol_info$S), Nb_out)) %>%
    arrange(desc(TI))
  
  return(df_res)
}
