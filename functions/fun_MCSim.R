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

f_import_data_DEB <- function(l_Experiments, Adults_alone = FALSE){
  
  if(!Adults_alone){
    df_data <- read.csv(here::here("data/Data_DEB_mean.csv")) |> 
      filter(((ID_experiment %in% c("Bart2019_XP1_2", "Bart2019_XP4_2"))|!(Status == "A" & Nb_ind == 1)) & ID_experiment %in% l_Experiments)
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

## DEB -----

f_write_DEB_models <- function() {
  
  l_DEB_models <- c(
    "A_DEB_FLYw",
    "B_DEB_FLYe",
    "C_DEB_CWYe",
    "D_DEB_CWKe",
    "E_DEB_CLYe",
    "F_DEB_CLKe",
    "G_DEB_FLKe",
    "H_DEB_FWKe"
  )
  
  l_RLWType <- c(
    "w",
    "e",
    "e",
    "e",
    "e",
    "e",
    "e",
    "e"
  )
  
  l_fType <- c(
    "Y",
    "Y",
    "Y",
    "K",
    "Y",
    "K",
    "K",
    "K"
  )
  
  l_KapChange <- c(
    "F",
    "F",
    "C",
    "C",
    "C",
    "C",
    "F",
    "F"
  )
  
  l_FmType <- c(
    "L",
    "L",
    "W",
    "W",
    "L",
    "L",
    "L",
    "W"
  )
  
  df_DEBmodels <- data.frame(
    Name = l_DEB_models,
    KapChange = l_KapChange,
    FmType = l_FmType,
    fType = l_fType,
    RLWType = l_RLWType
  )
  
  Text_Start <- "# DEB Apporectodea caliginosa
# -----------------------------------------------------------------------

States = {E,         # Energy in the reserve (J)
          L,         # Structural volumetric length (cm)
          Eh,        # Maturity (J)
          R,         # Cumulated reproduction (# of cocoons produced)
          OM_soil,   # OM soil quantity in cosm (for Dens earthworms) (g)
          OM_horse,  # OM horse dung quantity in cosm (for Dens earthworms) (g)
          z};        # Threshold for kappa change (-)


Outputs = {Energy,             # Energy in the reserve (J)
          Size,                # Size of the earthworm (cm)
          Length_struct,       # Structural volumetric length (cm)
          Maturity,            # Maturity (J)
          Reproduction,        # Cumulated reproduction (# of cocoons produced)
		      Organic_matter,
		      Weight,
		      kapreal,
          freal,
          Q_soil,
          Q_horse,
          Qf_soil,
          Qf_horse,
          Q_mattermax,
          Q_matter,
          OM_soil_out,
          OM_horse_out,
          OMsoil,
          OMhorse,
          X};             # Wet weight of the earthworm (g)
		   
Inputs = {Dens, # Number of earthworm in the cosm (#)
          Init,
          Soil,
          Winit,
          Texp,
          WeightSoilCosm,
          OM_horse_t0,
          OM_soil_t0,
          times_replace,
          OM_horse_replace,
          OM_horse_add,
          OM_soil_replace};              

# Error parameters
Sigma_W; # Error on Weight
Sigma_R; # Error on Reproduction


# PARAMETERS ---------------------------------------------------------------

# Environment
muOM         = 13000;            # Energy in 1 g of organic matter from horse dung (J/g) 0.9 * 13000
rOM_ClxHorse = 0.18;             # Ratio of energy in 1 g of OM from soil vs 1 g of OM from horse dung

# Energy assimilation
Fm           = 0.5;                # [F_m], max spec searching rate (g/d.g)
kapX         = 0.258;              # digestion efficiency of food to reserve (-)
pAm          = 1712.41;            # {p_Am}, spec assimilation flux (J/d.cm^2)
r_pAm_pM     = 1.5;                # Ratio between pAm and pM (pM = pAm * r_pAm_pM)
v            = 0.18505;            # energy conductance (cm/d)

K            = 0.1;                # Slope of functional threshold (K model)

kap          = 0.4373;             # allocation fraction to soma (-)
kapJ;                              # allocation fraction to soma for juveniles (-) (model A)
kapA;                              # allocation fraction to soma for grouped adults (-) (model A)
tau;                               # transition time between kapJ and kapA (d) (model A)

# Maintenance and growth
pM           = 1680.1;               # [p_M], vol-spec somatic maint (J/d.cm^3)
pT           = 0;                    # {p_T}, surf-spec somatic maint (J/d.cm^2)
Eg           = 4183;                 # 6880.33;     # [E_G], spec cost for structure (J/cm^3)
Shape        = 0.066193;             # shape coefficient (-)
dv           = 1;                    # Structure density (g/cm^3) (models w & e)
wE;                                  # Energy density of reserve (g/J) (model w)
de           = 27.24;                # Contribution of reserve to body weight (g/cm^3) (model e)
wv           = 6.72;                 # Contribution of water replacement to body fresh weight (g/cm^3) (model w)

# Maturity
kJ           = 0.002793;          # maturity maint rate coefficient (1/d)
Ehb          = 0.5;               # maturity at birth (J)
Ehp          = 100;               # maturity at puberty (J)

# Reproduction
L_coc        = 0.23;       # Structural length of a cocoon (cm)
E_coc        = 467;        # Energy in a cocoon (J)
UE0          = 0.078;      # Scaled cost of an egg (cm^2/d)
kapR         = 0.475;      # reproduction efficiency (-)

# Metabolic response to T
TAH          = 28750;       # Arrhenius temperature for upper boundary (K)
TH           = 293.2;       # upper boundary (K)
TA           = 7976;        # Arrhenius temperature (K)
Tref         = 293.15;      # Reference temperature (K)


Init_R = 1e-6;
"
  
  Text_End <- "CalcOutputs{
    Energy         = E;
    Size           = L*1/Shape;
    Length_struct  = L;
    Maturity       = (Eh <= 0 ? 1E-9 : Eh); 
    Reproduction   = (R < (Init_R*10) ? 1e-6 : R);
    Organic_matter = OM_soil + OM_horse;
    Weight         = Weight;
    Q_soil         = Q_soil;
    Q_horse        = Q_horse;
    Qf_soil        = Qf_soil; 
    Qf_horse       = Qf_horse;
    OMsoil         = OM_soil;
    OMhorse        = OM_horse;
    OM_soil_out    = OM_soil_out;
    OM_horse_out   = OM_horse_out;
    Q_mattermax    = Q_mattermax;
    Q_matter       = Q_matter;
    X              = X;
    kapreal        = kapreal;
    freal          = freal;
}

End."
  
  for (i in 1:length(df_DEBmodels$Name)) {
    
    Name_i <- df_DEBmodels$Name[i]
    KapChange_i <- df_DEBmodels$KapChange[i]
    FmType_i <- df_DEBmodels$FmType[i]
    fType_i <- df_DEBmodels$fType[i]
    RLWType_i <- df_DEBmodels$RLWType[i]
    
    if (KapChange_i == "F") {
      Text_Kappa <- "      kapreal = kap; # Model F"
      Text_Kappa_Dyn <- "      kapreal = kap; # Model F
            z = 0;"
    } else if (KapChange_i == "C") {
      Text_Kappa <- "      kapreal = (Eh < Ehp ? kapJ : (Dens > 1 ? kapA : kapJ)); # Model A"
      Text_Kappa_Dyn <- "      dt(z) = (Eh >= Ehp ? 1 : 0);                                 #(Model A)
            x = (z/tau > 1 ? 1 : z/tau);                                 # Model A
            kapAreal = x * kapA + (1-x) * kapJ;                          # Model A
            kapreal  = (Eh < Ehp ? kapJ : (Dens > 1 ? kapAreal : kapJ)); # Model A"
    }
    
    if (RLWType_i == "w") {
      Text_Linit <-  "        L = pow(((Winit)/(dv+wE)),(1.0/3.0)); # Model w"
      Text_Weight <- "      Weight = pow(L,3)*(dv + E/Em*wE); # Model w"
    } else if (RLWType_i == "e") {
      Text_Linit <- "        L = pow(((Winit-de*Em)/(dv)),(1.0/3.0)); # Model e"
      Text_Weight <- "      Weight = pow(L,3)*dv + E*de; # Model e"
    }
    
    if (fType_i == "K") {
      Text_f <- "     # Functional response of the earthworm to organic matter - K model
          X = Qf/(Qfmax);
          freal = X/(K*(1-X)+X);"
    } else if (fType_i == "Y") {
      Text_f <- "     # Functional response of the earthworm to organic matter - Y model
           X = Qf/(2*Qfmax);
          freal = X/(1+X);"
    }
    
    if (FmType_i == "W") {
      Text_Fm <- "      Type_Fm = Weight_struct; # Model W"
    } else if (FmType_i == "L") {
      Text_Fm <- "      Type_Fm = pow(L,2); # Model L"
    }

    
    Text_Init <- paste(
      "Initialize{
        
        ",
      Text_Kappa,
        "
        
        pM = r_pAm_pM * pAm;
        
        # Coefficient for temperature correction (Texp in °C)
        sA  = exp(TA/Tref-TA/(Texp+273.15));
        srH = (1 + exp(TAH/TH-TAH/Tref))/(1+exp(TAH/TH-TAH/(Texp+273.15)));
        Tc  = sA*((Texp+273.15>=Tref)*srH + (Texp+273.15<Tref));
        
        # Correction of coefficients depending on temperature
        pAm_t  = pAm * Tc;        
        v_t    = v   * Tc;
        kJ_t   = kJ * Tc;
        
        # Primary parameters          
        Km   = pM / Eg ;
        Km_t = Km * Tc;
        Em   = pAm / v ;
        g    = Eg / (kapreal * Em);
        Lm   = (v / (Km*g));
        
        # Starting points
        E = Em;",
        Text_Linit,
        "
        Weight = Winit;
        Weight_struct = dv*pow(L,3);
        Eh = Ehb;
        R = Init_R;
        freal = 1;
        z = 0;
        
        OM_soil_out = 0;
        OM_horse_out = 0;
        
      } # End of Initialize
      ",
      sep = "\n"
    )
    
    Text_Dyn <- paste(
      "Dynamics{

    # Coefficient for temperature correction (Texp in °C)
    	sA  = exp(TA/Tref-TA/(Texp+273.15));
      srH = (1 + exp(TAH/TH-TAH/Tref))/(1+exp(TAH/TH-TAH/(Texp+273.15)));
    	Tc  = sA*((Texp+273.15>=Tref)*srH + (Texp+273.15<Tref));
    	
    # Correction of coefficients depending on temperature
      pAm_t  = pAm * Tc;        
      v_t    = v   * Tc;
      kJ_t   = kJ  * Tc;
     ",
      Text_Kappa_Dyn,
      "

      pM = r_pAm_pM * pAm;
      
    # Primary parameters          
      Km   = pM / Eg ;
      Km_t = Km * Tc;
      Em   = pAm / v ;
      g    = Eg / (kapreal * Em);
      Lm   = (v / (Km*g));
      
    
    # Energie assimilation
      
      W_cosm = WeightSoilCosm + OM_horse/0.9;                        
      
      Weight_struct = dv*pow(L,3);
      ",
      Text_Fm,
      "
      
      Qfmax = pAm_t * pow(L,2)/kapX; # (J/d)
      
      Q_mattermax = (Qfmax * W_cosm) / (OM_horse * muOM + OM_soil * muOM * rOM_ClxHorse); # (g matter/d)
      
      # Quantity of matter ingested per earthworm per day (g/d)
      Q_matter_tmp = Fm * Type_Fm ;
      Q_matter = (Q_matter_tmp > Q_mattermax ? Q_mattermax : Q_matter_tmp);              # (g matter/d)
      
      Q_soil_tmp  = Q_matter * (OM_soil /W_cosm);
      Q_soil = ( (Q_soil_tmp * Dens) > OM_soil ? (OM_soil/Dens) : Q_soil_tmp);
      Q_horse_tmp = Q_matter * (OM_horse/W_cosm);
      Q_horse = ( (Q_horse_tmp * Dens) > OM_horse ? (OM_horse/Dens) : Q_horse_tmp);
      
      # Quantity of ingested energy per earthworm per day (J/d)
      Qf_soil  = Q_soil * muOM * rOM_ClxHorse;
      Qf_horse = Q_horse * muOM;
      Qf = Qf_soil + Qf_horse;
      
      # Quantity of organic matter consummed by all earthworms in one cosm per day (g/d)
      OM_soil_out  = Q_soil  * Dens; # gOM/cosm/d  
      OM_horse_out = Q_horse * Dens; # gOM/cosm/d   
      
      # Dynamics of organic matter in the cosm (which cannot be negative)
    
      dt(OM_soil)  = - OM_soil_out;
      dt(OM_horse) = - OM_horse_out;
      
      OM_soil  = (OM_soil  < 1e-9 ? 1e-9 : OM_soil);
      OM_horse = (OM_horse < 1e-9 ? 1e-9 : OM_horse);
      ",
      Text_f,
      "
      
    # Energy in the reserve (which cannot be negative)
    
      tmp_E = (pAm_t/L)*(freal-E/Em);
      dt(E) = (E < 1e-12 ? (-E) : tmp_E);

    # Structural length
    
      dt(L) = (v_t / (3 * ( ( (E/Em) + g) ) ) ) * ((E/Em) - (L/Lm));
      ",
      Text_Weight,
      "

    # Maturity
    
      pC       = ( (g * E ) / (g + (E/Em)) ) * ((v_t*pow(L,2)) + (Km_t*pow(L,3)));
      tmp_Eh   = ((1-kapreal) * pC) - (kJ_t* Eh);
      dt(Eh)   = (Eh < Ehp ? tmp_Eh : 0);

    # Reproduction
    
      tmp_R    = kapR * ( ((1-kapreal)*pC) - (kJ_t*Ehp));
  	  dt(R)    = (Eh >= Ehp ? tmp_R/(E_coc) : 0);

} # End of Dynamics
",
      sep = "\n"
    )
    
    
    Text_full_i <- paste(
      Text_Start,
      Text_Init,
      Text_Dyn,
      Text_End, 
      sep="\n"
      )
    
    path_model_i <- here::here("mod", paste0(Name_i, "/DEBcali.model"))
    writeLines(
      Text_full_i,
      path_model_i
    )
    
  }
  
}

f_In_experiments <- function(l_Experiments, Adult_alone){
  
  df_DEB_mean <- f_import_data_DEB(l_Experiments, Adults_alone) |> 
    mutate(
      OM_soil = signif(OM_soil, 3),
      OM_horse = signif(OM_horse, 3)
    )  |> 
    mutate(across(where(is.numeric), ~ replace_na(.x, -1))) 
  
  char_final <- ""
  
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
       arrange(Time) #|>                    
      # filter(row_number() == 1 | 
      #          Reproduction != lag(Reproduction) |
      #          row_number() == n())
    
    char_tR <- df_R_i |> pull(Time)         |>  paste(collapse = ", ") # Time in days
    char_R  <- df_R_i |> pull(Reproduction) |>  paste(collapse = ", ") # Cumulated reproduction (#)
    
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
    
    text_bol_XP3 <- paste("    bol_XP3=", bol_XP3, ";")
    
    text_OM_final <- paste(
      text_OM_soil_Replace,
      text_OM_horse_Replace,
      text_OM_horse_Add,
     # text_bol_XP3,
      sep="\n"
    )
    
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

f_In_tot <- function(file_path, l_Experiments, Adults_alone = FALSE, text_priors, text_likelihood, NbIter, seeds, RTOL, ATOL) {
  
  text_Level_global <- 
    'Level{ # Global
    '
  
  text_Level_exp <- '
Level{

  ############## Experiments ###################
  '
  
  text_experiment <- f_In_experiments(l_Experiments, Adult_alone)
  
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

SetPoints("DEB_setpoint_predobs_', compteur_exp, '.out", "tab_setpoint.out", 0,', l_param_name|>  paste(collapse = ", "),');

Integrate(Lsodes, 1E-4, 1E-5, 1);

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
    char_tR <- df_R_i |> pull(Time)         |>  paste(collapse = ", ") # Time in days
    
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
      mutate(No_sim = j) |> 
      mutate(across(everything(), as.numeric))
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

Integrate(Lsodes, 1E-4, 1E-5, 1);

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

f_read_Setpoint_full <- function(file_path, l_Experiments, Adults_alone = FALSE, Nb_experiment, select_times, step_sim, bol_Eh = FALSE, bol_OM = FALSE){
  
  df_pred_full <- NULL
  
  df_data_obs <- f_import_data_DEB(l_Experiments, Adults_alone)
  
  i <- 0
  
  for (ID_i in unique(df_data_obs$ID_experiment)){
    
    i <- i + 1
    
    Sim.Res.Exp_i <- readr::read_tsv(
      file.path(file_path, paste0("DEB_setpoint_full_", i, ".out"))
    ) |>
      dplyr::mutate(No_sim = i) |> 
      dplyr::select(where(~ any(!is.na(.x))))
    
    if (bol_Eh) {
      Nom_Endpoints = c("Weight", "Maturity", "Energy")
      Nsortie = 3
    } else if (bol_OM) {
      Nom_Endpoints = c("freal",
                        "OMhorse",
                        "OMsoil",
                        "OM_horse_out",
                        "OM_soil_out",
                        "Q_horse",
                        "Q_soil",
                        "Qf_horse",
                        "Qf_soil",
                        "Q_matter",
                        "Q_mattermax",
                        "Weight")
      Nsortie = 12

    } else {
      Nom_Endpoints = c("Weight", "Reproduction")
      Nsortie = 2
    }
  
    df_data_obs_i <- subset(df_data_obs, ID_experiment == ID_i)
    
    Time_pred <- seq(0, max(df_data_obs_i$Time, na.rm=TRUE), step_sim) # !!!

    
    MPV = IC_min = IC_max = Endpoint = index = NULL
    
    for (j in 1:length(Nom_Endpoints)){
      
      index = grep(paste0("^",Nom_Endpoints[j],"_") , colnames(Sim.Res.Exp_i), fixed=FALSE)
      
      
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


f_AS_Sobol_design <- function(Sobol_type, seed, Nb_sim, CV, Processors, Nboot, l_params_median, Temperature_bol){
  
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

# f_AS_Sobol_index <- function(df_AS_sim, Sobol_info){
#   
#   df_AS_sim_mean <- df_AS_sim |> 
#     summarise(
#       across(
#         .cols = is.numeric, 
#         .fns = list(Mean = mean),
#         na.rm = TRUE,
#         .names = "{col}"
#       )
#     )
#   
#   df_AS_sim <- df_AS_sim |>
#     mutate(
#       Weight_Juv = as.numeric(Weight_Juv),
#       across(
#         c(
#           Weight_Juv, Weight_Adult,
#           Energy_Juv, Energy_Adult,
#           Maturity_Juv, Maturity_Adult,
#           Reproduction, Organic_matter
#         ),
#         ~ coalesce(.x, df_AS_sim_mean[[cur_column()]][1])
#       )
#     )
#   
#   # Sobol calculations
#   
#   Nb_out <- length(df_AS_sim)
#   df_res <- data.frame()  # liste vide
#   
#   for (i in seq_along(df_AS_sim)) {
#     
#     Y <- dplyr::pull(df_AS_sim, i)
#     tell(x = Sobol_info, y = Y, nboot =  Nboot, conf = 0.95)
#     
#     df_res_i <- data.frame(
#       Variable = colnames(df_AS_sim)[i],
#       FOI = Sobol_info$S[,1],
#       FOI.borninf = Sobol_info$S[,4],
#       FOI.bornsup = Sobol_info$S[,5],
#       TI = Sobol_info$T[,1],
#       TI.borninf = Sobol_info$T[,4],
#       TI.bornsup = Sobol_info$T[,5]
#     )
#     
#     df_res <- rbind(df_res, df_res_i)
#   }
#   
#   df_res <- df_res %>%
#     mutate(Parameter = rep(rownames(Sobol_info$S), Nb_out)) %>%
#     arrange(desc(TI))
#   
#   return(df_res)
# }

f_AS_Sobol_index <- function(
    df_AS_sim,
    Sobol_info,
    l_sim_fail = integer()
) {
  
  # Dimensions du plan Jansen
  n <- nrow(Sobol_info$X1)
  p <- ncol(Sobol_info$X1)
  
  n_expected <- n * (p + 2L)
  
  # Numéros des simulations échouées
  l_sim_fail <- sort(unique(as.integer(l_sim_fail)))
  
  if (anyNA(l_sim_fail)) {
    stop("l_sim_fail contains NA(s).")
  }
  
  if (any(l_sim_fail < 1L | l_sim_fail > n_expected)) {
    stop(
      "Numbers of l_sim_fail must be in the range 1 to ",
      n_expected, "."
    )
  }
  
  if (nrow(df_AS_sim) + length(l_sim_fail) != n_expected) {
    stop(
      "Incoherence between outputs and l_sim_fail :\n",
      "- Outputs expected : ", n_expected, "\n",
      "- Outputs present : ", nrow(df_AS_sim), "\n",
      "- Declared failed simulations : ", length(l_sim_fail)
    )
  }
  
  # Conversion des sorties en numérique
  df_AS_sim <- df_AS_sim |>
    mutate(
      across(
        everything(),
        ~ suppressWarnings(as.numeric(.x))
      )
    )
  
  # Reconstruction du tableau complet avec NA aux positions échouées
  Y_full <- matrix(
    NA_real_,
    nrow = n_expected,
    ncol = ncol(df_AS_sim),
    dimnames = list(NULL, names(df_AS_sim))
  )
  
  positions_ok <- setdiff(
    seq_len(n_expected),
    l_sim_fail
  )
  
  Y_full[positions_ok, ] <- as.matrix(df_AS_sim)
  
  # Ligne de base correspondant à chaque simulation échouée
  # Exemple avec n = 30000 :
  # 1, 30001, 60001... correspondent à la même ligne de base
  base_rows_fail <- if (length(l_sim_fail) > 0L) {
    sort(unique((l_sim_fail - 1L) %% n + 1L))
  } else {
    integer()
  }
  
  base_rows_keep <- setdiff(
    seq_len(n),
    base_rows_fail
  )
  
  # Positions à conserver dans chacun des p + 2 blocs
  positions_keep <- unlist(
    lapply(
      0:(p + 1L),
      function(block) block * n + base_rows_keep
    ),
    use.names = FALSE
  )
  
  Y_reduced <- Y_full[
    positions_keep,
    ,
    drop = FALSE
  ]
  
  # Éventuels NA restants dans des sorties présentes
  # Ils sont remplacés par la moyenne de la sortie concernée
  for (j in seq_len(ncol(Y_reduced))) {
    
    bad <- !is.finite(Y_reduced[, j])
    
    if (all(bad)) {
      stop(
        "No finite value for the output ",
        colnames(Y_reduced)[j], "."
      )
    }
    
    if (any(bad)) {
      Y_reduced[bad, j] <- mean(
        Y_reduced[!bad, j],
        na.rm = TRUE
      )
    }
  }
  
  # Reconstruction du plan Sobol avec les lignes complètes uniquement
  Sobol_reduced <- sensitivity::soboljansen(
    model = NULL,
    X1 = Sobol_info$X1[
      base_rows_keep,
      ,
      drop = FALSE
    ],
    X2 = Sobol_info$X2[
      base_rows_keep,
      ,
      drop = FALSE
    ],
    nboot = Sobol_info$nboot,
    conf = Sobol_info$conf
  )
  
  message(
    length(l_sim_fail), " Failed simulations ; ",
    length(base_rows_fail), " Sobol lines deleted ; ",
    length(base_rows_keep), " lines kept out of ", n, "."
  )
  
  # Calcul des indices
  Nb_out <- ncol(Y_reduced)
  l_res <- vector("list", Nb_out)
  
  for (i in seq_len(Nb_out)) {
    
    Y <- Y_reduced[, i]
    
    # Copie indépendante pour chaque sortie
    Sobol_i <- Sobol_reduced
    
    tell(
      Sobol_i,
      y = Y
    )
    
    l_res[[i]] <- data.frame(
      Variable = colnames(Y_reduced)[i],
      Parameter = rownames(Sobol_i$S),
      
      FOI = Sobol_i$S[, 1],
      FOI.borninf = Sobol_i$S[, 4],
      FOI.bornsup = Sobol_i$S[, 5],
      
      TI = Sobol_i$T[, 1],
      TI.borninf = Sobol_i$T[, 4],
      TI.bornsup = Sobol_i$T[, 5]
    )
  }
  
  res <- dplyr::bind_rows(l_res) |>
    arrange(Variable, desc(TI))
  
  return(res)
}

## DEB-TKTD ----

f_write_DEBTKTD_models <- function(KapChange, FmType, fType, RLWType) {
  
  l_DEB_models <- c(
    "A_EPX_A",
    "B_EPX_M",
    "C_EPX_G",
    "D_EPX_R",
    "E_EPX_AM",
    "F_EPX_AG",
    "G_EPX_AR",
    "H_EPX_MG",
    "I_EPX_MR",
    "J_EPX_GR",
    
    "A_IMD_A",
    "B_IMD_M",
    "C_IMD_G",
    "D_IMD_R",
    "E_IMD_AM",
    "F_IMD_AG",
    "G_IMD_AR",
    "H_IMD_MG",
    "I_IMD_MR",
    "J_IMD_GR"
  )
  
  l_bol_A <- rep(c(
    1,
    0,
    0,
    0,
    1,
    1,
    1,
    0,
    0,
    0
  ),2)
  
  l_bol_M <- rep(c(
    0,
    1,
    0,
    0,
    1,
    0,
    0,
    1,
    1,
    0
  ), 2)
  
  l_bol_G <- rep(c(
    0,
    0,
    1,
    0,
    0,
    1,
    0,
    1,
    0,
    1
  ), 2)
  
  l_bol_R <- rep(c(
    0,
    0,
    0,
    1,
    0,
    0,
    1,
    0,
    1,
    1
  ), 2)
  
  l_Molecule <- c(
    rep("EPX", 10),
    rep("IMD", 10)
  )
  
  df_DEBmodels <- data.frame(
    Name = l_DEB_models,
    Molecule = l_Molecule,
    bol_A = l_bol_A,
    bol_M = l_bol_M,
    bol_G = l_bol_G,
    bol_R = l_bol_R
  )
  
  Text_End <- "CalcOutputs{
    Energy         = E;
    Size           = L*1/Shape;
    Length_struct  = L;
    Maturity       = (Eh <= 0 ? 1E-9 : Eh); 
    Reproduction   = (R < (Init_R*10) ? 1e-6 : R);
    Organic_matter = OM_soil + OM_horse;
    Weight         = Weight;
    Q_soil         = Q_soil;
    Q_horse        = Q_horse;
    Qf_soil        = Qf_soil; 
    Qf_horse       = Qf_horse;
    OMsoil         = OM_soil;
    OMhorse        = OM_horse;
    OM_soil_out    = OM_soil_out;
    OM_horse_out   = OM_horse_out;
    Q_mattermax    = Q_mattermax;
    Q_matter       = Q_matter;
    X              = X;
    kapreal        = kapreal;
    freal          = freal;
    Ci             = Ci;
    Ce             = Cext;
    Stress         = b * ((D - nec) > 0 ? (D - nec) : 0);
    Damage         = D;
}

End."

for (i in 1:length(df_DEBmodels$Name)) {
  Name_i <- df_DEBmodels$Name[i]
  Molecule_name <- df_DEBmodels$Molecule[i]
  bol_A <- df_DEBmodels$bol_A[i]
  bol_M <- df_DEBmodels$bol_M[i]
  bol_G <- df_DEBmodels$bol_G[i]
  bol_R <- df_DEBmodels$bol_R[i]
  
  Text_Start <- paste0("# DEB-TKTD ", Molecule_name, " - Apporectodea caliginosa
# -----------------------------------------------------------------------

States = {E,         # Energy in the reserve (J)
          L,         # Structural volumetric length (cm)
          Eh,        # Maturity (J)
          R,         # Cumulated reproduction (# of cocoons produced)
          OM_soil,   # OM soil quantity in cosm (for Dens earthworms) (g)
          OM_horse,  # OM horse dung quantity in cosm (for Dens earthworms) (g)
          z,
          Ci_cent,
          Ci_periph,
          Cext,
          D};        # Threshold for kappa change (-)


Outputs = {Energy,             # Energy in the reserve (J)
          Size,                # Size of the earthworm (cm)
          Length_struct,       # Structural volumetric length (cm)
          Maturity,            # Maturity (J)
          Reproduction,        # Cumulated reproduction (# of cocoons produced)
		      Organic_matter,
		      Weight,
		      kapreal,
          freal,
          Q_soil,
          Q_horse,
          Qf_soil,
          Qf_horse,
          Q_mattermax,
          Q_matter,
          OM_soil_out,
          OM_horse_out,
          OMsoil,
          OMhorse,
          X,
          Ci,
          Ce,
          Stress,
          Damage};        
		   
Inputs = {Dens, # Number of earthworm in the cosm (#)
          Init,
          Soil,
          Winit,
          Texp,
          WeightSoilCosm,
          OM_horse_t0,
          OM_soil_t0,
          times_replace,
          OM_horse_replace,
          OM_horse_add,
          OM_soil_replace};              

# Error parameters
Sigma_W; # Error on Weight
Sigma_R; # Error on Reproduction


# PARAMETERS ---------------------------------------------------------------

# Environment
muOM         = 13000;            # Energy in 1 g of organic matter from horse dung (J/g) 0.9 * 13000
rOM_ClxHorse = 0.0294325;        # Ratio of energy in 1 g of OM from soil vs 1 g of OM from horse dung

# Energy assimilation
Fm           = 3.87576;            # [F_m], max spec searching rate (g/d.g)
kapX         = 0.258;              # digestion efficiency of food to reserve (-)
pAm          = 341.857;            # {p_Am}, spec assimilation flux (J/d.cm^2)
r_pAm_pM     = 0.70081;            # Ratio between pAm and pM (pM = pAm * r_pAm_pM)
v            = 0.230918;           # energy conductance (cm/d)

K            = 0.0546564;          # Slope of functional threshold (K model)

kap          = 0.4373;             # allocation fraction to soma (-)
kapJ         = 0.983041;           # allocation fraction to soma for juveniles (-) (model A)
kapA         = 0.728915;           # allocation fraction to soma for grouped adults (-) (model A)
tau          = 11.9528;            # transition time between kapJ and kapA (d) (model A)

# Maintenance and growth
pM           = 1680.1;               # [p_M], vol-spec somatic maint (J/d.cm^3)
pT           = 0;                    # {p_T}, surf-spec somatic maint (J/d.cm^2)
Eg           = 4183;                 # [E_G], spec cost for structure (J/cm^3)
Shape        = 0.066193;             # shape coefficient (-)
dv           = 1.2135;               # Structure density (g/cm^3) (models w & e)
wE;                                  # Energy density of reserve (g/J) (model w)
de           = 1.130890e-09;         # Contribution of reserve to body weight (g/cm^3) (model e)
wv           = 6.72;                 # Contribution of water replacement to body fresh weight (g/cm^3) (model w)

# Maturity
kJ           = 0.002793;          # maturity maint rate coefficient (1/d)
Ehb          = 0.5;               # maturity at birth (J)
Ehp          = 73.1002;           # maturity at puberty (J)

# Reproduction
L_coc        = 0.23;       # Structural length of a cocoon (cm)
E_coc        = 74.966;    # Energy in a cocoon (J)
UE0          = 0.078;      # Scaled cost of an egg (cm^2/d)
kapR         = 0.475;      # reproduction efficiency (-)
kapH         = 1;          # Maturity efficiency (-)

# Metabolic response to T
TAH          = 28750;       # Arrhenius temperature for upper boundary (K)
TH           = 293.2;       # upper boundary (K)
TA           = 7976;        # Arrhenius temperature (K)
Tref         = 293.15;      # Reference temperature (K)


Init_R = 1e-6;
")
  
if (Molecule_name == "EPX") {
  
  Text_Param_TKTD <- "
# DEB-TKTD PARAMETERS ----------------------------------------------------------

# Initialisation of the environment
Ci0 = 5.84; # Initial concentration of toxicant in earthworm (ng/g) (5.84 for EPX | 6.2 for IMD)
Ce0; # Initial concentration of toxicant in soil (ng/g)

# Kinetics in earthworm
ku = 1.5;         # Uptake rate of toxicant (g(soil).g(worm)-1.d-1) | 0.851 x 1.8 (bc here earthworm weight with gut content)
ke = 1.8;         # Elimination rate of toxicant (d-1) 
kperiph = 0.0001; # Second compartment dynamics (d-1)

# Kinetics in soil 
DT50 = 353.5;      # EPX Half-life in soil (days)

# Effects
kr = 1;  # Damage reparation constant (d-1)
nec; # No-effect concentration for sublethal effects (ng/g)
b;   # Tolerance concentration (ng/g)
  "
  
  Text_TK <- "
      dt(Ci_periph) = (kperiph*Ci_cent - kperiph*Ci_periph) - (3/L)*dt(L)*Ci_periph;
      dt(Ci_cent) = (ku*Cext - ke*Ci_cent + kperiph*Ci_periph - kperiph*Ci_cent) - (3/L)*dt(L)*Ci_cent; # dt(L3) = 3 L2 dt(L)
      Ci = Ci_periph + Ci_cent;
  "
  
} else if (Molecule_name == "IMD") {
  
  Text_Param_TKTD <- "
# DEB-TKTD PARAMETERS ----------------------------------------------------------

# Initialisation of the environment
Ci0 = 6.2; # Initial concentration of toxicant in earthworm (ng/g) (5.84 for EPX | 6.2 for IMD)
Ce0; # Initial concentration of toxicant in soil (ng/g)

# Kinetics in earthworm
ku = 2.0;         # Uptake rate of toxicant (g(soil).g(worm)-1.d-1) | 0.851 x 2.4 (bc here earthworm weight with gut content)
ke = 0.12;         # Elimination rate of toxicant (d-1) 
kperiph = 0.043; # Second compartment dynamics (d-1)
C_sat = 1676; # Internal saturation concentration (ng/g)

# Kinetics in soil 
DT50 = 187;      # IMD Half-life in soil (days)

# Effects
kr = 1;  # Damage reparation constant (d-1)
nec; # No-effect concentration for sublethal effects (ng/g)
b;   # Tolerance concentration (ng/g)"
  
  Text_TK <- "
      tmp_Ci_periph = (kperiph*Ci_cent - kperiph*Ci_periph) - (3/L)*dt(L)*Ci_periph;
      tmp_Ci_cent = (ku*Cext - ke*Ci_cent + kperiph*Ci_periph - kperiph*Ci_cent) - (3/L)*dt(L)*Ci_cent; # dt(L3) = 3 L2 dt(L)
      
      dt(Ci_periph) = (Ci > C_sat ? 0 : tmp_Ci_periph);
      dt(Ci_cent) = (Ci > C_sat ? 0 : tmp_Ci_cent);
      
      Ci = ((Ci_periph + Ci_cent) > C_sat ? C_sat : (Ci_periph + Ci_cent));
  "
  
}
  
  Text_pMoA <- paste(
    "# pMoA taken into account",
    paste0("bol_A = ", bol_A,";"),
    paste0("bol_M = ", bol_M,";"),
    paste0("bol_G = ", bol_G,";"),
    paste0("bol_R = ", bol_R,";
           "),
    sep = "\n"
  )
  
  if (KapChange == "F") {
    Text_Kappa <- "      kapreal = kap; # Model F"
    Text_Kappa_Dyn <- "      kapreal = kap; # Model F
            z = 0;"
  } else if (KapChange == "C") {
    Text_Kappa <- "      kapreal = (Init == 2 ? kapA : (Eh < Ehp ? kapJ : (Dens > 1 ? kapA : kapJ))); # Model A"
    Text_Kappa_Dyn <- "      dt(z) = (Eh >= Ehp ? 1 : 0);                                # Model A
      x = (z/tau > 1 ? 1 : z/tau);                                # Model A
      kapAreal = x * kapA + (1-x) * kapJ;                         # Model A
      kapreal = (Init == 2 ? kapA : (Eh < Ehp ? kapJ : (Dens > 1 ? kapAreal : kapJ))); # Model A"
  }
  
  if (RLWType == "w") {
    Text_Linit <-  "        L = pow(((Winit)/(dv+wE)),(1.0/3.0)); # Model w"
    Text_Weight <- "      Weight = pow(L,3)*(dv + E/Em*wE); # Model w"
  } else if (RLWType == "e") {
    Text_Linit <- "        L = pow(((Winit-de*Em)/(dv)),(1.0/3.0)); # Model e"
    Text_Weight <- "      Weight = pow(L,3)*dv + E*de; # Model e"
  }
  
  if (fType == "K") {
    Text_f <- "     # Functional response of the earthworm to organic matter - K model
          X = Qf/(Qfmax);
          freal = X/(K*(1-X)+X) * ((1-s_f*bol_A) > 0 ? (1-s_f*bol_A) : 0);"
  } else if (fType == "Y") {
    Text_f <- "     # Functional response of the earthworm to organic matter - Y model
           X = Qf/(2*Qfmax);
          freal = X/(1+X) * ((1-s_f*bol_A) > 0 ? (1-s_f*bol_A) : 0);"
  }
  
  if (FmType == "W") {
    Text_Fm <- "      Type_Fm = Weight_struct; # Model W"
  } else if (FmType == "L") {
    Text_Fm <- "      Type_Fm = pow(L,2); # Model L"
  }
  
  
  Text_Init <- paste(
    "Initialize{
    
        Weight = Winit;
        Weight_struct = dv*pow(L,3);
        
        Ehjuv = 58.501 * pow(Weight,2) + 80.538 * Weight - 0.514; # Relationship determined from the predictions of Gollot2026_dens_D5
      
        Eh = (Init == 0 ? Ehb :
              (Init == 1 ? Ehjuv :
              (Init == 2 ? Ehp : Ehb)));
        
        ",
    Text_Kappa,
    "
        
        pM = r_pAm_pM * pAm;
        
        # Coefficient for temperature correction (Texp in °C)
        sA  = exp(TA/Tref-TA/(Texp+273.15));
        srH = (1 + exp(TAH/TH-TAH/Tref))/(1+exp(TAH/TH-TAH/(Texp+273.15)));
        Tc  = sA*((Texp+273.15>=Tref)*srH + (Texp+273.15<Tref));
        
        # Correction of coefficients depending on temperature
        pAm_t  = pAm * Tc;        
        v_t    = v   * Tc;
        kJ_t   = kJ * Tc;
        
        # Primary parameters          
        Km   = pM / Eg ;
        Km_t = Km * Tc;
        Em   = pAm / v ;
        g    = Eg / (kapreal * Em);
        Lm   = (v / (Km*g));
        
        # Starting points
        E = Em;",
    Text_Linit,
    "
        R = Init_R;
        freal = 1;
        z = 0;
        
        OM_soil_out = 0;
        OM_horse_out = 0;
        
        D = 0;
        Ci_cent = Ci0 * kperiph/(kperiph + kperiph);
        Ci_periph = Ci0 - Ci_cent;
        Cext = Ce0;
        
      } # End of Initialize
      ",
    sep = "\n"
  )
  
  Text_Dyn <- paste(
    "Dynamics{
    
    # Toxicokinetics
      dt(Cext) = -log(2)/DT50 * Cext;
      
      ",
    Text_TK,
      "
      dt(D) = kr*(Ci-D); # @jager_2020
      
    # Toxicodynamics
    
      s_f = b * ((D - nec) > 0 ? (D - nec) : 0);

    # Coefficient for temperature correction (Texp in °C)
    	sA  = exp(TA/Tref-TA/(Texp+273.15));
      srH = (1 + exp(TAH/TH-TAH/Tref))/(1+exp(TAH/TH-TAH/(Texp+273.15)));
    	Tc  = sA*((Texp+273.15>=Tref)*srH + (Texp+273.15<Tref));
    	
    # Correction of coefficients depending on temperature
      pAm_t  = pAm * Tc;        
      v_t    = v   * Tc;
      kJ_t   = kJ  * Tc * (1+s_f*bol_M);
     ",
    Text_Kappa_Dyn,
    "

      pM = r_pAm_pM * pAm * (1+s_f*bol_M);
      
    # Primary parameters          
      Km   = pM / (Eg * (1+s_f*bol_G)) ;
      Km_t = Km * Tc;
      Em   = pAm / v ;
      g    = (Eg * (1+s_f*bol_G)) / (kapreal * Em);
      Lm   = (v / (Km*g));
      
    
    # Energie assimilation
      
      W_cosm = WeightSoilCosm + OM_horse/0.9;                        
      
      Weight_struct = dv*pow(L,3);
      ",
    Text_Fm,
    "
      
      Qfmax = pAm_t * pow(L,2)/kapX; # (J/d)
      
      Q_mattermax = (Qfmax * W_cosm) / (OM_horse * muOM + OM_soil * muOM * rOM_ClxHorse); # (g matter/d)
      
      # Quantity of matter ingested per earthworm per day (g/d)
      Q_matter_tmp = Fm * Type_Fm ;
      Q_matter = (Q_matter_tmp > Q_mattermax ? Q_mattermax : Q_matter_tmp);              # (g matter/d)
      
      Q_soil_tmp  = Q_matter * (OM_soil /W_cosm);
      Q_soil = ( (Q_soil_tmp * Dens) > OM_soil ? (OM_soil/Dens) : Q_soil_tmp);
      Q_horse_tmp = Q_matter * (OM_horse/W_cosm);
      Q_horse = ( (Q_horse_tmp * Dens) > OM_horse ? (OM_horse/Dens) : Q_horse_tmp);
      
      # Quantity of ingested energy per earthworm per day (J/d)
      Qf_soil  = Q_soil * muOM * rOM_ClxHorse;
      Qf_horse = Q_horse * muOM;
      Qf = Qf_soil + Qf_horse;
      
      # Quantity of organic matter consummed by all earthworms in one cosm per day (g/d)
      OM_soil_out  = Q_soil  * Dens; # gOM/cosm/d  
      OM_horse_out = Q_horse * Dens; # gOM/cosm/d   
      
      # Dynamics of organic matter in the cosm (which cannot be negative)
    
      dt(OM_soil)  = - OM_soil_out;
      dt(OM_horse) = - OM_horse_out;
      
      OM_soil  = (OM_soil  < 1e-9 ? 1e-9 : OM_soil);
      OM_horse = (OM_horse < 1e-9 ? 1e-9 : OM_horse);
      ",
    Text_f,
    "
      
    # Energy in the reserve (which cannot be negative)
    
      tmp_E = (pAm_t/L)*(freal-E/Em);
      dt(E) = (E < 1e-12 ? (-E) : tmp_E);

    # Structural length
    
      dt(L) = (v_t / (3 * ( ( (E/Em) + g) ) ) ) * ((E/Em) - (L/Lm));
      ",
    Text_Weight,
    "

    # Maturity
    
      pC       = ( (g * E ) / (g + (E/Em)) ) * ((v_t*pow(L,2)) + (Km_t*pow(L,3)));
      tmp_Eh   = ((1-kapreal) * pC) - (kJ_t* Eh);
      kapH_real = kapH/(1+s_f*bol_G);
      dt(Eh)   = (Eh < Ehp ? tmp_Eh : 0) * kapH_real;

    # Reproduction
      kapR_real = kapR/(1+s_f*bol_R);
      tmp_R     = kapR_real * ( ((1-kapreal)*pC) - (kJ_t*Ehp));
  	  dt(R)     = (Eh >= Ehp ? tmp_R/(E_coc) : 0);

} # End of Dynamics
",
    sep = "\n"
  )
  
  
  Text_full_i <- paste(
    Text_Start,
    Text_Param_TKTD,
    Text_pMoA,
    Text_Init,
    Text_Dyn,
    Text_End, 
    sep="\n"
  )
  
  path_model_i <- here::here("mod/Z_DEB-TKTD", paste0(Molecule_name, "/",Name_i, "/DEB-", Molecule_name,".model"))
  writeLines(
    Text_full_i,
    path_model_i
  )
  
}

}


f_import_data_DEBTKTD <- function() {
  
  # EC50 Growth
  
  df_data_EC50_growth_raw <- read_excel(here::here("data/Data_expe_raw/Data_EC50.xlsx"), sheet = "Growth") |> 
    mutate(
      Time = t,
      Weight = w/1000, # conversion in g
      Dose = Dose * 1000, # mg/kg to ng/g
    ) |> 
    dplyr::select(-c(t,w, Date, Comments, ID_camp, ID))|> 
    mutate(Weight = ifelse(Weight == Inf, 0, Weight))
  
  df_data_EC50_growth_IMD_raw <- df_data_EC50_growth_raw |> 
    filter(Molec %in% c("IMD", "Ctrl_1", "Ctrl_2")) 
  
  df_data_EC50_growth_IMD <- df_data_EC50_growth_IMD_raw |> 
    group_by(Condition, Time, Dose, Soil_w) |> 
    summarise(
      Weight = mean(Weight, na.rm = TRUE),
      .groups = "drop"
    ) |> 
    arrange(Dose, Time) |> 
    mutate(
      Reproduction = 0,
      ID_experiment = paste0("IMD_EC50_growth_", Condition),
      Environment = case_when(
        Time %in% c(0,28) ~ "Soil_change",
        Time == 14 ~ "Food"
      ),
      OM_soil = 6.52, # 0.0326 gOM/gsoil x 200 gsoil
      OM_horse = 2.7, # 0.9 x 3 g/ind/14d
      Density = 1,
      Texp = 18,
      No_sim = match(Dose, unique(Dose)),
      Molecule = "IMD"
    )
  
  df_data_EC50_growth_EPX_raw <- df_data_EC50_growth_raw |> 
    filter(Molec %in% c("EPX", "Ctrl_1", "Ctrl_2")) 
  
  df_data_EC50_growth_EPX <- df_data_EC50_growth_EPX_raw |> 
    group_by(Condition, Time, Dose, Soil_w) |> 
    summarise(
      Weight = mean(Weight, na.rm = TRUE),
      .groups = "drop"
    ) |> 
    arrange(Dose, Time) |> 
    mutate(
      Reproduction = 0,
      ID_experiment = paste0("EPX_EC50_growth_", Condition),
      Environment = case_when(
        Time %in% c(0,28) ~ "Soil_change",
        Time == 14 ~ "Food"
      ),
      OM_soil = 6.52, # 0.0326 gOM/gsoil x 200 gsoil
      OM_horse = 2.7, # 0.9 x 3 g/ind/14d
      Density = 1,
      Texp = 18,
      No_sim = match(Dose, unique(Dose)),
      Molecule = "EPX"
    )
  
  df_data_EC50_growth <- rbind(df_data_EC50_growth_EPX, df_data_EC50_growth_IMD) |> 
    mutate(Experiment_type = "Growth")
  
  # EC50 Reproduction
  
  df_data_EC50_repro_raw <- read_excel(here::here("data/Data_expe_raw/Data_EC50.xlsx"), sheet = "Repro") |> 
    mutate(
      Time = t,
      Weight = w_tot/Nb_ind/1000, # conversion in g
      Dose = Dose * 1000, # mg/kg to ng/g
      Reproduction = case_when(
        Time == 0 ~ 0,
        Time == 28 ~ Nb_cocoons/2
      )
    ) |> 
    dplyr::select(-c(t,w_tot, Date, Comments, ID_camp, ID_cosm))|> 
    mutate(Weight = ifelse(Weight == Inf, 0, Weight))
  
  df_time14 <- df_data_EC50_repro_raw |> 
    distinct(Molec, Dose, Condition) |> 
    mutate(
      Time = 14,
      Nb_ind = 2,
      Status = "A",
      Weight = NA_real_,
      Reproduction = NA_real_,
      Soil_w = "200"
    )
  
  df_data_EC50_repro_raw <- bind_rows(
    df_data_EC50_repro_raw,
    df_time14
  ) |> 
    arrange(Dose, Condition, Time)
  
  df_data_EC50_repro_IMD_raw <- df_data_EC50_repro_raw |> 
    filter(Molec %in% c("IMD", "Ctrl_1", "Ctrl_2"))
  
  df_data_EC50_repro_IMD <- df_data_EC50_repro_IMD_raw |> 
    group_by(Condition, Time, Dose, Soil_w) |> 
    summarise(
      Weight = mean(Weight, na.rm = TRUE),
      Reproduction = mean(Reproduction, na.rm = TRUE),
      .groups = "drop"
    ) |> 
    arrange(Dose, Time) |> 
    mutate(
      ID_experiment = paste0("IMD_EC50_repro_", Condition),
      Environment = case_when(
        Time %in% c(0,28) ~ "Soil_change",
        Time == 14 ~ "Food"
      ),
      OM_soil = 6.52, # 0.0326 gOM/gsoil x 200 gsoil
      OM_horse = 5.4, # 0.9 x 3 g/14d x 2 ind
      Density = 2,
      Texp = 18,
      No_sim = match(Dose, unique(Dose)),
      Molecule = "IMD",
      Soil_w = 200
    )
  
  df_data_EC50_repro_EPX_raw <- df_data_EC50_repro_raw |> 
    filter(Molec %in% c("EPX", "Ctrl_1", "Ctrl_2"))
  
  df_data_EC50_repro_EPX <- df_data_EC50_repro_EPX_raw |> 
    group_by(Condition, Time, Dose, Soil_w) |> 
    summarise(
      Weight = mean(Weight, na.rm = TRUE),
      Reproduction = mean(Reproduction, na.rm = TRUE),
      .groups = "drop"
    ) |> 
    arrange(Dose, Time) |> 
    mutate(
      ID_experiment = paste0("EPX_EC50_repro_", Condition),
      Environment = case_when(
        Time %in% c(0,28) ~ "Soil_change",
        Time == 14 ~ "Food"
      ),
      OM_soil = 6.52, # 0.0326 gOM/gsoil x 200 gsoil
      OM_horse = 5.4, # 0.9 x 3 g/14d x 2 ind
      Density = 2,
      Texp = 18,
      No_sim = match(Dose, unique(Dose)),
      Molecule = "EPX",
      Soil_w = 200
    )
  
  df_data_EC50_repro <- rbind(df_data_EC50_repro_EPX, df_data_EC50_repro_IMD) |> 
    mutate(Experiment_type = "Repro")
  
  
  df_results <- rbind(df_data_EC50_growth, df_data_EC50_repro)
  
  return(df_results)
  
}

f_import_data_DEBTKTD_MIX <- function() {
  
  # EC50 Growth
  
  df_data_MIX_growth_raw <- read_excel(here::here("data/Data_expe_raw/Data_Mix.xlsx"), sheet = "ReproSpring2025weight") |> 
    mutate(
      Time = t,
      Weight = as.numeric(w)/1000, # conversion in g
      Dose_IMD = Dose_IMD * 1000, # mg/kg to ng/g
      Dose_EPX = Dose_EPX * 1000 # mg/kg to ng/g
    ) |> 
    dplyr::select(-c(t,w, Date, ID_cosm, ID, Color, Lot, Nb_rep, No_vdt, Box, ID_camp))|> 
    mutate(Weight = ifelse(Weight == Inf, 0, Weight))
  
  df_data_MIX_growth <- df_data_MIX_growth_raw |> 
    group_by(Ratio, Line, Time, Dose_IMD, Dose_EPX) |>
    summarise(
      Weight = mean(Weight, na.rm = TRUE),
      .groups = "drop"
    ) |> 
    arrange(Ratio, Line, Time) |> 
    mutate(
      Condition = paste0(Ratio, Line),
      ID_experiment = paste0("MIX_", Condition),
      Environment = case_when(
        Time %in% c(0,28) ~ "Soil_change",
        Time == 14 ~ "Food"
      ),
      OM_soil = 6.52, # 0.0326 gOM/gsoil x 200 gsoil
      OM_horse = 2.7*2, # 0.9 x 3 g/ind/14d
      Density = 2,
      Texp = 18,
      Soil_w = 200,
      No_sim = match(Condition, Condition),
      Reproduction = 0
    )
  
  # MIX Reproduction
  
  df_data_MIX_repro_raw <- read_excel(here::here("data/Data_expe_raw/Data_Mix.xlsx"), sheet = "ReproSpring2025coc") |> 
    mutate(
      Time = t,
      Reproduction_28 = case_when(
        Time == 0 ~ 0,
        Time == 28 ~ Nb_cocoons/2
      ),
      Condition = paste0(Ratio, Line)
    ) |> 
    dplyr::select(Time, Reproduction_28, Condition) 
  
  df_data_MIX_repro <- df_data_MIX_repro_raw |> 
    group_by(Condition, Time) |> 
    summarise(
      Reproduction_28 = mean(Reproduction_28, na.rm=TRUE),
      .groups = "drop"
    )
  
  df_time14 <- df_data_MIX_growth |> 
    distinct(Condition) |> 
    mutate(
      Time = 14,
      Nb_ind = 2,
      Status = "A",
      Weight = NA_real_,
      Reproduction = NA_real_,
      OM_soil = 6.52,
      OM_horse = 5.4,
      Texp = 18,
      Environment = "Food",
      Soil_w = 200, 
      ID_experiment = paste0("MIX_", Condition)
    )
  
  df_tot <- df_data_MIX_growth |> 
    left_join(df_data_MIX_repro) |> 
    mutate(
      Reproduction = case_when(
        Time == 28 ~ Reproduction_28,
        .default = Reproduction
      )
      ) |> 
    dplyr::select(-c(Reproduction_28)) 
  
  df_results <- bind_rows(
      df_tot,
      df_time14
    ) 
  
  return(df_results)
  
}


f_In_experiments_TKTD <- function(Molec){
  
  df_data <- f_import_data_DEBTKTD()
  df_data_EC50 <- df_data |> 
    filter(Molecule == Molec) |> 
    mutate(across(where(is.numeric), ~ replace_na(.x, -1)))
  
  char_final <- ""
  
  for (ID_i in unique(df_data_EC50$ID_experiment)) {
    
    
    df_data_EC50_i <- df_data_EC50 |> filter(ID_experiment == ID_i)
    
    # t0 = début manip 
    # W = W0 manip
    # Eh = Estimation par rapport à W0 (faire graph a partir DEB normal)
    # R = 0
    
    ################### Data on juvenile growth 
    
    # Weight data (remove Weigth(t=0))
    df_Ww_i <- df_data_EC50_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Weight != lag(Weight))
    
    char_tW <- df_Ww_i |> pull(Time)   |>  paste(collapse = ", ") # Time in days
    char_Ww <- df_Ww_i |> pull(Weight) |>  paste(collapse = ", ") # Weight (g)
    
    # Reproduction data
    df_R_i <- df_data_EC50_i |> 
      arrange(Time) 
    
    char_tR <- df_R_i |> pull(Time)         |>  paste(collapse = ", ") # Time in days
    char_R  <- df_R_i |> pull(Reproduction) |>  paste(collapse = ", ") # Cumulated reproduction (#)
    
    # Soil OM data
    df_OM_soil_i <- df_data_EC50_i |> 
      arrange(Time) |>  
      dplyr::select(Time, OM_soil, Environment) |> 
      filter(Environment == "Soil_change") 
    
    char_Replace_tOM_soil <- df_OM_soil_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_Replace_OM_soil  <- df_OM_soil_i |> pull(OM_soil) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Replace_soil <- paste(rep("Replace, ", length(df_OM_soil_i$Time)), collapse = "")
    
    # Horse OM data
    df_OM_horse_Replace_i <- df_data_EC50_i |> 
      arrange(Time) |>  
      dplyr::select(Time, OM_horse, Environment) |> 
      filter(Environment == "Soil_change") 
    
    char_Replace_tOM_horse <- df_OM_horse_Replace_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_Replace_OM_horse  <- df_OM_horse_Replace_i |> pull(OM_horse) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Replace_horse <- paste(rep("Replace, ", length(df_OM_horse_Replace_i$Time)), collapse = "")
    
    df_OM_horse_Add_i <- df_data_EC50_i |> 
      arrange(Time) |>  
      dplyr::select(Time, OM_horse, Environment) |> 
      filter(Environment == "Food") 
    
    char_Add_tOM_horse <- df_OM_horse_Add_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_Add_OM_horse  <- df_OM_horse_Add_i |> pull(OM_horse) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Add_horse <- paste(rep("Add, ", length(df_OM_horse_Add_i$Time)), collapse = "")
    
    # Soil quantity in cosm (g)
    df_Soilw_i <- df_data_EC50 |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Soil_w != lag(Soil_w))
    
    char_tsoilw <- df_Soilw_i |> pull(Time)   |>  paste(collapse = ", ") # Time in days
    char_soilw  <- df_Soilw_i |> pull(Soil_w) |>  paste(collapse = ", ") # Nb individual in a cosm
   
    # Density data
    df_density_i <- df_data_EC50_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Density != lag(Density))
    
    char_tdensity <- df_density_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_density  <- df_density_i |> pull(Density) |>  paste(collapse = ", ") # Nb individual in a cosm
    
    # MCSim syntaxe
    
    if (df_data_EC50_i$Experiment_type[1] == "Growth") {
      Init_i <- 1
    } else {
      Init_i <- 2
    }
    
    text_init <- paste("    Init = ", Init_i, ";", sep="") # Init = 1 for juveniles and 2 for adults
    text_Ce0 <- paste("    Ce0=", subset(df_data_EC50_i, Time==0)$Dose, ";", sep="")
    
    text_Winit <- paste("    Winit=", subset(df_data_EC50_i, Time==0)$Weight, ";", sep="")
    
    text_OM_soil_init <- paste("    OM_soil_t0=", subset(df_data_EC50_i, Time==0)$OM_soil, ";", sep="")
    text_OM_horse_init <- paste("    OM_horse_t0=", subset(df_data_EC50_i, Time==0)$OM_horse, ";", sep="")
    
    text_Ww <- paste(
      paste("    Print(Weight,", char_tW, ");" , sep=""),
      paste("    Data(Weight,", char_Ww,  ");" , sep=""), 
      sep = "\n"
    )
    
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
    

    text_OM_final <- paste(
      text_OM_soil_Replace,
      text_OM_horse_Replace,
      text_OM_horse_Add,
      sep="\n"
    )
    
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
    
    text_Texp <- paste("    Texp=",subset(df_data_EC50_i, Time==0)$Texp, ";", sep="")
    
    
    char_i <- paste(
      paste("Experiment { #", ID_i), 
      
      text_init,
      text_Ce0,
      text_Winit,
      text_WeightSoilCosm,
      text_OM_soil_init,
      text_OM_horse_init,
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

f_In_experiments_TKTD_MIX <- function(){
  
  df_data <- f_import_data_DEBTKTD_MIX()
  
  char_final <- ""
  
  for (ID_i in unique(df_data$ID_experiment)) {
    
    
    df_data_i <- df_data |> filter(ID_experiment == ID_i)
    
    # t0 = début manip 
    # W = W0 manip
    # Eh = Estimation par rapport à W0 (faire graph a partir DEB normal)
    # R = 0
    
    ################### Data on juvenile growth 
    
    # Weight data (remove Weigth(t=0))
    df_Ww_i <- df_data_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Weight != lag(Weight))
    
    char_tW <- df_Ww_i |> pull(Time)   |>  paste(collapse = ", ") # Time in days
    char_Ww <- df_Ww_i |> pull(Weight) |>  paste(collapse = ", ") # Weight (g)
    
    # Reproduction data
    df_R_i <- df_data_i |> 
      arrange(Time) 
    
    char_tR <- df_R_i |> pull(Time)         |>  paste(collapse = ", ") # Time in days
    char_R  <- df_R_i |> pull(Reproduction) |>  paste(collapse = ", ") # Cumulated reproduction (#)
    
    # Soil OM data
    df_OM_soil_i <- df_data_i |> 
      arrange(Time) |>  
      dplyr::select(Time, OM_soil, Environment) |> 
      filter(Environment == "Soil_change") 
    
    char_Replace_tOM_soil <- df_OM_soil_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_Replace_OM_soil  <- df_OM_soil_i |> pull(OM_soil) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Replace_soil <- paste(rep("Replace, ", length(df_OM_soil_i$Time)), collapse = "")
    
    # Horse OM data
    df_OM_horse_Replace_i <- df_data_i |> 
      arrange(Time) |>  
      dplyr::select(Time, OM_horse, Environment) |> 
      filter(Environment == "Soil_change") 
    
    char_Replace_tOM_horse <- df_OM_horse_Replace_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_Replace_OM_horse  <- df_OM_horse_Replace_i |> pull(OM_horse) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Replace_horse <- paste(rep("Replace, ", length(df_OM_horse_Replace_i$Time)), collapse = "")
    
    df_OM_horse_Add_i <- df_data_i |> 
      arrange(Time) |>  
      dplyr::select(Time, OM_horse, Environment) |> 
      filter(Environment == "Food") 
    
    char_Add_tOM_horse <- df_OM_horse_Add_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_Add_OM_horse  <- df_OM_horse_Add_i |> pull(OM_horse) |>  paste(collapse = ", ") # OM in horse available for one earthworm (%)
    char_Add_horse <- paste(rep("Add, ", length(df_OM_horse_Add_i$Time)), collapse = "")
    
    # Soil quantity in cosm (g)
    df_Soilw_i <- df_data_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Soil_w != lag(Soil_w))
    
    char_tsoilw <- df_Soilw_i |> pull(Time)   |>  paste(collapse = ", ") # Time in days
    char_soilw  <- df_Soilw_i |> pull(Soil_w) |>  paste(collapse = ", ") # Nb individual in a cosm
    
    # Density data
    df_density_i <- df_data_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Density != lag(Density))
    
    char_tdensity <- df_density_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_density  <- df_density_i |> pull(Density) |>  paste(collapse = ", ") # Nb individual in a cosm
    
    # MCSim syntaxe
    

      Init_i <- 2
    
    text_init <- paste("    Init = ", Init_i, ";", sep="") # Init = 1 for juveniles and 2 for adults
    text_Ce0 <- paste("    Ce0=", subset(df_data_i, Time==0)$Dose, ";", sep="")
    
    text_Winit <- paste("    Winit=", subset(df_data_i, Time==0)$Weight, ";", sep="")
    
    text_OM_soil_init <- paste("    OM_soil_t0=", subset(df_data_i, Time==0)$OM_soil, ";", sep="")
    text_OM_horse_init <- paste("    OM_horse_t0=", subset(df_data_i, Time==0)$OM_horse, ";", sep="")
    
    text_Ww <- paste(
      paste("    Print(Weight,", char_tW, ");" , sep=""),
      paste("    Data(Weight,", char_Ww,  ");" , sep=""), 
      sep = "\n"
    )
    
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
    
    
    text_OM_final <- paste(
      text_OM_soil_Replace,
      text_OM_horse_Replace,
      text_OM_horse_Add,
      sep="\n"
    )
    
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
    
    text_Texp <- paste("    Texp=",subset(df_data_i, Time==0)$Texp, ";", sep="")
    
    
    char_i <- paste(
      paste("Experiment { #", ID_i), 
      
      text_init,
      text_Ce0,
      text_Winit,
      text_WeightSoilCosm,
      text_OM_soil_init,
      text_OM_horse_init,
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

f_In_TKTD_tot <- function(file_path, Molecule, text_priors, text_likelihood, NbIter, seeds, RTOL, ATOL) {
    
    text_Level_global <- 
      'Level{ # Global
    '
    
    text_Level_exp <- '
Level{

  ############## Experiments ###################
  '
    
    text_experiment <- f_In_experiments_TKTD(Molecule)
    
    text_end <- '
} # End
} # End global

End.'
    
    # .in generation and saving
    for (id in names(seeds)) {
      text_start <- f_create_mcmc_block(paste0("DEB_TKTD_", id), seeds[[id]], NbIter, RTOL, ATOL)
      
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
        here::here(file_path, paste0("DEB_TKTD_", id, ".in"))
      )
    }
}

f_In_TKTD_tot <- function(file_path, Molecule, text_priors, text_likelihood, NbIter, seeds, RTOL, ATOL) {
    
    text_Level_global <- 
      'Level{ # Global
    '
    
    text_Level_exp <- '
Level{

  ############## Experiments ###################
  '
    
    text_experiment <- f_In_experiments_TKTD(Molecule)
    
    text_end <- '
} # End
} # End global

End.'
    
    # .in generation and saving
    for (id in names(seeds)) {
      text_start <- f_create_mcmc_block(paste0("DEB_TKTD_", id), seeds[[id]], NbIter, RTOL, ATOL)
      
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
        here::here(file_path, paste0("DEB_TKTD_", id, ".in"))
      )
    }
}


f_In_predobs_TKTD <- function(file_path, Molec, l_param_name){
  
  
  df_data <- f_import_data_DEBTKTD()
  df_data_EC50 <- df_data |> 
    filter(Molecule == Molec) |> 
    mutate(across(where(is.numeric), ~ replace_na(.x, -1)))
  
  compteur_exp <- 0 
  
  for(ID_experiment_i in unique(df_data_EC50$ID_experiment)){
    text_full_i <- ""
    char_i <- ""
    compteur_exp <- compteur_exp + 1
    
    text_start <- paste0('#### DEB TKTD A. caliginosa
#===============================================

SetPoints("DEB_setpoint_predobs_', compteur_exp, '.out", "tab_setpoint.out", 0,', l_param_name|>  paste(collapse = ", "),');

Integrate(Lsodes, 1E-6, 1E-6, 1);

########## Experiments ################################################
')
    
    text_end <- "
  End."
    
    # Data corresponding to this experiment
    df_DEB_i <- subset(df_data_EC50, ID_experiment == ID_experiment_i)
    
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
    
    # Soil quantity in cosm (g)
    df_Soilw_i <- df_DEB_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Soil_w != lag(Soil_w))
    
    char_tsoilw <- df_Soilw_i |> pull(Time)   |>  paste(collapse = ", ") # Time in days
    char_soilw  <- df_Soilw_i |> pull(Soil_w) |>  paste(collapse = ", ") # Nb individual in a cosm
    
    # Density data
    df_density_i <- df_DEB_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Density != lag(Density))
    
    char_tdensity <- df_density_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_density  <- df_density_i |> pull(Density) |>  paste(collapse = ", ") # Nb individual in a cosm
    
    # MCSim syntaxe
    
    if (df_DEB_i$Experiment_type[1] == "Growth") {
      Init_i <- 1
    } else {
      Init_i <- 2
    }
    
    text_init <- paste("    Init = ", Init_i, ";", sep="") # Init = 1 for juveniles and 2 for adults
    text_Ce0 <- paste("    Ce0=", subset(df_DEB_i, Time==0)$Dose, ";", sep="")
    
    text_Winit <- paste("    Winit=", subset(df_DEB_i, Time==0)$Weight, ";", sep="")
    
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
    
    
    

      text_OM_final <- paste(
        text_OM_soil_Replace,
        text_OM_horse_Replace,
        text_OM_horse_Add,
        sep="\n"
      )

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
    
    
    # Weight data times
    df_Ww_i <- df_DEB_i |> 
      arrange(Time)             
    char_tW <- df_Ww_i |> pull(Time)   |>  paste(collapse = ", ") # Time in days
    # Reproduction data times
    df_R_i <- df_DEB_i |> 
      arrange(Time)
    char_tR <- df_R_i |> pull(Time)         |>  paste(collapse = ", ") # Time in days
    
    text_print <- paste0(
      "    Print(Weight,", char_tW, ");
    Print(Reproduction,", char_tR,");"
    )
    
    
    char_i <- paste(
      paste("Simulation { #", ID_experiment_i), 
      
      text_init,
      text_Ce0,
      text_Winit,
      text_WeightSoilCosm,
      text_OM_soil_init,
      text_OM_horse_init,
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

f_In_predobs_TKTD_MIX <- function(file_path, l_param_name){
  
  
  df_data <- f_import_data_DEBTKTD_MIX() |> 
    mutate(across(where(is.numeric), ~ replace_na(.x, -1)))
  
  compteur_exp <- 0 
  
  for(ID_experiment_i in unique(df_data$ID_experiment)){
    text_full_i <- ""
    char_i <- ""
    compteur_exp <- compteur_exp + 1
    
    text_start <- paste0('#### DEB TKTD A. caliginosa
#===============================================

SetPoints("DEB_setpoint_predobs_', compteur_exp, '.out", "tab_setpoint.out", 0,', l_param_name|>  paste(collapse = ", "),');

Integrate(Lsodes, 1E-6, 1E-6, 1);

########## Experiments ################################################
')
    
    text_end <- "
  End."
    
    # Data corresponding to this experiment
    df_DEB_i <- subset(df_data, ID_experiment == ID_experiment_i)
    
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
    
    # Soil quantity in cosm (g)
    df_Soilw_i <- df_DEB_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Soil_w != lag(Soil_w))
    
    char_tsoilw <- df_Soilw_i |> pull(Time)   |>  paste(collapse = ", ") # Time in days
    char_soilw  <- df_Soilw_i |> pull(Soil_w) |>  paste(collapse = ", ") # Nb individual in a cosm
    
    # Density data
    df_density_i <- df_DEB_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Density != lag(Density))
    
    char_tdensity <- df_density_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_density  <- df_density_i |> pull(Density) |>  paste(collapse = ", ") # Nb individual in a cosm
    
    # MCSim syntaxe
    

      Init_i <- 2
    
    text_init <- paste("    Init = ", Init_i, ";", sep="") # Init = 1 for juveniles and 2 for adults
    text_Ce0_IMD <- paste("    Ce0_IMD=", subset(df_DEB_i, Time==0)$Dose_IMD, ";", sep="")
    text_Ce0_EPX <- paste("    Ce0_EPX=", subset(df_DEB_i, Time==0)$Dose_EPX, ";", sep="")
    
    text_Winit <- paste("    Winit=", subset(df_DEB_i, Time==0)$Weight, ";", sep="")
    
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
    
    
    
    
    text_OM_final <- paste(
      text_OM_soil_Replace,
      text_OM_horse_Replace,
      text_OM_horse_Add,
      sep="\n"
    )
    
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
    
    
    # Weight data times
    df_Ww_i <- df_DEB_i |> 
      arrange(Time)             
    char_tW <- df_Ww_i |> pull(Time)   |>  paste(collapse = ", ") # Time in days
    # Reproduction data times
    df_R_i <- df_DEB_i |> 
      arrange(Time)
    char_tR <- df_R_i |> pull(Time)         |>  paste(collapse = ", ") # Time in days
    
    text_print <- paste0(
      "    Print(Weight,", char_tW, ");
    Print(Reproduction,", char_tR,");"
    )
    
    
    char_i <- paste(
      paste("Simulation { #", ID_experiment_i), 
      
      text_init,
      text_Ce0_IMD,
      text_Ce0_EPX,
      text_Winit,
      text_WeightSoilCosm,
      text_OM_soil_init,
      text_OM_horse_init,
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


f_In_Setpoint_full_TKTD <- function(file_path, Molec, l_param_name, step_sim, Endpoints_print){
  
  
  df_data <- f_import_data_DEBTKTD()
  df_data_EC50 <- df_data |> 
    filter(Molecule == Molec) |> 
    mutate(across(where(is.numeric), ~ replace_na(.x, -1)))
  
  compteur_exp <- 0 
  
  for(ID_experiment_i in unique(df_data_EC50$ID_experiment)){
    text_full_i <- ""
    char_i <- ""
    compteur_exp <- compteur_exp + 1
    
    text_start <- paste0('#### DEB TKTD A. caliginosa
#===============================================

SetPoints("DEB_setpoint_full_', compteur_exp, '.out", "tab_setpoint.out", 0,', l_param_name|>  paste(collapse = ", "),');

Integrate(Lsodes, 1E-9, 1E-10, 1);

########## Experiments ################################################
')
    
    text_end <- "
  End."
    
    # Data corresponding to this experiment
    df_DEB_i <- subset(df_data_EC50, ID_experiment == ID_experiment_i)
    
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
    
    # Soil quantity in cosm (g)
    df_Soilw_i <- df_DEB_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Soil_w != lag(Soil_w))
    
    char_tsoilw <- df_Soilw_i |> pull(Time)   |>  paste(collapse = ", ") # Time in days
    char_soilw  <- df_Soilw_i |> pull(Soil_w) |>  paste(collapse = ", ") # Nb individual in a cosm
    
    # Density data
    df_density_i <- df_DEB_i |> 
      arrange(Time) |>                    
      filter(row_number() == 1 | Density != lag(Density))
    
    char_tdensity <- df_density_i |> pull(Time)    |>  paste(collapse = ", ") # Time in days
    char_density  <- df_density_i |> pull(Density) |>  paste(collapse = ", ") # Nb individual in a cosm

    # MCSim syntaxe
    
    if (df_DEB_i$Experiment_type[1] == "Growth") {
      Init_i <- 1
    } else {
      Init_i <- 2
    }
    
    text_init <- paste("    Init = ", Init_i, ";", sep="") # Init = 1 for juveniles and 2 for adults
    text_Ce0 <- paste("    Ce0=", subset(df_DEB_i, Time==0)$Dose, ";", sep="")
    
    text_Winit <- paste("    Winit=", subset(df_DEB_i, Time==0)$Weight, ";", sep="")
    
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

      text_OM_final <- paste(
        text_OM_soil_Replace,
        text_OM_horse_Replace,
        text_OM_horse_Add,
        sep="\n"
      )
    
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
      text_Ce0,
      text_Winit,
      text_WeightSoilCosm,
      text_OM_soil_init,
      text_OM_horse_init,
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

f_read_Setpoint_full_TKTD <- function(file_path, Molec, select_times, step_sim, bol_Eh = FALSE){
  
  df_pred_full <- NULL
  
  df_data <- f_import_data_DEBTKTD()
  df_data_obs <- df_data |> 
    filter(Molecule == Molec) 
  
  i <- 0
  
  for (ID_i in unique(df_data_obs$ID_experiment)){
    
    i <- i + 1
    
    Sim.Res.Exp_i <- readr::read_tsv(
      file.path(file_path, paste0("DEB_setpoint_full_", i, ".out"))
    ) |>
      dplyr::mutate(No_sim = i) |> 
      dplyr::select(where(~ !all(is.na(.x))))
    
    if (bol_Eh) {
      Nom_Endpoints = c("Weight", "Maturity", "Energy")
      Nsortie = 3
    } else {
      Nom_Endpoints = c("Weight", "Reproduction", "Ci")
      Nsortie = 3
    }
    
    df_data_obs_i <- subset(df_data_obs, ID_experiment == ID_i)
    
    Time_pred <- seq(0, max(df_data_obs_i$Time, na.rm=TRUE), step_sim) # !!!
    
    
    MPV = IC_min = IC_max = Endpoint = index = NULL
    
    for (j in 1:length(Nom_Endpoints)){
      
      index = grep(paste0("^",Nom_Endpoints[j],"_") , colnames(Sim.Res.Exp_i), fixed=FALSE)
      
      
      Data_sim = as.matrix( Sim.Res.Exp_i[ , index] )
      
      MPV    = c(   MPV,  as.numeric( Data_sim[nrow(Data_sim), ] ) )
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

f_get_temperatures_Closeaux <- function(Year, Start_day, Sim_duration){
  
  token_NOAA <- "ZFyQJPeWOXcEbvwIDWBvVJmfyoVPoTll"
  library(worldmet)
  library(lubridate)
  
  # Meadow location
  lat <- 48.8
  lon <- 2.083333
  
  # 1) Find closest station
  stations <- import_isd_stations(
    lat = lat,
    lon = lon,
    n_max = 10
  )
  
  # 2) Take closest station
  station_code <- stations$code[1]
  station_name <- stations$station[1]
  
  print(paste0("Station : ", station_name, " (Code : ", station_code, ")"))
  
  # 3) Get hourly temperature measurement for (a) given year(s)
  meteo <- import_isd_hourly(
    code = station_code,
    year = Year # Year can be a list of year
  )
  
  # 4) Keep only temperature
  temp <- meteo %>%
    dplyr::select(date, station = code, site = station, latitude, longitude, air_temp) %>%
    filter(!is.na(air_temp)) |> 
    mutate(month_day = format(date, "%y-%m-%d"))
  
  df_temp_sum <- temp |> 
    group_by(month_day, .drop = TRUE) |> 
    mutate(
      air_temp_mean_day = round(mean(air_temp, na.rm=TRUE),1)
    ) |> 
    ungroup() |> 
    dplyr::select(c(month_day, air_temp_mean_day)) |> 
    unique() |> 
    mutate(
      Time = round(as.numeric(difftime(month_day, Start_day, units = "days")))
    ) |> 
    filter((Time >= 0) & (Time <= Sim_duration))
  
  return(df_temp_sum)
}

f_In_sim_predictions <- function(file_path, Year, Start_day, Sim_duration, l_concentrations, bol_degradation, step_sim, No_scenario){
  
  # bol_degradation : 0 if no degradtion ; 1 if degradation
  
  WeightSoilCosm <- 10000
  Per_OM_soil <- 0.0326
  Per_OM_horse <- 1/400 # As if there is 1g of horse dung in a 400 g cosm
  
  df_temp_sim <- f_get_temperatures_Closeaux(Year, Start_day, Sim_duration)
  
  text_end <- '
End.'
  
  text_init <- "    Init = 1;"
  text_Winit <- "    Winit=0.0124;" # Initial weight of the controls in chap 5 (12.4 mg)
  text_dens <- "    Dens = 3;"
  text_Ci0 <- "     Ci0 = 0;"
  text_WeightSoilCosm <- paste0("   WeightSoilCosm =", WeightSoilCosm, ";")
  OMsoil0 <- WeightSoilCosm * Per_OM_soil
  OMhorse0 <- WeightSoilCosm * Per_OM_horse
  text_OMsoil_init<- paste0("   OM_soil_t0=", OMsoil0, ";")
  text_OMhorse_init <- paste0("   OM_horse_t0=", OMhorse0, ";")
  text_bol_degradation <- paste0("    bol_degradation = ", bol_degradation, ";") 
  
  char_Temp_t <- df_temp_sim |> pull(Time) |>  paste(collapse = ", ")
  char_Temp <- df_temp_sim |> pull(air_temp_mean_day) |>  paste(collapse = ", ")
  
  text_OM_temp <- paste("    Texp=NDoses(", 
                        length(df_temp_sim$Time), ",\n                          ", 
                        char_Temp,",\n                          ", 
                        char_Temp_t, ");", sep="")
  text_basis <- paste(
    text_init,
    text_Winit,
    text_OMhorse_init,
    text_OMsoil_init,
    text_WeightSoilCosm,
    text_dens,
    text_Ci0,
    text_bol_degradation,
    text_OM_temp,
    sep="\n"
  )
  
  text_sim <- paste0("    PrintStep(Weight, Reproduction, Organic_matter, Stress, Ce,0,", Sim_duration, ", ", step_sim, ");")
  
  for (i in 1:length(l_concentrations)) {
    
    text_start <- paste0('#### DEB TKTD A. caliginosa
#===============================================

SetPoints("DEB_IMD_prediction_', No_scenario, "_", i, '.out", "tab_setpoint.out", 0, nec, b, Sigma_W);

Integrate(Lsodes, 1E-10, 1E-8, 1);

########## Experiments ################################################')
    
    text_Ce0 <- paste0("    Ce0 = ", l_concentrations[i], ";")
    
    char_i <- paste(
      paste("Simulation { #", i), 
      
      text_basis,
      text_Ce0,
      text_sim,
      
      paste("}"),
      sep="\n"
      
    )
    
    Text_final_i <- paste(text_start, char_i, text_end, sep="\n")
    
    writeLines(
      Text_final_i,
      here::here(file_path, paste0("DEB_IMD_prediction_", No_scenario, "_", i, ".in"))
    )
  }
  
}


f_Read_fichier <- function(fichier){
  contenu_raw <- readBin(
    fichier,
    what = "raw",
    n = file.info(fichier)$size
  )
  
  contenu <- rawToChar(contenu_raw)
  contenu <- gsub("\r\n|\r", "\n", contenu)
  
  fichier_temp <- tempfile(fileext = ".out")
  on.exit(unlink(fichier_temp))
  
  writeBin(
    charToRaw(contenu),
    fichier_temp
  )
  
  res <- readr::read_tsv(
    fichier_temp,
    show_col_types = FALSE,
    progress = FALSE
  )
  return(res)
}

f_Read_fichier2 <- function(fichier) {
  contenu_raw <- readBin(
    fichier,
    what = "raw",
    n = file.info(fichier)$size
  )
  
  contenu <- rawToChar(contenu_raw)
  contenu <- gsub("\r\n|\r", "\n", contenu)
  
  # Ajout du nom de la première colonne
  contenu <- sub(
    "^",
    "Iter ",
    contenu
  )
  
  fichier_temp <- tempfile(fileext = ".out")
  on.exit(unlink(fichier_temp), add = TRUE)
  
  writeBin(
    charToRaw(contenu),
    fichier_temp
  )
  
  res <- readr::read_table(
    fichier_temp,
    show_col_types = FALSE,
    progress = FALSE
  )
  return(res)
}

f_Read_predictions <- function(file_path, No_scenario, l_concentrations, Sim_duration, step_sim){
  
  df_pred_full <- NULL
  
  for (i in 1:length(l_concentrations)){
    
    fichier <- file.path(file_path, paste0("DEB_IMD_prediction_", No_scenario,"_", i, ".out"))
  
    Sim.Res.Exp_i <- f_Read_fichier(fichier)
    
      Nom_Endpoints = c("Weight", "Reproduction", "Organic_matter", "Stress", "Ce")
      Nsortie = 5
    
    Time_pred <- seq(0, Sim_duration, step_sim) # !!!
    
    MPV = IC_min = IC_max = Endpoint = index = NULL
    
    for (j in 1:length(Nom_Endpoints)){
      
      index = grep(paste0("^",Nom_Endpoints[j],"_") , colnames(Sim.Res.Exp_i), fixed=FALSE)
      
      Data_sim = as.matrix( Sim.Res.Exp_i[ , index] )
      
      MPV    = c(   MPV,  as.numeric( Data_sim[nrow(Data_sim), ] ) )
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
      mutate(
        predict.endpoint = case_when(
          ( predict.endpoint <= 5e-6 & Endpt == "Reproduction") ~ 0,
          .default = predict.endpoint
        ),
        Ce = l_concentrations[i]
      )
    
    df_pred_full <- rbind(df_pred_full, df_pred_i)
  }
  
  saveRDS(df_pred_full, file.path(file_path, paste0("df_sim_", No_scenario, ".rds")))
  
  return(df_pred_full)
  
}


f_Read_setpoint_DEB_refit <- function(file_path, step_sim, max_time, DEB_type) {

  
  df_pred_full <- NULL
  
  for (i in c(1,2)){
    
    Sim.Res.Exp_i <- readr::read_tsv(
      file.path(file_path, paste0("DEB_setpoint_full_", DEB_type, "_", i, ".out"))
    ) |>
      dplyr::mutate(No_sim = i) |> 
      dplyr::select(where(~ !all(is.na(.x))))
    
    Nom_Endpoints = c("Weight", "Reproduction")
    Nsortie = 2
    
    
    Time_pred <- seq(0, max_time, step_sim) # !!!
    
    
    MPV = IC_min = IC_max = Endpoint = index = NULL
    
    for (j in 1:length(Nom_Endpoints)){
      
      index = grep(paste0("^",Nom_Endpoints[j],"_") , colnames(Sim.Res.Exp_i), fixed=FALSE)
      
      
      Data_sim = as.matrix( Sim.Res.Exp_i[ , index] )
      
      MPV    = c(   MPV,  as.numeric( Data_sim[nrow(Data_sim), ] ) )
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
      mutate(
        predict.endpoint = case_when(
          ( predict.endpoint <= 5e-6 & Endpt == "Reproduction") ~ 0,
          .default = predict.endpoint
        ),
        No_sim = i
      )
    
    df_pred_full <- rbind(df_pred_full, df_pred_i)
  }

saveRDS(df_pred_full, file.path(file_path, paste0("Sim.Res.Exp.full_", DEB_type, ".rds")))
        
return(df_pred_full)
}


# Plots ----

f_Plot_MCMC_Chains <- function(
  df, 
  bol_log = TRUE
) {
  
  pal_chains <- c(Nord_aurora, Nord_frost, Nord_polar)
  
  p_chains <-  
    ggplot(
      data = df$Chains |> 
        filter(Iteration !=1), 
      aes(
        x = Iteration, 
        y = value, 
        color = as.factor(Chain), 
        group = as.factor(Chain)
      )
    )+
    geom_line(alpha=0.7)+
    facet_wrap(~Parameter, scales = "free", ncol = 2)+
    
    scale_color_manual(name = "Chains", values = pal_chains)+
    
    theme_bw()+
    theme(
      legend.position = "right", 
      title=element_text(size=12, face="plain"), 
      axis.title.x = element_text(face="plain"),
      strip.background = element_rect(fill="white")
    )
  
  if(bol_log == TRUE){
    p_chains <- p_chains + scale_y_log10()
  }
  
  return(p_chains)
  
}

f_Plots_MCMC_Likelihood <- function(df) {
  
  pal_chains <- c(Nord_aurora, Nord_frost, Nord_polar)
  lab_chains <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
  
  p1 <- ggplot(
    df$df_LnPost_start |> 
      filter(iteration !=1), 
    aes(
      x = iteration, 
      y = LnPosterior, 
      color = Chain
    )
  ) +
    geom_line() +
    scale_color_manual(name = "Chains", label = lab_chains, values = pal_chains) +
    labs(
      x = "First half iterations", 
      y = "LnPosterior"
    ) +
    theme_minimal()

  p2 <- ggplot(
    df$df_LnPost |> 
      filter(iteration !=1), 
    aes(
      x = iteration, 
      y = LnPosterior, 
      color = Chain
    )
  ) +
    geom_line() +
    scale_color_manual(name = "Chains", label = lab_chains, values = pal_chains) +
    labs(
      x = paste(df$Nb_iter_kept, "last iterations"), 
      y = "LnPosterior"
    ) +
    theme_minimal()
  

  p3 <- ggplot(
    df$df_Devc |> 
      filter(iteration !=1), 
    aes(
      x = iteration, 
      y = RatioDeviance, 
      color = Chain
    )
  ) +
    geom_line() +
    scale_color_manual(name = "Chains", label = lab_chains, values = pal_chains) +
    labs(
      x = paste(df$Nb_iter_kept, "last iterations"), 
      y = "Ratio Deviance"
    ) +
    theme_minimal()
  
  p_LL <- p1 + p2 + p3 +
    plot_layout(ncol=3, guides = "collect")
  
  return(p_LL)
  
}


f_Plot_Convergence <- function(df){
  # Large format per chain
  mcmc_list_df <- lapply(
    split(
      df$Chains, 
      df$Chains$Chain
    ),
    function(d) {
      wide <- d %>%
        dplyr::select(Iteration, Parameter, value) %>%
        pivot_wider(names_from = Parameter, values_from = value) %>%
        arrange(Iteration)
      mcmc(as.matrix(wide[,-1]))
    })
  
  mcmc_list_df <- as.mcmc.list(mcmc_list_df)
  
  par(mfrow=c(1,3))
  
  gelman.plot(mcmc_list_df,
              bin.width = 10,      # Number of observations per segment,
              # first segment --> at least 50 iterations.
              max.bins = 50,       # Maximum number of bins, excluding the last one
              confidence = 0.95,   # Coverage probability of confidence interval.
              transform = FALSE,   # improve the normality of the distribution
              autoburnin=TRUE)     # Remove first half of sequence
  
  gelman <- gelman.diag(
    mcmc_list_df, 
    confidence = 0.95, 
    transform=FALSE, 
    autoburnin=TRUE, 
    multivariate=TRUE
  )
  
  return(gelman)
}

f_mcmc_list_corr <- function(df) {
  
  mcmc_list_df <- lapply(
    split(
      df$Chains, 
      df$Chains$Chain
    ),
    function(d) {
      wide <- d %>%
        dplyr::select(Iteration, Parameter, value) %>%
        pivot_wider(names_from = Parameter, values_from = value) %>%
        arrange(Iteration)
      mcmc(as.matrix(wide[,-1]))
    })
  
  mcmc_list_df <- as.mcmc.list(mcmc_list_df)
  
  mcmc_list_corr_df <- mcmc.list(
    lapply(mcmc_list_df, function(x) {
      mcmc(as.matrix(x)[-1, , drop = FALSE])  # <- reconvertir en mcmc
    })
  )
  return(mcmc_list_corr_df)
}

f_Plot_ipairs <- function(df, file_rds, size_label = 0.5) {
  
  mcmc_list_corr_df <- f_mcmc_list_corr(df)
  
  colramp_aurora <- colorRampPalette(rev(Nord_aurora[1:4]))
  
  if (!file.exists(file_rds)) {
    
    png(file_rds, width = 3000, height = 3000, res = 300)
    
    ipairs(
      as.matrix(mcmc_list_corr_df), 
      colramp = colramp_aurora, 
      pixs=0.25, 
      cex.diag = size_label
    )
    dev.off()
  }
  
  ipairs(
    as.matrix(mcmc_list_corr_df), 
    colramp = colramp_aurora, 
    pixs=0.25, 
    cex.diag = 0.5
  )
}

f_Get_Prior_distrib <- function(df){
  
  n_samples <- 1e5 # Simulation of the distribution
  
  df_priors_dist <- df$Priors |> 
    dplyr::select(Nom, Distribution, P1, P2, P3, P4) |> 
    rowwise() |> 
    mutate(samples = case_when(
      Distribution == "Uniform" ~ list(runif(n_samples, min = P1, max = P2)),
      Distribution == "Normal" ~ list(rnorm(n_samples, mean = P1, sd = P2)),
      Distribution == "LogUniform" ~ list(exp(runif(n_samples, min = log(P1), max = log(P2)))),
      Distribution == "TruncNormal" ~ list(rtruncnorm(n_samples, a = P3, b = P4, mean = P1, sd = P2)),     Distribution == "TruncNormal_cv" ~ list(rtruncnorm(n_samples, a = P3, b = P4, mean = P1, sd = P2*P1))
    )) |> 
    unnest(samples)
  
  return(df_priors_dist)
}

f_Plot_Posteriors <- function(df, bol_log = TRUE){
  
  L_parameters_df <- unique(df$Priors$Nom)
  
  df_priors_dist <- f_Get_Prior_distrib(df)
  
  alpha_prior <- 0.5
  alpha_post <- 0.5
  col_prior <- Nord_snow[1]
  col_post <- col_IMD
  
  df$Chains_plot <- df$Chains |> 
    mutate(Nom = str_remove(Parameter, "\\.\\d+\\.$")) 
  
  p_post <- ggplot() +
    
    # Priors
    stat_halfeye(
      data = df_priors_dist,
      aes(
        x = samples,
      ),
      alpha = alpha_prior,
      fill = col_prior,
      color = col_prior,
      normalize = "xy"
    ) +
    
    # Posterior
    stat_halfeye(
      data = subset(df$Chains_plot, Nom %in% L_parameters_df),
      mapping = aes(
        x=value
      ),
      fill = col_post,
      color = col_post,
      alpha = alpha_post,
      normalize = "xy"
    ) +
    
    facet_wrap(~Nom, scales = "free")+
    labs(
      title = "Posterior and prior distributions", 
      x = "Parameter value", 
      y = "Density"
    ) +
    theme_minimal()+
    theme(
      title=element_text(face="bold")
    )
  
  if (bol_log ==TRUE) {
    p_post <- p_post + scale_x_log10()
  }
  
  return(p_post)
  
}