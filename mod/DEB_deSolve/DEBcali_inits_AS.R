# 0. sub functions ----

f_events <- function(l_Experiments, No_experiment){
  
  df_data <- f_import_data_DEB(l_Experiments, Adults_alone = TRUE) %>% 
    filter(No_sim == No_experiment)
  
  # Event OM_soil
  
  df_data_OM_soil <- df_data |> 
    arrange(Time) |>  
    dplyr::select(Time, OM_soil, Environment) |> 
    filter(Environment == "Soil_change")
  
  Nb_times_OM_soil <- length(df_data_OM_soil$Time)
  
  df_events_OM_soil <- df_data_OM_soil %>% 
    mutate(
      var = rep("OM_soil", Nb_times_OM_soil),
      time = Time,
      value = OM_soil,
      method = rep("replace", Nb_times_OM_soil)
    ) %>% 
    dplyr::select(-c(Time, OM_soil, Environment))
  
  # Event OM_horse
  
  df_data_OM_horse <- df_data |> 
    arrange(Time) |>  
    dplyr::select(Time, OM_horse, Environment) %>%
    filter(Environment %in% c("Soil_change", "Food", "End"))
  
  Nb_times_OM_horse <- length(df_data_OM_horse$Time)
  
  df_events_OM_horse <- df_data_OM_horse |> 
    mutate(
      var = rep("OM_horse", Nb_times_OM_horse),
      time = Time,
      value = OM_horse,
      method =  case_when(
        Environment == "Soil_change" ~ "replace",
        Environment == "Food" ~ "add",
        Environment == "End" ~ "replace"
      )
    ) %>% 
    dplyr::select(-c(Time, OM_horse, Environment))
  
  res <- rbind(df_events_OM_soil,df_events_OM_horse) |> 
    arrange(time)
  return(res)
}

f_simulation_times <- function(l_Experiments, No_experiment, step){
  
  df_data <- f_import_data_DEB(l_Experiments, Adults_alone = TRUE) %>% 
    filter(No_sim == No_experiment)
  
  Max_times <- max(df_data$Time)
  
  res <- seq(0, Max_times, step)
  
  return(res)
  
}

f_forcing_inputs <- function(l_Experiments, No_experiment){
  
  df_data_Dens <- f_import_data_DEB(l_Experiments, Adults_alone = TRUE) %>% 
    filter(No_sim == No_experiment) %>% 
    dplyr::select(c(Time, Density)) %>% 
    arrange(Time) |>                    
    filter(
      Density != lag(Density, default = first(Density)) |
        Density != lead(Density, default = last(Density)) |
        row_number() == 1 |
        row_number() == n()
    )
  
  res <- list(df_data_Dens)
  
}

f_sobol_l_parms <- function(Sobol_type, seed, Nb_sim, Processors, Nboot, Nout, l_params_median, Temperature_bol){
  
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
  
  l_Lower <- l_params_median - l_params_median*0.1
  l_Upper <- l_params_median + l_params_median*0.1
  l_Lower <- pmax(l_Lower, 0)
  
  # Special cases
  l_Lower["kap"] <- 0.1
  l_Upper["kap"] <- 0.9
  l_Lower["rOM_ClxHorse"] <- 0
  l_Upper["rOM_ClxHorse"] <- 1
  
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

# 1. Initialisation ----

f_initParms <- function(l_parms, No_experiment, newParms = NULL){
  
  df_data_t0       <- subset(f_import_data_DEB(l_Experiments, Adults_alone = TRUE), (No_sim == No_experiment & Time == 0))
  Winit_i          <- df_data_t0$Weight[1]
  Texp_i           <- df_data_t0$Texp[1]
  WeightSoilCosm_i <- df_data_t0$Soil_w[1]
  
  parms <- c(
    # Paramètres 
    l_parms,
    # Inputs MCSim
    Winit = Winit_i,                   # Masse à t0 (g)
    Texp = Texp_i,                     # Température de l'expérience (°C)
    WeightSoilCosm = WeightSoilCosm_i  # Masse de sol dans le cosme (g)
  )
  
  parms_order <- c(
    "muOM", "rOM_ClxHorse", "Fm", "kapX", "pAm", "v", "kap",
    "pM", "pT", "Eg", "Shape", "w", "kJ", "Ehb", "Ehp",
    "L_coc", "E_coc", "kapR", "TAH", "TH", "TA", "Tref",
    "Winit", "Texp", "WeightSoilCosm"
  )
  
  parms <- as.numeric(parms[parms_order])
  names(parms) <- parms_order
  
  if (!is.null(newParms)) {
    if (!all(names(newParms) %in% c(names(parms)))) {
      stop("illegal parameter name")
    }
    parms[names(newParms)] <- newParms
  }
  
  out <- .C("getParms",  as.double(parms),
            out=double(length(parms)),
            as.integer(length(parms)))$out
  names(out) <- names(parms)
  return(out)
}

f_initStates <- function(parms, newStates = NULL){
  
  Y <- c(
    E = as.numeric(parms["pAm"]) / as.numeric(parms["v"]),
    L = (as.numeric(parms["Winit"])/(1+1*as.numeric(parms["w"])))^(1/3),   # pow((Winit/(1+1*w)),(1.0/3.0))
    Eh = as.numeric(parms["Ehb"]),
    R = 1e-6,
    OM_soil = 0.0,
    OM_horse = 0.0,
    W = as.numeric(parms["Winit"]),
    OM = 0
  )
  
  if (!is.null(newStates)) {
    if (!all(names(newStates) %in% c(names(Y)))) {
      stop("illegal state variable name in newStates")
    }
    Y[names(newStates)] <- newStates
  }
  
  .C("initState", as.double(Y));
  
  return(Y)
}

# 2. AS ----


f_AS_DEB <- function(Sobol_type, l_Experiments, No_experience, Nb_sim, Processors, Nboot, Nout, l_params_median, Temperature_bol){
  
  set.seed(1212)
  seed <- 1212
  Step_sim     <- 0.1
  
  # Model outputs
  Outputs <- c(
    "Energy",
    "Maturity",
    "Reproduction",
    "Weight",
    "Organic_matter"
  )
  
  Sobol_tot <- f_sobol_l_parms(Sobol_type, seed, Nb_sim, Processors, Nboot, Nout, l_params_median, Temperature_bol)
  l_parms <- Sobol_tot$Sobol_param
  Sobol_info <- Sobol_tot$Sobol
  
  NbIter <- length(l_parms)
  print(paste("Expérience :", l_Experiments[No_experience]))
  
  times_sim <- f_simulation_times(l_Experiments, No_experience, Step_sim)
  df_events <- f_events(l_Experiments, No_experience)
  forced_inputs <- f_forcing_inputs(l_Experiments, No_experience)
  
  df_AS_sim <- data.frame()
  
  for (i in seq(1,NbIter, 1)){
    
    l_parms_i <- l_parms[[i]]
    if (!Temperature_bol){
      l_parms_i <- c(l_parms_i, l_params_median[(length(l_params_median)-4):length(l_params_median)]) # Add TAH TH TA Tref back
    } 
    parms_deb_i <- f_initParms(l_parms_i, No_experience)
    Y_deb_i <- f_initStates(parms_deb_i)
    
    res <- tryCatch({
      
      # Résolution de l'ODE
      result_deb_i <- deSolve::ode(
        y = Y_deb_i,
        times_sim,
        func = "derivs",
        parms = parms_deb_i,
        rtol = 1e-5,
        atol = 1e-5,
        maxsteps = 10000,
        dllname = Model.name,
        initfunc = "initmod",
        initforc = "initforc",
        forcings = forced_inputs,
        fcontrol = list(method = "constant"),
        events = list(data = df_events),
        nout = length(Outputs),
        outnames = Outputs
      )
      
      # Conversion en dataframe et mise en forme
      df_result_deb_i <- as.data.frame(result_deb_i) |>
        pivot_longer(
          cols = -time,
          names_to = "Variable",
          values_to = "Value"
        ) |> 
        mutate(
          No_experiment = No_experience,
          Iter = i
        )
      
      get_val <- function(df, t, var){
        x <- subset(df, time == t & Variable == var)$Value
        if(length(x) == 0) return(NA)
        unique(x)[1]
      }
      
      df_result_deb_i_Endpoints <- data.frame(
        Weight_Juv     = get_val(df_result_deb_i, 50,  "Weight"),
        Weight_Adult   = get_val(df_result_deb_i, 130, "Weight"),
        Energy_Juv     = get_val(df_result_deb_i, 50,  "Energy"),
        Energy_Adult   = get_val(df_result_deb_i, 130, "Energy"),
        Maturity_Juv   = get_val(df_result_deb_i, 50,  "Maturity"),
        Maturity_Adult = get_val(df_result_deb_i, 130, "Maturity"),
        Reproduction   = get_val(df_result_deb_i, 130, "Reproduction"),
        Organic_matter = get_val(df_result_deb_i, 50,  "Organic_matter")
      )
      
      df_AS_sim <- rbind(df_AS_sim, df_result_deb_i_Endpoints)
      df_AS_sim  # retour normal
      
    }, error = function(e) {
      message(sprintf("Erreur à l'expérience %d, itération %d : %s", No_experience, i, e$message))
      
      df_result_deb_i_Endpoints <- data.frame(
        Weight_Juv     = get_val(df_result_deb_i, 50,  "Weight"),
        Weight_Adult   = get_val(df_result_deb_i, 130, "Weight"),
        Energy_Juv     = get_val(df_result_deb_i, 50,  "Energy"),
        Energy_Adult   = get_val(df_result_deb_i, 130, "Energy"),
        Maturity_Juv   = get_val(df_result_deb_i, 50,  "Maturity"),
        Maturity_Adult = get_val(df_result_deb_i, 130, "Maturity"),
        Reproduction   = get_val(df_result_deb_i, 130, "Reproduction"),
        Organic_matter = get_val(df_result_deb_i, 50,  "Organic_matter")
      )
      df_AS_sim <- rbind(df_AS_sim, df_result_deb_i_Endpoints)
      
      return(df_AS_sim)
    })
    
    print(i)
    
  } # End itération 
  
  #save(df_AS_sim, file = here::here("mod/DEB_deSolve/Resultat_AS_Sobolo_df_AS_sim.RData"))
  Nb_fail <- sum(rowSums(is.na(df_AS_sim)) > 0)
  print(paste0("Number of failed integration : ", Nb_fail))
  # Resolution des NA -> Moyenne du reste
  
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
    
    Y <- df_AS_sim[,i]
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
  
  res <- list(
    df_res = df_res,
    Nb_fail = Nb_fail
  )
  
  return(res)
}
