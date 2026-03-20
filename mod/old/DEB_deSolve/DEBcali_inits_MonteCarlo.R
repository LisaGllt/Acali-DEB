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

f_MonteCarlo_l_parms <- function(seed, Nb_iter, l_params_median, CV_parms, Parms_fixed, Parms_unif, Limits_parms_unif){
  
  # Nb_sim : Number of simulations
  
  library(sensitivity)
  set.seed(seed)
  
  l_param_names <- as.character(names(l_params_median))
  Nb_param_AS <- length(l_param_names)
  
  # Design parameters
  
  df_MonteCarlo_draws <- as.data.frame(lapply(names(l_params_median), function(p) {
    if (p %in% Parms_fixed) {
      rep(l_params_median[p], Nb_iter)
    } else if (p %in% Parms_unif) {
      runif(Nb_iter, Limits_parms_unif[[p]][1], Limits_parms_unif[[p]][2])
    } else {
      rnorm(Nb_iter, l_params_median[p], CV_parms * l_params_median[p])
    }
  }))
  names(df_MonteCarlo_draws) <- names(l_params_median)

  
  res <- apply(df_MonteCarlo_draws, 1, function(x) {
    as.list(x)
  })
  
  return(res)
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
    E = as.numeric(parms["pAm"]) /as.numeric( parms["v"][1]),
    L = (as.numeric(parms["Winit"])/(1+1*as.numeric(parms["w"])))^(1/3),   # pow((Winit/(1+1*w)),(1.0/3.0))
    Eh = as.numeric(parms["Ehb"]),
    R = 1e-6,
    OM_soil = 1e-6,
    OM_horse = 1e-6,
    W = as.numeric(parms["Winit"]),
    OM = 1e-6
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

# 2. Monte Carlo ----


f_MonteCarlo_DEB <- function(
    seed, Nb_iter, l_Experiments, 
    l_params_median, CV_parms, 
    Parms_fixed, Parms_unif, Limits_parms_unif
){
  
  set.seed(seed)
  Step_sim <- 0.1
  
  Nb_expe <- length(l_Experiments)
  
  Outputs <- c(
    "Energy",
    "Maturity",
    "Reproduction",
    "Weight",
    "Organic_matter"
  )
  
  l_parms <- f_MonteCarlo_l_parms(seed, Nb_iter, l_params_median, CV_parms, Parms_fixed, Parms_unif, Limits_parms_unif)
  
  df_MonteCarlo <- data.frame()
  
  for (expe_i in seq(1, Nb_expe)) {
    
    print(paste("Expérience :",expe_i))
    
    times_sim <- f_simulation_times(l_Experiments, expe_i, Step_sim)
    df_events <- f_events(l_Experiments, expe_i)
    forced_inputs <- f_forcing_inputs(l_Experiments, expe_i)
    
    for (i in seq(1,Nb_iter, 1)){
      
      l_parms_i <- l_parms[[i]]
      parms_deb_i <- f_initParms(l_parms_i, expe_i)
      Y_deb_i <- f_initStates(parms_deb_i)
      
      res <- tryCatch({
        
        # Résolution de l'ODE
        result_deb_i <- deSolve::ode(
          y = Y_deb_i,
          times_sim,
          func = "derivs",
          parms = parms_deb_i,
          rtol = 1e-7,
          atol = 1e-7,
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
        df_result_i <- as.data.frame(result_deb_i) |>
          pivot_longer(
            cols = -time,
            names_to = "Variable",
            values_to = "Value"
          ) |> 
          mutate(
            No_experiment = expe_i,
            Iter = i
          )
        
        # Ajout au dataframe final
        df_MonteCarlo <- rbind(df_MonteCarlo, df_result_i)
        
        df_MonteCarlo  # retour normal
        
      }, error = function(e) {
        message(sprintf("Erreur à l'expérience %d, itération %d : %s", expe_i, i, e$message))
        
        # On peut retourner NA ou df_MonteCarlo inchangé
        return(df_MonteCarlo)
      })
      
      print(i)
      
    } # End itération 
  } # End Expérience
  
  return(
    list(
      df_sim = df_MonteCarlo,
      df_parms = l_parms
      )
  )
}

f_RSS_calc <- function(df_data, df_MonteCarlo_sim, l_MonteCarlo_draws){
  
  df_data_simple <- df_data |> 
    dplyr::select(No_sim, Time, Reproduction, Weight) |> 
    mutate(
      No_experiment = as.numeric(No_sim),
      Time = as.numeric(Time)
      )
  
  df_MonteCarlo_sim_wide <- df_MonteCarlo_sim |> 
    filter(time %in% seq(0, 400, 0.5)) |> 
    mutate(
      Time = time,
    ) |> 
    distinct() |> 
    filter(Variable %in% c("Weight", "Reproduction")) |> 
    pivot_wider(
      id_cols = c(No_experiment, Time, Iter),
      names_from = Variable, 
      values_from = Value
    )
  
  df_params <- l_MonteCarlo_draws %>%
    map_dfr(as.list, .id = "Iter") %>%   # une ligne = un tirage
    mutate(Iter = as.integer(Iter))
  

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
      SSE_Weight = sum(err_weight, na.rm = TRUE),
      SSE_Reproduction = sum(err_repro, na.rm = TRUE),
      .by = Iter
    ) |> 
    left_join(df_params, by = "Iter")
  
  return(df_error)
}

