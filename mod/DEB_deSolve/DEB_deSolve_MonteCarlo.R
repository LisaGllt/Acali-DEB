library(here)
source(file = here::here("functions/fun.R"))
source(file = here::here("functions/fun_MCSim.R"))
f_load_libraries_colors()
col_Molecule <- rev(col_Molecule)
pal_chains <- c(Nord_aurora, Nord_frost, Nord_polar)

# 1. Compilation du modèle ----
Model.name <- "DEBcali"

mod_bin <- "/Users/lisagollot/mcsim-6.2.0/mod/mod"
mod_dir <- "/Users/lisagollot/Library/CloudStorage/OneDrive-Personnel/Documents/0_These/0_RepoGit/ew-deb-tktd/mod/DEB_deSolve"
mod_file <- "DEBcali.model"
c_file <- "DEBcali.c"

system(
  paste(
    mod_bin, "-R",
    file.path(mod_dir, mod_file),
    file.path(mod_dir, c_file)
  )
)

# A ajouter en haut du fichier .c : #include <math.h>

getLoadedDLLs()[["DEBcali"]]
dyn.unload(file.path(mod_dir, paste0(Model.name, .Platform$dynlib.ext))) # A faire sinon ça update pas bien
system(paste("R CMD SHLIB", file.path(mod_dir, paste0(Model.name, ".c"))))
dyn.load(file.path(mod_dir, paste0(Model.name, .Platform$dynlib.ext)))


source(file = here::here("mod/DEB_deSolve/DEBcali_inits_MonteCarlo.R"))

# 2. Data ----

l_Experiments <- c(
  "Bart2019_XP1_1",  
  "Bart2019_XP1_2",  
  #"Bart2019_XP2",
  "Bart2019_XP3",   
  "Bart2019_XP4_1",
  "Bart2019_XP4_2",
  "Bart2019_XP5", 
  "Bart2020",
  #"Gollot2026_lufa"
  "Gollot2026_dens_D1",
  "Gollot2026_dens_D2",
  "Gollot2026_dens_D2b",
  "Gollot2026_dens_D5",
  "Gollot2026_dens_D7"
)
Nb_expe <- length(l_Experiments)

df_data <- f_import_data_DEB(l_Experiments, TRUE)

df_data_long <- df_data |>
  dplyr::select(c(Time, Weight, Reproduction, No_sim)) |> 
  pivot_longer(
    cols = -c("Time", "No_sim"),
    names_to = "Variable",
    values_to = "Value"
  ) |> 
  mutate(No_experiment = No_sim)

# 3. Monte Carlo ----

l_params_median <- c(
  # Environment
  muOM = 11700, rOM_ClxHorse=0.3,
  # Energy assimilation
  Fm = 0.5, kapX = 0.258,  pAm = 1712.41, v = 0.018505, kap  = 0.4373,  
  # Maintenance and growth
  pM = 1680.1, pT = 0, Eg = 6880.33, Shape = 0.066193, w = 27.24,        
  # Maturity
  kJ = 0.002793, Ehb = 0.5, Ehp = 100,
  # Reproduction
  L_coc = 0.23, E_coc = 467, kapR = 0.475,   
  # Metabolic response to T
  TAH = 28750, TH = 293.2, TA = 7976, Tref = 293.15   
)

Parms_fixed <- c("kapX", "Eg","Shape","kapR", "TAH", "TH", "TA", "Tref")
Parms_unif <- list("rOM_ClxHorse", "Fm", "kap")
Limits_parms_unif <- list(
  rOM_ClxHorse = c(0,1),
  Fm = c(0.2, 0.5),
  kap = c(0.1, 0.9)
)

seed <- 121212
Nb_iter <- 100
CV_parms <- 0.3

l_MonteCarlo_tot <- f_MonteCarlo_DEB(
    seed, Nb_iter, l_Experiments, 
    l_params_median, CV_parms, 
    Parms_fixed, Parms_unif, Limits_parms_unif
)

df_succes <- df_MonteCarlo %>%
  group_by(No_experiment) %>%
  summarise(
    n_Iter = n_distinct(Iter),
    .groups = "drop"
  )

saveRDS(l_MonteCarlo_tot, "mod/DEB_deSolve/MonteCarlo_deSolve_20260202_100iter2.rds")
df_MonteCarlo <- readRDS("mod/DEB_deSolve/MonteCarlo_deSolve_20260202_100iter2.rds")
df_MonteCarlo <- l_MonteCarlo_tot$df_sim
l_MonteCarlo_draws <- l_MonteCarlo_tot$df_parms

# 4. Selection des meilleures combinaisons ----

df_errors <- f_RSS_calc(df_data, df_MonteCarlo_sim, l_MonteCarlo_draws)

df_errors_best10_SSE_Weight <- (df_errors |> 
  arrange(SSE_Weight))[1:10,]
df_errors_best10_SSE_Reproduction <- (df_errors |> 
  arrange(SSE_Reproduction))[1:10,]
  

df_errors_long <- df_errors |> 
  pivot_longer(
    cols = -c(Iter, SSE_Weight, SSE_Reproduction),
    names_to = "Parameter", 
    values_to = "Value"
  )

df_errors_best10_SSE_Weight_long <- df_errors_best10_SSE_Weight |> 
  pivot_longer(
    cols = -c(Iter, SSE_Weight, SSE_Reproduction),
    names_to = "Parameter", 
    values_to = "Value"
  )

df_errors_best10_SSE_Reproduction_long <- df_errors_best10_SSE_Reproduction |> 
  pivot_longer(
    cols = -c(Iter, SSE_Weight, SSE_Reproduction),
    names_to = "Parameter", 
    values_to = "Value"
  )

plot_Param <- c("E_coc", "Ehb", "Ehp", "Fm", "kap", "kJ", "L_coc", "muOM", "pAm", "pM", "rOM_ClxHorse", "v", "w")
#plot_Param <- c("Fm")

p_SSE_weight <- ggplot()+
  geom_point(
    data = df_errors_long |> filter(Parameter %in% plot_Param), 
    aes(
      x = Value, 
      y = SSE_Weight
    ),
    color = Nord_polar[4],
    alpha = 0.7
  )+
  geom_point(
    data = df_errors_best10_SSE_Weight_long |> filter(Parameter %in% plot_Param), 
    aes(
      x = Value, 
      y = SSE_Weight
    ),
    color = Nord_aurora[4]
  )+
  facet_wrap(~Parameter, scales = "free")+
  theme_minimal()
p_SSE_weight

p_SSE_repro <- ggplot()+
  geom_point(
    data = df_errors_long |> filter(Parameter %in% plot_Param), 
    aes(
      x = Value, 
      y = SSE_Reproduction
    ),
    color = Nord_polar[4],
    alpha = 0.7
  )+
  geom_point(
    data = df_errors_best10_SSE_Reproduction_long |> filter(Parameter %in% plot_Param), 
    aes(
      x = Value, 
      y = SSE_Reproduction
    ),
    color = Nord_aurora[4]
  )+
  facet_wrap(~Parameter, scales = "free")+
  theme_minimal()
p_SSE_repro

p_SSE_weight+p_SSE_repro


l_Best_SSE_weight <- unique(df_errors_best10_SSE_Weight$Iter)
l_Best_SSE_repro <- unique(df_errors_best10_SSE_Reproduction$Iter)
l_Best_SSE <- c(l_Best_SSE_weight, l_Best_SSE_repro)

Outputs_print <- c("Energy", "Maturity", "Weight", "Reproduction", "Organic_matter")

l_Draws <- sample(seq(1,100,1), 10, replace = FALSE)

df_deb_sim <- df_MonteCarlo |>
  filter(Variable %in% Outputs_print) |> 
  filter(time %in% seq(0,400, 0.5))




df_plot <- df_deb_sim |> 
  filter(Iter %in% c(l_Best_SSE))

xmax <- 365
p <- ggplot()+
  geom_point(
    data = df_data_long,
    mapping = aes(
      x = Time,
      y = Value,
      group = Variable
    ),
    color = Nord_polar[4],
    alpha = 0.7,
    size = 0.5
  )+
  geom_line(
    data = df_plot,
    mapping = aes(
      x = time,
      y = Value,
      group = Iter,
      color = Variable
      )
    )+
  
  facet_grid(rows = vars(Variable), cols = vars(No_experiment), scales = "free")+
  #facet_wrap(~Variable+No_experiment, scales = "free")+
  labs(x = "Time", y = "Value", title = "Simulation DEB") +
  scale_color_manual(values = c(Nord_aurora, Nord_frost))+
  scale_x_continuous(
    breaks = seq(0, xmax, by = 28),
    minor_breaks = seq(0, xmax, by = 14)
  ) +
  theme_minimal()+
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 6)
  )
p

# 
# 
# 
# 
# {
#   p <- ggplot(
#     df_deb_sim,
#     aes(
#       x = time,
#       y = Value,
#       color = Variable
#     )
#   ) +
#     geom_line(size = 1) +
#     facet_wrap(
#       ~Variable,
#       scales = "free"
#     ) +
#     scale_x_continuous(
#       limits = c(0, xmax),       
#       breaks = seq(0, xmax, by = 28),
#       minor_breaks = seq(0, xmax, by = 14) 
#     ) +
#     labs(x = "Time", y = "Value", title = "Simulation DEB") +
#     theme_minimal()
#   p
# }
# 
# 
