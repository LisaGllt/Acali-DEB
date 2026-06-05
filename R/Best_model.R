library(here)
source(file = here::here("functions/fun.R"))
source(file = here::here("functions/fun_MCSim.R"))
f_load_libraries_colors()
col_Molecule <- rev(col_Molecule)
pal_chains <- c(Nord_aurora, Nord_frost, Nord_polar)
lab_chains <- c("1", "2", "3", "4", "5", "6", "7", "8")


path_D_CWKe <- here::here("mod/W_D_DEB_CWKe")

l_Experiments_fit <- c(
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

Adults_alone <- FALSE

df_data_fit <- f_import_data_DEB(l_Experiments_fit, Adults_alone) 

Experiments_carac <- f_get_experiments()
l_Experiments_tot <- Experiments_carac$l_Experiments_tot
l_carac <- Experiments_carac$l_carac
df_carc_exp <- Experiments_carac$df_carc_exp


NbIter <- 100000
OM_diff <- TRUE
RTOL <- 1e-7
ATOL <- 1e-6

n_data_Weight <- length(na.omit(df_data_fit$Weight))
n_data_Reproduction <- length(na.omit(df_data_fit$Reproduction))
n_data_fit <- n_data_Weight + n_data_Reproduction


MCMC_out_D_CWKe <- f_MCSim(path_D_CWKe)

NC_D_CWKe <- MCMC_out_D_CWKe$NC
Nb_param_D_CWKe <- MCMC_out_D_CWKe$Nb_param
Nb_experiment_D_CWKe <- MCMC_out_D_CWKe$Nb_experiment
Nb_iter_kept_D_CWKe <- MCMC_out_D_CWKe$Nb_iter_kept
N_iter_setpoint_D_CWKe <- MCMC_out_D_CWKe$N_iter_setpoint
L_parameters_D_CWKe <- unique(MCMC_out_D_CWKe$Priors$Nom)

LL_D_CWKe <- f_get_mode(MCMC_out_D_CWKe$df_LnData$LnData)
k_D_CWKe <- Nb_param_D_CWKe - 1 # On enlève le Sigma_W

BIC_D_CWKe <- -2*LL_D_CWKe+k_D_CWKe*log(n_data_fit)

Summary_res_D_CWKe <- 
  as.data.frame(MCMC_out_D_CWKe$Summary_res)[1:length(L_parameters_D_CWKe)] |> 
  mutate(Model = "D_CWKe") 



# Chains ----

p_chains <-  
  ggplot(
    data = MCMC_out_D_CWKe$Chains |> 
      filter(Iteration !=1), 
    aes(
      x = Iteration, 
      y = value, 
      color = as.factor(Chain), 
      group = as.factor(Chain)
    )
  )+
  geom_line(alpha=0.7)+
  scale_y_log10()+
  facet_wrap(~Parameter, scales = "free", ncol = 2)+
  
  scale_color_manual(name = "Chains", values = pal_chains)+
  
  theme_bw()+
  theme(
    legend.position = "right", 
    title=element_text(size=12, face="plain"), 
    axis.title.x = element_text(face="plain"),
    strip.background = element_rect(fill="white")
  )


p_chains



## Correlations ----

# Large format per chain
mcmc_list_D_CWKe <- lapply(
  split(
    MCMC_out_D_CWKe$Chains, 
    MCMC_out_D_CWKe$Chains$Chain
  ),
  function(d) {
    wide <- d %>%
      dplyr::select(Iteration, Parameter, value) %>%
      pivot_wider(names_from = Parameter, values_from = value) %>%
      arrange(Iteration)
    mcmc(as.matrix(wide[,-1]))
  })

mcmc_list_D_CWKe <- as.mcmc.list(mcmc_list_D_CWKe)

mcmc_list_corr_D_CWKe <- mcmc.list(
  lapply(mcmc_list_D_CWKe, function(x) {
    mcmc(as.matrix(x)[-1, , drop = FALSE])  # <- reconvertir en mcmc
  })
)

colramp_aurora <- colorRampPalette(rev(Nord_aurora[1:4]))

ipairs(
  as.matrix(mcmc_list_corr_D_CWKe), 
  colramp = colramp_aurora, 
  pixs=0.25, 
  cex.diag = 0.5
)

Nb_experiment <- 12
step_sim <- 0.05
select_times <- seq(0,400, 0.1)

Endpoints_print <- "Weight,
    Reproduction"


l_param_name_D_CWKe <- names(Summary_res_D_CWKe) |>
  head(-1) |>
  stringr::str_remove_all("\\.1\\.")

f_In_Setpoint_full(
  path_D_CWKe, l_Experiments_fit, 
  Adults_alone, OM_diff, 
  l_param_name_D_CWKe, step_sim, Endpoints_print
)

file_rds_D_CWKe <- file.path(path_D_CWKe, paste0("Sim.Res.Exp.full.rds"))

if (file.exists(file_rds_D_CWKe)) {
  df_pred_full_D_CWKe <- readRDS(file_rds_D_CWKe)
} else {
  df_pred_full_D_CWKe <- f_read_Setpoint_full(path_D_CWKe, l_Experiments_fit, Adults_alone = FALSE, Nb_experiment, select_times, step_sim)
}

Nb_experiment <- 12

df_data_obs <- df_data_fit |> 
  mutate(
    Weight_obs = Weight,
    Reproduction_obs = Reproduction
  ) |> 
  dplyr::select(-c(Weight, Reproduction))


# Plot
alpha_ribbon <- 0.15
line_width <- 1.1
point_stroke <- 0.9
alpha_line <- 0.7
wrap_nb_col <- 2

col_experiment <- c(
  Nord_aurora[1], Nord_aurora[1], Nord_polar[1],
  Nord_aurora[2],
  Nord_aurora[2], Nord_aurora[3],
  Nord_aurora[4], 
  Nord_aurora[5], rev(Nord_frost)
)

pW <- ggplot()+
  
  geom_ribbon(
    data = subset(df_pred_full_D_CWKe, Endpt == "Weight"),
    aes(
      x = Time,
      ymax = up,
      ymin = low,
      fill = ID_experiment
    ),
    alpha = alpha_ribbon
  )+
  
  geom_line(
    data = subset(df_pred_full_D_CWKe, Endpt == "Weight"),
    aes(
      x = Time,
      y = predict.endpoint,
      color = ID_experiment
    ),
    linewidth = line_width,
    alpha = alpha_line
  )+
  
  geom_point(
    data = df_data_obs,    
    mapping = 
      aes(  
        x = Time, 
        y = Weight_obs, 
        color = ID_experiment  
      ),
    alpha = 0.7,
    stroke = point_stroke
  )+
  
  scale_color_manual(
    name = "Experiments",
    values = col_experiment
  )+
  scale_fill_manual(
    name = "Experiments",
    values = col_experiment
  )+
  facet_wrap(
    ~ID_experiment,
    ncol = wrap_nb_col,
    scales = "free"
  )+
  labs(
    x="Time (days)", 
    y = "Weight (g)",
    title ="Earthworm growth"
  )+
  theme_minimal()+
  theme(
    plot.title = element_text(face = "bold")
  )

pR <- ggplot()+
  
  geom_ribbon(
    data = subset(df_pred_full_D_CWKe, Endpt == "Reproduction"),
    aes(
      x = Time,
      ymax = up,
      ymin = low,
      fill = ID_experiment
    ),
    alpha = alpha_ribbon
  )+
  
  geom_line(
    data = subset(df_pred_full_D_CWKe, Endpt == "Reproduction"),
    aes(
      x = Time,
      y = predict.endpoint,
      color = ID_experiment
    ),
    linewidth = line_width,
    alpha = alpha_line
  )+
  
  geom_point(
    data = df_data_obs,    
    mapping = 
      aes(  
        x = Time, 
        y = Reproduction_obs, 
        color = ID_experiment  
      ),
    alpha = 0.8,
    stroke = point_stroke
  )+
  
  scale_color_manual(
    name = "Experiments",
    values = col_experiment
  )+
  scale_fill_manual(
    name = "Experiments",
    values = col_experiment
  )+
  facet_wrap(
    ~ID_experiment, 
    ncol=wrap_nb_col,
    scales = "free"
  )+
  labs(
    x="Time (days)", 
    y = "Reproduction (#)",
    title ="Earthworm reproduction"
  )+
  theme_minimal()+
  theme(
    plot.title = element_text(face = "bold")
  )

p <- pW + pR +
  plot_layout(guides = "collect") &
  guides(color = guide_legend("Experiments")) &
  theme(legend.position = "none")

p
