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

source(file = here::here("mod/DEB_deSolve/DEBcali_inits_AS.R"))

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

df_data <- f_import_data_DEB(l_Experiments, Adults_alone = TRUE)

# 3. AS ----

Sobol_type <- "Jansen"
Nb_sim     <- 1000
processors <- 2
Nboot      <- 100 # nb replicat
NOut       <- 3    # number of outputs
Temperature_bol <- TRUE

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

No_experience <- which(l_Experiments == "Gollot2026_dens_D2")

l_Sobol_tot <- f_AS_DEB(Sobol_type, l_Experiments, No_experience, Nb_sim, Processors, Nboot, Nout, l_params_median, Temperature_bol)
df_Sobol_results <- l_Sobol_tot$df_res
Nb_fail <- l_Sobol_tot$Nb_fail
save(df_Sobol_results, file = here::here(paste0("mod/DEB_deSolve/Resultat_AS_Sobol_",Sobol_type, Nb_sim, "F", Nb_fail,"_res.RData")))

# 4. Plots ----


Nb_sim <- 1000
Sobol_type <- "Jansen"

load(here::here(paste0("mod/DEB_deSolve/Resultat_AS_Sobol_Jansen1000F_res.RData")))

path_fig <- here::here("fig/AS")

## Plot Indice Totaux empilés ----
p_TI_tot <- ggplot(
  df_Sobol_results, 
  aes(
    x = reorder(Parameter, TI), 
    y = TI
    )
  ) +
  geom_col(
    aes(
      fill = Variable
      ),
    alpha = 0.7
    ) +
  coord_flip() +
  scale_fill_manual(
    values = rev(c(Nord_aurora, rev(Nord_frost))),
    name = "Endpoints"
    )+
  labs(
    title = "Sobol Sensitivity Analysis - Total Index",
    x = "Parameter",
    y = "Total Index",
    caption = paste0("n = ", Nb_sim)
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
    )

p_TI_tot

ggsave(
  filename = paste0("AS_Sobol",Sobol_type, Nb_sim,"_TotalIndex_perparam.png"),
  plot = p_TI_tot,
  width = 8, height = 8,
  path = path_fig
  )

## Plot Indice de sensibilité totaux par variable ----
{ 
  Color_AS <- rev(c(Nord_aurora, rev(Nord_frost)))
  size_title_facet <- 10
  size_axis_text_y <- 7
  for (i in 1:length(unique(df_Sobol_results$Variable))){
    
    Variable_i <- unique(df_Sobol_results$Variable)[i]
    
    if (i == 1){
      p1 <- ggplot(
        subset(df_Sobol_results, Variable == Variable_i), 
        aes(
          x = reorder(Parameter, TI), 
          y = TI
        )
      ) +
        geom_col(
          fill = Color_AS[i],
          alpha = 0.7
        ) +
        coord_flip() +
        labs(
          title = as.character(Variable_i),
          x = "Parameter",
          y = "Total Index"#,
          #caption = "Error bars = CI 95% | n = 1000"
        ) +
        theme_minimal() +
        theme(
          legend.position = "right",
          axis.text.y = element_text(size = size_axis_text_y),
          plot.title = element_text(face = "bold", size = size_title_facet),
          legend.title = element_text(face = "bold")
        )
    } else {
      p <- ggplot(
        subset(df_Sobol_results, Variable == Variable_i), 
        aes(
          x = reorder(Parameter, TI), 
          y = TI
        )
      ) +
        geom_col(
          fill = Color_AS[i],
          alpha = 0.7
        ) +
        coord_flip() +
        labs(
          title = as.character(Variable_i),
          x = "Parameter",
          y = "Total Index"#,
          #caption = "Error bars = CI 95% | n = 1000"
        ) +
        theme_minimal() +
        theme(
          legend.position = "right",
          axis.text.y = element_text(size = size_axis_text_y),
          plot.title = element_text(face = "bold", size = size_title_facet),
          legend.title = element_text(face = "bold")
        )
      p1 <- p1 + p
    }
  }
  
  p_TI_pervar <- p1 + plot_layout(ncol = 3) + 
    plot_annotation(
      title = "Sobol Sensitivity Analysis - Total Index",
      caption = paste0("n = ", Nb_sim),
      theme = theme(
        plot.title = element_text(size = 16, face = "bold"),
      )
    )
  p_TI_pervar
}

ggsave(
  filename = paste0("AS_Sobol",Sobol_type, Nb_sim,"TotalIndex_pervariable.png"),
  plot = p_TI_pervar,
  width = 12, height = 10,
  path = path_fig
)

## FOI vs. TI ----

df_Sobol_results_long <- df_Sobol_results %>%
  pivot_longer(
    cols = c(FOI, TI),
    names_to = "Index_Type",
    values_to = "Index_Value"
  ) %>%
  mutate(
    # Ajouter les IC correspondants
    IC_inf = ifelse(Index_Type == "FOI", FOI.borninf, TI.borninf),
    IC_sup = ifelse(Index_Type == "FOI", FOI.bornsup, TI.bornsup)
  )


# Graphique comparatif
p2 <- ggplot(
  df_Sobol_results_long, 
  aes(
    x = reorder(Parameter, Index_Value), 
    y = Index_Value, 
    fill = Index_Type)
  ) +
  geom_col(
    position = "dodge", 
    alpha = 0.8
    ) +
  geom_errorbar(
    aes(
      ymin = IC_inf, 
      ymax = IC_sup
      ),
    position = position_dodge(width = 0.8),
    width = 0.25,
    linewidth = 0.4,
    color = Nord_polar[4]
  ) +
  coord_flip() +
  facet_wrap(
    ~Variable, 
    ncol = 4,
    scales = "free"
    ) +
  scale_fill_manual(
    name = "Type d'indice",
    values = c("FOI" = "#5E81AC", "TI" = "#BF616A"),
    labels = c("FOI" = "Premier ordre", "TI" = "Total")
  ) +
  labs(
    title = "Comparaison des indices de sensibilité",
    x = "Paramètre",
    y = "Valeur de l'indice",
    caption = paste0("Error bars = CI 95% | n = ", Nb_sim)
  ) +
  theme_minimal()

p2

ggsave(
  filename = paste0("AS_Sobol",Sobol_type, Nb_sim,"TotalIndex_FOI.png"),
  plot = p2,
  width = 12, height = 10,
  path = path_fig
)

## Heatmap ----

p5 <- ggplot(
  df_Sobol_results, 
  aes(
    x = Variable, 
    y = Parameter, 
    fill = TI
    )
  ) +
  geom_tile(
    color = "black"
    ) +
  geom_text(
    aes(
      label = round(TI, 2)), 
    color = "white", 
    size = 3
    ) +
  scale_fill_gradient2(
    low = "#ECEFF4",
    mid = "#5E81AC",
    high = "#BF616A",
    midpoint = 0.5,
    name = "Indice total"
  ) +
  labs(
    title = "Total Index Heatmap",
    x = "Variable",
    y = "Parameter",
    caption = paste0("n = ", Nb_sim)
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p5

ggsave(
  filename = paste0("AS_Sobol",Sobol_type, Nb_sim,"TotalIndex_HeatMap.png"),
  plot = p5,
  width = 12, height = 10,
  path = path_fig
)

## Interactions ----

p6 <- ggplot(
  df_Sobol_results, 
  aes(
    x = FOI, 
    y = TI
    )
  ) +
  geom_abline(
    slope = 1, 
    intercept = 0, 
    linetype = "dashed", 
    color = "gray50"
    ) +
  geom_point(
    aes(
      color = Variable
      ), 
    size = 3, 
    alpha = 0.7
    ) +
  geom_text_repel(
    aes(
      label = Parameter
      ), 
    size = 3
    ) +
  facet_wrap(~Variable)+
  scale_color_manual(
    values = rev(c(Nord_aurora, rev(Nord_frost))),
    name = "Endpoints"
  )+
  labs(
    title = "Premier ordre vs Total",
    x = "Indice de premier ordre (FOI)",
    y = "Indice total (TI)",
    caption = "Points au-dessus de la diagonale = effets d'interaction importants"
  ) +
  theme_minimal() +
  coord_fixed()

p6



p3 <- ggplot(
  df_Sobol_results, 
  aes(
    x = TI, 
    y = reorder(Parameter, TI))
) +
  geom_errorbarh(
    aes(
      xmin = TI.borninf, 
      xmax = TI.bornsup
    ),
    height = 0.2,
    linewidth = 0.8,
    color = Nord_polar[4]
  ) +
  geom_point(
    size = 2, 
    color = "#5E81AC"
  ) +
  geom_vline(
    xintercept = 0, 
    linetype = "dashed", 
    color = "gray50"
  ) +
  facet_wrap(
    ~Variable, 
    scales = "free_y"
  ) +
  labs(
    title = "Indices de sensibilité totaux avec IC 95%",
    x = "Indice de sensibilité total",
    y = "Paramètre"
  ) +
  theme_minimal()

p3

## Radar plot ----

library(fmsb)

# Pour une seule variable
df_radar <- df_Sobol_results %>%
  filter(Variable == unique(df_Sobol_results$Variable)[1]) %>%
  dplyr::select(Parameter, TI) %>%
  pivot_wider(names_from = Parameter, values_from = TI)

# Ajouter les valeurs min/max pour fmsb
df_radar <- rbind(
  rep(1, ncol(df_radar)),  # Max
  rep(0, ncol(df_radar)),  # Min
  df_radar
)

# Tracer
radarchart(
  df_radar,
  axistype = 1,
  pcol = "#5E81AC",
  pfcol = scales::alpha("#5E81AC", 0.3),
  plwd = 2,
  cglcol = "grey",
  cglty = 1,
  axislabcol = "grey",
  caxislabels = seq(0, 1, 0.25),
  title = paste("Indices de sensibilité totaux -", unique(df_Sobol_results$Variable)[1])
)


