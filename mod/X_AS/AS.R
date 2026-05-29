library(here)
source(here::here("functions/fun.R"))
source(here::here("functions/fun_MCSim.R"))
f_load_libraries_colors()


# #. Paramètre de l'analyse de sensibilité ----

seed       <- 361224
Sobol_type <- "Jansen"
Nb_sim     <- 50000 # (p + 2) * Nb_sim simulations
CV         <- 0.1
processors <- 2
Nboot      <- 1000 # nb replicat
NOut       <- 3    # number of outputs
Temperature_bol <- TRUE

l_params_median <- c(
  # Environment
  muOM = 11700, rOM_ClxHorse = 0.0087741,
  # Energy assimilation
  Fm = 0.0006748850, kapX = 0.258, K = 0.00611915, pAm = 160.85200, v = 0.212075, 
  kapJ  = 0.913533, kapA = 0.494191, tau = 12.9371,
  # Maintenance and growth
  r_pAm_pM = 0.75549700, Eg = 4183, dv = 4.47722, de = 1.849940e-09,         
  # Maturity
  kJ = 0.002793, Ehb = 0.5, Ehp = 97.5249,
  # Reproduction
  E_coc = 27.4525, kapR = 0.475,   
  # Metabolic response to T
  TAH = 28750, TH = 293.2, TA = 7976, Tref = 293.15   
)

Sobol_design <- f_AS_Sobol_design(Sobol_type, seed, Nb_sim, CV, Processors, Nboot, Nout, l_params_median, Temperature_bol)

l_Sobol_param <- Sobol_design$Sobol_param
Sobol_info <- Sobol_design$Sobol

df_Sobol_param <- do.call(rbind, lapply(l_Sobol_param, as.data.frame)) 

df_Sobol_param <- cbind(id = seq_len(nrow(df_Sobol_param)), df_Sobol_param)

write.table(df_Sobol_param, file = paste0("mod/FINAL/Z_AS/Param_sobol_", Sobol_type,".txt"), sep = "\t", row.names = FALSE, quote = FALSE)


#./mcsim.DEBcali DEB_AS.in

# 2. Index calculations ----

Names_AS_outputs <- c("Weight_Juv", "Weight_Adult", "Energy_Juv", "Energy_Adult", "Maturity_Juv", "Maturity_Adult", "Reproduction", "Organic_matter")

df_AS_sim <- read_tsv(file.path(paste0("mod/Correction_X/Z_AS/Sobol_", Sobol_type,".out"))) |> 
  dplyr::select(-(1:(length(l_params_median)+1-(Temperature_bol)*4))) |> 
  setNames(Names_AS_outputs)

df_Sobol_results <- f_AS_Sobol_index(df_AS_sim, Sobol_info)

# 3. Plots ----

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
  scale_y_log10()+
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

  