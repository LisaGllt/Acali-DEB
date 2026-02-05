
# General setup ----

f_load_libraries <- function(){
  # 📦 Data Manipulation
  library(tidyverse)   # Collection de packages pour manipulation et visualisation des données
  library(here)        # Gestion des chemins de fichiers
  library(readxl)      # Lecture des fichiers Excel
  library(reshape2)    # Restructuration des données
  library(DT)          # Génération de tables interactives
  library(knitr)       # Génération de rapports dynamiques en RMarkdown/Quarto
  library(stringr)
  library(glue)
  
  # 🎨 Visualization
  library(ggplot2)     # Visualisation de données
  library(ggthemes)    # Thèmes pour ggplot2
  library(ggdist)      # Distribution et incertitude
  library(ggsci)       # Palettes de couleurs scientifiques
  library(viridis)     # Palette de couleurs perceptuellement uniforme
  library(wesanderson) # Palette de couleurs artistiques
  library(RColorBrewer)# Palettes de couleurs prédéfinies
  library(nord)        # Palettes de couleurs inspirées du Nord
  library(plotly)      # Graphiques interactifs
  #library(ggiraph)     # Graphiques interactifs pour ggplot2
  library(ggrepel)     # Étiquettes non chevauchantes sur ggplot2
  library(patchwork)   # Combinaison de plusieurs ggplots
  library(gridExtra)   # Arrangements de graphiques en grille
  library(grid)        # Outils de mise en page graphique
  library(ggbreak)     # Briser les axes dans ggplot2
  library(ggtext)      # Formatage avancé de texte dans ggplot2
  library(kableExtra)  # Mise en forme avancée des tables
  library(flextable)
  library(gt) # Tables rendering
  library(processx)
  library(metR)
  library(rnaturalearth)
  library(sf)
  library(scales)
  library(colorspace)
  library(IDPmisc)
  library(ggforce)
  library(GGally)
  
  # 📊 Statistical Modeling & Bayesian Analysis
  library(brms)        # Modélisation bayésienne avec Stan
  library(rstan)       # Interface R pour Stan
  #library(cmdstanr)    # Interface alternative pour Stan (CmdStan)
  library(tidybayes)   # Manipulation et visualisation des résultats bayésiens
  library(ggmcmc)      # Diagnostics des chaînes MCMC
 # library(rethinking)  # Modélisation bayésienne avancée
  library(priorsense)  # Analyse de sensibilité des priors
  library(coda)
  library(ggmcmc)
  library(bayesnec)
  library(truncnorm)
  library(minpack.lm)
  library(nimble)
  
  
  # 🔬 Regression & Hypothesis Testing
  library(car)         # Tests statistiques et régressions avancées
  library(nlstools)    # Outils pour modèles non linéaires
  library(lsmeans)     # Comparaisons post-hoc
  library(ggpubr)      # Outils pour publications scientifiques
  library(marginaleffects) # Effets marginaux des modèles
  library(brglm2)      # Régressions logistiques biais-réduits
  library(multcomp)
  
  # ⚙️ Computational Tools & Parallelization
  library(parallel)    # Calcul parallèle
  library(deSolve)     # Équations différentielles
  library(tmvtnorm)    # Distribution normale tronquée multivariée
  library(fdrtool)     # Faux taux de découverte (FDR)
  library(drc)         # Modélisation de réponses aux doses
  
  # 🛠️ Model Evaluation & Performance
  library(easystats)   # Outils pour statistiques et modèles
  library(performance) # Diagnostics et évaluation de modèles
  library(modelsummary)# Résumé des modèles statistiques
  library(plan)        # Planification de l'exécution des tâches
  
  # 🖋️ Math & LaTeX Support
  library(latex2exp)   # Expressions LaTeX dans ggplot2
  library(extrafont)   # Gestion des polices pour ggplot2
  
}


f_load_colors <- function(){
  Nord_frost <<- nord(palette = "frost")
  Nord_polar <<- nord(palette = "polarnight")
  Nord_aurora <<- nord(palette = "aurora")
  sizetitle <<- 12
  
  col_blue <<- "#5E81AC"
  col_red <<- "#f42404"
  
  pal_blue <<- c("#5E81AC", "#7F9DC4", "#A0C1D9", "#DCE9F2")
  pal_red <<- c("#f42404", "#F65E4B", "#F6876D", "#FBD3D0")
  
  Nord_frost <<- nord(palette = "frost")
  Nord_aurora <<- nord(palette = "aurora")
  Nord_polar <<- nord(palette = "polarnight")
  Nord_snow <<- nord(palette = "snowstorm")
  
  pal_col <<- c(Nord_aurora[1], Nord_frost[4])
  col_EPX <<- Nord_frost[2]
  col_IMD <<- Nord_frost[4]
  col_elim <<- Nord_aurora[4]
  col_uptake <<- Nord_aurora[1]
  shape_IMD <<- 16
  shape_EPX <<- 15
  col_Molecule <<- c(col_EPX, col_IMD)
  sizetitle <<- 12
  
  shape_Molecule <<- c(shape_EPX, shape_IMD)
  
  set.seed(121212)
  
}

f_load_libraries_colors <- function(){
  f_load_libraries()
  f_load_colors()
}


# Add My Pet parameter recuperation ----

getDEB.species <- function() {
  require(rvest)
  library(rvest)
  url <- "https://www.bio.vu.nl/thb/deb/deblab/add_my_pet/species_list.html"
  d1 <- read_html(url)
  
  phylum <- d1 %>% html_nodes("td:nth-child(1)") %>% html_text()
  class <- d1 %>% html_nodes("td:nth-child(2)") %>% html_text()
  order <- d1 %>% html_nodes("td:nth-child(3)") %>% html_text()
  family <- d1 %>% html_nodes("td:nth-child(4)") %>% html_text()
  species <- d1 %>% html_nodes("td:nth-child(5)") %>% html_text()
  common <- d1 %>% html_nodes("td:nth-child(6)") %>% html_text()
  type <- d1 %>% html_nodes("td:nth-child(7)") %>% html_text()
  mre <- d1 %>% html_nodes("td:nth-child(8)") %>% html_text()
  smre <- d1 %>% html_nodes("td:nth-child(9)") %>% html_text()
  complete <- d1 %>% html_nodes("td:nth-child(10)") %>% html_text()
  all.species <- as.data.frame(cbind(phylum, class, order, 
                                     family, species, common, type, mre, smre, complete), 
                               stringsAsFactors = FALSE)
  all.species$species <- gsub(" ", "_", all.species$species)
  all.species$mre <- as.numeric(mre)
  all.species$smre <- as.numeric(smre)
  all.species$complete <- as.numeric(complete)
  return(all.species)
}

getDEB.pars <- function(species) {
  require(rvest)
  library(rvest)
  baseurl <- "https://www.bio.vu.nl/thb/deb/deblab/add_my_pet/entries_web/"
  d1 <- read_html(paste0(baseurl, species, "/", species, "_par.html"))
  symbol1 <- d1 %>% html_nodes("td:nth-child(1)") %>% html_text()
  
  value1 <- d1 %>% html_nodes("td:nth-child(2)") %>% html_text()
  
  units1 <- d1 %>% html_nodes("td:nth-child(3)") %>% html_text()
  
  description1 <- d1 %>% html_nodes("td:nth-child(4)") %>% 
    html_text()
  
  extra1 <- d1 %>% html_nodes("td:nth-child(5)") %>% html_text()
  
  extra2 <- d1 %>% html_nodes("td:nth-child(6)") %>% html_text()
  end <- which(symbol1 == "T_ref")
  symbol <- symbol1[1:end]
  value <- value1[1:end]
  units <- units1[1:end]
  description <- description1[1:end]
  
  pars <- as.data.frame(cbind(symbol, value, units, description))
  pars$symbol <- as.character(symbol)
  pars$value <- as.numeric(value)
  pars$units <- as.character(units)
  pars$description <- as.character(description)
  
  chempot <- c(value1[end + 1], units1[end + 1], description1[end + 
                                                                1], extra1[1])
  dens <- c(value1[end + 2], units1[end + 2], description1[end + 
                                                             2], extra1[2])
  org.C <- c(units1[end + 3], description1[end + 3], extra1[3], 
             extra2[1])
  org.H <- c(value1[end + 4], units1[end + 4], description1[end + 
                                                              4], extra1[4])
  org.O <- c(value1[end + 5], units1[end + 5], description1[end + 
                                                              5], extra1[5])
  org.N <- c(value1[end + 6], units1[end + 6], description1[end + 
                                                              6], extra1[6])
  min.C <- c(units1[end + 7], description1[end + 7], extra1[7], 
             extra2[2])
  min.H <- c(value1[end + 8], units1[end + 8], description1[end + 
                                                              8], extra1[8])
  min.O <- c(value1[end + 9], units1[end + 9], description1[end + 
                                                              9], extra1[9])
  min.N <- c(value1[end + 10], units1[end + 10], description1[end + 
                                                                10], extra1[10])
  
  organics <- rbind(org.C, org.H, org.O, org.N)
  minerals <- rbind(min.C, min.H, min.O, min.N)
  colnames(organics) <- c("X", "V", "E", "P")
  colnames(minerals) <- c("CO2", "H2O", "O2", "N-waste")
  rownames(organics) <- c("C", "H", "O", "N")
  rownames(minerals) <- c("C", "H", "O", "N")
  class(chempot) <- "numeric"
  class(dens) <- "numeric"
  class(organics) <- "numeric"
  class(minerals) <- "numeric"
  
  return(list(pars = pars, chempot = chempot, dens = dens, 
              organics = organics, minerals = minerals))
}

getDEB.implied <- function(species) {
  require(rvest)
  library(rvest)
  baseurl <- "https://www.bio.vu.nl/thb/deb/deblab/add_my_pet/entries_web/"
  d1 <- read_html(paste0(baseurl, species, "/", species, "_stat.html"))
  symbol <- d1 %>% html_nodes("td:nth-child(1)") %>% html_text()
  
  value <- d1 %>% html_nodes("td:nth-child(2)") %>% html_text()
  
  units <- d1 %>% html_nodes("td:nth-child(3)") %>% html_text()
  
  description <- d1 %>% html_nodes("td:nth-child(4)") %>% html_text()
  
  final <- as.data.frame(cbind(symbol, value, units, description))
  final$symbol <- as.character(symbol)
  final$value <- as.numeric(value)
  final$units <- as.character(units)
  final$description <- as.character(description)
  return(final)
}

f_GetParsAddMyPet <- function(){
  Species <- "Aporrectodea_caliginosa"
  allpars = getDEB.pars(Species)
  df_pars = allpars$pars
  df_pars = df_pars[-1,] # remove duplicate T_A symbol
  
  df_rename <- data.frame(
    symbol = c("p_Am", "F_m", "kap_X", "kap_P", "v", "kap", "kap_R", "p_M", "p_T", 
               "k_J", "E_G", "E_Hb", "E_Hp", "h_a", "s_G", "L0", "T_AH", "T_H",
               "Wwg", "bw", "del_M", "f", "f_Bart_high", "f_Bart_high_b", 
               "f_Bart_low", "f_Bart_low_b", "f_Bart_medium", "m_0", "max_r_mb",
               "mu_OM", "mu_c", "r_mb", "t_0", "wV", "T_A", "T_ref", "Wd_0"),
    name =   c("pAm",  "Fm",  "kapX",  "kapP",  "v", "kap", "kapR",  "pM",  "pT",  
               "kJ",  "Eg",  "Ehb",  "Ehp",  "ha",  "sG",  "L0", "TAH",  "TH",
               "Wwg", "bw", "Shape", "f", "fBartHigh", "fBartHighb", 
               "fBartLow",   "fBartLowb",    "fBartMedium",   "m0",  "maxRmb", 
               "muOM",  "muC",  "Rmb",  "t0",  "wV", "TA",  "Tref",  "Wd0")
  )
  
  df_pars <- df_pars |> left_join(df_rename)
  
  df_pars = rbind(df_pars, data.frame(symbol = "UE0",
                                      value = 0.078, # from Gergs et al., 2022
                                      units = "cm^2/d",
                                      description = "Scaled cost of an egg",
                                      name = "UE0"))
  df_pars = rbind(df_pars, data.frame(symbol = "E0",
                                      value = 134.988, # from AddMyPet website
                                      units = "J",
                                      description = "Scaled cost of an egg",
                                      name = "E0"))
  df_pars = rbind(df_pars, data.frame(symbol = "w",
                                      value = 27.24, # from Gergs et al., 2022
                                      units = "g/cm^3",
                                      description = "Contribution of reserve to bw",
                                      name = "w"))
  
  return(df_pars)
}

f_AddMyPet_Tcorrection <- function(df_pars, Texp){
  
  # Correction des valeurs par la température (Gergs et al. 2022)

  df_pars$corr_value <- df_pars$value

  TA   <- subset(df_pars, name =="TA")$corr_value
  TAH  <- subset(df_pars, name =="TAH")$corr_value
  TH   <- subset(df_pars, name =="TH")$corr_value
  Tref <- subset(df_pars, name =="Tref")$corr_value

  sA   <- exp(TA/Tref-TA/(Texp+273.15))
  srH  <- (1 + exp(TAH/TH-TAH/Tref))/(1+exp(TAH/TH-TAH/(Texp+273.15)))
  Ft   <- sA*((Texp+273.15>=Tref)*srH + (Texp+273.15<Tref))

  df_pars_add <- df_pars |>
    filter(name %in% c("pAm", "pM", "v", "kJ")) |>
    mutate(corr_value = value*Ft)|>
    mutate(
      symbol = paste(symbol,"_t", sep=""),
      name = paste(name,"_t", sep="")
      )

  df_pars <- rbind(df_pars, df_pars_add)
  
  return(df_pars)
}

f_vec_pars <- function(df_pars){
  # Liste des paramètres pour le modèle
  vec_pars = df_pars$corr_value
  names(vec_pars) = df_pars$name
  vec_pars = append(vec_pars, c(kapH = 1)) # set the maturation efficiency to 1
  return(vec_pars)
}

f_GetParsAddMyPet_asglobal <- function(df_pars){
  
  pAm      <<- subset(df_pars, name =="pAm")$corr_value
  pAm_t    <<- subset(df_pars, name =="pAm_t")$corr_value
  Fm       <<- subset(df_pars, name =="Fm")$corr_value
  kapX     <<- subset(df_pars, name =="kapX")$corr_value
  kapP     <<- subset(df_pars, name =="kapP")$corr_value
  v        <<- subset(df_pars, name =="v")$corr_value
  v_t      <<- subset(df_pars, name =="v_t")$corr_value
  kap      <<- subset(df_pars, name =="kap")$corr_value
  kapR     <<- subset(df_pars, name =="kapR")$corr_value
  pM       <<- subset(df_pars, name =="pM")$corr_value
  pM_t     <<- subset(df_pars, name =="pM_t")$corr_value
  pT       <<- subset(df_pars, name =="pT")$corr_value
  kJ       <<- subset(df_pars, name =="kJ")$corr_value
  Eg       <<- subset(df_pars, name =="Eg")$corr_value
  Ehb      <<- subset(df_pars, name =="Ehb")$corr_value
  Ehp      <<- subset(df_pars, name =="Ehp")$corr_value
  ha       <<- subset(df_pars, name =="ha")$corr_value
  sG       <<- subset(df_pars, name =="sG")$corr_value
  L0       <<- subset(df_pars, name =="L0")$corr_value
  Wwg      <<- subset(df_pars, name =="Wwg")$corr_value
  bw       <<- subset(df_pars, name =="bw")$corr_value
  Shape    <<- subset(df_pars, name =="Shape")$corr_value
  m0       <<- subset(df_pars, name =="m0")$corr_value
  maxRmb   <<- subset(df_pars, name =="maxRmb")$corr_value
  muOM     <<- subset(df_pars, name =="muOM")$corr_value
  muC      <<- subset(df_pars, name =="muC")$corr_value
  Rmb      <<- subset(df_pars, name =="Rmb")$corr_value
  t0       <<- subset(df_pars, name =="t0")$corr_value
  wV       <<- subset(df_pars, name =="wV")$corr_value
  Wd0      <<- subset(df_pars, name =="Wd0")$corr_value
  E0       <<- subset(df_pars, name =="E0")$corr_value
  w        <<- subset(df_pars, name =="w")$corr_value
  
  TH        <<- subset(df_pars, name =="TH")$corr_value
  TAH        <<- subset(df_pars, name =="TAH")$corr_value
  TA        <<- subset(df_pars, name =="TA")$corr_value
  Tref        <<- subset(df_pars, name =="Tref")$corr_value
}

