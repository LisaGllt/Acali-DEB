library(here)
source(file = here::here("functions/fun.R"))
source(file = here::here("functions/fun_MCSim.R"))
f_load_libraries_colors()
col_Molecule <- rev(col_Molecule)
pal_chains <- c(Nord_aurora, Nord_frost, Nord_polar)
lab_chains <- c("1", "2", "3", "4", "5", "6", "7", "8")



Molec <- "IMD"
path_fig <- here::here("fig/DEBTKTD_Models")

path_E_IMD_A <- here::here("mod/Z_DEB-TKTD/IMD/Test_IMD")

df_data_fit_IMD <- f_import_data_DEBTKTD() |> 
  filter(Molecule == Molec) |> 
  mutate(No_sim = case_when(
    Experiment_type == "Repro" ~ No_sim+11,
    .default = No_sim
  ))


{
df_sim <- f_MCSim_read_sim(file.path(path_E_IMD_A, "sim.out"))

pW <- ggplot()+
  geom_line(
    data = df_sim,
    aes(
      x = Time, 
      y = Weight
    ),
    color = Nord_frost[2]
  )+
  geom_point(
    data = df_data_fit_IMD,
    aes(
      x = Time, 
      y = Weight
    )
  )+
  facet_wrap(~No_sim)+
  theme_minimal()

pR <- ggplot()+
  geom_line(
    data = df_sim,
    aes(
      x = Time, 
      y = Reproduction
    ),
    color = Nord_frost[2]
  )+
  geom_point(
    data = df_data_fit_IMD,
    aes(
      x = Time, 
      y = Reproduction
    )
  )+
  facet_wrap(~No_sim)+
  theme_minimal()


pS <- ggplot()+
  geom_line(
    data = df_sim,
    aes(
      x = Time,
      y = Stress
    ),
    color = Nord_frost[2]
  )+
  facet_wrap(~No_sim)+
  theme_minimal()
p <- pW + pR + pS
p
  }




Molec <- "EPX"
path_fig <- here::here("fig/DEBTKTD_Models")

path_E_EPX_A <- here::here("mod/Z_DEB-TKTD/EPX/Test_EPX")

df_data_fit_EPX <- f_import_data_DEBTKTD() |> 
  filter(Molecule == Molec) |> 
  mutate(No_sim = case_when(
    Experiment_type == "Repro" ~ No_sim+9,
    .default = No_sim
  ))

{
  df_sim <- f_MCSim_read_sim(file.path(path_E_EPX_A, "sim.out"))
  
  pW <- ggplot()+
    geom_line(
      data = df_sim,
      aes(
        x = Time, 
        y = Weight
      ),
      color = Nord_frost[2]
    )+
    geom_point(
      data = df_data_fit_EPX,
      aes(
        x = Time, 
        y = Weight
      )
    )+
    facet_wrap(~No_sim)+
    theme_minimal()
  
  pR <- ggplot()+
    geom_line(
      data = df_sim,
      aes(
        x = Time, 
        y = Reproduction
      ),
      color = Nord_frost[2]
    )+
    geom_point(
      data = df_data_fit_EPX,
      aes(
        x = Time, 
        y = Reproduction
      )
    )+
    facet_wrap(~No_sim)+
    theme_minimal()

  
  pS <- ggplot()+
    geom_line(
      data = df_sim,
      aes(
        x = Time, 
        y = Stress
      ),
      color = Nord_frost[2]
    )+
    facet_wrap(~No_sim)+
    theme_minimal()
p <- pW + pR + pS
p
  }



# EPX -----

Molec <- "EPX"
path_fig <- here::here("fig/DEBTKTD_Models")

path_A_EPX_A <- here::here("mod/Z_DEB-TKTD/Test_EPX")

df_data_fit_EPX <- f_import_data_DEBTKTD() |> 
  filter(Molecule == Molec) |> 
  mutate(No_sim = case_when(
    Experiment_type == "Repro" ~ No_sim+9,
    .default = No_sim
  ))

# ./mcsim.DEB-EPX_A DEB_TKTD_sim_ctrl.in sim_ctrl_A.out
# ./mcsim.DEB-EPX_A DEB_TKTD_sim10_0.01.in sim10_0.01_A.out
# ./mcsim.DEB-EPX_A DEB_TKTD_sim100_0.01.in sim100_0.01_A.out
# ./mcsim.DEB-EPX_A DEB_TKTD_sim100_0.1.in sim100_0.1_A.out


# ./mcsim.DEB-EPX_B DEB_TKTD_sim_ctrl.in sim_ctrl_B.out
# ./mcsim.DEB-EPX_B DEB_TKTD_sim10_0.01.in sim10_0.01_B.out
# ./mcsim.DEB-EPX_B DEB_TKTD_sim100_0.01.in sim100_0.01_B.out
# ./mcsim.DEB-EPX_B DEB_TKTD_sim100_0.1.in sim100_0.1_B.out


# ./mcsim.DEB-EPX_C DEB_TKTD_sim_ctrl.in sim_ctrl_C.out
# ./mcsim.DEB-EPX_C DEB_TKTD_sim10_0.01.in sim10_0.01_C.out
# ./mcsim.DEB-EPX_C DEB_TKTD_sim100_0.01.in sim100_0.01_C.out
# ./mcsim.DEB-EPX_C DEB_TKTD_sim100_0.1.in sim100_0.1_C.out

# ./mcsim.DEB-EPX_D DEB_TKTD_sim_ctrl.in sim_ctrl_D.out
# ./mcsim.DEB-EPX_D DEB_TKTD_sim10_0.01.in sim10_0.01_D.out
# ./mcsim.DEB-EPX_D DEB_TKTD_sim100_0.01.in sim100_0.01_D.out
# ./mcsim.DEB-EPX_D DEB_TKTD_sim100_0.1.in sim100_0.1_D.out

{
  df_sim_ctrl_A <- f_MCSim_read_sim(file.path(path_A_EPX_A, "sim_ctrl_A.out"))
  df_sim10_0.01_A <- f_MCSim_read_sim(file.path(path_A_EPX_A, "sim10_0.01_A.out"))
  df_sim100_0.01_A <- f_MCSim_read_sim(file.path(path_A_EPX_A, "sim100_0.01_A.out"))
  df_sim100_0.1_A <- f_MCSim_read_sim(file.path(path_A_EPX_A, "sim100_0.1_A.out"))
  
  pW <- ggplot()+
    geom_point(
      data = df_data_fit_EPX |> filter(ID_experiment == "EPX_EC50_growth_C8"),
      aes(
        x = Time,
        y = Weight
      )
    )+
    geom_line(
      data = df_sim_ctrl_A,
      aes(
        x = Time,
        y = Weight
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_A,
      aes(
        x = Time,
        y = Weight
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_A,
      aes(
        x = Time,
        y = Weight
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_A,
      aes(
        x = Time,
        y = Weight
      ),
      color = "red"
    )+
    theme_minimal()
  
  pR <- ggplot()+
    geom_line(
      data = df_sim_ctrl_A,
      aes(
        x = Time,
        y = Reproduction
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_A,
      aes(
        x = Time,
        y = Reproduction
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_A,
      aes(
        x = Time,
        y = Reproduction
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_A,
      aes(
        x = Time,
        y = Reproduction
      ),
      color = "red"
    )+
    theme_minimal()
  
  pCe <- ggplot()+
    geom_line(
      data = df_sim_ctrl_A,
      aes(
        x = Time,
        y = Ce
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_A,
      aes(
        x = Time,
        y = Ce
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_A,
      aes(
        x = Time,
        y = Ce
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_A,
      aes(
        x = Time,
        y = Ce
      ),
      color = "red"
    )+
    theme_minimal()
  
  pCi <- ggplot()+
    geom_line(
      data = df_sim_ctrl_A,
      aes(
        x = Time,
        y = Ci
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_A,
      aes(
        x = Time,
        y = Ci
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_A,
      aes(
        x = Time,
        y = Ci
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_A,
      aes(
        x = Time,
        y = Ci
      ),
      color = "red"
    )+
    theme_minimal()
  
  pD <- ggplot()+
    geom_line(
      data = df_sim_ctrl_A,
      aes(
        x = Time,
        y = Damage
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_A,
      aes(
        x = Time,
        y = Damage
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_A,
      aes(
        x = Time,
        y = Damage
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_A,
      aes(
        x = Time,
        y = Damage
      ),
      color = "red"
    )+
    theme_minimal()
  
  pS <- ggplot()+
    geom_line(
      data = df_sim_ctrl_A,
      aes(
        x = Time,
        y = Stress
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_A,
      aes(
        x = Time,
        y = Stress
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_A,
      aes(
        x = Time,
        y = Stress
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_A,
      aes(
        x = Time,
        y = Stress
      ),
      color = "red"
    )+
    theme_minimal()
  
  pfreal <- ggplot()+
    geom_line(
      data = df_sim_ctrl_A,
      aes(
        x = Time,
        y = freal
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_A,
      aes(
        x = Time,
        y = freal
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_A,
      aes(
        x = Time,
        y = freal
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_A,
      aes(
        x = Time,
        y = freal
      ),
      color = "red"
    )+
    theme_minimal()
  
  
  ppM <- ggplot()+
    geom_line(
      data = df_sim_ctrl_A,
      aes(
        x = Time,
        y = pM
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_A,
      aes(
        x = Time,
        y = pM
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_A,
      aes(
        x = Time,
        y = pM
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_A,
      aes(
        x = Time,
        y = pM
      ),
      color = "red"
    )+
    theme_minimal()
  
  pg <- ggplot()+
    geom_line(
      data = df_sim_ctrl_A,
      aes(
        x = Time,
        y = g
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_A,
      aes(
        x = Time,
        y = g
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_A,
      aes(
        x = Time,
        y = g
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_A,
      aes(
        x = Time,
        y = g
      ),
      color = "red"
    )+
    theme_minimal()
  
  pkapR_real <- ggplot()+
    geom_line(
      data = df_sim_ctrl_A,
      aes(
        x = Time,
        y = kapR_real
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_A,
      aes(
        x = Time,
        y = kapR_real
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_A,
      aes(
        x = Time,
        y = kapR_real
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_A,
      aes(
        x = Time,
        y = kapR_real
      ),
      color = "red"
    )+
    theme_minimal()
  
  
  p <- pW + pR + pfreal + ppM + pg + pkapR_real + pCe + pCi + pD + pS 
  p
  }

# B ----

{
  df_sim_ctrl_B <- f_MCSim_read_sim(file.path(path_A_EPX_A, "sim_ctrl_B.out"))
  df_sim10_0.01_B <- f_MCSim_read_sim(file.path(path_A_EPX_A, "sim10_0.01_B.out"))
  df_sim100_0.01_B <- f_MCSim_read_sim(file.path(path_A_EPX_A, "sim100_0.01_B.out"))
  df_sim100_0.1_B <- f_MCSim_read_sim(file.path(path_A_EPX_A, "sim100_0.1_B.out"))
  
  pW <- ggplot()+
    geom_point(
      data = df_data_fit_EPX |> filter(ID_experiment == "EPX_EC50_growth_C8"),
      aes(
        x = Time,
        y = Weight
      )
    )+
    geom_line(
      data = df_sim_ctrl_B,
      aes(
        x = Time,
        y = Weight
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_B,
      aes(
        x = Time,
        y = Weight
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_B,
      aes(
        x = Time,
        y = Weight
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_B,
      aes(
        x = Time,
        y = Weight
      ),
      color = "red"
    )+
    theme_minimal()
  
  pR <- ggplot()+
    geom_line(
      data = df_sim_ctrl_B,
      aes(
        x = Time,
        y = Reproduction
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_B,
      aes(
        x = Time,
        y = Reproduction
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_B,
      aes(
        x = Time,
        y = Reproduction
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_B,
      aes(
        x = Time,
        y = Reproduction
      ),
      color = "red"
    )+
    theme_minimal()
  
  pCe <- ggplot()+
    geom_line(
      data = df_sim_ctrl_B,
      aes(
        x = Time,
        y = Ce
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_B,
      aes(
        x = Time,
        y = Ce
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_B,
      aes(
        x = Time,
        y = Ce
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_B,
      aes(
        x = Time,
        y = Ce
      ),
      color = "red"
    )+
    theme_minimal()
  
  pCi <- ggplot()+
    geom_line(
      data = df_sim_ctrl_B,
      aes(
        x = Time,
        y = Ci
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_B,
      aes(
        x = Time,
        y = Ci
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_B,
      aes(
        x = Time,
        y = Ci
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_B,
      aes(
        x = Time,
        y = Ci
      ),
      color = "red"
    )+
    theme_minimal()
  
  pD <- ggplot()+
    geom_line(
      data = df_sim_ctrl_B,
      aes(
        x = Time,
        y = Damage
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_B,
      aes(
        x = Time,
        y = Damage
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_B,
      aes(
        x = Time,
        y = Damage
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_B,
      aes(
        x = Time,
        y = Damage
      ),
      color = "red"
    )+
    theme_minimal()
  
  pS <- ggplot()+
    geom_line(
      data = df_sim_ctrl_B,
      aes(
        x = Time,
        y = Stress
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_B,
      aes(
        x = Time,
        y = Stress
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_B,
      aes(
        x = Time,
        y = Stress
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_B,
      aes(
        x = Time,
        y = Stress
      ),
      color = "red"
    )+
    theme_minimal()
  
  pfreal <- ggplot()+
    geom_line(
      data = df_sim_ctrl_B,
      aes(
        x = Time,
        y = freal
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_B,
      aes(
        x = Time,
        y = freal
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_B,
      aes(
        x = Time,
        y = freal
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_B,
      aes(
        x = Time,
        y = freal
      ),
      color = "red"
    )+
    theme_minimal()
  
  
  ppM <- ggplot()+
    geom_line(
      data = df_sim_ctrl_B,
      aes(
        x = Time,
        y = pM
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_B,
      aes(
        x = Time,
        y = pM
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_B,
      aes(
        x = Time,
        y = pM
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_B,
      aes(
        x = Time,
        y = pM
      ),
      color = "red"
    )+
    theme_minimal()
  
  pg <- ggplot()+
    geom_line(
      data = df_sim_ctrl_B,
      aes(
        x = Time,
        y = g
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_B,
      aes(
        x = Time,
        y = g
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_B,
      aes(
        x = Time,
        y = g
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_B,
      aes(
        x = Time,
        y = g
      ),
      color = "red"
    )+
    theme_minimal()
  
  pkapR_real <- ggplot()+
    geom_line(
      data = df_sim_ctrl_B,
      aes(
        x = Time,
        y = kapR_real
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_B,
      aes(
        x = Time,
        y = kapR_real
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_B,
      aes(
        x = Time,
        y = kapR_real
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_B,
      aes(
        x = Time,
        y = kapR_real
      ),
      color = "red"
    )+
    theme_minimal()
  
  
  p <- pW + pR + pfreal + ppM + pg + pkapR_real + pCe + pCi + pD + pS 
  p
}






# C ----

{
  df_sim_ctrl_C <- f_MCSim_read_sim(file.path(path_A_EPX_A, "sim_ctrl_C.out"))
  df_sim10_0.01_C <- f_MCSim_read_sim(file.path(path_A_EPX_A, "sim10_0.01_C.out"))
  df_sim100_0.01_C <- f_MCSim_read_sim(file.path(path_A_EPX_A, "sim100_0.01_C.out"))
  df_sim100_0.1_C <- f_MCSim_read_sim(file.path(path_A_EPX_A, "sim100_0.1_C.out"))
  
  pW <- ggplot()+
    geom_point(
      data = df_data_fit_EPX |> filter(ID_experiment == "EPX_EC50_growth_C8"),
      aes(
        x = Time,
        y = Weight
      )
    )+
    geom_line(
      data = df_sim_ctrl_C,
      aes(
        x = Time,
        y = Weight
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_C,
      aes(
        x = Time,
        y = Weight
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_C,
      aes(
        x = Time,
        y = Weight
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_C,
      aes(
        x = Time,
        y = Weight
      ),
      color = "red"
    )+
    theme_minimal()
  
  pR <- ggplot()+
    geom_line(
      data = df_sim_ctrl_C,
      aes(
        x = Time,
        y = Reproduction
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_C,
      aes(
        x = Time,
        y = Reproduction
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_C,
      aes(
        x = Time,
        y = Reproduction
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_C,
      aes(
        x = Time,
        y = Reproduction
      ),
      color = "red"
    )+
    theme_minimal()
  
  pCe <- ggplot()+
    geom_line(
      data = df_sim_ctrl_C,
      aes(
        x = Time,
        y = Ce
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_C,
      aes(
        x = Time,
        y = Ce
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_C,
      aes(
        x = Time,
        y = Ce
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_C,
      aes(
        x = Time,
        y = Ce
      ),
      color = "red"
    )+
    theme_minimal()
  
  pCi <- ggplot()+
    geom_line(
      data = df_sim_ctrl_C,
      aes(
        x = Time,
        y = Ci
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_C,
      aes(
        x = Time,
        y = Ci
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_C,
      aes(
        x = Time,
        y = Ci
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_C,
      aes(
        x = Time,
        y = Ci
      ),
      color = "red"
    )+
    theme_minimal()
  
  pD <- ggplot()+
    geom_line(
      data = df_sim_ctrl_C,
      aes(
        x = Time,
        y = Damage
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_C,
      aes(
        x = Time,
        y = Damage
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_C,
      aes(
        x = Time,
        y = Damage
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_C,
      aes(
        x = Time,
        y = Damage
      ),
      color = "red"
    )+
    theme_minimal()
  
  pS <- ggplot()+
    geom_line(
      data = df_sim_ctrl_C,
      aes(
        x = Time,
        y = Stress
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_C,
      aes(
        x = Time,
        y = Stress
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_C,
      aes(
        x = Time,
        y = Stress
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_C,
      aes(
        x = Time,
        y = Stress
      ),
      color = "red"
    )+
    theme_minimal()
  
  pfreal <- ggplot()+
    geom_line(
      data = df_sim_ctrl_C,
      aes(
        x = Time,
        y = freal
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_C,
      aes(
        x = Time,
        y = freal
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_C,
      aes(
        x = Time,
        y = freal
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_C,
      aes(
        x = Time,
        y = freal
      ),
      color = "red"
    )+
    theme_minimal()
  
  
  ppM <- ggplot()+
    geom_line(
      data = df_sim_ctrl_C,
      aes(
        x = Time,
        y = pM
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_C,
      aes(
        x = Time,
        y = pM
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_C,
      aes(
        x = Time,
        y = pM
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_C,
      aes(
        x = Time,
        y = pM
      ),
      color = "red"
    )+
    theme_minimal()
  
  pg <- ggplot()+
    geom_line(
      data = df_sim_ctrl_C,
      aes(
        x = Time,
        y = g
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_C,
      aes(
        x = Time,
        y = g
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_C,
      aes(
        x = Time,
        y = g
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_C,
      aes(
        x = Time,
        y = g
      ),
      color = "red"
    )+
    theme_minimal()
  
  pkapR_real <- ggplot()+
    geom_line(
      data = df_sim_ctrl_C,
      aes(
        x = Time,
        y = kapR_real
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_C,
      aes(
        x = Time,
        y = kapR_real
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_C,
      aes(
        x = Time,
        y = kapR_real
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_C,
      aes(
        x = Time,
        y = kapR_real
      ),
      color = "red"
    )+
    theme_minimal()
  
  
  p <- pW + pR + pfreal + ppM + pg + pkapR_real + pCe + pCi + pD + pS 
  p
}




# D ----



{
  df_sim_ctrl_D <- f_MCSim_read_sim(file.path(path_A_EPX_A, "sim_ctrl_D.out"))
  df_sim10_0.01_D <- f_MCSim_read_sim(file.path(path_A_EPX_A, "sim10_0.01_D.out"))
  df_sim100_0.01_D <- f_MCSim_read_sim(file.path(path_A_EPX_A, "sim100_0.01_D.out"))
  df_sim100_0.1_D <- f_MCSim_read_sim(file.path(path_A_EPX_A, "sim100_0.1_D.out"))
  
  pW <- ggplot()+
    geom_point(
      data = df_data_fit_EPX |> filter(ID_experiment == "EPX_EC50_growth_C8"),
      aes(
        x = Time,
        y = Weight
      )
    )+
    geom_line(
      data = df_sim_ctrl_D,
      aes(
        x = Time,
        y = Weight
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_D,
      aes(
        x = Time,
        y = Weight
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_D,
      aes(
        x = Time,
        y = Weight
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_D,
      aes(
        x = Time,
        y = Weight
      ),
      color = "red"
    )+
    theme_minimal()
  
  pR <- ggplot()+
    geom_line(
      data = df_sim_ctrl_D,
      aes(
        x = Time,
        y = Reproduction
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_D,
      aes(
        x = Time,
        y = Reproduction
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_D,
      aes(
        x = Time,
        y = Reproduction
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_D,
      aes(
        x = Time,
        y = Reproduction
      ),
      color = "red"
    )+
    theme_minimal()
  
  pCe <- ggplot()+
    geom_line(
      data = df_sim_ctrl_D,
      aes(
        x = Time,
        y = Ce
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_D,
      aes(
        x = Time,
        y = Ce
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_D,
      aes(
        x = Time,
        y = Ce
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_D,
      aes(
        x = Time,
        y = Ce
      ),
      color = "red"
    )+
    theme_minimal()
  
  pCi <- ggplot()+
    geom_line(
      data = df_sim_ctrl_D,
      aes(
        x = Time,
        y = Ci
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_D,
      aes(
        x = Time,
        y = Ci
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_D,
      aes(
        x = Time,
        y = Ci
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_D,
      aes(
        x = Time,
        y = Ci
      ),
      color = "red"
    )+
    theme_minimal()
  
  pD <- ggplot()+
    geom_line(
      data = df_sim_ctrl_D,
      aes(
        x = Time,
        y = Damage
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_D,
      aes(
        x = Time,
        y = Damage
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_D,
      aes(
        x = Time,
        y = Damage
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_D,
      aes(
        x = Time,
        y = Damage
      ),
      color = "red"
    )+
    theme_minimal()
  
  pS <- ggplot()+
    geom_line(
      data = df_sim_ctrl_D,
      aes(
        x = Time,
        y = Stress
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_D,
      aes(
        x = Time,
        y = Stress
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_D,
      aes(
        x = Time,
        y = Stress
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_D,
      aes(
        x = Time,
        y = Stress
      ),
      color = "red"
    )+
    theme_minimal()
  
  pfreal <- ggplot()+
    geom_line(
      data = df_sim_ctrl_D,
      aes(
        x = Time,
        y = freal
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_D,
      aes(
        x = Time,
        y = freal
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_D,
      aes(
        x = Time,
        y = freal
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_D,
      aes(
        x = Time,
        y = freal
      ),
      color = "red"
    )+
    theme_minimal()
  
  
  ppM <- ggplot()+
    geom_line(
      data = df_sim_ctrl_D,
      aes(
        x = Time,
        y = pM
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_D,
      aes(
        x = Time,
        y = pM
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_D,
      aes(
        x = Time,
        y = pM
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_D,
      aes(
        x = Time,
        y = pM
      ),
      color = "red"
    )+
    theme_minimal()
  
  pg <- ggplot()+
    geom_line(
      data = df_sim_ctrl_D,
      aes(
        x = Time,
        y = g
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_D,
      aes(
        x = Time,
        y = g
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_D,
      aes(
        x = Time,
        y = g
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_D,
      aes(
        x = Time,
        y = g
      ),
      color = "red"
    )+
    theme_minimal()
  
  pkapR_real <- ggplot()+
    geom_line(
      data = df_sim_ctrl_D,
      aes(
        x = Time,
        y = kapR_real
      ),
      color = "green"
    )+
    geom_line(
      data = df_sim10_0.01_D,
      aes(
        x = Time,
        y = kapR_real
      ),
      color = "yellow"
    )+
    geom_line(
      data = df_sim100_0.01_D,
      aes(
        x = Time,
        y = kapR_real
      ),
      color = "orange"
    )+
    geom_line(
      data = df_sim100_0.1_D,
      aes(
        x = Time,
        y = kapR_real
      ),
      color = "red"
    )+
    theme_minimal()
  
  
  p <- pW + pR + pfreal + ppM + pg + pkapR_real + pCe + pCi + pD + pS 
  p
}






# ABCD ----


pW_ABCD <- ggplot()+
  geom_point(
    data = df_data_fit_EPX |> filter(ID_experiment == "EPX_EC50_growth_C8"),
    aes(
      x = Time,
      y = Weight
    )
  )+
  geom_line(
    data = df_sim100_0.01_A,
    aes(
      x = Time,
      y = Weight
    ),
    color = "red"
  )+
  geom_line(
    data = df_sim100_0.01_B,
    aes(
      x = Time,
      y = Weight
    ),
    color = "orange"
  )+
  geom_line(
    data = df_sim100_0.01_C,
    aes(
      x = Time,
      y = Weight
    ),
    color = "yellow"
  )+
  geom_line(
    data = df_sim100_0.01_D,
    aes(
      x = Time,
      y = Weight
    ),
    color = "green"
  )+
  theme_minimal()


pW_ABCD

















{
  df_sim <- f_MCSim_read_sim(file.path(path_A_EPX_A, "sim.out"))
  
  pW <- ggplot()+
    geom_line(
      data = df_sim,
      aes(
        x = Time,
        y = Weight
      )
    )+
    geom_point(
      data = df_data_fit_EPX,
      aes(
        x = Time,
        y = Weight
      )
    )+
    facet_wrap(~No_sim)+
    theme_minimal()
  
  pW
  
}
  
  pR <- ggplot()+
    geom_line(
      data = df_sim,
      aes(
        x = Time,
        y = Reproduction
      )
    )+
    geom_point(
      data = df_data_fit_EPX,
      aes(
        x = Time,
        y = Reproduction
      )
    )+
    facet_wrap(~No_sim)+
  theme_minimal()
  
  pCe <- ggplot()+
    geom_line(
      data = df_sim,
      aes(
        x = Time,
        y = Ce
      )
    )+
    facet_wrap(~No_sim)+
  theme_minimal()
  
  pCi <- ggplot()+
    geom_line(
      data = df_sim,
      aes(
        x = Time,
        y = Ci
      )
    )+
    facet_wrap(~No_sim)+
  theme_minimal()
  
  pD <- ggplot()+
    geom_line(
      data = df_sim,
      aes(
        x = Time,
        y = Damage
      )
    )+
    facet_wrap(~No_sim)+
  theme_minimal()
  
  pS <- ggplot()+
    geom_line(
      data = df_sim,
      aes(
        x = Time,
        y = Stress
      )
    )+
    facet_wrap(~No_sim)+
  theme_minimal()
  
  pfreal <- ggplot()+
    geom_line(
      data = df_sim,
      aes(
        x = Time,
        y = freal
      )
    )+
    facet_wrap(~No_sim)+
    theme_minimal()
  
  
  pW
  pR
  pfreal
  pCe
  pCi
  pD
  pS
}

