library(here)
source(file=here::here("functions/fun.R"))

f_load_libraries_colors()

w_deb <- 27.24
Shape_deb <- 0.066193;
f <- 1


# L_Ww_deb = (Ww/(1+w_deb*f))^(1/3),
# L_Lw_deb = Lw*Shape_deb,
# 
# Lw_Shape = (Ww/(1+w_deb))^(1/3)/Shape_deb#,
#Lw_w = Weight^(1/3)/alpha

df_rlw <- read_excel(here::here("data/Data_expe_raw/Data_RLW.xlsx")) |> 
  filter(Person == "Sylvain" | (Person=="Laura" & t==0)) |> 
  mutate(
    Ww=w/1000, # mg to g
    Cub_Lw = Lw^3
    ) 

# Ww = L^3 x [Shape^3(1+w)] = L^3 x alpha --> On ne peut estimer que alpha
lm.rlw <- lm(Ww~Cub_Lw+0, data=df_rlw)
summary(lm.rlw)
sum.lm.rlw <- summary(lm.rlw)
alpha <- lm.rlw$coefficients[1]
R2 <- sum.lm.rlw$r.squared

p <- ggplot(
  data = df_rlw
)+
  geom_point(
    alpha = 0.5,
    aes(
      y=Ww,
      x=Cub_Lw,
      color = Status
    )
  )+
  geom_abline(
    slope = alpha, 
    intercept = 0,
    color = Nord_polar[4]
      )+
  labs(
    x="Cubic measured length, Lw^3 (cm^3)",
    y="Measured wet weight, Ww"
  )+
  annotate(
    "text",
    label = paste0("Ww = ", round(alpha, 5), " x Lw3 with R2 = ", round(R2,2)),
    x = 300, 
    y = 0.12
  )+
  xlim(0,NA)+
  ylim(0,NA)+
  scale_color_manual(values = c(Nord_aurora[1], Nord_aurora[4]))+
  theme_minimal()
p

# Donc Shape^3 x (1+w) = 0.00205

Alpha_deb <- Shape_deb^3 * (1+w_deb)
Alpha_deb # 0.0082
