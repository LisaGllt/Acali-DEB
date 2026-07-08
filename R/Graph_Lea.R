df_data <- data.frame(
  Molecule = c("A", "A", "B", "B", "B", "C", "C"), # A éventuellement réordonner
  Dose = c(10, 100, 10, 100, 1000, 10, 100),       # idem
  Reponse = c(20, 10, 20, 10, 5, 20, 10)
) |> 
  mutate(
    Condition = paste(Molecule, Dose, sep=":")
  )

# Version facet_wrap ----

# pcq je crois que c'est ce qu'on avait fait ? Je ne sais plus

# Plot problématique

ggplot(
  data = df_data, 
  aes(
    x = as.factor(Dose),
    y = Reponse
  )
)+
  geom_col()+
  facet_wrap(
    ~Molecule
    )+
  theme_minimal()+
  theme(
    strip.text = element_blank()
  )

# Mais du coup est ce que la version du haut n'est pas plus simple que ce qu'on avait fait ? 


# On enlève la répétition de 1000 pour ceux où il n'y a rien
ggplot(
  data = df_data, 
  aes(
    x = as.factor(Dose),
    y = Reponse
  )
)+
  geom_col()+
  facet_wrap(
    ~Molecule,
    scales = "free_x" # Avec ça
  )+
  theme_minimal()+
  theme(
    strip.text = element_blank()
  )

# Problème : la taille des facets sont les mêmes 
# donc les largeurs des colonnes ne le sont pas

# TADAAAAA
ggplot(
  data = df_data, 
  aes(
    x = as.factor(Dose),
    y = Reponse
  )
)+
  geom_col()+
  facet_wrap(
    ~Molecule,
    scales = "free_x", 
    space = "free_x" # Avec ça, je connaissais pas du tout mdrr merci !
  )+
  theme_minimal()+
  theme(
    strip.text = element_blank()
  )


# Version sans facet wrap si jamais  ----

# Plot problématique initial simple
ggplot(
  data = df_data, 
  aes(
    x = Condition,
    y = Reponse
  )
) +
  geom_col()

# Plot corrigé avec réécriture de l'axe des x
ggplot(
  data = df_data, 
  aes(
    x = Condition,
    y = Reponse
  )
) +
  geom_col()+
  scale_x_discrete( # C'est bien discrete et pas continuous sinon il crie
    breaks = c("A:10", "A:100", "B:10", "B:100", "B:1000", "C:10", "C:100"),
    labels = c("10", "100", "10", "100", "1000", "10", "100")
  )
