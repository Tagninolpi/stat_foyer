library(coda)
library(rjags)

setwd("C:/Users/manch//Desktop/M2 Bioinfo/Stats/PROJET")


## Jeu de données


df <- read.table("2021_05_14_full_database_2doses.txt")
df

colnames(df) <- c(
  "individual_id",
  "house_id",
  "household_size",
  "date_symptom",
  "infection_status",
  "vaccination_status",
  "end_follow_up",
  "age",
  "index_cases",
  "isolation_status"
)

dim(df)
names(df)
head(df)
tail(df)

df$infected <- ifelse(df$infection_status == 0, 0, 1)
df$vaccinated <- ifelse(df$vaccination_status == 0, 0, 1)
df
dim(df)
names(df)

N <- nrow(df)
N

################################
##### Modele 0, vaccinated ou pas:
###################################

# 95 % des valeurs de b sont entre:
# OR = ex[]
# qnorm(c(0.025, 0.975), mean = 0, sd = 2.5) une seule personne supplémentaire dans le foyer multiplie les odds d'infection par le exp(res)

model_M0 <- "model {
  a ~ dnorm(0, 0.16)
  b ~ dnorm(0, 0.16)

  for (i in 1:N) {
    logit(p[i]) <- a + b * vaccinated[i]
    infected[i] ~ dbern(p[i])
  }
}"


data_M0 <- list(N = N, infected = df$infected, vaccinated = df$vaccinated)

# Creation du sampler
model_jags0 <- jags.model(
  textConnection(model_M0),
  data = data_M0,
  n.chains = 3
)

# le burnin pour que les valeurs convergent: les premieres 
# iterations de la chaine MCMC ne sont pas conservées
update(model_jags0, 2000)

# mcmc sample
samples0 <- coda.samples(
  model_jags0,
  variable.names = c("a", "b"),
  n.iter = 10000
)


# Bonne convergence?


summary(samples0)

gelman.diag(samples0) 
# Point est = 1 estun bon signe que les chaînes ont convergé vers la même distribution

effectiveSize(samples0)

plot(samples0)
# le graphique est cohérent avec une bonne convergence

autocorr.plot(samples0)



# Distribution a posteriori des parametres

post0 <- as.matrix(samples0)

hist(
  post0[, "a"],
  breaks = 50,
  main = "M0 - Posterior de a",
  xlab = "a",
  col = "lightblue"
)

hist(
  post0[, "b"],
  breaks = 50,
  main = "M0 - Posterior de b",
  xlab = "b",
  col = "lightblue"
)


## Probabilites d'infection a posteriori

# Personnes non protegees
p0_non_vaccinated <- plogis(post0[, "a"]) # transforme le logit en proba
mean(p0_non_vaccinated)
quantile(p0_non_vaccinated, c(0.025, 0.5, 0.975))

# Personnes protegees
p0_vaccinated <- plogis(post0[, "a"] + post0[, "b"])
mean(p0_vaccinated)
quantile(p0_vaccinated, c(0.025, 0.5, 0.975))


## Distribution a posteriori des probabilites

hist(
  p0_non_vaccinated,
  breaks = 50,
  main = "M0 - Probabilite d'infection : non proteges",
  xlab = "Probabilite d'infection",
  col = "lightblue"
)

hist(
  p0_vaccinated,
  breaks = 50,
  main = "M0 - Probabilite d'infection : proteges",
  xlab = "Probabilite d'infection",
  col = "lightblue"
)


## Comparaison obs vs predit

# Probabilites observees

# grouper les ind selon leur status vaccin 
# et calculer la proportion d'infectes dans chaque groupe

observed0 <- aggregate( 
  infected ~ vaccinated,
  data = df,
  FUN = mean
)

print(observed0)
names(observed0)[2] <- "observed_p"

# Probabilites predite

profiles0 <- data.frame(
  vaccinated = c(0, 1),
  mean_p = c(
    mean(p0_non_vaccinated),
    mean(p0_vaccinated)
  ),
  lower_95 = c(
    quantile(p0_non_vaccinated, 0.025),
    quantile(p0_vaccinated, 0.025)
  ),
  upper_95 = c(
    quantile(p0_non_vaccinated, 0.975),
    quantile(p0_vaccinated, 0.975)
  )
)


# Comparaison: on va merge selon leur statut de vaccination

comparison0 <- merge(
  observed0,
  profiles0,
  by = "vaccinated"
)

comparison0


## Graphique observe vs predit

# beside pour afficher les barres d'un même groupe côte à côte
barplot(
  rbind(comparison0$observed_p, comparison0$mean_p),
  beside = TRUE,
  names.arg = c("Non proteges", "Proteges"),
  col = c("blue", "orange"),
  ylim = c(0, 1),
  ylab = "Probabilite d'infection",
  main = "M0 - Observe vs predit"
)

legend(
  "topright",
  legend = c("Observe", "Predit"),
  fill = c("blue", "orange")
)

# On observe que les probabilites predites sont tres proches
# des probabilites obs pour les deux groupes
# Les personnes protegees ont une proba d'infection
# plus faible  22 %) que les non protegees (60 %)
# Le modele reproduit donc bien les proportions obs,
# mais cela ne suffit pas a valider sa qualite predictive vu qu'il utilise
# les meme données pour estimer les paramètres et comparer les res



## Verification du modele par simulation

observed <- sum(df$infected)
simulated0 <- numeric(1000)

for (j in 1:1000) {
  
  a <- post0[j, "a"]
  b <- post0[j, "b"]
  
  # Proba des 902 individus
  p <- plogis(a + b * df$vaccinated)
  
  # Simuler les infections pour chaque p
  infected_sim <- rbinom(N, size = 1, prob = p)
  
  # Nombre total d'infectes
  simulated0[j] <- sum(infected_sim)
}

hist(
  simulated0,
  breaks = 30,
  main = "M0 - PPC",
  xlab = "Nombre d'infectes simules",
  col = "lightblue"
)

abline(v = observed, col = "red", lwd = 2)

# La plupart des simulations donnent environ 460 à 510 infectes!!

# Faut refaire le meme travail sur les autres modeles...

################################
##### Modele 1, taille du foyer:
#####################################



model_M1 <- "model {
  a ~ dnorm(0, 0.16)
  b ~ dnorm(0, 0.16)

  for (i in 1:N) {
    logit(p[i]) <- a + b * size[i]
    infected[i] ~ dbern(p[i])
  }
}"


data_M1 <- list(N = N, infected = df$infected, size = df$household_size)

# Creation du sampler
model_jags1 <- jags.model(
  textConnection(model_M1),
  data = data_M1,
  n.chains = 3
)

# le burnin pour que les valeurs convergent: les premieres 
# iterations de la chaine MCMC ne sont pas conservées
update(model_jags1, 2000)

# mcmc sample
samples1 <- coda.samples(
  model_jags1,
  variable.names = c("a", "b"),
  n.iter = 10000
)


# Bonne convergence?


summary(samples1)

gelman.diag(samples1) 
# Point est = 1 estun bon signe que les chaînes ont convergé vers la même distribution

effectiveSize(samples1)

plot(samples1)
# le graphique est cohérent avec une bonne convergence

autocorr.plot(samples1)



# Distribution a posteriori des parametres

post1 <- as.matrix(samples1)

hist(
  post1[, "a"],
  breaks = 50,
  main = "M1 - Posterior de a",
  xlab = "a",
  col = "lightblue"
)

hist(
  post1[, "b"],
  breaks = 50,
  main = "M1 - Posterior de b",
  xlab = "b",
  col = "lightblue"
)

## Probabilites d'infection a posteriori

# Les differentes tailles de foyer
house_sizes <- sort(unique(df$household_size))
house_sizes

# Tableau pour stocker les resultats
profiles1 <- data.frame(
  household_size = house_sizes,
  mean_p = NA,
  median_p = NA,
  lower_95 = NA,
  upper_95 = NA
)

# Calcul du posterior de p pour chaque taille du foyer
for (j in 1:nrow(profiles1)) {
  
  L_post <- post1[, "a"] +
    post1[, "b"] * profiles1$household_size[j]
  
  p_post <- plogis(L_post)
  
  profiles1$mean_p[j] <- mean(p_post)
  profiles1$median_p[j] <- median(p_post)
  profiles1$lower_95[j] <- quantile(p_post, 0.025)
  profiles1$upper_95[j] <- quantile(p_post, 0.975)
}

profiles1


## Distribution a posteriori des probabilites

# Exemple pour un foyer de 3 personnes
p1_size3 <- plogis(post1[, "a"] + post1[, "b"] * 3)

hist(
  p1_size3,
  breaks = 50,
  main = "M1 - Probabilite d'infection : foyer de 3 personnes",
  xlab = "Probabilite d'infection",
  col = "lightblue"
)



## Comparaison observe vs predit

# Probabilites observees
observed1 <- aggregate(
  infected ~ household_size,
  data = df,
  FUN = mean
)

print(observed1)

names(observed1)[2] <- "observed_p"


# Comparaison
comparison1 <- merge(
  observed1,
  profiles1,
  by = "household_size"
)

comparison1

comparison1[, c("household_size", "observed_p", "mean_p")]


## Graphique observe vs predit

barplot(
  rbind(comparison1$observed_p, comparison1$mean_p),
  beside = TRUE,
  names.arg = comparison1$household_size,
  col = c("steelblue", "orange"),
  ylim = c(0, 1),
  xlab = "Taille du foyer",
  ylab = "Probabilite d'infection",
  main = "M1 - Observe vs predit"
)

legend(
  "topright",
  legend = c("Observe", "Predit"),
  fill = c("steelblue", "orange")
)
