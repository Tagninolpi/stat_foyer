library(coda)
library(rjags)

setwd("C:/Users/manch//Desktop/M2 Bioinfo/Stats/PROJET")


### Jeu de données


df <- read.table("2021_05_14_full_database_2doses.txt")

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

head(df)
dim(df)
names(df)
table(df$index_cases)

# Filtrer les cas contacts
df_model <- df[df$index_cases == 0, ] # on perd 215 ind index

# Garder uniquement les statuts vaccinaux 0 et 1
df_model <- df_model[df_model$vaccination_status %in% c(0, 1), ] #on perd 20 ind infecte précédemment 

# Créer la variable infected
df_model$infected <- ifelse(df_model$infection_status == 0, 0, 1)

dim(df)
dim(df_model)

N <- nrow(df_model)
N



## Param check

plogis(qnorm(c(0.025, 0.975), mean = 0, sd = 1))
exp(qnorm(c(0.025, 0.975), mean = 0, sd = 1))


### Modele 1 A, vaccinated ou pas:


model_M1_A <- "model {
  a ~ dnorm(0, 1)
  b ~ dnorm(0, 1)

  for (i in 1:N) {
    logit(p[i]) <- a + b * vaccinated[i]
    infected[i] ~ dbern(p[i])
  }
}"

data_M1_A <- list(N = N, infected = df_model$infected, vaccinated = df_model$vaccination_status)

model_jagsA <- jags.model(
  textConnection(model_M1_A),
  data = data_M1_A,
  n.chains = 3
)

update(model_jagsA, 2000)

samplesA <- coda.samples(
  model_jagsA,
  variable.names = c("a", "b"),
  n.iter = 10000
)

# Bonne convergence: oui

summary(samplesA)

gelman.diag(samplesA) 
# Point est = 1 est un bon signe que les chaînes ont convergé vers 
# la même distribution

effectiveSize(samplesA)

plot(samplesA)
# le graphique est cohérent avec une bonne convergence

autocorr.plot(samplesA)


### Posterior:

## Distribution a posteriori des parametres

postA <- as.matrix(samplesA)

# OR a posteriori : effet de la protection
quantile(exp(postA[, "b"]), c(0.025, 0.5, 0.975))


hist(
  postA[, "a"],
  breaks = 50,
  main = "Posterior de a (M1 vaccin)",
  xlab = "a",
  col = "lightblue"
)

hist(
  postA[, "b"],
  breaks = 50,
  main = "Posterior de b (M1 vaccin)",
  xlab = "b",
  col = "lightblue"
)



### Proba d'infection a posteriori

# proba a posteriori
pA_non_vaccinated <- plogis(postA[, "a"])
pA_vaccinated <- plogis(postA[, "a"] + postA[, "b"])

# tableau recap
tableau_A <- data.frame(
  Statut_vaccinal = c("Non vacciné", "Vacciné"),
  
  Probabilite_moyenne = round(100 * c(
    mean(pA_non_vaccinated),
    mean(pA_vaccinated)
  ), 1),
  
  IC_95 = c(
    paste0("[",
           round(100 * quantile(pA_non_vaccinated, 0.025), 1),
           " % ; ",
           round(100 * quantile(pA_non_vaccinated, 0.975), 1),
           " %]"),
    
    paste0("[",
           round(100 * quantile(pA_vaccinated, 0.025), 1),
           " % ; ",
           round(100 * quantile(pA_vaccinated, 0.975), 1),
           " %]")
  ),
  
  Effectif = c(
    sum(df_model$vaccination_status == 0),
    sum(df_model$vaccination_status == 1)
  )
)

print(tableau_A, row.names = FALSE)



## Distribution a posteriori des probabilites d'infection

hist(
  pA_non_vaccinated,
  breaks = 50,
  main = "Probabilite d'infection : non vacciné",
  xlab = "Probabilite d'infection",
  col = "lightblue"
)

hist(
  pA_vaccinated,
  breaks = 50,
  main = "Probabilite d'infection : vacciné",
  xlab = "Probabilite d'infection",
  col = "lightblue"
)


### PPC

# Nombre total d'infectes obs
observed <- sum(df_model$infected)

# Taux d'infection obs
observed_non_vaccinated <- mean(
  df_model$infected[df_model$vaccination_status == 0]
)

observed_vaccinated <- mean(
  df_model$infected[df_model$vaccination_status == 1]
)

# Stocker les simulations
simulatedA <- numeric(1000)
simulatedA_non_vaccinated <- numeric(1000)
simulatedA_vaccinated <- numeric(1000)


for (j in 1:1000) {
  
  a <- postA[j, "a"]
  b <- postA[j, "b"]
  
  # Proba d'infection de chaque contact
  p <- plogis(a + b * df_model$vaccination_status)
  
  # Simuler les infections
  infected_sim <- rbinom(N, size = 1, prob = p)
  
  # Nombre total d'infectes
  simulatedA[j] <- sum(infected_sim)
  
  #taux d'infection chez les no vaccines
  simulatedA_non_vaccinated[j] <- mean(
    infected_sim[df_model$vaccination_status == 0]
  )
  
  #taux d'infection chez les vaccines
  simulatedA_vaccinated[j] <- mean(
    infected_sim[df_model$vaccination_status == 1]
  )
}


# Nombre total d'infectés
quantile(simulatedA, c(0.025, 0.5, 0.975))
observed

hist(simulatedA, breaks = 30,
     main = "PPC Nombre d'infectés",
     xlab = "Nombre d'infectés simulés")
abline(v = observed, col = "red", lwd = 2)


# Taux d'infection chez les non-vaccinés
quantile(simulatedA_non_vaccinated, c(0.025, 0.5, 0.975))
observed_non_vaccinated

hist(simulatedA_non_vaccinated, breaks = 30,
     main = "PPC Non-vaccinés",
     xlab = "Taux d'infection simulé")
abline(v = observed_non_vaccinated, col = "red", lwd = 2)


# Taux d'infection chez les vaccinés
quantile(simulatedA_vaccinated, c(0.025, 0.5, 0.975))
observed_vaccinated

hist(simulatedA_vaccinated, breaks = 30,
     main = "PPC Vaccinés",
     xlab = "Taux d'infection simulé")
abline(v = observed_vaccinated, col = "red", lwd = 2)



### Sensibilite aux priors


# Nouveau prior pour b avec un prior plus large ecart type = 0.5


model_M1_A_sensi <- "model {
  a ~ dnorm(0, 1)
  b ~ dnorm(-1, 1)

  for (i in 1:N) {
    logit(p[i]) <- a + b * vaccinated[i]
    infected[i] ~ dbern(p[i])
  }
}"

model_jagsA_sensi <- jags.model(
  textConnection(model_M1_A_sensi),
  data = data_M1_A,
  n.chains = 3
)

update(model_jagsA_sensi, 2000)

samplesA_sensi <- coda.samples(
  model_jagsA_sensi,
  variable.names = c("a", "b"),
  n.iter = 10000
)

postA_sensi <- as.matrix(samplesA_sensi)


gelman.diag(samplesA_sensi)
effectiveSize(samplesA_sensi)


# Comparer les odds ratio a posteriori
OR_initial <- exp(postA[, "b"])
OR_sensi <- exp(postA_sensi[, "b"])


comparaison_OR <- data.frame(
  Prior = c("b ~ dnorm(0,1)", "b ~ dnorm(0,4)"),
  OR_median = c(
    median(OR_initial),
    median(OR_sensi)
  ),
  IC_inf = c(
    quantile(OR_initial, 0.025),
    quantile(OR_sensi, 0.025)
  ),
  IC_sup = c(
    quantile(OR_initial, 0.975),
    quantile(OR_sensi, 0.975)
  )
)

print(comparaison_OR, row.names = FALSE)

# Comparer les proba d'infection

p_vaccinated_sensi <- plogis(
  postA_sensi[, "a"] + postA_sensi[, "b"]
)

mean(pA_vaccinated)
mean(p_vaccinated_sensi)


### DIC

load.module("dic")

# Modele A

dic_M1_A <- dic.samples(
  model_jagsA,
  n.iter = 10000,
  type = "pD"
)

sum(dic_M1_A$deviance) + sum(dic_M1_A$penalty)



dic_M1_A_sensi <- dic.samples(
  model_jagsA_sensi,
  n.iter = 10000,
  type = "pD"
)

sum(dic_M1_A_sensi$deviance) + sum(dic_M1_A_sensi$penalty)


