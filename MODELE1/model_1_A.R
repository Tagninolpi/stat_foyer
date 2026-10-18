library(coda)
library(rjags)

setwd("C:/Users/manch//Desktop/M2 Bioinfo/Stats/PROJET")


## Jeu de données


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

# Bonne convergence:


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
  main = "M0 - Posterior de a",
  xlab = "a",
  col = "lightblue"
)

hist(
  postA[, "b"],
  breaks = 50,
  main = "M0 - Posterior de b",
  xlab = "b",
  col = "lightblue"
)


## Probabilites d'infection a posteriori

# Personnes non vaccinées
pA_non_vaccinated <- plogis(postA[, "a"]) # transforme le logit en proba
mean(pA_non_vaccinated)
quantile(pA_non_vaccinated, c(0.025, 0.5, 0.975))

# Personnes protegees
pA_vaccinated <- plogis(postA[, "a"] + postA[, "b"])
mean(pA_vaccinated)
quantile(pA_vaccinated, c(0.025, 0.5, 0.975))


## Distribution a posteriori des probabilites

hist(
  pA_non_vaccinated,
  breaks = 50,
  main = "M0 - Probabilite d'infection : non vacciné",
  xlab = "Probabilite d'infection",
  col = "lightblue"
)

hist(
  pA_vaccinated,
  breaks = 50,
  main = "M0 - Probabilite d'infection : vacciné",
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

N <- nrow(df_model)

for (j in 1:1000) {
  
  a <- postA[j, "a"]
  b <- postA[j, "b"]
  
  #Probabilite d'infection de chaque contact
  p <- plogis(a + b * df_model$vaccination_status)
  
  #Simuler les infections
  infected_sim <- rbinom(N, size = 1, prob = p)
  
  #Nombre total d'infectes
  simulatedA[j] <- sum(infected_sim)
  
  # taux d'infection chez les no vaccines
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
     main = "PPC - Nombre d'infectés",
     xlab = "Nombre d'infectés simulés")
abline(v = observed, col = "red", lwd = 2)

# Taux d'infection chez les non-vaccinés
quantile(simulatedA_non_vaccinated, c(0.025, 0.5, 0.975))
observed_non_vaccinated

hist(simulatedA_non_vaccinated, breaks = 30,
     main = "PPC - Non-vaccinés",
     xlab = "Taux d'infection simulé")
abline(v = observed_non_vaccinated, col = "red", lwd = 2)

# Taux d'infection chez les vaccinés
quantile(simulatedA_vaccinated, c(0.025, 0.5, 0.975))
observed_vaccinated

hist(simulatedA_vaccinated, breaks = 30,
     main = "PPC - Vaccinés",
     xlab = "Taux d'infection simulé")
abline(v = observed_vaccinated, col = "red", lwd = 2)






