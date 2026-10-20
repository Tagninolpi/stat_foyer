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

df_model <- df[df$index_cases == 0, ] 
df_model <- df_model[df_model$vaccination_status %in% c(0, 1), ] 
df_model$infected <- ifelse(df_model$infection_status == 0, 0, 1)

N <- nrow(df_model)
N

## Param check

# a
plogis(qnorm(c(0.025, 0.975), mean = 0, sd = 1))

# c initial
exp(qnorm(c(0.025, 0.975), mean = 0, sd = 0.5))

# c sensi
exp(qnorm(c(0.025, 0.975), mean = 0, sd = 2.5))


### Modele 1 B, taille du foyer:


model_M1_B <- "model {
  a ~ dnorm(0, 1)
  c ~ dnorm(0, 4)

  for (i in 1:N) {
    logit(p[i]) <- a + c * (size[i] - 2) 
    infected[i] ~ dbern(p[i])
  }
}"

# On fait size[i] - 2 car l'intercept a représente la proba d'infection quand la taille 
# du foyer est de 2. En plus la taille minimale d'un foyer dans notre jeux de données est de 2

# A check le householdsize
data_M1_B <- list(N = N, infected = df_model$infected, size = df_model$household_size)

model_jagsB <- jags.model(
  textConnection(model_M1_B),
  data = data_M1_B,
  n.chains = 3
)

update(model_jagsB, 2000)

samplesB <- coda.samples(
  model_jagsB,
  variable.names = c("a", "c"),
  n.iter = 10000
)

summary(samplesB)
gelman.diag(samplesB) 
effectiveSize(samplesB)
plot(samplesB)
autocorr.plot(samplesB)



### Posterior:

## Distribution a posteriori des parametres

postB <- as.matrix(samplesB)

# OR a posteriori : effet de la taille du foyer
quantile(exp(postB[, "c"]), c(0.025, 0.5, 0.975))


hist(
  postB[, "a"],
  breaks = 50,
  main = "Posterior de a (M1 taille du foyer)",
  xlab = "a",
  col = "lightblue"
)

hist(
  postB[, "c"],
  breaks = 50,
  main = "Posterior de c (M1 taille du foyer)",
  xlab = "c",
  col = "lightblue"
)



### Proba d'infection a posteriori

# proba a posteriori
tailles <- sort(unique(df_model$household_size))

pB <- sapply(tailles, function(taille) {
  plogis(postB[, "a"] + postB[, "c"] * (taille - 2))
})

# tableau recap
tableau_B <- data.frame(
  Taille_foyer = tailles,
  
  Probabilite_moyenne = round(100 * colMeans(pB), 1),
  
  IC_95 = sapply(1:ncol(pB), function(i) {
    paste0("[",
           round(100 * quantile(pB[, i], 0.025), 1),
           " % ; ",
           round(100 * quantile(pB[, i], 0.975), 1),
           " %]")
  }),
  
  Effectif = as.integer(table(
    factor(df_model$household_size, levels = tailles)
  ))
)

print(tableau_B, row.names = FALSE)



## Distribution a posteriori des probabilites d'infection

for (i in seq_along(tailles)) {
  hist(
    pB[, i],
    breaks = 50,
    main = paste("Probabilite d'infection : foyer de", tailles[i], "personnes"),
    xlab = "Probabilite d'infection",
    col = "lightblue"
  )
}



### PPC

# Nombre total d'infectes obs
observed <- sum(df_model$infected)

# Taux d'infection obs par taille
observed_B <- sapply(tailles, function(taille) {
  mean(df_model$infected[df_model$household_size == taille])
})

# Stocker les simulations
simulatedB <- numeric(1000)
simulatedB_tailles <- matrix(NA, nrow = 1000, ncol = length(tailles))

colnames(simulatedB_tailles) <- tailles

indices <- sample(1:nrow(postB), 1000)

for (j in 1:1000) {
  
  a <- postB[indices[j], "a"]
  c <- postB[indices[j], "c"]
  
  # Proba d'infection de chaque contact
  p <- plogis(a + c * (df_model$household_size - 2))
  
  # Simuler les infections
  infected_sim <- rbinom(N, size = 1, prob = p)
  
  # Nombre total d'infectes
  simulatedB[j] <- sum(infected_sim)
  
  # taux d'infection selon la taille du foyer
  simulatedB_tailles[j, ] <- sapply(tailles, function(taille) {
    mean(infected_sim[df_model$household_size == taille])
  })
}


# Nombre total d'infectés
quantile(simulatedB, c(0.025, 0.5, 0.975))
observed

hist(simulatedB, breaks = 30,
     main = "PPC Nombre d'infectés",
     xlab = "Nombre d'infectés simulés")
abline(v = observed, col = "red", lwd = 2)


# Taux d'infection selon la taille du foyer

for (i in seq_along(tailles)) {
  
  cat("Taille :", tailles[i], "\n")
  print(quantile(simulatedB_tailles[, i], c(0.025, 0.5, 0.975)))
  print(observed_B[i])
  
  hist(simulatedB_tailles[, i], breaks = 30,
       main = paste("PPC Foyer de", tailles[i], "personnes"),
       xlab = "Taux d'infection simulé")
  abline(v = observed_B[i], col = "red", lwd = 2)
}


### Sensibilite aux priors


# Prior alternatif pour c : ecart type = 2.5
# contre 0.5 dans le modèle initial

model_M1_B_sensi <- "model {
  a ~ dnorm(0, 1)
  c ~ dnorm(0, 0.16)

  for (i in 1:N) {
    logit(p[i]) <- a + c * (size[i] - 2)
    infected[i] ~ dbern(p[i])
  }
}"

model_jagsB_sensi <- jags.model(
  textConnection(model_M1_B_sensi),
  data = data_M1_B,
  n.chains = 3
)

update(model_jagsB_sensi, 2000)

samplesB_sensi <- coda.samples(
  model_jagsB_sensi,
  variable.names = c("a", "c"),
  n.iter = 10000
)

postB_sensi <- as.matrix(samplesB_sensi)


gelman.diag(samplesB_sensi)
effectiveSize(samplesB_sensi)


# Comparer les odds ratio a posteriori
OR_initial <- exp(postB[, "c"])
OR_sensi <- exp(postB_sensi[, "c"])


comparaison_OR <- data.frame(
  Prior = c("c ~ dnorm(0,4)", "c ~ dnorm(0,0.16)"),
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

# Comparer les proba d'infection

pB_sensi <- sapply(tailles, function(taille) {
  plogis(postB_sensi[, "a"] + postB_sensi[, "c"] * (taille - 2))
})

comparaison_proba_B <- data.frame(
  Taille_foyer = tailles,
  
  "Proba_initiale (dnorm(0,4))" = round(100 * colMeans(pB), 1),
  
  "Proba_sensi (dnorm(0,0.16))" = round(100 * colMeans(pB_sensi), 1),
  
  check.names = FALSE
)

print(comparaison_proba_B, row.names = FALSE)



### DIC

load.module("dic")

# Modele B

dic_M1_B <- dic.samples(
  model_jagsB,
  n.iter = 10000,
  type = "pD"
)

sum(dic_M1_B$deviance) + sum(dic_M1_B$penalty)



dic_M1_B_sensi <- dic.samples(
  model_jagsB_sensi,
  n.iter = 10000,
  type = "pD"
)

sum(dic_M1_B_sensi$deviance) + sum(dic_M1_B_sensi$penalty)


summary(samplesB)
print(tableau_B, row.names = FALSE)

par(mfrow=c(1,1))

