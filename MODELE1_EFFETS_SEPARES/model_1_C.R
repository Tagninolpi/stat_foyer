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

table(df_model$isolation_status)

## Param check

# a initial
plogis(qunif(c(0.025, 0.975), min = 0, max = 10))
# d et e initial
exp(qunif(c(0.025, 0.975), min = -10, max = 0))
# a sensi
plogis(qnorm(c(0.025, 0.975), mean = 0, sd = 1))
# d sensi
exp(qnorm(c(0.025, 0.975), mean = -0.5, sd = 1))
# e sensi
exp(qnorm(c(0.025, 0.975), mean = -1, sd = 1))


### Modele 1 C, isolement:


model_M1_C <- "model {
  a ~ dunif(0, 10)
  d ~ dunif(-10,0)
  e ~ dunif(-10, 0)

  for (i in 1:N) {
    logit(p[i]) <- a + 
                   d * (isolation[i] == 1) + 
                   e * (isolation[i] == 2)
    infected[i] ~ dbern(p[i])
  }
}"

data_M1_C <- list(N = N, infected = df_model$infected, isolation = df_model$isolation_status)

model_jagsC <- jags.model(
  textConnection(model_M1_C),
  data = data_M1_C,
  n.chains = 3
)

update(model_jagsC, 2000)

samplesC <- coda.samples(
  model_jagsC,
  variable.names = c("a", "d", "e"),
  n.iter = 10000
)

summary(samplesC)
gelman.diag(samplesC)
effectiveSize(samplesC)
plot(samplesC)
autocorr.plot(samplesC)



### Posterior:

## Distribution a posteriori des parametres

postC <- as.matrix(samplesC)

quantile(postC[, "a"], c(0.025, 0.5, 0.975))
quantile(postC[, "d"], c(0.025, 0.5, 0.975))
quantile(postC[, "e"], c(0.025, 0.5, 0.975))

# OR a posteriori : effet de l'isolement partiel et complet
quantile(exp(postC[, "d"]), c(0.025, 0.5, 0.975))
quantile(exp(postC[, "e"]), c(0.025, 0.5, 0.975))
#quantile(exp(postC[, "e"] - postC[, "d"]), c(0.025, 0.5, 0.975))

hist(
  postC[, "a"],
  breaks = 50,
  main = "Posterior de a (M1 isolement)",
  xlab = "a",
  col = "lightblue"
)

hist(
  postC[, "d"],
  breaks = 50,
  main = "Posterior de d (M1 isolement)",
  xlab = "d",
  col = "lightblue"
)

hist(
  postC[, "e"],
  breaks = 50,
  main = "Posterior de e (M1 isolement)",
  xlab = "e",
  col = "lightblue"
)



### Proba d'infection a posteriori

# proba a posteriori
isolations <- sort(unique(df_model$isolation_status))

pC <- sapply(isolations, function(isolation) {
  plogis(postC[, "a"] +
           postC[, "d"] * (isolation == 1) +
           postC[, "e"] * (isolation == 2))
})

# tableau recap
tableau_C <- data.frame(
  Statut_isolement = isolations,
  
  Probabilite_moyenne = round(100 * colMeans(pC), 1),
  
  IC_95 = sapply(1:ncol(pC), function(i) {
    paste0("[",
           round(100 * quantile(pC[, i], 0.025), 1),
           " % ; ",
           round(100 * quantile(pC[, i], 0.975), 1),
           " %]")
  }),
  
  Effectif = as.integer(table(
    factor(df_model$isolation_status, levels = isolations)
  ))
)

print(tableau_C, row.names = FALSE)



## Distribution a posteriori des probabilites d'infection

# Taux d'infection obs par isolement
observed_C <- sapply(isolations, function(isolation) {
  mean(df_model$infected[df_model$isolation_status == isolation])
})


for (i in seq_along(isolations)) {
  hist(
    pC[, i],
    breaks = 50,
    main = paste("Probabilite d'infection : isolement", isolations[i]),
    xlab = "Probabilite d'infection",
    col = "lightblue"
  )
  
  abline(v = observed_C[i], col = "red", lwd = 2)
}



### PPC

# Nombre total d'infectes obs
observed <- sum(df_model$infected)

# Stocker les simulations
simulatedC <- numeric(1000)
simulatedC_isolations <- matrix(NA, nrow = 1000, ncol = length(isolations))

colnames(simulatedC_isolations) <- isolations

indices <- sample(1:nrow(postC), 1000)

for (j in 1:1000) {
  
  a <- postC[indices[j], "a"]
  d <- postC[indices[j], "d"]
  e <- postC[indices[j], "e"]
  
  # Proba d'infection de chaque contact
  p <- plogis(a +
                d * (df_model$isolation_status == 1) +
                e * (df_model$isolation_status == 2))
  
  # Simuler les infections
  infected_sim <- rbinom(N, size = 1, prob = p)
  
  # Nombre total d'infectes
  simulatedC[j] <- sum(infected_sim)
  
  # taux d'infection selon l'isolement
  simulatedC_isolations[j, ] <- sapply(isolations, function(isolation) {
    mean(infected_sim[df_model$isolation_status == isolation])
  })
}


# Nombre total d'infectés
quantile(simulatedC, c(0.025, 0.5, 0.975))
observed

hist(simulatedC, breaks = 30,
     main = "PPC Nombre d'infectés",
     xlab = "Nombre d'infectés simulés")
abline(v = observed, col = "red", lwd = 2)


# Taux d'infection selon l'isolement

for (i in seq_along(isolations)) {
  
  cat("Isolement :", isolations[i], "\n")
  print(quantile(simulatedC_isolations[, i], c(0.025, 0.5, 0.975)))
  print(observed_C[i])
  
  hist(simulatedC_isolations[, i], breaks = 30,
       main = paste("PPC Isolement", isolations[i]),
       xlab = "Taux d'infection simulé")
  abline(v = observed_C[i], col = "red", lwd = 2)
}


### Sensibilite aux priors


# Priors alternatifs pour d et e :
# lois normales d'ecart type 1 autorisant les deux signes

model_M1_C_sensi <- "model {
  a ~ dnorm(0, 1)
  d ~ dnorm(-0.5, 1)
  e ~ dnorm(-1, 1)

  for (i in 1:N) {
    logit(p[i]) <- a + d * (isolation[i] == 1) + e * (isolation[i] == 2)
    infected[i] ~ dbern(p[i])
  }
}"

model_jagsC_sensi <- jags.model(
  textConnection(model_M1_C_sensi),
  data = data_M1_C,
  n.chains = 3
)

update(model_jagsC_sensi, 2000)

samplesC_sensi <- coda.samples(
  model_jagsC_sensi,
  variable.names = c("a", "d", "e"),
  n.iter = 10000
)


postC_sensi <- as.matrix(samplesC_sensi)

quantile(postC_sensi[, "a"], c(0.025, 0.5, 0.975))
quantile(postC_sensi[, "d"], c(0.025, 0.5, 0.975))
quantile(postC_sensi[, "e"], c(0.025, 0.5, 0.975))

gelman.diag(samplesC_sensi)
effectiveSize(samplesC_sensi)


# Comparer les odds ratio a posteriori
OR_initial_d <- exp(postC[, "d"])
OR_sensi_d <- exp(postC_sensi[, "d"])

OR_initial_e <- exp(postC[, "e"])
OR_sensi_e <- exp(postC_sensi[, "e"])


comparaison_OR <- data.frame(
  Prior = c("d ~ dunif(-10,0)", "d ~ dnorm(-0.5, 1)",
            "e ~ dunif(-10,0)", "e ~ dnorm(-1, 1)"),
  OR_median = c(
    median(OR_initial_d),
    median(OR_sensi_d),
    median(OR_initial_e),
    median(OR_sensi_e)
  ),
  IC_inf = c(
    quantile(OR_initial_d, 0.025),
    quantile(OR_sensi_d, 0.025),
    quantile(OR_initial_e, 0.025),
    quantile(OR_sensi_e, 0.025)
  ),
  IC_sup = c(
    quantile(OR_initial_d, 0.975),
    quantile(OR_sensi_d, 0.975),
    quantile(OR_initial_e, 0.975),
    quantile(OR_sensi_e, 0.975)
  )
)

print(comparaison_OR, row.names = FALSE)


# Comparer les proba d'infection

pC_sensi <- sapply(isolations, function(isolation) {
  plogis(postC_sensi[, "a"] +
           postC_sensi[, "d"] * (isolation == 1) +
           postC_sensi[, "e"] * (isolation == 2))
})

comparaison_proba_C <- data.frame(
  Statut_isolement = isolations,
  
  "Proba_initiale (d,e ~ dunif(-10,0))" = round(100 * colMeans(pC), 1),
  
  "Proba_sensi (d ~ dnorm(-0.5,1) ; e ~ dnorm(-1,1))" =
    round(100 * colMeans(pC_sensi), 1),
  
  check.names = FALSE
)

print(comparaison_proba_C, row.names = FALSE)


### DIC

load.module("dic")

# Modele C

dic_M1_C <- dic.samples(
  model_jagsC,
  n.iter = 10000,
  type = "pD"
)

sum(dic_M1_C$deviance) + sum(dic_M1_C$penalty)



dic_M1_C_sensi <- dic.samples(
  model_jagsC_sensi,
  n.iter = 10000,
  type = "pD"
)

sum(dic_M1_C_sensi$deviance) + sum(dic_M1_C_sensi$penalty)


summary(samplesC)
