


################################
##### Modele 0, isolement:
#####################################


# Modèle JAGS
model_M0_C <- "model {
  a ~ dnorm(0, 1)
  d ~ dnorm(0, 1)
  e ~ dnorm(-1, 1)

  for (i in 1:N) {
    logit(p[i]) <- a + d * (isolation[i] == 1) + e * (isolation[i] == 2)

    infected[i] ~ dbern(p[i])
  }
}"

# Données
data_M0_C <- list(
  N = N,
  infected = df$infected,
  isolation = df$isolation_status
)

# Création du modèle
model_jags2 <- jags.model(
  textConnection(model_M0_C),
  data = data_M0_C,
  n.chains = 3
)

# Burn-in
update(model_jags2, 2000)

# Échantillonnage
samples2 <- coda.samples(
  model_jags2,
  variable.names = c("a", "d", "e"),
  n.iter = 10000
)

#######################################
##### CONVERGENCE
#######################################

summary(samples2)
gelman.diag(samples2)
effectiveSize(samples2)

plot(samples2)
autocorr.plot(samples2)

#######################################
##### DISTRIBUTIONS A POSTERIORI
#######################################

post2 <- as.matrix(samples2)

hist(
  post2[, "a"],
  breaks = 50,
  main = "M0-C - Posterior de a",
  xlab = "a",
  col = "lightblue"
)

hist(
  post2[, "d"],
  breaks = 50,
  main = "M0-C - Posterior de d",
  xlab = "d",
  col = "lightblue"
)

hist(
  post2[, "e"],
  breaks = 50,
  main = "M0-C - Posterior de e",
  xlab = "e",
  col = "lightblue"
)

#######################################
##### ODDS RATIOS A POSTERIORI
#######################################

# Isolement 1 par rapport à 0
quantile(
  exp(post2[, "d"]),
  c(0.025, 0.5, 0.975)
)

# Isolement 2 par rapport à 0
quantile(
  exp(post2[, "e"]),
  c(0.025, 0.5, 0.975)
)

#######################################
##### PROBABILITES A POSTERIORI
#######################################

# Catégorie 0
p2_iso0 <- plogis(post2[, "a"])

# Catégorie 1
p2_iso1 <- plogis(post2[, "a"] + post2[, "d"])

# Catégorie 2
p2_iso2 <- plogis(post2[, "a"] + post2[, "e"])

# Résumés des probabilités
profiles2 <- data.frame(
  isolation_status = c(0, 1, 2),
  
  mean_p = c(
    mean(p2_iso0),
    mean(p2_iso1),
    mean(p2_iso2)
  ),
  
  median_p = c(
    median(p2_iso0),
    median(p2_iso1),
    median(p2_iso2)
  ),
  
  lower_95 = c(
    quantile(p2_iso0, 0.025),
    quantile(p2_iso1, 0.025),
    quantile(p2_iso2, 0.025)
  ),
  
  upper_95 = c(
    quantile(p2_iso0, 0.975),
    quantile(p2_iso1, 0.975),
    quantile(p2_iso2, 0.975)
  )
)

profiles2

#######################################
##### HISTOGRAMMES DES PROBABILITES
#######################################

hist(
  p2_iso0,
  breaks = 50,
  main = "M0-C - Probabilite : isolement 0",
  xlab = "Probabilite d'infection",
  col = "lightblue"
)

hist(
  p2_iso1,
  breaks = 50,
  main = "M0-C - Probabilite : isolement 1",
  xlab = "Probabilite d'infection",
  col = "lightblue"
)

hist(
  p2_iso2,
  breaks = 50,
  main = "M0-C - Probabilite : isolement 2",
  xlab = "Probabilite d'infection",
  col = "lightblue"
)

#######################################
##### OBSERVE VS PREDIT
#######################################

observed2 <- aggregate(
  infected ~ isolation_status,
  data = df,
  FUN = mean
)

names(observed2)[2] <- "observed_p"

comparison2 <- merge(
  observed2,
  profiles2,
  by = "isolation_status"
)

comparison2

barplot(
  rbind(
    comparison2$observed_p,
    comparison2$mean_p
  ),
  beside = TRUE,
  names.arg = comparison2$isolation_status,
  col = c("steelblue", "orange"),
  ylim = c(0, 1),
  xlab = "Categorie d'isolement",
  ylab = "Probabilite d'infection",
  main = "M0-C - Observe vs predit"
)

legend(
  "topright",
  legend = c("Observe", "Predit"),
  fill = c("steelblue", "orange")
)

#######################################
##### POSTERIOR PREDICTIVE CHECK
#######################################

observed <- sum(df$infected)

simulated2 <- numeric(1000)

for (j in 1:1000) {
  
  a <- post2[j, "a"]
  d <- post2[j, "d"]
  e <- post2[j, "e"]
  
  # Probabilité pour chaque contact
  p <- plogis(
    a +
      d * (df$isolation_status == 1) +
      e * (df$isolation_status == 2)
  )
  
  # Simulation des infections
  infected_sim <- rbinom(
    N,
    size = 1,
    prob = p
  )
  
  simulated2[j] <- sum(infected_sim)
}

hist(
  simulated2,
  breaks = 30,
  main = "M0-C - PPC",
  xlab = "Nombre d'infectes simules",
  col = "lightblue"
)

abline(
  v = observed,
  col = "red",
  lwd = 2
)

mean(simulated2)
quantile(simulated2, c(0.025, 0.5, 0.975))
