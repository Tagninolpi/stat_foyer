library(coda)
library(rjags)

df <- read.table("DATA/2021_05_14_full_database_2doses.txt")

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

# on regroupe 1,2et 3 =>1 pour avoir seulement 0 et 1
# loi de bernoulli n'accepte que 0 et 1
df$infected <- ifelse(df$infection_status == 0, 0, 1) 

# comme dans l'étude de l'article on garde que les cas contacts (index_cases == 0) pour ne pas biaiser
#les résultats
df <- subset(df, index_cases == 0 & vaccination_status != -1)
dim(df)   # 667 lignes (comme dans l'étude de l'article)

# Isolation : 0 = aucun ; 1 = partiel ; 2 = complet
# On regroupe 1 et 2 (isolé), et 0 reste seul (non isolé)
df$isolated <- ifelse(df$isolation_status == 0, 0, 1)

N <- nrow(df)

##### Modele isolement, isolé ou pas:

# b = effet de l'isolement sur le logit de la probabilite d'infection
# prior dnorm(0, 0.16) : precision 0.16, donc ecart-type 2.5 (meme prior que les autres modeles)

model_M_iso <- "model {
  a ~ dnorm(0, 0.16)
  b ~ dnorm(0, 0.16)

  for (i in 1:N) {
    logit(p[i]) <- a + b * isolated[i]
    infected[i] ~ dbern(p[i])
  }
}"


data_M_iso <- list(N = N, infected = df$infected, isolated = df$isolated)

# points de depart eloignes pour les 3 chaines
init_iso <- list(
  list(a = -3, b = -3),
  list(a =  0, b =  0),
  list(a =  3, b =  3)
)
# Creation du sampler
model_jags_iso <- jags.model(
  textConnection(model_M_iso),
  data = data_M_iso,
  inits = init_iso,
  n.chains = 3
)

# le burnin pour que les valeurs convergent: les premieres
#iterations de la chaine MCMC ne sont pas conservees
update(model_jags_iso, 2000)

# mcmc sample
samples_iso <- coda.samples(
  model_jags_iso,
  variable.names = c("a", "b"),
  n.iter = 10000
)

summary(samples_iso)
"""
1. Empirical mean and standard deviation for each variable,
   plus standard error of the mean:

     Mean     SD  Naive SE Time-series SE
a  0.5483 0.1296 0.0007482       0.001814
b -1.5792 0.1713 0.0009890       0.002377

2. Quantiles for each variable:

    2.5%     25%     50%     75%   97.5%
a  0.294  0.4614  0.5486  0.6344  0.8045
b -1.918 -1.6936 -1.5794 -1.4647 -1.2447
"""

gelman.diag(samples_iso)
#Gelman-Rubin	1 pour a, b et en multivarié =>	les chaînes se confondent
effectiveSize(samples_iso)
#Taille effective	environ 5100 (a), 5200 (b)	=> assez de tirages indépendants
plot(samples_iso)
#Traces	superposées, stables =>	pas de dérive ni de blocage

autocorr.plot(samples_iso)
# diminution rapide de l'autocorelation => pas d'autoccorélation


# Distribution a posteriori des parametres
post_iso <- as.matrix(samples_iso)

hist(
  post_iso[, "a"],
  breaks = 50,
  main = "M_iso - Posterior de a",
  xlab = "a",
  col = "lightblue"
)

hist(
  post_iso[, "b"],
  breaks = 50,
  main = "M_iso - Posterior de b",
  xlab = "b",
  col = "lightblue"
)

# Probabilite a posteriori que l'isolement reduise le risque (b < 0)
mean(post_iso[, "b"] < 0)
# Odds ratio de l'isolement (mediane et intervalle credible a 95 %)
exp(quantile(post_iso[, "b"], c(0.025, 0.5, 0.975)))

## Probabilites d'infection a posteriori

# Personnes non isolees
p_iso_non_isolated <- plogis(post_iso[, "a"]) # transforme le logit en proba
mean(p_iso_non_isolated)
quantile(p_iso_non_isolated, c(0.025, 0.5, 0.975))

# Personnes isolees
p_iso_isolated <- plogis(post_iso[, "a"] + post_iso[, "b"])
mean(p_iso_isolated)
quantile(p_iso_isolated, c(0.025, 0.5, 0.975))


## Distribution a posteriori des probabilites

hist(
  p_iso_non_isolated,
  breaks = 50,
  main = "M_iso - Probabilite d'infection : non isoles",
  xlab = "Probabilite d'infection",
  col = "lightblue"
)

hist(
  p_iso_isolated,
  breaks = 50,
  main = "M_iso - Probabilite d'infection : isoles",
  xlab = "Probabilite d'infection",
  col = "lightblue"
)


## Comparaison obs vs predit

# Probabilites observees

# grouper les ind selon leur statut d'isolement
# et calculer la proportion d'infectes dans chaque groupe

observed_iso <- aggregate(
  infected ~ isolated,
  data = df,
  FUN = mean
)

print(observed_iso)
names(observed_iso)[2] <- "observed_p"

# Probabilites predites

profiles_iso <- data.frame(
  isolated = c(0, 1),
  mean_p = c(
    mean(p_iso_non_isolated),
    mean(p_iso_isolated)
  ),
  lower_95 = c(
    quantile(p_iso_non_isolated, 0.025),
    quantile(p_iso_isolated, 0.025)
  ),
  upper_95 = c(
    quantile(p_iso_non_isolated, 0.975),
    quantile(p_iso_isolated, 0.975)
  )
)


# Comparaison: on merge selon le statut d'isolement

comparison_iso <- merge(
  observed_iso,
  profiles_iso,
  by = "isolated"
)

comparison_iso


## Graphique observe vs predit

# beside pour afficher les barres d'un meme groupe cote a cote
barplot(
  rbind(comparison_iso$observed_p, comparison_iso$mean_p),
  beside = TRUE,
  names.arg = c("Non isoles", "Isoles"),
  col = c("blue", "orange"),
  ylim = c(0, 1),
  ylab = "Probabilite d'infection",
  main = "M_iso - Observe vs predit"
)

legend(
  "topright",
  legend = c("Observe", "Predit"),
  fill = c("blue", "orange")
)

# Le modele reproduit-il les proportions observees ?
# Attention : cela ne suffit pas a valider sa qualite predictive, vu qu'il utilise
# les memes donnees pour estimer les parametres et comparer les resultats



## Verification du modele par simulation

observed_total <- sum(df$infected)
simulated_iso <- numeric(1000)

for (j in 1:1000) {
  
  a <- post_iso[j, "a"]
  b <- post_iso[j, "b"]
  
  # Proba des N individus
  p <- plogis(a + b * df$isolated)
  
  # Simuler les infections pour chaque p
  infected_sim <- rbinom(N, size = 1, prob = p)
  
  # Nombre total d'infectes
  simulated_iso[j] <- sum(infected_sim)
}

hist(
  simulated_iso,
  breaks = 30,
  main = "M_iso - PPC",
  xlab = "Nombre d'infectes simules",
  col = "lightblue"
)

abline(v = observed_total, col = "red", lwd = 2)

