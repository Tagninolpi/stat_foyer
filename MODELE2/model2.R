donnees <- read.table(file.choose(), header = FALSE)

### PRÉPARATION DES DONNÉES ###

# Donner un nom à chaque des 10 colonnes 
names(donnees) <- c(
  "id", "foyer", "taille_foyer", "jour_infection",
  "statut_infection", "statut_vaccinal", "fin_suivi",
  "classe_age", "index_contact", "isolement"
)

# vérif les colonnes 
head(donnees)

# regrouper le status d'infecton (nouvelle colonne)
donnees$infecte <- ifelse(donnees$statut_infection == 0, 0, 1)

# verif
table(donnees$statut_infection, donnees$infecte)

# voir infections antérieures
contacts <- subset(donnees, index_contact == 0)

table(contacts$statut_vaccinal)

# Données du premier modèle 

analyse <- subset(contacts, statut_vaccinal != -1)
with(analyse, table(statut_vaccinal, infecte))

## donc dans notre échantillon, les vaccinés ont moins souvent été infectés

### PRÉPARER LES INFOS A TRANSMETTRE A JAGS ###

data_jags <- list(
  N = nrow(analyse),
  infecte = analyse$infecte,
  vacc = analyse$statut_vaccinal,
  taille_foyer = analyse$taille_foyer
)

### DÉFINITION DU MODÈLE INTERMÉDIAIRE. ###

modele <- "model {

  a ~ dnorm(0,0.16)
  b ~ dnorm(0,0.16)
  c ~ dnorm(0,0.16)

  for (i in 1:N) {
    logit(p[i]) <- a + b*vacc[i] + c*taille_foyer[i]
    infecte[i] ~ dbern(p[i])
  }

}"

### CHARGER LES OUTILS ET DÉFNIR POINTS DE DÉPART DES CHAINES ###

library(rjags)
library(coda)

# points de départs différents

init <- list(
  list(a = 0,  b = 0,  c = 0),
  list(a = -1, b = -1, c = 0.1),
  list(a = 1,  b = 1,  c = -0.1)
)




### CRÉER DES INFECTIONS SIMULÉES POUR TESTER LE MODÈLE ###

# Fixer générateur
set.seed(123)

# Coefficients fictifs choisis pour le test
a_vrai <- -0.5
b_vrai <- -1
c_vrai <- 0.1

# Probabilité d'infection pour chaque contact
p_sim <- plogis(
  a_vrai +
    b_vrai * analyse$statut_vaccinal +
    c_vrai * analyse$taille_foyer
)

# Simuler le statut infection de chaque contact
infecte_sim <- rbinom(
  n = nrow(analyse),
  size = 1,
  prob = p_sim
)

table(infecte_sim)




### ESTIMER LES COEFFICIENTS À PARTIR DES DONNÉES SIMULÉES ###

# Préparer les données simulées
data_sim <- list(
  N = nrow(analyse),
  infecte = infecte_sim,
  vacc = analyse$statut_vaccinal,
  taille_foyer = analyse$taille_foyer
)

# Ajuster le modèle aux infections simulées
modele_sim <- jags.model(
  file = textConnection(modele),
  data = data_sim,
  inits = init,
  n.chains = 3,
  n.adapt = 1000
)

# Écarter les 5k premières itérations
update(modele_sim, n.iter = 5000)

# Conserver les tirages des coefficients
tirages_sim <- coda.samples(
  model = modele_sim,
  variable.names = c("a", "b", "c"),
  n.iter = 10000
)

# Examiner les chaînes et les estimations
plot(tirages_sim)
gelman.diag(tirages_sim, multivariate = FALSE)
effectiveSize(tirages_sim)
summary(tirages_sim)




### AJUSTER LE MODÈLE AUX VRAIES DONNÉES ###

modele_jags <- jags.model(
  file = textConnection(modele),
  data = data_jags,
  inits = init,
  n.chains = 3,
  n.adapt = 1000
)

# burn in 

update(modele_jags, n.iter = 5000)

# on conserve 10k tirage par chaine

tirages <- coda.samples(
  model = modele_jags,
  variable.names = c("a", "b", "c"),
  n.iter = 10000
)

# vérif convergence

plot(tirages)
gelman.diag(tirages, multivariate = FALSE)
effectiveSize(tirages)

# veif à posteriori 

summary(tirages) #Résumer le posterior de chaque coefficient
echantillons <- as.matrix(tirages) #Regrouper les tirages des trois chaînes dans une matrice
mean(echantillons[, "b"] < 0) #Estimer la probabilité a posteriori que b soit négatif : vaccination associée à une probabilité d'infection plus faible, à taille de foyer identique
mean(echantillons[, "c"] > 0) # Estimer la probabilité a posteriori que c soit positif : taille du foyer associée à une probabilité d'infection plus élevée, à statut vaccinal identique





### CALCUL DU DIC ###

load.module("dic")

dic_M2 <- dic.samples(
  modele_jags,
  n.iter = 10000,
  type = "pD"
)

sum(dic_M2$deviance) + sum(dic_M2$penalty)




### TEST AVEC PRIOR PLUS LARGE ###

modele_large <- "model {

  a ~ dnorm(0, 0.01)
  b ~ dnorm(0, 0.01)
  c ~ dnorm(0, 0.01)

  for (i in 1:N) {
    logit(p[i]) <- a + b*vacc[i] + c*taille_foyer[i]
    infecte[i] ~ dbern(p[i])
  }
}"

# Ajustement aux mêmes données réelles
jags_large <- jags.model(
  file = textConnection(modele_large),
  data = data_jags,
  inits = init,
  n.chains = 3,
  n.adapt = 1000
)

update(jags_large, n.iter = 5000)

tirages_large <- coda.samples(
  model = jags_large,
  variable.names = c("a", "b", "c"),
  n.iter = 10000
)

# Résultats
plot(tirages_large)
gelman.diag(tirages_large, multivariate = FALSE)
effectiveSize(tirages_large)
summary(tirages_large)




### CALCULER LES PROBABILITÉS D’INFECTION PAR PROFIL ###

# Regrouper les tirages des trois chaînes dans une matrice
post <- as.matrix(tirages)

profils <- expand.grid(
  vacc = c(0, 1),
  taille_foyer = sort(unique(analyse$taille_foyer))
)

profils$p_moyenne <- NA_real_
profils$borne_inf_95 <- NA_real_
profils$borne_sup_95 <- NA_real_
profils$effectif <- NA_integer_

for (j in 1:nrow(profils)) {
  
  # Une probabilité pour chaque tirage du posterior
  p <- plogis(
    post[, "a"] +
      post[, "b"] * profils$vacc[j] +
      post[, "c"] * profils$taille_foyer[j]
  )
  
  profils$p_moyenne[j] <- mean(p)
  profils$borne_inf_95[j] <- quantile(p, 0.025)
  profils$borne_sup_95[j] <- quantile(p, 0.975)
  
  # Nombre de contacts réellement observés pour ce profil
  profils$effectif[j] <- sum(
    analyse$statut_vaccinal == profils$vacc[j] &
      analyse$taille_foyer == profils$taille_foyer[j]
  )
}

# résultat
print(profils)




### VÉRIFICATION PRÉDICTIVE PAR STATUT VACCINAL ###

set.seed(456)

# Choisir au hasard 2000 tirages du posterior
indices <- sample(seq_len(nrow(post)), size = 2000)

#Repérer les contacts appartenant à chaque groupe
non_vaccines <- analyse$statut_vaccinal == 0
vaccines <- analyse$statut_vaccinal == 1

#matrice
proportions_sim <- matrix(
  NA_real_, nrow = 2000, ncol = 2,
  dimnames = list(NULL, c("Non vaccines", "Vaccines"))
)

#Réaliser 2000 simulations.
for (s in seq_along(indices)) {
  
  k <- indices[s]
  
  #Probabilité de chaque contact pour ce tirage
  p_rep <- plogis(
    post[k, "a"] +
      post[k, "b"] * analyse$statut_vaccinal +
      post[k, "c"] * analyse$taille_foyer
  )
  
  #Simuler les infections des 667 contacts
  y_rep <- rbinom(nrow(analyse), size = 1, prob = p_rep)
  
  proportions_sim[s, 1] <- mean(y_rep[non_vaccines])
  proportions_sim[s, 2] <- mean(y_rep[vaccines])
}

#Construire le tableau comparant observations et simulations
comparaison_ppc <- data.frame(
  groupe = colnames(proportions_sim),
  proportion_observee = c(
    mean(analyse$infecte[non_vaccines]),
    mean(analyse$infecte[vaccines])
  ),
  moyenne_simulee = colMeans(proportions_sim),
  borne_inf_95 = apply(proportions_sim, 2, quantile, 0.025),
  borne_sup_95 = apply(proportions_sim, 2, quantile, 0.975),
  row.names = NULL
)

print(comparaison_ppc)




### VÉRIFICATION PRÉDICTIVE PAR FOYER ###

set.seed(789)

# Nombre de contacts analysés et d'infectés dans chaque foyer
n_contacts <- as.numeric(table(analyse$foyer))

infectes_obs <- as.numeric(
  tapply(analyse$infecte, analyse$foyer, sum)
)

# Deux caractéristiques observées
observations <- c(
  sum(infectes_obs == 0),
  sum(infectes_obs == n_contacts)
)

# Simuler ces mêmes caractéristiques avec le modèle
stats_foyers <- matrix(
  NA_real_, nrow = length(indices), ncol = 2
)

for (s in seq_along(indices)) {
  
  k <- indices[s]
  
  p_rep <- plogis(
    post[k, "a"] +
      post[k, "b"] * analyse$statut_vaccinal +
      post[k, "c"] * analyse$taille_foyer
  )
  
  y_rep <- rbinom(nrow(analyse), size = 1, prob = p_rep)
  
  infectes_rep <- as.numeric(
    tapply(y_rep, analyse$foyer, sum)
  )
  
  stats_foyers[s, 1] <- sum(infectes_rep == 0)
  stats_foyers[s, 2] <- sum(infectes_rep == n_contacts)
}

# Tableau de comparaison
comparaison_foyers <- data.frame(
  situation = c(
    "Aucun contact infecte",
    "Tous les contacts infectes"
  ),
  nombre_foyers_observe = observations,
  moyenne_simulee = colMeans(stats_foyers),
  borne_inf_95 = apply(stats_foyers, 2, quantile, 0.025),
  borne_sup_95 = apply(stats_foyers, 2, quantile, 0.975)
)

print(comparaison_foyers)

