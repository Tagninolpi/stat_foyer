donnees <- read.table(file.choose(), header = FALSE)

# Nomer les colonnes 
names(donnees) <- c(
  "id", "foyer", "taille_foyer", "jour_infection",
  "statut_infection", "statut_vaccinal", "fin_suivi",
  "classe_age", "index_contact", "isolement"
)

# vérif
head(donnees)

# regrouper le status d'infecton
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

# jags

data_jags <- list(
  N = nrow(analyse),
  infecte = analyse$infecte,
  vacc = analyse$statut_vaccinal,
  taille_foyer = analyse$taille_foyer,
  classe_age =  analyse$classe_age
)


modele <- "model {

  a ~ dnorm(0,0.16)
  b ~ dnorm(0,0.16)
  c ~ dnorm(0,0.16)
  d ~ dnorm(0,0.16)


  for (i in 1:N) {
    logit(p[i]) <- a + b*vacc[i] + c*taille_foyer[i] + d*classe_age[i]
    infecte[i] ~ dbern(p[i])
  }

}"

# on lance 3 chaines MCMC

library(rjags)
library(coda)

# points de départs dif

init <- list(
  list(a = 0,  b = 0,  c = 0,    d = 0),
  list(a = -1, b = -1, c = 0.1,  d = -1),
  list(a = 1,  b = 1,  c = -0.1, d = 1)
)

#compilation

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
  variable.names = c("a", "b", "c", "d"),
  n.iter = 10000
)

# vérif convergence

plot(tirages)
gelman.diag(tirages, multivariate = FALSE)
effectiveSize(tirages)

# veif à posteriori 

summary(tirages)
echantillons <- as.matrix(tirages)
mean(echantillons[, "b"] < 0)
echantillons <- as.matrix(tirages)
mean(echantillons[, "d"] > 0)

