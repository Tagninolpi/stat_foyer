############################################################
# MODELE COMPLET : vaccination + taille du foyer + isolement
############################################################
# Question : l'isolement apporte-t-il de l'information EN PLUS de la
# vaccination et de la taille du foyer ?
# Pour y repondre on regarde (1) si l'effet de l'isolement reste fort une fois
# les deux autres variables prises en compte, (2) si le modele s'ajuste mieux
# (DIC), (3) si les predictions collent aux donnees (PPC).

library(coda)
library(rjags)
load.module("dic")  # necessaire pour calculer le DIC avec dic.samples()


################################
##### 1. Preparation des donnees
################################

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

# Variable a expliquer : infecte (1) ou non (0).
# dbern (loi de Bernoulli) n'accepte que 0 ou 1, donc on regroupe tous les
# statuts d'infection non nuls en 1.
df$infected <- ifelse(df$infection_status == 0, 0, 1)

# On ne garde que les contacts (index_cases == 0), car les cas index sont
# infectes par definition et ont presque tous isolation_status = 0 : les garder
# ferait croire que "ne pas s'isoler" cause l'infection.
# On retire aussi les contacts deja infectes l'annee precedente
# (vaccination_status == -1), comme dans l'article.
df <- subset(df, index_cases == 0 & vaccination_status != -1)
dim(df)   # 667 lignes attendues

# Vaccination : 0 = non vaccine, 1 = vaccine (apres le filtre, il ne reste que 0 et 1)
df$vaccinated <- df$vaccination_status

# Isolement : 0 = aucun ; 1 = partiel ; 2 = complet
# On regroupe 1 et 2 (isole) et on laisse 0 seul (non isole), comme dans M_iso.
df$isolated <- ifelse(df$isolation_status == 0, 0, 1)

# Taille du foyer CENTREE : on retire la moyenne.
# Pourquoi : (1) a devient la proba d'un non vaccine, non isole, dans un foyer
# de taille moyenne (au lieu d'un foyer de taille 0, qui n'existe pas),
# (2) a et b_size sont moins correles, donc les chaines melangent mieux.
# Les coefficients b_vacc, b_size et b_iso ne changent pas : seul a change de sens.
df$size_c <- df$household_size - mean(df$household_size)
mean(df$household_size)

N <- nrow(df)
N

# Verifications avant de modeliser
table(df$vaccination_status)               # uniquement 0 et 1
table(df$isolation_status, df$isolated)    # regroupement correct
table(df$vaccinated, df$isolated)
# -> si presque tous les vaccines sont isoles (ou l'inverse), les effets de la
#    vaccination et de l'isolement seront difficiles a separer (grande incertitude)

aggregate(infected ~ vaccinated + isolated, data = df, FUN = mean)
# -> proportions observees d'infectes dans chaque groupe : on les compare plus
#    tard aux predictions du modele


################################
##### 2. Le modele
################################

# Regression logistique bayesienne :
# logit(p) = a + b_vacc * vaccine + b_size * taille_centree + b_iso * isole
# - le logit garantit une proba entre 0 et 1
# - exp(b_vacc), exp(b_size), exp(b_iso) sont des odds ratios, AJUSTES : l'effet
#   de chaque variable quand les deux autres restent identiques
# - prior dnorm(0, 0.16) : precision 0.16, donc ecart-type 2.5 (meme prior que
#   les autres modeles) : large, centre sur 0, donc on ne presuppose aucun effet

model_M3 <- "model {
  a ~ dnorm(0, 0.16)
  b_vacc ~ dnorm(0, 0.16)
  b_size ~ dnorm(0, 0.16)
  b_iso ~ dnorm(0, 0.16)

  for (i in 1:N) {
    logit(p[i]) <- a + b_vacc * vaccinated[i] + b_size * size_c[i] + b_iso * isolated[i]
    infected[i] ~ dbern(p[i])
  }
}"

data_M3 <- list(
  N = N,
  infected = df$infected,
  vaccinated = df$vaccinated,
  size_c = df$size_c,
  isolated = df$isolated
)


################################
##### 3. Estimation par MCMC
################################

# 3 chaines qui partent de points eloignes : si elles finissent par explorer la
# meme zone, c'est un signe de convergence (c'est ce que mesure gelman.diag).
# .RNG.seed differente par chaine : resultats reproductibles.
# La taille du foyer a des valeurs de depart plus petites (+/- 1) car son
# coefficient s'applique a des ecarts jusqu'a plusieurs personnes.
init_M3 <- list(
  list(a = -2, b_vacc = -2, b_size = -1, b_iso = -2,
       .RNG.name = "base::Mersenne-Twister", .RNG.seed = 1),
  list(a =  0, b_vacc =  0, b_size =  0, b_iso =  0,
       .RNG.name = "base::Mersenne-Twister", .RNG.seed = 2),
  list(a =  2, b_vacc =  2, b_size =  1, b_iso =  2,
       .RNG.name = "base::Mersenne-Twister", .RNG.seed = 3)
)

# Creation du sampler
model_jags_M3 <- jags.model(
  textConnection(model_M3),
  data = data_M3,
  inits = init_M3,
  n.chains = 3
)

# burn-in : les premieres iterations ne sont pas conservees, le temps que les
# chaines quittent leur point de depart et rejoignent la zone stable
update(model_jags_M3, 2000)

# echantillon a posteriori : 10000 tirages conserves par chaine
samples_M3 <- coda.samples(
  model_jags_M3,
  variable.names = c("a", "b_vacc", "b_size", "b_iso"),
  n.iter = 10000
)


################################
##### 4. Convergence
################################
# A regarder AVANT d'interpreter quoi que ce soit : des chaines non convergees
# donnent des estimations sans valeur.

summary(samples_M3)

gelman.diag(samples_M3)
# Point est. proche de 1 (< 1.05) = les chaines ont converge vers la meme distribution

effectiveSize(samples_M3)
# nombre de tirages "independants" (les tirages MCMC sont correles) : viser
# plusieurs centaines au minimum

plot(samples_M3)
# traces superposees, sans derive = bon signe

autocorr.plot(samples_M3)
# l'autocorrelation doit retomber rapidement vers 0
# (appuyer sur Entree entre deux graphiques)


################################
##### 5. Distribution a posteriori des parametres
################################

post_M3 <- as.matrix(samples_M3)

for (par in c("a", "b_vacc", "b_size", "b_iso")) {
  hist(
    post_M3[, par],
    breaks = 50,
    main = paste("M3 - Posterior de", par),
    xlab = par,
    col = "lightblue"
  )
}

# Correlation entre les parametres a posteriori.
# Pourquoi : si b_vacc et b_iso sont fortement correles negativement, les donnees
# n'arrivent pas a separer les deux effets (vaccines et isoles se recouvrent).
round(cor(post_M3), 2)


################################
##### 6. Effets ajustes : l'isolement garde-t-il son effet ?
################################

# Probabilites a posteriori du sens de chaque effet
# (b < 0 = reduit le risque ; b_size > 0 = augmente le risque)
mean(post_M3[, "b_vacc"] < 0)
mean(post_M3[, "b_iso"] < 0)
mean(post_M3[, "b_size"] > 0)

# Odds ratios et intervalles credibles a 95 % (lignes : 2.5 %, 50 %, 97.5 %)
# exp(b_size) = effet d'UNE personne de plus dans le foyer
exp(apply(post_M3[, c("b_vacc", "b_size", "b_iso")], 2, quantile, c(0.025, 0.5, 0.975)))

# COMPARAISON AVEC LE MODELE ISOLEMENT SEUL :
# dans M_iso, b = -1.58 et OR = 0.21 (0.15 - 0.29), effet BRUT.
# Si l'OR de b_iso reste proche de 0.21 ici, l'effet de l'isolement tient une fois
# la vaccination et la taille du foyer prises en compte.
# S'il se rapproche de 1, une partie de l'effet brut venait des autres variables
# (confusion).

# Chevauchement des intervalles : l'effet de l'isolement est-il different de 0 ?
quantile(post_M3[, "b_iso"], c(0.025, 0.975))


################################
##### 7. Probabilites d'infection a posteriori par profil
################################

# On applique plogis a CHAQUE tirage (a, b_vacc, b_iso) : cela conserve toute
# l'incertitude et la correlation entre les parametres.
# Taille du foyer fixee a la moyenne (size_c = 0).

profiles_M3 <- expand.grid(vaccinated = c(0, 1), isolated = c(0, 1))

for (k in seq_len(nrow(profiles_M3))) {
  p_k <- plogis(
    post_M3[, "a"] +
      post_M3[, "b_vacc"] * profiles_M3$vaccinated[k] +
      post_M3[, "b_iso"] * profiles_M3$isolated[k]
  )
  profiles_M3$mean_p[k]   <- mean(p_k)
  profiles_M3$lower_95[k] <- quantile(p_k, 0.025)
  profiles_M3$upper_95[k] <- quantile(p_k, 0.975)
}

profiles_M3
# Les profils "non vaccine + isole" et "vaccine + non isole" sont ceux qui montrent
# l'effet propre de chaque protection.


################################
##### 8. Observe vs predit par groupe
################################

# On utilise 2000 tirages choisis au hasard (et non les 2000 premiers, qui
# viendraient tous de la chaine 1)
set.seed(123)
idx <- sample(nrow(post_M3), 2000)

group <- factor(paste0("vacc=", df$vaccinated, " iso=", df$isolated))

# Pour chaque tirage : proba predite de chaque individu, moyenne dans chaque groupe
pred_group <- sapply(idx, function(j) {
  p <- plogis(
    post_M3[j, "a"] +
      post_M3[j, "b_vacc"] * df$vaccinated +
      post_M3[j, "b_size"] * df$size_c +
      post_M3[j, "b_iso"] * df$isolated
  )
  tapply(p, group, mean)
})

comparison_M3 <- data.frame(
  group = levels(group),
  n = as.vector(table(group)),
  observed_p = as.vector(tapply(df$infected, group, mean)),
  mean_p = rowMeans(pred_group),
  lower_95 = apply(pred_group, 1, quantile, 0.025),
  upper_95 = apply(pred_group, 1, quantile, 0.975)
)
comparison_M3
# Attention aux groupes avec peu d'individus (colonne n) : l'observe y est tres
# variable et les comparaisons sont peu informatives.

barplot(
  rbind(comparison_M3$observed_p, comparison_M3$mean_p),
  beside = TRUE,
  names.arg = comparison_M3$group,
  col = c("blue", "orange"),
  ylim = c(0, 1),
  ylab = "Probabilite d'infection",
  main = "M3 - Observe vs predit",
  cex.names = 0.8
)

legend(
  "topright",
  legend = c("Observe", "Predit"),
  fill = c("blue", "orange")
)

# Les memes donnees servent a estimer et a comparer : un bon accord ici ne suffit
# pas a valider le modele, d'ou les verifications par simulation ci-dessous.


################################
##### 9. Verification par simulation (PPC)
################################

# Principe : on simule des infections avec le modele et on regarde si les donnees
# reelles restent dans ce que le modele produit.

idx_ppc <- sample(nrow(post_M3), 1000)

sim_total <- numeric(1000)
sim_none  <- numeric(1000)   # nombre de foyers sans aucun contact infecte
sim_all   <- numeric(1000)   # nombre de foyers ou tous les contacts sont infectes

# Observe, par foyer (limite connue du modele : il traite les contacts comme
# independants alors qu'ils partagent la meme exposition dans le foyer)
n_contacts <- tapply(df$infected, df$house_id, length)
obs_inf    <- tapply(df$infected, df$house_id, sum)

observed_total <- sum(df$infected)
observed_none  <- sum(obs_inf == 0)
observed_all   <- sum(obs_inf == n_contacts)

for (k in 1:1000) {
  j <- idx_ppc[k]
  
  p <- plogis(
    post_M3[j, "a"] +
      post_M3[j, "b_vacc"] * df$vaccinated +
      post_M3[j, "b_size"] * df$size_c +
      post_M3[j, "b_iso"] * df$isolated
  )
  
  infected_sim <- rbinom(N, size = 1, prob = p)
  inf_sim <- tapply(infected_sim, df$house_id, sum)
  
  sim_total[k] <- sum(infected_sim)
  sim_none[k]  <- sum(inf_sim == 0)
  sim_all[k]   <- sum(inf_sim == n_contacts)
}

# 9a. Nombre total d'infectes
hist(
  sim_total,
  breaks = 30,
  main = "M3 - PPC : nombre total d'infectes",
  xlab = "Nombre d'infectes simules",
  col = "lightblue"
)
abline(v = observed_total, col = "red", lwd = 2)

# 9b. Foyers sans aucun infecte
hist(
  sim_none,
  breaks = 30,
  main = "M3 - PPC : foyers sans aucun contact infecte",
  xlab = "Nombre de foyers simules",
  col = "lightblue"
)
abline(v = observed_none, col = "red", lwd = 2)

# 9c. Foyers ou tous les contacts sont infectes
hist(
  sim_all,
  breaks = 30,
  main = "M3 - PPC : foyers ou tous les contacts sont infectes",
  xlab = "Nombre de foyers simules",
  col = "lightblue"
)
abline(v = observed_all, col = "red", lwd = 2)

# Resume : observe, moyenne simulee, et proportion de simulations >= observe
# (valeur proche de 0 ou de 1 = le modele reproduit mal cette caracteristique)
data.frame(
  caracteristique = c("infectes (total)", "foyers sans infecte", "foyers tous infectes"),
  observe = c(observed_total, observed_none, observed_all),
  moyenne_simulee = c(mean(sim_total), mean(sim_none), mean(sim_all)),
  prop_sim_sup_ou_egal = c(
    mean(sim_total >= observed_total),
    mean(sim_none >= observed_none),
    mean(sim_all >= observed_all)
  )
)


################################
##### 10. DIC (ajustement penalise par la complexite)
################################

# DIC = deviance moyenne + penalite (nombre de parametres effectifs).
# Plus il est bas, mieux le modele equilibre ajustement et complexite.
dic_M3 <- dic.samples(model_jags_M3, n.iter = 10000, type = "pD")
dic_M3
sum(dic_M3$deviance) + sum(dic_M3$penalty)

# ATTENTION : le DIC ne se compare qu'entre modeles ajustes sur EXACTEMENT les
# memes lignes. Les modeles vaccination et taille du foyer doivent donc etre
# ajustes sur ces 667 contacts (meme filtre subset) pour etre compares a celui-ci.
# Un ecart d'environ 5 ou plus est generalement considere comme net.