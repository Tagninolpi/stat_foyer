library(coda)
library(rjags)
load.module("dic")

#1. Préparation des données

df <- read.table("DATA/2021_05_14_full_database_2doses.txt")

colnames(df) <- c(
  "individual_id", "house_id", "household_size", "date_symptom",
  "infection_status", "vaccination_status", "end_follow_up",
  "age", "index_cases", "isolation_status"
)

# on regroupe tous les statuts infectés (1,2,3) et on laisse non-infecté (0) seul
# => infected => 0 ou 1 (pour Bernoulli)
df$infected <- ifelse(df$infection_status == 0, 0, 1)

# Comme dans l'étude :
# on ne garde que les cas contact (les cas index sont infectés par définition => biais)
# on retire les 20 individus déjà infectés (immunité non liée au vaccin)
df <- subset(df, index_cases == 0 & vaccination_status != -1)
dim(df)   # 667 lignes (comme l'étude)

df$vaccinated <- df$vaccination_status   # 0 = non vacciné ; 1 = vacciné (0 ou 1 pour Bernoulli)

# Isolement : 0 = aucun ; 1 = partiel ; 2 = complet
# l'étude regroupe isolation partiel et complet
# nous ne les regroupons pas car on veut voir s'ils ont des effets différents
df$isolated <- df$isolation_status

N <- nrow(df)
N


# 2. Modèle

# Le plus petit foyer est composé d'au moins 2 personnes (cas contact et cas index (testé positif))
# pour que l'intercept ait une valeur interprétable on centre sur 2 (household_size -2)

# a      : logit d'un contact non vacciné, non isolé, dans un foyer de 2 personnes
# b_vacc : effet de la vaccination
# c_size : effet d'une personne de plus dans le foyer
# d_iso  : effet de l'isolement partiel
# e_iso  : effet de l'isolement complet

model_M3 <- "model {
  a      ~ dnorm(0, 1)
  b_vacc ~ dnorm(0, 1)
  c_size ~ dnorm(0, 4)
  d_iso  ~ dnorm(0, 1)
  e_iso  ~ dnorm(-1, 1)

  for (i in 1:N) {
    logit(p[i]) <- a + b_vacc * vaccinated[i] + c_size * (household_size[i] - 2) +
                   d_iso * (isolated[i] == 1) + e_iso * (isolated[i] == 2)
    infected[i] ~ dbern(p[i])
  }
}"

data_M3 <- list(
  N = N,
  infected = df$infected,
  vaccinated = df$vaccinated,
  household_size = df$household_size,
  isolated = df$isolated
)

# Points de départ éloignés pour les 3 chaînes : 
# une convergence depuis des départs proches peut être un faux positif
init_M3 <- list(
  list(a = -2, b_vacc = -2, c_size = -1, d_iso = -2, e_iso =  0),
  list(a =  0, b_vacc =  0, c_size =  0, d_iso =  0, e_iso =  2),
  list(a =  2, b_vacc =  2, c_size =  1, d_iso =  2, e_iso = -2)
)

model_jags_M3 <- jags.model(
  textConnection(model_M3), data = data_M3, inits = init_M3, n.chains = 3
)

update(model_jags_M3, 2000)   # burn-in

samples_M3 <- coda.samples(
  model_jags_M3,
  variable.names = c("a", "b_vacc", "c_size", "d_iso", "e_iso"),
  n.iter = 10000
)

params <- c("a", "b_vacc", "c_size", "d_iso", "e_iso")


# 3. Convergence et autocorrélation

# Traces
for (par_name in params) {
  png(paste0("trace_", par_name, ".png"), width = 1000, height = 500)
  plot(samples_M3[, par_name], main = par_name)
  dev.off()
}

# Gelman-Rubin (valeurs proches de 1 = chaînes cohérentes entre elles)
gelman.diag(samples_M3)

for (par_name in params) {
  png(paste0("gelman_", par_name, ".png"), width = 800, height = 500)
  gelman.plot(samples_M3[, par_name], main = par_name)
  dev.off()
  
  # Autocorrélation des 3 chaînes (une par panneau)
  png(paste0("autocorr_", par_name, ".png"), width = 1500, height = 500)
  par(mfrow = c(1, 3))
  for (ch in 1:3) {
    autocorr.plot(samples_M3[[ch]][, par_name],
                  main = paste(par_name, "- chaîne", ch), auto.layout = FALSE)
  }
  par(mfrow = c(1, 1))
  dev.off()
}

# Taille d'échantillon effective (nombre de tirages "indépendants" équivalents)
effectiveSize(samples_M3)

# Résumé a posteriori : moyennes, écarts-types, quantiles
summary(samples_M3)

# Erreur de Monte Carlo relative = Time-series SE / SD a posteriori
stats_M3 <- summary(samples_M3)$statistics
round(stats_M3[, "Time-series SE"] / stats_M3[, "SD"], 3)

# Corrélations a posteriori entre paramètres
round(crosscorr(samples_M3), 2)



# 4. Apport des données : posterior vs prior

# priors : (moyenne, précision JAGS) -> sd = 1 / sqrt(précision)
priors <- list(
  a      = c(mean =  0, prec = 1),
  b_vacc = c(mean =  0, prec = 1),
  c_size = c(mean =  0, prec = 4),
  d_iso  = c(mean =  0, prec = 1),
  e_iso  = c(mean = -1, prec = 1)
)

post_M3 <- as.matrix(samples_M3)

xlims <- list(a = c(-5, 5), b_vacc = c(-5, 5), c_size = c(-1.5, 1.5),
              d_iso = c(-5, 5), e_iso = c(-5, 5))

for (par_name in names(priors)) {
  m        <- unname(priors[[par_name]]["mean"])
  sd_prior <- 1 / sqrt(unname(priors[[par_name]]["prec"]))
  
  png(paste0("prior_posterior_", par_name, ".png"), width = 1000, height = 600, res = 120)
  hist(post_M3[, par_name], breaks = 50, prob = TRUE, xlim = xlims[[par_name]],
       main = paste("Posterior de", par_name), xlab = par_name,
       col = "lightblue", border = "black")
  curve(dnorm(x, m, sd_prior), add = TRUE, lwd = 3, col = "steelblue3")   # prior
  legend("topright",
         legend = c("Posterior", paste0("Prior N(", m, " ; ", sd_prior, ")")),
         fill = c("lightblue", NA), border = c("black", NA),
         lty = c(NA, 1), lwd = c(NA, 3), col = c(NA, "steelblue3"), bty = "n")
  dev.off()
}



# 5. Effets ajustés


# Tableau des coefficients (échelle logit) : moyenne et intervalle crédible à 95 %
coefs <- t(sapply(params, function(p) {
  c(moyenne = mean(post_M3[, p]),
    inf = unname(quantile(post_M3[, p], 0.025)),
    sup = unname(quantile(post_M3[, p], 0.975)))
}))
round(coefs, 3)

# Sens de l'effet : probabilité a posteriori que le coefficient soit du côté attendu
mean(post_M3[, "b_vacc"] < 0)   # vaccination protectrice
mean(post_M3[, "d_iso"]  < 0)   # isolement partiel protecteur
mean(post_M3[, "e_iso"]  < 0)   # isolement complet protecteur
mean(post_M3[, "c_size"] > 0)   # foyer plus grand => risque plus élevé

# Taille de l'effet : odds ratios (médiane [IC crédible à 95 %])
exp(apply(post_M3[, c("b_vacc", "c_size", "d_iso", "e_iso")], 2,
          quantile, c(0.025, 0.5, 0.975)))

# Isolement complet vs partiel
mean(post_M3[, "e_iso"] < post_M3[, "d_iso"])
exp(quantile(post_M3[, "e_iso"] - post_M3[, "d_iso"], c(0.025, 0.5, 0.975)))



# 6. Impact de l'isolement à taille de foyer et statut vaccinal égaux

# Dans des foyers de même taille, et à statut vaccinal égal,
#l'isolement (partiel ou complet) change-t-il la probabilité d'être infecté ?

# On cherche la taille de foyer avec le plus d'effectif
table(df$household_size)
# 4 et 5 ont le plus d'effectif, on compare les proportions de vaccinés et d'isolement
# (on cherche des proportions les plus égales possible)
sub <- subset(df, household_size %in% c(4, 5))

# Proportions de vaccination
table(taille = sub$household_size, vacc = sub$vaccinated)
round(prop.table(table(taille = sub$household_size, vacc = sub$vaccinated), 1), 2)
# 4 : 70-30
# 5 : 80-20

# Proportions d'isolement
table(taille = sub$household_size, iso = sub$isolated)
round(prop.table(table(taille = sub$household_size, iso = sub$isolated), 1), 2)
# 35, 15, 50 pour les 2

# Les proportions de vaccination et d'isolement sont globalement similaires,
# on prend donc une taille de foyer de 5 pour avoir un plus grand effectif

# Probabilités prédites : isolement variable, taille 5 et vaccination fixées
taille_ref <- 5

resume <- function(x) c(moy = mean(x),
                        inf = unname(quantile(x, 0.025)),
                        sup = unname(quantile(x, 0.975)))

lignes <- list()
for (v in c(0, 1)) {
  lp0 <- post_M3[, "a"] + post_M3[, "b_vacc"] * v + post_M3[, "c_size"] * (taille_ref - 2)
  p_aucun   <- plogis(lp0)
  p_partiel <- plogis(lp0 + post_M3[, "d_iso"])
  p_complet <- plogis(lp0 + post_M3[, "e_iso"])
  lignes[[length(lignes) + 1]] <- c(
    vacc = v, taille = taille_ref,
    p_aucun = resume(p_aucun),
    p_partiel = resume(p_partiel),
    p_complet = resume(p_complet),
    baisse_partiel = resume(p_aucun - p_partiel),
    baisse_complet = resume(p_aucun - p_complet)
  )
}
resultats_iso <- as.data.frame(do.call(rbind, lignes))
round(resultats_iso, 2)

# 7. Posterior predictive check (PPC) par foyer

# Principe : simuler des jeux de données à partir du modèle ajusté.
# Comparer des caractéristiques simulées à celles observées.
# Question : le modèle reproduit-il la répartition des infections entre foyers ?

set.seed(456)
idx <- sample(nrow(post_M3), 2000)   # tirages au hasard parmi les 3 chaînes

n_contacts <- tapply(df$infected, df$house_id, length)   # contacts par foyer
obs_inf    <- tapply(df$infected, df$house_id, sum)      # infectés par foyer

obs_foyer <- c(
  total = sum(df$infected),
  none  = sum(obs_inf == 0),             # foyers sans contact infecté
  all   = sum(obs_inf == n_contacts),    # foyers où tous les contacts sont infectés
  sd    = sd(obs_inf)                    
)

sim_foyer <- matrix(NA, nrow = length(idx), ncol = 4,
                    dimnames = list(NULL, names(obs_foyer)))

for (k in seq_along(idx)) {
  j <- idx[k]
  p_sim <- plogis(post_M3[j, "a"] + post_M3[j, "b_vacc"] * df$vaccinated +
                    post_M3[j, "c_size"] * (df$household_size - 2) +
                    post_M3[j, "d_iso"] * (df$isolated == 1) +
                    post_M3[j, "e_iso"] * (df$isolated == 2))
  y_sim   <- rbinom(N, size = 1, prob = p_sim)
  inf_sim <- tapply(y_sim, df$house_id, sum)
  sim_foyer[k, ] <- c(sum(y_sim), sum(inf_sim == 0),
                      sum(inf_sim == n_contacts), sd(inf_sim))
}

ppc_M3 <- data.frame(
  caracteristique = c("infectes (total)", "foyers sans infecte",
                      "foyers tous infectes", "sd infectes par foyer"),
  observe = as.numeric(obs_foyer),
  moyenne_simulee = colMeans(sim_foyer),
  borne_inf_95 = apply(sim_foyer, 2, quantile, 0.025),
  borne_sup_95 = apply(sim_foyer, 2, quantile, 0.975),
  p_bayesien = sapply(1:4, function(i) mean(sim_foyer[, i] >= obs_foyer[i])),
  row.names = NULL
)
ppc_M3



# 8. DIC :
dic_M3 <- dic.samples(model_jags_M3, n.iter = 10000, type = "pD")
dic_M3
