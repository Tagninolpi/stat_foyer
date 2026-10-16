library(coda)
library(rjags)
load.module("dic")

# 1. preparation des données
df <- read.table("DATA/2021_05_14_full_database_2doses.txt")

colnames(df) <- c(
  "individual_id", "house_id", "household_size", "date_symptom",
  "infection_status", "vaccination_status", "end_follow_up",
  "age", "index_cases", "isolation_status"
)

# infecté = 1, non infecté = 0 (dbern accepte que 0 et 1) (1,2,3 = 1) (0=0)
df$infected <- ifelse(df$infection_status == 0, 0, 1)

# comme dans l'étude on étudie que les cas contact (pas de biais)
# comme dans l'étude on retire aussi les 20 individus déja infectés (immunité non liée au vaccin) ???
df <- subset(df, index_cases == 0 & vaccination_status != -1)
dim(df)

df$vaccinated <- df$vaccination_status # (on as enlevé -1)

# comme dans l'étude, isolement partiel et complet sont regroupé
df$isolated   <- ifelse(df$isolation_status == 0, 0, 1)

# on transform la taille du foyer en centrant sur 0 pour que dans a taille= 0 correspond as
# un foyer de taille moyenne et pas un foyer avec 0 personnes
df$size_c     <- df$household_size - mean(df$household_size)

N <- nrow(df)

# 2. Le modele (3 variables)

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

# on fixe les poitns de départ eloigné possible
# une convergence avec des départs aléatoire proche peux etre un faux positif 
init_M3 <- list(
  list(a = -2, b_vacc = -2, b_size = -1, b_iso = -2,),
  list(a =  0, b_vacc =  0, b_size =  0, b_iso =  0,),
  list(a =  2, b_vacc =  2, b_size =  1, b_iso =  2,)
)

model_jags_M3 <- jags.model(
  textConnection(model_M3), data = data_M3, inits = init_M3, n.chains = 3
)

update(model_jags_M3, 2000)   # burn-in

samples_M3 <- coda.samples(
  model_jags_M3,
  variable.names = c("a", "b_vacc", "b_size", "b_iso"),
  n.iter = 10000
)


# 3. Convergence et autocorrelation
for (par_name in c("a", "b_vacc", "b_size", "b_iso")) {
  png(paste0("trace_", par_name, ".png"), width = 1000, height = 500)
  plot(samples_M3[, par_name], main = par_name)
  dev.off()
}

gelman.diag(samples_M3)

for (par_name in c("a", "b_vacc", "b_size", "b_iso")) {
  
  # Gelman-Rubin : un fichier par parametre
  png(paste0("gelman_", par_name, ".png"), width = 800, height = 500)
  gelman.plot(samples_M3[, par_name], main = par_name)
  dev.off()
  
  # Autocorrelation (chaine 1) : un fichier par parametre
  png(paste0("autocorr_", par_name, ".png"), width = 800, height = 500)
  autocorr.plot(samples_M3[, par_name][[1]], main = par_name)
  dev.off()
}
effectiveSize(samples_M3)
"""
gelman diag donne 1 et sur les plots retombe vers 1 au bout de ~4000 ittéretions =>  bonne convergence
les traces des 4 paramètres se superposes => bonne convergence

tailles effectives : 
        a     b_iso    b_size    b_vacc 
 4743.345  4935.909 19267.365 13367.161 

size et vacc retombent à 0 en 2 à 4 pas (pas autocorélation)
iso et a mettent plus de temps (8 à 10 pas => quand meme absence d'autocorélation)
ce qui s'explique par taille effective plus faible
"""

#correlation entre parametres:
round(crosscorr(samples_M3), 2)
# corrélations a posteriori proche de 0 (faible corrélation entre les paramètres)
# sauf entre a et iso -0.75  (normal (on pas un nombre de personnes similaire dans chaque groupe
#cette différence crée une corrélation entre a et iso), pas de corrélation entre variables)
"""       a  b _iso  b_size b_vacc
a       1.00 -0.75  -0.03  -0.32
b_iso  -0.75  1.00  -0.02   0.12
b_size -0.03 -0.02   1.00   0.09
b_vacc -0.32  0.12   0.09   1.00
"""

# 4. Apport des donnees : posterior vs prior

post_M3 <- as.matrix(samples_M3)

for (par_name in c("a", "b_vacc", "b_size", "b_iso")) {
  png(paste0("prior_posterior_", par_name, ".png"), width = 1000, height = 600, res = 120)
  hist(post_M3[, par_name], breaks = 50, prob = TRUE, xlim = c(-5, 5),
       main = paste("Posterior de", par_name), xlab = par_name,
       col = "lightblue", border = "white")
  curve(dnorm(x, 0, 2.5), add = TRUE, lwd = 3, col = "steelblue3")   # prior
  legend("topright", legend = c("Posterior", "Prior N(0 ; 2.5)"),
         fill = c("lightblue", NA), border = c("white", NA),
         lty = c(NA, 1), lwd = c(NA, 3), col = c(NA, "steelblue3"), bty = "n")
  dev.off()
}
# les distributions des postérior (données) sont 9 à 55 fois plus grand que les 
# priors, et etroite et centré => les estimations proviennent des données 
# et non pas du prior.

# 5. Effets ajustes
# sens de l'effet
mean(post_M3[, "b_vacc"] < 0) # =1 => vaccination réduit risque
mean(post_M3[, "b_iso"] < 0) # =1 => isolement réduit risque
mean(post_M3[, "b_size"] > 0) # ~1 => plus grand foyer augmente risque

# taille de l'effet 
exp(apply(post_M3[, c("b_vacc", "b_size", "b_iso")], 2, quantile, c(0.025, 0.5, 0.975)))
# vaccination: odds d’infection sont divisés par environ 4,4 (1 / 0,23= 4,4)
# isolement: odds d’infection sont divisés par environ 5,14 (1 / 0,19= 5,1)
# odds d’infection sont augmenté d'environ 9% par personne suplémentaire
quantile(post_M3[, "b_vacc"], c(0.025, 0.5, 0.975))
'''2.5%        50%      97.5% 
-2.0357354 -1.4845661 -0.9718969 
'''
quantile(post_M3[, "b_iso"],  c(0.025, 0.5, 0.975))
'''2.5%       50%     97.5% 
-1.990271 -1.637842 -1.292367 
''' 
# à comparer avec modèles précédents

# 6. Posterior predictive check

set.seed(123)
idx <- sample(nrow(post_M3), 1000)

group <- factor(paste0("vacc=", df$vaccinated, " iso=", df$isolated))

n_contacts <- tapply(df$infected, df$house_id, length)
obs_inf    <- tapply(df$infected, df$house_id, sum)
k_max      <- max(n_contacts)

obs_total <- sum(df$infected)
obs_none  <- sum(obs_inf == 0)
obs_all   <- sum(obs_inf == n_contacts)
obs_sd    <- sd(obs_inf)
obs_freq  <- table(factor(obs_inf, levels = 0:k_max))

sim_total <- sim_none <- sim_all <- sim_sd <- numeric(length(idx))
sim_freq  <- matrix(NA, nrow = length(idx), ncol = k_max + 1)
sim_group <- matrix(NA, nrow = length(idx), ncol = nlevels(group))

for (k in seq_along(idx)) {
  j <- idx[k]
  p <- plogis(post_M3[j, "a"] + post_M3[j, "b_vacc"] * df$vaccinated +
                post_M3[j, "b_size"] * df$size_c + post_M3[j, "b_iso"] * df$isolated)
  y_sim   <- rbinom(N, size = 1, prob = p)
  inf_sim <- tapply(y_sim, df$house_id, sum)
  
  sim_total[k]  <- sum(y_sim)
  sim_none[k]   <- sum(inf_sim == 0)
  sim_all[k]    <- sum(inf_sim == n_contacts)
  sim_sd[k]     <- sd(inf_sim)
  sim_freq[k, ] <- table(factor(inf_sim, levels = 0:k_max))
  sim_group[k, ] <- tapply(y_sim, group, mean)
}

comparison_group <- data.frame(
  group = levels(group),
  n = as.vector(table(group)),
  observed_p = as.vector(tapply(df$infected, group, mean)),
  mean_p = colMeans(sim_group),
  lower_95 = apply(sim_group, 2, quantile, 0.025),
  upper_95 = apply(sim_group, 2, quantile, 0.975)
)
comparison_group
'''
  group   n observed_p     mean_p   lower_95  upper_95
1 vacc=0 iso=0 210  0.7047619 0.69315714 0.60952381 0.7761905
2 vacc=0 iso=1 332  0.3012048 0.30564759 0.24096386 0.3734940
3 vacc=1 iso=0  42  0.2857143 0.33347619 0.16666667 0.5238095
4 vacc=1 iso=1  83  0.1084337 0.08972289 0.02409639 0.1686747
'''


data.frame(
  caracteristique = c("infectes (total)", "foyers sans infecte",
                      "foyers tous infectes", "sd infectes par foyer"),
  observe = c(obs_total, obs_none, obs_all, obs_sd),
  moyenne_simulee = c(mean(sim_total), mean(sim_none), mean(sim_all), mean(sim_sd)),
  borne_inf_95 = c(quantile(sim_total, .025), quantile(sim_none, .025),
                   quantile(sim_all, .025), quantile(sim_sd, .025)),
  borne_sup_95 = c(quantile(sim_total, .975), quantile(sim_none, .975),
                   quantile(sim_all, .975), quantile(sim_sd, .975)),
  p_bayesien = c(mean(sim_total >= obs_total), mean(sim_none >= obs_none),
                 mean(sim_all >= obs_all), mean(sim_sd >= obs_sd))
)
'''
 caracteristique    observe moyenne_simulee borne_inf_95 borne_sup_95 p_bayesien
1      infectes (total) 269.000000      268.491000   237.000000   301.000000      0.512
2   foyers sans infecte  97.000000       71.738000    59.000000    86.000000      0.000
3  foyers tous infectes  43.000000       27.014000    17.000000    37.000000      0.001
4 sd infectes par foyer   1.615393        1.379335     1.226613     1.548888      0.004
'''

# bonne description moyenne du modèle, mais mauvaisse explication de la répartition
# des infections dans les foyers (le modèle ne prend pas  en compte
# l'états des autres personnes dans le foyer et autres facteurs (effet de l'age))

# 7. DIC
dic_M3 <- dic.samples(model_jags_M3, n.iter = 10000, type = "pD")
dic_M3

'''
Mean deviance:  769.5 
penalty 3.958 
Penalized deviance: 773.4 
'''
