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
# comme dans l'étude on retire aussi les 20 individus déja infectés (immunité non liée au vaccin)
df <- subset(df, index_cases == 0 & vaccination_status != -1)
dim(df)

df$vaccinated <- df$vaccination_status # (on as enlevé -1)

# comme dans l'étude, on ne regroupe pas les partiel et complet (on met 2 paramètre)
df$isolated   <- df$isolation_status

N <- nrow(df)
N
# 2. Le modele (3 variables)

model_M3 <- "model {
  a ~ dnorm(0, 0.16)
  b_vacc ~ dnorm(0, 0.16)
  c_size ~ dnorm(0, 0.16)
  d_iso ~ dnorm(0, 0.16)
  e_iso ~ dnorm(0, 0.16)

  for (i in 1:N) {
    logit(p[i]) <- a + b_vacc * vaccinated[i] + c_size * (household_size[i]-1) + d_iso * (isolated[i]==1) + e_iso * (isolated[i]==2)
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

# on fixe les points de départ eloigné possible
# une convergence avec des départs aléatoire proche peux etre un faux positif 
init_M3 <- list(
  list(a = -2, b_vacc = -2, c_size = -1, d_iso = -2, e_iso = 0),
  list(a =  0, b_vacc =  0, c_size =  0, d_iso =  0, e_iso = 2),
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


# 3. Convergence et autocorrelation
for (par_name in c("a", "b_vacc", "c_size", "d_iso", "e_iso")) {
  png(paste0("trace_", par_name, ".png"), width = 1000, height = 500)
  plot(samples_M3[, par_name], main = par_name)
  dev.off()
}

gelman.diag(samples_M3)

for (par_name in c("a", "b_vacc", "c_size", "d_iso", "e_iso")) {
  
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
  a       b_vacc     c_size     d_iso     e_iso 
1525.888 12036.266  1859.648  5151.273 20251.230 

size et vacc retombent à 0 en 2 à 4 pas (pas autocorélation)
iso et a mettent plus de temps (8 à 10 pas => quand meme absence d'autocorélation)
ce qui s'explique par taille effective plus faible
"""
summary(samples_M3)


#correlation entre parametres:
round(crosscorr(samples_M3), 2)
# corrélations a posteriori proche de 0 (faible corrélation entre les paramètres)
# sauf entre a et iso -0.75  (normal (on pas un nombre de personnes similaire dans chaque groupe
#cette différence crée une corrélation entre a et iso), pas de corrélation entre variables)
"""
           a b_vacc c_size d_iso e_iso
a       1.00  -0.27  -0.80 -0.45  0.00
b_vacc -0.27   1.00   0.09  0.12 -0.01
c_size -0.80   0.09   1.00  0.00  0.00
d_iso  -0.45   0.12   0.00  1.00  0.00
e_iso   0.00  -0.01   0.00  0.00  1.00
"""

# 4. Apport des donnees : posterior vs prior

post_M3 <- as.matrix(samples_M3)

for (par_name in c("a", "b_vacc", "c_size", "d_iso", "e_iso")) {
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
mean(post_M3[, "d_iso"] < 0) # =1 => isolement réduit risque
mean(post_M3[, "e_iso"] < 0)# =0,498 => isolement réduit risque
mean(post_M3[, "c_size"] > 0) # ~9,8 => plus grand foyer augmente risque

# taille de l'effet 
exp(apply(post_M3[, c("b_vacc", "c_size", "d_iso", "e_iso")], 2, quantile, c(0.025, 0.5, 0.975)))
# vaccination: odds d’infection sont divisés par environ 4,4 (1 / 0,23= 4,4)
# isolement: odds d’infection sont divisés par environ 5,14 (1 / 0,19= 5,1)
# odds d’infection sont augmenté d'environ 9% par personne suplémentaire
quantile(post_M3[, "b_vacc"], c(0.025, 0.5, 0.975))
'''2.5%        50%      97.5% 
-2.0336178 -1.4856816 -0.9665277 
'''
quantile(post_M3[, "d_iso"],  c(0.025, 0.5, 0.975))
'''     2.5%       50%     97.5% 
-1.990509 -1.632518 -1.283018 
''' 
quantile(post_M3[, "e_iso"],  c(0.025, 0.5, 0.975))
'''       2.5%         50%       97.5% 
-4.83664397  0.01465433  4.91245252 
''' 

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
