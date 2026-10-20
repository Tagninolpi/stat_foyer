library(coda)
library(rjags)

# 1. Preparation des donnees

df <- read.table("DATA/2021_05_14_full_database_2doses.txt")

colnames(df) <- c(
  "individual_id", "house_id", "household_size", "date_symptom",
  "infection_status", "vaccination_status", "end_follow_up",
  "age", "index_cases", "isolation_status"
)

# infecté = 1, non infecté = 0 (dbern accepte que 0 et 1) (1, 2, 3 = 1) (0 = 0)
df$infected <- ifelse(df$infection_status == 0, 0, 1)

# comme dans l'étude, on étudie que les cas contact (les index sont infectés par définition => biais)
# comme dans l'étude, on retire aussi les 20 individus déjà infectés (immunité non liée au vaccin)
df <- subset(df, index_cases == 0 & vaccination_status != -1)
dim(df)   # 667 lignes (comme l'étude)

df$vaccinated <- df$vaccination_status   # (on a enlevé -1)

# isolement : 0 = aucun ; 1 = partiel ; 2 = complet
# contrairement à l'étude (qui regroupe partiel et complet), on garde 2 paramètres
# pour voir si l'effet est différent entre isolement partiel et complet
df$isolated <- df$isolation_status

N <- nrow(df)
N   # 667


# 2. Le modèle

# taille du foyer - 1 = nombre d'autres personnes dans le foyer
# => a = logit d'un non vacciné, non isolé, dans un foyer d'une seule personne
# b_vacc = vaccination ; c_size = 1 personne de plus ; d_iso = isolement partiel ; e_iso = isolement complet
model_M3 <- "model {
  a ~ dnorm(0, 0.16)
  b_vacc ~ dnorm(0, 0.16)
  c_size ~ dnorm(0, 0.16)
  d_iso ~ dnorm(0, 0.16)
  e_iso ~ dnorm(0, 0.16)

  for (i in 1:N) {
    logit(p[i]) <- a + b_vacc * vaccinated[i] + c_size * (household_size[i] - 1) + d_iso * (isolated[i] == 1) + e_iso * (isolated[i] == 2)
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

# on fixe les points de départ éloignés
# une convergence avec des départs aléatoires proches peut être un faux positif
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


# 3. Convergence et autocorrelation

for (par_name in c("a", "b_vacc", "c_size", "d_iso", "e_iso")) {
  png(paste0("trace_", par_name, ".png"), width = 1000, height = 500)
  plot(samples_M3[, par_name], main = par_name)
  dev.off()
}

gelman.diag(samples_M3)
# gelman diag = 1 pour les 5 paramètres (multivarié = 1) => bonne convergence
# traces : les 3 chaînes se superposent sans dérive

for (par_name in c("a", "b_vacc", "c_size", "d_iso", "e_iso")) {
  
  # Gelman-Rubin (les 3 chaînes sont déjà comparées) :
  png(paste0("gelman_", par_name, ".png"), width = 800, height = 500)
  gelman.plot(samples_M3[, par_name], main = par_name)
  dev.off()
  
  # Autocorrélation des 3 chaînes (une par panneau) :
  png(paste0("autocorr_", par_name, ".png"), width = 1500, height = 500)
  par(mfrow = c(1, 3))
  for (ch in 1:3) {
    autocorr.plot(samples_M3[[ch]][, par_name], main = paste(par_name, "- chaîne", ch), auto.layout = FALSE)
  }
  par(mfrow = c(1, 1))
  dev.off()
}
# gelman.plot :  pour tous les paramètres le facteur retombe vers 1 rapidement (burn-in suffisant)
# autocorrelation : 
# toutes les courbes retombe vers 0 (a et size sont plus lent que les autres)

effectiveSize(samples_M3)
# tailles effectives (sur 30 000 tirages) :
#         a    b_vacc    c_size     d_iso     e_iso 
#   1740.268 12173.179  1955.836  6893.075  6482.179 

# a et size sont plus lent que les autres

summary(samples_M3)
# Mean      SD  Naive SE Time-series SE
#a       0.44172 0.23615 0.0013634       0.005660
#b_vacc -1.48073 0.27474 0.0015862       0.002504
#c_size  0.08671 0.04484 0.0002589       0.001015
#d_iso  -1.10478 0.22699 0.0013105       0.002740
#e_iso  -1.96627 0.20495 0.0011833       0.002548

#2. Quantiles for each variable:
  
#  2.5%      25%      50%     75%   97.5%
#a      -0.0170363  0.28185  0.43957  0.6020  0.9101
#b_vacc -2.0234627 -1.66557 -1.47879 -1.2924 -0.9518
#c_size -0.0006352  0.05697  0.08676  0.1163  0.1741
#d_iso  -1.5541456 -1.25591 -1.10309 -0.9508 -0.6636
#e_iso  -2.3681491 -2.10215 -1.96489 -1.8280 -1.5687

# erreur de Monte Carlo reste petite (2,4 % de l'écart-type a posteriori)

# moyennes => précises, les ESS sont supérieru à plusieurs centaines => nombre d'itération suffisant

# corrélation entre paramètres :
round(crosscorr(samples_M3), 2)
#           a b_vacc c_size d_iso e_iso
#a       1.00  -0.25  -0.79 -0.33 -0.34
#b_vacc -0.25   1.00   0.08  0.08  0.10
#c_size -0.79   0.08   1.00 -0.04 -0.07
#d_iso  -0.33   0.08  -0.04  1.00  0.41
#e_iso  -0.34   0.10  -0.07  0.41  1.00

# pas de corrélation entre variables (sauf pour d_iso et e_iso (0.41))
# une corrélation existe entre a et les différents paramètres (normale car a est le groupe de référence)

# 4. Apport des données : posterior vs prior

post_M3 <- as.matrix(samples_M3)

for (par_name in c("a", "b_vacc", "c_size", "d_iso", "e_iso")) {
  png(paste0("prior_posterior_", par_name, ".png"), width = 1000, height = 600, res = 120)
  hist(post_M3[, par_name], breaks = 50, prob = TRUE, xlim = c(-5, 5),
       main = paste("Posterior de", par_name), xlab = par_name,
       col = "lightblue", border = "black")
  curve(dnorm(x, 0, 2.5), add = TRUE, lwd = 3, col = "steelblue3")   # prior
  legend("topright", legend = c("Posterior", "Prior N(0 ; 2.5)"),
         fill = c("lightblue", NA), border = c("black", NA),
         lty = c(NA, 1), lwd = c(NA, 3), col = c(NA, "steelblue3"), bty = "n")
  dev.off()
}

# a postériori est beaucoup plus étroit que le prior pour tous
# => les données dominent (les estimations ne viennent pas du prior).


# 5. Effets ajustés

# sens de l'effet
mean(post_M3[, "b_vacc"] < 0)   # = 1 => vaccination réduit le risque
mean(post_M3[, "d_iso"] < 0)    # = 1 => isolement partiel réduit le risque
mean(post_M3[, "e_iso"] < 0)    # = 1 => isolement complet réduit le risque
mean(post_M3[, "c_size"] > 0)   # = 0,98 => un foyer plus grand augmente le risque

# taille de l'effet (odds ratios)
exp(apply(post_M3[, c("b_vacc", "c_size", "d_iso", "e_iso")], 2, quantile, c(0.025, 0.5, 0.975)))
#        b_vacc   c_size     d_iso      e_iso
# 2.5%  0.1314927 1.002023 0.2113478 0.09305805
# 50%   0.2293449 1.090021 0.3338711 0.14058670
# 97.5% 0.3845266 1.190233 0.5201868 0.21074101

# Variable	               OR médian	Odds divisés par
#Vaccination	              0,23	        4,4 (1/0,23=4,4)
#Isolement partiel	        0,33	        3,0
#Isolement complet	        0,14	        7,1
#+1 personne dans le foyer	1,09	     augmenté 9 %

mean(post_M3[, "e_iso"] < post_M3[, "d_iso"])
# 0.9998333 l'isolement complet est plus protecteur que l'isolement partiel
exp(quantile(post_M3[, "e_iso"] - post_M3[, "d_iso"], c(0.025, 0.5, 0.975)))
#  2.5%       50%     97.5% 
#0.2644623 0.4216535 0.6706395 


# 6. Posterior predictive check

set.seed(456)
idx <- sample(nrow(post_M3), 2000)   # tirages au hasard (pas seulement la chaîne 1)

n_contacts <- tapply(df$infected, df$house_id, length)   # contacts par foyer
obs_inf    <- tapply(df$infected, df$house_id, sum)      # infectés par foyer

obs <- c(
  total = sum(df$infected),
  none  = sum(obs_inf == 0),             # foyers sans contact infecté
  all   = sum(obs_inf == n_contacts),    # foyers où tous les contacts sont infectés
  sd    = sd(obs_inf)                    # écart-type du nombre d'infectés par foyer
)

sim <- matrix(NA, nrow = length(idx), ncol = 4,
              dimnames = list(NULL, names(obs)))

for (k in seq_along(idx)) {
  j <- idx[k]
  p <- plogis(post_M3[j, "a"] + post_M3[j, "b_vacc"] * df$vaccinated +
                post_M3[j, "c_size"] * (df$household_size - 1) +
                post_M3[j, "d_iso"] * (df$isolated == 1) +
                post_M3[j, "e_iso"] * (df$isolated == 2))
  y_sim   <- rbinom(N, size = 1, prob = p)
  inf_sim <- tapply(y_sim, df$house_id, sum)
  sim[k, ] <- c(sum(y_sim), sum(inf_sim == 0), sum(inf_sim == n_contacts), sd(inf_sim))
}

ppc_M3 <- data.frame(
  caracteristique = c("infectes (total)", "foyers sans infecte",
                      "foyers tous infectes", "sd infectes par foyer"),
  observe = as.numeric(obs),
  moyenne_simulee = colMeans(sim),
  borne_inf_95 = apply(sim, 2, quantile, 0.025),
  borne_sup_95 = apply(sim, 2, quantile, 0.975),
  p_bayesien = sapply(1:4, function(i) mean(sim[, i] >= obs[i])),
  row.names = NULL
)
ppc_M3

#         caracteristique    observe moyenne_simulee borne_inf_95 borne_sup_95 p_bayesien
#1      infectes (total) 269.000000      269.552500   239.000000   302.000000     0.5210
#2   foyers sans infecte  97.000000       72.744500    59.000000    87.000000     0.0000
#3  foyers tous infectes  43.000000       27.819500    18.000000    38.000000     0.0020
#4 sd infectes par foyer   1.615393        1.391265     1.233034     1.544284     0.0025

# 7. DIC
load.module("dic")
dic_M3 <- dic.samples(model_jags_M3, n.iter = 10000, type = "pD")
dic_M3
#Mean deviance:  757.3 
#penalty 4.932 
#Penalized deviance: 762.2 