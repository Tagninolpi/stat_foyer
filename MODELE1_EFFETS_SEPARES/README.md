
## Objectif

On cherche à étudier le risque d'infection des cas contacts selon trois facteurs :

- **Modèle A** : la vaccination ;
- **Modèle B** : la taille du foyer ;
- **Modèle C** : l'isolement.

Chaque facteur est étudié dans un script différent.

## 1. Fichiers et bibliothèques

Les analyses sont faites sur R avec les bibliothèques suivantes :

```r
library(coda)
library(rjags)
```

Le fichier de données utilisé est `2021_05_14_full_database_2doses.txt`.

Avant de lancer les scripts, il faut modifier `setwd()` pour indiquer le dossier où se trouve le fichier.

## 2. Préparation des données

On commence par importer les données avec `read.table()` et renommer les colonnes.

Ensuite, on conserve uniquement les cas contacts (`index_cases == 0`) et les individus dont le statut vaccinal vaut `0` ou `1`.

```r
df_model <- df[df$index_cases == 0, ]
df_model <- df_model[df_model$vaccination_status %in% c(0, 1), ]
df_model$infected <- ifelse(df_model$infection_status == 0, 0, 1)
```

La variable `infected` vaut `0` si la personne n'est pas infectée et `1` sinon.

## 3. Les trois modèles

Les trois modèles utilisent une régression logistique bayésienne. Pour chaque personne, `p[i]` représente la probabilité d'être infectée.

### Modèle A — Vaccination

```r
logit(p[i]) <- a + b * vaccinated[i]
```

- `a` : valeur de référence pour les non-vaccinés, sur l'échelle du logit.
- `b` : effet de la vaccination sur le logit de la probabilité d'infection.

Priors : `a ~ dnorm(0,1)` et `b ~ dunif(-10,0)`.

### Modèle B — Taille du foyer

```r
logit(p[i]) <- a + c * (size[i] - 2)
```

- `a` : valeur de référence pour un foyer de 2 personnes, sur l'échelle du logit.
- `c` : effet d'une personne supplémentaire dans le foyer.

On soustrait `2` à la taille du foyer pour prendre un foyer de 2 personnes comme référence.

Priors : `a ~ dnorm(0,1)` et `c ~ dnorm(0,4)`.

### Modèle C — Isolement

```r
logit(p[i]) <- a + d * (isolation[i] == 1) + e * (isolation[i] == 2)
```

- `a` : valeur de référence sans isolement, sur l'échelle du logit.
- `d` : effet de l'isolement partiel par rapport à l'absence d'isolement.
- `e` : effet de l'isolement complet par rapport à l'absence d'isolement.

Priors : `a ~ dunif(0,10)`, `d ~ dunif(-10,0)` et `e ~ dunif(-10,0)`.


## 4. Estimation et convergence

Les modèles sont lancés avec `jags.model()`, en utilisant **3 chaînes MCMC**.

On effectue d'abord **2 000 itérations de chauffe**, puis **10 000 itérations par chaîne** avec `coda.samples()`.

Les fonctions utilisées pour vérifier les chaînes sont :

- `summary()` : résumé des paramètres estimés ;
- `gelman.diag()` : vérification de la convergence des chaînes ;
- `effectiveSize()` : nombre effectif d'échantillons ;
- `plot()` : graphiques des chaînes ;
- `autocorr.plot()` : autocorrélation des échantillons.

## 5. Résultats

Après l'estimation, on récupère les paramètres avec `as.matrix()`.

On utilise ensuite :

- `exp()` pour calculer les **odds ratios** à partir des coefficients ;
- `plogis()` pour obtenir les **probabilités d'infection** ;
- `mean()` et `quantile()` pour calculer les moyennes et les intervalles de crédibilité à 95 %.

Pour chaque modèle, un tableau donne les probabilités d'infection selon les différents groupes. Des histogrammes permettent aussi de voir la distribution des paramètres et des probabilités.

## 6. Vérifications prédictives (PPC)

On simule **1 000 fois** les infections avec `rbinom()` à partir des paramètres estimés.

On compare les simulations aux données réelles, pour le nombre total d'infectés et pour le taux d'infection dans chaque groupe.

Sur les histogrammes, la ligne rouge correspond à la valeur observée.

## 7. Analyse de sensibilité

On relance chaque modèle avec d'autres priors pour voir si les résultats changent beaucoup.

| Modèle | Prior initial | Prior alternatif |
|---|---|---|
| A | `b ~ dunif(-10,0)` | `b ~ dnorm(0,1)` |
| B | `c ~ dnorm(0,4)` | `c ~ dnorm(0,0.16)` |
| C | `a ~ dunif(0,10)` ; `d,e ~ dunif(-10,0)` | `a ~ dnorm(0,1)` ; `d ~ dnorm(-0.5,1)` ; `e ~ dnorm(-1,1)` |

On compare ensuite les probabilités d'infection et les odds ratios obtenus avec les deux versions.

## 8. DIC

À la fin de chaque script, on calcule le **DIC** avec `dic.samples()` pour le modèle initial et le modèle avec les nouveaux priors.
