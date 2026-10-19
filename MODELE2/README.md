## 1. Préparation des données

Simplificatins : 

- La variable `infecte` regroupe les statuts d’infection en deux catégories : **0 = non infecté ; 1 = infecté**.

- Nous conservons uniquement les contacts (`index_contact == 0`) et excluons les personnes ayant une infection durant l’année précédente (`statut_vaccinal == -1`). 

L’objectif est d’étudier l’infection chez les contacts selon leur vaccination et la taille de leur foyer.

| Statut vaccinal | Non infectés | Infectés | Total | Proportion infectée |
|---|---:|---:|---:|---:|
| Non-vaccinés | 294 | 248 | 542 | 45,8 % |
| Vaccinés | 104 | 21 | 125 | 16,8 % |

## 2. Définition du modèle et initialisation

Nous cherchons à estimer la probabilité d’infection d’un contact en prenant en compte son statut vaccinal et la taille de son foyer.

Pour chaque contact `i`, l’infection suit une loi Bernoulli de probabilité `p[i]`, avec :

logit(p[i]) = a + b × vacc[i] + c × (taille_foyer[i]-1)

- `a` : constante sur l’échelle logit.
- `b` : association avec la vaccination à taille de foyer identique.
- `c` : association avec la taille du foyer à statut vaccinal identique.

Les priots pour a et b sont `dnorm(0, 1)`. Pour c, nous utilisons `dnorm(0, 4)`. 

Les bibliothèques `rjags` et `coda` sont chargées pour réaliser l’estimation et examiner les chaînes MCMC. Trois jeux de valeurs initiales différents sont préparés pour les trois chaînes.

## 3. Création des données simulées

Pour tester l’estimation, nous générons des infections fictives avec des coefficients connus : **a = −0,5 ; b = −1 ; c = 0,1**.

 La fonction `plogis()` calcule leurs probabilités d’infection à partir des coefficients choisis, puis `rbinom(size = 1)` génère un résultat **0 ou 1** pour chaque contact.

| Statut simulé | Nombre de contacts |
|---|---:|
| Non infecté | 371 |
| Infecté | 296 |
| **Total** | **667** |


## 4. Estimation sur les données simulées

Nous ajustons le modèle aux infections simulées pour vérifier s’il retrouve des coefficients compatibles avec ceux utilisés pour les générer.

Trois chaînes sont lancées avec 1000 itérations d’adaptation, puis 2000 itérations de burn-in. Nous conservons ensuite **10000 tirages par chaîne**, soit 30000 tirages au total.

### Diagnostics des chaînes

![Traces des chaînes et distributions a posteriori sur les données simulées](PLOT/PLOT1.png)

Les trois chaînes se superposent sans dérive. 

Les diagnostics de Gelman–Rubin sont arrondis à 1 pour les trois coefficients. Les tailles effectives sont d’environ 1940 pour a, 14398 pour b et 1976 pour c.

Les graphiques suivants montrent l’autocorrélation des trois paramètres dans chacune des trois chaînes.

![Autocorrélation des tirages sur les données simulées — première partie](PLOT/AUTOCOR1A.png)

![Autocorrélation des tirages sur les données simulées — deuxième partie](PLOT/AUTOCOR1B.png)

![Autocorrélation des tirages sur les données simulées — troisième partie](PLOT/AUTOCOR1C.png)

L’axe Lag représente le nombre d’itérations séparant les tirages comparés. L’autocorrélation diminue rapidement vers zéro pour `b`, mais plus lentement pour `a` et `c`. Cela explique leurs tailles effectives plus faibles.

### Comparaison aux coefficients connus

| Coefficient | Valeur utilisée pour simuler | Moyenne estimée | Intervalle de crédibilité à 95 % |
|---|---:|---:|---|
| a | −0,5 | −0,382 | [−0,770 ; 0,009] |
| b | −1 | −1,227 | [−1,696 ; −0,779] |
| c | 0,1 | 0,084 | [0,004 ; 0,165] |

Les trois valeurs connues sont comprises dans leurs intervalles de crédibilité. Les écarts entre valeurs choisies et estimées sont compatibles avec la variabilité des données simulées.

Ce test est rassurant pour le fonctionnement de la procédure d’estimation. 

## 5. Estimation sur les données réelles

### Diagnostics des chaînes

![Traces des chaînes et distributions a posteriori sur les données réelles](PLOT/PLOT2.png)

Les chaînes se superposent sans dérive . Les diagnostics de Gelman Rubin sont arrondis à 1. Les tailles effectives sont d’environ 1806 pour a, 13377 pour b et 1914 pour c (nous avons une bonne convergence).

![Autocorrélation sur les données réelles — première partie](PLOT/AUTOCOR2A.png)

![Autocorrélation sur les données réelles — deuxième partie](PLOT/AUTOCOR2B.png)

![Autocorrélation sur les données réelles — troisième partie](PLOT/AUTOCOR2C.png)

### Estimations et interprétation

| Coefficient | Moyenne estimée | Intervalle de crédibilité à 95 % |
|---|---:|---|
| a - constante | −0,558 | [−0,936 ; −0,180] |
| b - vaccination | −1,318 | [−1,811 ; −0,846] |
| c - taille du foyer | 0,088 | [0,011 ; 0,167] |

- **Vaccination :** à taille de foyer identique, elle est associée à une probabilité d’infection plus faible. Tous les tirages conservés de `b` sont négatifs : la probabilité a posteriori estimée de `b < 0` est proche de 1, sans signifier une certitude absolue.

- **Taille du foyer :** à statut vaccinal identique, une taille plus élevée est associée à une probabilité d’infection plus grande. La probabilité a posteriori estimée de `c > 0` est de **98,92%**.

- **Constante :**  `a` correspond au logit de la probabilité pour un contact non vacciné dans un foyer d’une personne. Cette référence reste théorique, car les foyers analysés comptent au moins deux personnes.

## 6. Comparaison entre prior et posterior

Nous comparons pour chaque coefficient la distribution choisi avant l’analyse. le prior correspond à la courbe bleue et le posterior à l'histogramme en bleu clair.

![Prior et posterior du coefficient a](PLOT/PRIOR1.png)

Les distributions a posteriori sont beaucoup plus concentrées que les priors --> les données ont réduit l’incertitude sur les coefficients.

Le posterior de `b` se concentre sur des valeurs négatives, tandis que celui de `c` se situe principalement sur des valeurs positives. Cela correspind aux résultats précédents concernant la vaccination et la taille du foyer.

## 7. Corrélations a posteriori entre paramètres

Nous utilisons `crosscorr(tirages)` pour examiner comment les coefficients varient ensemble dans les tirages a posteriori.

Une corrélation proche de 1 indique qu’ils tendent à augmenter ensemble,  proche de −1, qu’ils évoluent en sens opposés, proche de 0, qu’ils présentent peu de liaison linéaire.

| Paramètres | Corrélation |
|---|---:|
| a et b | −0,22 |
| a et c | −0,90 |
| b et c | 0,08 |

La constante `a` et le coefficient de taille du foyer `c` présentent une forte corrélation négative. Une augmentation de `c` peut être partiellement compensée par une diminution de `a`, tout en conservant des probabilités d’infection similaires pour les profils observés.

Le coefficient de vaccination `b` présente une faible corrélation avec `a` et une corrélation presque nulle avec `c`.

## 8. Sensibilité au choix des priors

Nous ajustons le même modèle aux mêmes données, en élargissant les priors de `a`, `b` et `c`. L’objectif est de vérifier si les résultats dépendent fortement de ce choix.

Nous conservons les mêmes réglages.

### Diagnostics des chaînes

![Traces des chaînes et distributions a posteriori avec les priors élargis](PLOT/PLOT3.png)

Nous avons une bonne convergence.

### Comparaison des résultats

| Coefficient | Moyenne, prior initial | Moyenne, prior élargi | Intervalle crédible à 95 % avec le prior élargi |
|---|---:|---:|---|
| `a` : constante | −0,558 | −0,571 | [−0,969 ; −0,188] |
| `b` : vaccination | −1,318 | −1,391 | [−1,920 ; −0,902] |
| `c` : taille du foyer | 0,088 | 0,092 | [0,012 ; 0,175] |

Les estimations et les intervalles crédibles restent très proches de ceux obtenus avec les priors initiaux.

## 9. Probabilités d’infection par profil

Nous utilisons les tirages a posteriori du modèle avec les priors initiaux pour répondre à une question concrète : **quelle probabilité d’infection le modèle estime-t-il pour un contact, selon sa vaccination et la taille de son foyer ?**

Pour chaque profil, chaque tirage des coefficients donne une probabilité d’infection. Nous résumons ces probabilités par leur moyenne et leur intervalle crédible à 95 %. L’effectif indique combien de contacts correspondent à ce profil dans les données.

### Résultats pour quelques profils

| Taille du foyer | Statut vaccinal | Probabilité moyenne | Intervalle crédible à 95 % | Effectif |
|---|---|---:|---|---:|
| 2 | Non vacciné | 38,5 % | [31,5 % ; 46,0 %] | 24 |
| 2 | Vacciné | 14,6 % | [9,0 % ; 21,5 %] | 8 |
| 4 | Non vacciné | 42,7 % | [37,9 % ; 47,6 %] | 95 |
| 4 | Vacciné | 16,9 % | [11,1 % ; 23,6 %] | 45 |
| 8 | Non vacciné | 51,4 % | [44,7 % ; 58,2 %] | 55 |
| 8 | Vacciné | 22,4 % | [14,4 % ; 32,0 %] | 8 |
| 12 | Non vacciné | 59,9 % | [46,5 % ; 72,4 %] | 20 |
| 12 | Vacciné | 29,3 % | [16,3 % ; 45,3 %] | 0 |

### Interprétation

À taille de foyer identique, le modèle estime une probabilité d’infection plus faible chez les vaccinés. Par exemple, dans un foyer de quatre personnes, elle est de **16,9 % pour un contact vacciné**, contre **42,7 % pour un contact non vacciné**.

Les probabilités moyennes augmentent avec la taille du foyer dans les deux groupes. 

Le profil « vacciné dans un foyer de 12 personnes » ne correspond à aucun contact observé. Le modèle peut néanmoins calculer une probabilité pour ce profil : il s’agit donc d’une extrapolation.

## 10. Vérification prédictive par statut vaccinal

Nous vérifions si le modèle peut générer des données ressemblant aux observations, concernant la proportion d’infectes chez les vaccinés et les non-vaccinés.

### Démarche

1. Nous choisissons un tirage des coefficients estimés à partir des données réelles.
2. Nous calculons la probabilité d’infection de chaque contact, en conservant sa vaccination et la taille de son foyer.
3. Nous tirons au hasard son statut d’infection, selon cette probabilité.
4. Nous calculons la proportion d’infectés chez les vaccinés et les non-vaccinés.

Nous comparons ensuite ces proportions simulées aux proportions réellement observées.

### Résultats

| Groupe | Proportion observée | Moyenne simulée | Intervalle prédictif à 95 % |
|---|---:|---:|---|
| Non-vaccinés | 45,8 % | 45,6 % | [39,9 % ; 51,5 %] |
| Vaccinés | 16,8 % | 18 % | [9,6 % ; 28 %] |

L’intervalle prédictif contient les 95 % centraux des proportions simulées. Il tient compte de l’incertitude sur les coefficients et du hasard dans les infections générées.

### Interprétation

Les proportions observées sont proches des moyennes simulées et se trouvent dans les intervalles prédictifs. **Le modèle reproduit donc bien les proportions globales d’infection dans ces deux groupes.**


## 11. Vérification prédictive par foyer

Le modèle peut retrouver le bon nombre total d’infectés sans reproduire leur répartition entre les foyers. C'est ce que nous allons vérifier 

### Démarche

Nous simulons 2000 versions des infections des 667 contacts, à partir des coefficients a posteriori. Nous conservons leurs caractéristiques et leur appartenance aux foyers.

Pour chaque simulation, nous comptons les infectés dans chaque foyer, puis calculons :

- Le nombre total d’infectés.
- Le nombre de foyers sans contact infecté.
- Le nombre de foyers où tous les contacts sont infectés.
- L’écart-type du nombre d’infectés par foyer (à quel point ces nombres diffèrent entre les foyers).

### Résultats

| Caractéristique | Valeur observée | Moyenne simulée | Intervalle prédictif à 95 % |
|---|---:|---:|---|
| Total d’infectés | 269 | 269,6 | [235 ; 303] |
| Foyers sans contact infecté | 97 | 62,4 | [48 ; 77] |
| Foyers avec tous les contacts analysés infectés | 43 | 22,3 | [14 ; 32] |
| Écart-type du nombre d’infectés par foyer | 1,62 | 1,26 | [1,10 ; 1,43] |

### Interprétation

Le modèle reproduit bien le total d’infectés : 269 observés contre 269,6 en moyenne dans les simulations.

En revanche, il produit trop peu de foyers sans infection et de foyers où tous les contacts analysés sont infectés. Les nombres observés, respectivement **97 et 43**, dépassent les bornes supérieures des intervalles de confiance prédictifs.

Le modèle répartit les infections de manière trop uniforme entre les foyers. L’écart-type observé, plus élevé que celui généralement simulé, confirme qu’il sous-estime les différences entre foyers.


### Conclusion

**Le modèle 2 reproduit les proportions globales d’infection, mais représente mal leur répartition entre les foyers.**

Il ne tient compte que de la vaccination et de la taille du foyer. Des facteurs absents du modèle pourraient contribuer à ces écarts.

La comparaison avec le modèle 3 permettra d’examiner si l’ajout de l’isolement améliore ces résultats.