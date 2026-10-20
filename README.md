# Statistique bayésienne : modélisation de l'infection COVID

## Modélisation bayésienne du risque d'infection par le SARS-CoV-2 chez les contacts de foyers israéliens : rôle de la vaccination, de la taille du foyer et de l'isolement

> **Problématique :** dans quelle mesure la vaccination, la taille du foyer et l'isolement sont-ils associés au risque d'infection des contacts ?

Projet de statistiques bayésiennes (M2) réalisé avec JAGS à partir des données de l'étude de Layan et al. (2022), *Impact of BNT162b2 Vaccination and Isolation on SARS-CoV-2 Transmission in Israeli Households: An Observational Study*, American Journal of Epidemiology.

## Répartition des tâches

| Modèle | Dossier | Variables | Responsable |
|---|---|---|---|
| M0 | `MODELE1_EFFETS_SEPARES/` | Vaccination seule | Ahmad JARKAS |
| M1 | `MODELE1_EFFETS_SEPARES/` | Taille du foyer seule | Ahmad JARKAS |
| M2 | `MODELE2_VACC+FOYER/` | Vaccination + taille du foyer | Duzguncan MERAL |
| M3 | `MODELE3_VACC+FOYER+ISOLEMENT/` | Vaccination + taille du foyer + isolement (modèle final) | Colin PALANDRE |
