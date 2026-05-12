# SRC_MPCD_projection_matlab

Prototype MATLAB pour le développement, le diagnostic et la validation de méthodes SRC/MPCD quasi-incompressibles fondées sur des projections de vitesse et de flux de masse.

Ce dépôt contient une chaîne expérimentale autour de méthodes SRC/MPCD/SRD pour écoulements périodiques ou en canal, avec un accent actuel sur la réduction des modes compressifs de densité sans imposer une redistribution particulaire dure.

État de référence courant :

```text
Branche de travail : feature/q8-mass-flux-projection
Commit de base Q9 : 14c84fe
Méthode validée actuelle : Q9
Cas validé : Poiseuille, 30000 steps, Nx=32, Ny=16, gamma=20, seed=11
```

---

## 1. Objectif scientifique

L’objectif est de construire une méthode SRC/MPCD hybride qui conserve les avantages du SRC classique — particules, fluctuations, transport mésoscopique — tout en supprimant les modes compressifs incompatibles avec un comportement liquide ou quasi-incompressible.

La contrainte incompressible de base est :

```math
\nabla \cdot u = 0
```

Dans une méthode particulaire à occupation cellulaire fluctuante, cette contrainte seule ne garantit pas que les grandes structures de densité soient maîtrisées. La suite de développement a donc introduit une correction du flux de masse :

```math
\nabla \cdot (N u)
```

où `N` désigne l’occupation particulaire par cellule et `u` la vitesse moyenne cellulaire.

---

## 2. Évolution des méthodes

### Q6 — projection de vitesse

Q6 applique une projection de vitesse sur grille après une étape SRC/MPCD classique.

Schéma conceptuel :

```text
SRC/MPCD classique
-> dépôt particules-grille
-> projection de vitesse : div(u) ≈ 0
-> interpolation de la correction vers les particules
-> thermostat optionnel
```

Q6 stabilise très bien le cas Poiseuille et donne une bonne réduction du transport compressif. En revanche, certains modes basse fréquence de densité restent présents.

---

### Q7 — réparation densitaire / virielle faible

Q7 explore une correction positionnelle et virielle faible afin de réduire les défauts de densité tout en restaurant le champ de vitesse. Cette approche reste intéressante pour une fermeture liquide plus physique, mais elle est plus intrusive et plus délicate à stabiliser.

Q7c a fourni un candidat liquide exploitable, mais Q6 reste la base la plus propre pour une projection cinématique.

---

### Q8 — projection du flux de masse

Q8 introduit une projection du flux de masse `N u`, avec l’idée de contrôler directement la divergence du transport de population.

Formellement, la correction agit sur le flux :

```math
M = N u
```

et cherche une correction de la forme :

```math
M_{\mathrm{new}} = M - D^\top \lambda
```

avec :

```math
D M_{\mathrm{new}} = \mathrm{target}
```

où `D` est l’opérateur de divergence discret.

Cette approche est plus directement liée au transport de densité, mais une correction trop globale peut devenir intrusive et dégrader le champ hydrodynamique.

---

### Q9 — référence actuelle

Q9 est la méthode de référence actuelle.

Elle combine :

```text
Q6 : projection de vitesse div(u) ≈ 0
+
projection low-k du flux de masse div(Nu)
```

Le choix important est que Q9 ne cherche pas à corriger tout le bruit cellulaire de densité. Elle cible seulement les grandes longueurs d’onde de densité, via un filtrage low-k.

Paramètres de référence :

```matlab
params.projectionStrength = 1.0;

params.massFluxProjectionMode = 'relax_to_uniform_lowk';
params.massFluxProjectionStrength = 1.0;
params.massFluxDensityRelaxationBeta = 0.002;

params.massFluxApplyAfterVelocityProjection = true;
params.massFluxTargetFilter = 'lowpass_fft';
params.massFluxLowKMaxIndex = 2;
```

Schéma conceptuel :

```text
SRC/MPCD classique
-> projection de vitesse Q6
-> calcul du défaut de densité N - gamma
-> filtrage low-k du défaut
-> relaxation du flux de masse sur les modes low-k
-> interpolation de la correction vers les particules
-> thermostat optionnel
```

---

## 3. Interprétation de Q9

Q9 doit être interprété comme une méthode SRC/MPCD hybride incompressible effective.

Elle ne prétend pas préserver exactement les coefficients de transport du SRC classique ou de Q6. La modification de `nu_eff` est donc acceptable, à condition que :

1. la solution reste stable ;
2. le profil hydrodynamique reste cohérent ;
3. les modes compressifs basse fréquence soient fortement réduits ;
4. les structures fines, notamment vorticitaires, ne soient pas détruites.

Le point clé est que Q9 agit principalement sur les grandes structures de densité, et non sur le bruit local d’occupation cellulaire.

---

## 4. Résultat de référence Q9 — Poiseuille 30000 steps

Configuration :

```text
Nx = 32
Ny = 16
gamma = 20
seed = 11
initialPopulationMode = exact_per_cell
nSteps = 30000
sampleEvery = 100
bodyForceX = 0.02
wallModeY = thermalize
thermostatAfterProjection = true
```

Comparaison avec Q6 projection seule :

| Métrique | Q9 | Q6 projection seule | Ratio Q9/Q6 |
|---|---:|---:|---:|
| mean std(N) | 3.52220 | 3.54127 | 0.9946 |
| mean out-band | 0.19936 | 0.20255 | 0.9843 |
| mean low-k energy | 4.03686e-5 | 1.74986e-4 | 0.2307 |
| mean rho transport | 9.97653e-3 | 9.44305e-3 | 1.0565 |
| time-avg rel RMS | 0.0221224 | 0.0232799 | 0.9503 |
| nu_eff | 0.0126258 | 0.0201185 | 0.6276 |
| R2 | 0.95636 | 0.96418 | 0.9919 |

Conclusion :

```text
Q9 conserve la qualité locale de densité de Q6,
réduit fortement l'énergie basse fréquence de densité,
et conserve un profil Poiseuille stable.
```

La baisse de `nu_eff` est interprétée comme une modification effective des propriétés de transport du fluide hybride, et non comme un échec de la méthode.

---

## 5. Fichiers principaux

### Scripts de référence

| Fichier | Rôle |
|---|---|
| `run_q9_reference_lowk_mass_flux_30000.m` | Script de référence Q9 validé, Poiseuille 30000 steps |
| `run_compare_density_projection_lowk_mass_flux_poiseuille.m` | Comparaison classic / Q6 projection / Q9 low-k mass-flux |
| `run_compare_density_projection_mass_flux_poiseuille.m` | Comparaison avec projection de flux de masse générale |
| `run_projection_poiseuille_demo.m` | Démonstration Poiseuille courte |
| `run_projection_poiseuille_long_demo.m` | Run Poiseuille long avec analyse viscosité |
| `test_density_homogeneity_lowk_mass_flux_compare_short.m` | Test court de non-régression Q9 |
| `test_projection_mass_flux_grid_only.m` | Test grille seule de la projection flux de masse |

---

### Noyaux de simulation

| Fichier | Rôle |
|---|---|
| `mpcd_step_classic_poiseuille.m` | Étape SRC/MPCD classique en canal |
| `mpcd_step_projection_poiseuille.m` | Étape SRC/MPCD + projection vitesse + option Q9 |
| `mpcd_apply_wall_bc_y.m` | Conditions limites en `y` |
| `projection_initialize_particles.m` | Initialisation particulaire, dont `exact_per_cell` |

---

### Projections

| Fichier | Rôle |
|---|---|
| `projection_project_grid_periodic_x_neumann_y.m` | Projection de vitesse sur domaine périodique en `x`, borné en `y` |
| `projection_project_grid_periodic_fft.m` | Projection périodique FFT |
| `projection_project_mass_flux_periodic_x_neumann_y.m` | Projection du flux de masse `N u`, utilisée par Q9 |
| `projection_interpolate_grid_delta_to_particles.m` | Interpolation de la correction grille vers particules |
| `projection_deposit_particles_to_grid.m` | Dépôt particules-grille |

---

### Diagnostics

| Fichier | Rôle |
|---|---|
| `projection_density_homogeneity_metrics.m` | Métriques de densité : `std(N)`, out-band, low-k energy, RMS |
| `projection_population_diagnostics.m` | Diagnostics d’occupation cellulaire |
| `projection_population_transport_diagnostics.m` | Transport de population |
| `projection_density_transport_continuous_diagnostics.m` | Transport continu de densité |
| `projection_thermal_diagnostics.m` | Diagnostics thermiques |
| `projection_apply_cell_thermostat.m` | Thermostat cellulaire |
| `analyze_projection_poiseuille_viscosity.m` | Estimation de `nu_eff` à partir du profil Poiseuille |
| `plot_density_homogeneity_comparison.m` | Visualisation des diagnostics de densité |

---

### Ancienne chaîne incompressible / liquid closure

Ces fichiers sont conservés comme base historique et pour les futurs tests piston/liquide :

| Fichier | Rôle |
|---|---|
| `mpcd_incompressible_simulator.m` | Ancien simulateur incompressible |
| `mpcd_incompressible_simulator_pistondiag.m` | Ancien simulateur avec diagnostics piston |
| `mpcd_incompressible_step_liquidclosure.m` | Étape avec fermeture liquide |
| `mpcd_apply_liquid_closure_step.m` | Fermeture liquide / réparation |
| `reconstruct_grid_fields_from_runout_piston.m` | Reconstruction pour cas piston |
| `reconstruct_grid_fields_from_runout_liquidclosure*.m` | Reconstructions de champs liquid-closure |
| `run_mpcd_incompressible_main.m` | Ancien script principal incompressible |

Ces fichiers ne sont pas encore la chaîne Q9 proprement dite. Ils serviront pour construire un futur cas piston Q9.

---

## 6. Reproduire le run Q9 de référence

Depuis MATLAB, dans le répertoire du dépôt :

```matlab
clear
clc
close all

run_q9_reference_lowk_mass_flux_30000
```

Le script crée un dossier de sortie de la forme :

```text
q9_reference_lowk_mass_flux_30000_YYYYMMDD_HHMMSS/
```

Ce dossier contient typiquement :

```text
console_log.txt
q9_reference_lowk_mass_flux_30000.mat
q9_reference_lowk_mass_flux_30000_summary.csv
q9_reference_lowk_mass_flux_30000_summary.txt
```

Le fichier `.mat` sauvegarde :

```text
params
outQ9
metricsQ9
summary
refQ6
```

---

## 7. Paramètres essentiels

### Géométrie et population

```matlab
params.Nx = 32;
params.Ny = 16;
params.gamma = 20;
params.initialPopulationMode = 'exact_per_cell';
```

`exact_per_cell` est recommandé pour les validations de référence, car il élimine les différences initiales de population entre runs.

---

### Forçage et murs

```matlab
params.bodyForceX = 0.02;
params.wallModeY = 'thermalize';
```

Le cas Poiseuille courant est un canal :

```text
périodique en x
borné en y
forçage volumique en x
```

---

### Projection de vitesse Q6

```matlab
params.projectionStrength = 1.0;
params.projectedStrength = 1.0;
```

La projection de vitesse est appliquée après l’étape SRC/MPCD classique.

---

### Correction Q9 du flux de masse

```matlab
params.massFluxProjectionMode = 'relax_to_uniform_lowk';
params.massFluxProjectionStrength = 1.0;
params.massFluxDensityRelaxationBeta = 0.002;
params.massFluxApplyAfterVelocityProjection = true;
params.massFluxTargetFilter = 'lowpass_fft';
params.massFluxLowKMaxIndex = 2;
```

Interprétation :

| Paramètre | Sens |
|---|---|
| `massFluxProjectionMode` | Mode de correction du flux de masse |
| `relax_to_uniform_lowk` | Relaxation des modes basse fréquence de `N - gamma` |
| `massFluxProjectionStrength` | Amplitude de la correction interpolée aux particules |
| `massFluxDensityRelaxationBeta` | Intensité de relaxation densitaire |
| `massFluxApplyAfterVelocityProjection` | Applique Q9 après la projection de vitesse Q6 |
| `massFluxTargetFilter` | Filtre appliqué à la cible densitaire |
| `massFluxLowKMaxIndex` | Taille du domaine spectral corrigé |

---

### Diagnostics de densité

```matlab
params.storeDensityMaps = true;
params.densityBandFraction = 0.20;
params.densityExcludeWallCells = 0;
params.lowKMaxIndex = 2;
```

Métriques principales :

| Métrique | Signification |
|---|---|
| `meanStdN` | Écart-type moyen de l’occupation cellulaire |
| `meanOutBandFraction` | Fraction de cellules hors bande `gamma(1 ± densityBandFraction)` |
| `meanLowKEnergy` | Énergie spectrale basse fréquence de `N/gamma - 1` |
| `timeAvgRelRms` | RMS relatif du champ de densité moyenné en temps |
| `meanDensityTransportProjectedRms` | Mesure RMS du transport compressif de densité |
| `nuEff` | Viscosité effective estimée par fit Poiseuille |
| `R2` | Qualité du fit parabolique |

---

## 8. Tests utiles

### Test grille seule de la projection flux de masse

```matlab
test_projection_mass_flux_grid_only
```

Objectif : vérifier l’algèbre de projection du flux `N u` indépendamment du bruit particulaire.

---

### Test court Q9

```matlab
test_density_homogeneity_lowk_mass_flux_compare_short
```

Objectif : vérifier rapidement que la chaîne classic / Q6 / Q9 fonctionne et que les cartes de densité sont bien stockées.

---

### Comparaison Poiseuille classic / Q6 / Q9

```matlab
params = struct();
params.Nx = 32;
params.Ny = 16;
params.gamma = 20;
params.nSteps = 30000;
params.sampleEvery = 100;
params.seed = 11;
params.initialPopulationMode = 'exact_per_cell';

params.projectionStrength = 1.0;

params.massFluxProjectionMode = 'relax_to_uniform_lowk';
params.massFluxProjectionStrength = 1.0;
params.massFluxDensityRelaxationBeta = 0.002;
params.massFluxApplyAfterVelocityProjection = true;
params.massFluxTargetFilter = 'lowpass_fft';
params.massFluxLowKMaxIndex = 2;

params.storeDensityMaps = true;
params.lowKMaxIndex = 2;

cmp = run_compare_density_projection_lowk_mass_flux_poiseuille(params);
```

---

## 9. Points d’attention

### Q9 ne remplace pas encore une fermeture liquide complète

Q9 contrôle des modes cinématiques et densitaires, mais ne définit pas encore une loi d’état liquide explicite.

En particulier, Q9 ne doit pas encore être interprété comme imposant directement :

```math
P = P_{\mathrm{kin}} + P_{\mathrm{vir}}
```

ou un module de compressibilité réaliste.

Pour le futur cas piston, il faudra donc distinguer :

```text
comportement quasi-incompressible cinématique
```

et :

```text
réponse liquide thermodynamique / pression effective
```

---

### La viscosité effective peut changer

La modification de `nu_eff` observée avec Q9 n’est pas considérée comme problématique en soi.

Elle signifie que la méthode hybride modifie les propriétés de transport effectives. La validation doit donc se faire sur l’ensemble :

```text
stabilité
profil moyen
réduction des modes compressifs
préservation des structures fines
diagnostics thermiques
```

et non uniquement sur la conservation de la viscosité Q6.

---

### Le filtrage low-k est central

Le choix `massFluxLowKMaxIndex = 2` est important.

Une correction trop large en nombre d’onde risquerait de lisser le bruit local et de perturber les structures fines. Une correction trop faible laisserait les modes compressifs de grande échelle se développer.

Le compromis actuel est :

```text
corriger les modes cohérents de densité
ne pas corriger brutalement le bruit cellulaire
```

---

## 10. Roadmap

### Étape 1 — Stabiliser Q9 comme référence

- conserver `run_q9_reference_lowk_mass_flux_30000.m` comme script de référence ;
- documenter les résultats Poiseuille 30000 steps ;
- utiliser `seed=11` et `exact_per_cell` pour les comparaisons reproductibles.

---

### Étape 2 — Refactoriser Q9 en brique générique

Objectif futur :

```matlab
[stateOut, diag] = mpcd_apply_q9_projection_channel(stateClassic, params)
```

Cette extraction permettra d’utiliser Q9 en dehors du seul fichier :

```text
mpcd_step_projection_poiseuille.m
```

Elle facilitera ensuite les cas piston et cylindre.

---

### Étape 3 — Cas piston

Objectif : tester le comportement liquide ou quasi-incompressible en compression.

Premiers diagnostics recommandés :

```text
masse totale
std(N)
out-band
low-k density energy
div(u)
div(Nu)
densité moyenne active
énergie cinétique
pression cinétique Pkin
```

Le premier cas piston doit être lent et modéré, par exemple une compression de 5 % à 10 %, afin de tester la réponse quasi-statique avant de passer à des compressions plus fortes.

Point d’interprétation :

```text
Q9 peut réduire les modes compressifs sans encore fournir une vraie EOS liquide.
```

Si la pression piston est un objectif central, il faudra probablement réintroduire une fermeture de pression ou virielle compatible avec Q9.

---

### Étape 4 — Cas cylindre / von Kármán

Objectif : vérifier que Q9 ne détruit pas les structures fines et la vorticité.

Cas minimal recommandé :

```text
canal périodique en x
murs en y
obstacle circulaire
forçage uniforme
comparaison Q6 vs Q9
```

Diagnostics recommandés :

```text
champ de vorticité
enstrophie moyenne
spectre temporel de vorticité derrière le cylindre
présence d’un pic de shedding
énergie low-k de densité
std(N)
out-band
débit moyen
dérive thermique
```

Question physique centrale :

```text
Q9 réduit-il les modes de densité sans effacer la dynamique vorticitaires ?
```

---

## 11. Convention de validation

Pour chaque nouvelle méthode ou nouveau cas, conserver :

```text
1 script de référence clairement nommé
1 dossier de sortie daté
1 fichier .mat complet
1 fichier .csv compact
1 fichier .txt lisible
1 log console
```

Les comparaisons doivent toujours indiquer :

```text
seed
Nx, Ny
gamma
nSteps
sampleEvery
initialPopulationMode
projectionStrength
massFluxProjectionMode
massFluxDensityRelaxationBeta
massFluxLowKMaxIndex
```

---

## 12. Résumé court

Q9 est la référence actuelle du dépôt.

Elle correspond à :

```text
SRC/MPCD classique
+ projection de vitesse Q6
+ relaxation low-k du flux de masse
```

Elle est validée sur Poiseuille 30000 steps avec :

```text
densité locale comparable à Q6
low-k density energy réduite d’environ 4.3x
profil Poiseuille stable
R2 ≈ 0.956
```

La prochaine validation doit porter sur :

```text
1. piston : comportement quasi-incompressible / liquide
2. cylindre von Kármán : préservation des structures fines et de la vorticité
```