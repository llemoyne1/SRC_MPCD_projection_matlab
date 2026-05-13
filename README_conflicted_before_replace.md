# SRC_MPCD_projection_matlab

Prototype MATLAB pour le développement, le diagnostic et la validation de méthodes SRC/MPCD quasi-incompressibles fondées sur des projections de vitesse et de flux de masse.

Ce dépôt contient une chaîne expérimentale autour de méthodes SRC/MPCD/SRD pour écoulements périodiques, en canal, avec piston, obstacle ou géométrie interne. L’objectif actuel est de réduire les modes compressifs de densité sans imposer une redistribution particulaire dure, tout en conservant les structures hydrodynamiques pertinentes.

État de référence courant :

```text
Branche de travail : feature/q8-mass-flux-projection
Méthode actuelle : Q9 low-k mass-flux projection
Base historique Q9 : 14c84fe
Statut : prototype MATLAB validé qualitativement et quantitativement sur plusieurs cas tests
Prochaine étape logique : gel du prototype MATLAB, puis portage C++/OpenMP pour runs plus longs et mieux résolus
```

---

## 1. Objectif scientifique

Le but est de construire une méthode SRC/MPCD hybride qui conserve les propriétés utiles du SRC classique :

```text
particules,
fluctuations mésoscopiques,
transport effectif,
conditions aux limites cinétiques,
```

tout en supprimant les modes compressifs incompatibles avec un comportement liquide ou quasi-incompressible.

La projection incompressible de base impose :

```math
\nabla \cdot u \simeq 0.
```

Cependant, dans une méthode particulaire avec occupation cellulaire fluctuante, cette seule contrainte ne suffit pas à contrôler les grandes structures de densité. La correction Q9 agit donc aussi sur le flux de masse discret :

```math
\nabla \cdot (N u),
```

où `N` est l’occupation particulaire par cellule et `u` la vitesse moyenne cellulaire.

---

## 2. Évolution des méthodes

### Q6 — projection de vitesse

Q6 applique une projection de vitesse après une étape SRC/MPCD classique :

```text
SRC/MPCD classique
-> dépôt particules-grille
-> projection de vitesse : div(u) ≈ 0
-> interpolation de la correction vers les particules
-> thermostat optionnel
```

Q6 stabilise bien Poiseuille et donne une référence cinématique propre, mais certains modes basse fréquence de densité restent présents.

---

### Q7 — réparation densitaire / virielle faible

Q7 explore une correction positionnelle et virielle faible afin de réduire les défauts de densité et de réintroduire une fermeture liquide effective.

Cette approche peut produire une raideur de type loi d’état via un kick de vitesse lié au défaut de densité. Elle est physiquement intéressante, mais plus intrusive et plus délicate à stabiliser.

---

### Q8 — projection du flux de masse

Q8 introduit une projection du flux de masse :

```math
M = N u.
```

L’idée est de corriger directement la divergence du transport de population :

```math
D M_{new} = target.
```

Cette approche est plus directement liée au transport de densité, mais une correction trop globale peut perturber le champ hydrodynamique.

---

### Q9 — référence actuelle

Q9 est la méthode de référence actuelle.

Elle combine :

```text
Q6 : projection de vitesse div(u) ≈ 0
+
relaxation low-k du flux de masse div(Nu)
```

L’idée essentielle est de **ne pas corriger brutalement tout le bruit cellulaire**. Q9 cible surtout les grandes longueurs d’onde de densité, via un filtrage spectral low-k.

Paramètres de référence génériques :

```matlab
params.projectionStrength = 1.0;

params.massFluxProjectionMode = 'relax_to_uniform_lowk';
params.massFluxProjectionStrength = 1.0;
params.massFluxApplyAfterVelocityProjection = true;
params.massFluxTargetFilter = 'lowpass_fft';
params.massFluxLowKMaxIndex = 2;

params.massFluxFinalVelocityProjectionCleanup = true;
params.massFluxFinalVelocityProjectionStrength = 0.5;
```

Le paramètre `massFluxDensityRelaxationBeta` dépend du cas :

```text
Poiseuille / piston EOS : beta typiquement autour de 0.002
Marche / géométries à coin : beta optimal observé autour de 5e-4
```

Schéma conceptuel :

```text
SRC/MPCD classique
-> projection de vitesse Q6
-> calcul du défaut N - gamma
-> filtrage low-k du défaut
-> relaxation du flux de masse sur les modes low-k
-> interpolation de la correction vers les particules
-> cleanup de divergence optionnel
-> thermostat optionnel
```

---

## 3. Interprétation de Q9

Q9 doit être interprété comme une méthode SRC/MPCD hybride quasi-incompressible effective.

Elle ne conserve pas nécessairement les coefficients de transport du SRC classique. Une variation de `nu_eff` n’est donc pas automatiquement un échec. La validation doit porter sur l’ensemble :

```text
stabilité,
profil moyen,
réduction des modes compressifs,
préservation des structures fines,
réponse mécanique effective,
diagnostics thermiques.
```

La méthode agit principalement sur les structures cohérentes de densité, et non sur le bruit local d’occupation cellulaire.

---

## 4. Résultats de validation

### 4.1 Poiseuille — validation de référence

Configuration représentative :

```text
Nx = 32
Ny = 16
gamma = 20
seed = 11
initialPopulationMode = exact_per_cell
nSteps = 30000
bodyForceX = 0.02
wallModeY = thermalize
```

Comparaison Q9 / Q6 projection seule :

| Métrique | Ratio Q9/Q6 | Lecture |
|---|---:|---|
| mean std(N) | ≈ 0.995 | densité locale comparable |
| mean out-band | ≈ 0.984 | légèrement meilleur |
| mean low-k energy | ≈ 0.231 | forte réduction low-k |
| mean rho transport | ≈ 1.057 | transport compressif similaire |
| time-avg rel RMS | ≈ 0.950 | légèrement meilleur |
| nu_eff | ≈ 0.628 | viscosité effective modifiée |
| R2 profil | ≈ 0.992 | profil conservé |

Conclusion : Q9 conserve un profil Poiseuille stable tout en réduisant fortement les modes basse fréquence de densité.

---

### 4.2 Piston — compression quasi-incompressible

Le cas piston a été utilisé pour vérifier la réponse sous compression géométrique lente.

Résultat long représentatif :

```text
nSteps = 30000
compression ≈ 10 %
cleanupStrength = 0.5
```

Comparaison Q9 / Q6 :

| Métrique | Ratio Q9/Q6 | Lecture |
|---|---:|---|
| mean std(N) | ≈ 0.991 | densité locale légèrement meilleure |
| mean out-band | ≈ 0.981 | légèrement meilleur |
| mean low-k energy | ≈ 0.163 | réduction forte des modes low-k |
| mean rho transport | ≈ 1.006 | quasi inchangé |
| mean div(u) après particules | résiduel ≈ 0.026 | compromis cleanup acceptable |

Conclusion : Q9 contrôle fortement les modes compressifs de grande échelle sous compression sans modifier fortement la pression cinétique moyenne.

---

### 4.3 Taylor–Green — préservation d’un vortex lisse

Configuration validante courte :

```text
Nx = Ny = 32
gamma = 20
nSteps = 500
kBT = 0.01
TG amplitude = 0.30
```

Résultat Q9 / Q6 :

| Métrique | Ratio Q9/Q6 | Lecture |
|---|---:|---|
| low-k density energy | ≈ 0.707 | réduction low-k |
| final TG amplitude | ≈ 1.015 | amplitude conservée |
| final TG mode energy | ≈ 1.031 | énergie modale conservée |
| final TG coherence | ≈ 1.038 | cohérence conservée |
| final high-k vel fraction | ≈ 0.979 | pas d’augmentation high-k |

Conclusion : Q9 réduit les modes compressifs sans détruire une structure vorticitiaire lisse.

---

### 4.4 Marche / backward-facing step — séparation et recirculation

Le cas marche est devenu le test structural principal après les difficultés rencontrées sur le cylindre.

Configuration validée :

```text
Nx = 48
Ny = 24
gamma = 20
nSteps = 10000
stepWallMode = bounceback
initialMeanVelocityX = 0.04
meanFlowControlMode = relax_to_target
massFluxDensityRelaxationBeta = 5e-4
massFluxLowKMaxIndex = 2
cleanupStrength = 0.5
```

Résultat Q9 / Q6 :

| Métrique | Ratio Q9/Q6 | Lecture |
|---|---:|---|
| std(N) | ≈ 0.991 | densité locale légèrement meilleure |
| out-band | ≈ 0.977 | meilleur |
| low-k density energy | ≈ 0.655 | réduction nette low-k |
| rho transport | ≈ 1.001 | inchangé |
| div(u) après particules | ≈ 0.840 | meilleure divergence |
| mean enstrophy | ≈ 1.030 | vorticité conservée |
| shear omega RMS | ≈ 0.996 | couche de cisaillement conservée |
| recirc length | ≈ 1.001 | recirculation conservée |
| recirc area | ≈ 1.048 | proche |

Conclusion : Q9 réduit les modes compressifs basse fréquence tout en préservant la recirculation, la couche de cisaillement et la vorticité produites par la marche.

---

### 4.5 Cylindre / von Kármán — non validé à ce stade

Un banc cylindre/von Kármán a été mis en place avec obstacle circulaire, diagnostics de vorticité et estimation de Strouhal.

Les runs courts ont montré que Q9 ne lissait pas le sillage et pouvait réduire les modes low-k. En revanche, les runs longs n’ont pas pu être stabilisés de manière satisfaisante à résolution MATLAB raisonnable.

Problèmes identifiés :

```text
obstacle circulaire sous-résolu,
forte sensibilité au traitement solide/fluide,
instabilité locale près du cylindre,
coût prohibitif des raffinements en MATLAB,
extraction Strouhal non fiable sur fenêtres courtes.
```

Conclusion : le cylindre est mis en pause comme benchmark validant. Il devra être repris dans une version C++/OpenMP plus résolue, avec traitement solide plus robuste.

---

## 5. Diagnostic EOS effective par piston

Un protocole de diagnostic de loi d’état effective a été ajouté.

Trois pressions sont distinguées :

```text
Pkin    = pression cinétique volumique moyenne
Pwall   = pression mécanique au piston par impulsion pariétale
Pexcess = Pwall - Pkin
```

La variable de compression est :

```math
\chi = \frac{\rho}{\rho_0} = \frac{y_{top,0}}{y_{top}}.
```

La forme EOS locale retenue est :

```math
P(\rho) = P(\rho_0) + K_{eff}\ln\left(\frac{\rho}{\rho_0}\right).
```

Pour faibles compressions :

```math
P(\rho) \simeq P(\rho_0) + K_{eff}\left(\frac{\rho}{\rho_0}-1\right).
```

Résultat représentatif du protocole staircase résolu :

```text
rhoRatioList = [1.00 1.05 1.10]
stepsPerPlateau = 6000
sampleEvery = 50
plateauDiscardFraction = 0.5
beta = 0.002
```

Valeurs typiques :

```text
Kkin  ≈ 5.66e3
Kwall ≈ 6.79e3
Kextra = Kwall - Kkin ≈ 1.13e3
```

Interprétation : Q9 induit une raideur mécanique effective mesurable au piston, analogue dans sa fonction au kick viriel, mais obtenue par relaxation du flux de masse et projection.

---

## 6. Sweep de raideur en beta

Un sweep en `massFluxDensityRelaxationBeta` a été introduit pour tester si la compressibilité mécanique effective peut être réglée.

Scripts associés :

```text
run_q9_piston_eos_beta_sweep.m
run_q9_piston_eos_beta_sweep_quick.m
```

Résultat qualitatif :

```text
beta influence fortement Kwall,
mais Kwall(beta) n’est pas une fonction monotone simple sur les runs actuels.
```

Un sweep résolu sur :

```text
betaList = [5e-4, 1e-3, 2e-3]
rhoRatioList = [1.00, 1.05, 1.10]
stepsPerPlateau = 6000
```

a montré que :

```text
la raideur mécanique peut dépasser la raideur cinétique,
les cas beta = 5e-4 et beta = 2e-3 donnent une réponse mécanique nette,
le cas beta = 1e-3 est plus ambigu,
la calibration doit être faite statistiquement par graines et par plages stables.
```

Conclusion : la compressibilité peut être influencée, et probablement calibrée, par les paramètres de relaxation, mais pas encore décrite par une loi monotone simple `K = K(beta)`.

---

## 7. Fichiers principaux

### Références Poiseuille / Q9

| Fichier | Rôle |
|---|---|
| `run_q9_reference_lowk_mass_flux_30000.m` | Référence Q9 Poiseuille longue |
| `run_compare_density_projection_lowk_mass_flux_poiseuille.m` | Comparaison classic / Q6 / Q9 |
| `run_projection_poiseuille_demo.m` | Démonstration Poiseuille courte |
| `run_projection_poiseuille_long_demo.m` | Run long avec analyse viscosité |

### Piston / EOS

| Fichier | Rôle |
|---|---|
| `run_q9_reference_piston.m` | Piston Q9/Q6 long de référence |
| `run_q9_piston_eos_cycle.m` | Cycle piston lent initial |
| `run_q9_piston_eos_cycle_long.m` | Variante longue, coûteuse |
| `run_q9_piston_eos_staircase.m` | Protocole plateaux quasi-statiques |
| `run_q9_piston_eos_staircase_cycle.m` | Variante compression/décompression par plateaux |
| `run_q9_piston_eos_staircase_resolved3.m` | Trois plateaux mieux résolus |
| `run_q9_piston_eos_beta_sweep.m` | Sweep EOS en beta |
| `run_q9_piston_eos_beta_sweep_quick.m` | Sweep rapide |

### Taylor–Green

| Fichier | Rôle |
|---|---|
| `run_q9_reference_taylor_green_short.m` | Validation courte Taylor–Green |
| `run_compare_projection_taylor_green_q6_q9.m` | Comparaison Q6/Q9 |
| `mpcd_apply_q9_projection_periodic.m` | Projection Q9 périodique |
| `projection_taylor_green_diagnostics.m` | Diagnostics vortex/cohérence |

### Marche / step-channel

| Fichier | Rôle |
|---|---|
| `run_q9_reference_step_short.m` | Référence marche ; le nom historique contient encore `short` |
| `run_projection_step_channel_demo.m` | Démo marche Q6 ou Q9 |
| `run_compare_projection_step_q6_q9.m` | Comparaison marche Q6/Q9 |
| `run_q9_reference_step_visual_demo.m` | Démo visuelle avec particules/champs |
| `projection_step_visualize_frame.m` | Affichage frame marche |
| `projection_step_replay_fields.m` | Relecture de champs stockés |
| `projection_step_diagnostics.m` | Diagnostics recirculation/vorticité |

### Cylindre / von Kármán

| Fichier | Rôle |
|---|---|
| `run_q9_reference_vonkarman_short.m` | Banc cylindre court |
| `run_q9_vonkarman_re_st_sweep.m` | Sweep Re/St expérimental |
| `projection_wake_shedding_diagnostics.m` | Diagnostics shedding |
| `projection_vorticity_diagnostics.m` | Diagnostics vorticité |

---

## 8. Visualisation

Une visualisation du cas marche a été réintroduite afin d’illustrer qualitativement les effets de Q9.

Script principal :

```matlab
run_q9_reference_step_visual_demo
```

Affichage :

```text
particules colorées par vitesse,
densité relative N/gamma - 1,
champ de vitesse avec quiver,
vorticité omega.
```

Paramètres utiles :

```matlab
params.visualEvery = 25;
params.visualMaxParticles = 3500;
params.visualPause = 0.01;
params.visualSaveFrames = false;
```

---

## 9. Roadmap immédiate

Avant le passage C++ :

```text
1. Mettre à jour le dépôt GitHub avec un état propre.
2. Compléter le rapport SRC/MPCD avec Q6/Q7/Q8/Q9, validations et EOS effective.
3. Ajouter un cas dam-break MATLAB avec et sans Q9 pour démonstration qualitative.
```

Ensuite :

```text
4. Geler le prototype MATLAB.
5. Porter Q9 vers C++/OpenMP.
6. Reprendre les cas plus résolus : marche affinée, obstacle carré, cylindre/von Kármán.
7. Évaluer la calibration EOS sur runs plus longs et plusieurs graines.
```

---

## 10. Conclusion actuelle

À ce stade, Q9 fournit une régularisation quasi-incompressible efficace du SRC/MPCD :

```text
réduction des modes compressifs low-k,
préservation des structures vorticitaires et séparées,
réponse mécanique effective mesurable au piston,
meilleure parallélisabilité que les méthodes à redistribution dure.
```

La compressibilité mécanique effective peut être influencée par les paramètres de projection, notamment `massFluxDensityRelaxationBeta`, `massFluxLowKMaxIndex` et `cleanupStrength`. Toutefois, la relation entre ces paramètres et le module effectif `Kwall` n’est pas encore monotone ni universelle. Elle doit être calibrée statistiquement.

La prochaine étape technique importante est donc le passage à un code parallèle pour permettre des runs plus longs, mieux résolus et statistiquement plus robustes.
