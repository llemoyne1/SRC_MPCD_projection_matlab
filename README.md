# SRC_MPCD_projection_matlab

Prototype MATLAB pour le développement, le diagnostic et la validation de méthodes SRC/MPCD quasi-incompressibles fondées sur des projections de vitesse et de flux de masse.

Ce dépôt contient une chaîne expérimentale autour de méthodes SRC/MPCD/SRD pour écoulements périodiques ou en canal, avec un accent actuel sur la réduction des modes compressifs de densité sans imposer une redistribution particulaire dure.

État de référence courant :

```text
Branche de travail historique : feature/q8-mass-flux-projection
Commit de base Q9 : 14c84fe
Méthode de référence actuelle : Q9
Validation acquise :
  - Poiseuille 30000 steps
  - piston quasi-incompressible 30000 steps
  - Taylor-Green vortex court 500 steps
  - marche / step-channel 10000 steps
Validation non stabilisée à ce stade :
  - cylindre / von Kármán long
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

Le critère physique n’est donc pas uniquement de réduire `div(u)`, mais de limiter le transport compressif cohérent de population tout en préservant les structures de vitesse.

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
+
reprojection finale optionnelle de vitesse
```

Le choix important est que Q9 ne cherche pas à corriger tout le bruit cellulaire de densité. Elle cible seulement les grandes longueurs d’onde de densité, via un filtrage low-k.

Paramètres génériques de référence :
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

Pour les cas piston, Taylor-Green et les cas structurés récents, le compromis courant utilise en plus :

```matlab
params.massFluxFinalVelocityProjectionCleanup = true;
params.massFluxFinalVelocityProjectionStrength = 0.5;
```

Schéma conceptuel :

```text
SRC/MPCD classique
-> projection de vitesse Q6
-> calcul du défaut de densité N - gamma
-> filtrage low-k du défaut
-> relaxation du flux de masse sur les modes low-k
-> interpolation de la correction vers les particules
-> reprojection finale partielle éventuelle
-> thermostat optionnel
```

La reprojection finale partielle est un compromis : elle réduit la divergence finale réintroduite par la correction de flux de masse, sans annuler entièrement le bénéfice low-k.

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

Pour les cas avec structures, la validation ne se limite donc pas à `std(N)` ou au low-k densitaire. Elle doit aussi inclure :

```text
amplitude de structure cohérente
cohérence modale
enstrophie
fraction haute fréquence de vitesse
vorticité RMS
longueur de recirculation
stabilité temporelle
```

---

## 4. Résultats de validation Q9

### 4.1 Poiseuille 30000 steps — validé
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

### 4.2 Piston 30000 steps — validé

Le cas piston teste le comportement quasi-incompressible sous compression lente. Il ne constitue pas encore une vraie fermeture liquide thermodynamique avec loi d’état virielle, mais il vérifie que Q9 limite les modes compressifs de densité pendant une compression géométrique.

Configuration de référence :

```text
Nx = 32
Ny = 16
gamma = 20
seed = 11
initialPopulationMode = exact_per_cell
nSteps = 30000
sampleEvery = 50
compression finale = 10 %
piston y0/ymin = 1 / 0.9
piston vy = -0.00333333333333
```

Paramètres Q9 :

```matlab
params.projectionStrength = 1.0;
params.massFluxProjectionMode = 'relax_to_uniform_lowk';
params.massFluxProjectionStrength = 1.0;
params.massFluxDensityRelaxationBeta = 0.002;
params.massFluxApplyAfterVelocityProjection = true;
params.massFluxTargetFilter = 'lowpass_fft';
params.massFluxLowKMaxIndex = 2;

params.massFluxFinalVelocityProjectionCleanup = true;
params.massFluxFinalVelocityProjectionStrength = 0.5;
```

Résultat long de référence :

| Métrique | Q9 | Q6 projection seule | Ratio Q9/Q6 |
|---|---:|---:|---:|
| mean std(N) | 3.58696 | 3.61962 | 0.9910 |
| mean out-band | 0.20791 | 0.21202 | 0.9806 |
| mean low-k energy | 5.24213e-5 | 3.21293e-4 | 0.1632 |
| time-avg rel RMS | 0.03296 | 0.03594 | 0.9171 |
| mean rho transport | 9.48122e-3 | 9.42365e-3 | 1.0061 |
| mean div(u) after particles | 2.5614e-2 | 1.0348e-14 | — |
| mean Pkin | 5678.3982 | 5678.4084 | ≈ 1 |
| final Pkin | 5988.9230 | 5988.9147 | ≈ 1 |

Conclusion :

```text
Q9 piston long validé :
- densité locale légèrement meilleure que Q6 ;
- low-k density energy réduite d’un facteur ≈ 6.1 ;
- transport de densité quasi inchangé ;
- pression cinétique conservée ;
- divergence finale résiduelle modérée, conforme au compromis cleanupStrength=0.5.
```

Cette validation soutient l’interprétation de Q9 comme méthode quasi-incompressible effective.

---

### 4.3 Taylor-Green vortex 500 steps — validé

Le cas Taylor-Green est un test propre de préservation des structures vorticitaires lisses. Il ne contient ni murs, ni obstacle, ni contrôle de débit, ni géométrie solide. Il sert donc à isoler l’effet de Q9 sur un mode cohérent.

Champ initial :

```math
u_x = U_0 \sin(k_x x)\cos(k_y y)
```

```math
u_y = -U_0 \cos(k_x x)\sin(k_y y)
```

Configuration validée :

```text
Nx = 32
Ny = 32
gamma = 20
seed = 11
initialPopulationMode = exact_per_cell
nSteps = 500
sampleEvery = 25
dt = 0.002
kBT = 0.01
alphaDeg = 90
Taylor-Green amplitude = 0.3
mode = (1,1)
```

Paramètres Q9 :

```matlab
params.massFluxProjectionMode = 'relax_to_uniform_lowk';
params.massFluxProjectionStrength = 1.0;
params.massFluxDensityRelaxationBeta = 0.002;
params.massFluxApplyAfterVelocityProjection = true;
params.massFluxFinalVelocityProjectionCleanup = true;
params.massFluxFinalVelocityProjectionStrength = 0.5;
params.massFluxTargetFilter = 'lowpass_fft';
params.massFluxLowKMaxIndex = 2;
```

Résultat de référence :

| Métrique | Q9 | Q6 projection seule | Ratio Q9/Q6 |
|---|---:|---:|---:|
| mean std(N) | 3.22595 | 3.24774 | 0.9933 |
| mean out-band | 0.17434 | 0.17322 | 1.0064 |
| mean low-k density energy | 6.63111e-6 | 9.38327e-6 | 0.7067 |
| time-avg rel RMS | 0.09960 | 0.10110 | 0.9851 |
| mean div(u) after particles | 3.20513e-3 | 1.19491e-15 | — |
| final TG amplitude | 0.014872 | 0.014649 | 1.0153 |
| final TG mode energy | 5.52974e-5 | 5.36486e-5 | 1.0307 |
| final TG coherence | 0.20383 | 0.19641 | 1.0378 |
| final enstrophy | 1.18841 | 1.21564 | 0.9776 |
| final high-k velocity fraction | 0.59954 | 0.61231 | 0.9791 |

Conclusion :

```text
Taylor-Green Q9 court validé :
Q9 réduit les modes low-k de densité sans détruire un vortex cohérent lisse.
```

Les runs plus longs `3000` ou `10000` steps avec ce type de champ ne doivent pas être interprétés à partir des valeurs finales, car le mode Taylor-Green finit par tomber au niveau du bruit high-k. Le protocole validant est donc le run court, où le mode cohérent reste mesurable.

---

### 4.4 Cylindre / von Kármán — non stabilisé à ce stade

Un premier banc cylindre a été ajouté pour tester un sillage de type von Kármán :

```text
canal périodique en x
murs en y
obstacle circulaire fixe
conditions cylindre bounceback ou specular
diagnostics de vorticité et sondes de sillage
```

Le premier run court Q6/Q9 a montré que Q9 ne détruisait pas immédiatement la vorticité :

| Métrique | Ratio Q9/Q6 |
|---|---:|
| mean low-k density energy | 0.522 |
| mean std(N) | 1.003 |
| mean out-band | 0.994 |
| mean enstrophy | 1.022 |
| wake omega RMS | 1.028 |
| wake Uy RMS | 1.039 |

Ce résultat court est encourageant comme test de non-destruction instantanée de structures.

En revanche, les runs longs destinés à mesurer un Strouhal n’ont pas pu être stabilisés. Les essais avec forçage volumique, contrôle de vitesse moyenne, paroi `bounceback`, paroi `specular`, réduction de vitesse, réduction du rayon et augmentation de résolution ont montré des décrochements tardifs mais brutaux :

```text
croissance rapide de l’enstrophie
dérive de Ux final
augmentation forte de div(u) reconstruit
échec de l’estimation Strouhal
```

Diagnostic courant :

```text
Le problème n’est pas identifié comme un échec direct de Q9.
Il semble plutôt lié à la combinaison :
- obstacle interne circulaire sous-résolu ;
- traitement solide/fluide encore rudimentaire ;
- projection non strictement masquée sur domaine fluide ;
- domaine périodique fermé ;
- contrôle de vitesse moyenne ;
- injection locale de vorticité près de l’obstacle.
```

Statut :

```text
Le cas cylindre/von Kármán est mis en pause.
Il nécessitera probablement une projection fluide/solide masquée,
un traitement obstacle plus robuste,
ou un domaine plus adapté avant de pouvoir mesurer Re/St.
```

---

### 4.5 Marche / step-channel 10000 steps — validé

Le cas marche est introduit comme test intermédiaire plus robuste que le cylindre :

```text
canal périodique en x
murs en y
bloc solide rectangulaire attaché au mur bas
géométrie alignée grille
écoulement moyen contrôlé
comparaison Q6 vs Q9
```

Objectif :

```text
Tester si Q9 réduit les modes low-k de densité
sans effacer la couche de cisaillement,
la recirculation et la vorticité produites par la marche.
```

Ce cas est moins exigeant que von Kármán car il évite la paroi courbe sous-résolue. Il garde néanmoins des structures hydrodynamiques pertinentes :

```text
séparation
recirculation
couche de cisaillement
vorticité localisée
zone de wake derrière une discontinuité géométrique
```

Configuration validée :

```text
Nx = 48
Ny = 24
gamma = 20
seed = 11
initialPopulationMode = exact_per_fluid_cell
nSteps = 10000
sampleEvery = 50
initialMeanVelocityX = 0.04
meanFlowControlMode = relax_to_target
targetMeanVelocityX = 0.04
meanFlowRelaxationTau = 0.5
dt = 0.001
kBT = 0.02
alphaDeg = 90
step x0/x1/height = 0.3 / 0.7 / 0.25
stepWallMode = bounceback
```

Paramètres Q9 validés pour ce cas :

```matlab
params.projectionStrength = 1.0;
params.massFluxProjectionMode = 'relax_to_uniform_lowk';
params.massFluxProjectionStrength = 1.0;
params.massFluxDensityRelaxationBeta = 5e-4;
params.massFluxApplyAfterVelocityProjection = true;
params.massFluxTargetFilter = 'lowpass_fft';
params.massFluxLowKMaxIndex = 2;
params.massFluxFinalVelocityProjectionCleanup = true;
params.massFluxFinalVelocityProjectionStrength = 0.5;
```

Résultat long de référence :

| Métrique | Q9 | Q6 projection seule | Ratio Q9/Q6 |
|---|---:|---:|---:|
| mean std(N) | 3.08693 | 3.11508 | 0.9910 |
| mean out-band | 0.14474 | 0.14816 | 0.9769 |
| mean low-k density energy | 1.19140e-5 | 1.81760e-5 | 0.6555 |
| time-avg rel RMS | 0.06258 | 0.06177 | 1.0132 |
| mean rho transport | 2.78117e-3 | 2.77753e-3 | 1.0013 |
| mean div(u) after particles | 1.73360e-2 | 2.06486e-2 | 0.8396 |
| mean enstrophy | 0.20414 | 0.19811 | 1.0304 |
| mean shear omega RMS | 0.63974 | 0.64253 | 0.9957 |
| mean recirculation length | 0.98024 | 0.97900 | 1.0013 |
| mean recirculation area | 0.12968 | 0.12370 | 1.0483 |
| probe omega RMS | 0.57106 | 0.59820 | 0.9546 |
| probe Uy RMS | 0.01761 | 0.01732 | 1.0166 |

Conclusion :

```text
Step-channel Q9 long validé :
- stabilité jusqu’à 10000 steps ;
- low-k density energy réduite d’environ 34.5 % ;
- std(N), out-band et div(u) meilleurs que Q6 ;
- transport de densité quasi inchangé ;
- couche de cisaillement et recirculation conservées ;
- aucune destruction visible des structures séparées.
```

Cette validation complète Taylor-Green : Q9 préserve non seulement un vortex lisse, mais aussi une structure de séparation/recirculation générée par une géométrie alignée grille.

---

## 5. Fichiers principaux

### Scripts de référence

| Fichier | Rôle |
|---|---|
| `run_q9_reference_lowk_mass_flux_30000.m` | Script de référence Q9 validé, Poiseuille 30000 steps |
| `run_q9_reference_piston.m` | Script de référence piston Q9, runs courts ou longs |
| `run_q9_reference_taylor_green_short.m` | Script de référence Taylor-Green Q6/Q9 |
| `run_q9_reference_step_short.m` | Script de référence marche / step-channel validé ; peut être renommé en version long |
| `run_q9_reference_vonkarman_short.m` | Script cylindre/von Kármán court ; non stabilisé en long |
| `run_q9_vonkarman_re_st_sweep.m` | Sweep exploratoire Q9-only pour Re/St ; non validé à ce stade |
| `run_compare_density_projection_lowk_mass_flux_poiseuille.m` | Comparaison classic / Q6 projection / Q9 low-k mass-flux |
| `run_compare_projection_piston_q6_q9.m` | Comparaison piston Q6/Q9 |
| `run_compare_projection_taylor_green_q6_q9.m` | Comparaison Taylor-Green Q6/Q9 |
| `run_compare_projection_step_q6_q9.m` | Comparaison marche Q6/Q9 |
| `run_compare_projection_cylinder_q6_q9.m` | Comparaison cylindre Q6/Q9 |
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
| `mpcd_step_classic_piston.m` | Étape classique pour piston mobile |
| `mpcd_step_projection_piston.m` | Étape piston avec projection Q6/Q9 |
| `mpcd_step_projection_periodic_q9.m` | Étape périodique pour Taylor-Green |
| `mpcd_step_classic_cylinder.m` | Étape classique avec obstacle circulaire |
| `mpcd_step_projection_cylinder.m` | Étape cylindre avec projection Q6/Q9 |
| `mpcd_step_classic_step_channel.m` | Étape classique pour marche / step-channel |
| `mpcd_step_projection_step_channel.m` | Étape marche avec projection Q6/Q9 |
| `mpcd_apply_q9_projection_channel.m` | Brique Q9 pour domaines canal / masques fluides |
| `mpcd_apply_q9_projection_periodic.m` | Brique Q9 pour domaine périodique |
| `mpcd_apply_wall_bc_y.m` | Conditions limites en `y` |
| `mpcd_cylinder_mask.m` | Masque solide cylindre |
| `mpcd_step_mask.m` | Masque solide marche |
| `projection_initialize_particles.m` | Initialisation particulaire générique |
| `projection_initialize_particles_cylinder.m` | Initialisation hors cylindre |
| `projection_initialize_particles_step_channel.m` | Initialisation fluide pour marche |
| `projection_initialize_particles_taylor_green.m` | Initialisation Taylor-Green |
| `mpcd_apply_wall_bc_y.m` | Conditions limites en `y` |
| `projection_initialize_particles.m` | Initialisation particulaire, dont `exact_per_cell` |

---

### Projections

| Fichier | Rôle |
|---|---|
| `projection_project_grid_periodic_x_neumann_y.m` | Projection de vitesse sur domaine périodique en `x`, borné en `y` |
| `projection_project_grid_periodic_fft.m` | Projection périodique FFT |
| `projection_project_mass_flux_periodic_x_neumann_y.m` | Projection du flux de masse `N u`, utilisée par Q9 en canal |
| `projection_project_mass_flux_periodic_fft.m` | Projection du flux de masse en domaine périodique |
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
| `projection_vorticity_diagnostics.m` | Vorticité, enstrophie et diagnostics de sillage |
| `projection_wake_shedding_diagnostics.m` | Sondes temporelles et estimation de fréquence de shedding |
| `projection_taylor_green_diagnostics.m` | Amplitude, énergie modale et cohérence Taylor-Green |
| `projection_step_diagnostics.m` | Recirculation, vorticité et couche de cisaillement pour marche |
| `plot_density_homogeneity_comparison.m` | Visualisation des diagnostics de densité |

---

### Ancienne chaîne incompressible / liquid closure

Ces fichiers sont conservés comme base historique et pour de futures fermetures liquides plus physiques :
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

---

## 6. Reproduire les runs de référence

### Poiseuille Q9 30000 steps
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

Dossier typique :
Le script crée un dossier de sortie de la forme :

```text
q9_reference_lowk_mass_flux_30000_YYYYMMDD_HHMMSS/
```

---

### Piston Q9

```matlab
clear
clc
close all

run_q9_reference_piston
```

Pour le jalon long, utiliser :

```matlab
params.nSteps = 30000;
params.sampleEvery = 50;
params.massFluxFinalVelocityProjectionCleanup = true;
params.massFluxFinalVelocityProjectionStrength = 0.5;
```

---

### Taylor-Green Q9

```matlab
clear
clc
close all

run_q9_reference_taylor_green_short
```

Paramètres de référence :

```matlab
params.nSteps = 500;
params.sampleEvery = 25;
params.taylorGreenAmplitude = 0.30;
params.kBT = 0.01;
params.dt = 0.002;
```

---

### Marche / step-channel

```matlab
clear
clc
close all

run_q9_reference_step_short
```

Ce cas est en cours de validation.

---

### Cylindre / von Kármán

```matlab
clear
clc
close all

run_q9_reference_vonkarman_short
```

ou, pour les essais Q9-only :

```matlab
run_q9_vonkarman_re_st_sweep
```

Ce cas n’est pas encore validé en long. Les résultats Strouhal ne doivent pas être interprétés tant que la stabilité temporelle du banc n’est pas obtenue.

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

Pour les géométries solides, utiliser de préférence :

```matlab
params.initialPopulationMode = 'exact_per_fluid_cell';
```

afin de ne peupler que les cellules fluides.
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

### Reprojection finale partielle

Pour les cas structurés actuels :

```matlab
params.massFluxFinalVelocityProjectionCleanup = true;
params.massFluxFinalVelocityProjectionStrength = 0.5;
```

Interprétation :

```text
cleanup = false :
  meilleur contrôle low-k, mais div(u) finale plus élevée.

cleanup = true, strength = 1.0 :
  div(u) finale très faible, mais perte partielle du bénéfice low-k.

cleanup = true, strength = 0.5 :
  compromis actuel retenu.
```

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

Pour les futurs cas piston thermodynamiques, il faudra donc distinguer :
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

### Les valeurs finales doivent parfois être interprétées avec prudence

Dans les cas de vortex libres, comme Taylor-Green, le mode cohérent peut décroître jusqu’au plancher de bruit. Les métriques finales de structure ne sont alors plus significatives.

Pour ces cas, privilégier :

```text
runs courts où le mode reste mesurable
moyennes sur fenêtre cohérente
cohérence modale
énergie du mode cible
fraction high-k
```

---

### Le cas von Kármán n’est pas encore un benchmark validant

Les essais cylindre actuels sont utiles pour le développement, mais ne doivent pas encore être utilisés pour conclure sur un Reynolds ou un Strouhal effectif.

Avant une validation Re/St, il faudra probablement :

```text
améliorer le traitement fluide/solide ;
implémenter une projection réellement masquée ;
stabiliser le contrôle de vitesse moyenne ;
allonger le domaine aval ;
calibrer nu_eff dans les mêmes conditions que le cylindre.
```

---

## 10. Roadmap

### Étape 1 — Stabiliser Q9 comme référence

Statut : fait pour Poiseuille.

- conserver `run_q9_reference_lowk_mass_flux_30000.m` comme script de référence ;
- documenter les résultats Poiseuille 30000 steps ;
- utiliser `seed=11` et `exact_per_cell` pour les comparaisons reproductibles.

---

### Étape 2 — Refactoriser Q9 en brique générique

Statut : partiellement fait.

Fonctions cibles :

```matlab
[stateOut, diag] = mpcd_apply_q9_projection_channel(stateClassic, params)
[stateOut, diag] = mpcd_apply_q9_projection_periodic(stateClassic, params)
```

Objectif :

```text
utiliser la même correction Q9 sur Poiseuille, piston,
Taylor-Green, marche et futurs cas structurés.
```

---

### Étape 3 — Cas piston

Statut : validé en version quasi-incompressible cinématique.

Conclusion :

```text
Q9 piston long réduit fortement le low-k sous compression
sans dérive de pression cinétique.
```

Travail futur :

```text
réintroduire éventuellement une fermeture EOS/virielle
si l’objectif devient une réponse liquide thermodynamique.
```

---

### Étape 4 — Taylor-Green

Statut : validé en run court.

Conclusion :

```text
Q9 préserve une structure vorticitiaire lisse tout en réduisant le low-k densitaire.
```

Travail futur :

```text
ajouter un sweep de graines ;
ajouter des métriques sur fenêtre cohérente ;
éventuellement tester plusieurs modes (1,1), (2,1), (2,2).
```

---

### Étape 5 — Marche / step-channel

Statut : validé en run long 10000 steps.

Conclusion :

```text
Q9 préserve la séparation, la recirculation et la couche de cisaillement,
tout en réduisant le low-k densitaire et la divergence reconstruite.
```

Travail futur :

```text
renommer les scripts/sorties `short` en référence générique ou `long` ;
ajouter éventuellement un sweep de graines ;
conserver ce cas comme benchmark de non-destruction des structures séparées.
```

---

### Étape 6 — Cylindre / von Kármán

Statut : en pause.

Objectif final :

```text
vérifier le comportement sur sillage instationnaire
et, à terme, mesurer un Strouhal.
```

Condition de reprise :

```text
stabiliser le traitement obstacle/projection/contrôle
ou disposer d’un domaine mieux adapté.
```

---
### Étape 7 — Passage vers C++ / parallélisation

Statut : recommandé après gel des validations MATLAB actuelles.

Les validations MATLAB disponibles couvrent maintenant quatre niveaux complémentaires :

```text
Poiseuille    : profil moyen et viscosité effective.
Piston        : quasi-incompressibilité sous compression.
Taylor-Green  : préservation d’un vortex lisse.
Marche        : préservation d’une séparation/recirculation alignée grille.
```

Avant de passer en C++, une petite passe intermédiaire est recommandée :

```text
1. figer les scripts de référence et les paramètres ;
2. renommer les scripts/sorties `short` devenus des références longues ;
3. ajouter un petit tableau de benchmarks dans le README ;
4. éventuellement lancer un sweep de quelques graines sur Taylor-Green et marche ;
5. isoler clairement les briques Q9 canal/périodique pour faciliter le portage.
```

Le passage C++ devient pertinent pour traiter des cas plus résolus, en particulier :

```text
marche plus longue ou plus résolue ;
obstacle carré aligné grille ;
cylindre avec meilleure résolution ;
éventuellement projection fluide/solide masquée.
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
massFluxFinalVelocityProjectionCleanup
massFluxFinalVelocityProjectionStrength
```

Pour les cas à structure, ajouter autant que possible :

```text
enstrophy
vorticity RMS
mode amplitude / coherence
recirculation length
recirculation area
high-k velocity fraction
actualLastStep
stoppedEarly
```

---

## 12. Résumé court

Q9 est la référence actuelle du dépôt.

Elle correspond à :

```text
SRC/MPCD classique
+ projection de vitesse Q6
+ relaxation low-k du flux de masse
+ reprojection finale partielle optionnelle
```

Elle est validée sur :

```text
Poiseuille :
  profil moyen stable,
  low-k density energy réduite d’environ 4.3x.

Piston :
  compression quasi-incompressible stable,
  low-k density energy réduite d’environ 6.1x,
  densité locale légèrement meilleure que Q6.

Taylor-Green :
  structure vorticitiaire lisse préservée,
  low-k density energy réduite d’environ 29 %,
  amplitude/cohérence du mode non dégradées.

Marche / step-channel :
  séparation, recirculation et couche de cisaillement conservées,
  low-k density energy réduite d’environ 34.5 %,
  divergence reconstruite améliorée,
  run long 10000 steps stable.
```

Le cas cylindre/von Kármán n’a pas encore pu être stabilisé en run long. Il est mis en pause, car il teste simultanément la projection Q9, le traitement solide/fluide, le contrôle de débit et une géométrie courbe sous-résolue.
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
