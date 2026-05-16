# SRC_MPCD_projection_matlab — Q9 forced Taylor–Green validation

Prototype MATLAB minimal pour la validation d’une méthode SRC/MPCD quasi-incompressible par projection de vitesse et projection low-k du flux de masse.

Cette branche est volontairement centrée sur un cas test unique et propre : le vortex de Taylor–Green forcé en domaine entièrement périodique. Elle sert à figer une base Q6/Q9 bulk-only corrigée, avec thermostat commun, correction globale exacte de quantité de mouvement, visualisation et fit de viscosité effective.

```text
Branche recommandée : q9-mass-flux-projection-momentum-correction-parallel
Méthode principale  : Q9 = Q6 + projection low-k du flux de masse N u
Cas validant        : Taylor–Green forcé périodique
Statut              : validation bulk positive sur structure tourbillonnaire forcée
Hors périmètre      : piston, marche, cylindre, surface libre, Q10/interface
```

---

## 1. Objectif scientifique

Le SRC/MPCD classique possède des fluctuations mésoscopiques utiles, mais il laisse apparaître des modes compressifs de densité importants dans des configurations où l’on souhaite approcher un comportement liquide ou quasi-incompressible.

L’objectif de cette branche est de valider une méthode de projection qui réduit ces modes compressifs sans détruire une structure tourbillonnaire organisée.

La projection incompressible de base vise :

```math
\nabla \cdot u \simeq 0.
```

La correction Q9 agit en plus sur le flux de masse discret :

```math
M = N u,
```

avec `N` l’occupation particulaire cellulaire et `u` la vitesse moyenne cellulaire. L’objectif est de réduire les grandes longueurs d’onde de densité en corrigeant préférentiellement les modes low-k de :

```math
\nabla \cdot (N u).
```

---

## 2. Méthodes comparées

Le script principal compare trois méthodes.

### CLASSIC

SRC/MPCD périodique classique avec forçage Taylor–Green et thermostat post-step commun.

```text
forcing Taylor–Green
-> streaming périodique
-> collision SRC/MPCD
-> thermostat cellulaire
-> diagnostics
```

### Q6

SRC/MPCD classique suivi d’une projection de vitesse :

```text
classic step
-> projection de vitesse div(u) ≈ 0
-> correction globale exacte de quantité de mouvement
-> thermostat cellulaire
-> diagnostics
```

### Q9

Q9 ajoute à Q6 une relaxation low-k du flux de masse :

```text
classic step
-> projection de vitesse Q6
-> correction globale exacte de quantité de mouvement Q6
-> correction low-k du flux de masse N u
-> correction globale exacte de quantité de mouvement Q9
-> cleanup final éventuel de div(u)
-> correction globale exacte de quantité de mouvement cleanup
-> thermostat cellulaire
-> diagnostics
```

Paramètres Q9 de référence :

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

---

## 3. Correction globale exacte de quantité de mouvement

Une projection de pression modifie les vitesses particulaires :

```math
v_i^{new} = v_i + \delta v_i.
```

Elle peut donc injecter une quantité de mouvement globale numérique :

```math
\Delta P_{proj} = \sum_i m_i \delta v_i.
```

La branche intègre une correction globale exacte appliquée aux corrections numériques de projection :

```math
\delta v_i^{corr}
= \delta v_i
- \frac{\sum_j m_j \delta v_j}{\sum_j m_j}.
```

Pour masses unitaires :

```math
\delta v_i^{corr}
= \delta v_i - \langle \delta v \rangle.
```

Ainsi :

```math
\sum_i \delta v_i^{corr} = 0.
```

Cette correction est appliquée aux sous-étapes numériques Q6, Q9 et cleanup. Elle ne retire pas les impulsions physiques éventuelles associées au forcing, aux collisions ou aux conditions aux limites.

Diagnostics associés :

```text
meanMomentumRawDV       : norme moyenne du kick global brut par particule
meanMomentumResidualDV  : résidu après correction globale exacte
```

Dans le run de référence, le résidu corrigé est de l’ordre de la précision machine.

---

## 4. Thermostat commun et diagnostic thermique corrigé

Le cas Taylor–Green forcé utilise un thermostat cellulaire commun aux trois méthodes. La comparaison CLASSIC/Q6/Q9 serait biaisée si le thermostat n’était appliqué qu’après projection.

Ordre final retenu :

```text
CLASSIC :
  forcing -> streaming -> collision -> thermostat -> diagnostics

Q6/Q9 :
  forcing -> streaming -> collision -> projection(s)
  -> correction globale du moment
  -> thermostat -> diagnostics
```

Le thermostat cellulaire recale les vitesses relatives autour de la vitesse moyenne de cellule :

```math
v_i \leftarrow u_c + s_c (v_i - u_c),
```

ce qui conserve la quantité de mouvement de chaque cellule à l’arrondi près.

Le diagnostic thermique a été corrigé pour utiliser la même convention d’indexation cellulaire que le thermostat. Le champ suivi `kBT cell` correspond maintenant correctement à la température imposée.

---

## 5. Cas test : Taylor–Green forcé périodique

Le cas validant est un vortex de Taylor–Green forcé en domaine périodique 2D.

Champ de vitesse modal :

```math
u_x = A \sin(k_x x) \cos(k_y y),
```

```math
u_y = -A \frac{k_x}{k_y} \cos(k_x x) \sin(k_y y).
```

Forçage divergence-free de même forme :

```math
f_x = F \sin(k_x x) \cos(k_y y),
```

```math
f_y = -F \frac{k_x}{k_y} \cos(k_x x) \sin(k_y y).
```

Le cas est intéressant parce que la structure tourbillonnaire attendue est connue et mesurable. La validation ne repose donc pas seulement sur des diagnostics indirects comme une longueur de recirculation.

Diagnostics principaux :

```text
amplitude du mode Taylor–Green
cohérence modale
énergie du mode
vorticité / enstrophie
fraction high-k de l’énergie de vitesse
low-k density energy
density relative RMS
kBT cell
viscosité effective forcée
correction globale de moment
```

---

## 6. Lancer le cas de validation

Depuis MATLAB, à la racine du dépôt :

```matlab
run_q9_forced_taylor_green_validation
```

Le script compare par défaut :

```text
CLASSIC
Q6
Q9
```

Configuration de référence :

```text
Nx = Ny = 64
gamma = 20
nSteps = 10000
sampleEvery = 100
visualEvery = 250
dt = 0.001
kBT = 0.01
TG initial amplitude = 0.10
TG force amplitude   = 0.08
TG mode              = (1,1)
```

Le run génère un dossier de sortie du type :

```text
q9_forced_taylor_green_medium64_compare_classic_q6_q9_YYYYMMDD_HHMMSS/
```

avec notamment :

```text
forced_tg_summary.txt
forced_tg_summary.csv
forced_tg_classic_timeseries.csv
forced_tg_q6_timeseries.csv
forced_tg_q9_timeseries.csv
forced_tg_timeseries.png
forced_tg_viscosity_fit.png
forced_tg_validation.mat
```

Ces sorties ne doivent pas être versionnées.

---

## 7. Visualisation

Le cas Taylor–Green forcé inclut une visualisation live pour suivre les structures, et non seulement les scalaires de diagnostic.

La visualisation affiche :

```text
particules,
densité relative N/gamma - 1,
vitesse |u| avec quiver,
vorticité omega.
```

Des champs moyens peuvent également être accumulés selon les paramètres du run.

Paramètres utiles :

```matlab
params.visualEnable = true;
params.visualEvery = 250;
params.visualSaveFrames = false;
```

Pour un run plus long ou plus coûteux, augmenter `visualEvery`.

---

## 8. Fit de viscosité effective

Le cas forcé permet d’estimer une viscosité effective à partir de l’amplitude modale.

Le modèle utilisé est :

```math
\frac{dA}{dt} = F - \nu_{eff} (k_x^2 + k_y^2) A.
```

D’où :

```math
\nu_{eff} = \frac{F - dA/dt}{(k_x^2 + k_y^2) A}.
```

Le code fournit plusieurs estimations :

```text
meanNuEffForced
finalNuEffForced
nuEffFitPreferred
nuEffFitDerivativeKnownF
nuEffFitExpKnownF
nuEffFitPlateau
viscosityFitR2
```

Le fit est surtout interprétable lorsque l’amplitude atteint une fenêtre quasi stationnaire ou suit une relaxation exponentielle suffisamment nette. Dans le run de référence, les différentes estimations donnent des viscosités effectives cohérentes entre CLASSIC, Q6 et Q9.

---

## 9. Résultat de référence actuel

Run validant :

```text
methods = CLASSIC, Q6, Q9
Nx = Ny = 64
gamma = 20
nSteps = 10000
dt = 0.001
kBT = 0.01
TG initial/force amplitude = 0.10 / 0.08
Q9 beta/lowK/cleanup = 5e-4 / 2 / 0.5
```

Synthèse :

| Métrique | CLASSIC | Q6 | Q9 |
|---|---:|---:|---:|
| final amplitude | 0.05294 | 0.04916 | 0.05046 |
| final coherence | 0.557 | 0.739 | 0.755 |
| final high-k fraction | 0.358 | 0.203 | 0.188 |
| final low-k density | 5.40e-3 | 5.83e-6 | 5.24e-7 |
| mean density rel RMS | 0.241 | 0.166 | 0.164 |
| final kBT cell | 0.010 | 0.010 | 0.010 |
| nu_eff fit preferred | 0.01916 | 0.02058 | 0.02045 |
| mean raw momentum kick / particle | NaN | 4.51e-5 | 3.83e-5 |
| mean residual momentum kick / particle | NaN | 8.12e-19 | 1.77e-18 |

Lecture :

```text
CLASSIC conserve légèrement plus d’amplitude brute,
mais avec beaucoup plus de bruit high-k et de modes compressifs low-k.

Q6 améliore fortement la cohérence modale et réduit les modes compressifs.

Q9 donne le meilleur contrôle de densité low-k,
la meilleure cohérence finale,
et une viscosité effective très proche de Q6.
```

Conclusion du run : Q9 rend le fluide nettement moins compressible sans supprimer la structure tourbillonnaire imposée. La méthode augmente légèrement la viscosité effective, ce qui doit être interprété comme une propriété effective du fluide projeté.

---

## 10. Fichiers nécessaires dans cette branche minimale

### Script principal

| Fichier | Rôle |
|---|---|
| `run_q9_forced_taylor_green_validation.m` | Lance la comparaison CLASSIC/Q6/Q9 et écrit les résumés |
| `run_projection_forced_taylor_green_demo.m` | Moteur d’un run Taylor–Green forcé pour une méthode donnée |

### Initialisation, forcing, diagnostics et visualisation

| Fichier | Rôle |
|---|---|
| `projection_initialize_particles_taylor_green_forced.m` | Initialise les particules et le champ TG |
| `projection_taylor_green_mode_at_points.m` | Évalue le mode TG aux positions particulaires |
| `projection_apply_taylor_green_forcing.m` | Applique le forçage TG divergence-free |
| `projection_taylor_green_diagnostics.m` | Calcule amplitude, cohérence, enstrophie, high-k, low-k |
| `projection_taylor_green_visualize_frame.m` | Visualisation instantanée des champs |
| `projection_fit_forced_taylor_green_viscosity.m` | Fit de viscosité effective |

### Steps SRC/MPCD et projections

| Fichier | Rôle |
|---|---|
| `mpcd_step_classic_periodic_forced.m` | Step SRC/MPCD périodique forcé avec thermostat commun |
| `mpcd_step_projection_periodic_q9_forced.m` | Step Q6/Q9 périodique forcé |
| `mpcd_apply_q9_projection_periodic.m` | Projection périodique Q6/Q9, cleanup, correction de moment |
| `projection_project_grid_periodic_fft.m` | Projection de vitesse périodique par FFT |
| `projection_project_mass_flux_periodic_fft.m` | Projection low-k du flux de masse périodique |

### Utilitaires communs indispensables

| Fichier | Rôle |
|---|---|
| `projection_deposit_particles_to_grid.m` | Dépôt particules -> grille |
| `projection_interpolate_grid_delta_to_particles.m` | Interpolation correction grille -> particules |
| `projection_population_diagnostics.m` | Diagnostics d’occupation cellulaire |
| `projection_apply_global_momentum_correction.m` | Correction globale exacte de quantité de mouvement |
| `projection_apply_cell_thermostat.m` | Thermostat cellulaire vectorisé |
| `projection_thermal_diagnostics.m` | Diagnostic thermique corrigé |

Les anciens scripts `piston`, `cylinder`, `step`, `staircase`, `surface_leveling`, `inclined_layer`, `dam_break` et `Q10` ne sont pas nécessaires dans cette branche minimale.

---

## 11. Fichiers générés à ignorer

Ajouter ou conserver dans `.gitignore` :

```gitignore
# Generated validation outputs
q9_*/
*.mat
*.asv
forced_tg_*_timeseries.csv
forced_tg_summary.csv
forced_tg_summary.txt
forced_tg_timeseries.png
forced_tg_viscosity_fit.png
```

---

## 12. Limites actuelles

Cette branche ne prétend pas valider tous les écoulements SRC/MPCD possibles.

Limites assumées :

```text
pas de surface libre,
pas de contact liquide/air,
pas de piston,
pas d’obstacle solide courbe,
pas de von Kármán,
pas de validation locale complète du moment q = rho u.
```

La correction de quantité de mouvement actuelle est globale. Elle empêche toute dérive du moment total due à la projection, mais ne garantit pas encore une reconstruction locale conservative de :

```math
q = \rho u.
```

Une étape ultérieure devra ajouter des diagnostics locaux de `delta(rho u)` et éventuellement une reconstruction conservative locale du champ de moment.

---

## 13. Prochaines étapes recommandées

1. Confirmer le run Taylor–Green forcé avec quelques amplitudes de forcing :

```text
tgForceAmplitude = 0.04, 0.08, 0.12, 0.16
```

2. Vérifier la stabilité de `nu_eff` estimé.

3. Ajouter des diagnostics locaux du moment :

```text
Delta q_c = (N u)_after_projection - (N u)_before_projection
spectre de Delta q
norme RMS bulk
corrélation avec grad(phi)
```

4. Seulement ensuite, reprendre des cas plus complexes : dipôle vortex, Kelvin–Helmholtz, marche raffinée, cylindre.

5. Porter la méthode Q9 bulk vers C++/OpenMP lorsque la branche MATLAB minimale est figée.

---

## 14. Conclusion

Cette branche fournit un état propre et minimal de la méthode Q9 :

```text
projection div(u),
projection low-k du flux de masse N u,
thermostat commun corrigé,
correction globale exacte de quantité de mouvement,
visualisation des structures,
fit de viscosité effective.
```

Le cas Taylor–Green forcé montre que Q9 réduit les modes compressifs de densité de plusieurs ordres de grandeur tout en conservant une structure tourbillonnaire organisée. La méthode modifie légèrement la viscosité effective, mais ce comportement est mesurable, stable et compatible avec l’interprétation de Q9 comme fluide SRC/MPCD quasi-incompressible effectif.
