# SRC/MPCD quasi-incompressible — prototype MATLAB `weighted-resampling-core`

Cette branche MATLAB est un prototype de formulation SRC/MPCD quasi-incompressible par **particules pondérées**, **recyclage local de population**, **masque wet/dry**, **particules latentes**, **projection Q6** et **particules virtuelles de paroi wallVP**.

Elle a été développée pour traiter un problème non résolu par les sécurités Q6/Q9/viriel/guards précédentes : la robustesse locale de population lorsque des cellules deviennent pauvres, vides ou fortement surpeuplées. Le principe retenu est de ne plus corriger les cellules pauvres par une force eulérienne attractive, mais de corriger directement la **mesure particulaire locale**.

## État validé de la branche

Jalon actuel : **prototype MATLAB resampling pondéré v1**.

Fonctionnalités présentes :

- masses variables par particule : `state.m` ;
- rôle particulaire : fluide, latent, pool libre ;
- pool préalloué pour éviter les réallocations dynamiques ;
- dépôt pondéré masse/moment/vitesse ;
- collision SRC pondérée ;
- projection Q6 pondérée ;
- insertion de particules dans les cellules sous-peuplées ;
- extraction/recyclage depuis les cellules surpeuplées ;
- remap local conservatif masse + moment ;
- sécurité masses bornées avec correction de vitesse ;
- renormalisation optionnelle préservant localement masse, vitesse moyenne et énergie thermique ;
- masque de cellules mouillées/non mouillées ;
- particules latentes pour domaines non mouillés ou futurs solides mobiles ;
- wallVP pondéré pour valider Poiseuille avec parois ;
- validateurs Taylor--Green et Poiseuille.

Résultat principal : la méthode conserve un comportement fluide crédible sur Taylor--Green forcé et Poiseuille wallVP, tout en contrôlant la densité cellule par cellule.

## Pourquoi cette branche existe

Les développements antérieurs Q6/Q9/viriel amélioraient fortement le comportement quasi-incompressible, mais restaient fragiles face aux cellules pauvres :

- `Q6` impose essentiellement \(\nabla\cdot u=0\) ;
- `Q9` corrige les modes bas de masse-flux \(\nabla\cdot(\rho u)\) ;
- le viriel / pressure guard limite les accumulations ;
- les guards attractifs dense-vers-pauvre testés ont créé des artefacts : cellules vides agissant comme obstacles mous ou attracteurs eulériens.

Conclusion méthodologique : une cellule pauvre ne doit pas créer une pression négative attractive. Le trou de population doit être traité comme un défaut de **support particulaire**, pas comme une force physique.

## Formulation résumée

Chaque particule fluide porte :

\[
(x_p,v_p,m_p),
\]

et les grandeurs cellulaires sont déposées sous forme pondérée :

\[
M_c = \sum_{p\in c} m_p,
\qquad
P_c = \sum_{p\in c} m_p v_p,
\qquad
U_c = \frac{P_c}{M_c}.
\]

La collision SRC pondérée utilise le centre de masse local :

\[
v'_p = U_c + R_\alpha (v_p-U_c).
\]

La boucle de population est :

```text
extraction cellules surpeuplées
-> pool libre
-> insertion cellules sous-peuplées
-> remap local masse/moment
-> thermostat pondéré optionnel
-> renormalisation ponctuelle/conditionnelle des masses si nécessaire
```

Les seuils usuels sont :

```matlab
NMin    = 14;
NTarget = 20;
NMax    = 26;
```

pour `gamma=20`.

## Rôles particulaires

Le champ `state.particleRole` encode :

```text
0 : inactive / free pool
1 : fluid particle
2 : latent particle
```

Les particules latentes sont stockées dans les tableaux mais exclues des dépôts, collisions, projections, remaps, thermostats et diagnostics fluides. Cette couche est la base pour les domaines mouillés évolutifs, les injections, les futurs solides mobiles et les surfaces libres.

## Wet/dry cells

Le masque `state.cellWetMask` distingue :

```text
wet : cellule fluide active, soumise au remap et au recyclage
dry : cellule non mouillée, non forcée vers Mtarget
```

Important : le masque wet/dry seul ne rend pas automatiquement les particules de la zone dry latentes. Pour représenter un vrai domaine vide ou inactif, utiliser également la couche `particleRole`.

## WallVP pondéré

Poiseuille ne doit pas être jugé avec des parois minimales. La validation exploitable utilise des particules virtuelles de paroi pondérées :

```matlab
wallVirtualParticlesEnable = true;
wallVirtualParticlesGeometryMode = 'shifted_solid_fraction';
wallVirtualParticlesDensityFactor = 1.0;
wallVirtualParticlesThermal = true;
```

Les VP sont des contributions collisionnelles agrégées. Elles ne sont pas transportées, ne participent pas au pool et ne modifient pas directement la masse fluide réelle.

## Principaux scripts

### Infrastructure et smokes de base

```text
run_resamp_weighted_mass_smoke.m
run_resamp_local_remap_smoke.m
run_resamp_pool_insertion_smoke.m
run_resamp_pool_insertion_visual_demo.m
```

### Recyclage et tests conçus

```text
run_resamp_recycling_pockets_demo.m
run_resamp_population_heterogeneity_demo.m
run_resamp_population_heterogeneity_sweep.m
```

### Wet/dry, latent, injection

```text
run_resamp_injection_fill_demo.m
```

### Validation fluide

```text
run_resamp_tg_physics_validation.m
run_resamp_poiseuille_physics_validation.m
run_resamp_poiseuille_wallvp_validation.m
run_resamp_poiseuille_wallvp_mass_regularization_once.m
run_resamp_fluid_physics_validation_suite.m
```

### Analyse et visualisation Poiseuille

```text
resamp_poiseuille_visualize_frame.m
resamp_poiseuille_fit_profile.m
analyze_resamp_poiseuille_viscosity.m
```

## Ordre conseillé pour retrouver le développement

### 1. Smoke masse pondérée

```matlab
out = run_resamp_weighted_mass_smoke( ...
    'method','weighted_q6', ...
    'steps',500, ...
    'summaryEvery',50);
```

Objectif : vérifier que `state.m`, le dépôt pondéré, la collision pondérée et Q6 fonctionnent en masse unitaire.

### 2. Remap local masse/moment

```matlab
out = run_resamp_local_remap_smoke( ...
    'method','weighted_q6', ...
    'remapMethod','scale_preserve_velocity', ...
    'steps',200, ...
    'summaryEvery',20);
```

Critère attendu : `MRelRms` et résidus de remap au roundoff.

### 3. Test conçu poches vide/surpeuplée

```matlab
outP = run_resamp_recycling_pockets_demo( ...
    'steps',300, ...
    'visualEvery',10, ...
    'summaryEvery',25, ...
    'rngSeed',12345);
```

Critère attendu :

```text
extractedParticlesCumulative = insertedParticlesCumulative = 500
N[min,max] revient près de [19,22]
MRelRms ~ 1e-16
```

### 4. Hétérogénéité distribuée de population

```matlab
outSweep = run_resamp_population_heterogeneity_sweep( ...
    'PopulationHeterogeneityStdList',[0 4 8 16], ...
    'steps',300, ...
    'summaryEvery',25, ...
    'NMin',14, ...
    'NTarget',20, ...
    'NMax',26, ...
    'RemapMassSafetyEnable',true, ...
    'rngSeeds',12345);
```

Critère attendu : population ramenée dans `[14,26]`, masse cellulaire au roundoff, masses bornées par la sécurité.

### 5. Remplissage/injection latent/wet

```matlab
outFill = run_resamp_injection_fill_demo( ...
    'method','weighted_classic', ...
    'steps',200, ...
    'visualEvery',5, ...
    'summaryEvery',5, ...
    'Nx',64, ...
    'Ny',32, ...
    'gamma',20, ...
    'dt',0.05, ...
    'inletVelocityX',2.0, ...
    'inletVelocityY',0.0, ...
    'injectionParticlesPerCell',20, ...
    'NMin',10, ...
    'NTarget',20, ...
    'NMax',30, ...
    'ThermostatAfterRemap',true, ...
    'ThermostatStrength',0.25, ...
    'ThermostatTargetKBT',0.01, ...
    'rngSeed',12345);
```

Critère attendu : `nWetCells` augmente progressivement, les cellules dry ne sont pas forcées vers `Mtarget`, et `MRelRms` reste petit sur les cellules wet.

## Validation fluide

### Suite rapide

```matlab
outSmoke = run_resamp_fluid_physics_validation_suite( ...
    'quick',true, ...
    'visualEvery',0);
```

### Taylor--Green forcé

```matlab
outTG = run_resamp_tg_physics_validation( ...
    'steps',3000, ...
    'summaryEvery',50, ...
    'visualEvery',0, ...
    'Nx',64, ...
    'Ny',64, ...
    'gamma',20, ...
    'NMin',14, ...
    'NTarget',20, ...
    'NMax',26, ...
    'TaylorGreenForcingAmplitude',0.12, ...
    'ThermostatStrength',0.25, ...
    'rngSeed',12345);
```

Résultat observé à ce jalon :

```text
q6_resampled:
finalRmsDiv ~ 4e-17
finalMRelRms ~ 1e-16
finalKBT ~ 0.011
TG coherence ~ 0.89
N[min,max] = [14,26]
```

### Poiseuille wallVP : run décisif

```matlab
outW = run_resamp_poiseuille_wallvp_validation( ...
    'steps',5000, ...
    'sampleEvery',50, ...
    'summaryEvery',100, ...
    'visualEvery',200, ...
    'saveFinalFigures',true, ...
    'Nx',64, ...
    'Ny',32, ...
    'gamma',20, ...
    'dt',0.005, ...
    'bodyForceX',0.005, ...
    'NMin',14, ...
    'NTarget',20, ...
    'NMax',26, ...
    'ThermostatStrength',0.25, ...
    'wallVirtualParticlesEnable',true, ...
    'wallVirtualParticlesGeometryMode','shifted_solid_fraction', ...
    'wallVirtualParticlesDensityFactor',1.0, ...
    'wallVirtualParticlesThermal',true, ...
    'poiseuilleFitModel','slip', ...
    'excludeWallCells',2, ...
    'fitStartFraction',0.5, ...
    'rngSeed',12345);
```

Résultat observé :

```text
classic_wallvp_reference:
nuEffSlip ~ 16.99
R2Slip ~ 0.9959
MRelRMS ~ 0.15
N[min,max] ~ [8,32]

q6_resampled_wallvp:
nuEffSlip ~ 17.65
R2Slip ~ 0.9979
MRelRMS ~ 1e-16
N[min,max] = [14,26]
finalMassRelStd ~ 0.234
```

Interprétation : le SRC classic wallVP fournit la référence physique ; Q6/resampling conserve essentiellement la viscosité et le profil tout en contrôlant la densité.

### Lancer uniquement le SRC classic wallVP pur

```matlab
outClassicOnly = run_resamp_poiseuille_wallvp_validation( ...
    'cases',{'src_classic_wallvp_only'}, ...
    'steps',5000, ...
    'sampleEvery',50, ...
    'summaryEvery',100, ...
    'visualEvery',200, ...
    'Nx',64, ...
    'Ny',32, ...
    'gamma',20, ...
    'dt',0.005, ...
    'bodyForceX',0.005, ...
    'wallVirtualParticlesEnable',true, ...
    'wallVirtualParticlesGeometryMode','shifted_solid_fraction', ...
    'wallVirtualParticlesDensityFactor',1.0, ...
    'wallVirtualParticlesThermal',true, ...
    'poiseuilleFitModel','slip', ...
    'excludeWallCells',2, ...
    'fitStartFraction',0.5, ...
    'rngSeed',12345);
```

## Renormalisation des masses

La renormalisation locale préserve :

\[
M_c,\qquad U_c,\qquad E_{\mathrm{th},c}.
\]

Elle rapproche les masses de `Mcell/Ncell`, puis corrige les vitesses par translation et rescaling thermique.

### Smoke court : renormalisation au step 300

```matlab
outRegShort = run_resamp_poiseuille_wallvp_mass_regularization_once( ...
    'steps',500, ...
    'massRegularizationStep',300, ...
    'sampleEvery',25, ...
    'summaryEvery',25, ...
    'visualEvery',0, ...
    'saveFinalFigures',true, ...
    'Nx',64, ...
    'Ny',32, ...
    'gamma',20, ...
    'dt',0.005, ...
    'bodyForceX',0.005, ...
    'NMin',14, ...
    'NTarget',20, ...
    'NMax',26, ...
    'ThermostatStrength',0.25, ...
    'wallVirtualParticlesEnable',true, ...
    'wallVirtualParticlesGeometryMode','shifted_solid_fraction', ...
    'wallVirtualParticlesDensityFactor',1.0, ...
    'wallVirtualParticlesThermal',true, ...
    'poiseuilleFitModel','slip', ...
    'excludeWallCells',2, ...
    'fitStartFraction',0.5, ...
    'rngSeed',12345);
```

Contrôle :

```matlab
T = outRegShort.cases.q6_resampled_wallvp.summary;
T(T.step >= 250 & T.step <= 350, {'step','mParticleRelStd', ...
    'massRegCells','massRegParticles', ...
    'massRegRelStdBefore','massRegRelStdAfter', ...
    'MRelRms','kBTWeighted'})
```

Résultat observé : la renormalisation au step 300 réduit `mParticleRelStd` sans rupture visible de `MRelRms`, `kBT` ou de la dynamique court terme.

### Run long : renormalisation au step 3000

```matlab
outReg = run_resamp_poiseuille_wallvp_mass_regularization_once( ...
    'steps',5000, ...
    'massRegularizationStep',3000, ...
    'sampleEvery',50, ...
    'summaryEvery',100, ...
    'visualEvery',200, ...
    'saveFinalFigures',true, ...
    'Nx',64, ...
    'Ny',32, ...
    'gamma',20, ...
    'dt',0.005, ...
    'bodyForceX',0.005, ...
    'NMin',14, ...
    'NTarget',20, ...
    'NMax',26, ...
    'ThermostatStrength',0.25, ...
    'wallVirtualParticlesEnable',true, ...
    'wallVirtualParticlesGeometryMode','shifted_solid_fraction', ...
    'wallVirtualParticlesDensityFactor',1.0, ...
    'wallVirtualParticlesThermal',true, ...
    'poiseuilleFitModel','slip', ...
    'excludeWallCells',2, ...
    'fitStartFraction',0.5, ...
    'rngSeed',12345);
```

Résultat observé :

```text
q6_resampled_wallvp + renormalisation au step 3000:
finalMassRelStd ~ 0.200
nuEffSlip ~ 18.01
R2Slip ~ 0.9978
MRelRms ~ 1e-16
kBT ~ 0.0099
N[min,max] = [14,26]
```

Interprétation : une renormalisation ponctuelle réduit l’hétérogénéité des masses sans dégrader de manière significative le profil Poiseuille. Elle doit être considérée comme sécurité ponctuelle ou conditionnelle, plutôt qu’application systématique par défaut.

## Sorties produites

Les scripts écrivent dans `runs/`, notamment :

```text
runs/resamp_pool_insertion_smoke/
runs/resamp_pool_insertion_visual/
runs/resamp_population_heterogeneity_sweep/
runs/resamp_tg_physics_validation/
runs/resamp_poiseuille_wallvp_validation/
```

Fichiers utiles :

```text
summary.csv
timeseries.csv
profile_late_mean.csv
*_spatial.png
*_profile.png
```

## Diagnostics à regarder systématiquement

```text
MRelRms
mParticleRelStd
mParticleMin / mParticleMax
NMin / NMax
NpActive / Nfree
insertedParticlesCumulative
extractedParticlesCumulative
kBTWeighted
rmsDivParticleAfter
R2Slip / R2NoSlip
nuEffSlip / nuEffNoSlip
slipRatio / wallSlipRatioObserved
```

## Points de vigilance

1. **Ne pas juger Poiseuille sans wallVP.** Les parois minimales produisent trop de glissement et des fits non exploitables.
2. **La masse cellulaire est contrôlée au prix d’une hétérogénéité des masses particulaires.** La renormalisation M,u,E est la sécurité prévue contre une dérive longue.
3. **Les cellules dry ne sont pas des solides physiques à elles seules.** Pour un vrai vide ou solide, utiliser aussi les rôles particulaires latents/inactifs.
4. **La projection Q6 masquée pour domaines partiellement mouillés n’est pas encore la validation principale.** Les validations actuelles concernent surtout domaine plein, injection démonstrative et canal wallVP.
5. **Le passage OpenMP doit reproduire TG et Poiseuille wallVP avant solides mobiles ou surface libre.**

## Documentation complémentaire

Le dossier `doc/` contient les étapes détaillées du développement :

```text
doc/README_RESAMPLING_CORE.md
doc/README_RESAMPLING_REMAP.md
doc/README_RESAMPLING_POOL_INSERTION.md
doc/README_RESAMPLING_EXTRACTION_RECYCLING.md
doc/README_RESAMPLING_HETEROGENEOUS_POPULATION.md
doc/README_RESAMPLING_MASS_SAFETY.md
doc/README_RESAMPLING_WET_CELLS.md
doc/README_RESAMPLING_LATENT_INJECTION.md
doc/README_RESAMPLING_WALLVP_POISEUILLE.md
doc/README_RESAMPLING_POISEUILLE_DECISIVE_VISUAL.md
doc/README_RESAMPLING_MASS_REGULARIZATION.md
```

Le rapport complet du développement est disponible sous forme LaTeX/PDF dans les fichiers produits pendant la clôture de la branche.

## Suite recommandée

La suite logique du projet est le portage vers `SRC_incompressible_openMP` sur une branche dédiée, par exemple :

```text
feature/weighted-resampling-openmp
```

Objectif initial du portage : reproduire en OpenMP les deux validateurs MATLAB :

```text
Taylor--Green forcé périodique
Poiseuille wallVP q6_resampled
```

avant d’ajouter :

```text
solides immergés fixes
solides mobiles
inlet/outlet
surface libre / tension de surface
```
