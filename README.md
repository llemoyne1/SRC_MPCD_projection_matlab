# SRC_MPCD_projection_matlab

Prototype MATLAB pour le développement, le diagnostic et la validation d’une méthode SRC/MPCD quasi-incompressible fondée sur :

```text
SRC/MPCD classique
→ projection de vitesse Q6 : div(u) ≈ 0
→ correction Q9 low-k du flux de masse div(Nu)
→ thermostat cellulaire corrigé
→ fermeture virielle optionnelle pour réponse pression–densité
```

Le dépôt sert actuellement à valider un modèle particulaire mésoscopique capable de conserver les avantages du SRC/MPCD — particules, fluctuations, transport mésoscopique — tout en limitant les modes compressifs cohérents de densité et en produisant une réponse pression–densité exploitable sous compression piston.

---

## 1. État courant

### Modèle courant

La chaîne de référence n’est plus une comparaison de variantes `classic`, `Q6`, `Q9`, mais le modèle complet :

```text
q9_full = SRC/MPCD classic + Q6 + Q9 + thermostat corrigé
```

La fermeture virielle est ajoutée comme fermeture EOS au-dessus du modèle `q9_full`, et non comme méthode concurrente.

### Jalons validés

État courant de la branche piston :

```text
- Q6/Q9 complet piston : fonctionnel
- opérateur general_bc : stabilisé par gauge-fix du mode constant
- thermostat cellulaire corrigé : stable
- piston full-model : conserve la compression globale
- EOS virielle bulk : fonctionnelle
- Ptot = Pkin + Pvir : validé
- Pdrive = Ptot : validé comme diagnostic EOS
- kick viriel actif faible : stable
- correction globale du moment du kick : validée
- limiteur viriel : implémenté et exposé, non déclenché dans les runs doux
- staircase EOS : fonctionnel
- wallVP pressure : algébriquement cohérent, mais physiquement non validé
```

### Dernier run staircase EOS de référence

Configuration représentative :

```text
Nx = 32, Ny = 32, gamma = 20
Np = 20480
dt = 0.001
kBT = 0.01
compression levels = [0 0.0025 0.005 0.01 0.02 0.035 0.05]
Kvirial = 0.05
virialBeta = 0.1
virialMaxDuFractionThermal = 0.02
massFluxProjectionOperator = general_bc
massFluxTargetFilter = elliptic_lowpass
```

Résultats principaux :

| Quantité | Keff | R² | Interprétation |
|---|---:|---:|---|
| `PkinMean` | ≈ 219.66 | 1.0 | pression cinétique isotherme cohérente |
| `PvirMean` | ≈ 1042.6 | 1.0 | contribution virielle linéaire en densité |
| `PtotMean` | ≈ 1262.2 | 1.0 | pression totale = pression cinétique + virielle |
| `PdriveMean` | ≈ 1262.2 | 1.0 | diagnostic du champ de drive aligné avec `Ptot` |

À compression finale 5 % :

```text
final rho physical mean = 21557.8947368
final PkinMean          ≈ 226.9456
final PvirMean          ≈ 53.8947
final PtotMean          ≈ 280.8403
final kBTCell           = 0.01
final rho rel. error    = 0
final std(N)            ≈ 3.47
final low-k density     ≈ 1.0e-7
```

---

## 2. Objectif scientifique

L’objectif est de construire un fluide SRC/MPCD quasi-incompressible qui vérifie simultanément :

1. conservation de la masse particulaire et de la compression globale ;
2. homogénéisation des grandes longueurs d’onde de densité ;
3. stabilité thermique ;
4. conservation d’une dynamique particulaire mésoscopique ;
5. réponse pression–densité contrôlable sous piston.

Le test piston n’est donc pas seulement un test de réduction de `std(N)` ou de `lowKDensity`. Il doit vérifier que le code conserve la compression globale tout en augmentant la pression lorsque la densité augmente.

---

## 3. Méthode numérique

### 3.1 Q6 — projection de vitesse

Q6 applique une projection de vitesse sur grille après l’étape SRC/MPCD classique :

```text
SRC/MPCD classic
→ dépôt particules-grille
→ projection de vitesse : div(u) ≈ 0
→ interpolation de la correction vers particules
```

Q6 contrôle la composante incompressible de la vitesse, mais ne suffit pas toujours à contrôler les grandes structures de densité advectées par `Nu`.

### 3.2 Q9 — correction low-k du flux de masse

Q9 ajoute une correction du flux de masse :

```text
M = N u
```

afin de réduire le transport compressif cohérent :

```text
div(Nu)
```

La correction cible les modes basse fréquence de densité. Elle ne cherche pas à lisser brutalement le bruit local d’occupation cellulaire.

Paramètres représentatifs :

```matlab
params.massFluxProjectionMode = 'relax_to_uniform_lowk';
params.massFluxProjectionStrength = 1.0;
params.massFluxDensityRelaxationBeta = 5e-4;
params.massFluxApplyAfterVelocityProjection = true;
params.massFluxProjectionOperator = 'general_bc';
params.massFluxTargetFilter = 'elliptic_lowpass';
params.massFluxLowKMaxIndex = 2;
params.massFluxFinalVelocityProjectionCleanup = false;
```

### 3.3 Opérateur `general_bc` et gauge-fix

L’opérateur elliptique général `projection_project_mass_flux_general_bc.m` est utilisé pour les géométries périodiques/bornées.

Pour les cas fermés de type périodique/Neumann, le potentiel elliptique possède un mode constant nul. Le solveur utilise désormais un gauge-fix du mode 0 plutôt qu’une régularisation artificielle :

```text
mode0GaugeActive = true
effectiveRegularization = 0
```

Le gauge-fix ne change pas la physique : il fixe seulement la constante additive du potentiel, qui ne modifie pas le gradient.

### 3.4 Fermeture virielle EOS

La fermeture virielle ajoute une pression excédentaire :

```matlab
Pvir = Kvirial * (rho - rhoEOSRef);
Ptot = Pkin + Pvir;
Pdrive = Ptot;
```

Le kick actif applique une correction de vitesse proportionnelle au gradient de pression :

```matlab
du = - virialBeta * dt ./ rhoKick .* grad(Pdrive);
```

Le kick est suivi d’une correction globale exacte du moment, puis d’un thermostat cellulaire final.

### 3.5 Limiteur de sécurité du kick

Le kick viriel est borné par :

```matlab
duMax = virialMaxDuFractionThermal * sqrt(kBT);
```

Diagnostics associés :

```text
virialDuRms
virialDuOverThermalRms
virialLimitedCellFraction
virialLimitedCellCount
virialLimiterDuMax
virialResidualMomentumKickNorm
```

Dans le run staircase de référence, le limiteur est exposé mais non déclenché :

```text
maxVirialLimitedCellFraction = 0
virialResidualMomentumKickNorm ~ 1e-15
```

---

## 4. Scripts principaux

### Validation piston / EOS

| Fichier | Rôle |
|---|---|
| `run_q9_piston_staircase_eos_validation.m` | Validation EOS par paliers de compression |
| `run_q9_piston_full_model_validation.m` | Campagne piston full-model |
| `run_projection_piston_wallvp_demo.m` | Driver piston général avec Q6/Q9, viriel, wallVP |
| `analyze_piston_compressibility.m` | Fits `Pkin`, `Pvir`, `Ptot`, `Pdrive` vs densité |
| `mpcd_step_projection_piston.m` | Étape piston avec Q6/Q9 complet |
| `mpcd_step_classic_piston.m` | Étape SRC/MPCD piston + diagnostics wallVP |

### Noyaux Q6/Q9

| Fichier | Rôle |
|---|---|
| `mpcd_apply_q9_projection_channel.m` | Application Q6/Q9 en canal / domaine borné |
| `projection_project_mass_flux_general_bc.m` | Opérateur elliptique général pour flux de masse |
| `projection_project_mass_flux_periodic_x_neumann_y.m` | Opérateur historique périodique-x / Neumann-y |
| `projection_project_grid_periodic_x_neumann_y.m` | Projection de vitesse en canal |
| `projection_interpolate_grid_delta_to_particles.m` | Interpolation grille → particules |
| `projection_deposit_particles_to_grid.m` | Dépôt particules → grille |
| `projection_apply_cell_thermostat.m` | Thermostat cellulaire corrigé |

### Viriel et wallVP

| Fichier | Rôle |
|---|---|
| `mpcd_apply_virial_pressure_kick_channel.m` | Diagnostic et kick viriel actif |
| `mpcd_srd_collision_channel_virtual_walls.m` | Collision SRD avec particules virtuelles de paroi |
| `projection_piston_visualize_frame.m` | Visualisation piston |
| `reconstruct_grid_fields_from_runout_piston.m` | Reconstruction post-traitement |

### Autres validations

| Fichier | Rôle |
|---|---|
| `run_q9_reference_taylor_green_short.m` | Taylor-Green périodique |
| `run_q9_reference_lowk_mass_flux_30000.m` | Poiseuille Q9 historique |
| `run_q9_reference_step_short.m` | Step-channel / marche |
| `run_q9_reference_vonkarman_short.m` | Von Kármán court, non stabilisé en long |

---

## 5. Lancer un staircase EOS

Exemple de run long de validation :

```matlab
clear functions;
rehash;
close all;
clc;

runTag = datestr(now, 'yyyymmdd_HHMMSS');
outRoot = fullfile(pwd, ['overnight_staircase_eos_' runTag]);
mkdir(outRoot);

diary(fullfile(outRoot, 'console_overnight_staircase_eos.txt'));
diary on;

disp(which('run_q9_piston_staircase_eos_validation', '-all'));
disp(which('run_projection_piston_wallvp_demo', '-all'));
disp(which('mpcd_step_projection_piston', '-all'));
disp(which('mpcd_step_classic_piston', '-all'));
disp(which('mpcd_srd_collision_channel_virtual_walls', '-all'));
disp(which('mpcd_apply_virial_pressure_kick_channel', '-all'));
disp(which('projection_project_mass_flux_general_bc', '-all'));

lastwarn('');

out = run_q9_piston_staircase_eos_validation( ...
    'outputRoot', outRoot, ...
    'compressionLevels', [0 0.0025 0.005 0.01 0.02 0.035 0.05], ...
    'initialHoldSteps', 1000, ...
    'moveSteps', 150, ...
    'holdSteps', 1500, ...
    'averageHoldFraction', 0.5, ...
    'fitMinSamples', 3, ...
    'excludeInvalidPlateauxFromFit', true, ...
    'method', 'q9_full', ...
    'Nx', 32, ...
    'Ny', 32, ...
    'gamma', 20, ...
    'densityFactor', 0.8, ...
    'massFluxProjectionOperator', 'general_bc', ...
    'massFluxTargetFilter', 'elliptic_lowpass', ...
    'massFluxDensityRelaxationBeta', 5e-4, ...
    'massFluxFinalVelocityProjectionCleanup', false, ...
    'pistonQ9TargetMode', 'reference_gamma', ...
    'pistonQ9FreezeEllipticMetric', false, ...
    'virialDiagnosticsEnable', true, ...
    'virialKickEnable', true, ...
    'Kvirial', 0.05, ...
    'virialBeta', 0.1, ...
    'virialLimiterEnable', true, ...
    'virialMaxDuFractionThermal', 0.02, ...
    'sampleEvery', 50, ...
    'progressEvery', 250, ...
    'visualEnable', false, ...
    'storeDensityMaps', true);

[msg, wid] = lastwarn;

save(fullfile(outRoot, 'overnight_staircase_eos_output.mat'), ...
     'out', 'msg', 'wid', '-v7.3');

disp(wid);
disp(msg);
disp(out.summary);
disp(out.staircase.plateauTable);
disp(out.staircase.fitTable);

diary off;
```

Sorties typiques :

```text
console_overnight_staircase_eos.txt
staircase_plateaux.csv
staircase_fits.csv
staircase_output.mat
```

---

## 6. Diagnostics à examiner

### Bulk EOS

Diagnostics principaux :

```text
PkinMean
PvirMean
PtotMean
PdriveMean
PkinIdealRatio
PtotIdealRatio
KeffPkin
KeffPvir
KeffPtot
KeffPdrive
R2Ptot
R2Pdrive
```

Critères attendus :

```text
PtotMean = PkinMean + PvirMean
PdriveMean = PtotMean
KeffPtot ≈ KeffPkin + KeffPvir
R2Ptot ≈ 1
R2Pdrive ≈ 1
```

### Quasi-incompressibilité

Diagnostics principaux :

```text
rhoPhysicalRelError
stdN
lowKDensity
massFluxBefore
massFluxAfter
massFluxResidual
```

Critères attendus :

```text
rhoPhysicalRelError = 0
stdN borné
lowKDensity faible
pas de dérive thermique
```

### Kick viriel

Diagnostics principaux :

```text
virialDuRms
virialDuOverThermalRms
virialLimitedCellFraction
virialResidualMomentumKickNorm
```

Critères attendus :

```text
virialDuOverThermalRms << 1
virialResidualMomentumKickNorm ~ 1e-15
kBTCell ≈ 0.01
```

### WallVP

Les diagnostics wallVP sont maintenant algébriquement cohérents :

```text
pressureTopWallTotal = pressureTopWallImpact + pressureTopWallVPPositiveCompression
pressureTopWallTotalConsistencyResidual = 0
```

Mais leur interprétation physique n’est pas encore validée. Les valeurs instantanées et même certaines moyennes de plateau peuvent être fortement bruitées ou changer de signe.

Pour l’instant, ne pas utiliser `topWallPressureTotal` comme validation principale de l’EOS. La validation pression–densité doit reposer sur le bulk :

```text
PkinMean, PvirMean, PtotMean, PdriveMean
```

Le chantier wallVP restant consiste à isoler un cas plus simple, idéalement un gaz à l’équilibre entre deux murs sans piston ni viriel, afin de comparer la pression mécanique mesurée à `rho*kBT`.

---

## 7. Reproductibilité

### Ne pas modifier le path pendant un run

MATLAB peut recharger une fonction `.m` modifiée entre deux appels d’une boucle temporelle. Il ne faut donc pas copier de fichiers ou modifier une fonction située sur le path actif pendant un run long.

Avant chaque validation :

```matlab
clear functions;
rehash;

which run_q9_piston_staircase_eos_validation -all
which run_projection_piston_wallvp_demo -all
which mpcd_step_projection_piston -all
which mpcd_step_classic_piston -all
which mpcd_srd_collision_channel_virtual_walls -all
which mpcd_apply_virial_pressure_kick_channel -all
which projection_project_mass_flux_general_bc -all
```

Pour des runs longs, utiliser de préférence un snapshot de code dédié.

### Initialisation

Pour les validations piston et canal :

```matlab
params.initialPopulationMode = 'exact_per_cell';
```

Pour les domaines avec obstacles ou masques solides :

```matlab
params.initialPopulationMode = 'exact_per_fluid_cell';
```

---

## 8. Limitations connues

1. `wallVP pressure` est algébriquement cohérent mais pas encore physiquement validé.
2. Le cas cylindre / von Kármán reste non stabilisé en long.
3. `general_bc + elliptic_lowpass` est stable, mais le diagnostic `massFluxResidual` ne doit pas être interprété seul comme mesure physique.
4. La méthode Q6/Q9 modifie les propriétés de transport effectives ; ce n’est pas un échec si `nu_eff` diffère du SRC classique.
5. Le limiteur viriel est implémenté mais n’a pas été déclenché dans les runs de référence doux.

---

## 9. Prochaines étapes

Priorités recommandées :

1. figer la branche piston comme version bulk EOS fonctionnelle ;
2. ajouter un test dédié d’activation du limiteur avec seuil volontairement bas ;
3. construire un cas wallVP simple sans piston ni viriel pour valider la pression mécanique de paroi ;
4. reprendre ensuite les cas structurés : Poiseuille wallVP, step-channel, puis seulement von Kármán ;
5. documenter les commits de référence dans `BRANCH_REF.txt`.

---

## 10. Résumé opérationnel

La méthode actuelle permet déjà de tester un fluide SRC/MPCD quasi-incompressible avec fermeture pression–densité :

```text
compression globale conservée
densité low-k contrôlée
Pkin croît avec rho
Pvir ajoute une raideur réglable
Ptot et Pdrive sont cohérents
kick viriel actif stable
thermostat stable
```

Le verrou principal restant est la mesure mécanique de pression wallVP, qui doit être traitée comme diagnostic séparé et non comme critère bloquant pour l’EOS bulk.
