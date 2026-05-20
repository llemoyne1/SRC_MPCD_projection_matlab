# SRC_MPCD_projection_matlab — Q6/Q9 quasi-incompressible validation suite

Prototype MATLAB pour le développement et la validation d'une méthode SRC/MPCD quasi-incompressible en deux dimensions.

La branche de référence actuelle est :

```text
validation/integrated-q9-suite
```

État de référence au moment de l'intégration :

```text
201fc13 Integrate piston ramp wallVP diagnostics into clean Q6/Q9 validation base
```

Cette branche réunit :

```text
SRC/MPCD classique
→ projection de vitesse Q6 : div(u) ≈ 0
→ correction Q9 low-k du flux de masse div(Nu)
→ thermostat cellulaire corrigé
→ conditions limites solides avec particules virtuelles de paroi, wallVP
→ fermeture EOS virielle optionnelle : Ptot = Pkin + Pvir
```

Le dépôt sert maintenant à valider un modèle particulaire mésoscopique capable de conserver les avantages du SRC/MPCD — particules, fluctuations, transport mésoscopique — tout en réduisant les modes compressifs cohérents de densité et en produisant une réponse pression–densité de type liquide sous compression piston.

---

## 1. État courant

### 1.1 Modèle fluide de référence

Le modèle final testé n'est plus une simple comparaison de variantes indépendantes, mais la chaîne complète :

```text
q9_full = SRC/MPCD classic + Q6 + Q9 + thermostat corrigé
```

La fermeture virielle est ajoutée au-dessus de `q9_full` comme équation d'état effective, et non comme une méthode concurrente.

Les trois familles de validation utilisées sont :

| Cas | Conditions limites | Ce qui est validé |
|---|---|---|
| Taylor--Green forcé | périodique / sans mur | dynamique bulk, projection Q6, correction Q9, conservation des structures tourbillonnaires |
| Poiseuille | canal périodique-x / borné-y, wallVP | conditions limites solides, profil parabolique, viscosité effective, correction de grille en y |
| Piston rampe | piston mobile, wallVP, viriel | EOS liquide bulk, compression quasi-statique, diagnostic mécanique wall/bulk |

### 1.2 Jalons validés

État actuel de la méthode :

```text
- Q6/Q9 complet : fonctionnel en périodique, canal et piston
- thermostat cellulaire corrigé : stable
- opérateur general_bc : stabilisé par gauge-fix du mode constant
- Q9 low-k mass flux : réduction nette des modes compressifs cohérents
- Poiseuille wallVP : conditions de paroi améliorées
- correction de viscosité en taille de grille Ny : documentée
- piston smooth-ramp : compression quasi-statique fonctionnelle
- EOS virielle bulk : Pvir = Kvirial*(rho-rhoEOSRef)
- Ptot = Pkin + Pvir : validé côté bulk
- Pdrive = Ptot : validé côté bulk
- kick viriel actif faible : stable
- correction globale du moment du kick : résidu ~1e-15
- limiteur viriel : implémenté et diagnostiqué, non déclenché en régime doux
- diagnostic wallVP : cohérent à quelques pourcents près en régime post-rampe
```

### 1.3 Point encore à qualifier finement

Le diagnostic mécanique de pression murale `Pwall` est maintenant algébriquement cohérent et utilisable, mais il reste un diagnostic bruité. Sur le run piston-rampe court `32 x 32`, `9500` steps, les fenêtres strictement post-rampe donnent typiquement :

```text
Pwall_sym      ≈ 227.9
Ptot_layer_sym ≈ 221.0
Pwall-Ptot     ≈ +6.8, soit ≈ 3 %
```

Ce résultat est excellent au regard du bruit impulsionnel de paroi, mais il ne permet pas encore de séparer finement la contribution virielle locale lorsque celle-ci est seulement de l'ordre de `2-3` unités de pression. La validation fine de `Pwall` comme mesure indépendante de `Ptot` nécessitera un meilleur rapport signal/bruit : run plus long, grille plus résolue, compression/viriel plus marqué, ou moyennes indépendantes.

---

## 2. Méthode numérique

### 2.1 SRC/MPCD classique

Le fluide est représenté par des particules ponctuelles alternant :

```text
advection balistique + forçage
→ application des conditions limites
→ collision SRD/MPCD par rotation des vitesses relatives cellule par cellule
→ thermostat cellulaire optionnel/corrigé
```

La collision SRD conserve la quantité de mouvement cellulaire. Le thermostat cellulaire corrigé renormalise uniquement les fluctuations relatives pour imposer la température cible `kBT`, sans changer la vitesse moyenne de cellule.

### 2.2 Q6 — projection de vitesse

Q6 applique une projection de vitesse sur grille après l'étape SRC/MPCD classique :

```text
dépôt particules → grille
projection de vitesse : div(u) ≈ 0
interpolation de la correction grille → particules
correction globale du moment
```

Q6 contrôle la composante compressive du champ de vitesse, mais ne suffit pas toujours à contrôler les grandes structures de densité advectées par `Nu`.

### 2.3 Q9 — correction low-k du flux de masse

Q9 ajoute une correction du flux de masse :

```text
M = N u
```

afin de réduire le transport compressif cohérent :

```text
div(Nu)
```

La correction cible les modes basse fréquence de densité. Elle ne cherche pas à supprimer le bruit local d'occupation particulaire, qui fait partie du modèle mésoscopique.

Paramètres représentatifs :

```matlab
params.massFluxProjectionMode = 'relax_to_uniform_lowk';
params.massFluxProjectionStrength = 1.0;
params.massFluxDensityRelaxationBeta = 5.0e-4;
params.massFluxApplyAfterVelocityProjection = true;
params.massFluxProjectionOperator = 'general_bc';
params.massFluxTargetFilter = 'elliptic_lowpass';
params.massFluxLowKMaxIndex = 2;
params.massFluxFinalVelocityProjectionCleanup = false;
```

### 2.4 Opérateur `general_bc` et gauge-fix

L'opérateur `projection_project_mass_flux_general_bc.m` traite les géométries périodiques/bornées. Pour les cas fermés de type périodique/Neumann, le potentiel elliptique possède un mode constant arbitraire. Le solveur fixe désormais explicitement ce mode constant par gauge-fix plutôt que par une régularisation artificielle.

Cette correction ne change pas le champ physique, car seul le gradient du potentiel intervient.

### 2.5 Fermeture virielle EOS

La fermeture EOS utilise :

```matlab
Pvir   = Kvirial * (rho - rhoEOSRef);
Ptot   = Pkin + Pvir;
Pdrive = Ptot;
```

Le kick actif applique une correction de vitesse proportionnelle au gradient de pression :

```matlab
du = - virialBeta * dt ./ rhoKick .* grad(Pdrive);
```

Le kick est suivi d'une correction globale exacte du moment, puis d'un thermostat cellulaire final.

Diagnostics principaux :

```text
PkinMean
PvirMean
PtotMean
PdriveMean
virialDuRms
virialDuOverThermalRms
virialLimitedCellFraction
virialLimitedCellCount
virialLimiterDuMax
virialResidualMomentumKickNorm
```

En régime doux validé, le limiteur est exposé mais non déclenché :

```text
virialLimitedCellFraction = 0
virialResidualMomentumKickNorm ~ 1e-15
```

---

## 3. Conditions limites et wallVP

### 3.1 Conditions limites classiques

Les scripts utilisent selon les cas :

```text
périodique-x / périodique-y  : Taylor--Green
périodique-x / murs-y       : Poiseuille
périodique-x / piston-y     : piston compression
```

Les murs peuvent utiliser des réflexions de type `bounceback`, des parois thermalisées, et des particules virtuelles de paroi.

### 3.2 Particules virtuelles de paroi

Les particules virtuelles de paroi, wallVP, participent à la collision SRD dans les cellules proches du mur afin d'améliorer le couplage fluide/solide :

```text
particules fluides + particules virtuelles solides
→ collision SRD commune
→ reconstruction de l'impulsion échangée avec le mur
```

Le fichier central est :

```text
mpcd_srd_collision_channel_virtual_walls.m
```

Il reconstruit séparément :

```text
impulsion impact géométrique
impulsion wallVP collisionnelle
impulsion totale top/bottom
pressions signées
pressions positives de compression
pressions cumulées
```

### 3.3 Diagnostic mécanique Pwall

La pression murale ne doit pas être lue à partir d'un échantillon instantané. Elle est un flux d'impulsion :

```text
Pwall = DeltaI / (dt * Lx)
```

Avec `dt = 1e-3`, les valeurs instantanées peuvent fluctuer de plusieurs centaines ou milliers d'unités. Le diagnostic robuste utilise donc les impulsions cumulées et des fenêtres longues :

```text
PtopWall(Tw)    =  (Itop(t)-Itop(t-Tw))       / (Tw*Lx)
PbottomWall(Tw) = -(Ibottom(t)-Ibottom(t-Tw)) / (Tw*Lx)
```

Puis compare aux pressions bulk locales près des murs :

```text
Rtop    = PtopWall(Tw)    - <Pkin+Pvir>topLayer(Tw)
Rbottom = PbottomWall(Tw) - <Pkin+Pvir>bottomLayer(Tw)
Rsym    = 0.5*(Rtop + Rbottom)
Ranti   = 0.5*(Rtop - Rbottom)
```

Interprétation :

```text
Rsym  : biais moyen wall/bulk
Ranti : asymétrie top/bottom, donc transitoire ou incohérence de signe
```

Le protocole préféré pour valider wallVP est maintenant une rampe lisse suivie d'un hold final. Les fenêtres utilisées pour la conclusion doivent être entièrement situées dans le hold final.

---

## 4. Scripts principaux

### 4.1 Taylor--Green

| Fichier | Rôle |
|---|---|
| `run_projection_forced_taylor_green_demo.m` | Driver fonctionnel TG forcé |
| `run_q9_forced_taylor_green_validation.m` | Script de validation comparative classic/Q6/Q9 |
| `mpcd_step_classic_periodic_forced.m` | Étape SRC/MPCD périodique forcée |
| `mpcd_step_projection_periodic_q9_forced.m` | Étape périodique Q6/Q9 |
| `projection_taylor_green_diagnostics.m` | Diagnostics TG |
| `projection_fit_forced_taylor_green_viscosity.m` | Fit de viscosité TG |
| `projection_taylor_green_visualize_frame.m` | Visualisation TG |

### 4.2 Poiseuille

| Fichier | Rôle |
|---|---|
| `run_projection_poiseuille_medium64_demo.m` | Driver Poiseuille principal |
| `run_q9_poiseuille_wall_virtual_particles_validation.m` | Validation wallVP Poiseuille |
| `run_q9_poiseuille_wallvp_v2_q9_height_campaign.m` | Campagne hauteur/Ny wallVP + Q9 |
| `mpcd_step_classic_poiseuille.m` | Étape SRC/MPCD canal |
| `mpcd_step_projection_poiseuille_q9.m` | Étape canal Q6/Q9 |
| `analyze_projection_poiseuille_viscosity.m` | Analyse viscosité/profil |
| `postprocess_poiseuille_campaign_scaling.m` | Post-traitement correction Ny |
| `projection_poiseuille_visualize_frame.m` | Visualisation Poiseuille |

### 4.3 Piston, viriel et wallVP

| Fichier | Rôle |
|---|---|
| `run_projection_piston_wallvp_demo.m` | Driver piston général Q6/Q9, viriel, wallVP |
| `run_q9_piston_smooth_wallbulk_validation.m` | Validation piston rampe lisse + diagnostic wall/bulk rolling |
| `run_q9_piston_staircase_eos_validation.m` | Validation EOS par paliers de compression, surtout utile pour fits bulk |
| `run_wallvp_equilibrium_pressure_validation.m` | Test gaz stationnaire entre deux murs, wallVP off/on |
| `run_q9_piston_full_model_validation.m` | Campagne piston full-model historique |
| `mpcd_step_projection_piston.m` | Étape piston avec Q6/Q9 complet |
| `mpcd_step_classic_piston.m` | Étape SRC/MPCD piston + diagnostics wallVP |
| `mpcd_apply_virial_pressure_kick_channel.m` | Diagnostic et kick viriel actif |
| `analyze_piston_compressibility.m` | Fits `Pkin`, `Pvir`, `Ptot`, `Pdrive` vs densité |
| `projection_piston_visualize_frame.m` | Visualisation piston |

---

## 5. Lancer les validations principales

Avant un run long :

```matlab
clear functions;
rehash;
close all;
clc;
addpath(genpath(pwd));
```

Vérifier les fonctions utilisées :

```matlab
which run_projection_forced_taylor_green_demo -all
which run_projection_poiseuille_medium64_demo -all
which run_q9_piston_smooth_wallbulk_validation -all
which mpcd_srd_collision_channel_virtual_walls -all
which projection_project_mass_flux_general_bc -all
```

Ne pas modifier les fichiers `.m` pendant un run MATLAB long : MATLAB peut recharger une fonction modifiée entre deux appels.

### 5.1 Taylor--Green Q9

Exemple court ou moyen :

```matlab
params = struct();
params.method = 'q9';
params.Nx = 64;
params.Ny = 64;
params.gamma = 20;
params.dt = 1.0e-3;
params.kBT = 0.01;
params.nSteps = 30000;
params.sampleEvery = 100;
params.progressEvery = 1000;
params.taylorGreenInitialAmplitude = 0.10;
params.taylorGreenForceAmplitude = 0.12;
params.massFluxDensityRelaxationBeta = 5.0e-4;
params.massFluxLowKMaxIndex = 2;
params.visualEnable = false;
params.meanFieldEnable = true;

outTG = run_projection_forced_taylor_green_demo(params);
save('tg_q9_validation_output.mat','outTG','-v7.3');
```

Diagnostics à examiner :

```text
summary.finalAmplitude
summary.meanCoherence
summary.nuEffForcedFit
summary.finalKBTCell
series.lowKDensity
series.massFluxResidual
```

### 5.2 Poiseuille Q9 + wallVP

Exemple :

```matlab
resultsPois = run_q9_poiseuille_wall_virtual_particles_validation( ...
    'methods', {'q9'}, ...
    'Nx', 64, ...
    'Ny', 64, ...
    'gamma', 20, ...
    'nSteps', 30000, ...
    'sampleEvery', 100, ...
    'progressEvery', 1000, ...
    'visualEnable', false, ...
    'bodyForceX', 0.005, ...
    'wallModeY', 'bounceback', ...
    'geometryMode', 'shifted_solid_fraction', ...
    'densityFactor', 1.0, ...
    'thermal', true, ...
    'stochasticCount', true, ...
    'forceRandomShiftY', true, ...
    'initialPoiseuille', true, ...
    'nuGuess', 0.021, ...
    'nuReferenceNy', 64, ...
    'nuScalingEnable', true);

save('poiseuille_q9_wallvp_validation_output.mat','resultsPois','-v7.3');
```

Diagnostics à examiner :

```text
R2 du fit parabolique
nu_eff brut
nu_eff corrigé Ny
centerMinusWall
kBTCell
stdN
lowKDensity
massFluxBefore / massFluxAfter / massFluxResidual
```

### 5.3 Piston rampe lisse + wall/bulk

Exemple :

```matlab
outPiston = run_q9_piston_smooth_wallbulk_validation( ...
    'compressionFinal', 0.01, ...
    'initialHoldSteps', 3000, ...
    'rampSteps', 8000, ...
    'finalHoldSteps', 14000, ...
    'wallBulkResidualWindowSteps', 4000, ...
    'wallBulkResidualLayerCells', 3, ...
    'method', 'q9_full', ...
    'Nx', 32, ...
    'Ny', 32, ...
    'gamma', 20, ...
    'densityFactor', 0.8, ...
    'massFluxProjectionOperator', 'general_bc', ...
    'massFluxTargetFilter', 'elliptic_lowpass', ...
    'massFluxDensityRelaxationBeta', 5.0e-4, ...
    'massFluxFinalVelocityProjectionCleanup', false, ...
    'pistonQ9TargetMode', 'reference_gamma', ...
    'pistonQ9FreezeEllipticMetric', false, ...
    'virialDiagnosticsEnable', true, ...
    'virialKickEnable', true, ...
    'Kvirial', 0.01, ...
    'virialBeta', 0.02, ...
    'virialLimiterEnable', true, ...
    'virialMaxDuFractionThermal', 0.025, ...
    'sampleEvery', 100, ...
    'progressEvery', 1000, ...
    'visualEnable', false, ...
    'saveOutput', true, ...
    'plotEnable', true);

save('piston_q9_smooth_wallbulk_validation_output.mat','outPiston','-v7.3');
```

Sorties typiques :

```text
wall_bulk_timeseries.csv
wall_bulk_summary.csv
smooth_wallbulk_output.mat
wall_bulk_residual_vs_compression.png
wall_vs_bulk_pressure_timeseries.png
```

Diagnostics à examiner :

```text
PkinMean
PvirMean
PtotMean
PdriveMean
kBTCell
virialLimitedCellFraction
virialResidualMomentumKickNorm
pressureWallSymRollingAverage
PtotLayerSymRollingAverage
wallBulkResidualSymPtotRollingAverage
wallBulkResidualAntiPtotRollingAverage
```

---

## 6. Recalcul post-run du diagnostic wall/bulk

Le diagnostic `Pwall` doit être analysé sur des fenêtres strictement post-rampe. Une fenêtre est acceptable si :

```text
step - windowSteps >= rampEndStep
```

où :

```text
rampEndStep = initialHoldSteps + rampSteps
```

Exemple de principe :

```matlab
T = readtable(fullfile(outPiston.outputRoot,'wall_bulk_timeseries.csv'));
R = wall_bulk_recomputed_windows(T, [1500 2500 3500 5000]); % si helper disponible
```

Si aucun helper n'est disponible, recalculer les pressions murales par différence de cumul :

```text
I/Lx = pressureWallTotalCumulativeMean * wallImpulseTimeCumulative
Pwall_window = [I(t)-I(t-Tw)]/[tcum(t)-tcum(t-Tw)]
```

Les conclusions doivent privilégier :

```text
fenêtres longues mais entièrement dans le hold final
résidu symétrique Rsym
résidu antisymétrique Ranti
comparaison relative (Pwall-Ptot)/Ptot
```

---

## 7. Correction de viscosité Poiseuille en taille de grille

Les campagnes Poiseuille multi-résolutions ont montré que la viscosité brute obtenue par fit parabolique dépend fortement de `Ny` lorsque le domaine physique est conservé mais que la discrétisation transverse change.

Le diagnostic de référence conserve donc à la fois :

```text
nu_raw  : viscosité issue directement du fit
nu_corr : viscosité corrigée par la résolution transverse
```

Correction utilisée :

```text
nu_corr = nu_raw * (Ny / Ny_ref)^2
```

avec typiquement :

```text
Ny_ref = 64
```

Cette correction est essentielle pour comparer des runs à résolution transverse différente. Les valeurs corrigées se regroupent autour de la référence issue de Taylor--Green, alors que les valeurs brutes varient avec `Ny`.

---

## 8. Campagne intégrée longue recommandée

Pour vérifier qu'aucun développement récent n'a cassé les validations précédentes, lancer les trois tests dans le même snapshot de code :

```text
TG Q9 périodique
Poiseuille Q9 + wallVP
Piston Q9 + viriel + smooth-ramp wall/bulk
```

Budget indicatif pour une allocation de 15 h :

| Test | Steps indicatifs | Objectif |
|---|---:|---|
| Taylor--Green | 20000--30000 | dynamique bulk, projection, structures |
| Poiseuille | 25000--30000 | profil, viscosité, wallVP canal |
| Piston rampe | 25000--30000 | EOS liquide, wall/bulk rolling |

Désactiver la visualisation pour les runs de mesure :

```matlab
visualEnable = false;
```

ou utiliser un `visualEvery` très grand. La visualisation peut fortement ralentir les runs.

Une campagne autonome future doit écrire un dossier par test et un fichier `DONE.txt` après chaque cas afin de ne pas perdre les résultats déjà acquis en cas de coupure.

---

## 9. Reproductibilité Git

Vérifier la branche :

```powershell
git branch --show-current
git status
git log --oneline -5
```

État attendu :

```text
validation/integrated-q9-suite
working tree clean
```

Après modification du README :

```powershell
git add README.md
git commit -m "Update README for integrated Q6Q9 validation suite"
git push
```

---

## 10. Limitations connues

1. `Pwall` est un diagnostic impulsionnel bruité. Il doit être analysé par cumul/fenêtre, pas par valeurs instantanées.
2. Le diagnostic wall/bulk est cohérent à quelques pourcents près sur le cas piston-rampe court, mais ne sépare pas encore finement `Pkin` et `Pvir` lorsque `Pvir` est plus petit que l'incertitude effective de paroi.
3. Les validations TG, Poiseuille et piston utilisent des conditions limites différentes ; c'est un avantage pour la non-régression, mais impose de garder des scripts séparés.
4. La méthode Q6/Q9 modifie les propriétés de transport effectives ; une différence de viscosité avec SRC classique n'est pas automatiquement un échec.
5. Les diagnostics `massFluxResidual` et `lowKDensity` doivent être lus avec les autres observables : température, profils, pression, conservation globale.
6. Le cas von Karman/cylindre n'est pas inclus dans la validation finale de cette branche.

---

## 11. Résumé opérationnel

La branche `validation/integrated-q9-suite` fournit une base cohérente pour la validation finale MATLAB de la méthode :

```text
Q6/Q9 quasi-incompressible : validé côté bulk
thermostat corrigé        : stable
Poiseuille wallVP         : exploitable pour viscosité et profils
EOS virielle piston       : validée côté bulk
Pwall rolling             : cohérent avec Ptot local à quelques pourcents près
```

Le verrou restant n'est plus la méthode Q6/Q9/viriel elle-même, mais la métrologie fine de la pression mécanique de paroi lorsque l'on souhaite distinguer séparément `Pkin` et `Pvir` dans un signal wallVP bruité.
