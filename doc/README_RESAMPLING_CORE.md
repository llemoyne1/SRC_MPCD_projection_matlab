# Weighted resampling core — étape 1 : masses particulaires variables

Cette étape ajoute une infrastructure pondérée minimale au prototype MATLAB, sans introduire encore de création, suppression, split ou fusion de particules.

## Principe

L'état particulaire devient :

```matlab
state.x   % positions Np x 2
state.v   % vitesses  Np x 2
state.m   % masses    Np x 1
```

Par défaut `state.m(:)=1`. Dans ce cas, les dépôts pondérés doivent reproduire l'interprétation historique à particules de masse constante.

Le dépôt hydrodynamique devient :

```text
M_c  = sum_p m_p
P_c  = sum_p m_p v_p
U_c  = P_c / M_c
rho_c = M_c / (dx dy)
```

La collision SRC est effectuée autour du centre de masse pondéré :

```text
v'_p = U_c + R(alpha) (v_p - U_c)
```

ce qui conserve le moment pondéré cellule par cellule, à l'arrondi près.

## Fichiers ajoutés

- `resamp_initialize_particles_taylor_green_forced.m`
- `resamp_cell_ids_periodic.m`
- `resamp_deposit_weighted_to_grid.m`
- `resamp_step_classic_periodic_weighted.m`
- `resamp_apply_q6_periodic_weighted.m`
- `resamp_step_projection_periodic_weighted.m`
- `resamp_apply_global_momentum_correction_weighted.m`
- `resamp_apply_cell_thermostat_weighted.m`
- `resamp_population_mass_diagnostics.m`
- `run_resamp_weighted_mass_smoke.m`

## Exécution MATLAB

Depuis la racine du dépôt :

```matlab
out = run_resamp_weighted_mass_smoke('method','weighted_classic','steps',500);
out = run_resamp_weighted_mass_smoke('method','weighted_q6','steps',500);
```

Sortie CSV par défaut :

```text
runs/resamp_weighted_mass_smoke/resamp_weighted_mass_smoke_weighted_classic.csv
runs/resamp_weighted_mass_smoke/resamp_weighted_mass_smoke_weighted_q6.csv
```

## Limites assumées de cette étape

Cette étape ne fait pas encore de resampling :

- pas de création de particules ;
- pas de suppression/fusion de particules ;
- pas de réallocation locale des masses ;
- pas de contrainte stricte `M_c = M_target` ;
- pas de contrainte stricte `P_c = P_target`.

Elle sert uniquement à établir le socle pondéré nécessaire avant le resampling conservatif local.
