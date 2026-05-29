# Weighted-resampling fluid-physics validation

This patch adds two hydrodynamic validation cases for the weighted/resampled
SRC/MPCD prototype.

## 1. Forced periodic Taylor--Green

`run_resamp_tg_physics_validation.m` compares:

- `classic_reference`: weighted classic SRC/MPCD, no Q6, no remap/recycling;
- `q6_resampled`: weighted Q6 with extraction/insertion/remap, bounded-mass safety and post-remap thermostat.

The domain is fully wet and periodic.  It measures TG mode amplitude,
coherence, enstrophy, density fluctuation, divergence, particle-mass spread
and an effective TG viscosity estimate from the forced amplitude equation.

Example:

```matlab
outTG = run_resamp_tg_physics_validation( ...
    'steps',3000, ...
    'Nx',64, 'Ny',64, 'gamma',20, ...
    'NMin',14, 'NTarget',20, 'NMax',26, ...
    'TaylorGreenForcingAmplitude',0.12, ...
    'ThermostatStrength',0.25);
```

Outputs are written under:

```text
runs/resamp_tg_physics_validation/
```

## 2. Forced Poiseuille channel

`run_resamp_poiseuille_physics_validation.m` compares the same two families in
a periodic-x / bounded-y channel with `bodyForceX` and bounceback walls.  This
is intentionally a minimal weighted channel benchmark; it does not yet include
virtual wall particles.

Example:

```matlab
outP = run_resamp_poiseuille_physics_validation( ...
    'steps',5000, ...
    'Nx',64, 'Ny',32, 'gamma',20, ...
    'bodyForceX',0.005, ...
    'NMin',14, 'NTarget',20, 'NMax',26, ...
    'ThermostatStrength',0.25);
```

The runner stores weighted y-profiles and fits a quadratic Poiseuille profile
with `analyze_resamp_poiseuille_viscosity.m`.

## 3. Suite wrapper

For a quick smoke of both cases:

```matlab
out = run_resamp_fluid_physics_validation_suite('quick',true);
```

For longer runs:

```matlab
out = run_resamp_fluid_physics_validation_suite( ...
    'quick',false, ...
    'tgSteps',3000, ...
    'poiseuilleSteps',5000);
```

## Interpretation

The target is not exact equality with the fixed-mass SRC reference.  The goal
is to determine whether the weighted/resampled formulation gives a stable,
calibratable fluid with controlled density, bounded particle masses,
controlled temperature and reproducible effective transport coefficients.

The next step after this patch should be driven by these two questions:

1. Does the TG case preserve coherent vortical structures with stable effective viscosity?
2. Does the channel case produce a stable parabolic profile with a measurable `nu_eff`?
