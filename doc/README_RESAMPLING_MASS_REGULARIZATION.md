# Mass regularization preserving cell M, U and thermal energy

This patch adds a local particle-mass regularizer for the weighted-resampled
SRC/MPCD prototype.

The regularizer acts cell by cell on fluid-active particles only.  For each
wet cell, it stores the pre-regularization hydrodynamic quantities

- `M = sum(m)`;
- `U = sum(m*v)/M`;
- `Eth = 1/2 sum(m*|v-U|^2)`.

It then relaxes particle masses toward the uniform cell value `M/N`, projects
them into a bounded interval, restores the weighted mean velocity by a uniform
velocity shift, and optionally rescales fluctuations to restore the original
scalar thermal energy.  The target is therefore to preserve, up to roundoff,
cell mass, cell velocity and scalar thermal energy while reducing the spread of
particle masses.

## New files

- `resamp_regularize_particle_masses_preserve_mue.m`
- `run_resamp_poiseuille_wallvp_mass_regularization_once.m`

## Modified runner

`run_resamp_poiseuille_wallvp_validation.m` now accepts:

```matlab
'MassRegularizationEnable', true/false
'MassRegularizationStep', 3000
'MassRegularizationSteps', [3000 6000]
'MassRegularizationEvery', 0
'MassRegularizationStrength', 1.0
'MassRegularizationTriggerMode', 'always' | 'bounds_or_relstd' | 'bounds' | 'relstd'
'MassRegularizationRelStdTrigger', 0.20
'MassRegularizationMinFactor', 0.25
'MassRegularizationMaxFactor', 4.0
'MassRegularizationPreserveEnergy', true
```

Default is disabled, so previous runs are unchanged.

## Decisive comparison requested

To compare against the previous `q6_resampled_wallvp` 5000-step run while
changing only the addition of one regularization event at step 3000:

```matlab
outReg = run_resamp_poiseuille_wallvp_mass_regularization_once( ...
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

Inspect:

```matlab
outReg.cases.q6_resampled_wallvp.summary(:, {'step','mParticleRelStd', ...
    'massRegCells','massRegParticles','massRegRelStdBefore', ...
    'massRegRelStdAfter','massRegEnergyResidualRelRms', ...
    'massRegMomentumResidualRms','MRelRms','kBTWeighted'})

outReg.summary
outReg.cases.q6_resampled_wallvp.viscosity
```

Expected immediate effect at step 3000: lower particle-mass relative standard
deviation, with mass/momentum/energy residuals near roundoff.  The final
5000-step comparison indicates whether the improvement persists and whether the
Poiseuille viscosity/profile are affected.
