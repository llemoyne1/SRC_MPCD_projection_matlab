# Weighted resampling Poiseuille with wall virtual particles

This patch adds weighted virtual wall particles (wallVP) to the channel SRD collision used by the weighted-resampling Poiseuille validation.

The goal is not to add a new mass-management mechanism.  WallVPs are collision-only aggregates: they are not transported, not stored in the particle pool, not inserted/extracted, not remapped, and not counted as fluid mass.  They only modify the SRD collision-cell mean near the top and bottom walls so that the wall collision statistics impose a stronger no-slip behavior.

## Added files

- `resamp_srd_collision_channel_virtual_walls_weighted.m`
- `run_resamp_poiseuille_wallvp_validation.m`

## Modified file

- `resamp_step_classic_poiseuille_weighted.m`

The weighted channel step now calls `resamp_srd_collision_channel_virtual_walls_weighted`, which also handles the no-wallVP case when `wallVirtualParticlesEnable=false`.

## Collision model

For each collision cell, the real weighted fluid aggregates are

```text
M_real = sum_p m_p
P_real = sum_p m_p v_p
```

Near walls, virtual aggregate mass and momentum are added:

```text
M_total = M_real + M_virtual
P_total = P_real + P_virtual
U_collision = P_total / M_total
```

Only real fluid particles are then rotated around `U_collision`.

For thermal wallVPs, the aggregate virtual wall momentum is sampled with variance proportional to `M_virtual*kBT`, consistent with weighted momentum fluctuations.

## Nominal run

```matlab
outW = run_resamp_poiseuille_wallvp_validation( ...
    'steps',5000, ...
    'sampleEvery',50, ...
    'summaryEvery',100, ...
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
    'rngSeed',12345);
```

The runner compares:

- `classic_wallvp_reference`
- `q6_resampled_wallvp`

Outputs are written to:

```text
runs/resamp_poiseuille_wallvp_validation/
```

## Diagnostics to inspect

```matlab
outW.summary
outW.cases.q6_resampled_wallvp.viscosity
outW.cases.q6_resampled_wallvp.profileTable
outW.cases.q6_resampled_wallvp.summary(:, {'step','NMin','NMax','MRelRms', ...
    'mParticleRelStd','kBTWeighted','centerMinusWall', ...
    'wallMeanVelocity','nVirtualWallParticlesEquivalent'})
```

The expected improvement over the previous minimal-wall Poiseuille case is a much lower wall slip and a more parabolic late-time profile.
