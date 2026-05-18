# SRC/MPCD Q9 — Taylor–Green & Poiseuille clean validation

This repository contains the clean MATLAB version of the SRC/MPCD projection framework focused on two validated benchmark cases:

1. **Periodic forced Taylor–Green flow** for bulk viscosity and Q9 non-regression.
2. **Poiseuille channel flow with virtual wall particles** for wall/no-slip validation and comparison with the Taylor–Green bulk viscosity.

The repository is intentionally restricted to the scripts required for these two cases. Obsolete free-surface, dam-break, inclined-layer, step-channel and intermediate debug scripts have been removed from this clean branch.

---

## 1. Numerical methods included

### Classic SRC/MPCD

The baseline method is standard stochastic rotation dynamics / MPCD:

- particle streaming,
- optional forcing,
- random grid shift,
- cell-wise SRC rotation,
- thermal control,
- diagnostics and visualization.

### Q6 velocity projection

Q6 adds an incompressibility projection of the grid velocity field:

\[
\nabla \cdot \mathbf{u} \simeq 0.
\]

It is used mostly as a reference between classic SRC and the Q9 method.

### Q9 mass-flux projection

Q9 adds a low-wavenumber correction of the mass flux:

\[
\mathbf{J} = N \mathbf{u},
\]

designed to reduce slow density/compressibility modes while preserving global momentum.

In periodic Taylor–Green, Q9 uses FFT-based periodic operators.  
In Poiseuille channel flow, Q9 uses the current elliptic/general boundary-condition operator.

### Wall virtual particles, wallVP-v2

The Poiseuille channel uses a shifted-grid virtual wall particle model:

```matlab
params.wallVirtualParticlesEnable = true;
params.wallVirtualParticlesGeometryMode = 'shifted_solid_fraction';
params.wallVirtualParticlesDensityFactor = 0.8;
```

Virtual particles:

- participate in the SRC collision average,
- have wall mean velocity,
- have thermal fluctuations at `kBT`,
- are not stored in the real particle state,
- are not advected,
- are discarded after collision.

This restores the momentum balance at the wall and removes the previous excessive slip observed in Poiseuille.

---

## 2. Main entry points

### Taylor–Green periodic benchmark

```matlab
run_projection_forced_taylor_green_demo.m
```

Optional validation wrapper:

```matlab
run_q9_forced_taylor_green_validation.m
```

### Poiseuille channel benchmark

```matlab
run_projection_poiseuille_medium64_demo.m
```

Wall virtual particle validation:

```matlab
run_q9_poiseuille_wall_virtual_particles_validation.m
```

Q9-inclusive height convergence campaign:

```matlab
run_q9_poiseuille_wallvp_v2_q9_height_campaign.m
```

Poiseuille post-processing with viscosity scaling:

```matlab
postprocess_poiseuille_campaign_scaling.m
```

---

## 3. Taylor–Green non-regression

The current reference Q9 Taylor–Green run is:

```matlab
params = struct();

params.method = 'q9';

params.Lx = 1.0;
params.Ly = 1.0;
params.Nx = 64;
params.Ny = 64;
params.gamma = 20;
params.seed = 11;

params.dt = 1.0e-3;
params.kBT = 0.01;
params.alphaDeg = 90;

params.nSteps = 30000;
params.sampleEvery = 100;
params.progressEvery = 1000;

params.taylorGreenInitialAmplitude = 0.10;
params.taylorGreenAmplitude = 0.10;
params.taylorGreenThermalNoise = true;
params.taylorGreenForceEnable = true;
params.taylorGreenForceAmplitude = 0.12;
params.taylorGreenForceZeroMeanKick = true;
params.taylorGreenModeX = 1;
params.taylorGreenModeY = 1;

params.initialPopulationMode = 'exact_per_cell';
params.initialVelocityZeroGlobalMean = true;
params.useRandomGridShift = true;

params.thermostatAfterStep = true;
params.thermostatAfterProjection = true;

params.projectionEnable = true;
params.projectionStrength = 1.0;
params.projectionInterpolationMethod = 'nearest';
params.projectionMomentumCorrectionEnable = true;
params.projectionMomentumCorrectionMode = 'particle_global_exact';

params.massFluxProjectionMode = 'relax_to_uniform_lowk';
params.massFluxProjectionStrength = 1.0;
params.massFluxDensityRelaxationBeta = 5.0e-4;
params.massFluxApplyAfterVelocityProjection = true;
params.massFluxLowKMaxIndex = 2;
params.lowKMaxIndex = 2;
params.massFluxFinalVelocityProjectionCleanup = true;
params.massFluxFinalVelocityProjectionStrength = 0.5;

params.taylorGreenViscosityFitFraction = 0.60;

params.visualEnable = true;
params.visualEvery = 500;
params.visualPause = 0.001;
params.visualMaxParticles = 8000;

out = run_projection_forced_taylor_green_demo(params);
```

### Current reference result

For `Nx=Ny=64`, `gamma=20`, `dt=1e-3`, `kBT=0.01`, `F_TG=0.12`:

```text
Q9 Taylor–Green:
    nu_eff preferred      ≈ 0.02050
    nu_eff derivative     ≈ 0.02052
    nu_eff plateau        ≈ 0.02052
    kBT                   ≈ 0.01
    final low-k density   ≈ 1.1e-6
    final coherence       ≈ 0.87
    final high-k fraction ≈ 0.10
```

This is the current bulk-viscosity non-regression reference for Q9.

---

## 4. Poiseuille wallVP-v2 validation

The preferred Poiseuille wall model is:

```matlab
params.wallModeY = 'bounceback';

params.wallVirtualParticlesEnable = true;
params.wallVirtualParticlesGeometryMode = 'shifted_solid_fraction';
params.wallVirtualParticlesDensityFactor = 0.8;
params.wallVirtualParticlesThermal = true;
```

The main validation campaign is:

```matlab
results = run_q9_poiseuille_wallvp_v2_q9_height_campaign();
```

Default campaign:

```text
densityFactor = 0.8
geometryMode  = shifted_solid_fraction
bodyForceX    = 0.005
nuGuess       = 0.032
visualEnable  = true

Runs:
    CLASSIC 32 x 48
    Q9      32 x 48
    CLASSIC 48 x 64
    Q9      48 x 64
    CLASSIC 64 x 96
    Q9      64 x 96
```

The campaign writes:

```text
console_log.txt
summary.csv
summary.mat
campaign_results.mat
one .mat file per run
```

---

## 5. Poiseuille viscosity scaling

The raw Poiseuille fit gives:

\[
\nu_\mathrm{raw}
=
-\frac{f_x}{2a_2},
\]

where \(a_2\) is the quadratic coefficient of the fitted profile.

When comparing different `Ny`, the viscosity must be scaled consistently with the wall-normal cell size. The current post-processing reports:

```matlab
nuEffRaw
nuEffScaledNyRef
nuEffCellY
dy
nuScaleFactorNyRef
```

with:

```matlab
nuEffScaledNyRef = nuEffRaw * (Ny / NyRef)^2;
```

By default:

```matlab
NyRef = 64;
```

To reprocess a campaign:

```matlab
summaryScaled = postprocess_poiseuille_campaign_scaling( ...
    'E:\GitHub\SRC_MPCD_projection_matlab_par\q9_wallvp_v2_q9_height_20260518_084250\campaign_results.mat', ...
    'nuReferenceNy', 64, ...
    'fitWindowTime', 5, ...
    'accelerationWindowTime', 5);
```

This writes:

```text
summary_scaled.csv
```

---

## 6. Current Poiseuille validation state

With `densityFactor=0.8`, `bodyForceX=0.005`, wallVP-v2 and `NyRef=64`, the current scaled viscosity results are:

| Method | Grid | `nuEffRaw` | `nuEffScaledNyRef` | `R2` | `ratioAccelRecent` |
|---|---:|---:|---:|---:|---:|
| classic | 32×48 | 0.03489 | 0.01963 | 0.9706 | -0.057 |
| Q9 | 32×48 | 0.03614 | 0.02033 | 0.9853 | -0.085 |
| classic | 48×64 | 0.01960 | 0.01960 | 0.9950 | -0.020 |
| Q9 | 48×64 | 0.02114 | 0.02114 | 0.9941 | -0.037 |
| classic | 64×96 | 0.00948 | 0.02133 | 0.9991 | +0.067 |
| Q9 | 64×96 | 0.01038 | 0.02336 | 0.9979 | +0.078 |

The most reliable current Poiseuille reference point is:

```text
Q9, 48 x 64:
    nuEffScaledNyRef ≈ 0.02114
    R2               ≈ 0.994
    ratioAccelRecent ≈ -0.037
    final low-k      ≈ 2.9e-6
```

This is consistent with the current Taylor–Green Q9 bulk viscosity:

```text
nu_TG,Q9 ≈ 0.0205
```

The `64 x 96` runs still show a small positive acceleration ratio, so they should be considered not fully converged yet. They are stable and physically consistent, but a longer run is recommended for final confirmation.

---

## 7. Important interpretation

The old Poiseuille wall model was not valid as a no-slip benchmark: classic, Q6 and Q9 all showed excessive global acceleration, with approximately:

```text
d<ux>/dt / bodyForceX ≈ 0.8
```

After wallVP-v2, the acceleration ratio is reduced to a few percent. Therefore:

```text
Poiseuille before wallVP-v2 was not a reliable viscosity benchmark.
Poiseuille with wallVP-v2 is now consistent with Taylor–Green after Ny scaling.
```

Taylor–Green remains the cleanest bulk-viscosity benchmark because it is periodic and does not involve walls.

Poiseuille now validates:

1. the wall model,
2. the no-slip momentum balance,
3. the compatibility of Q9 with bounded domains,
4. the consistency of wall-bounded viscosity with periodic bulk viscosity.

---

## 8. Visualization policy

Visualizations are enabled by default in validation scripts.

This is intentional: live plots of

```text
meanUx
Ucenter - Uwall
Umax
acceleration ratio
low-k density
kBT
Poiseuille profile
Taylor–Green amplitude/coherence
```

are essential to detect slow drifts, non-stationarity and wall problems early.

To disable visualization explicitly:

```matlab
params.visualEnable = false;
```

---

## 9. Minimal clean script set

The clean repository should contain only the scripts required for Taylor–Green and Poiseuille validation.

### Main scripts

```text
run_projection_forced_taylor_green_demo.m
run_projection_poiseuille_medium64_demo.m
run_q9_forced_taylor_green_validation.m
run_q9_poiseuille_wall_virtual_particles_validation.m
run_q9_poiseuille_wallvp_v2_q9_height_campaign.m
postprocess_poiseuille_campaign_scaling.m
```

### Poiseuille-specific scripts

```text
projection_initialize_particles_poiseuille.m
mpcd_step_classic_poiseuille.m
mpcd_step_projection_poiseuille_q9.m
mpcd_apply_q9_projection_channel.m
mpcd_apply_wall_bc_y.m
mpcd_srd_collision_channel_virtual_walls.m
projection_project_grid_periodic_x_neumann_y.m
projection_project_mass_flux_channel_operator.m
projection_project_mass_flux_general_bc.m
projection_project_mass_flux_periodic_x_neumann_y.m
analyze_projection_poiseuille_viscosity.m
poiseuille_acceleration_diagnostics.m
projection_poiseuille_visualize_frame.m
```

### Taylor–Green-specific scripts

```text
projection_initialize_particles_taylor_green_forced.m
mpcd_step_classic_periodic_forced.m
mpcd_step_projection_periodic_q9_forced.m
mpcd_apply_q9_projection_periodic.m
projection_apply_taylor_green_forcing.m
projection_project_grid_periodic_fft.m
projection_project_mass_flux_periodic_fft.m
projection_fit_forced_taylor_green_viscosity.m
projection_taylor_green_diagnostics.m
projection_taylor_green_mode_at_points.m
projection_taylor_green_visualize_frame.m
```

### Shared scripts

```text
projection_apply_cell_thermostat.m
projection_apply_global_momentum_correction.m
projection_deposit_particles_to_grid.m
projection_interpolate_grid_delta_to_particles.m
projection_population_diagnostics.m
projection_thermal_diagnostics.m
```

---

## 10. Recommended next validation steps

1. Keep Q9 Taylor–Green as the bulk reference:

```text
nu_TG,Q9 ≈ 0.0205
```

2. Use Q9 Poiseuille `48 x 64` wallVP-v2 as the current wall-bounded reference:

```text
nu_scaled ≈ 0.0211
```

3. Later, relaunch a longer Q9 `64 x 96` Poiseuille run to confirm asymptotic convergence.

4. When moving toward C++/OpenMP, preserve the following design principles:

```text
- virtual wall particles must enter the collision average only;
- they must not be stored as real particles;
- Q9 must preserve global momentum correction;
- wall/boundary handling should remain parameter-driven;
- visual diagnostics should remain available by default during validation.
```
