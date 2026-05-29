# Weighted resampling periodic-step Kolmogorov forcing

This patch adds a mean-free Kolmogorov body force to the periodic-step weighted-resampling demo.

The goal is not to introduce wall-driven Poiseuille physics. The periodic-step case is fully periodic, so a uniform `bodyForceX` would mostly accelerate the box mean momentum. Instead, the new forcing is spatially periodic and has zero continuum mean:

```text
a_x(y) = A sin(2*pi*k*y/Ly + phase),   a_y = 0
```

The finite-particle kick can also be corrected to have exactly zero mass-weighted mean momentum increment. This is enabled by default.

## New files

- `resamp_apply_kolmogorov_body_force.m`

## Modified files

- `run_resamp_periodic_step_visual_demo.m`
- `resamp_periodic_step_visualize_frame.m`

## New options for `run_resamp_periodic_step_visual_demo`

```matlab
'ForceMode','kolmogorov_x'
'KolmogorovForceAmplitude',0.02
'KolmogorovForceWaveNumber',1
'KolmogorovForcePhase',0.0
'KolmogorovForceDirection','x'
'KolmogorovForceZeroMeanKick',true
'KolmogorovForceApplyDt',true
```

Aliases are accepted for some options, for example `ForceAmplitude` and `KolmogorovAmplitude`.

## Step ordering

For `weighted_q6`, the runner now splits the step so that the forcing is applied after SRC streaming/collision and before the Q6 projection:

```text
SRC weighted streaming/collision
Kolmogorov body-force kick
Q6 projection
insertion
local mass/momentum remap
post-remap weighted thermostat
```

For `weighted_classic`, the force is applied after the weighted SRC step. There is no Q6 projection in that mode.

## Example

```matlab
outK = run_resamp_periodic_step_visual_demo( ...
    'method','weighted_q6', ...
    'steps',1000, ...
    'visualEvery',10, ...
    'summaryEvery',25, ...
    'Nx',64, ...
    'Ny',32, ...
    'gamma',20, ...
    'dt',0.01, ...
    'ULower',0.20, ...
    'UUpper',-0.05, ...
    'ForceMode','kolmogorov_x', ...
    'KolmogorovForceAmplitude',0.02, ...
    'KolmogorovForceWaveNumber',1, ...
    'KolmogorovForceZeroMeanKick',true, ...
    'NMin',10, ...
    'NTarget',20, ...
    'ThermostatAfterRemap',true, ...
    'ThermostatStrength',0.25, ...
    'ThermostatTargetKBT',0.01, ...
    'rngSeed',12345);
```

A more aggressive forcing can use `KolmogorovForceAmplitude = 0.05` or `0.1`, but this should be interpreted as a stress test of the weighted resampling machinery.

## Diagnostics added to CSV

- `forceEnabled`
- `forceMode`
- `kolmogorovForceAmplitude`
- `kolmogorovForceWaveNumber`
- `kolmogorovKickRms`
- `kolmogorovKickMax`
- `kolmogorovMomentumDeltaNorm`
- `kolmogorovMeanKickX`
- `kolmogorovMeanKickY`

For a zero-mean kick, `kolmogorovMomentumDeltaNorm` should remain near roundoff.
