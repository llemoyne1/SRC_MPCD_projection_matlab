# Weighted resampling: overpopulated-cell extraction / recycling

This patch closes the first population-management loop for the weighted
resampling prototype.

Before this patch, the pool could only activate free slots in underpopulated
cells.  This repaired poor cells, but it could only increase `Nactive`.  It
also left overpopulated cells untouched, so local cell-mass constraints could
lead to progressively broad particle masses.

The new step adds extraction from overpopulated cells:

```text
SRC / Q6 step
save pre-edit weighted grid velocity U_target
extract particles from cells with N > NMax
insert particles into cells with N < NMin
local mass/moment remap toward U_target
post-remap weighted thermostat, optional
```

The extraction itself does not conserve local mass or momentum.  It is only a
population-support edit.  The following remap is responsible for restoring, in
each non-empty cell,

```text
sum_p m_p = M_target
sum_p m_p v_p = M_target U_target
```

where `U_target` is the weighted velocity deposited before extraction and
insertion.  This avoids making the population edit act as an artificial force.

## New file

- `resamp_extract_overpopulated_particles.m`

It deactivates particles in cells where `N > NMax`, down to `NTarget`.
The storage slots remain in the preallocated pool and can later be reused by
`resamp_insert_underpopulated_particles.m`.

Supported selection modes:

- `closest_to_cell_mean` nominal default; remove particles whose velocities are closest to the cell mean.
- `random` neutral diagnostic mode.
- `lightest` experimental.
- `heaviest` experimental.

## Runner options

`run_resamp_pool_insertion_visual_demo.m` now accepts:

```matlab
'ExtractEvery',1
'ExtractSelectionMode','closest_to_cell_mean'
'PreservePreEditVelocity',true
```

`PreservePreEditVelocity=true` means that, when extraction or insertion has actually modified the population at the current step, the following remap uses the grid velocity saved before the population edit.  If no population edit occurred, the normal remap mode is preserved.  This avoids changing the cheap nominal scale remap when the pool is idle.

## Suggested Taylor-Green smoke

```matlab
outR = run_resamp_pool_insertion_visual_demo( ...
    'method','weighted_q6', ...
    'steps',500, ...
    'visualEvery',10, ...
    'summaryEvery',25, ...
    'initialDepletion','none', ...
    'NMin',10, ...
    'NTarget',20, ...
    'NMax',30, ...
    'ExtractEvery',1, ...
    'InsertEvery',1, ...
    'ExtractSelectionMode','closest_to_cell_mean', ...
    'ThermostatAfterRemap',true, ...
    'ThermostatStrength',0.25, ...
    'ThermostatTargetKBT',0.01, ...
    'rngSeed',12345);
```

More aggressive population regularization:

```matlab
outR14 = run_resamp_pool_insertion_visual_demo( ...
    'method','weighted_q6', ...
    'steps',500, ...
    'visualEvery',10, ...
    'summaryEvery',25, ...
    'initialDepletion','none', ...
    'NMin',14, ...
    'NTarget',20, ...
    'NMax',26, ...
    'ExtractEvery',1, ...
    'InsertEvery',1, ...
    'ThermostatAfterRemap',true, ...
    'ThermostatStrength',0.25, ...
    'ThermostatTargetKBT',0.01, ...
    'rngSeed',12345);
```

## Diagnostics to inspect

The summary table includes:

```text
extractedParticles
extractCells
poorCellsAfterExtract
emptyCellsAfterExtract
overCellsAfterExtract
extractedParticlesCumulative
extractCellsCumulative
maxExtractedParticlesPerStep
lastExtractionStep
insertedParticlesCumulative
Nactive / Ncapacity / Nfree
NMin / NMax
mParticleRelStd
MRelRMS
kBTWeighted
```

The desired behavior for a closed recycling loop is not necessarily
`extractedParticles == insertedParticles` at every step, but rather:

```text
Nactive remains bounded
Nfree remains bounded
N[min,max] is controlled by [NMin,NMax]
mParticleRelStd stops drifting upward
MRelRMS stays near roundoff
kBTWeighted remains near target when the thermostat is enabled
```
