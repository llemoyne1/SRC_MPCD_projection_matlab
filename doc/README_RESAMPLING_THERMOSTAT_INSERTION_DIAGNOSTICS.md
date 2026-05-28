# Weighted resampling: insertion counters and post-remap thermal reinjection

This note documents patch step 5 for the weighted-resampling MATLAB branch.

## Purpose

The previous visualization displayed only the number of particles inserted at the current visualized step. This was misleading for runs where the initial depleted patch is repaired immediately at step 1 and no further insertion occurs. This patch adds cumulative insertion counters and displays the last insertion step.

The patch also adds an explicit post-insertion/post-remap thermal reinjection option. This is separate from the historical `thermostatAfterProjection` flag used inside the Q6 step. The new option is intended for weighted resampling runs, where the preferred ordering is:

```text
weighted SRC/Q6 step
underpopulated-cell insertion
local mass/momentum remap
optional weighted thermostat after remap
```

## New visual-demo options

```matlab
'ThermostatAfterRemap', true/false
'ThermostatTargetKBT', value
'ThermostatStrength', value
```

The thermostat is weighted and cell-local. It rescales fluctuations around the weighted cell mean velocity and therefore preserves the weighted cell momentum.

## New diagnostics in the CSV summary

```text
insertedParticlesCumulative
insertCellsCumulative
maxInsertedParticlesPerStep
lastInsertionStep
thermostatAfterRemapEnabled
thermostatAfterRemapMeanKBTBefore
thermostatAfterRemapMeanKBTAfter
thermostatAfterRemapMeanScale
thermostatAfterRemapRmsVelocityChange
```

## Visualization changes

The title now reports:

```text
insertedNow, cum, last, thermAfter
```

The population-class panel now uses a fixed color axis and explicit labels:

```text
-2 empty, -1 poor, 0 ok, 1 over, 2 inserted
```

This makes comparisons across `NMin`, `insertEvery`, and thermostat settings less ambiguous.

## Suggested commands

Baseline without thermal reinjection:

```matlab
out0 = run_resamp_pool_insertion_visual_demo( ...
    'method','weighted_q6', ...
    'steps',300, ...
    'visualEvery',5, ...
    'summaryEvery',25, ...
    'initialDepletion','patch', ...
    'NMin',10, ...
    'NTarget',20, ...
    'ThermostatAfterRemap',false, ...
    'rngSeed',12345);
```

Moderate post-remap reinjection:

```matlab
outT25 = run_resamp_pool_insertion_visual_demo( ...
    'method','weighted_q6', ...
    'steps',300, ...
    'visualEvery',5, ...
    'summaryEvery',25, ...
    'initialDepletion','patch', ...
    'NMin',10, ...
    'NTarget',20, ...
    'ThermostatAfterRemap',true, ...
    'ThermostatStrength',0.25, ...
    'ThermostatTargetKBT',0.01, ...
    'rngSeed',12345);
```

Stronger reinjection:

```matlab
outT50 = run_resamp_pool_insertion_visual_demo( ...
    'method','weighted_q6', ...
    'steps',300, ...
    'visualEvery',5, ...
    'summaryEvery',25, ...
    'initialDepletion','patch', ...
    'NMin',10, ...
    'NTarget',20, ...
    'ThermostatAfterRemap',true, ...
    'ThermostatStrength',0.5, ...
    'ThermostatTargetKBT',0.01, ...
    'rngSeed',12345);
```

## Expected interpretation

A useful regime should keep:

```text
MrelRMS near machine precision
N[min,max] near the target band
mRelStd moderate
coherent Q6 vortical structures
kBT not collapsed to the unthermostatted Q6 value
```

If the thermostat strength is too high, the reinjected fluctuations may obscure coherent low-wavenumber structures. In that case, reduce `ThermostatStrength` or use a smaller target value.
