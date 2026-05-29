# Designed empty/rich pocket test for weighted resampling recycling

This step adds a deliberately constructed initial population edit for the
weighted resampling branch.  The goal is not to add another soft/hard
parameter sweep.  The goal is to create, from the start, the conditions that
must exercise the whole population-management method:

```text
empty pocket + overpopulated pocket
    -> extraction from rich cells
    -> insertion into empty cells
    -> local conservative mass/moment remap
    -> optional weighted thermostat
```

## New initial condition mode

The runner `run_resamp_pool_insertion_visual_demo` now accepts:

```matlab
'initialDepletion','paired_pockets'
```

This mode moves all active particles from one patch into another patch.  It
therefore keeps the total active particle count unchanged at `t=0`, but it
creates:

- a genuinely empty patch;
- a genuinely overpopulated patch.

The legacy mode remains available:

```matlab
'initialDepletion','patch'
```

which only deactivates particles in the empty patch.

## Dedicated runner

A short, intentionally diagnostic runner is provided:

```matlab
out = run_resamp_recycling_pockets_demo( ...
    'steps',20, ...
    'visualEvery',1, ...
    'summaryEvery',1);
```

The expected first population-edit event is:

```text
extractedParticles > 0
insertedParticles  > 0
Nactive approximately unchanged
MRelRMS near roundoff after remap
emptyCellsAfterInsert = 0
poorCellsAfterInsert  = 0
```

With the default `gamma=20`, `depletionPatchSize=[6 6]`, `NTarget=20`,
`NMin=10`, `NMax=30`, the empty patch contains 36 cells and the moved
particles make the rich patch roughly twice as populated.  The extraction
should therefore feed the free pool, and insertion should immediately refill
the empty pocket.

## Suggested command

```matlab
outP = run_resamp_recycling_pockets_demo( ...
    'steps',20, ...
    'visualEvery',1, ...
    'summaryEvery',1, ...
    'rngSeed',12345);

outP.summary(:, {'step','NpActive','Nfree','NMin','NMax', ...
    'extractedParticles','extractedParticlesCumulative', ...
    'insertedParticles','insertedParticlesCumulative', ...
    'poorCellsAfterInsert','emptyCellsAfterInsert', ...
    'overCellsAfterInsert','mParticleRelStd','MRelRms','kBTWeighted'})
```

This is the designed functional test for the closed pool loop.  More physical
cases should be assessed only after this test shows that the mechanism works
algebraically and visually.
