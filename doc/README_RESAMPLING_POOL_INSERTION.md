# Resampling step 3 — preallocated pool and underpopulated-cell insertion

This patch adds the first memory-oriented population repair layer for the weighted SRC/MPCD prototype.

## Scope

The patch implements:

- a preallocated particle pool through `state.active`, `state.Nactive`, `state.Ncapacity`;
- active-mask support in weighted deposit, diagnostics, weighted SRC collision, Q6 projection, thermostat and mass/moment remap;
- insertion of nominal-mass particles into cells with `N < NMin`;
- velocity memory fields `state.uMemUx`, `state.uMemUy`, `state.uMemValid` for empty or strongly underpopulated cells;
- a smoke runner with an optional initial population hole.

The patch deliberately does **not** implement overpopulation fusion/removal. Overpopulated cells are diagnosed only.

## New files

- `resamp_active_mask.m`
- `resamp_enable_particle_pool.m`
- `resamp_particle_pool_info.m`
- `resamp_insert_underpopulated_particles.m`
- `run_resamp_pool_insertion_smoke.m`
- `doc/README_RESAMPLING_POOL_INSERTION.md`

## Modified files

- `resamp_deposit_weighted_to_grid.m`
- `resamp_population_mass_diagnostics.m`
- `resamp_step_classic_periodic_weighted.m`
- `resamp_apply_q6_periodic_weighted.m`
- `resamp_apply_cell_thermostat_weighted.m`
- `resamp_local_mass_moment_remap.m`

## Default policy

For `gamma = 20`, the default policy is equivalent to:

```matlab
NTarget = gamma;
NMin    = ceil(0.5*gamma);
NMax    = ceil(1.5*gamma);
capacityFactor = 2.0;
```

If a cell has `N < NMin`, the insertion layer tries to activate enough free pool slots to reach `NTarget`.

New particles receive:

- random positions inside the target cell;
- nominal mass `m0 = params.resampParticleMass`;
- velocity equal to the current reliable cell velocity when available, otherwise the stored memory velocity;
- zero-mean thermal fluctuations around that velocity.

The subsequent remap stage imposes the exact target cell mass and preserves the selected cell velocity.

## MATLAB smoke commands

Classic weighted path with an initial 6x6 depleted patch:

```matlab
outC = run_resamp_pool_insertion_smoke( ...
    'method','weighted_classic', ...
    'steps',200, ...
    'summaryEvery',20, ...
    'initialDepletion','patch');
```

Q6 weighted path:

```matlab
outQ = run_resamp_pool_insertion_smoke( ...
    'method','weighted_q6', ...
    'steps',200, ...
    'summaryEvery',20, ...
    'initialDepletion','patch');
```

No artificial initial depletion, useful to monitor natural underpopulation and capacity use:

```matlab
outN = run_resamp_pool_insertion_smoke( ...
    'method','weighted_q6', ...
    'steps',500, ...
    'summaryEvery',50, ...
    'initialDepletion','none');
```

CSV output goes to:

```text
runs/resamp_pool_insertion_smoke/
```

## Main diagnostics

The smoke CSV includes:

- `NpActive`, `Ncapacity`, `Nfree`, `activeFraction`;
- `NMin`, `NMax`, `NStd`, `nEmptyCells`;
- `insertedParticles`, `insertCells`, `poorCellsAfterInsert`, `capacityHit`;
- `MRelRms`, `remapMassResidualRelRms`, `remapMomentumResidualRms`;
- `mParticleRelStd`, `mParticleMin`, `mParticleMax`.

The acceptance target for this patch is not physical validation yet.  It is only to verify that the pool/insertion/remap chain is algebraically coherent and that the active particle count remains bounded by the preallocated capacity.
