# Weighted resampling — local mass/moment remap, step 2

This patch adds the first conservative remapping layer for the MATLAB
weighted SRC/MPCD prototype.

## Scope

The patch only changes particle masses.  It does **not** create, delete,
split or merge particles.

Added files:

- `resamp_solve_cell_mass_moment_weights.m`
- `resamp_local_mass_moment_remap.m`
- `run_resamp_local_remap_smoke.m`

The historical Q6/Q9 MATLAB files are not modified.

## Principle

For each non-empty cell, the remap seeks particle masses such that

```text
sum_p m_p = M_target
sum_p m_p v_p = P_target
```

The default target is deliberately conservative:

```text
P_target = M_target * U_old
U_old    = sum_p m_old v_p / sum_p m_old
```

Thus the remap enforces the target cell mass while preserving the current
weighted cell velocity.  No velocity kick and no artificial pressure force
is introduced by this step.

## Methods

### `scale_preserve_velocity`

Default and most robust mode.  For the default target velocity, all masses
inside a cell are multiplied by the same factor:

```text
m_new = m_old * M_target / M_old
```

This exactly enforces the cell target mass and preserves the weighted cell
velocity.  If bounds are violated, the code falls back to the bounded
minimum-change solver.

### `min_change`

Uses `resamp_solve_cell_mass_moment_weights` to solve a local bounded
least-change problem:

```text
minimize ||m_new - m_old||_2
subject to mass and momentum constraints
         massMin <= m_new <= massMax
```

No Optimization Toolbox is required.

## Smoke runs

From MATLAB at repository root:

```matlab
outC = run_resamp_local_remap_smoke( ...
    'method','weighted_classic', ...
    'steps',500, ...
    'summaryEvery',50);
```

```matlab
outQ = run_resamp_local_remap_smoke( ...
    'method','weighted_q6', ...
    'steps',500, ...
    'summaryEvery',50);
```

The script writes CSV summaries under:

```text
runs/resamp_local_remap_smoke/
```

Key columns:

- `MRelRms`: RMS cell mass error relative to target cell mass.
- `mParticleRelStd`: relative standard deviation of particle masses.
- `remapMassResidualRelRms`: remap constraint residual for mass.
- `remapMomentumResidualRms`: remap constraint residual for momentum.
- `remapSuccessFractionNonEmpty`: fraction of non-empty cells satisfying the constraints.
- `remapCellsEmpty`: cells that cannot be remapped because no particle is present.

## Interpretation

This step tests the algebraic core of weighted resampling only.  It is not
yet a complete solution for poor/empty cells, because empty cells cannot be
repaired without particle creation.  That is the next layer.
