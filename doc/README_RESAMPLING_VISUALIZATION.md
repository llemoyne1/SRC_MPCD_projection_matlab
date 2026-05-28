# Weighted resampling visualization

This note documents the first visualization layer for the weighted-resampling
prototype.

The historical Taylor-Green visualizer is useful for velocity and modal
content, but it assumes that particle count is also the density proxy.  In the
weighted formulation, this is no longer true.  The diagnostic view must show,
side by side:

- active particles colored by individual mass `m_p`;
- cell population `N/gamma - 1`;
- weighted cell mass `M/M_target - 1`;
- weighted velocity field and quiver;
- vorticity from the weighted velocity field;
- poor/nominal/overpopulated/inserted cell class.

## New files

- `resamp_pool_visualize_frame.m`
- `run_resamp_pool_insertion_visual_demo.m`

The visual demo mirrors `run_resamp_pool_insertion_smoke.m`, but calls the new
visualizer at step 0 and every `visualEvery` steps.  It does not change the
resampling, insertion, remap or Q6 algorithms.

## Basic MATLAB commands

```matlab
outV = run_resamp_pool_insertion_visual_demo( ...
    'method','weighted_q6', ...
    'steps',300, ...
    'visualEvery',5, ...
    'summaryEvery',25, ...
    'initialDepletion','patch');
```

Classic weighted reference:

```matlab
outC = run_resamp_pool_insertion_visual_demo( ...
    'method','weighted_classic', ...
    'steps',300, ...
    'visualEvery',5, ...
    'summaryEvery',25, ...
    'initialDepletion','patch');
```

Debug view and PNG export:

```matlab
outD = run_resamp_pool_insertion_visual_demo( ...
    'method','weighted_q6', ...
    'steps',200, ...
    'visualEvery',10, ...
    'showDebugFigure',true, ...
    'saveFrames',true, ...
    'frameDir',fullfile('runs','resamp_pool_insertion_visual','frames'));
```

## Main figure interpretation

Panel 1 shows active particles colored by `m_p/m0`.  This is the first place to
look for mass outliers.

Panel 2 shows `N/gamma - 1`.  It may remain non-zero even if the method works,
because population and mass are no longer identical.

Panel 3 shows `M/M_target - 1`.  After the local remap, this should remain close
to zero except in genuinely unresolved/capacity-hit cells.

Panel 6 uses the following codes:

- `-2`: empty cell;
- `-1`: poor cell;
- `0`: nominal cell;
- `1`: overpopulated cell;
- `2`: cell that received inserted particles at the current step.

## Debug figure

The optional debug figure adds:

- mean particle mass per cell;
- relative standard deviation of particle masses per cell;
- current insertion count per cell;
- `|u-u_mem|` where memory is available;
- active particle mass histogram;
- pool counters.

This view is intended to detect progressive pool consumption, persistent local
insertion regions and pathological mass dispersion.
