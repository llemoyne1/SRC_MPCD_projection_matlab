# Weighted resampling: wet/non-wet cell layer

This patch introduces a first cell-level wet/non-wet tag for the weighted
resampling branch.

Motivation
----------

The earlier weighted formulation imposed a full target mass in every grid cell:

```text
M_target(c) = gamma * m0
```

This is correct only for a fully fluid domain.  It is not compatible with future
immersed solids, mobile solids, injection/filling fronts or any variable active
fluid region.  Non-wet cells must not be treated as fluid cells by the local mass
remap, otherwise they are artificially filled.

New fields and helpers
----------------------

A cell mask is represented as an `Nx`-by-`Ny` logical grid:

```matlab
state.cellWetMask
params.cellWetMask
```

`true` means the cell belongs to the active fluid domain for population and mass
remapping.  `false` means the cell is currently non-wet/non-fluid for the mass
management layer.

New helper functions:

```matlab
resamp_cell_wet_mask.m
resamp_update_cell_wet_mask.m
```

Default behavior is all-wet, so existing TG/pocket/heterogeneity tests are
unchanged unless a wet mask is explicitly supplied or updated.

Mass-management behavior
------------------------

The following routines now honor the wet mask:

```matlab
resamp_insert_underpopulated_particles.m
resamp_extract_overpopulated_particles.m
resamp_local_mass_moment_remap.m
resamp_population_mass_diagnostics.m
resamp_deposit_weighted_to_grid.m
resamp_pool_visualize_frame.m
```

The rule is:

```text
wet cell     -> normal NMin/NMax/NTarget and M_target remap
non-wet cell -> no insertion, no extraction, no local M_target constraint
```

In `resamp_local_mass_moment_remap`, dry cells are skipped.  Their target mass is
not interpreted as `gamma*m0`; they are excluded from mass residual statistics.
This is the essential safeguard for future solid or empty regions.

Mask evolution
--------------

The mask can be fixed or updated through:

```matlab
[state, wetInfo] = resamp_update_cell_wet_mask(state, params, ...
    'mode','mass_hysteresis', ...
    'wetMassOnThreshold',0.5*params.resampTargetCellMass, ...
    'wetMassOffThreshold',0.05*params.resampTargetCellMass);
```

Supported update modes include:

```text
none / fixed            keep the current mask, or initialize all-wet
all                     force all cells wet
manual                  use an explicit cellWetMask
mass_threshold          wet if M >= threshold
mass_hysteresis         wet/dry with separate on/off mass thresholds
particle_threshold      wet if N >= threshold
particle_hysteresis     wet/dry with separate on/off population thresholds
```

Runner options
--------------

`run_resamp_pool_insertion_visual_demo` now accepts:

```matlab
'InitialWetMaskMode','all'      % default, all cells wet
'WetMaskUpdateMode','none'      % default, fixed mask
'WetMaskUpdateEvery',0          % 0 = no dynamic update
'CellWetMask',mask              % explicit Nx-by-Ny logical mask
'WetMassOnThreshold',value
'WetMassOffThreshold',value
'WetParticleOnThreshold',value
'WetParticleOffThreshold',value
```

A future injection/filling prototype can therefore start from dry cells and turn
them wet when deposited mass or particles exceed a threshold.

Important limitation
--------------------

This patch only protects the *mass/population management layer*.  The Q6
projection remains the current periodic Q6 operator and is not yet a masked
fluid-domain projection.  For immersed/mobile solids, the next layer will have to
combine this wet-cell mask with a proper projection/collision boundary treatment.
