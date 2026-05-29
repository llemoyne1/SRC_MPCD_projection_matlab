# Remap mass safety: bounded masses + velocity shift

This note documents the safety fallback added to the weighted-resampling MATLAB prototype.

## Motivation

The closed-loop population mechanism

```text
overpopulated cells -> extraction -> pool -> insertion -> underpopulated cells -> remap
```

works algebraically, but the previous remap could satisfy local mass and momentum by creating very small or very large particle masses when the population edit was strong.  This is acceptable as a diagnostic of the algebra, but it is not acceptable as a robust SRC/MPCD model: transport properties should remain close to a quasi-constant-mass particle system whenever the population support is kept near `gamma`.

## New fallback

The normal remap is still attempted first.  If its candidate masses exceed a configurable safety band, the cell is remapped using:

```matlab
m_p = M_target / N_cell
v_p <- v_p + (U_target - U_current)
```

where `U_current` is computed after assigning the uniform cell mass.  This exactly preserves:

```text
sum_p m_p = M_target
sum_p m_p v_p = M_target U_target
```

while preserving the relative velocity fluctuations inside the cell.  The fallback is therefore a local bounded-mass safety mechanism, not a new population force.

## Options

The main runner forwards these options to `resamp_local_mass_moment_remap`:

```matlab
'RemapMassSafetyEnable',true
'RemapMassSafetyMode','uniform_mass_velocity_shift'
'RemapMassSafetyMinFactor',0.25
'RemapMassSafetyMaxFactor',4.0
```

The safety is conditional by default: it is applied only in cells where the normal remap would produce a candidate mass outside the safety band:

```text
m_p / m0 < RemapMassSafetyMinFactor
m_p / m0 > RemapMassSafetyMaxFactor
```

To disable the safety path explicitly:

```matlab
'RemapMassSafetyEnable',false
```

To apply the uniform-mass/velocity-shift path to all remapped non-empty cells:

```matlab
'RemapMassSafetyEnable',true, ...
'RemapMassSafetyMode','always_uniform_mass_velocity_shift'
```

The always mode is a diagnostic mode; it is not the recommended default.

## Diagnostics

The summary table now includes:

```text
remapMassSafetyCells
remapMassSafetyParticles
remapMassSafetyInfeasibleCells
remapMassSafetyVelocityShiftRms
remapMassSafetyVelocityShiftMax
remapMassSafetyCandidateMinFactor
remapMassSafetyCandidateMaxFactor
```

The figure title displays the number of safety cells and the RMS velocity shift for the current visualized step.

## Suggested check

A distributed heterogeneity sweep can be rerun with the same command as before.  The expected improvement is not necessarily a lower `mParticleRelStd` for mild cases, because the safety should remain inactive there.  The expected improvement is the disappearance of extreme mass tails in strong cases, especially in `mParticleMinFinal` and `mParticleMaxFinal`.

