# Hotfix latent injection source-cell selection

Fixes `run_resamp_injection_fill_demo` source-cell construction.  The previous
periodic index expression could return an empty source list when
`injectionPatchSize=[1 1]`, so no latent particles were activated and the smoke
remained fully dry.

Also writes the time column as `step*dt` in the injection-fill CSV.

Expected smoke after this hotfix: `activatedParticles > 0`, `NpActive > 0`, and
`nWetCells > 0` from the first injection step.
