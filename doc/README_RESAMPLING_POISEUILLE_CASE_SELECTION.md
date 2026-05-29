# Poiseuille wallVP case selection

This hotfix adds a `cases` option to `run_resamp_poiseuille_wallvp_validation`.

Supported values:

- `{'classic_wallvp_reference'}` or alias `{'src_classic_wallvp_only'}`
- `{'q6_resampled_wallvp'}`
- `{'all'}` or `{'both'}`

The alias `src_classic_wallvp_only` maps to the existing pure SRC baseline:
classic weighted wallVP, body force, no Q6 projection, no insertion/extraction,
no remap, and no mass safety.
