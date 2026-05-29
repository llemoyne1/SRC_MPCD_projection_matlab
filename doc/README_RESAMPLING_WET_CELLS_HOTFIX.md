# Wet-cell mask hotfix

This hotfix prevents an empty forwarded `cellWetMask` option from forcing
`resamp_cell_wet_mask` into explicit/manual mode.

Before this fix, initialization paths such as
`resamp_initialize_particles_taylor_green_forced -> resamp_deposit_weighted_to_grid`
passed an empty optional mask while the state did not yet contain
`state.cellWetMask`.  The empty mask was interpreted as an explicit mask and
raised:

```text
cellWetMask is empty.
```

The corrected behavior is:

- empty forwarded mask: keep auto resolution, then fall back to state/params/all-wet;
- non-empty forwarded mask: treat as explicit and validate dimensions.
