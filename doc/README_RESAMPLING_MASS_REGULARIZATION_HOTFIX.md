# Hotfix: activeMask support in cell id helper

This hotfix makes `resamp_cell_ids_periodic` accept the `activeMask` option used by
`resamp_regularize_particle_masses_preserve_mue`.

The regularization patch calls the cell-id helper on the full pooled particle
arrays and passes `activeMask` so that inactive/latent particles are ignored.  The
previous helper only accepted `periodicX` and `periodicY`, which caused:

```text
Unknown option: activemask
```

The new behavior is backward compatible: if `activeMask` is omitted, all supplied
positions are mapped exactly as before.  If it is provided, inactive entries are
assigned cell id `0` by default and active entries are mapped normally.
