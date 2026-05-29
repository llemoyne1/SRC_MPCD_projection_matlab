# Latent particles and injection/fill prototype

This patch adds a particle-level role layer on top of the wet-cell mask.

Particle roles:

- `0`: inactive/free pool slot;
- `1`: fluid particle, included in deposit, collision, Q6, remap and thermostat;
- `2`: latent marker, stored in the particle arrays but hydrodynamically inert.

The helper `resamp_active_mask` now returns only role-1 fluid particles when
`state.particleRole` exists.  This keeps legacy states unchanged while making
latent states safe by default.

New helpers:

- `resamp_particle_role_mask.m` summarizes role masks;
- `resamp_set_dry_cell_particles_latent.m` converts particles in dry cells to
  latent markers (`active=false`, `particleRole=2`, `m=0`, `v=0`);
- `resamp_activate_particles_in_cell.m` activates latent/free slots as fluid in
  selected cells;
- `run_resamp_injection_fill_demo.m` starts from a fully latent domain and
  injects fluid support at a source cell.

Minimal smoke:

```matlab
outFill = run_resamp_injection_fill_demo( ...
    'method','weighted_classic', ...
    'steps',100, ...
    'visualEvery',5, ...
    'summaryEvery',5, ...
    'Nx',64, 'Ny',32, 'gamma',20, ...
    'inletVelocityX',0.15, ...
    'injectionParticlesPerCell',20, ...
    'NMin',10, 'NTarget',20, 'NMax',30, ...
    'ThermostatAfterRemap',true, ...
    'ThermostatStrength',0.25, ...
    'rngSeed',12345);
```

This is not yet a physical inlet/outlet validation.  It is a stress test of the
new active-domain layer: dry cells are not remapped to full mass, and latent
particles do not create pressure, collision or projection contributions.
