# SRC_MPCD_projection_matlab

Prototype MATLAB pour le développement d’une méthode SRC/MPCD incompressible par projection de pression.

Ce dépôt part d’une version MATLAB existante du solveur SRC/MPCD incompressible fondée sur :
- collision SRC/MPCD classique ;
- redistribution particulaire pour contrôler l’occupation cellulaire ;
- fermeture liquide par réparation de vitesse et pression motrice ;
- diagnostics Poiseuille, diffusion, onde de densité, compression et piston.

L’objectif de cette nouvelle branche de travail est d’évaluer une formulation alternative de l’incompressibilité fondée sur la contrainte :

```math
\nabla \cdot u = 0

## Q9 reference: velocity projection with low-k mass-flux density relaxation

Q9 is the current reference hybrid incompressible SRC/MPCD method. It combines the velocity projection used in Q6 with a low-wavenumber mass-flux correction inspired by Q8.

The objective is to preserve the robust Poiseuille behavior of the velocity projection while further reducing large-scale density modes without acting directly on the cell-scale particle-number noise.

### Method

Starting from a standard SRC/MPCD step, Q9 applies:

```text
SRC/MPCD step
-> velocity projection: div(u) ≈ 0
-> density defect: N - gamma
-> low-k filtering of the density defect
-> mass-flux relaxation on the filtered defect
-> thermostat

Q9 preserves the local density quality of Q6 while reducing the low-k density energy by a factor of about 4.3.