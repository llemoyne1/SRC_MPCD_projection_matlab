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
