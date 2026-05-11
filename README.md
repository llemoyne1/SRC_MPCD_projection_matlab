# SRC\_MPCD\_projection\_matlab

Prototype MATLAB pour le développement d’une méthode SRC/MPCD incompressible par projection de pression.

Ce dépôt part d’une version MATLAB existante du solveur SRC/MPCD incompressible fondée sur :

* collision SRC/MPCD classique ;
* redistribution particulaire pour contrôler l’occupation cellulaire ;
* fermeture liquide par réparation de vitesse et pression motrice ;
* diagnostics Poiseuille, diffusion, onde de densité, compression et piston.

L’objectif de cette nouvelle branche de travail est d’évaluer une formulation alternative de l’incompressibilité fondée sur la contrainte :

```math
\\nabla \\cdot u = 0



@'



\## Q7c reference: weak virial position repair with velocity restoration



This milestone extends the Q6 pressure-projection method with a weak liquid-like density repair step.



The validated Q7c reference configuration is:



```matlab

params.useVirialDensityRepair = true;

params.virialDensityRepairStrength = 0.001;

params.virialDensityRepairK = 1.0;

params.virialDensityRepairSmoothPasses = 1;

params.virialDensityRepairMaxDisplacementFraction = 0.01;

params.virialDensityRepairRestoreVelocity = true;

