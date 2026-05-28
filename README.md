# SRC/MPCD weighted resampling prototype

Branche MATLAB minimale dédiée à l'exploration d'une formulation SRC/MPCD pondérée avec resampling local conservatif.

Objectif :
- maintenir un support particulaire local suffisant ;
- autoriser des masses particulaires variables mais proches de m0 ;
- imposer localement masse et moment ;
- utiliser Q6 comme projection incompressible principale ;
- exclure provisoirement Q9, viriel, piston, cylindre et parois complexes.

État initial :
- noyau périodique classic/Q6 conservé comme référence ;
- nouvelles fonctions resamp_* à ajouter progressivement ;
- pas de dépendance aux anciens scripts de redistribution/liquid closure.
