Nous poursuivons le développement MATLAB SRC/MPCD/Q6/Q9 pour fluide quasi-incompressible.

Règle impérative : toute modification de code doit partir uniquement des fichiers que je fournis dans CE nouveau chat, pas de versions mémorisées ou de patchs précédents. Si un fichier manque, il faut me le demander ou me dire de le récupérer depuis le snap historique.

Contexte actuel :
- Les validations Taylor–Green et Poiseuille ont été faites avec le modèle Q9 complet : SRC/MPCD classic → projection Q6 de vitesse → correction Q9 low-k du flux de masse → thermostat corrigé, avec wallVP-v2 pour Poiseuille.
- Le piston ne doit pas valider une variante partielle q9_mass_only ; il doit tester le même modèle Q6/Q9 complet que TG/Poiseuille.
- Le mode q9_mass_only peut rester comme diagnostic de sensibilité, mais pas comme validation principale.
- Classic et Q6 piston passent déjà correctement sur 32x32, gamma=20, compression 5 %, wallVP-v2, thermostat corrigé.
- Résultats observés :
  classic final stdN ≈ 4.60, final lowK ≈ 1.335e-6, final Pkin ≈ 226.97, rhoRelError=0, kBT=0.01.
  Q6 final stdN ≈ 3.52, final lowK ≈ 3.05e-8, final Pkin ≈ 226.95, rhoRelError=0, kBT=0.01.
- Le Q6 homogénéise donc déjà fortement la densité tout en conservant la compression globale.

Problèmes actuels :
1. Le run Q9_FULL avec `massFluxProjectionOperator='general_bc'` rencontre un problème de matrice mal conditionnée dans `projection_project_mass_flux_general_bc`, probablement dû au mode constant nul d’un problème elliptique périodique/Neumann. Il faut corriger cela proprement par un gauge-fix du potentiel, sans changer la physique.
Le warning n’est pas causé par Q9 lui-même.
Il est causé par le solveur elliptique general_bc appliqué à un problème périodique/Neumann dont le potentiel possède un mode constant nul.
La bonne correction est un gauge-fix dans general_bc, sans changer le cœur physique Q6/Q9 ni spécialiser le piston.
2. Il faut aussi vérifier que le diagnostic de pression wallVP piston n’a pas régressé : dans un run level3 précédent, la pression top wall totale était dominée par wallVP et cohérente avec Pkin ; dans le run full-model récent, la contribution wallVP ressortait à 0, ce qui suggère soit un fichier shadowed, soit une régression du diagnostic.

Objectif immédiat :
- Reprendre proprement le modèle piston Q6/Q9 complet à partir des fichiers fournis.
- Tester d’abord `q9_historical` avec `periodic_x_neumann_y` + `lowpass_fft`, pour valider le modèle Q6/Q9 complet sans être bloqué par `general_bc`.
- Ensuite seulement corriger `general_bc` par gauge-fix et comparer `q9_full` elliptique à `q9_historical`.

Les résultats de ces deux runs sont joints

Je vais joindre aussi :
- le snap courant du répertoire piston,
- le snap historique Q6/Q9 contentant les tests piston avant la détection du bug thermostat, l'opérateur fft historiques et la parallélisation partielle du code
- le snap des scripts sans le bug thermostat, semi-parallélisés et avec opérateur elliptique 
Commence par inspecter les fichiers fournis, identifier précisément les fonctions manquantes et les différences de chemin modèle, puis analyse les résultats des trois versions classic/q6/q9 pour les deux opérateurs ans le cas piston.


