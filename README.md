Notre premier ensemble d'objectifs de référence se compose de six tâches individuelles inspirées de la conception de photovoltaïques organiques (OPV). Le développement de cellules solaires organiques (OSC) présente un grand intérêt en raison de leur potentiel à remplacer les dispositifs inorganiques actuellement prédominants et à élargir leur domaine d'application.

Ce projet vise à calculer le PCE (Pouvoir de Conversion de Puissance) des cellules photovoltaïques organiques (OPV) en utilisant différentes méthodes de calcul combinées.

*  Notre objectif de recherche consiste à concevoir des petites molécules dotées de propriétés électroniques spécifiques, notamment des molécules capables d'effectuer la séparation des charges, inspirées de la conception photovoltaïque organique (OPV). Nous avons deux tâches individuelles :

- Recherche d'une molécule donneuse organique à utiliser avec l'ester méthylique de l'acide l'ester méthylique de
    l'acide [6,6]­phényl­C61­butyrique(PCBM).
- Recherche d'une molécule accepteuse à utiliser dans des dispositifs basés sur le $ poly[N­90­heptadécanyl­2,7­carbazole­alt­5,5­(40,70­di­2­thiényl­20,10,30­benzothiaMachine  (PCDTBT)).
Le code accepte une structure proposée sous forme de chaîne SMILES, génère des coordonnées cartésiennes initiales avec rdkit et effectue une recherche de conformateur et une optimisation de la géométrie avec crest et xtb, respectivement. Enfin, un calcul en un seul point au niveau théorique GFN2-xTB fournit les propriétés d'intérêt, notamment les énergies HOMO et LUMO, le gap HOMO-LUMO et le moment dipolaire moléculaire. L'efficacité de conversion de puissance (PCE) est calculée à partir de ces propriétés simulées sur la base du modèle Scharber.
Calcul du PCE des OPV


Méthodes utilisées
Calcul des niveaux HOMO, LUMO et du gap via RDKit, Crest et XTB :
Utilisation de RDKit pour la manipulation moléculaire
Optimisation de la géométrie avec Crest
Calcul des niveaux HOMO, LUMO et du gap avec XTB
Modèle de Scharber:
Mise en œuvre du modèle de Scharber pour estimer le PCE à partir des niveaux HOMO, LUMO et du gap
Calculs DFT, HF et MDA avec RDKit, Crest, XTB et PySCF :
Utilisation de RDKit, Crest et XTB pour la préparation des molécules
Calculs DFT, HF et MDA avec PySCF
Analyse et comparaison des résultats
Intégration avec l'algorithme de synthèse SAScores :
Utilisation de l'algorithme SAScores pour estimer la facilité de synthèse des composés
Combinaison des résultats de PCE et de synthèse pour une évaluation complète des molécules
Objectifs
L'objectif principal de ce projet est de développer une méthodologie complète pour évaluer le potentiel des matériaux OPV en termes de PCE et de facilité de synthèse. Les résultats obtenus pourront aider à identifier les composés les plus prometteurs pour le développement de cellules photovoltaïques organiques performantes.

Utilisation
Le code Python et les notebooks Jupyter associés permettent de réaliser les différentes étapes de calcul et d'analyse décrites ci-dessus. Les instructions d'installation et d'utilisation seront fournies dans le README du projet
