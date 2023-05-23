# Bureau d'Etude : Résolution numérique d’une équation aux dérivées partielles d’advection-diffusion par une méthode de type Volumes Finis

Pour visualiser la mise en page de ce Readme.md : https://github.com/GregHsr/BE_Volumes_Finis/blob/master/README.md

L'objectif de ce projet est de simuler numériquement le mélange de deux composé chimiques à l'aide de la méthode des différences finies. Nous nous trouvons ici dans le cas d'une résolution numérique d’équations aux dérivées partielles d’advection-diffusion sur __Fortran90__.

Ce projet est composé d'un code principal **prog.f90** faisant appel aux différentes fonctions définies dans les fichiers annexes. 
Il est aussi composé d'un code python permettant de récupérer les données calculées et de les traiter différemment. 

#### Liste des fichiers contenus dans ce répertoire :
- prog.f90 : code principal
- subroutines.f90 : code contenant les subroutines de calcul
- verif.f90 : code permettant d'effectuer les vérifications des calculs effectués
- m_type.f90 : définition du module utilisé
- VTSWriter.f90 : subroutine permettant d'enregistrer les données
- Makefile : fichier de compilation
- analyse.py : code d'analyse python
- data.txt : paramètres du problèmes * à modifier *


#### Comment faire tourner ce code (sous Linux)

- Placer tous les fichiers listés ci-dessus dans un même répertoire (en dézippant le dossier rendu, ou à l'aide de "git clone" : ``` git clone https://github.com/GregHsr/BE_Volumes_Finis ``` pour télécharger tous les fichiers dans le répertoire courant).
- Compiler le programme : dans un terminale, écrire ``` make ```
- Executer le programme : toujours dans le terminale, écrire ```./prog.exe```
- Finalement, des fichiers du type __sol_00*__ et un fichier __sol.pvd__ ont du être générés et peuvent-être ouverts à l'aide de __Paraview__ ou du __code Python__.

