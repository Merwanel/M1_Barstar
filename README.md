# Dynamique moléculaire; Étude de la barstar

Le script 'projet.py' permet de créer des figures et fichiers utiles à l'analyse de la
dynamique moléculaire de la barstar ou de toute autre protéine.


## Prérequis

	python3 avec les packages suivants installés: matplotlib , numpy , math et string

	Le fichier 'PDBTools.py' doit ȇtre placé dans un dossier 'fun'

	
## PDBTools.py

Le fichier 'PDBTools.py' vient des travaux pratiques effectués au cours de l'année.
Le fichier 'PDBTools.py' contient toutes les fonctions utilisés dans le script 'projet.py'

Les fonctions qui ont été crée et ajouté à 'PDBTools.py' sont :  
rmsd, rmsdExtractor, figuerRmsdMoy, figureRmsdMoyETParRes_rmsd, rmsdCroisee, plotRmsdCroisee, computesProteinCenter, enfouissementExtractor,
giration, plotGiration, figureEnfouissementRmsdToutLesPoints, figureEnfouissementRmsdMoyParRes, 
makeFigurePourUnResidu, Rij, Norme, getXYZ, Env_i, sumVect, CVi, FigureCV_res,
writeRmsdMoy_Giration, writeRmsdEnfouissementMoyen_ParRes, writeRmsdGiration_ParRes_ParModel et usage

La fonction parsePDBMultiChains a également été modifié de manière à pouvoir récupérer plusieurs structures dans le mȇme fichier.


## Usage


python3 projet.py -file <FILE_NAME> -res <numéros des résidus à analyser> -seuil 3 -rc 8 -c -w





-file:	nom du fichier pdb. C'est le seul paramètre obligatoire. Son absence provoque l'arrȇt du script.

-res:	les noms des résidus doivent être séparés par une virgule, sans espace.

-seuil:	pour définir le seuil, les résisdus dont le rmsd > rmsd_moyen + seuil   sont mis en évidence. Par défaut seuil = 3 Å

-rc:	seuil en Å, à partir duquel la variance circulaire est calculée. Par défaut rc = 8 Å

-c:		à écrire si l'utilisateur veut obtenir un plot représentant le rmsd entre chaque modèle.
		ATTENTION : l'algorithme pour calculer la matrice à l'origine de cet figure est 
                    d'une compléxité n², avec n le nombre de modèle, si il y a 2000 modèles, 
                    le calcul peut durer plus d'une heure.
                         
-w:		à entrer si l'utilisateur veut écrire les fichiers textes.
		
		
		Entrez : <NOM_DU_SCRIPT_PYTHON> -file md_prot_only_skip100.pdb -res 76,80,39,35 -c -w && <NOM_DU_SCRIPT_PYTHON> -file md_prot_only_skip10.pdb -res 76,80,39,35 -w -c
				pour créer les mêmes fichiers que dans le rapport

#### Exemples d'utilisation

 
Dans un terminal:	

python3 projet.py -file md_prot_only_skip100.pdb          					crée les figures 2, 6, 8 et 10
python3 projet.py -file md_prot_only_skip100.pdb  -res 76,35,39,80			crée les figures 2, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34 et 36
python3 projet.py -file md_prot_only_skip100.pdb	-c							crée les figures 2, 4, 6, 8 et 10
python3 projet.py -file md_prot_only_skip100.pdb	-w							crée les figures 2, 6, 8, 10 et les 3 fichiers textes



#### Commentaires sur les calculs d'enfouissement

Pour l'analyse globale, les calculs d'enfouissement ont tous été approximé par la moyenne des distances euclidiennes
entre un résidu et le centre de la protéine.
La variance circulaire prend en compte toutes les distances des atomes proches, le temps de calcul est important,
donc la variance circulaire n'est pas calculée pour l'analyse globale.

Pour l'analyse des résidus demandée avec la commande '-res', la variance circulaire ET la distance euclidienne
entre un résidu et le centre de la protéine sont calculées.
La variance circulaire est calculée pour le barycentre du résidu et en prenant en compte tous les atomes autours situés à moins
d'un certain seuil fixé par la commande '-rc' (par défaut rc = 8 Å )


## Auteur

* EL ASRI Merwan

