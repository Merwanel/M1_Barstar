#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@purpose: SUJET DYNAMIQUE MOLECULAIRE ETUDE DE LA BARSTAR
@author: EL ASRI Merwan
"""

import matplotlib.pyplot as plt
from fun.PDBTools import *
import numpy as np

if "-help" in sys.argv:
    usage()

try :
    fileIndex = sys.argv.index("-file") + 1
    fileName = sys.argv[fileIndex]
except:
    print("\nLe nom du fichier pdb n'a pas été trouvé\n\n")
    usage()
            
try : 
    seuilIndex = sys.argv.index("-seuil") + 1
    seuil = sys.argv[seuilIndex] 
except:
    seuil = 3
    
try : 
    resIndex = sys.argv.index("-res") + 1
    res_pour_analyse_locale = sys.argv[resIndex].split(",")
except :
    res_pour_analyse_locale = []

try : 
    cvIndex = sys.argv.index("-rc") + 1
    rc = sys.argv[cvIndex]
except :
    rc = 8

try : 
    sys.argv.index("-c")
    croise = True
except :
    croise = False

try:
	sys.argv.index("-w")
	write = True
except:
	write = False

try:
    TabPDB = parsePDBMultiChains(fileName)
except:
    print("\nLe nom du fichier pdb n'a pas été trouvé\n\n")
    usage()







tab_rmsd,tab_rmsd_moy,tab_rmsd_point,rmsd_moyen_par_res , nb_res = rmsdExtractor(TabPDB)
nModel = list(range(0, len(tab_rmsd_moy)))
nb_res = len(rmsd_moyen_par_res) # nombre de résidu


figuerRmsdMoy(tab_rmsd_moy )
figureRmsdMoyETParRes(tab_rmsd ,tab_rmsd_point ,  tab_rmsd_moy , seuil , nb_res )


#enfouissement

enfouissementEuclidien ,tab_enfouissement_point , enfouissement_moyen_par_res = enfouissementExtractor(TabPDB)

figureEnfouissementRmsdMoyParRes(rmsd_moyen_par_res , enfouissement_moyen_par_res , res_pour_analyse_locale , len(tab_rmsd_moy))
figureEnfouissementRmsdToutLesPoints(tab_rmsd_point, tab_enfouissement_point , len(tab_rmsd_moy))

#rayon de giration
giration = giration(enfouissementEuclidien)
plotGiration(giration)


##graph locaux

for res in res_pour_analyse_locale:
    makeFigurePourUnResidu( res , tab_rmsd, tab_rmsd_moy  , 'rmsd')
    makeFigurePourUnResidu( res , enfouissementEuclidien, giration  , 'distance euclidienne du centre')
    FigureCV_res(TabPDB , res , rc)


#### Ecriture de fichiers
if write == True :
	writeRmsdMoy_Giration(tab_rmsd_moy ,giration )
	writeRmsdEnfouissementMoyen_ParRes(rmsd_moyen_par_res ,enfouissement_moyen_par_res, len(tab_rmsd))
	writeRmsdGiration_ParRes_ParModel( tab_rmsd , enfouissementEuclidien ) 


### croisé :
if croise == True:
    rmsdCroise = rmsdCroisee(TabPDB)
    plotRmsdCroisee(rmsdCroise)
