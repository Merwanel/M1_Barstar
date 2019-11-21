#!/usr/bin/env python3

import string, sys, math
import numpy as np
import matplotlib.colors as colors
import matplotlib.pyplot as plt


######################################
#
#       Parsers & Writers
#
########################################

def parsePDBMultiChains(infile) :
    """ propos : parser un fichier PDB
        input : un fichier PDB pouvant contenir plusieurs structures
                séparées par "ENDMDL"
        output : un tableau de dicos de type dPDB;
                le i-ème élément du tableau correspond à la i-ème structure
                présent dans le fichier PDB
    """


    # lecture du fichier PDB
    #-----------------------
    f = open(infile, "r")
    lines = f.readlines()
    f.close()


    TabPDB = [] #initialisation du tableau de dictionnaires
    i = -1 # initialisation du numéro de la ligne

    while( i < len(lines) ) :
        i+= 1

        # var ini
        #---------
        chaine = True
        firstline = True
        prevres = None
        dPDB = {}
        #dPDB["reslist"] = []
        dPDB["chains"] = []


        # parcoure le PDB
        #-----------------

        # A chaque fois que "ENDMDL est rencontré, cela tarduit la fin de la
        # structure
        while( (i < len(lines)) and ( "ENDMDL" not in lines[i]) ) :
            #print("i= ", i , "   ",lines[i][0:3])
            line = lines[i]
            i+= 1
            #print(line[0:3] == "TER", " len(lines)= ",len(lines))
            if line[0:4] == "ATOM" :

                # on recupere l'info de la chaine
                chain = line[21]

                # si la chaine n'existe pas, on cree la cle correspondante et on ajoute la chaine a la liste des chaines
                if not chain in dPDB["chains"] :
                    dPDB["chains"].append(chain) # ajout de "chain" a la liste des chaines
                    dPDB[chain] = {} # creation du sous-dico pour la chaine
                    # on prepare la structure de donnees pour cette chaine
                    dPDB[chain]["reslist"] = []

                # on recupere l'info du residu
                curres = "%s"%(line[22:26]).strip()

                # si le residu pour cette chaine "chain" n'existe pas, on cree la cle correspondante et on ajoute le res a la liste des res
                if not curres in dPDB[chain]["reslist"] :
                    dPDB[chain]["reslist"].append(curres)
                    dPDB[chain][curres] = {}
                    # on prepare la structure de donnees pour ce residu
                    dPDB[chain][curres]["atomlist"] = []
                    # on recupere l'info du residu
                    dPDB[chain][curres]["resname"] = line[17:20].strip()

                # on recupere les info pour l'atome de ce res de cette chaine (type atomique + coords x, y, z)
                atomtype = line[12:16].strip()
                dPDB[chain][curres]["atomlist"].append(atomtype)
                dPDB[chain][curres][atomtype] = {}
                dPDB[chain][curres][atomtype]["x"] = float(line[30:38])
                dPDB[chain][curres][atomtype]["y"] = float(line[38:46])
                dPDB[chain][curres][atomtype]["z"] = float(line[46:54])
                dPDB[chain][curres][atomtype]["id"] = line[6:11].strip()
                dPDB[chain][curres][atomtype]["bfactor"] = float(line[60:67].strip())

        # Si la une structure a fini d'être parsé, elle est alors ajoutée au tableau
        if( len(dPDB["chains"]) > 0 ) :
            TabPDB.append(dPDB)

    return TabPDB

#################################################
#           WRITING TOOLS
#################################################


def writePDB(dPDB, filout = "out.pdb", bfactor = False) :
    """purpose: according to the coordinates in dPDB, writes the corresponding PDB file.
       If bfactor = True, writes also the information corresponding to the key bfactor
       of each atom in dPDB.
       input: a dico with the dPDB format
       output: PDB file.
    """

    fout = open(filout, "w")

    for chain in dPDB["chains"]:
        for res in dPDB[chain]["reslist"] :
            for atom in dPDB[chain][res]["atomlist"] :
                if bfactor :
                    fout.write("ATOM  %5s  %-4s%3s %s%4s    %8.3f%8.3f%8.3f  1.00%7.3f X X\n"%(dPDB[chain][res][atom]["id"], atom, dPDB[chain][res]["resname"],chain, res,dPDB[chain][res][atom]["x"], dPDB[chain][res][atom]["y"],dPDB[chain][res][atom]["z"],dPDB[chain][res][atom]["bfactor"] ))
                else:
                    fout.write("ATOM  %5s  %-4s%3s %s%4s    %8.3f%8.3f%8.3f  1.00  1.00 X X\n"%(dPDB[chain][res][atom]["id"], atom, dPDB[chain][res]["resname"],chain, res,dPDB[chain][res][atom]["x"], dPDB[chain][res][atom]["y"],dPDB[chain][res][atom]["z"] ))

    fout.close()


######################################
#
#           extract info from PDB
#
########################################



# ceci est une possibilite, vous pouviez aussi renvoyer la sortie avec une liste de listes
# de type [[chainname, nb_aa1],[chainname, nb_aa2],etc] --> [["A", 23],["B", 29],["C", 11]]
def getNbAA(dPDB, aa) :
    """ Purpose: counts the number of aa per chain of a protein
        Input: dPDB (a dico PDB) and the type of aa (string) in 3letter code (i.e. ALA, CYS...)
        Output: a dico containing for each chain, the nb of residues corresponding to "aa"
    """

    # init var
    d_aaPerChain = {}

    for chaini in dPDB["chains"] :
        # initiation de la cle chaini. La cle chaini est pointe vers un compteur initialise a zero
        d_aaPerChain[chaini] = 0

        for resi in dPDB[chaini]["reslist"] :
            # compte nb aa pr la chaine i
            if dPDB[chaini][resi]["resname"] == aa :
                d_aaPerChain[chaini]+= 1

    return d_aaPerChain

######################################
#
#           3D manipulations
#
########################################


def distancePoints(L1,L2):
    """Computes the distance between the two sets of coordinates
       input: 2 tuples with the corresponding coordinates
       output: distance"""

    x = L1[0]-L2[0]
    y = L1[1]-L2[1]
    z = L1[2]-L2[2]
    return math.sqrt(x*x+y*y+z*z)
    

def getCentroid(d_res):
    """Purpose: Calculates the center of mass of a residue
       Input: a dico residue
       Output: coords x, y, z of the centroid (tuple format)
    """

    x = y = z = 0.0

    # loop over all atoms
    for atom in d_res["atomlist"] :
        x +=d_res[atom]["x"]
        y +=d_res[atom]["y"]
        z +=d_res[atom]["z"]

    Xcen = float(x)/len(d_res["atomlist"])
    Ycen = float(y)/len(d_res["atomlist"])
    Zcen = float(z)/len(d_res["atomlist"])

    return (Xcen, Ycen, Zcen)




def compDistance(d_res1, d_res2, mode = "minvalue") :
    """Purpose: computes the distance between 2 residues (res1 and res2)
                Distance can be calculated either as the min distance between all the pairwise distances (atom-atom)
                or between the 2 centroids of the residues

       Inputs:  d_res1, d_res2 which are dico corresponding to res 1 and res 2 respectively
                mode: either "minvalue" or "centroid" (string)

       Output:  distance (float) """

    # mode minvalue: min distance atom-atom
    #--------------------------------------
    if mode == "minvalue" :
        distance = 1000000

        #*** computes all pairwise distances between atoms of res1, res2

        # loop over all atoms from res1
        for atom1 in d_res1["atomlist"] :
            coord1 = [d_res1[atom1]["x"], d_res1[atom1]["y"], d_res1[atom1]["z"]]
            # lloop over all atoms from res2
            for atom2 in d_res2["atomlist"] :
                coord2 = [d_res2[atom2]["x"], d_res2[atom2]["y"], d_res2[atom2]["z"]]
                # computes distance between atom1 from res1 and atom2 from res2
                dist_tmp = distancePoints((coord1[0], coord1[1], coord1[2]),(coord2[0],coord2[1], coord2[2]))
                if distance > dist_tmp :
                    distance = dist_tmp

    # mode centroid: distance between the centroids of the 2 given residues
    #----------------------------------------------------------------------
    elif mode == "centroid" :
        cent1 = getCentroid(d_res1)
        cent2 = getCentroid(d_res2)
        #print(cent1, cent2)

        distance = distancePoints(cent1, cent2)

    return distance


def rmsd( L1,L2 ) :
    """purpose : calcul le rmsd entre deux points
        input : 2 vecteur contenant trois dimensions
        output : rmsd
    """
    return math.sqrt(distancePoints(L1,L2)**2)


def rmsdExtractor(TabPDB) :
    """
    propos : calculer le rmsd et créer plusieurs tableau et dictionnaire qui vont 
    servir pour les graphiques à plotter
    input : un tableau de dicos de type dPDB;
            le i-ème élément du tableau correspond à la i-ème structure
            présent dans le fichier PDB 
    output :
        tab_rmsd = un tableau de dictionnaire, dont le i-ème élément est le
                   dictionnaire correspondant au i-ème modèle. Chaque dictionnaire contient
                   le rmsd moyen du modèle et les rmsd par résidu
        tab_rmsd_moy = tableau dont le i-ème élément est le rmsd moyen du i-ème modèle
        tab_rmsd_point = tableau contenant tous les rmsd de tous les résidus à la suite
        rmsd_moyen_par_res = dictionnaire contenant pour chaque résidu son rmsd moyen
        nb_res = nombre de résidu dans une structure
    """
    dPDB_t0 = TabPDB[0] #première structure du fichier, utilisée comme référence
    
    tab_rmsd = []
    tab_rmsd_moy = []
    tab_rmsd_point = [] 
    rmsd_moyen_par_res = {}
    for dPDB in TabPDB  :
    
        #dico_rmsd est un dictionnaire qui contient le rmsd moyen du modèle qui 
        # est en train d'être lu et les rmsd par résidu
        dico_rmsd = {}
        dico_rmsd["rmsd_moyen"] = 0
        nb_res = 0 # nombre de résidu dans la structure
        for chain in dPDB["chains"] :
            for res in dPDB[chain]["reslist"] :
                center = getCentroid( dPDB[chain][res] )# barycentre du résidu sur le modèle qui est en train d'être lu
                center_ref = getCentroid( dPDB_t0[chain][res] )# barycentre du résidu sur le modèle de référence
                dico_rmsd["rmsd_moyen"] += rmsd(center , center_ref)
                dico_rmsd[res] = rmsd(center , center_ref)
    
                tab_rmsd_point.append(rmsd(center , center_ref))
                nb_res += 1

                if(res not in rmsd_moyen_par_res):
                    rmsd_moyen_par_res[res] = rmsd(center , center_ref)
                else:
                    rmsd_moyen_par_res[res] += rmsd(center , center_ref)
            
        dico_rmsd["rmsd_moyen"] /= nb_res
        tab_rmsd_moy.append( dico_rmsd["rmsd_moyen"] )
        tab_rmsd.append(dico_rmsd)
    
    for res in rmsd_moyen_par_res.keys() :
        rmsd_moyen_par_res[res] /= (len(tab_rmsd_moy)-1)

    return tab_rmsd,tab_rmsd_moy,tab_rmsd_point,rmsd_moyen_par_res , nb_res



def figuerRmsdMoy(tab_rmsd_moy ):
    """
    propos : plot de rmsd
    input :
        tab_rmsd = un tableau de dictionnaire, dont le i-ème élément est le
               dictionnaire correspondant au i-ème modèle. Chaque dictionnaire contient
               le rmsd moyen du modèle et les rmsd par résidu
    output :
            un plot représentant le rmsd moyen, dans le format .png
    """
    nModel = list(range(0, len(tab_rmsd_moy)))
    plt.figure(figsize=(20, 16))
    plt.plot(nModel, tab_rmsd_moy, ms=3.5, color="black", marker='o', linestyle='-', label='rmsd moyen')
    plt.xlabel('model', fontsize=18)
    plt.ylabel('rmsd moyen', fontsize=18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.savefig("rmsd_moyen_par_model." + str(len(tab_rmsd_moy)) + ".model.png")



    
def figureRmsdMoyETParRes(tab_rmsd ,tab_rmsd_point ,  tab_rmsd_moy , seuil , nb_res) :
    """
    propos : Plot de rmsd
    input : 
            tab_rmsd = un tableau de dictionnaire, dont le i-ème élément est le
                   dictionnaire correspondant au i-ème modèle. Chaque dictionnaire contient
                   le rmsd moyen du modèle et les rmsd par résidu
            tab_rmsd_moy = tableau dont le i-ème élément est le rmsd moyen du i-ème modèle
            tab_rmsd_point = tableau contenant tous les rmsd de tous les résidus à la suite
            seuil = pour un modèle et un résidu donnés :
                    si seuil + rmsd moyen > rmsd du residu, alors ce résidu est 
                    marqué d'un diamant coloré
            nb_res = nombre de résidu dans une structure
    output :
            un plot représentant le rmsd moyen , les rmsd par résidu et met en évidence 
            les rmsd au-dessus d'un seuil prédéfini, dans le format .png
    """
    nModel = list(range(0, len(tab_rmsd_moy))) # absicsi    
    x = [] # abscisse pour le plot des rmsd par résidu
    for i in range( len(tab_rmsd) ) :
        x += [i] * nb_res
    plt.figure(figsize=(20,20))
    plt.plot(x,tab_rmsd_point, "go",ms = 1 , alpha = 0.1 , label='rmsd par AA')
    plt.plot( nModel , tab_rmsd_moy ,ms= 3.5 , color= "black" , marker='o' , linestyle ='-',label='rmsd moyen' )
    colors_list = list(colors._colors_full_map.values()) #liste de couleur disponible 
    nColor = 10 # variable d'iteration sur la liste des couleurs (les couleurs avant 10 sont moins joli)
    handler = []
    AA = {} # AA est un dictionnaire spécifiant pour chaque résidu sa couleur attribuée
    for i in range(len(tab_rmsd)) :
        dico_rmsd = tab_rmsd[i]
        
        for res in dico_rmsd.keys() :
            if dico_rmsd[res] > (dico_rmsd["rmsd_moyen"] + seuil) :
                if res not in AA :
                    AA[res] = {}
                    AA[res]["count"] = 0
                    AA[res]["color"] = colors_list[nColor]
                    nColor+= 1
                AA[res]["count"] += 1
                plt.plot(i , dico_rmsd[res], color = AA[res]["color"],
                         marker='D', markersize = 4  )
    
    handler.append(plt.plot([],[], "go",ms = 1 ,color="green",label='rmsd par AA', mec=None)[0])
    for res in AA.keys() :
        handler.append( plt.plot([],[], marker="D", ms=3, ls="", mec=None,
                         color= AA[res]["color"],
                         label=res+" ("+str(AA[res]["count"])+" fois)" )[0])
    plt.xlabel('n°model' , fontsize = 18)
    plt.ylabel('rmsd' , fontsize=18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.title('Graphique des rmsd moyen et de chaque acide aminé pour chaque model,\nen prenant comme référence le model n°0',
              fontsize = 18)
    plt.legend(fontsize ='large' ,  markerscale = 2, handles=handler)
    plt.savefig("rmsd_moyen_et_parRes_"+str(len(tab_rmsd))+"_model.png")



def rmsdCroisee(TabPDB) :
    """
    propos : calculer la matrice contenant les rmsd entre chaque structure
             contenus dans TabPDB
    input : 
        TabPDB = un tableau de dicos de type dPDB;
                 le i-ème élément du tableau correspond à la i-ème structure
                 présent dans le fichier PDB
    output :
        rmsdCroise = matrice carré contenant les rmsd entre chaque structure
                     contenus dans TabPDB
    """
    
    rmsdCroise = np.zeros((len(TabPDB) , len(TabPDB)))
    for num1 in range(len(TabPDB))  :
        dPDB1 = TabPDB[num1]
        # Le calcul peu être long, en effet la complexité est de n² (en prenant
        # n x n les dimensions de la matrice carrée), donc l'utilisateur est prévenu
        # de l'avancée de la fonction
        print("The computation of the matrice containing each rmsd between every conformations is at the "+str(num1)+"th line out of "+str(len(TabPDB)-1) )
        for num2 in range(num1,len(TabPDB)) :
            nb_res = 0 # nombre de résidu dans la protéine
            dPDB2 = TabPDB[num2]
            rmsd_moyen = 0
            for chain in dPDB1["chains"] :
                for res in dPDB1[chain]["reslist"] :
                    nb_res += 1
                    center1 = getCentroid( dPDB1[chain][res] )
                    center2 = getCentroid( dPDB2[chain][res] )
                    rmsd_moyen += rmsd(center1 , center2)
            rmsd_moyen /= nb_res
            rmsdCroise[num1 , num2 ] = rmsd_moyen
            rmsdCroise[num2 , num1 ] = rmsd_moyen

    return rmsdCroise

def plotRmsdCroisee(rmsdCroise):
    """ 
    propos : Représenter la matrice de rmsd
    input : 
        rmsdCroise = matrice carré contenant les rmsd entre chaque structure
                     contenus dans TabPDB
    output : plot une de la matrice carrée avec une échelle de couleur, dans le format .png
    """
    plt.figure(figsize=(16,12))
    plt.pcolor(rmsdCroise, cmap='coolwarm' )
    
    plt.colorbar()
    plt.xlabel('n°model' , fontsize = 18)
    plt.ylabel('n°model' , fontsize = 18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.title('Rmsd moyen entre chaque modèle',
              fontsize = 18)
    plt.savefig("rmsd_croisee_"+str(len(rmsdCroise))+"_model.png")
    return
    








def computesProteinCenter(dPDB) :
    """ propos : Trouver le centre du structure protéique
        input : un fichier PDB
        output : Une liste contenant les coordonnées en 3D du centre
                 de la protéine
        principe : Le barycentre de chaque résidu est calculé,
                   le centre de la protéine est le barycentre de ceux-ci
    """

    centreProt = [0, 0, 0]
    nRes = 0
    for chain in dPDB["chains"] :
        for res in dPDB[chain]["reslist"] :
            centreRes = getCentroid( dPDB[chain][res] )
            centreProt[0] += centreRes[0]
            centreProt[1] += centreRes[1]
            centreProt[2] += centreRes[2]
            nRes += 1
    centreProt[0] /= nRes
    centreProt[1] /= nRes
    centreProt[2] /= nRes

    return centreProt







def enfouissementExtractor(TabPDB) :
    """
    propos : calculer les valeurs d'enfouissemnt euclidien et créer plusieurs tableau
            et dictionnaire qui vont servir pour les graphiques à plotter
    input :
        TabPDB = un tableau de dicos de type dPDB;
                 le i-ème élément du tableau correspond à la i-ème structure
                 présent dans le fichier PDB
    output :
        enfouissementEuclidien =  un tableau dont le i-ème élément est un dictionnaire
                                  associé au i-ème modèle. Chaque dico contient les distances
                                  entre le centre de la protéine et le centre de chaque résidu
        tab_enfouissement_point = tableau contenant toutes les valeurs d'enfouissement 
                                   de tous les résidus à la suite
        enfouissement_moyen_par_res = dictionnaire contenant, pour chaque résidu,
                                      son enfouissement moyen

    """
    enfouissementEuclidien = []
    tab_enfouissement_point = []
    enfouissement_moyen_par_res = {} 
    for dPDB in TabPDB  :
        #Calcul du centre de la protéine
        centreProt = computesProteinCenter(dPDB)
    
        #Calcul distance euclidienne centre-AA
        dico_enfouissement = {}
        for chain in dPDB["chains"] :
            for res in dPDB[chain]["reslist"] :
                #dico_enfouissement est un dico contenant les distances entre le centre de
                # la protéine et le centre de chaque résidu, pour le modèle en train
                # d'être lu
                dico_enfouissement[res] = distancePoints(centreProt , getCentroid(dPDB[chain][res]))
                tab_enfouissement_point.append(dico_enfouissement[res])
                
                if(res not in enfouissement_moyen_par_res):
                    enfouissement_moyen_par_res[res] = dico_enfouissement[res]
                else:
                    enfouissement_moyen_par_res[res] += dico_enfouissement[res]
        
        enfouissementEuclidien.append(dico_enfouissement)
        
    for res in enfouissement_moyen_par_res.keys() :
        enfouissement_moyen_par_res[res] /= len(enfouissementEuclidien )
        
    return enfouissementEuclidien ,tab_enfouissement_point , enfouissement_moyen_par_res
    



def giration(enfouissementEuclidien) :
    """
    propos : calculer les valeurs de giration pour chaque modèle
    input :
        enfouissementEuclidien =  un tableau dont le i-ème élément est un dictionnaire
                                  associé au i-ème modèle. Chaque dico contient les distances
                                  entre le centre de la protéine et le centre de chaque résidu
    ouput :
        giration = tableau dont le i-ème élément est la valeur maximal d'enfouissement du
                   i_ème model            
    """
    #pour chaque modèle, il faut rechercher l'enfouissement maximal
    giration = []
    for dico_enfouissement in enfouissementEuclidien :
        val_max = 0 ;
        for valeur in dico_enfouissement.values() :
            if(val_max < valeur):
                val_max = valeur
        giration.append(val_max)
    return giration




def plotGiration(giration) :
    """
    propos : représenter l'évolution du rayon de giration 
    input :
        giration = tableau dont le i-ème élément est la valeur maximal d'enfouissement du
                   i_ème model  
    output : fichier .png qui est un plot du rayon de giration selon le model
    """
    nModel = list(range(0, len(giration)))
    plt.figure(figsize=(16,16))
    plt.plot( nModel , giration ,ms= 3.5 , color= "black" , marker='o' , linestyle ='-',label='rmsd moyen' )
    plt.xlabel('n°model' , fontsize = 18)
    plt.ylabel('giration' , fontsize=18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.title('Graphique du rayon de giration pour chaque model',
              fontsize = 18)
    plt.savefig("giration_"+str(len(giration))+"_model.png")



def figureEnfouissementRmsdToutLesPoints(tab_rmsd_point, tab_enfouissement_point , nb_model):
    """
    propos: plot
    input:
        nb_model = nombre de résidu dans la protéine
        tab_enfouissement_point = tableau contenant toutes les valeurs d'enfouissement
                                   de tous les résidus à la suite
        tab_rmsd_point = tableau contenant tous les rmsd de tous les résidus à la suite
    output: plot de toutes les valeurs de rmsd et d'enfouissement de tout les modèles en format .png

    """
    plt.figure(figsize=(16, 16))
    plt.plot(tab_rmsd_point, tab_enfouissement_point,
             c='k', marker='o', markersize=2, alpha=0.4, linestyle='None')
    plt.title(
        'Graphique de la distance euclidienne avec le centre de la protéine\nen fonction du rmsd pour chaque résidu et pour chaque model',
        fontsize=18)
    plt.xlabel('rmsd', fontsize=18)
    plt.ylabel('distance du centre de la protéine', fontsize=18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.savefig("enfouissement.euclid_rmsd_parRes." + str(nb_model) + ".model.png")



def figureEnfouissementRmsdMoyParRes(rmsd_moyen_par_res , enfouissement_moyen_par_res , res_pour_analyse_locale , nb_model):
	"""
    propos: plot
    input:
        nb_model = nombre de modèle
		res_pour_analyse_locale = liste des résidus à mettre en évidence
        rmsd_moyen_par_res = dictionnaire contenant pour chaque résidu son rmsd moyen
        enfouissement_moyen_par_res = dictionnaire contenant, pour chaque résidu,
                                      son enfouissement moyen
    output: plot des rmsd moyens et des enfouissements moyens par résidu, dans le format .png
	"""
	abscisse_rmsd_enfouissement_moy = []
	ordonne_rmsd_moy = []
	ordonne_enfouissement_moy = []
	for res, val_rmsd in rmsd_moyen_par_res.items():
		abscisse_rmsd_enfouissement_moy.append(res)
		ordonne_rmsd_moy.append(val_rmsd)
	for enfouissement in enfouissement_moyen_par_res.values():
		ordonne_enfouissement_moy.append(enfouissement)
	plt.figure(figsize=(26, 26))
	plt.plot(ordonne_rmsd_moy, ordonne_enfouissement_moy, marker="o", linestyle="none")
	for res in res_pour_analyse_locale:
		plt.plot(rmsd_moyen_par_res[res], enfouissement_moyen_par_res[res], marker="o", c="r", linestyle="none" , label = "résidus\nmarqué par\nl'utilisateur")
	for res in rmsd_moyen_par_res.keys():
		plt.text(rmsd_moyen_par_res[res], enfouissement_moyen_par_res[res], res)
	plt.xlabel('rmsd moyen', fontsize=18)
	plt.ylabel('distance euclidienne moyenne avec le centre de masse de la protéine', fontsize=18)
	plt.xticks(fontsize=16)
	plt.yticks(fontsize=16)
	plt.title('Rmsd moyen enfouissement (euclidien) moyen par résidu',
		      fontsize=18)
	plt.savefig("rmsd_enfouis_moy_par_res_"+str(nb_model)+"_model.png")



def makeFigurePourUnResidu( res , tab , tabRef  , choix):
    """
    propos : Représenter, pour un résidu donnée , l'évolution du rmsd (par rapport 
             à la structure de référence) ou de son enfouissement
    input : res = numéro du résidu à représenter
            tab = un tableau de dictionnaire, dont le i-ème élément est le
                  dictionnaire correspondant au i-ème modèle.
                  Chaque dictionnaire contient les rmsd ou l'enfouissemnt par résidu
            choix = 'rmsd' ou 'distance euclidienne du centre'
            tabRef = rmsd moyen ou rayon de giration 
    output : plot du rmsd ou de l'enfouissemnt d'un résidu donnée selon le modèle
             dans le format .png
    """
    if choix == 'distance euclidienne du centre':
        greenLabel = "rayon de giration"
    else : 
        greenLabel = "rmsd moyen"
    nModel = list(range(0 , len(tab)))

    plt.figure(figsize=(18,16))
    ordonee = []
    for dico in tab:
	    ordonee.append(dico[res])
    plt.plot(nModel,ordonee ,c = "r", linestyle = "-" ,  linewidth = 0.3 ,
	         marker = 'x', markersize = 4 , label = res)
    plt.plot(nModel, tabRef , c = "g" , linewidth = 0.3 ,
	         marker = 'o', markersize = 2  , label = greenLabel)
    plt.xlabel('model' ,  fontsize =18)
    plt.ylabel(choix , fontsize =18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend()
    plt.savefig(choix+"_locale_residu_"+res+"_"+str(len(tab))+"_model.png")









def Rij ( i , j ):
    """
    input : i et j deux listes contenant les coordonnées de l'atome i, 
    respectivement j
    output : vecteur Rij 
    """
    x , y  , z = i[0]-j[0] , i[1]-j[1] , i[2]-j[2]
    return (x , y , z )
    


def Norme ( vec ): 
    """
    input :  vec, un vecteur
    output : norme de vec
    """
    return math.sqrt( vec[0]**2 + vec[1]**2 + vec[2]**2 ) 


def getXYZ (dPDB):
    """
    input : dPDB dico de type PDB
    output :  X,Y,Z contenant les coordonnées de tous les atomes de la protéine
    """
    X , Y , Z = [] , [] , []
    for chain in dPDB["chains"] :
        for res in dPDB[chain]["reslist"] : 
            for atom in dPDB[chain][res]["atomlist"] :
                j = dPDB[chain][res][atom]
                X.append(j["x"])
                Y.append(j["y"])
                Z.append(j["z"])  
    return X , Y , Z
    
    
    
    
def Env_i ( rc , i , X,Y,Z ):
    """
    input : rc le seuil et i, une liste des 
    coordonées de l'atome i se trouvant dans dPDB;
    X,Y,Z contenant les coordonnées de tous les atomes de la protéine
    output: 3 listes des coordonnées des atomes situés à une distance de i <= rc  
    """
    Xenv , Yenv , Zenv = [] , [] , []
    for k in range(len(X)):
        if i != [ X[k] , Y[k] , Z[k] ] :
            L1 = [i[0] , i[1] , i[2] ]
            L2 = [X[k] , Y[k] , Z[k] ]
            if distancePoints(L1,L2) <= rc :
                Xenv.append(X[k])
                Yenv.append(Y[k])
                Zenv.append(Z[k])  
    return Xenv , Yenv , Zenv

   

def sumVect(Vec) : 
    """
    input: liste Vec contenant les coordonnées des vecteurs à sommer
    output: S liste contenant les trois coordonnées du vecteur résultant de la
    somme
    """
    S = [0,0,0]
    for i in range(len(Vec)):
        for xyz in range(3) : 
            S[xyz]+= Vec[i][xyz]
    return S

def CVi ( rc , i , X,Y,Z) : 
    """
    input :rc le seuil 
    i, une liste comportant les coordonnées de l'atome i; 
    X,Y,Z contenant   les coordonnées de tous les atomes de la protéine
    output : variance circulaire de l'atome i
    """
    Xenv , Yenv , Zenv = Env_i ( rc , i , X,Y,Z )
    ri = [ Rij(i ,   [ Xenv[j] , Yenv[j] , Zenv[j]] ) for j in range(len(Xenv))]
    Normeri = [Norme(rij) for rij in ri ]
    
    riInf = [ri[j] for j in range(len(ri)) if Normeri[j] < rc ]# atomes dont la norme est infèrieur à rc
    NormeriInf = [normeij for normeij in Normeri if normeij < rc ]
    ni = len(riInf)
    
    riNormed = []    
    for j in range(len(riInf)) : 
        add = []
        for k in range(len(riInf[j])) :
            add.append(riInf[j][k] / NormeriInf[j]) 
        riNormed.append(add)
        
    return (1 - 1/ni * Norme(sumVect(riNormed)))
    


def FigureCV_res(TabPDB , res , rc):
    """
    propos: Représenter l'évolution de la variance ciculaire au cours du temps 
    input: TabPDB = un tableau de dicos de type dPDB;
                     le i-ème élément du tableau correspond à la i-ème structure
                     présent dans le fichier PDB 
            res = numéro du résidu dont la variance circulaire est suivi
            rc = seuil utilisé pour calculé la variance circulaire
    output: plot de la variance circulaire en fonction des modèles, en 
            format .png
    """

    tabCV_res = [] # tableau dont le i-ème élément est la variance circulaire d'un résidu
    for dPDB in TabPDB  :
        X , Y ,Z = getXYZ(dPDB)
        #la chaine où se trouve le résidu est recherché
        n = 0
        chain = dPDB["chains"][n]
        while res not in dPDB[chain]["reslist"] :
            n = n+1
            chain = dPDB["chains"][n]
        #la variance circulaire est calculé à partir du barycentre du résidu
        i = getCentroid( dPDB[chain][res] )
        tabCV_res.append(CVi(rc , i , X ,Y , Z))
    nModel = list(range(0, len(TabPDB)))
    plt.figure(figsize=(20, 16))
    plt.plot(nModel , tabCV_res, linestyle = "-" ,  linewidth = 0.3 ,
	         marker = 'x', markersize = 4 ,)
    plt.xlabel('model', fontsize=18)
    plt.ylabel('variance circulaire', fontsize=18)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    plt.savefig("VarianceCirculaire_res_"+res+"_"+str(len(TabPDB))+".model.png")












def writeRmsdMoy_Giration(tab_rmsd_moy ,giration ):
    """
    propos : écriture de fichier
    input :
        tab_rmsd_moy = tableau dont le i-ème élément est le rmsd moyen du i-ème modèle
        giration = tableau dont le i-ème élément est le rayon de giration du
                   i_ème model         
    output : fichier texte qui spécifie pour chaque modèle le rmsd moyen et le rayon de giration
    """
    outFile_rmsd_giration = open("rmsd_giration."+str(len(tab_rmsd_moy))+".model.txt", "w")
    outFile_rmsd_giration.write("model\trmsd\trayon de giration")
    nModel = list(range(0 , len(tab_rmsd_moy)))

    for n in nModel:
        outFile_rmsd_giration.write("\n"+str(n)+"\t"+str(tab_rmsd_moy[n])+"\t"+str(giration[n]))
    outFile_rmsd_giration.close()
    
    


def writeRmsdEnfouissementMoyen_ParRes(rmsd_moyen_par_res ,enfouissement_moyen_par_res , nbrDeModel) :
    """
    propos : écriture de fichier
    input :
        rmsd_moyen_par_res = dictionnaire contenant pour chaque résidu son rmsd moyen
        enfouissement_moyen_par_res = dictionnaire contenant, pour chaque résidu,
                                      son enfouissement moyen
             
    output : fichier texte qui spécifie, pour chaque résidu, le rmsd moyen 
             et l'enfouissement (euclidien) moyen
    """
    outFile_rmsd_enfouissement_moyen_par_res = open("rmsd_enfouissement_moyen_par_res."+str(nbrDeModel)+".model.txt", "w")
    outFile_rmsd_enfouissement_moyen_par_res.write("res\trmsd\tdistance euclidienne du centre")
    for res in rmsd_moyen_par_res.keys() :
        outFile_rmsd_enfouissement_moyen_par_res.write("\n"+res+"\t"+
                                                       str(rmsd_moyen_par_res[res])+
                                                       "\t"+
                                                       str(enfouissement_moyen_par_res[res]))
    outFile_rmsd_enfouissement_moyen_par_res.close()



def writeRmsdGiration_ParRes_ParModel( tab_rmsd , enfouissementEuclidien ) :
    """
    propos : écriture de fichier
    input :
        tab_rmsd = un tableau de dictionnaire, dont le i-ème élément est le
                   dictionnaire correspondant au i-ème modèle. Chaque dictionnaire
                   contient le rmsd moyen du modèle et les rmsd par résidu
        enfouissementEuclidien = un tableau dont le i-ème élément est un dictionnaire
                                  associé au i-ème modèle. Chaque dico contient les distances
                                  entre le centre de la protéine et le centre de chaque résidu
          
    output : fichier texte qui spécifie pour chaque modèle et pour chaque résidu
             le rmsd et l'enfouissement
    """
    outFile_rmsd_enfouissement_par_res = open("rmsd_enfouissement_par_res."+str(len(tab_rmsd))+".model.txt", "w")
    outFile_rmsd_enfouissement_par_res.write("res")
    for nRes in enfouissementEuclidien[0] :
        outFile_rmsd_enfouissement_par_res.write("\t"+str(nRes))
    outFile_rmsd_enfouissement_par_res.write("\nmodel")
    for model in range(len(tab_rmsd)) :
        outFile_rmsd_enfouissement_par_res.write("\n"+str(model))
        dico_enfouissement = enfouissementEuclidien[model]
        dico_rmsd = tab_rmsd[model]
        
        outFile_rmsd_enfouissement_par_res.write("\n\trmsd")
        for res  in dico_rmsd.keys() - {'rmsd_moyen'} :
            outFile_rmsd_enfouissement_par_res.write("\t"+str(dico_rmsd[res]))
            
        outFile_rmsd_enfouissement_par_res.write("\n\tdistance du centre")
        for res  in dico_enfouissement.keys() :
            outFile_rmsd_enfouissement_par_res.write("\t"+str(dico_enfouissement[res]))
    outFile_rmsd_enfouissement_par_res.close()
    
    
def usage():
    """
    propos : expliquer à l'utilisateur
    output : une explication de l'usage du script affichée dans la sortie standard
             et arrêt du script
    """
    print("\nUsage:"+
    	  "\n-file <FILE_NAME> -res <numéros des résidus à analyser> -seuil 3 -rc 8 -c -w"+
          "\n\n-file:\tnom du fichier, c'est le seul paramètre obligatoire, son absence provoque l'arrêt du script."+
          "\n-res:\tles noms des résidus doivent être séparés par une virgule, sans espace"+
          "\n-seuil:\tpour définir le seuil, les résisdus dont le rmsd > rmsd_moyen + seuil   sont mis en évidence. Par défaut seuil = 3"+
          "\n-rc:\tseuil en Å à partir duquel la variance circulaire est calculée. Par défaut rc = 8 Å"+
          "\n-c:\tà écrire si l'utilisateur veut obtenir un plot représentant le rmsd "+
                 "entre chaque modèle."+
          "\n\tATTENTION : l'algorithme pour calculer la matrice à l'origine de cet figure est "+
                             "d'une compléxité n², avec n le nombre de modèle,"+
                           "\n\t\tsi il y a 2000 modèles, le calcul peut durer plus d'une heure."+
		 "\n-w:\tà entrer si l'utilisateur veut écrire les fichiers textes"+
          "\n\nEntrez : <NOM_DU_SCRIPT_PYTHON> -file md_prot_only_skip100.pdb -res 76,80,39,35 -c -w "+
             "&& <NOM_DU_SCRIPT_PYTHON> -file md_prot_only_skip10.pdb -res 76,80,39,35 -w -c  "+
             "pour créer les mêmes fichiers que dans le rapport"+
          "\n\n\"-help\" pour revoir cette aide")
    sys.exit("\nArrêt du script\n")
