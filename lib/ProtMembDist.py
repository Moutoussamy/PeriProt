#!/usr/bin/env python3

import common
import MDAnalysis as mda
from scipy.spatial import distance
import matplotlib.pyplot as plt

"""
SET OF FUNCTION TO CALCULATE THE PROT-MEMB DISTANCE
"""


def WhichDist():
    """
    Which type of distance do you want
    :return: 1, 2 or 3 (distance type)
    """
    whichdist = input("""
    
    Which dist:
    1: Prot (COM) - Memb (COM) distance
    2: Prot (COM) - Phosphate Plane (COM) distance
    3: Custom distance between two atoms (one in the prot. and one on the memb.)
    
    Selection:
    """)

    if whichdist not in [1,2,3]: #recrusion, to get the right value
        WhichDist()

    return whichdist

def CustomSelection():
    """
    Custom selection for distance calculation (option 3)
    :return: index of atom for the protein and index of atom for the bilayer
    """
    selprot = input('Atom ID on the protein: ') #index for the prot.
    selmemb = input('Atom ID on the membrane: ') #index for the memb.

    return selprot,selmemb




def WhichSide(psf,segmemb,segprot,univers):
    """
    Find where is the protein: under or above the bilayer?

    :param psf: PSF file
    :param segmemb: segid of the membrane
    :param segprot: segid of the protein
    :param univers: mdAnalysis univers
    :return: a list containing the phosphate resid of the leaflet close to the protein
    """
    lipids = common.get_lipid_resid_list(psf,segmemb) # get the list of the lipids resid

    prot = univers.select_atoms("segid %s"%segprot)
    memb = univers.select_atoms("segid %s"%segmemb)

    #center of mass:
    prot = prot.center_of_mass()
    memb = memb.center_of_mass()

    good_lipids = []

    for lipid in lipids:
        lipid_coor = univers.select_atoms("segid %s and (resid %s)"%(segmemb,lipid))
        lipid_coor = lipid_coor.positions[2]

        if lipid_coor[2] > memb[2]: #protein above bilayer
            if lipid_coor[2] < prot[2]:
                good_lipids.append(lipid)

        elif lipid_coor[2] < memb[2]: #protein under bilayer
            if lipid_coor[2] > prot[2]:
                good_lipids.append(lipid)

    return good_lipids

def plot_results_dist(frames,dist,outname):
    """
    Plot the results
    :param frames: list contaning the frames number
    :param dist: list containing the dist
    :return: none
    """
    fig, ax = plt.subplots(figsize=(10,8 ))

    plt.rcParams["xtick.labelsize"] = 20
    plt.rcParams["ytick.labelsize"] = 20
    plt.rcParams["axes.labelsize"] = 24
    plt.rcParams["axes.labelsize"] = 24

    plt.plot(frames,dist,color = "black",linewidth = 1)
    plt.xlabel("Frame")
    plt.ylabel(r"Prot. - Memb. Distance ($\AA$)")
    plt.savefig("%s_prot_memb_distance.png"%outname, dpi=300, papertype="a4", orientation="landscape",\
                format="png", bbox_inches='tight')


def ComputeDistance(psf,dcd,segprot,segmemb,outname):
    """
    Calculate the distance between the protein and the bilayer.
    three type of disatnce are available:
    1: Prot (COM) - Memb (COM) distance
    2: Prot (COM) - Phosphate Plane (COM) distance
    3: Custom distance between two atoms (one in the prot. and one on the memb.)


    :param psf: PSF file
    :param dcd: trajectory filee
    :param segprot: segid of the protein
    :param segmemb: segid of the membrane
    :param outname: name of the output file
    :return: None
    """

    whichdist = WhichDist()
    univers = mda.Universe(psf, dcd)
    flag = 0
    output = open("%s_prot_memb_distance.dat"%outname,"w")
    output.write("#Frame,dist\n")


    frames = []
    dists = []

    for ts in univers.trajectory:

        frames.append(ts.frame)

        if whichdist == 1:
            ###Prot (COM) - Memb (COM) distance
            prot = univers.select_atoms("segid %s"%segprot)
            memb = univers.select_atoms("segid %s"%segmemb)

            prot = prot.center_of_mass()
            memb = memb.center_of_mass()
            dist = distance.euclidean(prot, memb)
            dists.append(dist)
            output.write("%i,%f\n"%(ts.frame,dist))

        elif whichdist == 2:
            ###Prot (COM) - Phosphate Plane (COM) distance
            if flag == 0:
                good_lipids = WhichSide(psf, segmemb,segprot,univers)
                flag = 1

            prot = univers.select_atoms("segid %s"%segprot)
            memb = univers.select_atoms("segid %s and (resid %s)"%(segmemb," ".join(good_lipids)))

            prot = prot.center_of_mass()
            memb = memb.center_of_mass()
            dist = distance.euclidean(prot, memb)
            dists.append(dist)
            output.write("%i,%f\n" % (ts.frame, dist))



        elif whichdist == 3:
            ###Custom distance between two atoms (one in the prot. and one on the memb.)

            if flag == 0:
                sel_prot,sel_memb = CustomSelection()
                flag = 1

            prot = univers.select_atoms("index %i"%sel_prot)
            memb = univers.select_atoms("index %i"%sel_memb)

            prot = prot.positions
            memb = memb.positions
            dist = distance.euclidean(prot, memb)
            dists.append(dist)
            output.write("%i,%f\n" % (ts.frame, dist))

    plot_results_dist(frames, dists, outname)
