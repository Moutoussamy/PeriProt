#!/usr/bin/env python3

import numpy as np
import  MDAnalysis as mda
import pandas as pd
import matplotlib.pyplot as plt

"""
SET OF FUNCTION TO EVALUATE THE DEPTH OF ANCHORING
the depth of anchoring correponds to the z-distance
between the CA atom anf the average phosphate plane
"""

def initializeDico(dico_resid,residues_list):
    """
    Initialize a dico to contain the distances
    :param dico_resid: dico to initialize
    :param residues_list: list of residues
    :return: a dico with residues ids as keys and empty list as value
    """

    for resid in residues_list:
        dico_resid[resid] = []

    return dico_resid

def center_of_mass_z(coord):
    """
    Calculate the center of mass in the z-direction
    :param coord: set od coordinates
    :return: the center of mass in the z-direction
    """
    mean_z = np.mean(coord[:,2])

    return mean_z

def Getphosphateids(psf,segidMEMB):
    """
    Get lipids ids of lipids cointaining phosphate atoms
    :param psf: PSF file
    :param segidMEMB: segid of the mebrane
    :return: a list containg the ids of lipids containg phosphate
    """
    phos_ids = []
    with open(psf) as inputfile:
        for line in inputfile:
            line = line.split()
            if len(line) == 9:
                if line[1] == segidMEMB:
                    if line[4] == "P":
                        phos_ids.append(line[2])
    return phos_ids


def WhereIsTheProtein(univers,segidMEMB,psf):
    """
    find where is the protein below or under the membrane
    :param univers: univers from MDAnalysis
    :param segidMEMB: membrane segid
    :param psf: PSF file
    :return: resid of lipids that containg phosphate and that are close to the PMP
    """

    good_phos_ids = []
    prot = univers.select_atoms("protein")
    memb = univers.select_atoms("segid %s" % (segidMEMB))
    prot_com = center_of_mass_z(prot.positions)
    memb_com = center_of_mass_z(memb.positions)

    resids_phos = Getphosphateids(psf,segidMEMB)

    if prot_com > memb_com:
        for phos in resids_phos:
            phosphate_coord = univers.select_atoms("segid %s and (name P) and (resid %s)" % (segidMEMB,phos))
            z_phos = phosphate_coord.positions[0][2]
            if z_phos > memb_com:
                good_phos_ids.append(phos)

    if prot_com > memb_com:
        for phos in resids_phos:
            phosphate_coord = univers.select_atoms("segid %s and (name P) and (resid %s)"\
                                                   % (segidMEMB,phos))

            z_phos = phosphate_coord.positions[0][2]
            if z_phos > memb_com:
                good_phos_ids.append(phos)

    return good_phos_ids

def write_results_depth(avg,sd,residues_list,outname):
    """
    Wirite the resulst in a file in a csv format
    :param dico_resid: dico conating the DOA
    :param residues_list: residues list
    :return: two list conating the avg. DOA and the associated sd
    """

    output = open("%s_depth_of_anchoring.dat"%outname,"w")
    output.write("#Resid,DOA_avg,DOA_sd\n")

    for i in range(len(residues_list)):

        output.write("%s,%f,%s\n"%(residues_list[i],avg[i],sd[i]))

    output.close()

def plot_results_depth(avg,sd, residues_list,outname):
    """
    Plot the results
    :param avg: average DOA
    :param sd: standard deviation od the DOA
    :param residues_list: list of residue
    :return: none
    """
    fig, ax = plt.subplots(figsize=(10,8 ))

    plt.rcParams["xtick.labelsize"] = 20
    plt.rcParams["ytick.labelsize"] = 20
    plt.rcParams["axes.labelsize"] = 24
    plt.rcParams["axes.labelsize"] = 24

    avg = np.array(avg)
    sd = np.array(sd)
    residues_list = [float(i) for i in residues_list]
    plt.plot(residues_list,avg,color = "black",linewidth = 1)
    plt.fill_between(residues_list, avg - sd, avg + sd, alpha=0.3, color="green")
    plt.xticks(np.arange(min(residues_list), max(residues_list) + 1, 10.0), rotation=45)
    plt.xlabel("Residue")
    plt.ylabel(r"Depth of Anchoring ($\AA$)")
    plt.hlines(0,min(residues_list),max(residues_list))
    plt.savefig("%s_depth_of_anchoring.png"%outname, dpi=300, papertype="a4", orientation="landscape",\
                format="png", bbox_inches='tight')

def GetAVGandSD(depht_data_frame,residues_list):
    """
    Get the average depth and its standard deviation
    :param depht_data_frame: pd dataframe containg all the depth
    :param residues_list: list of residues
    :return:
    """
    depth_avg = []
    depth_sd = []

    for resid in residues_list:
        depth_avg.append(np.mean(depht_data_frame[resid]))
        depth_sd.append(np.std(depht_data_frame[resid]))

    return depth_avg,depth_sd

def RunDepthOfAnchoring(psf,dcd,residues_list,segidMEMB,first_frame,last_frame,skip,outname):
    """
    RUN THE DEPTH OF ANCHORING (DOA) CALCULATION
    :param psf: PSF file
    :param dcd: trajectory file
    :param residues_list: residues list
    :param segidMEMB: membrane segid
    :return: none
    """
    univers = mda.Universe(psf, dcd)

    dico_resid = {}
    dico_resid = initializeDico(dico_resid,residues_list)
    flag = 0

    frame_to_read = first_frame

    frames = []

    for ts in univers.trajectory:
        if ts.frame == frame_to_read and frame_to_read < last_frame :
            frames.append(ts.frame)
            if flag == 0:
                good_phos_ids = WhereIsTheProtein(univers,segidMEMB,psf)
                phos_plane = univers.select_atoms("segid %s and (resid %s)"%(segidMEMB," ".join(good_phos_ids)))
                flag = 1


            prot = univers.select_atoms("(resid %s) and (name CA)" %(" ".join(residues_list)))
            z_prot = prot.positions[:,2] - np.mean(phos_plane.positions[:,2])

            if 'depht_data_frame' not in locals():
                depht_data_frame = z_prot

            else:
                depht_data_frame = np.vstack((depht_data_frame, z_prot))

            frame_to_read += skip

    depht_data_frame = pd.DataFrame(depht_data_frame, columns=residues_list, index=frames)

    depth_avg, depth_sd = GetAVGandSD(depht_data_frame,residues_list)

    write_results_depth(depth_avg,depth_sd,residues_list,outname)
    plot_results_depth(depth_avg,depth_sd,residues_list,outname)


