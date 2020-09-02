#!/usr/bin/env python3

import sys
import numpy as np
import  MDAnalysis as mda
from MDAnalysis.analysis import contacts

"""
SET OF FUNCTION TO EVALUATE HYDROPHOBIC CONTACT
"""


def filterCanditatesAtom(atom_name,charge,candidate_list):
    """
    select candidate atoms for hydroiphobic interaction according to their charge
    charge of carbon should be below 0.3 and absolute charge of hydrogen should be below 0.1
    :param atom_name: name of the atom
    :param charge: charge of the atom
    :param candidate_list: list to append
    :return: candidtae list append with the appropriate atoms
    """
    if atom_name[0] == 'C':
        if abs(charge) < 0.3:
            if atom_name not in candidate_list:
                candidate_list.append(atom_name)
    elif atom_name[0] == 'H':
        if abs(charge) < 0.1:
            if atom_name not in candidate_list:
                candidate_list.append(atom_name)

    return candidate_list

def HydroCandidateSelection(segidMEMB,segidPROT,psf):
    """
    Get the selected atoms for hydrophobic contact
    :param segidMEMB: segid of MEMB
    :param segidPROT:  segid of PROT
    :param psf: PSF file
    :return: list of candiates in the proteion and the membrane
    """

    prot_hydro_candidates = []
    memb_hydro_candidates = []

    with open(psf) as inputfile:
        for line in inputfile:
            line = line.split()
            if len(line) == 9:
                if line[1] == segidMEMB:
                    memb_hydro_candidates = filterCanditatesAtom(line[4],float(line[6]),memb_hydro_candidates)
                elif line[1] == segidPROT:
                    prot_hydro_candidates = filterCanditatesAtom(line[4], float(line[6]), prot_hydro_candidates)

    return prot_hydro_candidates, memb_hydro_candidates


def CountLowerThanCutoff(r,Cutoff):
    """
    Count the number of time that we have a distance lower that the cut-off
    :param r: ditance
    :param Cutoff: cutoff
    :return:number of time that we have a distance lower that the cut-off
    """
    count = 0

    for distance in r:
        if distance < Cutoff:
            count += 1

    return count

def is_any_closer(r, r0, dist=3.0):
    """
    find atom close enought to make hydrophobic contact
    :param r: distance between atoms
    :param dist: cut-off
    :return:
    """
    return CountLowerThanCutoff(r,dist)


def write_hydro_results(avgcontact, residues_list,psf):
    """
    write the results of the calculation in a file:

    :param avgcontact: average contact per frame
    :param residues_list: list of residues
    :param psf: PSF file
    :return:
    """
    output = open("average_hydrophobic_contact_per_frame.dat","w")
    output.write("#resid,avg_contact_per_frame\n")
    for resid in residues_list:
        output.write("%s,%f\n"%(resid,avgcontact[resid]))

    output.close()

def RunHydroAnalysis(segidMEMB,psf,dcd,residues_list,prot_candidates,memb_candidates):
    """
    Run calculation of hydrophobic contact using MDanalysis
    :param segidMEMB: segid of the membrane
    :param psf: PSF file
    :param dcd: trajectory file
    :param residues_list: liste of residues of interest
    :param prot_candidates: atom candidates in the prot
    :param memb_candidates: atom candidates in the memb
    :return:
    """
    univers = mda.Universe(psf, dcd)
    sel_memb = "(name %s) and (segid %s)" % (" ".join(memb_candidates),segidMEMB)
    memb = univers.select_atoms(sel_memb)

    avgcontact = {}

    for resid in residues_list:
        sel_prot = "(name %s) and (protein) and (resid %s)" % (" ".join(prot_candidates), str(resid))
        prot = univers.select_atoms(sel_prot)
        contacthydro = contacts.Contacts(univers, selection=(sel_memb, sel_prot), method=is_any_closer, refgroup=(memb, prot))
        contacthydro.run()
        avgcontact[resid] = np.mean(contacthydro.timeseries[:,1])

    write_hydro_results(avgcontact, residues_list,psf)
