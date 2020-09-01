#!/usr/bin/env python3

import sys
import numpy as np
import  MDAnalysis as mda

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