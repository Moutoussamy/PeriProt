#!/usr/bin/env python3

import common
import  MDAnalysis as mda
from MDAnalysis.analysis import contacts

"""
SET OF FUNCTION TO EVALUATE HYDROPHOBIC CONTACT
"""



def HydroCandidateSelection(psf_info):
    """
    Get the selected atoms for hydrophobic contact
    :param psf_info: class conating infos about the system (topology from ./lib/common.py)
    :return: list of candiates in the proteion and the membrane
    """

    prot_hydro_candidates = []
    memb_hydro_candidates = []


    # FIND CANDIATE ATOM IN THE PROTEIN
    for res_infos in psf_info.protCharge.keys():

        res_infos_decompose = res_infos.split("_")

        if res_infos_decompose[3][0] == "C":
            if abs(psf_info.protCharge[res_infos]) < 0.3:
                if res_infos_decompose[3] not in prot_hydro_candidates:
                    prot_hydro_candidates.append(res_infos_decompose[3])

        elif res_infos_decompose[3][0] == "H":
            if abs(psf_info.protCharge[res_infos]) < 0.1:
                if res_infos_decompose[3] not in prot_hydro_candidates:
                    prot_hydro_candidates.append(res_infos_decompose[3])

    # FIND CANDIATE ATOM IN THE MEMBRANE
    for res_infos in psf_info.MembCharge.keys():
        res_infos_decompose = res_infos.split("_")

        if res_infos_decompose[3][0] == "C":
            if abs(psf_info.MembCharge[res_infos]) < 0.3:
                if res_infos_decompose[3] not in memb_hydro_candidates:
                    memb_hydro_candidates.append(res_infos_decompose[3])

        elif res_infos_decompose[3][0] == "H":
            if abs(psf_info.MembCharge[res_infos]) < 0.1:
                if res_infos_decompose[3] not in memb_hydro_candidates:
                    memb_hydro_candidates.append(res_infos_decompose[3])

    return prot_hydro_candidates, memb_hydro_candidates


def GetNearAtom(close_res,NbAtomPerRes):
    """
    Get atom close to the bilayer based on the residues close to the bilayer
    :param close_res: list of residues close to the bilayer
    :param NbAtomPerRes: dico key = residues; value = [fist atom, last atom]
    :return: list of atoms close to the bilayer
    """

    near_atom = []
    for residues in close_res:
        for atomid in range(NbAtomPerRes[residues][0],NbAtomPerRes[residues][1]):
            near_atom.append(atomid)

    return near_atom

def selectionBasedOnCharge(charge_dico,near_atom_list):
    """
    Select atom based on the charge:
    if charge(carbon) < 0.3 or charge(hydrogen) < 0.1

    :param charge_dico: dico key = atom id ; value = charge
    :param near_atom_list: list of atoms close to the bilayer
    :return:
    """

    selected_atom = []

    dico_nb_atom_selected = {}

    for atom_key in charge_dico.keys():
        atomname = atom_key.split("_")[3]
        atomid = atom_key.split("_")[2]

        if atomname[0] == "C" and  abs(charge_dico[atom_key]) < 0.3:
            if int(atomid) in near_atom_list:
                selected_atom.append(int(atomid))

                if atom_key.split("_")[1] not in dico_nb_atom_selected.keys():
                    dico_nb_atom_selected[atom_key.split("_")[1]] = 1
                else:
                    dico_nb_atom_selected[atom_key.split("_")[1]] += 1

        elif atomname[0] == "H" and abs(charge_dico[atom_key]) < 0.1:
            if atomid in near_atom_list:
                selected_atom.append(int(atomid))

                if atom_key.split("_")[1] not in dico_nb_atom_selected.keys():
                    dico_nb_atom_selected[atom_key.split("_")[1]] = 1
                else:
                    dico_nb_atom_selected[atom_key.split("_")[1]] += 1


    return selected_atom


def GetCandidate(psf_info,close_lipid,close_prot):
    """
    Obtain candiate atom for hydrophobic contact analysis
    :param psf_info: class conating infos about the system (topology from ./lib/common.py)
    :param close_lipid: lipid close to the bilayer
    :param close_prot: amino acids close to the bilayer
    :return: list of candidates for the protein and lkist of candidates for the membrane
    """

    prot_atom_near = GetNearAtom(close_prot,psf_info.NbAtomPerProtRes)
    prot = selectionBasedOnCharge(psf_info.protCharge,prot_atom_near)
    prot = sorted(prot)

    memb_atom_near = GetNearAtom(close_lipid,psf_info.NbAtomPerLipRes)
    memb = selectionBasedOnCharge(psf_info.MembCharge,memb_atom_near)
    memb = sorted(memb)

    return prot,memb



def RunHydroAnalysis(psf, dcd, psf_info, close_lipid,close_prot,outname,segprot,pdb):
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

    prot_candidates, memb_candidates = GetCandidate(psf_info, close_lipid, close_prot)
    prot_candidates = common.convertIntToStr(prot_candidates)
    memb_candidates = common.convertIntToStr(memb_candidates)

    univers = mda.Universe(psf, dcd)

    counting = {}
    nbFrame = 0

    output_avg = open("%s_average_hydrophobic_contact_per_frame.csv"%outname,"w")
    output_avg.write("#Residue,AVG_nb_of_contact\n")

    output_raw = open("%s_hydrophobic_contact_raw_data.csv"%outname,"w")
    output_raw .write("#Frame,ProtResidue,MembRes\n")

    for ts in  univers.trajectory:
        nbFrame += 1
        sel_memb = "index %s" % (" ".join(memb_candidates))
        memb = univers.select_atoms(sel_memb)
        sel_prot = "index %s" % (" ".join(prot_candidates))
        prot = univers.select_atoms(sel_prot)
        dist_mat = mda.analysis.distances.distance_array(prot.positions, memb.positions)

        conctacts =[]

        for i in range(dist_mat.shape[0]):
            for j in range(dist_mat.shape[1]):
                if dist_mat[i,j] < 3:
                    conctacts.append([prot_candidates[i],memb_candidates[j]])


        for i in range(len(conctacts)):
            residues_prot = common.obtainResInfo(int(conctacts[i][0]), psf_info.NbAtomPerProtRes)
            residues_lipids =  common.obtainResInfo(int(conctacts[i][1]), psf_info.NbAtomPerLipRes)

            if residues_prot not in counting:
                counting[residues_prot] = 1

            else:
                counting[residues_prot] += 1

            output_raw.write("%i,%s,%s\n"%(ts.frame,residues_prot,residues_lipids))


    for residue in counting.keys():
        counting[residue] = float(counting[residue])/nbFrame

        output_avg.write("%s,%f\n"%(residue,counting[residue]))

    for residue in psf_info.protresid:
        if str(residue) not in counting.keys():
            counting[str(residue)] = 0

    if pdb != 'None':
        common.Mapped(pdb,counting,segprot,outname,"_average_hydrophobic_contact_per_frame")


    output_avg.close()
    output_raw.close()