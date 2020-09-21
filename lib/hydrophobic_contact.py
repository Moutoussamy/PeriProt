#!/usr/bin/env python3

import common
import  MDAnalysis as mda
from MDAnalysis.analysis import contacts

"""
SET OF FUNCTION TO EVALUATE HYDROPHOBIC CONTACTS 
"""


def readCandidates(Memb_res):
    """
    Read the atom list in the file: lib/hydrophobic_candidates_prot.dat
    :param Memb_res: list of lipids in the system
    :return: two dictionary (one for the membrane and one for the protein); index = resname
    and value = list of cand. atoms
    """

    prot_cand_atom = {}
    memb_cand_atom = {}

    with open("lib/hydrophobic_candidates_prot.dat") as inputfile:
        for line in inputfile:
            line = line.split()
            if len(line) == 2:
                prot_cand_atom[line[0]] = line[1].split(",")


    with open("lib/hydrophobic_candidates_lipid.dat") as inputfile:
        for line in inputfile:
            line = line.split()
            if line[0] in Memb_res:
                memb_cand_atom[line[0]] = line[1].split(",")

    return prot_cand_atom,memb_cand_atom


def HydroCandidateSelection(psf,segprot,segmemb,prot_cand_atom,memb_cand_atom):
    """
    Take the list of candidates atoms and select the corresponding atomid in the PSF
    :param psf: PSF file
    :param segprot: segid of the protein
    :param segmemb: segid of the membrane
    :param prot_cand_atom: candidates atoms names for the protein
    :param memb_cand_atom: candidates atoms names for the membrane
    :return:
    """

    atomid_cand_prot = []
    atomid_cand_memb = []


    with open(psf) as topology_file:

        for line in topology_file:
            line = line.split()
            if len(line) == 9:
                if line[1] == segprot:
                    if line[3] in prot_cand_atom.keys():
                        if line[4] in prot_cand_atom[line[3]]:
                            atomid_cand_prot.append(line[0])
                elif line[1] == segmemb:
                    if line[3] in memb_cand_atom.keys():
                        if line[4] in memb_cand_atom[line[3]]:
                            atomid_cand_memb.append(line[0])

    return  atomid_cand_prot, atomid_cand_memb



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



def EpureSelectedAtom(near_atoms,atomid_cand):
    """
    Take the list of candidates atom id and and keep only the one close to the bilayer
    :param near_atoms: atom near the bilayer
    :param atomid_cand: list of candidates atom id
    :return: list of candidates atom id with only atoms close to the bilayer/protein interface
    """

    selected = []

    for atom in atomid_cand:
        if int(atom) in near_atoms:
            selected.append(atom)

    return selected


def GetCandidate(psf_info,close_lipid,close_prot,atomid_cand_prot,atomid_cand_memb):
    """
    Obtain candiate atom for hydrophobic contact analysis
    :param psf_info: class conating infos about the system (topology from ./lib/common.py)
    :param close_lipid: lipid close to the bilayer
    :param close_prot: amino acids close to the bilayer
    :return: list of candidates for the protein and lkist of candidates for the membrane
    """

    prot_atom_near = GetNearAtom(close_prot,psf_info.NbAtomPerProtRes)
    prot = EpureSelectedAtom(prot_atom_near,atomid_cand_prot )
    prot = sorted(prot)

    memb_atom_near = GetNearAtom(close_lipid,psf_info.NbAtomPerLipRes)
    memb = EpureSelectedAtom(memb_atom_near,atomid_cand_memb)
    memb = sorted(memb)

    return prot,memb



def RunHydroAnalysis(psf, dcd, psf_info, close_lipid,close_prot,outname,segprot,segmemb,pdb\
                     ,first_frame,last_frame,skip):
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

    prot_cand_atom, memb_cand_atom = readCandidates(psf_info.lipids)
    atomid_cand_prot, atomid_cand_memb = HydroCandidateSelection(psf,segprot, segmemb, prot_cand_atom, memb_cand_atom)
    prot_candidates, memb_candidates = GetCandidate(psf_info,close_lipid,close_prot,atomid_cand_prot, atomid_cand_memb)
    prot_candidates = common.convertIntToStr(prot_candidates)
    memb_candidates = common.convertIntToStr(memb_candidates)



    univers = mda.Universe(psf, dcd)

    counting = {}
    nbFrame = 0

    output_avg = open("%s_average_hydrophobic_contact_per_frame.csv"%outname,"w")
    output_avg.write("#Residue,AVG_nb_of_contact\n")

    output_raw = open("%s_hydrophobic_contact_raw_data.csv"%outname,"w")
    output_raw .write("#Frame,ProtResidue,MembRes\n")

    frame_to_read = first_frame

    for ts in  univers.trajectory:
        if ts.frame == frame_to_read and frame_to_read < last_frame:
            nbFrame += 1
            sel_memb = "index %s" % (" ".join(memb_candidates))
            memb = univers.select_atoms(sel_memb)
            sel_prot = "index %s" % (" ".join(prot_candidates))
            prot = univers.select_atoms(sel_prot)
            dist_mat = mda.analysis.distances.distance_array(prot.positions, memb.positions)

            conctacts =[]

            #Contact will be fill with list of 2 elements: [atomid in prot,atomid in memb]
            for i in range(dist_mat.shape[0]):
                for j in range(dist_mat.shape[1]):
                    if dist_mat[i,j] < 3:
                        conctacts.append([prot_candidates[i],memb_candidates[j]])


            #parse contact and write the raw data
            for i in range(len(conctacts)):
                residues_prot = common.obtainResInfo(int(conctacts[i][0]), psf_info.NbAtomPerProtRes)
                residues_lipids =  common.obtainResInfo(int(conctacts[i][1]), psf_info.NbAtomPerLipRes)

                if residues_prot not in counting:
                    counting[residues_prot] = 1

                else:
                    counting[residues_prot] += 1

                output_raw.write("%i,%s,%s\n"%(ts.frame,residues_prot,residues_lipids))

            frame_to_read += skip



    #Get the average
    for residue in counting.keys():
        counting[residue] = float(counting[residue])/nbFrame

        output_avg.write("%s,%f\n"%(residue,counting[residue]))

    for residue in psf_info.protresid:
        if str(residue) not in counting.keys():
            counting[str(residue)] = 0

    if pdb != 'None':
        common.Mapped(pdb,counting,segprot,outname,"average_hydrophobic_contact_per_frame")


    output_avg.close()
    output_raw.close()