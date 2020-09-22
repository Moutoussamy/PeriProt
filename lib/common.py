#!/usr/bin/env python3

import numpy as np
import  MDAnalysis as mda

'''
SET OF COMMON FUNCTIONS FOR THE DIFFERENT ANALYSIS
'''


__author__ = "Emmanuel Edouard MOUTOUSSAMY"
__date__ = "2020/08"



AMINO_ACIDS = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}



#CLASS FOR TOPOLOGY: PSF FILE
class topology():
    """
    Class containg information about the PSF file:
    - Sequence
    - Membrane resid
    - Protein resid
    - Membrane index [start,last]
    - Protein index [start,last]
    - Charge {index:charge}
    - Mass {index:mass}
    """
    def __init__ (self,PSF, segidMEMB, segidPROT):
        """
        Inisitialize the structure
        :param PSF: PSF file
        :param segidMEMB: Membrane segid
        :param segidPROT: Protein segid
        """

        self.membresid = []
        self.lipids = []
        self.protresid = []
        self.protindex = [0,0]
        self.membindex = [0,0]
        self.protCharge = {}
        self.tyrosines = []
        self.tryptophanes = []
        self.NbAtomPerLipRes = {}
        self.NbAtomPerProtRes = {}

        self.fillTopoStruct(PSF, segidMEMB, segidPROT)

    def FillProtInfos(self,psfline):
        """
        Fill information regarding the protein
        :param psfline: a line of the PSF file
        :return:
        """

        if self.protindex[0] == 0:
            self.protindex[0] = int(psfline[0])

        if int(psfline[0]) < self.protindex[0]: #find the first atom number of the protein
            self.protindex[0] = int(psfline[0])

        if int(psfline[0])  > self.protindex[1]: #find the last atom number of the protein
            self.protindex[1] = int(psfline[0])

        key = "%s_%s"%(psfline[3],psfline[4])

        self.protCharge[key] = float(psfline[6]) #Collect the charge


        if int(psfline[2]) not in self.protresid: #Collect the resid
            self.protresid.append(int(psfline[2]))

        if int(psfline[2]) not in self.NbAtomPerProtRes.keys():
            self.NbAtomPerProtRes[int(psfline[2])] = [int(psfline[0]),int(psfline[0])]

        else:
            if int(psfline[0]) > self.NbAtomPerProtRes[int(psfline[2])][1]:
                self.NbAtomPerProtRes[int(psfline[2])][1] = int(psfline[0])

        if psfline[3] == "TYR":
            if psfline[2] not in self.tyrosines:
                self.tyrosines.append(psfline[2])

        elif psfline[3] == "TRP":
            if psfline[2] not in self.tryptophanes:
                self.tryptophanes.append(psfline[2])


    def FillMembInfos(self,psfline):
        """
        Fill information regarding the membrane
        :param psfline: a line of the PSF file
        :return:
        """

        if self.membindex[0] == 0:
            self.membindex[0] = int(psfline[0])

        if int(psfline[0]) < self.membindex[0]: #find the first atom number of the membrane
            self.membindex[0] = int(psfline[0])

        if int(psfline[0])  > self.membindex[1]: #find the last atom number of the membrane
            self.membindex[1] = int(psfline[0])

        key = "%s_%s_%s_%s"%(psfline[3],psfline[2],psfline[0],psfline[4])


        if int(psfline[2]) not in self.membresid: #Collect the resid
            self.membresid.append(int(psfline[2]))

        if int(psfline[2]) not in self.NbAtomPerLipRes.keys():
            self.NbAtomPerLipRes[int(psfline[2])] = [int(psfline[0]),int(psfline[0])]

        else:
            if int(psfline[0]) > self.NbAtomPerLipRes[int(psfline[2])][1]:
                self.NbAtomPerLipRes[int(psfline[2])][1] = int(psfline[0])

        if psfline[3] not in self.lipids:
            self.lipids.append(psfline[3])


    def Collectinfos(self,psf,segidMEMB,segidPROT):
        """
        Collect information from PSF file
        :param psf: PSF file
        :param segidMEMB: membrane segid
        :param segidPROT: protein segid
        :return:
        """

        with open(psf) as topology_file:
            for line in topology_file:
                line = line.split()
                if len(line) == 9:
                    if line[1] == segidPROT:
                        self.FillProtInfos(line)

                    elif line[1] == segidMEMB:
                        self.FillMembInfos(line)

    def fillTopoStruct(self,PSF,segidMEMB,segidPROT):
        """
        Fill the structure TOPOLOGY
        :param psf: PSF file
        :param segidMEMB: membrane segid
        :param segidPROT: protein segid
        :return:
        """
        self.Collectinfos(PSF,segidMEMB,segidPROT)



#FUNCTIONS

def Mapped(pdb,dico,segid_prot,outname,analysis):
    """
    Map the dist. on a pdb file
    :param pdb: a PDB file of a protein/membrane complex
    :param dico: dico containing data to change
    :param segid_prot: segid pf the protein
    :param outname: name of output
    :return: write a pdb file with the new dada as the B-factor
    """

    pdbout = "%s_%s.pdb"%(outname,analysis)

    output = open(pdbout,"w")

    with open(pdb) as inputfile:
        for line in inputfile:
            if "ATOM" in line:
                if segid_prot in line:
                    resid = str(int(line[22:26]))
                    extend = line[46:].replace(line[60:66],"%6.2f"%dico[resid])
                    line = "%s%s"%(line[0:46],extend)
                    output.write(line)

def get_residues_list(psf_file,segidPROT):
    """
    Obtain the list of residues that compose the protein
    :param psf_file: PSF file
    :param segidPROT: segid of the protein in the PSF filel
    :return:
    """

    residues_list = []
    with open(psf_file) as inputfile:
        for line in inputfile:
            line = line.split()
            if len(line) == 9:
                if line[1] == segidPROT:
                    if line[4].replace(" ","") == 'CA':
                        residues_list.append(line[2])

    return residues_list


def get_lipid_resid_list(psf_file,segidMEMB):
    """
    Obtain the list of residues that compose the protein
    :param psf_file: PSF file
    :param segidPROT: segid of the protein in the PSF filel
    :return:
    """

    residues_list = []
    with open(psf_file) as inputfile:
        for line in inputfile:
            line = line.split()
            if len(line) == 9:
                if line[1] == segidMEMB:
                    if line[2] not in residues_list:
                        residues_list.append(line[2])

    return residues_list



def ParseDistanceMatrix(dist_mat,psf_info,closed_lipids,close_amino_acids):
    """
    Parse the distance matric in order to select residues and lipids within 7 A
    :param dist_mat:
    :param psf_info:
    :param closed_lipids:
    :param close_amino_acids:
    :return:
    """
    start_memb = psf_info.membindex[0]

    for resid in psf_info.NbAtomPerLipRes.keys():
        index_low = psf_info.NbAtomPerLipRes[resid][0] - start_memb
        index_up = psf_info.NbAtomPerLipRes[resid][1] - start_memb


        if np.amin(dist_mat[:,index_low:index_up]) < 7:
            if resid not in closed_lipids:
                closed_lipids.append(resid)


    for resid in psf_info.NbAtomPerProtRes.keys():
        index_low = psf_info.NbAtomPerProtRes[resid][0]
        index_up = psf_info.NbAtomPerProtRes[resid][1]


        if np.amin(dist_mat[index_low:index_up,:]) < 7:
            if resid not in close_amino_acids:
                close_amino_acids.append(resid)


    return closed_lipids,close_amino_acids


def GetClosePartner(psf,dcd,segprot,segmemb,psf_info):
    """
    Get the lipid close to the prot and the aminoacid close to the bilayer during the simulation
    :param psf: PSF file
    :param dcd: trajectory
    :param segprot: segid of the protein
    :param segmemb: segid of the membrane
    :param psf_info: structure containing info about the system
    :return:
    """
    univers = mda.Universe(psf, dcd)

    closed_lipids = []
    close_amino_acids = []

    for ts in univers.trajectory:
        prot = univers.select_atoms("segid %s"%segprot)
        memb = univers.select_atoms("segid %s"%segmemb)
        dist_mat = mda.analysis.distances.distance_array(prot.positions,memb.positions)
        closed_lipids, close_amino_acids = ParseDistanceMatrix(dist_mat,psf_info,closed_lipids,close_amino_acids)


    return closed_lipids,close_amino_acids


def convertIntToStr(list_int):
    """
    Convert a list of int to a list os str
    :param list_int: list of int
    :return: a list of str
    """
    list_str = []
    for i in list_int:
        list_str.append(str(i))

    return list_str

def obtainResInfo(atomid, NbAtomPerRes):
    """
    Obetain the resid based on atom id
    :param atomid: atom id
    :param NbAtomPerRes:  dico made by topology structure (common)
    :return:
    """

    for residue in NbAtomPerRes.keys():
        if NbAtomPerRes[residue][0] < atomid < NbAtomPerRes[residue][1]:
            return str(residue)
            break