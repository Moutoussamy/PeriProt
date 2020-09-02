#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import  MDAnalysis as mda
sys.path.insert(0, './lib/')
import hbond as hb
import hydrophobic_contact as hydro
import depth as dp

"""
PeriProt 1.0
PeriProt is a software allowing to analysis Peripheral Membrane Protein (PMP) interaction with bilayer
during Molecular Dynamics (MD) simulation.

In this version the available tools are:
- Hbond Analysis: calculate the Hydrogen bond event between the proitein and the bilayer
- Hydrophobic contact analysis: evaluation of the hydrophobic contact between the protein and the bilayer
- Cation-pi int. analysis: evaluate cation interaction between the PMP and the membrane
- Anchoring depth: calculate the depth of anchoring of each residues (avg. during the simulation)
- electron density profile: calculate the electron density profile (avg. during the simulation)
"""

__author__ = "Emmanuel Edouard MOUTOUSSAMY"
__version__  = "1.0.0"
__date__ = "2020/08"
__copyright__ = "CC_by_SA"
__dependencies__ = "Numpy,MDAnalysis and argparse"


def GetArgs():
    """
    Get all the necessary arguments
    :return: args: a list:
    """

    parser = argparse.ArgumentParser(description='PeriProt Analysis:')

    parser.add_argument('-top', help= "psf file only")
    parser.add_argument('-traj', help="trajectories (DCD format)")
    parser.add_argument('-hydro',action='store_true', help="Hydrophobic contact analysis")
    parser.add_argument('-hbond', action='store_true', help="hbond analysis")
    parser.add_argument('-hbcand', type=str, default='None', help="hbond candidates")
    parser.add_argument('-catpi', action='store_true', help="cation pi analysis")
    parser.add_argument('-depth', action='store_true', help="depth of anchoring analysis")
    parser.add_argument('-segmemb', type=str, default='MEMB', help="segid for membrane")
    parser.add_argument('-segprot', type=str, default='PROA', help="segid for protein")
    parser.add_argument('-edp', action='store_true', help="electro density profile")
    parser.add_argument('-first', type=int, default=0, help="first frame to read")
    parser.add_argument('-last', type=int, default=10000000, help="last frame to read")
    parser.add_argument('-skip', type=int, default=1, help="first frame to read")

    args = parser.parse_args()

    return args

def check_extension(file,extensions):
    """
    Check files extensions
    :param file: argument
    :param extensions: expected extension
    :return: error if wrong extension
    """
    if file[-4:] not in extensions: # '.psf' extension
        sys.exit(
        """ERROR: argument not reconized
        REMINDER:
        top: only PSF ! NO PDB !s
        traj: only dcd file
        """)

def check_args(psf,dcd):
    """
    Check arguments for PSF and DCD files
    :param psf: topology file (PSF or PDB)
    :param dcd: Trajectories file (DCD)
    :return: none
    """
    check_extension(psf,[".psf"]) #check extension of topology file
    check_extension(dcd, [".dcd"]) #check traj file


def get_selection(hbcand,segid_protein):
    candidates = []
    resid = []
    with open(hbcand) as inputfile:
        for line in inputfile:
            if "#" not in line:
                line = line.split()
                if line != []:
                    candidates.append("%s_%s"%(line[0],line[1]))
                    resid.append(line[1])

    selection = "segid %s and (resid %s)"%(segid_protein," ".join(resid))

    return candidates,selection

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




def runHbonds(a):
    print("k")



if __name__ == '__main__':

    arguments = GetArgs()
    print(arguments)
    check_args(arguments.top,arguments.traj)

    ### HBOND ANALYSIS
    if arguments.hbond:
        print("Hydrogen bonds Analysis...")
        if arguments.hbcand != "None":
            residues_selection,selection = get_selection(arguments.hbcand,arguments.segprot)
            runHbonds(selection)
        else:
            selection = arguments.segprot
            runHbonds(selection)

    ### HYDROPHOBIC CONTACT ANALYSIS
    if arguments.hydro:
        print("Hydrophobic contact calculations...")
        residues_list = get_residues_list(arguments.top, arguments.segprot)
        prot_candidates, memb_candidates = hydro.HydroCandidateSelection(arguments.segmemb,\
                                                                         arguments.segprot,arguments.top)

        hydro.RunHydroAnalysis(arguments.segmemb,arguments.top,arguments.traj,\
                               residues_list,prot_candidates, memb_candidates)


    ### CATION-PI ANALYSIS

    ### DEPTH
    if arguments.depth:
        print("Depth of anchoring calculations...")
        residues_list = get_residues_list(arguments.top, arguments.segprot)
        dp.RunDepthOfAnchoring(arguments.top,arguments.traj,residues_list,arguments.segmemb,arguments.first,\
                               arguments.last,arguments.skip)

