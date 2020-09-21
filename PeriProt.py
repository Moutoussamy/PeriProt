#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import  MDAnalysis as mda
sys.path.insert(0, './lib/')
import hbond as hb
import common
import hydrophobic_contact as hydro
import depth as dp
import ProtMembDist as pmd
import cationpi as catpi
import Macrodipole
import BodyGuard

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
    parser.add_argument('-dist', action='store_true', help="Prot. - Memb distance")
    parser.add_argument('-segmemb', type=str, default='MEMB', help="segid for membrane")
    parser.add_argument('-segprot', type=str, default='PROA', help="segid for protein")
    parser.add_argument('-edp', action='store_true', help="electro density profile")
    parser.add_argument('-first', type=int, default=0, help="first frame to read")
    parser.add_argument('-last', type=int, default=10000000, help="last frame to read")
    parser.add_argument('-skip', type=int, default=1, help="first frame to read")
    parser.add_argument('-out', type=str, default="periprot", help="output name")
    parser.add_argument('-pdb', type=str, default="None", help="PDB file")
    parser.add_argument('-mdipole', action='store_true', help="Macrodipole calculation")

    args = parser.parse_args()

    return args



def done():
    print("...done")


if __name__ == '__main__':

    arguments = GetArgs()
    print(arguments)

    BodyGuard.checkSegID(arguments.top, arguments.segprot, arguments.segmemb) #Did you gice the right segids ?
    psf_info = common.topology(arguments.top, arguments.segmemb, arguments.segprot)


    if arguments.hydro or arguments.hbond or arguments.catpi :
        print("Look for lipids close to the protein...")
        #close_lipids, close_amino_acids = common.GetClosePartner(arguments.top,arguments.traj, arguments.segprot, arguments.segmemb, psf_info)

        close_lipids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 28, 29,\
                        31, 33, 40, 41, 103, 104, 113, 118, 120, 121, 132, 137, 36, 38, 39, 109, 130, 133, 24, 106,\
                        112, 134, 142, 50, 131, 135, 26, 235, 98, 94, 140, 155, 244, 115, 146, 42, 93, 25, 32, 250,\
                        49, 52, 45, 144, 158, 47, 129, 164]

        close_amino_acids = [32, 33, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 69,\
                             71, 73, 74, 75, 76, 77, 79, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 93, 115, 116, 117,\
                             118, 119, 120, 121, 122, 123, 124, 125, 180, 181, 199, 200, 201, 202, 203, 204, 205, 206,\
                             207, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251,\
                             252, 275, 277, 279, 280, 35, 36, 97, 126, 182, 234, 274, 281, 34, 57, 67, 163, 198, 209,\
                             232, 58, 72, 254, 276, 278, 178, 179, 255, 258, 177, 208]
        done()


    BodyGuard.PriminilaryCheck(arguments)

    ### HBOND ANALYSIS
    if arguments.hbond:
        print("Hbond network analysis...")
        hb.runHbond(arguments.top,arguments.traj,arguments.out,psf_info.protresid,arguments.segprot,arguments.pdb)
        done()

    ### HYDROPHOBIC CONTACT ANALYSIS
    if arguments.hydro:
        print("Hydrophobic contact calculations...")
        #prot_candidates, memb_candidates = hydro.HydroCandidateSelection(psf_info)
        hydro.RunHydroAnalysis(arguments.top,arguments.traj,psf_info,close_lipids,close_amino_acids,arguments.out,\
                               arguments.segprot,arguments.segmemb,arguments.pdb)
        done()



    ### CATION-PI ANALYSIS

    if arguments.catpi:
        print("Cation-Pi interaction calculation...")

        catpi.RunAromatic(arguments.top,arguments.traj,psf_info,arguments.segmemb,close_lipids, close_amino_acids,\
                          arguments.out,arguments.pdb,arguments.segprot)
        done()

    ### DEPTH
    if arguments.depth:
        print("Depth of anchoring calculations...")
        residues_list = common.get_residues_list(arguments.top, arguments.segprot)
        dp.RunDepthOfAnchoring(arguments.top,arguments.traj,residues_list,arguments.segmemb,arguments.segprot,\
                               arguments.first,arguments.last,arguments.skip,arguments.out,arguments.pdb)
        done()


    ### DISTANCE

    if arguments.dist:
        pmd.ComputeDistance(arguments.top,arguments.traj,arguments.segprot,arguments.segmemb,arguments.out,psf_info)
        done()


    if arguments.mdipole:
        Macrodipole.RunMacrodipoleCalculation(arguments.pdb,psf_info.protCharge,arguments.out,arguments.segprot)
