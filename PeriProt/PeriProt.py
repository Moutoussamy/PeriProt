#!/usr/bin/env python3

import sys
import argparse
import os

LIBPATH = '%s/../lib/'%(os.path.dirname(sys.argv[0])) #Path of Periprot library
sys.path.insert(0,LIBPATH)
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
1) Hbond Analysis: calculate the Hydrogen bond event between the proitein and the bilayer
2) Hydrophobic contact analysis: evaluation of the hydrophobic contact between the protein and the bilayer
3) Cation-pi int. analysis: evaluate cation interaction between the PMP and the membrane
4) Anchoring depth: calculate the depth of anchoring of each residues (avg. during the simulation)
5) Prot./Memb. distance: calculate the distance between the protein and the membrane during the simulation
6) Macrodipole: calculate the macrodipole of the protein


Usage:
  -h, --help        show this help message and exit
  -top TOP          psf file only
  -traj TRAJ        trajectories (DCD format)
  -hydro            Hydrophobic contact analysis
  -hbond            hbond analysis
  -catpi            cation pi analysis
  -depth            depth of anchoring analysis
  -dist             Prot. - Memb distance
  -segmemb SEGMEMB  segid for membrane
  -segprot SEGPROT  segid for protein
  -edp              electro density profile
  -first FIRST      first frame to read
  -last LAST        last frame to read
  -skip SKIP        first frame to read
  -out OUT          output name
  -pdb PDB          PDB file
  

1) Hbond Analysis:

python PeriProt.py -top mydata.psf -traj mydata.dcd -pdb mydata.pdb  -out example -hbond

if -pdb flag is mentioned, the b-factor column of the protein will be replace by the occupancies of each residue

2) Hydrophobic contact analysis:

python PeriProt.py -top mydata.psf -traj mydata.dcd -pdb mydata.pdb  -out example -hydro

if -pdb flag is mentioned, the b-factor column of the protein will be replace by the avg nb. of hydrophobic contact
of each residue with the bilayer

3) Cation-pi int. analysis:

python PeriProt.py -top mydata.psf -traj mydata.dcd -pdb mydata.pdb  -out example -catpi

if -pdb flag is mentioned, the b-factor column of the protein will be replace by the occupancies of each TYR or TRP

4) Anchoring depth:

python PeriProt.py -top mydata.psf -traj mydata.dcd -pdb mydata.pdb  -out example -depth
if -pdb flag is mentioned, the b-factor column of the protein will be replace by the occupancies of each TYR or TRP


5) Prot./Memb. distance:

python PeriProt.py -top mydata.psf -traj mydata.dcd -out example -dist

three possible choice of distance:
    Which distance ?:
    1: Prot (COM) - Memb (COM) distance
    2: Prot (COM) - Phosphate Plane (COM) distance
    3: Custom distance between two atoms (one in the prot. and one on the memb.)

6) Macrodipole:

python PeriProt.py -top mydata.psf -pdb mydata.pdb  -out example -mdipole

"""

__author__ = "Emmanuel Edouard MOUTOUSSAMY"
__date__ = "2020/08"

def GetArgs():
    """
    Get all the necessary arguments
    :return: a list with all arguments
    """

    parser = argparse.ArgumentParser(description='PeriProt Analysis:')

    parser.add_argument('-top', default='None', help= "psf file only") #toplogy file (PSF)
    parser.add_argument('-traj', help="trajectories (DCD format)") #Trajectory file (DCD)
    parser.add_argument('-hydro',action='store_true', help="Hydrophobic contact analysis")
    parser.add_argument('-hbond',action='store_true', help="hbond analysis")
    parser.add_argument('-catpi', action='store_true', help="cation pi analysis")
    parser.add_argument('-depth',action='store_true', help="depth of anchoring analysis")
    parser.add_argument('-dist', action='store_true', help="Prot. - Memb distance")
    parser.add_argument('-segmemb', type=str, default='MEMB', help="segid for membrane")
    parser.add_argument('-segprot', type=str, default='PROA', help="segid for protein")
    parser.add_argument('-first', type=int, default=0, help="first frame to read")
    parser.add_argument('-last', type=int, default=10000000, help="last frame to read")
    parser.add_argument('-skip', type=int, default=1, help="first frame to read")
    parser.add_argument('-out', type=str, default="periprot", help="output name")
    parser.add_argument('-pdb', type=str, default="None", help="PDB file")
    parser.add_argument('-mdipole',action='store_true', help="Macrodipole calculation")
    parser.add_argument('-tres', action='store_true', help="target  some residues")

    args = parser.parse_args()

    return args



def done():
    """
    Print "...done" when the job is done
    :return: none
    """
    print("...done")

def write_close_partner(close_lipids, close_amino_acids,logfile):
    """

    :param close_lipids:
    :param close_amino_acids:
    :param logfile:
    :return:
    """
    logfile.write("lipids close to the protein: %s\n" % (" ".join(str(close_lipids))))
    logfile.write("residues close to the bilayer: %s\n" % (" ".join(str(close_amino_acids))))

def DoYouWantToTargetRes(protresid,tres_flag):
    """
    check if targeting residues is requested
    :param protresid: protein resid
    :param tres_flag: tres flag
    :return:
    """

    if tres_flag:
        res_list = common.GetDesireResidueList()
        return res_list

    else:
        return protresid




if __name__ == '__main__':

    arguments = GetArgs()
    logfile = open("%s.log"%arguments.out,"w")
    logfile.write("command line: %s\n"%(" ".join(sys.argv)))

    BodyGuard.PriminilaryCheck(arguments) #check if input are correct
    BodyGuard.checkSegID(arguments.top, arguments.segprot, arguments.segmemb) #Did you gice the right segids ?

    psf_info = common.topology(arguments.top, arguments.segmemb, arguments.segprot)

    residues_list = DoYouWantToTargetRes(psf_info.protresid,arguments.tres)
    print(residues_list)

    ### RESIDUES AT THE PROT/MEMB INTERFACE
    if arguments.hydro or arguments.hbond or arguments.catpi :
        print("Look for lipids close to the protein...")
        close_lipids, close_amino_acids = common.GetClosePartner(arguments.top,arguments.traj, arguments.segprot,\
                                                                 arguments.segmemb, psf_info,arguments.first,\
                                                                 arguments.last,arguments.skip)

        write_close_partner(close_lipids, close_amino_acids,logfile)
        done()



    ### HBOND ANALYSIS
    if arguments.hbond:
        print("Hbond network analysis...")

        hb.runHbond(arguments.top,arguments.traj,arguments.out,residues_list,arguments.segprot,\
                    arguments.pdb,LIBPATH)
        done()

    ### HYDROPHOBIC CONTACT ANALYSIS
    if arguments.hydro:
        print("Hydrophobic contact calculations...")
        hydro.RunHydroAnalysis(arguments.top,arguments.traj,psf_info,close_lipids,close_amino_acids,arguments.out,\
                               arguments.segprot,arguments.segmemb,arguments.pdb,arguments.first,arguments.last,\
                               arguments.skip,LIBPATH)
        done()



    ### CATION-PI ANALYSIS

    if arguments.catpi:
        print("Cation-Pi interaction calculation...")

        catpi.RunAromatic(arguments.top,arguments.traj,psf_info,arguments.segmemb,close_lipids, close_amino_acids,\
                          arguments.out,arguments.pdb,arguments.segprot,arguments.first,arguments.last,arguments.skip)
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


    ### MACRODIPOLE
    if arguments.mdipole:
        print("Macrodipole calculation...")
        Macrodipole.RunMacrodipoleCalculation(arguments.pdb,psf_info.protCharge,arguments.out,arguments.segprot)
        done()


    logfile.close()