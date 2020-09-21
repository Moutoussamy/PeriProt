#!/usr/bin/env python3

import sys

"""
BodyGuard will check input to kick you out if you did something wrong.
"""

__author__ = "Emmanuel Edouard MOUTOUSSAMY"
__version__  = "1.0.0"
__date__ = "2020/08"
__copyright__ = "CC_by_SA"
__dependencies__ = "Numpy,MDAnalysis and argparse"




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
        pdb: PDB file
        """)

def ErroMessageSegID(segid,flag):

    message = \
        "\nsegid: {0} not found in the PSF ! You probably used an other segid. Please mention it with the flag: {1}\n"\
        .format(segid,flag)

    sys.exit(message)



def check_args(psf,dcd,pdb):
    """
    Check arguments for PSF and DCD files
    :param psf: topology file (PSF or PDB)
    :param dcd: Trajectories file (DCD)
    :return: none
    """
    check_extension(psf,[".psf"]) #check extension of topology file
    check_extension(dcd, [".dcd"]) #check traj file
    check_extension(pdb, [".pdb"])  # check traj file


def checkSegID(psf,segprot,segmemb):

    segids = []

    flag = 1

    with open(psf) as inputfile:
        for line in inputfile:
            if "!NBOND" in line:
                flag = 0

            if flag == 1:
                line = line.split()
                if len(line) == 9:
                    if line[1] not in segids:
                        segids.append(line[1])


    if segmemb not in segids:
        if segprot not in segids:
            ErroMessageSegID("%s and %s"%(segmemb,segprot), "-segmemb and  -segprot")
        else:
            ErroMessageSegID(segmemb, "-segmemb")


    if segprot not in segids:
        ErroMessageSegID(segprot, "-segprot")


def PriminilaryCheck(arguments):
    if arguments.hydro or arguments.hbond or arguments.catpi or arguments.depth or arguments.dist:
        if arguments.traj == None:
            sys.exit("\nNo DCD file given. Please use -traj flag to specify of dcd file\n")


def CheckinghydrophobicContact(lipids_list,our_list_of_lipids):

    for lipid in lipids_list:
        if lipids_list not in our_list_of_lipids:
            sys.exit("""\nlipid: %s not reconized. Please fill lib/hydrophobic_candidates_lipid.dat manually or use
            the script lib/ReadToppar.py\n""")

