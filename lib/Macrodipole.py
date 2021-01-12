#!/usr/bin/env python
# ~*~ coding:utf-8 ~*~


__author__ = "Emmanuel Edouard MOUTOUSSAMY"
__version__  = "1.0.0"
__copyright__ = "copyleft"
__dependencies__ = "os, sys, pandas, numpy, matplotlib"
__date__ = "2019/10"


import numpy as np
import pandas as pd


"""
MACRODIPOLE OF THE PROTEIN CALCULATION 
IT IS INTERESITING TO CALCULATE MACRODIPOLE FOR PMP
"""

__author__ = "Emmanuel Edouard MOUTOUSSAMY"
__date__ = "2020/08"


def WritePDB2PQR(DicoCharge,pdb,outname,segidprot):
    """
    Write a kind of PQRFILE
    :param DicoCharge: dico with key = RES_ATOMNAME: VALUE: atom charge.
    :param pdb: the pdb file
    :param pdbname: the name of PDB without extension
    :return:  write a pqr file
    """
    output = open("%s_charge_info.pqr"%(outname),"w")
    output.write("#ATOMNUM,RES,ATOMNAME,X,Y,Z,CHARGE\n")

    with open(pdb) as PDBFile:
        for line in PDBFile:
            if line[0:4] == "ATOM":
                if segidprot in line:
                    AtomeNum = line[6:11].replace(" ", "")
                    Residue = line[17:20].replace(" ", "")
                    AtomName = line[12:16].replace(" ","")
                    AtomCode = "%s_%s"%(Residue,AtomName)
                    charge = DicoCharge[AtomCode]
                    CoordX = float(line[30:38])
                    CoordY = float(line[38:46])
                    CoordZ = float(line[46:54])
                    output.write("{0},{1},{2},{3},{4},{5},{6}\n".format(AtomeNum,Residue,AtomName,CoordX,CoordY,CoordZ,charge))


def CalCenterofCharge(PqrFile):
    """
    Calculate the center of charge
    :param PqrFile: pqr file
    :return: a tuple with the center of mass
    """
    cx, cy, cz, SumAbsCharge = 0, 0, 0, 0

    for index, row in PqrFile.iterrows():
        cx += row["X"] * abs(row["CHARGE"])
        cy += row["Y"] * abs(row["CHARGE"])
        cz += row["Z"] * abs(row["CHARGE"])
        SumAbsCharge += abs(row["CHARGE"])

    cx = cx / SumAbsCharge
    cy = cy / SumAbsCharge
    cz = cz / SumAbsCharge

    return (cx,cy,cz)


def CalDipoleMoment(CenterOfCharge,PqrFile):
    """
    Calculate the dipole moment
    :param CenterOfCharge: center of charge
    :param PqrFile: pqr file from WritePDB2PQR function
    :return: dipole moment magnitude and dipole vecto in x,y anz
    """
    dp_x, dp_y, dp_z = 0, 0, 0

    for index, row in PqrFile.iterrows():
        dp_x += ( row["X"] - CenterOfCharge[0]) * row["CHARGE"]
        dp_y += (row["Y"] - CenterOfCharge[1]) * row["CHARGE"]
        dp_z += (row["Z"] - CenterOfCharge[2]) * row["CHARGE"]

    dpmag = np.sqrt(dp_x*dp_x + dp_y*dp_y + dp_z*dp_z) #dipole moment magnitude
    dpmag_SI = dpmag*1.6022e-19*1e-10
    dpmag_Deb = dpmag_SI/3.33564e-30

    return dpmag_Deb,dp_x,dp_y,dp_z


def DrawArrow(outname,dpmag_Deb, dp_x, dp_y, dp_z,CenterOfCharge):
    """
    Create a .bild file to visualize the macrodipole on chimera
    the dipole start from the center of charge of the protein

    :param outname:  outpout name
    :param dpmag_Deb: dipole moment magnitude
    :param dp_x: final point on X of the dipole vector
    :param dp_y: final point on Y of the dipole vector
    :param dp_z: final point on Z of the dipole vector
    :param CenterOfCharge:  center of charge
    :return:
    """
    bildname = "%s_macrodipole_arrow.bild" % outname
    output = open(bildname, "w")
    output.write("!  Dipole moment magnitude %8.3f D\n" % dpmag_Deb)
    output.write(".color blue\n")
    output.write(".arrow  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  1 4 0.9\n" % (CenterOfCharge[0], CenterOfCharge[1]\
    , CenterOfCharge[2], dp_x + CenterOfCharge[0], dp_y + CenterOfCharge[1], dp_z + CenterOfCharge[2]))
    output.close()



def RunMacrodipoleCalculation(pdb,charge_dico,outname,segidprot):
    """
    Run the macrodipole calculation
    :param pdb: PSD file
    :param charge_dico: dico of charge for topology structure
    :param outname: outputname
    :param segidprot: segid of the protein
    :return: none
    """
    WritePDB2PQR(charge_dico,pdb,outname,segidprot)
    PqrFile = pd.read_csv("%s_charge_info.pqr"%(outname))
    CenterOfCharge = CalCenterofCharge(PqrFile)
    dpmag_Deb, dp_x, dp_y, dp_z = CalDipoleMoment(CenterOfCharge,PqrFile)
    DrawArrow(outname,dpmag_Deb, dp_x, dp_y, dp_z,CenterOfCharge)
