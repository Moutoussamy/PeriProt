#!/usr/bin/env python3

import common
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import hbonds


"""
SET OF FUNCTION TO EVALUATE HYDROPHOBIC CONTACT
"""


__author__ = "Emmanuel Edouard MOUTOUSSAMY"
__date__ = "2020/08"



def CollectHbondCandInfo(libpath):
    """
    Collect hbond donors and acceptors from lib/hbond_acceptor_donnor_list_lipids.dat
    :return: donors,acceptors,head,glycerol,phos lists
    """

    donors = [] #store donor atoms
    acceptors = [] # store acceptor atoms
    head = [] #store head atoms
    glycerol = [] # store glycerol atoms
    phos = [] #store phosphate atoms

    with open("%s/hbond_acceptor_donnor_list_lipids.dat"%libpath) as input_file:
        for line in input_file:

            if "#acceptor" in line : #Acceptor list

                acceptors = acceptors +  line.split()[1:]

            elif "#donor" in line: #Donor list
                donors = donors + line.split()[1:]

            elif "#head" in line:  #Head list
                head = head + line.split()[1:]


            elif "#phosphate" in line:  #phosphate list
                phos = phos + line.split()[1:]


            elif "#glycerol" in line:  #Glycerol list
                glycerol = glycerol + line.split()[1:]


    return donors,acceptors,head,glycerol,phos


def GetTimeSeries(psf,dcd,outname,libpath):
    """
    Get time series from the MDanalysis
    :param psf: PSF file
    :param dcd: trajectory file (DCD)
    :param outname: output name
    :return: times series results
    """

    donor, acceptor, head, glycerol, phos = CollectHbondCandInfo(libpath)

    universe = mda.Universe(psf,dcd)
    hbonds = mda.analysis.hbonds.HydrogenBondAnalysis(universe,"segid MEMB","segid PROA",donors= donor, acceptors = acceptor )
    hbonds.run()
    hbonds.generate_table()
    hbonds_data = pd.DataFrame.from_records(hbonds.table)
    hbonds_data.to_csv("%s_hbonds_analysis_raw_data.csv"%outname)

    return hbonds_data


def Gettime(hbonds_data):
    """
    get a dictionary with the 'time' as keys
    :param hbonds_data: hbond data from GetTimeSeries
    :return: dictionary with the 'time' as keys
    """

    times = hbonds_data["time"].tolist()
    times = list(dict.fromkeys(times))

    time_series = {}

    for time in times:
        time_series[str(time)] = []

    return time_series

def ReadRwData(outname,times):
    """
    collect data from the file write by GetTimeSeries
    :param outname: outname
    :param times: dictionary with the 'time' as keys
    :return: the time dictinary fill with the info from the file write by GetTimeSeries
    """

    with open("%s_hbonds_analysis_raw_data.csv"%outname) as inputfile:
        for line in inputfile:
            line = line.split(",")

            if line[4] in common.AMINO_ACIDS.keys():
                times[line[1][0:10]].append(line[5])


            if line[7] in common.AMINO_ACIDS.keys():
                times[line[1][0:10]].append(line[8])

    return times

def Occupancies(times_series_data,nbFrame):
    """
    Caculate the occupancies of the hbonds
    :param times_series_data: time series data from Gettime function
    :param nbFrame: number of frame in the trajectory
    :return: a dictionary; index = resid and value = occupancy
    """
    occupancy = {}

    for time in times_series_data.keys():
        resids = list(dict.fromkeys(times_series_data[time]))
        for resid in resids:
            if resid not in occupancy:
                occupancy[resid] = 1
            else:
                occupancy[resid] += 1


    for resid in occupancy.keys():
        occupancy[resid] = (float(occupancy[resid])/float(nbFrame))*100

    return occupancy

def write_results(occupancies,outname):
    """
    Write the final results in a file
    :param occupancies: outout of the function Occupancies(times_series_data,nbFrame)
    :param outname: output name
    :return:
    """
    output = open("%s_hbond_occupancies.csv"%outname,"w")
    output.write("#Resid,occupancy(%)\n")

    for resid in list(dict.fromkeys(occupancies.keys())):
        output.write("%s,%f\n"%(resid,occupancies[resid]))

    output.close()

def runHbond(psf,dcd,outname,resid_list,segprot,pdb,libpath):
    """
    Run Hbond calculation
    :param psf: PSF file
    :param dcd: DCD file
    :param outname: output name
    :param resid_list: list of resid in the protein
    :param segprot: segif of the proetin
    :param pdb: PDB file
    :return:
    """
    hbonds_data = GetTimeSeries(psf, dcd, outname,libpath)
    times = Gettime(hbonds_data)
    times_series_data = ReadRwData(outname,times)
    occupancies = Occupancies(times_series_data,len(times.keys()))
    write_results(occupancies,outname)

    for residue in resid_list:
        if str(residue) not in occupancies.keys():
            occupancies[str(residue)] = 0

    if pdb != "None":
        common.Mapped(pdb, occupancies, segprot, outname, "hbond_occupancies")

