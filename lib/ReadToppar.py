#!/usr/bin/env python3

import sys

'''
READ TOPOLOGY FROM CHARMM AND MAKE A LIST OF POTENTIAL ATOM CANDIDATES FOR HYDROPHOBIC CONTACT CALCULATION
'''


__author__ = "Emmanuel Edouard MOUTOUSSAMY"
__date__ = "2020/08"



def ReadTopparCharmm(toppar):
    """
    read the topology and make list of potential atom candidates for hydrophobic contact calculation
    :param toppar: topology file from CHARMM FF
    :return: a list of hydrophobic candidate
    """
    residues = {}
    group = []

    with open(toppar) as inputfile:
        for line in inputfile:

            if line[0:4] == "RESI":
                if group != []:
                    flag = 0
                    for atom in group:
                        if atom[0] not in ["C","H"]:
                            flag = 1

                    if flag == 0:
                        residues[resi] = residues[resi] + group

                group = []
                line = line.split()
                residues[line[1]] = []
                resi = line[1]

            if line[0:5] == "GROUP":
                if group != []:
                    flag = 0
                    for atom in group:
                        if atom[0] not in ["C","H"]:
                            flag = 1

                    if flag == 0:
                        residues[resi] = residues[resi] + group

                group = []

            if line[0:5] == "ATOM ":
                line = line.split()
                group.append(line[1])

    return residues


def write_cand(residues):
    """
    write the potential atom candidates fro hydrophobic contact calculation in the folowing format:
    RESID [list of potential atom candidates]
    :param residues:
    :return:none
    """

    for key in residues.keys():
        atoms = ",".join(residues[key])
        print("{0}\t{1}".format(key,atoms))

if __name__ == '__main__':
    res = ReadTopparCharmm(sys.argv[1])
    write_cand(res)