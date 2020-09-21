#!/usr/bin/env python3

import common
import MDAnalysis as mda

"""
SET OF FUNCTION TO CALCULATE THE CATION-PI INTERACTION
"""

TYR = ["CG", "CZ", "CD1", "CD2", "CE2", "CE1"]
TRP = ["CE3", "CD2", "CE2", "CZ2", "CH2", "CZ3"]


def IsTheAromClose(close_amino_acids,aromatics):
    """
    Check if the aromatic residues are close to the bilayer
    :param close_amino_acids: list of residues close to the bilayer
    :param aromatics: list aromatics
    :return: a list with the aromatics close to the bilayer
    """

    aromatics_close  = []

    for arom in aromatics:
        if int(arom) in close_amino_acids:
            aromatics_close.append(arom)

    return aromatics_close


def ParseAromDist(dist_mat,arom_list,frame,output_raw,lipids_list,counting):
    """
    Parse the distance matrix calculated by Mdanalysis to find Cation-pi int.
    it also write the raw data
    :param dist_mat: distance matrix ( N atom from bilayer vs close aromatics)
    :param arom_list: list of aromatics
    :param frame: current frame
    :param output_raw: outputfile
    :param lipids_list: list of lipids
    :param counting: dico to count interaction per residues
    :return:
    """

    for i in range(dist_mat.shape[0]):
        index_low = 0
        index_up = 6
        for arom in arom_list:
            dist_mat_arom = dist_mat[i, index_low:index_up]

            index_low = index_up
            index_up = index_up + 6

            if max(dist_mat_arom) <= 7:
                if max(dist_mat_arom) - min(dist_mat_arom) < 1.5:
                    if arom not in counting.keys():
                        counting[arom] = 1
                    else:
                        counting[arom] +=1

                    output_raw.write("{0},{1},{2}\n".format(frame,arom,lipids_list[i]))


def RunAromatic(psf,dcd,psf_info,segmemb,close_lipids,close_amino_acids,outname,pdb,segprot,first_frame,\
                last_frame,skip):
    """
    Run Analysis of Cation-Pi interaction
    :param psf: PSF (topology) file
    :param dcd: trajectory file
    :param psf_info: structure conataing infos from PSF
    :param segmemb: membrane segid
    :param close_lipids: list of lipids close to the protein
    :param close_amino_acids: list of amino_acids close to the bilayer
    :param outname: name of the output
    :param pdb: pdb file
    :param segprot: segid of the protein
    :return:
    """
    univers = mda.Universe(psf,dcd)

    close_lipids = common.convertIntToStr(close_lipids)

    output_raw = open("%s_cation_pi_int_raw_data.csv"%outname,"w")
    output_raw.write("#frame,resProt,resMemb\n")
    output = open("%s_cation_pi_int.csv"%outname,"w")
    output.write("#Residue,occupancy\n")

    tyrosine_list = IsTheAromClose(close_amino_acids,psf_info.tyrosines)
    tryptophane_list = IsTheAromClose(close_amino_acids,psf_info.tryptophanes)

    counting = {}

    nb_frame = 0

    frame_to_read = first_frame

    for ts in univers.trajectory:
            nb_frame += 1
            all_nitrogen = univers.select_atoms("(name N) and (segid %s) and (resid %s)" %(segmemb," ".join(close_lipids)))
            tyrs = univers.select_atoms("(resid %s) and (name %s)"%(" ".join(tyrosine_list)," ".join(TYR)))
            trps = univers.select_atoms("(resid %s) and (name %s)"%(" ".join(tryptophane_list)," ".join(TRP)))

            dist_tyrosines = mda.analysis.distances.distance_array(all_nitrogen.positions, tyrs.positions,\
                                                                   box=None, result=None,backend='serial')


            dist_trp = mda.analysis.distances.distance_array(all_nitrogen.positions, trps.positions, box=None,\
                                                             result=None,backend='serial')

            ParseAromDist(dist_tyrosines,tyrosine_list,ts.frame,output_raw,close_lipids,counting)
            ParseAromDist(dist_trp,tryptophane_list,ts.frame,output_raw,close_lipids,counting)

            frame_to_read += skip

    occupancies = {}

    for amino_acid in counting.keys():
        output.write("{0},{1}\n".format(amino_acid,(float(counting[amino_acid])/float(nb_frame))*100))

    for residue in psf_info.protresid:
        if str(residue) not in occupancies.keys():
            occupancies[str(residue)] = 0

    if pdb != 'None':
        common.Mapped(pdb,occupancies,segprot,outname,"cation_pi_occupancies")