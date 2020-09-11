#!/usr/bin/env python3

'''
SET OF COMMON FUNCTION for different analysis
'''


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

        self.protseq = [] # sequence
        self.membresid = []  # sequence
        self.protresid = []
        self.protindex = [0,0]
        self.membindex = [0,0]
        self.protCharge = {}
        self.mass = {}

        self.fillTopoStruct(PSF, segidMEMB, segidPROT)

    def FillProtInfos(self,psfline):
        """
        Fill information regarding the protein
        :param psfline: a line of the PSF file
        :return:
        """

        if int(psfline[0]) < self.protindex[0]: #find the first atom number of the protein
            self.protindex[0] = int(psfline[0])

        if int(psfline[0])  > self.protindex[1]: #find the last atom number of the protein
            self.protindex[1] = int(psfline[0])

        self.protCharge[psfline[0]] = float(psfline[6]) #Collect the charge
        self.mass[psfline[0]] = float(psfline[7]) #Collect the mass

        if int(psfline[2]) not in self.protresid: #Collect the resid
            self.protresid.append(int(psfline[2]))


    def FillMembInfos(self,psfline):
        """
        Fill information regarding the membrane
        :param psfline: a line of the PSF file
        :return:
        """
        if int(psfline[0]) < self.membindex[0]: #find the first atom number of the membrane
            self.membindex[0] = int(psfline[0])

        if int(psfline[0])  > self.membindex[1]: #find the last atom number of the membrane
            self.membindex[1] = int(psfline[0])

        self.protCharge[psfline[0]] = float(psfline[6]) #Collect the charge
        self.mass[psfline[0]] = float(psfline[7]) #Collect the mass


        if int(psfline[2]) not in self.membresid: #Collect the resid
            self.membresid.append(int(psfline[2]))


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

def Mapped(pdb,dico,segid_prot,outname):
    """
    Map the dist. on a pdb file
    :param pdb: a PDB file of a protein/membrane complex
    :param dico: dico containing data to change
    :param segid_prot: segid pf the protein
    :param outname: name of output
    :return: write a pdb file with the new dada as the B-factor
    """

    pdbout = "%s_depth_of_anchoring.pdb"%outname

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