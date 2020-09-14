
![](images/logo_periprot.png "logo" )
# PeriProt

Periprot1.0 is a sofware built to analysis Peripheral Membrane Proteins (PMPs) structure or Molecular Dynamics (MD) simulation data. Periprot can only analyse date from NAMD or CHARMM simulation (PSF format for the topology and dcd format for the trajectory). 

The available analysis are:

 - Hydrogen bond network analysis (MD)
 - Hydrophobic contact analysis (MD)
 - Cation-Pi interaction analysis (MD)
 - Depth of Anchoring (MD)
 - Protein - Membrane Distances (MD)
 - Electron Density Profil (MD)
 - Macrodipole (PDB structure)

# Usage

optional arguments:
  -h, --help        show this help message and exit
  -top TOP          psf file only
  -traj TRAJ        trajectories (DCD format)
  -hydro            Hydrophobic contact analysis
  -hbond            hbond analysis
  -hbcand HBCAND    hbond candidates
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
 
 # Hydrophobic contact analysis
 
In Periprot, an hydrophobic contact are consider between to candidates atom (one from the protein and one from the bilayer) if there distance are less than 3 Å. The candidates atom can be picked automatiocally or the user can give his own candidates list.

usage: python PeriProt.py -top mydata.psf -traj mydata.dcd -pdb mydata.pdb -hydro -out example

The outputs will two file:
 - a csv file ("example_hydrophobic_contact_raw_data.csv") contating the raw data (all hydrophobic contact in all frame)
 - a csv file ("example_average_hydrophobic_contact_per_frame.csv") containg the average number of hydrophobic contact per residues
 - a pdb file ("example__average_hydrophobic_contact_per_frame.pdb"). If a PDB file is given on the command line, a new PDB file will be created with the average number of hydrophobic contact as the B-factor. 
 
Example of output:
 
 
 
