
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
 
 
