
![](images/logo_periprot.png "logo" )
## PeriProt

Periprot1.0 is a sofware built to analysis Peripheral Membrane Proteins (PMPs) structure or Molecular Dynamics (MD) simulation data. Periprot can only analyse date from NAMD or CHARMM simulation (PSF format for the topology and dcd format for the trajectory). 

The available analysis are:

 - (1) Depth of Anchoring (MD)
 - (2) Hydrophobic contact analysis (MD)
 - (3) Cation-Pi interaction analysis (MD)
 - (4) Hydrogen bond network analysis (MD)
 - (5) Electron Density Profil (MD)
 - (6) Protein - Membrane Distances (MD)
 - (7) Macrodipole (PDB structure)

## Dependencies

## Usage
# ARGUMENTS 
```text
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
```

 ## (1) Depth of anchoring
 

Periprot can be use to calculate the average depth of anchoring (DOA) of each residue in the bilayer. The DOA here represent the distance from the carbon alpha (CA) of the residue and the average phosphate plane. I the CA is below the phosphate plane the DOA will be negative.

usage: python PeriProt.py -top mydata.psf -traj mydata.dcd -pdb mydata.pdb -depth -out example

The output files:
 - a csv file ("example_depth_of_anchoring.csv") containg the average DOA and the associated standard deviation for each residue
 - A figure ("example_depth_of_anchoring.png") with the average DOA vs the residue number.
 - A PDB file ("example_depth_of_anchoring.pdb"). If a PDB file is given on the command line, a new PDB file will be created with the average average DOA during the MD as the B-factor. 

Example of output:
 ![](images/doa_out_example.png "DOA" )

 ## (2) Hydrophobic contact analysis
 
In Periprot, an hydrophobic contact are consider between to candidates atom (one from the protein and one from the bilayer) if there distance are less than 3 Å. The candidates atom can be picked automatiocally or the user can give his own candidates list.

usage: python PeriProt.py -top mydata.psf -traj mydata.dcd -pdb mydata.pdb -hydro -out example

The output files:
 - a csv file ("example_hydrophobic_contact_raw_data.csv") contating the raw data (all hydrophobic contact in all frame)
 - a csv file ("example_average_hydrophobic_contact_per_frame.csv") containg the average number of hydrophobic contact per residues
 - a PDB file ("example__average_hydrophobic_contact_per_frame.pdb"). If a PDB file is given on the command line, a new PDB file will be created with the average number of hydrophobic contact as the B-factor. 
 
Example of output:
![](images/hydrophobic_out_example.png "hydrophobes" )

## (3) Cation-Pi interaction analysis
The cation-pi interaction between tyrosines or tryptophane and the bilayer can be evaluated using Periprot. The cation-pi interaction is consider when the distances between each carbon of the aromatic cycle and the nitrogen atom of a PC lipid is less that 7 Å. Moreover, these distance should not differ by more than 1.5 Å.

/!\ Warning: the calculation will be done between PC lipid and the protein, other lipids (PS or PE) will be add later /!\

usage: python PeriProt.py -top mydata.psf -traj mydata.dcd -pdb mydata.pdb -catpi -out example

The output files:
 - a csv file ("example_cation_pi_int.csv") contating the occupencies of the cation-pi interactions.
 - a PDB file ("example_cation_pi_occupencies.pdb"). If a PDB file is given on the command line, a new PDB file will be created with the average number of hydrophobic contact as the B-factor. 
 
 Example of output:
 ![](images/cation_pi_out_example.png "catpi" )


## (4) Hydrogen bond network analysis

## (5) Electron density profil

## (6) Protein - Membrane Distance

The protein-Membrane distance can be calculated with Periprot. The user can be should between three different distances. When '-dist' the following dialogue will appears:
```text
    Which dist:
    1: Prot (COM) - Memb (COM) distance
    2: Prot (COM) - Phosphate Plane (COM) distance
    3: Custom distance between two atoms (one in the prot. and one on the memb.)
    
    Selection:
```

1: distance between the protein center of mass and the membrane center of mass
2: distance between the protein center of mass and the phosphate plane center of mass
3: Custom distance between two atoms (one in the prot. and one on the memb.)

usage: python PeriProt.py -top mydata.psf -traj mydata.dcd -pdb mydata.pdb -dist -out example

The output files:
 - a csv file ("example__prot_memb_distance.csv") contating the prot./memb. distance.
 - a figure ("example__prot_memb_distance.png"). Distance vs frame

Example of output:
 ![](images/distance_out_example.png "dist" )

## (7) Macrodipole

