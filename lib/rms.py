#!/usr/bin/env python3


import sys
import MDAnalysis
import MDAnalysis.analysis.rms
import matplotlib.pyplot as plt



PSF = sys.argv[1]
DCD = sys.argv[2]

u = MDAnalysis.Universe(PSF,DCD)
ref = MDAnalysis.Universe(PSF,DCD)     # reference closed AdK (1AKE) (with the default ref_frame=0)
#ref = MDAnalysis.Universe(PSF,CRD)    # reference open AdK (4AKE)


R = MDAnalysis.analysis.rms.RMSD(u, ref,
           select="backbone",             # superimpose on whole backbone of the whole protein
           groupselections=["backbone"])                                    # NMP
R.run()

rmsd = R.rmsd.T   # transpose makes it easier for plotting
time = rmsd[1]
fig = plt.figure(figsize=(4,4))
ax = fig.add_subplot(111)
ax.plot(time, rmsd[2], 'k-',  label="all")
ax.legend(loc="best")
ax.set_xlabel("time (ps)")
ax.set_ylabel(r"RMSD ($\AA$)")
fig.savefig("rmsd_all_CORE_LID_NMP_ref1AKE.pdf")


frame = 1
for i in range(len(time)):
    print("{0},{1}".format(frame,rmsd[2][i]))
    frame += 1