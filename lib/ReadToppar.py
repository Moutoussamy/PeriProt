#!/usr/bin/env python3

import sys
import numpy as np
import  MDAnalysis as mda

'''
SET OF COMMON FUNCTION for different analysis
'''



def ReadTopparCharmm(toppar):
    resi = {}
    with open(toppar) as inputfile:
        for line in inputfile:
            if line[0:4] == "RESI":
                line = line.split()
                resi[line[1]] = []

            if line[0:5] == "GROUP":
                


if __name__ == '__main__':
    ReadTopparCharmm(sys.argv[1])