#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import  MDAnalysis as mda

"""
PeriProt 1.0
Hbond Anlaysis
"""

__author__ = "Emmanuel Edouard MOUTOUSSAMY"
__version__  = "1.0.0"
__date__ = "2020/08"
__copyright__ = "CC_by_SA"
__dependencies__ = "Numpy,MDAnalysis and argparse"

def hbond_list(rawdata):
    """
    This function allow to get additinal ...
    :param rawdata:
    :return:
    """
    donor = []
    acceptor = []
    groups = {}

    with open(rawdata) as inputfile:
        for line in inputfile:
            splitted_line = line.split()

            if "!" in line:
                groups[splitted_line[0]] = splitted_line[1:]

            elif "#donor" in line:
                donor = splitted_line[1:]

            elif "#acceptor" in line:
                acceptor = splitted_line[1:]

    return acceptor,donor,groups

