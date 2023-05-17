#!/usr/bin/env python3
"""
read and plot o3 spectra from text file
"""

import os
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np

def read_o3(inFile):
    f = open(inFile)
    nhead = 72
    for i in range(nhead):
        hh = f.readline()
    
    l = f.readline()
    l = l.replace('k','').split(',')[1:]
    temp = np.array(l).astype(float)


    wav_o3, xsec_o3 = [],[]
    for l in f:
        a = np.array(l.split(',')).astype(float)
        wav_o3.append(a[0])
        xsec_o3.append(a[1:])

    wav_o3 = np.array(wav_o3)
    xsec_o3 = np.array(xsec_o3)

    return temp,wav_o3, xsec_o3
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    inRoot = "../o3_bremen"
    inFile = '{}/SerdyuchenkoGorshelevVersionJuly2013.dat.csv'.format(inRoot)

    temp,wav_o3,xsec_o3 = read_o3(inFile)

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.semilogy(wav_o3,xsec_o3)
    plt.show()

    i = wav_o3 == 350
    i = np.argmax(i)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(temp,xsec_o3[i,:])
    plt.show()
