#!/bin/python
import sys, os
import re
import matplotlib.pyplot as plt
import numpy as np


def readFile(inputFile):
    maxIdentities = []
    f = open(inputFile,'r')
    for line in f.readlines():
        id1,id2,identity = line.rstrip().split()
        maxIdentities.append(float(identity))
    f.close()
    return maxIdentities;

if __name__ == "__main__":
    if (len(sys.argv) < 3):
        print("USAGE:\n%s <max identities file> <titlte>\n"%(sys.argv[0]));
        exit()

    maxIdentitiesFile = sys.argv[1]
    title = sys.argv[2]
    maxIdentities = readFile(maxIdentitiesFile)

    plt.ylim(ymax=350)
    plt.title(title,fontsize=16)
    plt.ylabel("Count",fontsize=14)
    plt.xlabel("Max Global Alignment Identity",fontsize=14)
    plt.hist(maxIdentities, bins=10,range=(0,100),label=title)
    plt.legend()
    plt.savefig(maxIdentitiesFile + ".svg")
    plt.show()
