
# import
import os
import sys
sys.path.insert(1, "../../../src")
sys.path.insert(1, "../../../thirdParty")
import csv
import math
import numpy as np
import dill as pickle
import multiprocessing
import platypusModV4 as plat
from testFunctions import *
from CFDNNetAdaptV5 import *


# cost function evaluation 
def smpEvaluation(i):
    # evaluate the cases
    checkOut = [0] * (algorithm.nObjs)

    netPars = algorithm.population[i].variables[:]
    netOuts = algorithm.population[i].objectives[:]

    # rescale the pars
    for p in range(len(netPars)):
        netPars[p] = netPars[p]*(algorithm.smpMaxs[p] - algorithm.smpMins[p]) + algorithm.smpMins[p]

    checkOut = F7(netPars)

    # rescale for ANN
    for co in range(len(checkOut)):
        checkOut[co] = (checkOut[co] - algorithm.smpMins[algorithm.nPars+co])/(algorithm.smpMaxs[algorithm.nPars+co] - algorithm.smpMins[algorithm.nPars+co])

    # compute the error
    delta = 0
    delta += abs(netOuts[0] - checkOut[0])
    delta += abs(netOuts[1] - checkOut[1])

    algorithm.outFile.write("Doing CFD check no. " + str(algorithm.toCompare.index(i)) + " with parameters " + str(netPars) + "\n")
    algorithm.outFile.write("no. " + str(algorithm.toCompare.index(i)) + " ANN outs were " + str(netOuts) + "\n")
    algorithm.outFile.write("no. " + str(algorithm.toCompare.index(i)) + " CFD outs were " + str(checkOut) + " delta " + str(delta) + "\n")
    algorithm.outFile.write("CFD check no. " + str(algorithm.toCompare.index(i)) + " done\n")
    algorithm.outFile.flush()

    return delta

def dnnEvaluation(vars):
    """ function to return the costs for optimization """
    dummy = 1e6

    # prepare neural network
    netIn = np.array(vars)
    netIn = np.expand_dims(netIn, axis = 1)

    costOut = list()

    for i in range(len(algorithm.nets)):
        costOut.append(prn.NNOut(netIn, algorithm.nets[i]).squeeze())

    costOut = np.array(costOut)
    costOut = costOut.mean(axis = 0)

    return costOut

# declare CFDNNetAdapt
algorithm = CFDNNetAdapt()

# problem specification
algorithm.nPars = 10
algorithm.nObjs = 2
algorithm.nOuts = 2
algorithm.mainDir = "01_algoRuns/"
algorithm.smpDir = "00_prepData/"
algorithm.prbDir = ""
algorithm.dataNm = "10_platypusAllSolutions.dat"
algorithm.minMax = ""

# algorithm parameters
algorithm.nSam = 10000
algorithm.deltaNSams = [10000]
algorithm.nNN = 1
algorithm.minN = 3
algorithm.maxN = 30
algorithm.nHidLay = 3
algorithm.tol = 1e-6
algorithm.iMax = 3
algorithm.dRN = 1
algorithm.nComps = 16
algorithm.nSeeds = 1

# parameters for ANN training
algorithm.trainPro = 75
algorithm.valPro = 15
algorithm.testPro = 10
algorithm.kMax = 10000
algorithm.rEStop = 1e-5
algorithm.dnnVerbose = False

# parameters for MOP
algorithm.pMins = [0.0] + [0.0]*(algorithm.nPars - 1)
algorithm.pMaxs = [1.0]*algorithm.nPars
algorithm.offSize = 300
algorithm.popSize = 300
algorithm.nGens = 500
algorithm.archive = plat.AdaptiveGridArchive(
        algorithm.popSize*algorithm.nGens*0.01,
        algorithm.nObjs,
        algorithm.nPars**algorithm.nObjs
        )
algorithm.moeaVerbose = False

# evaluation funcs
algorithm.dnnEvalFunc = dnnEvaluation
algorithm.smpEvalFunc = smpEvaluation

# initialize
algorithm.initialize()

# run
algorithm.run()
