
# import
import os
import sys
sys.path.insert(1, "../../src")
sys.path.insert(1, "../../thirdParty")
import csv
import numpy as np
import dill as pickle
from matplotlib import cm
from testFunctions import *
from CFDNNetAdaptV3 import *
import matplotlib.pyplot as plt

# parameters
runDir = "01_algoRuns/run_01/"
xName = "f1"
yName = "f2"
logName = "log.out"
parName = "optimOut.plat"
moeaName = "31_igdValues.dat"

# prepare CFDNNetAdapt
algorithm = CFDNNetAdapt()

# problem specification
algorithm.nPars = 2
algorithm.nObjs = 2
algorithm.nOuts = 2
algorithm.mainDir = "01_algoRuns/"
algorithm.smpDir = "00_prepData/"
algorithm.prbDir = ""
algorithm.dataNm = "10_platypusAllSolutions.dat"
algorithm.minMax = ""

# prepare plot
fig = plt.figure(figsize = (16,9))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

# read scales
smpMins, smpMaxs = algorithm.getScalesFromFile(algorithm.smpDir + algorithm.prbDir, algorithm.dataNm)

# read nsgaii data
fileName = algorithm.smpDir + algorithm.prbDir + moeaName
with open(fileName, 'r') as file:
    reader = csv.reader(file)

    cols = next(reader)

    data = list()
    for line in reader:
        data.append([float(i) for i in line])

xI = cols.index("step")
yI = cols.index("igd")

data = np.array(data)
ax1.scatter(data[:,xI], data[:,yI], label = "NSGA-II", color = "black", marker = "x")

# prepare and plot optimal solution
optSols = optSolsZDT6(100, algorithm.nPars)
f1s = list()
f2s = list()
for i in range(len(optSols)):
    f1, f2 = ZDT6(optSols[i])
    f1s.append(f1)
    f2s.append(f2)

# read cfdnnetadapt log
fileName = runDir + logName
with open(fileName, 'r') as file:
    data = file.readlines()

# get the best DNNs from each iteration
sizes = list()
bestDNNs = list()
for line in data:
    if "Best DNN found " in line:
        bestDNNs.append(line.split()[-1])

    elif "Using" in line:
        sams = line.split()
        num = int(sams[1]) + int(sams[4]) + int(sams[8])
        sizes.append(num)

# go over steps and plot
igds = list()
for n in range(len(bestDNNs)):
    stepDir = "step_%04d/" %(n+1)
    fileName = runDir + stepDir + bestDNNs[n] + "/" + parName

    # read data from optimization
    with open(fileName, 'rb') as file:
        [population,result,name,problem] = pickle.load(file, encoding="latin1")

    # process data
    data = list()
    for i in range(len(result)):
        netPars = result[i].variables[:]
        netOuts = result[i].objectives[:]

        # concatenate and descale
        point = netPars + netOuts
        point = np.array(point)
        point = point*(smpMaxs - smpMins) + smpMins

        # add
        data.append(point)
    data = np.array(data)

    # compute inverted generational distance
    igd = 0.0
    dists = list()
    for i in range(len(optSols)):
        dist = 100.0
        for j in range(len(data)):
            nondom = np.array(optSols[i][:algorithm.nPars])
            netple = data[j][:algorithm.nPars]

            aux = np.linalg.norm(nondom - netple)
            if aux < dist:
                dist = aux

        dists.append(dist)
        igd += dist

    igd /= len(optSols)
    igds.append(igd)

ax1.scatter(sizes, igds, label = "CFDNNetAdapt")
ax1.legend()

ax1.set_yscale("symlog", base = 10, linthresh = 1e-3)
ax1.set_xlabel("number of samples")
ax1.set_ylabel("mean inverted generational distance (IGD)")

plt.savefig("./igdValues.png")
plt.close()
