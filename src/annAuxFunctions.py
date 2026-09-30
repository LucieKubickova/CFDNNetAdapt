
import os
import math
import numpy as np
import pyrennModV3 as prn
from postOptimizerV1 import *
from postProcesserV26 import *
from configureAndRunV15 import *

def draw_neural_net(ax, left, right, bottom, top, layer_sizes, layerNames=None):
    '''
    Draw a neural network cartoon using matplotilb.
    
    :usage:
        >>> fig = plt.figure(figsize=(12, 12))
        >>> draw_neural_net(fig.gca(), .1, .9, .1, .9, [4, 7, 2])
    
    :parameters:
        - ax : matplotlib.axes.AxesSubplot
            The axes on which to plot the cartoon (get e.g. by plt.gca())
        - left : float
            The center of the leftmost node(s) will be placed here
        - right : float
            The center of the rightmost node(s) will be placed here
        - bottom : float
            The center of the bottommost node(s) will be placed here
        - top : float
            The center of the topmost node(s) will be placed here
        - layer_sizes : list of int
            List of layer sizes, including input and output dimensionality
    '''
    n_layers = len(layer_sizes)
    v_spacing = (top - bottom)/float(max(layer_sizes))
    h_spacing = (right - left)/float(len(layer_sizes) - 1)
    # Nodes
    for n, layer_size in enumerate(layer_sizes):
        layer_top = v_spacing*(layer_size - 1)/2. + (top + bottom)/2.
        for m in range(layer_size):
            circle = plt.Circle((n*h_spacing + left, layer_top - m*v_spacing), v_spacing/4.,
                                color='w', ec='k', zorder=4)
            ax.add_artist(circle)
    # Edges
    for n, (layer_size_a, layer_size_b) in enumerate(zip(layer_sizes[:-1], layer_sizes[1:])):
        layer_top_a = v_spacing*(layer_size_a - 1)/2. + (top + bottom)/2.
        layer_top_b = v_spacing*(layer_size_b - 1)/2. + (top + bottom)/2.
        for m in range(layer_size_a):
            for o in range(layer_size_b):
                line = plt.Line2D([n*h_spacing + left, (n + 1)*h_spacing + left],
                                  [layer_top_a - m*v_spacing, layer_top_b - o*v_spacing], c='k')
                ax.add_artist(line)

    # Layer names
    if layerNames is not None:
        for n,layer_size in enumerate(layer_sizes):
            layer_top = v_spacing*(layer_size - 1)/2. + (top + bottom)/2.
            layer_bot = layer_top - (layer_size-0.5)*v_spacing

            if n > 0 and n < len(layer_sizes)-1:
                ax.text(
                    n*h_spacing+left,layer_bot,
                    layerNames[n-1],
                    horizontalalignment='center',
                    verticalalignment='center',
                )

def loadAndScaleData(dataDir,dataNm,nPars,nObjs):
    """ function to load CFD samples and scale then in <0,1> """

    # load samples
    with open(dataDir + dataNm,'r') as file:
        data = file.readlines()

    # remove annotation row
    data = data[1::]

    # convert the data to numpy array
    dataNum = []
    for line in data:
        lineSpl = line.split(',')
        row = []
        for value in lineSpl:
            row.append(float(value))
        dataNum.append(row)

    dataNum = np.array(dataNum)

    # scale the data
    colMins = np.min(dataNum,axis=0)
    colMaxs = np.max(dataNum,axis=0)
    for rowInd in range(dataNum.shape[0]):
        for colInd in range(dataNum.shape[1]):
            dataNum[rowInd,colInd] = (dataNum[rowInd,colInd]-colMins[colInd])/(colMaxs[colInd]-colMins[colInd])

    source = dataNum[:,:nPars].T
    target = dataNum[:,nPars:nPars+nObjs].T

    return source,target


def loadAndScaleDataExtScales(dataDir,dataNm,colMins,colMaxs,nPars,nObjs):
    """ function to load CFD samples and scale them based on supplied mins and maxs"""

    # load samples
    with open(dataDir + dataNm,'r') as file:
        data = file.readlines()

    # remove annotation row
    data = data[1::]

    # convert the data to numpy array
    dataNum = []
    for line in data:
        lineSpl = line.split(',')
        row = []
        for value in lineSpl:
            row.append(float(value))
        dataNum.append(row)

    dataNum = np.array(dataNum)

    # scale the data
    for rowInd in range(dataNum.shape[0]):
        for colInd in range(dataNum.shape[1]):
            dataNum[rowInd,colInd] = (dataNum[rowInd,colInd]-colMins[colInd])/(colMaxs[colInd]-colMins[colInd])

    source = dataNum[:,:nPars].T
    target = dataNum[:,-nObjs:].T

    return source,target

def sortData(source,target):
    """ function to sort the data according to source for postprocessing """

    sourceSorted = source
    targetSorted = []
    sourceSorted = sourceSorted.flatten().tolist()
    for ind in range(target.shape[0]):
        auxCol = list(target[ind,:])
        auxColSorted = [val for _,val in sorted(zip(sourceSorted,auxCol))]
        targetSorted.append(auxColSorted)

    targetSorted =  np.array(targetSorted)
    dummy = sourceSorted.sort()

    sourceSorted = np.array(sourceSorted)

    return sourceSorted,targetSorted

def getScalesFromFile(dataDir,dataNm):
    """ function to get scales from the given file """

    # load samples
    with open(dataDir + dataNm,'r') as file:
        data = file.readlines()

    # remove annotation row
    data = data[1::]

    # convert the data to numpy array
    dataNum = []
    for line in data:
        lineSpl = line.split(',')
        row = []
        for value in lineSpl:
            row.append(float(value))
        dataNum.append(row)

    dataNum = np.array(dataNum)

    # scale the data
    colMins = np.min(dataNum,axis=0)
    colMaxs = np.max(dataNum,axis=0)

    return colMins,colMaxs

def annOptim(vars,nets,nObjs,constr,cfdMaxs,cfdMins,dummy = 1e6):
    """ function to return the costs for optimization """

    # rescale the pars
    netPars = list()
    for p in range(len(vars)):
        netPars.append(vars[p]*(cfdMaxs[p] - cfdMins[p]) + cfdMins[p])

    convCPs = [[netPars[0], netPars[1]], [netPars[2], netPars[3]]]
    diffuserCPs = [[netPars[4], netPars[5]],[netPars[6], netPars[7]],[netPars[8],netPars[9]],[netPars[10],netPars[11]]]
    LConv = netPars[-2]
    LDiff = netPars[-1]

    isGood = True
    for j in range(len(convCPs)-1):
        if convCPs[j][0] >= convCPs[j+1][0]:
            isGood = False

    for j in range(len(diffuserCPs)-1):
        if diffuserCPs[j][0] >= diffuserCPs[j+1][0]:
            isGood = False

    if LConv+LDiff > constr[0]:
        isGood = False

    if isGood:
        netIn = np.array(vars)
        netIn = np.expand_dims(netIn,axis = 1)

        costOut = list()

        for i in range(len(nets)):
            costOut.append(prn.NNOut(netIn,nets[i]).squeeze())

        costOut = np.array(costOut)
        costOut = costOut.mean(axis = 0)

        return costOut

    else:
        return np.array([dummy]*nObjs)

def annOptimLen(vars,nets,nObjs,constr,cfdMaxs,cfdMins,dummy = 1e6):
    """ function to return the costs for optimization """

    # rescale the pars
    netPars = list()
    for p in range(len(vars)):
        netPars.append(vars[p]*(cfdMaxs[p] - cfdMins[p]) + cfdMins[p])

    convCPs = [[netPars[0], netPars[1]], [netPars[2], netPars[3]]]
    diffuserCPs = [[netPars[4], netPars[5]],[netPars[6], netPars[7]],[netPars[8],netPars[9]],[netPars[10],netPars[11]]]
    LConv = netPars[-2]
    LDiff = netPars[-1]

    isGood = True
    for j in range(len(convCPs)-1):
        if convCPs[j][0] >= convCPs[j+1][0]:
            isGood = False

    for j in range(len(diffuserCPs)-1):
        if diffuserCPs[j][0] >= diffuserCPs[j+1][0]:
            isGood = False

    if LConv+LDiff > constr[0]:
        isGood = False

    if isGood:
        netIn = np.array(vars)
        netIn = np.expand_dims(netIn,axis = 1)

        costOut = list()

        for i in range(len(nets)):
            costOut.append(prn.NNOut(netIn,nets[i]).squeeze())

        costOut = np.array(costOut)
        costOut = costOut.mean(axis = 0)

        return [costOut,vars[-2]+vars[-1]]

    else:
        return np.array([dummy]*nObjs)

def annOptimLenMinMax(vars,nets,nObjs,constr,cfdMaxs,cfdMins,dummy = 1e6):
    """ function to return the costs for optimization """

    # rescale the pars
    netPars = list()
    for p in range(len(vars)):
        netPars.append(vars[p]*(cfdMaxs[p] - cfdMins[p]) + cfdMins[p])

    convCPs = [[netPars[0], netPars[1]], [netPars[2], netPars[3]]]
    diffuserCPs = [[netPars[4], netPars[5]],[netPars[6], netPars[7]],[netPars[8],netPars[9]],[netPars[10],netPars[11]]]
    LConv = netPars[-2]
    LDiff = netPars[-1]

    isGood = True
    for j in range(len(convCPs)-1):
        if convCPs[j][0] >= convCPs[j+1][0]:
            isGood = False

    for j in range(len(diffuserCPs)-1):
        if diffuserCPs[j][0] >= diffuserCPs[j+1][0]:
            isGood = False

    if LConv+LDiff > constr[0]:
        isGood = False

    # check angles
    xC1 = netPars[0]
    yC1 = netPars[1]
    xC2 = netPars[2]
    yC2 = netPars[3]
    xD1 = netPars[4]
    yD1 = netPars[5]
    xD2 = netPars[6]
    yD2 = netPars[7]
    xD3 = netPars[8]
    yD3 = netPars[9]
    xD4 = netPars[10]
    yD4 = netPars[11]
    LConv = netPars[12]
    LDiff = netPars[13]

    WConv = 0.0175
    WMxT = 0.014*0.5
    WDiff = 0.0178

    # converging part
    dXC1 = (xC1 - 0)*LConv
    dXC2 = (xC2 - xC1)*LConv
    dXC3 = (1 - xC2)*LConv

    dYC1 = (1 - yC1)*(WConv - WMxT)
    dYC2 = (yC1 - yC2)*(WConv - WMxT)
    dYC3 = (yC2 - 0)*(WConv - WMxT)

    angC1 = math.atan(dYC1/dXC1)
    if angC1 < constr[1] or angC1 > constr[2]:
        isGood = False
    angC2 = math.atan(dYC2/dXC2)
    if angC2 < constr[3] or angC2 > constr[4]:
        isGood = False
    angC3 = math.atan(dYC3/dXC3)
    if angC3 < constr[5] or angC3 > constr[6]:
        isGood = False

    # diffuser
    dXD1 = (xD1 - 0)*LDiff
    dXD2 = (xD2 - xD1)*LDiff
    dXD3 = (xD3 - xD2)*LDiff
    dXD4 = (xD4 - xD3)*LDiff
    dXD5 = (1 - xD4)*LDiff

    dYD1 = (yD1 - 0)*(WDiff - WMxT)
    dYD2 = (yD2 - yD1)*(WDiff - WMxT)
    dYD3 = (yD3 - yD2)*(WDiff - WMxT)
    dYD4 = (yD4 - yD3)*(WDiff - WMxT)
    dYD5 = (1 - yD4)*(WDiff - WMxT)

    angD1 = math.atan(dYD1/dXD1)
    if angD1 < constr[7] or angD1 > constr[8]:
        isGood = False
    angD2 = math.atan(dYD2/dXD2)
    if angD2 < constr[9] or angD2 > constr[10]:
        isGood = False
    angD3 = math.atan(dYD3/dXD3)
    if angD3 < constr[11] or angD3 > constr[12]:
        isGood = False
    angD4 = math.atan(dYD4/dXD4)
    if angD4 < constr[13] or angD4 > constr[14]:
        isGood = False

    if isGood:
        netIn = np.array(vars)
        netIn = np.expand_dims(netIn,axis = 1)

        costOut = list()

        for i in range(len(nets)):
            costOut.append(prn.NNOut(netIn,nets[i]).squeeze())

        costOut = np.array(costOut)
        costOut = costOut.mean(axis = 0)

        return [costOut,vars[-2]+vars[-1]]

    else:
        return np.array([dummy]*nObjs)

def annOptimLenCheckMesh(vars,nets,nObjs,constr,cfdMaxs,cfdMins,consPars,Allrun,dummy = 1e6):
    """ function to return the costs for optimization """

    # rescale the pars
    netPars = list()
    for p in range(len(vars)):
        netPars.append(vars[p]*(cfdMaxs[p] - cfdMins[p]) + cfdMins[p])

    convCPs = [[netPars[0], netPars[1]], [netPars[2], netPars[3]]]
    diffuserCPs = [[netPars[4], netPars[5]],[netPars[6], netPars[7]],[netPars[8],netPars[9]],[netPars[10],netPars[11]]]
    LConv = netPars[-2]
    LDiff = netPars[-1]

    isGood = True
    for j in range(len(convCPs)-1):
        if convCPs[j][0] >= convCPs[j+1][0]:
            isGood = False

    for j in range(len(diffuserCPs)-1):
        if diffuserCPs[j][0] >= diffuserCPs[j+1][0]:
            isGood = False

    if LConv+LDiff > constr[0]:
        isGood = False

    # UGLYYYYYYYYYYYYYYYYY - for consPars
    QInLst = [0.3e-3, 0.4e-3, 0.5e-3]
    pSucLst = [81.325, 81.325, 81.325]
    DNz2 = 0.0043
    endTime = 5000
    edgeFunction = "polyLine"
    cwd = os.getcwd().split("caseConstructorEjectorHorizontal")[0] + "/"
    genDir = cwd + "caseConstructorEjectorHorizontal/"
    baseCase = genDir + "10_baseCaseFixedP/"
    baseDir = cwd + "Simulations/ANNCheck/"
    LMxT = 0.2
    caseConstructor = "caseConstructorFixedPHorizontalV2.py"
    caseID = os.getpid()
    
    # checkMesh
    if isGood:
        caseConsArgs = []
        pars = consPars[:]
        parValues = []
        for par in pars:
            parValues.append(eval(par))

        changeDict = {
                "changeParsInCaseConstructor":[pars, parValues]
                }

        simulations = configureAndRun(genDir, caseConstructor, baseCase, caseID = caseID)
        simulations.makeCase(changeDict) # create simulation

        caseDirList = glob.glob(baseDir+"*"+str(caseID))

        for caseDir in caseDirList:
            simulations.runSingleCase(caseDir + "/", Allrun) # run simulation

        # evaluate the cases
        for caseDir in caseDirList:
            case = postProcesser(caseDir + "/")
            if not case.meshCheck(checkMesh = False): # checkMesh checked during simulation evaluation
                isGood = False
                break

            if not case.isComputed(): # checkMesh check will show itself here
                isGood = False
                break

    if isGood:
        netIn = np.array(vars)
        netIn = np.expand_dims(netIn,axis = 1)

        costOut = list()

        for i in range(len(nets)):
            costOut.append(prn.NNOut(netIn,nets[i]).squeeze())

        costOut = np.array(costOut)
        costOut = costOut.mean(axis = 0)

        return [costOut,vars[-2]+vars[-1]]

    else:
        return np.array([dummy]*nObjs)

def annOptimWScales(vars,net,scales):
    """ scaled variant of annOptim, not used at the moment """

    varMins,varMaxs,critMins,critMaxs = scales

    # Note: net can work with inputs in <0,1>
    # Note: net gives outputs in <0,1>
    sVars = deepcopy(vars)
    # -- prepare variables to be taken in by the net
    for varInd in range(len(vars)):
        sVars[varInd] = (vars[varInd] - varMins[varInd])/(varMaxs[varInd]-varMins[varInd])

    # -- case specific
    W = sVars
    netIn = np.array([W])

    # -- evaluate the net
    costOut = prn.NNOut(netIn,net).squeeze()

    # -- scale the net outs back
    for critInd in range(len(costOut)):
        costOut[critInd] = costOut[critInd]*(critMaxs[critInd] - critMins[critInd]) + critMins[critInd]

    return costOut

def prepOutDir(outDir, dirLstMk = []):
    """ function to prepare the output directory """

    if not os.path.exists(outDir):                                      #check if outDir exists
        os.makedirs(outDir)

    for dr in dirLstMk:
        if not os.path.exists(outDir + dr):                             #if it does not already exist
            os.makedirs(outDir + dr)
