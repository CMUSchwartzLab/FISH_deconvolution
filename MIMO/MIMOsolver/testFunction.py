import numpy as np
import random
from scipy.optimize import linear_sum_assignment
import csv, sys
import os, errno

'''
    This function is used to arrange the inferred C matrix to best match the real single cells matrix used in simulation. 
    If CellsInCol is True, then the format is B - C F. Else, the format is B - F C.
'''

def arrangeC (C, Ctrue, threshold=0.001, CellsInCol=True):
    if (CellsInCol):
        cells = C.shape[1]
        accMatrix = np.zeros([cells, cells])
        for i in range(cells):
            for j in range(cells):
                accMatrix[i, j] = 1.0 * np.sum(abs(Ctrue[:, i] - C[:, j]) <= threshold)/C.shape[0]
        row, col = linear_sum_assignment(-accMatrix)
        return [C[:, col], col]
    else:
        cells = C.shape[0]
        accMatrix = np.zeros([cells, cells])
        for i in range(cells):
            for j in range(cells):
                accMatrix[i, j] = 1.0 * np.sum(abs(Ctrue[i, :] - C[j, :]) <= threshold)/C.shape[1]
        row, col = linear_sum_assignment(-accMatrix)
        return [C[col, :], col]

def calcAccuracyAfterArrange(C, Ctrue, threshold= 0.001, CellsInCol = True):
    if (CellsInCol):
        cells = C.shape[1]
        accuracy = np.zeros([cells, 1])
        for i in range(cells):
            accuracy[i, 0] = 1.0 * np.sum(abs(Ctrue[:, i] - C[:, i]) <= threshold)/C.shape[0]
    else:
        cells = C.shape[0]
        accuracy = np.zeros([cells, 1])
        for i in range(cells):
            accuracy[i, 0] = 1.0 * np.sum(abs(Ctrue[i, :] - C[i, :]) <= threshold)/C.shape[1]
    return accuracy

def calcAccuracy (C, Ctrue, threshold= 0.001, CellsInCol = True):
    if (CellsInCol):
        cells = C.shape[1]
        accuracy = np.zeros([cells, 1])
        accMatrix = np.zeros([cells, cells])
        for i in range(cells):
            for j in range(cells):
                accMatrix[i, j] = 1.0 * np.sum(abs(Ctrue[:, i] - C[:, j]) <= threshold)/C.shape[0]
        row, col = linear_sum_assignment(-accMatrix)
        accuracy[:, 0] = accMatrix[row, col]
    else:
        cells = C.shape[0]
        accuracy = np.zeros([cells, 1])
        accMatrix = np.zeros([cells, cells])
        for i in range(cells):
            for j in range(cells):
                accMatrix[i, j] = 1.0 * np.sum(abs(Ctrue[i, :] - C[j, :]) <= threshold)/C.shape[1]
        row, col = linear_sum_assignment(-accMatrix)
        accuracy[:, 0] = accMatrix[row, col]
    
    return accuracy

def calcAccuracyByRow(C, Ctrue, threshold=0.001, CellsInCol=True):
    right_row = 0
    if (CellsInCol):
        DD = np.abs(C - Ctrue)
        for i in range(C.shape[0]):
            if np.sum(DD[i, :]) <= threshold:
                right_row += 1
        right_row = right_row/float(C.shape[0])

    else:
        C = C.T
        Ctrue = Ctrue.T
        DD = np.abs(C - Ctrue)
        for i in range(C.shape[0]):
            if np.sum(DD[i, :]) <= threshold:
                right_row += 1
        right_row = right_row / float(C.shape[0])
        
    return right_row


'''
    This function is used to generate the Fraction Matrix with Dirichlet distribution.
    alpha is a vector with length same as the number of single cells used for simulation.
    Parameters: 
        cellsList: a list of the number of major single cells components in each of three
        regions.
        tumor_samples: a list of tumor samples in each of three regions.
        CellsNoiseList: a list of the number of minor single cells components in each of 
        three regions.
'''
def generateF(cellsList, cellsNoiseList, tumor_samples, alpha=[100, 1, 0.1], CellsInCol=True):
    '''Simulate a cell fraction matrix, choosing the mixture fractions
    from a Dirichlet distribution.  Return the fraction matrix, with
    each column representing a component, and the rows the fractional
    contribution of each single cell.  Also return the Dirichlet
    parameters actually used.

    alpha: Dirichlet parameters for simulation.  There are three
    components: 0) the parameter for the dominant clone; 1) the
    parameter for noise cells in the same region; 2) the parameter for
    noise cells in another region.

    cellsList: an array of length 3 with cellsList[i] equal to the
    number of dominant cells to choose from region i

    cellsNoiseList: an array of length 3 with cellsNoiseList[i] equal
    to the number of dominant cells to choose from region i

    tumorSamples: The number of bulk components to simulate for each
    region, usually 1, 2, or 3.

    CellsInCol: default True.  If false, transpose the simulated
    frequency matrix.

    '''

    # There are 3 sampled regions for each tumor.  We will choose
    # CSel[i] cells to have nonzero contribution to the simulated
    # tumors.
    CSel = np.array(cellsList) + np.array(cellsNoiseList)
    dirAList = []

    # cellList[0] dominant cells from region 0, cellsNoiseList[0]
    # noise cells from region 0, and noise cells from the other two
    # regions.
    dirAList.append([alpha[0]] * cellsList[0] + [alpha[1]] * cellsNoiseList[0] + [alpha[2]] * (CSel[1] + CSel[2]))

    F1 = np.random.dirichlet(dirAList[0], tumor_samples[0])

    # cellList[1] dominant cells from region 1, cellsNoiseList[1]
    # noise cells from region 1, and noise cells from the other two
    # regions.
    dirAList.append([alpha[2]] * CSel[0] + [alpha[0]] * cellsList[1] + [alpha[1]] * cellsNoiseList[1] + [alpha[2]] * CSel[2])

    F2 = np.random.dirichlet(dirAList[1], tumor_samples[1])

    dirAList.append([alpha[2]] * (CSel[0]+CSel[1]) + [alpha[0]] * cellsList[2] + [alpha[1]] * cellsNoiseList[2])

    # cellList[2] dominant cells from region 2, cellsNoiseList[2]
    # noise cells from region 2, and noise cells from the other two
    # regions.
    F3 = np.random.dirichlet(dirAList[2], tumor_samples[2])

    F = np.concatenate((F1, F2, F3), axis=0)
    
    dirA = np.zeros([0, len(dirAList[0])])
    for i in range(len(tumor_samples)):
        if (tumor_samples[i] != 0):
            dirA = np.concatenate((dirA, np.matrix([dirAList[i]] * tumor_samples[i])), axis=0)
    if CellsInCol:
        return [F.transpose(), dirA]
    else:
        return [F, dirA]
'''
    This function is used to read in the single cells data in .csv file.
'''
def readMatrix(filename, dtype=np.int32, capToTen=False):
    arr = []
    with open(filename, 'rb') as f:
        reader = csv.reader(f)
        try:
            for row in reader:
                arr.append(row)
        except csv.Error as e:
            sys.exit('file %s, line %d: %s' % (filename, reader.line_num, e))
    matrixOut = np.matrix(arr, dtype=np.float)
    matrixOut = np.matrix(matrixOut, dtype=dtype)
    if capToTen:
        matrixOut[matrixOut > 10] = 10
    else:
        pass
    return matrixOut


'''
    This function is used to randomly select single cells for simulating the bulk tumor
    B matrix in the simulation. CSel/CReferSel is a vector of the number of single cells
    selected from each region.
'''
def generateCandCRefer (allSC, CSel, CReferSel, tumorType='GBM07', dtype=np.float, CellsInCol=True):
    cells = np.sum(CSel)
    copyNum = allSC.shape[0]
    cellsOrig = np.sum(CReferSel)

    if (tumorType == 'GBM07'or tumorType=='simulated_GBM07'):
        regions = [list(range(0,57)), list(range(57, 132)), list(range(132, allSC.shape[1]))]
    elif (tumorType == 'GBM33' or tumorType =='simulated_GBM33'):
        regions = [list(range(0,59)), list(range(59,135)), list(range(135, allSC.shape[1]))]
    else:
        sys.exit('Please Use GBM07 or GBM33 (or corresponding data from simulated SCS) as tumorType!')


    CIndexList = [0] * cells
    CReferIndexList = [0] * cellsOrig
    
    index0 = 0
    index1 = 0
    for i in range(len(regions)):
        temp = random.sample(regions[i], (CSel[i]+CReferSel[i]))
        CIndexList[index0:(index0+CSel[i])] = temp[0:CSel[i]]
        CReferIndexList[index1:(index1+CReferSel[i])] = temp[CSel[i]:(CSel[i]+CReferSel[i])]
        index0 += CSel[i]
        index1 += CReferSel[i]

    C = allSC[:, CIndexList]
    CRefer = allSC[:, CReferIndexList]

    if (CellsInCol):
        return [C, CRefer, CIndexList, CReferIndexList]
    else:
        return [C.transpose(), CRefer.transpose(), CIndexList, CReferIndexList]


def generateFISHRefer(FISHdata, FiSel,FiReferSel, tumorType='GBM07'):
    '''
    This function is to pick FISH data from each region as Reference in the objective
    Input: FISH dataframe, number of Selected samples, number of Selected refer samples
    Output: selected n sample from each region
    '''
    FISH = FISHdata.iloc[:, 4:].values
    cells = np.sum(FiSel)
    cellsOrig = np.sum(FiReferSel)
       
    #every region contains 150 cell for every patient
    regions = [list(range(0, 150)), list(range(150, 300)),
               list(range(300, FISH.shape[1]))]
    CIndexList = [0] * cells
    CReferIndexList = [0] * cellsOrig

    index0 = 0
    index1 = 0
    for i in range(len(regions)):
        temp = random.sample(regions[i], (FiSel[i]+FiReferSel[i]))
        CIndexList[index0:(index0+FiSel[i])] = temp[0:FiSel[i]]
        CReferIndexList[index1:(index1+FiReferSel[i])
                        ] = temp[FiSel[i]:(FiSel[i]+FiReferSel[i])]
        index0 += FiSel[i]
        index1 += FiReferSel[i]

    C = FISH[:, CIndexList]
    CRefer = FISH[:, CReferIndexList]
    return [C, CRefer, CIndexList, CReferIndexList]


'''
    This function is used to get the initial C matrix.
    Hint: list.extend(list) can be used to merge two lists into one.
'''
def initialC (allSC, cells, method='sc', usedList=[], CellsInCol=True):
    '''Choose 'cells' cells to use as initial values for the copy
    numbers in the deconvolution iterative method(s).

    Usually these will be cells not already chosen for simulation or
    as reference cells, i.e. not the cells is usedList.

    If method == 'random', just choose cells randomly and ignore
    'usedList'.

    '''

    copyNum = allSC.shape[0]
    if (method == 'sc'):
        newList = list(range(allSC.shape[1]))
        for element in usedList:
            newList.remove(element)
        index = random.sample(newList, cells)
        C = allSC[:, index]
    elif (method == 'random'):
        C = np.round(np.random.random([copyNum, cells]) * 10)
    else:
        sys.exit('Use sc or random as method!')

    if CellsInCol:
        return C
    else:
        return C.transpose()

'''
    Calculate Root-mean-square deviation between two matrices.
'''
def calcRMSD(mat1, mat2):
    diff = mat1 - mat2
    sumVal = np.sum(np.multiply(diff, diff)) / (mat1.shape[0] * mat1.shape[1])
    return np.sqrt(sumVal)

'''
    Return the index of truth major single cell components in C matrix.
'''
def findMajorIndex(cellsList, cellsNoiseList):
    index = 0
    majorList = []
    for i in range(len(cellsList)):
        majorList.extend(list(range(index, index+cellsList[i])))
        index += (cellsList[i] + cellsNoiseList[i])
    return majorList


def calcRMSDBySwitch (mat1, mat2, byCol=True):
    if (byCol):
        cells = mat1.shape[1]
        rmsdMat = np.zeros([cells, cells])
        for i in range(cells):
            for j in range(cells):
                diff = mat2[:, i] - mat1[:, j]
                rmsdMat[i, j] = 1.0 * np.sum(np.multiply(diff, diff))
        row, col = linear_sum_assignment(rmsdMat)
        return np.sqrt(np.sum(rmsdMat[row, col])  / (mat1.shape[0] * mat1.shape[1]) )
    else:
        cells = mat1.shape[0]
        rmsdMat = np.zeros([cells, cells])
        for i in range(cells):
            for j in range(cells):
                diff = mat2[i, :] - mat1[j, :]
                rmsdMat[i, j] = 1.0 * np.sum(np.multiply(diff, diff))
        row, col = linear_sum_assignment(rmsdMat)
        return np.sqrt(np.sum(rmsdMat[row, col]) / (mat1.shape[0] * mat1.shape[1]) )

'''
This is calculate the RMS of copy number or fraction in each cell component
'''
def calcRMSInCell(mat1, mat2, Cell=True, CellsInCol=True):
    RMS = []
    if Cell:
        if CellsInCol:
            mat1 = mat1
            mat2 = mat2
        else:
            mat1 = mat1.T
            mat2 = mat2.T

        m, n = mat1.shape
        for i in range(n):
            diff = mat1[:, i] - mat2[:, i]
            sumVal = np.sum(np.multiply(diff, diff)) / float(m)
            RMS.append(np.sqrt(sumVal))
    else:
        if CellsInCol:
            mat1 = mat1.T
            mat2 = mat2.T
        else:
            mat1 = mat1
            mat2 = mat2
        m, n = mat1.shape
        for i in range(n):
            diff = mat1[:, i] - mat2[:, i]
            sumVal = np.sum(np.multiply(diff, diff)) / float(m)
            RMS.append(np.sqrt(sumVal))
    return np.array(RMS)

def CheckDirectory(directory):
    '''
    Check if the directory named 'directory' exists.
    if not, create it.
    '''
    if not os.path.exists(directory):
        try:
            os.makedirs(directory)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
    return

'''
Function to extract the Major Component 
'''
def extraMajorComponent(whole_mat, major_index):
    temp = []
    for i in major_index:
        temp.append(whole_mat[:, i])
    return np.array(temp)
##

'''
Function to add noise to the reference cell
'''
def addNoise(originalMat, NoiseLv):
    if NoiseLv == 0:
        return originalMat
    else:
        m, n = originalMat.shape
        k = int(m * NoiseLv)
        for i in range(n):
            pos = np.random.choice(m, k, replace=False)
            for j in pos:
                flag = np.random.rand()
                if flag < 0.5:
                    originalMat[j,i] = originalMat[j,i] + 1
                else:
                    if originalMat[j,i] == 0:
                        originalMat[j,i] = originalMat[j,i] + 1
                    else:
                        originalMat[j,i] = originalMat[j,i] - 1
        return originalMat

def matchSCSandFISH(SCS, FISH):
    '''
    This function is to match the coordinates of FISH and SCS
    get the index of FISH that is in SCS data
    Input: SCS with coordinate, FISH with coordinates
    Output: a dict that the key is the index in SCS and value is the gene probe
    '''
    SCScoor = SCS[:, 0:3]
    #FISHcoor = FISH[:, 0:3]
    m1 = SCScoor.shape[0]
    m2 = FISH.shape[0]
    ind = dict()
    for i in range(m1):
        for j in range(m2):
            if SCScoor[i][0] == FISH['chromosome'][j]:
                if SCScoor[i][1] <= FISH['begin'][j] and SCScoor[i][2] >= FISH['end'][j]:
                    ind[i] = FISH['probe'][j]
                if SCScoor[i][1] <= FISH['begin'][j] and FISH['begin'][j] <= SCScoor[i][2] and SCScoor[i+1][1] <= FISH['end'][j] and SCScoor[i+1][2] >= FISH['end'][j]:
                    ind[i] = FISH['probe'][j]
                    ind[i+1] = FISH['probe'][j]

    return ind