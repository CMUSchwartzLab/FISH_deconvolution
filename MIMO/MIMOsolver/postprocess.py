'''
post-prcess the results:
. rearrange the columns of cell
. rearrange the columns of frequencies
. calculate the overall accuracy of copy number
. calculate the accuracy of each cell component
. calculate the RMSD of copy number
. calculate the RMSD of frequencies
. calculate the pairwise RMSD of samples in cell
. calculate the pairwise RMSD of samples in frequencies
...

'''

import numpy as np
from scipy.optimize import linear_sum_assignment


class PProcess:
    def __init__(self, trueCells, inferredCells, trueFreqs, inferredFreqs, truePloidy, inferredPloidy):
        self.trueCells = trueCells
        self.inferredCells = inferredCells
        self.trueFreqs = trueFreqs
        self.inferredFreqs = inferredFreqs
        self.truePloidy = np.array(truePloidy)
        self.inferredPloidy = np.array(inferredPloidy)
        self.order = None

    '''
    re-arrange the order of cells according to the distance
    '''
    def arrangeCells(self):
        _, k = self.trueCells.shape
        distMat = np.zeros((k, k))
        for i in range(k):
            for j in range(k):
                distMat[i, j] = np.sum(abs(self.trueCells[:, i] - self.inferredCells[:, j]))
        _, self.order = linear_sum_assignment(distMat)
        self.inferredCells = self.inferredCells[:, self.order]

    '''
    re-arrange the order of frequencies according the new order from cells
    '''
    def arrangeFreqs(self):
        self.inferredFreqs = self.inferredFreqs[self.order, :]

    '''re-arrange the order of the ploidy'''
    def arrangePloidy(self):
        self.inferredPloidy = self.inferredPloidy[self.order]

    '''
    calculate the accuracy of copy number in each cell component and
    overall accuracy
    '''
    def calcAccuracy(self, threshold=0.5):
        m, k = self.trueCells.shape
        acc = np.zeros(k)
        for i in range(k):
            acc[i] = np.sum(abs(self.trueCells[:, i] -
                                      self.inferredCells[:, i]) <= threshold)/m
        return list(acc), [np.sum(acc)/k]

    '''
    calculate the RMSD in each cell
    '''
    def calcRMSDInCell(self):
        rms = []
        m, k = self.trueCells.shape
        for i in range(k):
            diff = self.trueCells[:, i] - self.inferredCells[:, i]
            sumVal = np.sum(np.multiply(diff, diff)) / float(m)
            rms.append(np.sqrt(sumVal))
        return rms

    '''
    calculate the RMSD in each frequencies
    '''
    def calcRMSDInFreqs(self):
        rms = []
        k, n = self.trueFreqs.shape
        for i in range(k):
            diff = self.trueFreqs[i, :] - self.inferredFreqs[i, :]
            sumVal = np.sum(np.multiply(diff, diff)) / float(n)
            rms.append(np.sqrt(sumVal))
        return rms

    '''
    calculate the overall RMSD of copy number
    '''
    def calcRMSD_C(self):
        return calcRMSD(self.trueCells, self.inferredCells)

    '''
    calculate the overall RMSD of frequencies
    '''
    def calcRMSD_F(self):
        return calcRMSD(self.trueFreqs, self.inferredFreqs)
    
    '''
    calculate the RMSD in each ploidy
    '''
    def calcRMSD_P(self):
        diff = self.truePloidy - self.inferredPloidy
        sumVal = np.sum([i*i for i in diff])/self.truePloidy.shape[0]
        return [np.sqrt(sumVal)]
    
    def calcRMSDInPloidy(self):
        rms = []
        k = self.truePloidy.shape[0]
        for i in range(k):
            diff = self.truePloidy[i] - self.inferredPloidy[i]
            sumVal = diff * diff
            rms.append(np.sqrt(sumVal))
        return rms

        


'''
calculate the overall RMSD of two matrix
'''


def calcRMSD(mat1, mat2):
    diff = mat1 - mat2
    sumVal = np.sum(np.multiply(diff, diff)) / \
        (mat1.shape[0] * mat1.shape[1])
    return [np.sqrt(sumVal)]


'''
write the result
'''
def yaml_dump_matrix(lst, outfile):
    ncol = len(lst)
    print(' - [ ', file=outfile, end='')
    for j in range(ncol):
        if j != 0:
            print(', ', file=outfile, end='')
        print(lst[j], file=outfile, end='')
    print(' ]', file=outfile)


def write_result(cl, output):
    cl.arrangeCells()
    cl.arrangeFreqs()
    print('--- overall copy number accuracy ---', file=output)
    acc, overAcc = cl.calcAccuracy()
    yaml_dump_matrix(overAcc, output)
    print('--- copy number accuracy in cell ---', file=output)
    yaml_dump_matrix(acc, output)

    print('--- overall RMSD of copy number ---', file=output)
    rmsdC = cl.calcRMSD_C()
    yaml_dump_matrix(rmsdC, output)

    print('--- RMSD of copy number in cells ---', file=output)
    rmsdInC = cl.calcRMSDInCell()
    yaml_dump_matrix(rmsdInC, output)

    print('--- overall RMSD of freqs ---', file=output)
    rmsdF = cl.calcRMSD_F()
    yaml_dump_matrix(rmsdF, output)

    print('--- RMSD of freqs in cells ---', file=output)
    rmsdInF = cl.calcRMSDInFreqs()
    yaml_dump_matrix(rmsdInF, output)

    print('--- overall RMSD of ploidy ---', file=output)
    rmsdP = cl.calcRMSD_P()
    yaml_dump_matrix(rmsdP, output)

    print('--- RMSD of ploidy in cells ---', file=output)
    rmsdInP = cl.calcRMSDInPloidy()
    yaml_dump_matrix(rmsdInP, output)


'''
main function to run all the methods above and save data
'''
def calculations(trueCells, inferredCells, trueFreqs, inferredFreqs, truePloidy, inferredPloidy, date, postfix):
    pi = PProcess(trueCells, inferredCells, trueFreqs, inferredFreqs, truePloidy, inferredPloidy)
    basename = '../results/' + date + '/stastics_output_' + postfix
    with open(basename + '.txt', 'wt') as outfile:
        write_result(cl=pi, output=outfile)
    outfile.close()
