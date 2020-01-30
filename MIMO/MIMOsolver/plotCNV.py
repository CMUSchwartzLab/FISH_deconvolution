import seaborn as sns
import numpy as np
import argparse
import glob
import pandas as pd
import re
import itertools
import collections
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import testFunction
from deconvsolver import Simulation
from scipy.optimize import linear_sum_assignment
import os


PROBLEM_INDEX = re.compile('simulate_data_(.*?).pkl')

'''
need to retrieve data from simulated and result
'''
def readResult(filename):
    dic = collections.defaultdict(list)
    level = -1
    label = ['cnv', 'freqs', 'ploidy', 'tree']
    try: 
        with open(filename) as f:
            for item in f:       
                if '---' in item:
                    level += 1
                else:
                    p = re.compile(r'\d+.\d+|\d+')
                    item = item.split(',')
                    levelList = []
                    for i in item:
                        b = p.findall(i)
                        if len(b) > 1:
                            b[0] = float(b[0]) * (10**-float(b.pop()))  #deal with e
                        levelList.append(float(b[0]))
                    tag = label[level]
                    dic[tag].append(levelList)
        return dic

    except:

        return None
        
    
'''
re-order based on ploidy
'''
def arrangeCellsOnP(mat1, mat2):
    k = mat1.shape[0]
    distMat = np.zeros((k, k))
    for i in range(k):
        for j in range(k):
            distMat[i, j] = np.sum(abs(mat1[i] - mat2[j]))
    _, order = linear_sum_assignment(distMat)
    return order

'''
re-order based on Cell CNV
'''
def arrangeCellsOnC(mat1, mat2):
    k = mat1.shape[1]
    distMat = np.zeros((k, k))
    for i in range(k):
        for j in range(k):
            distMat[i, j] = np.sum(abs(mat1[:,i] - mat2[:,j]))
    _, order = linear_sum_assignment(distMat)
    return order


def calcAccuracy(mat1, mat2, threshold=0.5):
    m, k = mat1.shape
    acc = np.zeros(k)
    for i in range(k):
        acc[i] = np.sum(abs(mat1[:, i] - mat2[:, i]) <= threshold)/m
    return list(acc), [np.sum(acc)/k]

def calcRMSDInCell(mat1, mat2):
    rms = []
    m, k = mat1.shape
    for i in range(k):
        diff = mat1[:, i] - mat2[:, i]
        sumVal = np.sum(np.multiply(diff, diff)) / float(m)
        rms.append(np.sqrt(sumVal))
    return rms

def calcRMSDInFreqs(mat1, mat2):
    rms = []
    k, n = mat1.shape
    for i in range(k):
        diff = mat1[i, :] - mat2[i, :]
        sumVal = np.sum(np.multiply(diff, diff)) / float(n)
        rms.append(np.sqrt(sumVal))
    return rms

def calcRMSD_P(mat1, mat2):
    diff = mat1 - mat2
    sumVal = np.sum([i*i for i in diff])/mat1.shape[0]
    return [np.sqrt(sumVal)]

def calcRMSDInPloidy(mat1, mat2):
    rms = []
    k = mat1.shape[0]
    for i in range(k):
        diff = mat1[i] - mat2[i]
        sumVal = diff * diff
        rms.append(np.sqrt(sumVal))
    return rms

        
def calcRMSD(mat1, mat2):
    diff = mat1 - mat2
    sumVal = np.sum(np.multiply(diff, diff)) / \
        (mat1.shape[0] * mat1.shape[1])
    return [np.sqrt(sumVal)]


def paras(parameters, val):
    '''
    get all the combinations of parameters
    '''
    n = len(parameters[0])
    paras = parameters
    for i in range(1,n+1):
        c = itertools.combinations(range(n), i)
        for item in c:
            temp = [0.0] * n
            for j in item:
                temp[j] = val
            paras.append(temp)
    return parameters

def loadSimulatedData(sample):
    simulation = Simulation.read(sample)
    true_cells, true_freqs = simulation.true_cells, simulation.true_frequencies
    true_ploidy = simulation.true_cells_ploidy
    return true_cells, true_freqs, np.array(true_ploidy)

def find_simulated_datafiles(dirs):
    simulated_datafiles = list(
        sorted(glob.glob('%s/*.pkl' % (dirs))))
    return simulated_datafiles

def calcAverge(paras, simulate_files, dirs, reorder='c'):
    totalData = []
    totalErr = []
    for i in range(len(paras)):
        dicData = collections.defaultdict(list)
        err = {}
        for j, sample in enumerate(simulate_files):
            prob_id = os.path.basename(sample)
            m = PROBLEM_INDEX.match(prob_id)
            if m:
                prob_id = m.group(1)
    
            input_name = '%ssample_output_%s' % (dirs, prob_id)
            for item in paras[i]:
                input_name += ('_' + str(item))
            input_name += '.txt'
            dic = readResult(input_name)
            true_cells, true_freqs, true_ploidy = loadSimulatedData(sample)
            if dic is not None:
                infercell = np.array(dic['cnv'])
                inferF = np.array(dic['freqs'])
                inferP = np.array(dic['ploidy'][0])
                if reorder == 'c':
                    order = arrangeCellsOnC(true_cells, infercell)
                else:
                    order = arrangeCellsOnP(true_ploidy, inferP)
        
                dicData['overall copy number accuracy'].append(calcAccuracy(true_cells, infercell[:, order])[1][0])
                dicData['overall RMSD of copy number'].append(calcRMSD(true_cells, infercell[:, order])[0])
                dicData['overall RMSD of freqs'].append(calcRMSD(true_freqs, inferF[order, :])[0])
                dicData['overall RMSD of ploidy'].append(calcRMSD_P(true_ploidy, inferP[order])[0])
                '''
                calculate the RMSD between unnormalized copy number
                '''
                diagPTrue = np.diag(true_ploidy)
                diagPInfer = np.diag(inferP[order])
                utrueCNV = np.dot(true_cells, diagPTrue) / 2
                uinferCNV = np.dot(infercell[:,order], diagPInfer) / 2
                dicData['overall RMSD of unnormalized copy number'].append(calcRMSD(utrueCNV, uinferCNV)[0])

        for key in dicData.keys():
            err[key] = np.std(dicData[key])
            dicData[key] = sum(dicData[key])/len(dicData[key])
        totalData.append(dicData)
        totalErr.append(err)
    return totalData, totalErr

def plot(data, err, labels, path, reorder):
    plt.style.use('classic')
    plt.rcParams['figure.facecolor'] = '1'
    keys = list(data.keys())
    ylabel = ['Accuracy', 'RMSD', 'RMSD', 'RMSD', 'RMSD']
    my_colors = ['red', 'lime', 'violet', 'lightgrey', 'blue', 'orange','cyan', 'coral']
    fig = plt.figure(figsize=(60, 6))
    #fig = plt.figure(figsize=(12, 12))
    for i in range(len(keys)):
        key = keys[i]
        subData = data[key]
        subErr =  err[key]
        ax = fig.add_subplot(1,6,i+1)
        #ax = fig.add_subplot(3,2,i+1)
        ax.bar(x=range(1, len(subData)+1), height=subData, yerr=subErr, color=my_colors)
        plt.axhline(y=0, color='k')
        if i ==0:
            g = 0.02
            y1 = subData[4]+0.012 + subErr[4]
            y2 = subData[5]+0.012 + subErr[5]
            y3 = subData[6]+0.012 + subErr[6]
            y4 = subData[7]+0.012 + subErr[7]
        
        elif i == 2:
            g = 0.01
            y1 = subData[4]+0.01 + subErr[4]
            y2 = subData[5]+0.01 + subErr[5]
            y3 = subData[6]+0.01 + subErr[6]
            y4 = subData[7]+0.01 + subErr[7]
        elif i == 1 or i == 3:
            g = 0.1
            y1 = subData[4]+0.055 + subErr[4]
            y2 = subData[5]+0.055 + subErr[5]
            y3 = subData[6]+0.055 + subErr[6]
            y4 = subData[7]+0.055 + subErr[7]
        else:
            g = 0.5
            y1 = subData[4]+0.5 + subErr[4]
            y2 = subData[5]+0.5 + subErr[5]
            y3 = subData[6]+0.5 + subErr[6]
            y4 = subData[7]+0.5 + subErr[7]

        # label the mean
        # ax.text(4+0.6, 0 - 4*g, "%.3f" % subData[4] , color='k', fontweight='bold', fontsize=9)
        # ax.text(5+0.6, 0 - 4*g, "%.3f" % subData[5], color='k', fontweight='bold', fontsize=9)
        # ax.text(6+0.6, 0 - 4*g, "%.3f" % subData[6], color='k', fontweight='bold', fontsize=9)
        # ax.text(7+0.6, 0 - 4*g, "%.3f" % subData[7], color='k', fontweight='bold', fontsize=9)
        ax.text(4+0.8, 0 - 4*g, "%.3f" % subData[4] , color='k', fontweight='bold', fontsize=9)
        ax.text(5+0.8, 0 - 4*g, "%.3f" % subData[5], color='k', fontweight='bold', fontsize=9)
        ax.text(6+0.8, 0 - 4*g, "%.3f" % subData[6], color='k', fontweight='bold', fontsize=9)
        ax.text(7+0.8, 0 - 4*g, "%.3f" % subData[7], color='k', fontweight='bold', fontsize=9)
        #label the std
        # ax.text(4+0.5, y1, "(%.3f)" % subErr[4], color='k', fontweight='bold', fontsize=8)
        # ax.text(5+0.5, y2, "(%.3f)" % subErr[5], color='k', fontweight='bold', fontsize=8)
        # ax.text(6+0.5, y3, "(%.3f)" % subErr[6], color='k', fontweight='bold', fontsize=8)
        # ax.text(7+0.5, y4, "(%.3f)" % subErr[7], color='k', fontweight='bold', fontsize=8)
        ax.text(4+0.7, y1, "(%.3f)" % subErr[4], color='k', fontweight='bold', fontsize=8)
        ax.text(5+0.7, y2, "(%.3f)" % subErr[5], color='k', fontweight='bold', fontsize=8)
        ax.text(6+0.7, y3, "(%.3f)" % subErr[6], color='k', fontweight='bold', fontsize=8)
        ax.text(7+0.7, y4, "(%.3f)" % subErr[7], color='k', fontweight='bold', fontsize=8)
        
        low = min(0.0, np.min(np.array(subData) - np.array(subErr)))
        low = low - 6 * g
        plt.ylim(low, np.max(np.array(subData) + np.array(subErr)) * 1.1)
        ax.title.set_text(keys[i])
        ax.set_ylabel(ylabel[i])
        ax.set_xlabel('model')
        ax.get_xaxis().set_visible(False)
        plt.grid()
    
    legend_elements = []
    af = '\\alpha_{f}'
    a1 = '\\alpha_{1}'
    a2 = '\\alpha_{2}'
    for i in range(len(labels)):
        legend_elements.append(Patch(facecolor=my_colors[i] ,label=', '.join([str(i) for i in labels[i]])))
    fig.legend(handles=legend_elements,  title="Parameters values: \n           $%s$,    $%s$,    $%s$" %(af, a1, a2), loc=(0.855, 0.38))  
    #fig.legend(handles=legend_elements,  title="Parameters values: \n           $%s$,    $%s$,    $%s$" %(af, a1, a2), loc=(0.65, 0.06))
    #fig.text(0.5, -0.04, 'models', size='15')
    plt.savefig('%s/overallAver_%s.eps' % (path, reorder), bbox_inches='tight', format='eps', dpi=100)

def reconstructTable(table):
    '''
    each key contain the result for each parameter combination
    '''
    storedKey = ['overall copy number accuracy',
                 'overall RMSD of copy number', 
                 'overall RMSD of freqs',
                 'overall RMSD of ploidy',
                 'overall RMSD of unnormalized copy number']

    cTable = collections.defaultdict(list)
    for key in storedKey:
        for item in table:
            cTable[key].append(item[key])
    
    return cTable

def parseParameter(inputDir):
    '''
    detect all the combinations of parameter from the file name
    '''
    seen = set()
    paras = []
    p = re.compile(r'\d*\.\d')
    paths = glob.glob('%sstastics_output_*' % inputDir)
    for path in paths:
        path = p.findall(path)
        para = [float(val) for val in path[-3:]]
        seen.add(tuple(para))
    for item in seen:
        paras.append(item)
    return paras


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='plot the result')
    parser.add_argument('simu_dir', type=str)
    #parser.add_argument('simu_date', type=str)
    parser.add_argument('input_dir')
    parser.add_argument('output_dir')
    parser.add_argument('out_date', type=str)
    parser.add_argument('tumor_name', type=str)
    parser.add_argument('cell_num', type=int)
    parser.add_argument('--parameter_num', type=int)
    parser.add_argument('--parameter_val', type=float)
    parser.add_argument('--reorder_option', type=str)
    return parser.parse_args()

def main():
    args = parse_args()
    #out_date = args.out_date
    #num = args.parameter_num
    #val = args.parameter_val
    simulate_path = args.simu_dir + '/' + args.tumor_name + '/' + str(args.cell_num) + '/'
    dataDir = args.input_dir 
    output =  args.output_dir 
    reorder= args.reorder_option
    testFunction.CheckDirectory(output)
    simulated_datafiles = find_simulated_datafiles(simulate_path)
    #parameters = [[0.0] * num]
    #parameters = paras(parameters, val)
    parameters = parseParameter(dataDir)
    parameters.sort(key=lambda x: sum(x))
    parameters.sort(key=lambda x: x[1])
    print(parameters)
    totalData, totalError = calcAverge(parameters, simulated_datafiles, dataDir, reorder)
    totalData = reconstructTable(totalData)
    totalError = reconstructTable(totalError)
    plot(totalData, totalError, parameters, output, reorder)

if __name__ == '__main__':
    main()