import numpy as np
import pandas as pd
import pickle
import testFunction
import sys
import argparse
import collections

class Simulation:
    '''
    To simulate bulk and FISH data
    '''

    def __init__(self, data, cellNum, tumorNum, alpha, divider, fishInfo, ploidy, nosie):
        self.scsdata = data[:,3:]            #copy number data of SCS
        self.scsdata[self.scsdata>10] = 10   #cap the copy to 10
        self.interval = data[:, 0:3]         #interval data of sequencing
        self.alpha = alpha                   #parameter for dirichlet distribution 
        self.divider = divider               #divider to separate cells in different region
        self.cellNum = cellNum               #the number of types of clones
        self.tumorNum = tumorNum             #the number of tumor samples
        self.FISH = fishInfo                 #the interval information of FISH data
        self.FISHposition = None             #the index of rows of FISH probes in SCS data
        self.freqs = None                    #simulated frequency
        self.ploidy = None                   #simulated ploidy
    
    
    def simulate_freqs(self):
        '''
        pick 25 cells from each region, 
        cellNum / 3 is the number of dominant cells from the region
        miniCell is the number of other cells excluding dominant cells from the region
        noiseCell is the number of cells from other regions.
        '''
        dormCell_n = [int(self.cellNum/3) ] * 3
        miniCell_n = [int(25 - self.cellNum/3)] * 3
        noiseCell_n = [25] * 3

        self.DirALists = [[self.alpha[0]] * dormCell_n[0] + [self.alpha[1]] *
                           miniCell_n[0], [self.alpha[2]] * noiseCell_n[1], [self.alpha[2]]* noiseCell_n[2]]
        self.Order = [[0, 1, 2], [1, 0, 2], [1, 2, 0]]
        self.freqs = np.empty((75, int(self.tumorNum)))
        self.OrderedDirAList = []        #for FISH simulation use

        for i in range(len(self.Order)):
            temp = []
            for j in range(len(self.DirALists)):
                temp += self.DirALists[self.Order[i][j]]
            self.OrderedDirAList.append(temp)
            self.freqs[:, i] = np.random.dirichlet(temp, int(self.tumorNum/3))


    def simulate_ploidy(self, ploidy):
        '''
        uniformally simulate a ploidy vector that only contains [1,...,8]
        '''
        k = self.scsdata.shape[1]
        if ploidy.isdigit():
            self.ploidy = np.ones(k) * int(ploidy)
        if ploidy == 'random':
            p = [0.1/6, 0.6, 0.1/6, 0.3, 0.1/6, 0.1/6, 0.1/6, 0.1/6]
            self.ploidy = np.random.choice(range(1,9), k, p=p)

    def fishProbePosition(self):
        m1, _ = self.scsdata.shape
        probe_data = self.FISH.values
        m2, _ = probe_data.shape
        #col = self.FISH['probe']
        self.FISHposition = []
        self.maker = []        
        '''
        0: the interval of FISH is included in the interval of SC
        1: the FISH interval cross two interval of SC
        2: the inverval of FISH contain the interval SC
        '''
        for i in range(m2):
            for j in range(m1):
                if self.interval[j, 0] == probe_data[i, 0] and self.interval[j, 1] <= probe_data[i, 1] and self.interval[j, 2] >= probe_data[i, 2]:
                    self.FISHposition.append(j)
                    self.maker.append(0)
                elif self.interval[j, 0] == probe_data[i, 0] and self.interval[j, 1] >= probe_data[i, 1] and self.interval[j, 2] <= probe_data[i, 2]:
                    self.FISHposition.append(j)
                    self.maker.append(2)
                elif self.interval[j, 0] == probe_data[i, 0] and self.interval[j, 1] <= probe_data[i, 1] and self.interval[j, 2] <= probe_data[i, 2] and self.interval[j, 2] > probe_data[i, 1]:
                    self.FISHposition.append(j)
                    self.maker.append(1)

    '''
    the number of columns of SCS data and FISH data, and the length of ploidy are the same (SCS sample number)
    '''

    def generateSelection(self):
        '''
        select the index for simulation 
        '''
        dormCell_n = [int(self.cellNum/3)] * 3
        miniCell_n = [int(25 - self.cellNum/3)] * 3
        noiseCell_n = [25] * 3

        regions = [list(range(self.divider[0])), list(range(self.divider[0], self.divider[1])),
                   list(range(self.divider[1], self.scsdata.shape[1]))]

        self.OriginCIndex = [0] * sum(noiseCell_n)

        index0 = 0
        for i in range(len(regions)):
            temp = np.random.choice(regions[i], noiseCell_n[i] + dormCell_n[i]*2, replace=False)
            self.OriginCIndex[index0:(
                index0+noiseCell_n[i])] = temp[0:noiseCell_n[i]]
        
            index0 += noiseCell_n[i]

        '''
        retrieve the index of dorminant clone from the OriginalIndex, 
        naturally, it is the first two of each region
        '''
        index4 = 0
        self.majorList = []
        for i in range(len(dormCell_n)):
            self.majorList.extend(list(range(index4, index4+dormCell_n[i])))
            index4 += (dormCell_n[i] + miniCell_n[i])

        self.OriginC = self.scsdata[:, self.OriginCIndex]
        self.TrueCIndex = [self.OriginCIndex[x] for x in self.majorList]
        self.TrueFreqs = self.freqs[self.majorList, :]
        self.TrueFreqs = self.TrueFreqs/np.sum(self.TrueFreqs, axis=0)
        # ploidy of the major clones in each region
        # ploidy of the samples used to simulate Bulk
        self.OriginPlodiy = [self.ploidy[x] for x in self.OriginCIndex]
        self.TruePloidy = [self.OriginPlodiy[x] for x in self.majorList]



    def simulateFISH(self):
        '''
        simulateFISH according to the simulated frequency and selected cells
        '''
        #--------------------
        #retrieve the copy number according to the selected real SCS
        _, k = self.OriginC.shape
        probe_data = self.FISH.values
        m, _ = probe_data.shape
        self.FISH = np.zeros((m, k))

        for i in range(m):
            for j in range(k):
                if self.maker[i] == 0:
                    self.FISH[i, j] = self.OriginC[self.FISHposition[i], j]
                if self.maker[i] == 1:
                    '''
                    take the copy to the portion of the length and then average:
                    p = ((x-y1) * p1 + (y2-x) * p2) /(y2-y1), 
                        x is the end position in the SC file for current row
                        y1 is the start position in the probe file for the current row
                        y2 is the end position in the proble file for the current row
                        p1 is the CNV of current row, p2 is the CNV of next row
                    '''
                    x = self.interval[self.FISHposition[i], 2]
                    y1 = probe_data[i, 1]
                    y2 = probe_data[i, 2]
                    p1 = self.OriginC[self.FISHposition[i], j]
                    p2 = self.OriginC[self.FISHposition[i] + 1, j]

                    self.FISH[i, j] = ((x-y1) * p1 + (y2-x) * p2) / (y2-y1)

        #--------------------
        #simulate the proportion (count) of each FISH sample according to frequency
        #sampling the FISH cell from multinomial distrition
        #dirichlet distribution is conjugate to Multinomial Distribution
        _, n = self.freqs.shape
        self.counts = []                              #number of sampling FISH cells
        self.OriginRefFreqs = np.empty((k, n))        #fraction of sampling FISH cells
        for i in range(n):
            prior = self.freqs[:,i]
            multi = np.random.multinomial(1000, prior) #sampling 1000 FISH cells for each region
            prob = [i/sum(multi) for i in multi]
            #posterior = np.random.dirichlet(list(multi+np.array(self.OrderedDirAList[i])), int(self.tumorNum/3))
            #self.OriginRefFreqs[:, i] = posterior
            self.OriginRefFreqs[:, i] = prob
            # FISH count for each FISH sample would be the reference frequency
            # the most frequent types would be selected for each region
            self.counts.append(list(multi))

        '''
        reference frequence come from the most frequent FISH cell in each region
        which is the index in majorList
        (*) we do not normalize the reference frequence 
        '''
        counts = np.array([x/sum(x) for x in self.counts])
        self.TrueRefFreqs = counts[:, self.majorList].T

    def simulateCell(self, noise):
        '''
        simulate SCS with diploidy for cells to be 
        selected as reference and initialization
        '''
        tempOriginC = self.scsdata[:, self.OriginCIndex]
        m, k = tempOriginC.shape
        _, n = self.freqs.shape
        # repeatTimes = 100 # repeat the selection 100 times and calculate the average
        #self.referC, self.initialC, self.referPloidy, self.initialPloidy = averageCPN(repeatTimes, k, self.cellNum//3, n, self.freqs, tempOriginC, self.OriginPlodiy)

        self.referC, self.initialC, self.referPloidy, self.initialPloidy = [], [], [], []
        for i in range(n):
            prior = self.freqs[:, i]
            multi = np.random.multinomial(1000, prior)
            p = multi / sum(multi)
            index = np.random.choice(range(k), self.cellNum//3*2, p=p)
            for j in range(2):
                self.referC.append(tempOriginC[:, index[2*j]])
                self.initialC.append(tempOriginC[:, index[2*j+1]])
                self.referPloidy.append(self.OriginPlodiy[index[2*j]])
                self.initialPloidy.append(self.OriginPlodiy[index[2*j+1]])

        self.referC = np.array(self.referC).T
        self.initialC = np.array(self.initialC).T

        #add noise to the copy number
        p = noise
        for i in range(m):
            for j in range(self.cellNum//3*2):
                if tempOriginC[i, j] < 1:
                    self.referC[i, j] += np.random.choice([0,1],1, p=[1-p, p])[0]
                    self.initialC[i, j] += np.random.choice([0,1],1, p=[1-p, p])[0]
                else:
                    self.referC[i, j] += np.random.choice([-1, 0, 1], 1, p=[p, 1-p*2, p])[0]
                    self.initialC[i, j] += np.random.choice([-1, 0, 1], 1, p=[p, 1-p*2, p])[0]

        '''
        cap the copy number to be 10 if copy number > 10
        '''
        self.referC[self.referC > 10] = 10
        self.initialC[self.initialC > 10] = 10


    def simulateReferFISH(self):
        '''
        (12.20) We should pick most frequent types for each region
        which are the first two clone of each region
        '''
        #k = self.FISH.shape[1]
        #n = self.TrueFreqs.shape[1]
        # the index is 0, 1, 25, 26, 50, 51
        idx = self.majorList
        counts = self.counts
        self.referFISHploidy = [self.OriginPlodiy[i] for i in idx]
        self.referFISH = self.FISH[:, idx]

        




    def addNOiseFISH(self, noise):
        self.NoiseFISH = np.zeros((self.referFISH.shape[0], self.referFISH.shape[1]))
        self.NoiseFISH[-1, :] = self.referFISH[-1, :]
        p = noise   #add the same level of noise to FISH
        for i in range(self.referFISH.shape[0]-1):
            for j in range(self.referFISH.shape[1]):
                if self.referFISH[i, j] == 0:
                    self.NoiseFISH[i, j] = self.referFISH[i, j] + np.random.choice([0,1],1, p=[1-p, p])[0]
                else:
                    self.NoiseFISH[i, j] = self.referFISH[i, j] + np.random.choice([-1, 0, 1], 1, p=[p, 1-p*2, p])[0]
        
        # cap the copy number of FISH to be 10 after adding noise
        self.NoiseFISH[self.NoiseFISH>10] = 10

    def corr_feature(self, threshold=0.95):
        ''' calculate the correlation matrix for the feature 
            make the features around FISH probe similar        
        '''
        corrMat = np.corrcoef(self.scsdata)
        self.link = []
        for p in self.FISHposition:
            t = []
            corrPos = corrMat[0:, p]
            i1, i2 = p, p
            while True:
                if i1 >=0 and corrPos[i1] >= threshold:
                    t.append(i1)
                    i1 -= 1  #go forward
                elif i2 >=0 and corrPos[i2+1] >= threshold:
                    t.append(i2+1)
                    i2 += 1 #go back ward
                else:
                    break
            self.link.append(sorted(t))
        
    def enrichReferFISHMat(self, noise):
        '''
        expand the referFISH matrix to the neighbouring positions
        '''
        self.enrichReferFISH = []
        for i, item in enumerate(self.link):
            for _ in range(len(item)):
                self.enrichReferFISH.append(self.referFISH[i, :])
        self.enrichReferFISH = np.array(self.enrichReferFISH)
        p = noise   #add the same level of noise to FISH
        for i in range(self.enrichReferFISH.shape[0]):
            for j in range(self.enrichReferFISH.shape[1]):
                if self.enrichReferFISH[i, j] == 0:
                    self.enrichReferFISH[i, j] += np.random.choice([0,1],1, p=[1-p, p])[0]
                else:
                    self.enrichReferFISH[i, j] += np.random.choice([-1, 0, 1], 1, p=[p, 1-p*2, p])[0]
        
        # cap the copy number of FISH to be 10 after adding noise
        self.enrichReferFISH[self.enrichReferFISH > 10] = 10



    def simulateBulk(self):
        '''
        simulate the bulk tumor from C, ploidy(p) and F
        ploidy is the rescale factor of the cell samples:
        B = C*diag(p)/2*F
        '''        
        p = np.diag(self.OriginPlodiy)
        self.bulkTumor = np.dot(np.dot(self.OriginC, p/2), self.freqs)

def mostFrequent(n, p, k, m):
    '''return the most frequent items for FISH selections'''
    counts = collections.Counter({})
    for _ in range(n):
        count = collections.Counter(np.random.choice(range(k), m, p))
        counts += count
    counts = dict(counts)
    val = list(counts.values())
    val.sort(reverse=True)
    idx = [key for i in range(2) for key in counts.keys() if counts[key]==val[i]]
    return idx

def averageCPN(n, k, m, r, Freq, CPNdata, ploidyData):
    '''repeat the selection and average he cpn and ploidy as the final values for reference'''
    refCPN, iniCPN = np.zeros((CPNdata.shape[0], m*r)), np.zeros((CPNdata.shape[0], m*r))
    ref_ploidy, init_ploidy = [0]*m*r, [0]*m*r
    for ri in range(r):
        prior = Freq[:, ri]
        multi = np.random.multinomial(1000, prior)
        p = multi / sum(multi)
        for _ in range(n):
            index = np.random.choice(range(k), m*2, p=p)
            for j in range(2):
                temprefC = CPNdata[:, index[2*j]]
                tempinitC = CPNdata[:, index[2*j+1]]
                temprefP = ploidyData[index[2*j]]
                tempinitP = ploidyData[index[2*j+1]]
                ref_ploidy[2*ri+j] += temprefP
                init_ploidy[2*ri+j] += tempinitP
                refCPN[:, [2*ri+j]] += (temprefC * temprefP / 2).reshape((-1,1))
                iniCPN[:, [2*ri+j]] += (tempinitC * tempinitP / 2).reshape((-1,1))
    
    ref_ploidy = [item / n for item in ref_ploidy]
    init_ploidy = [item / n for item in init_ploidy]
    for i in range(m*r):
        refCPN[:, i] = refCPN[:, i] / n * 2 / ref_ploidy[i]
        iniCPN[:, i] = iniCPN[:, i] / n * 2 / init_ploidy[i]
    return refCPN, iniCPN, ref_ploidy, init_ploidy

def readSCS(path):
    return np.loadtxt(path, delimiter=',')

def readFISHInterval(filename):
    fil = []
    gene = []
    with open(filename) as txtfile:
        for i in txtfile:
            i = i.strip('\n').split('\t')
            fil.append(i[1:])
            gene.append(i[0])
    col = ['chr', 'start', 'end']
    df = pd.DataFrame(np.array(fil[1:]).astype(int), columns=col)
    df['probe'] = gene[1:]
    return df

def saveData(obj, filename, i):
    output_file = '%s/simulate_data_%d.pkl' % (filename, i)
    pickle_data = {
        # the simulated bulk matrix (an np.array with genes as rows)
        'TumorSample': obj.bulkTumor,
        # copy number of the dominant cells only (transposed)
        'CTrue': obj.scsdata[:, obj.TrueCIndex],
        # several reference cells used in some solvers
        'CRefer': obj.referC,
        # Mixture fractions of the dominant cells only.
        'FTrue': obj.TrueFreqs,
        # Mixture fractions used to generated B.
        'FTrueAll': obj.freqs,
        'FRefer': obj.TrueRefFreqs,
        # a (possible) starting point for deconvolution algorithms
        'CInitial': obj.initialC,
        'COrigin': obj.scsdata[:, obj.OriginCIndex],  # Copy numbers used to generate B
        'PloidyOrigin': obj.OriginPlodiy, #ploidy of selected cell to simulate B
        'TruePloidy': obj.TruePloidy,     #ploidy of major clones in each region
        'ReferPloidy': obj.referPloidy,   #ploidy of reference clones
        'InitialPloidy': obj.initialPloidy, #ploidy of initialization
        #'FISHRefer': obj.FISH,   # serverl reference FISH probe used
        'FISHcounts':obj.counts, # count for each FISH cell in the region
        'SimulatedFISH': obj.FISH,  # simulated FISH corresponding to SCS samples
        'FISHRefer': obj.NoiseFISH, #FISH cell for reference (with noise)
        'referFISHploidy': obj.referFISHploidy, #ploidy for the reference FISH cell
        'ProbeIndex': obj.FISHposition,   # FISH probes's corresponding indexes in SCS
        #'ProbeName': np.asanyarray(ind_values,dtype=str) #FISH probes
        'EnrichReferFISH': obj.enrichReferFISH,   #expand FISH mat for reference      
        'correlatePositions': obj.link,    #
        'uFISHRefer': np.dot(obj.NoiseFISH, np.diag(obj.referFISHploidy)) / 2, 
        'uEnrichReferFISH': np.dot(obj.enrichReferFISH, np.diag(obj.referFISHploidy)) / 2


    }
    output = open(output_file, 'wb')
    pickle.dump(pickle_data, output)
    output.close()

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Simulate bulk tumor')
    parser.add_argument('scsdata_directory')
    parser.add_argument('fishdata_directory')
    parser.add_argument('save_directory')
    parser.add_argument('--date')
    parser.add_argument('--tumorName', type=str, default='GBM07')
    parser.add_argument('--cellNums', type=int, default=6)
    parser.add_argument('--tumorNums', type=int, default=3)
    parser.add_argument('--simuNums', type=int, default=20)
    parser.add_argument('--cellNoise', type=float, default=0.0)
    parser.add_argument('--ploidy', type=str, default=2)
    return parser.parse_args()


def simulateInstance(output_dir, scsData, cellNum, tumorNum, alpha, divider, fishInfo, ploidy, noise, i):
    s = Simulation(scsData, cellNum, tumorNum, alpha, divider, fishInfo, ploidy, noise)
    #simulate frequency, ploidy and FISH data
    s.simulate_freqs()
    s.simulate_ploidy(ploidy)
    s.fishProbePosition()
    s.generateSelection()
    s.simulateFISH()
    s.simulateCell(noise)
    s.simulateReferFISH()
    s.addNOiseFISH(noise)
    s.corr_feature()
    s.enrichReferFISHMat(noise)
    s.simulateBulk()
    saveData(s, output_dir, i)

def main():
    args = parse_args()
    scsData = readSCS(args.scsdata_directory)
    fishInfo = readFISHInterval(args.fishdata_directory)
    output = args.save_directory
    DateFolder = args.date
    TumorName = args.tumorName
    cellNum = args.cellNums
    tumorNum = args.tumorNums
    N = args.simuNums
    noise = args.cellNoise
    ploidy = args.ploidy
    output_dir = '%s/%s/%s/%s' % (output, DateFolder, TumorName, tumorNum)
    testFunction.CheckDirectory(output_dir)


    alpha = [100, 1, 0.1]
    if args.tumorName == 'GBM07':
        divider = [57, 132]
    elif args.tumorName == 'GBM33':
        divider = [59, 135]
    else:
        exit("Please choose GBM07 or GBM33")
        
    for i in range(N):
       simulateInstance(output_dir, scsData, cellNum, tumorNum,
                        alpha, divider, fishInfo, ploidy, noise, i)

if __name__ == '__main__':
    main()
    
