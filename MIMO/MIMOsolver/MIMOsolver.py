from gurobipy import GRB, Model
import numpy as np
import time
import random
import pickle
import testFunction as test
import os

'''
Update Fraction Matrix.
'''
def updateProportion(B, Ctotal, pr, FRefer, alphaf, cells=None, root=0, vType='C'):
    tumors = B.shape[0]
    cellsTotal = Ctotal.shape[0]
    copyNum = Ctotal.shape[1]

    if cells == None:
        cells = cellsTotal
    m = Model('phylo')
    m.setParam("OutputFlag", 0)
    F = m.addVars(tumors, cells, name='F', vtype=vType)
    bDelta = m.addVars(tumors, copyNum, name='bDelta', vtype='C')
    for p in range(tumors):
        for s in range(copyNum):
            expr = B[p, s]
            for k in range(cells):
                expr += -(F[p, k]*pr[k, k]/2*Ctotal[k, s])
            m.addConstr(bDelta[p, s]-(expr) >= 0,
        name='residu(%s,%s)' % (p, s))
            m.addConstr(bDelta[p, s]+(expr) >= 0,
        name='residl(%s,%s)' % (p, s))

    for p in range(tumors):
        for k in range(cells):
            m.addConstr(F[p, k] <= 1, name='Fupp(%s,%s)' % (p,k))
            m.addConstr(F[p, k] >= 10 ** (-4),
        name='Flow(%s,%s)' % (p,k))
    m.addConstrs((F.sum(p, '*') == 1 for p in range(tumors)),
    name='Fsum')

    # if or not impose the constraints on F
    if alphaf > 0:
        fDelta = m.addVars(tumors, cells, name='fDelta', vtype=vType)
        for p in range(tumors):
            for s in range(cells):
                expr = F[p, s] - FRefer[p, s]
                m.addConstr(
                    fDelta[p, s] - (expr) >= 0, name='fesu(%s,%s)' % (p, s))
                m.addConstr(
                    fDelta[p, s] + (expr) >= 0, name='fesl(%s,%s)' % (p, s))

        objExpr = bDelta.sum('*', '*') + alphaf * fDelta.sum("*", "*")
    else:
        objExpr = bDelta.sum('*', '*')
        
    m.setObjective(objExpr, GRB.MINIMIZE)
    m.write('update_proportion.lp')
    m.optimize()

    Fresult = np.zeros([tumors, cells], dtype=np.float)
    for p in range(tumors):
        for k in range(cells):
            Fresult[p, k] = F[p, k].getAttr(GRB.Attr.X)

    objVal = m.objVal
    m.reset()
    return Fresult, objVal


def getDistanceMatrix(Ctotal):
    cellsTotal = Ctotal.shape[0]
    distance = np.zeros((cellsTotal, cellsTotal))

    for i in range(cellsTotal):
        for j in range(cellsTotal):
            distance[i, j] = np.sum(abs(Ctotal[i, :] - Ctotal[j, :]))

    return distance


'''
Update topological tree structure of phylogenetic tree.
'''
def updateTree(Ctotal, cells=None, alpha1=0.1, root=0):
    cellsTotal = Ctotal.shape[0]

    if (cells == None):
        cells = cellsTotal

    m = Model('phylo')
    m.setParam("OutputFlag", 0)

    g = m.addVars(cellsTotal, cellsTotal, cellsTotal, vtype='B', name='g')
    S = m.addVars(cellsTotal, cellsTotal, vtype='B', name='S')

    for t in range(cellsTotal):
        for u in range(cellsTotal):
            if u != t and u != root:
                expr = 0
                for v in range(cellsTotal):
                    expr += (g[u, v, t] - g[v, u, t])
                m.addConstr(expr == 0, name='conserv(%s,%s)' % (u, t))

    for t in range(cellsTotal):
        if t != root:
            m.addConstr(g.sum('*', t, t) == 1, name='sink(%s)' % (t, ))
    m.addConstrs(
        (g.sum(t, '*', t) == 0 for t in range(cellsTotal)), name='leaf')
    m.addConstr(g.sum('*', root, '*') == 0, name='root')
    for t in range(cellsTotal):
        if t != root:
            m.addConstr(g.sum(root, '*', t) == 1, name='source(%s)' % (t, ))

    for t in range(cellsTotal):
        for u in range(cellsTotal):
            for v in range(cellsTotal):
                m.addConstr(
                    S[u, v] - g[u, v, t] >= 0,
                    name='present(%s,%s,%s)' % (u, v, t))

    wDelta = getDistanceMatrix(Ctotal)

    objExpr = S.sum('*', '*')
    for u in range(cellsTotal):
        for v in range(cellsTotal):
            objExpr += wDelta[u, v] * S[u, v]
    objExpr = alpha1 * objExpr

    m.setObjective(objExpr, GRB.MINIMIZE)
    m.optimize()

    if (m.status != 2):
        return None

    Sresult = np.zeros([cellsTotal, cellsTotal], dtype=np.int8)
    for u in range(cellsTotal):
        for v in range(cellsTotal):
            Sresult[u, v] = S[u, v].getAttr(GRB.Attr.X)

    objVal = m.objVal
    m.reset()
    return Sresult, objVal


def calcObjVal(B, F, C, metric='L1'):
    diff = B - np.matmul(F, C)
    if metric == 'L1':
        return(np.sum(abs(diff)))
    elif metric == 'L2':
        return(np.sum(np.multiply(diff, diff)))


'''
Update CNV in single-cell matrix.
'''
def updateCopyNum(B, F, S1, S2, CRefer, pr, cells, FISHRefer, FISHPositions, mask, alpha1=0.1, alpha2=0.1, vType='I', stopVal=None, Cap=10):
    tumors = B.shape[0]
    # CRefer = np.dot(ploidy, CRefer)
    cellsTotal = cells + CRefer.shape[0]
    fishTotal = cells +  FISHRefer.shape[0]
    copyNum = B.shape[1]
    maskedCopy = len(mask)
    m = Model('phylo')
    m.setParam("OutputFlag", 0)
    FISHnum, FISHlociNum = FISHRefer.shape
    if stopVal != None:
        m.setParam("BestObjStop", stopVal)
    C = m.addVars(cells, copyNum, name='C', vtype=vType)

    bDelta = m.addVars(tumors, copyNum, name='bDelta', vtype='C')
    for p in range(tumors):
        for s in range(copyNum):
            expr = B[p, s]
            for k in range(cells):
                expr += -(F[p, k] * pr[k, k]/2 * C[k, s])
            m.addConstr(
                bDelta[p, s] - (expr) >= 0, name='resu(%s,%s)' % (p, s))
            m.addConstr(
                bDelta[p, s] + (expr) >= 0, name='resl(%s,%s)' % (p, s))
    m.addConstrs(
        (C[k, s] >= 0 for k in range(cells) for s in range(copyNum)),
        name='pos')
        
    '''need turn this off if work with ploidy'''
    if Cap is not None:
        m.addConstrs(
            (C[k, s] <= Cap for k in range(cells) for s in range(copyNum)),
            name='cap')

    objExpr = bDelta.sum('*', '*')

    '''tree on cell'''
    if alpha1 > 0:      
        wDelta = m.addVars(
            cellsTotal, cellsTotal, copyNum, vtype='C', name='wDelta')

        for u in range(cellsTotal):
            for v in range(cellsTotal):
                for i in range(maskedCopy):
                    if (u < cells):
                        Cui = C[u, mask[i]]
                    else:
                        Cui = CRefer[u-cells, mask[i]]
                    if (v < cells):
                        Cvi = C[v, mask[i]]
                    else:
                        Cvi = CRefer[v - cells, mask[i]]
                    m.addConstr(
                        wDelta[u, v, mask[i]] - Cui + Cvi >= 0,
                        name='delu(%s,%s,%s)' % (u, v, i))
                    m.addConstr(
                        wDelta[u, v, mask[i]] - Cvi + Cui >= 0,
                        name='dell(%s,%s,%s)' % (u, v, i))
                objExpr += alpha1 * wDelta.sum(u, v, '*') * S1[u, v]
    '''
    constraints in FISH probes
    '''
    if alpha2 > 0:
        hDelta = m.addVars(FISHnum, FISHlociNum, name='hl', vtype='C')
        for i in range(FISHnum):
            for j in range(FISHlociNum):
                expr = FISHRefer[i, j] - C[i, FISHPositions[j]] * pr[i, i] / 2

                m.addConstr(hDelta[i, j] - (expr) >= 0, name='hesu(%s,%s)' % (i, j))
                m.addConstr(hDelta[i, j] + (expr) >= 0, name='hesl(%s,%s)' % (i, j))
                
        objExpr += alpha2 * hDelta.sum('*', '*')


    m.setObjective(objExpr, GRB.MINIMIZE)
    m.optimize()

    if (m.status != 2 and m.status != 15):
        return None

    Cresult = np.zeros([cells, copyNum], dtype=np.float)
    for k in range(cells):
        for s in range(copyNum):
            Cresult[k, s] = C[k, s].getAttr(GRB.Attr.X)

    objVal = m.objVal
    m.reset()
    return Cresult, objVal

def updatePloidy(B, C, F, H, FISHInfer, pr, alpha2, vType='C'):
    '''
    decompose the FISH to get the ploidy
    a_2 * ||H-P'HInfer||
    '''
    m = Model('phylo')
    m.setParam("OutputFlag", 0)  

    tumors = B.shape[0]
    cells= C.shape[0]
    copyNum = C.shape[1]
    FISHnum, FISHlociNum = H.shape
    P = m.addVars(FISHnum, FISHnum, name='p', vtype='C')
    bDelta = m.addVars(tumors, copyNum, name='bDelta', vtype='C')
    for i in range(tumors):
        for j in range(copyNum):
            expr = B[i, j]
            for k in range(cells):
                expr += -(F[i, k] * P[k, k]/2 * C[k, j])
            m.addConstr(
                bDelta[i, j] - (expr) >= 0, name='resu(%s,%s)' % (i, j))
            m.addConstr(
                bDelta[i, j] + (expr) >= 0, name='resl(%s,%s)' % (i, j))
    objExpr = bDelta.sum('*', '*')

    '''the ploidy value: 0 <= P <= 8 '''
    m.addConstrs(
        (P[k, s] >= 0 for k in range(FISHnum) for s in range(FISHnum)),
        name='pos1')
    
    m.addConstrs((P[k, s] <= 8.0 for k in range(FISHnum) for s in range(FISHnum)),
        name='lowP')
    
    '''make all the element not diagnal to be 0'''
    for i in range(FISHnum):
        for j in range(FISHnum):
            if i != j:
                m.addConstr(P[i,j] == 0, name='sparseP(%s,%s)' % (i, j))
    
    if alpha2 > 0:
        hDelta = m.addVars(FISHnum, FISHlociNum, name='hl', vtype=vType)
        for i in range(FISHnum):
            for j in range(FISHlociNum):
                expr = H[i, j]
                for k in range(FISHnum):
                    expr += -(P[i, k]*FISHInfer[k, j])/2
                m.addConstr(
                    hDelta[i, j] - (expr) >= 0, name='hesu(%s,%s)' % (i, j))
                m.addConstr(
                    hDelta[i, j] + (expr) >= 0, name='hesl(%s,%s)' % (i, j))

        objExpr += alpha2 * hDelta.sum('*', '*')

        
    m.setObjective(objExpr, GRB.MINIMIZE)
    m.write('update_proportion.lp')
    m.optimize()

    Presult = np.zeros((FISHnum, FISHnum), dtype=np.float)
    for i in range(FISHnum):
        for j in range(FISHnum):
            Presult[i, j] = P[i, j].getAttr(GRB.Attr.X)

    objVal = m.objVal
    m.reset()
    return Presult, objVal





def makeSurePath(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    return
