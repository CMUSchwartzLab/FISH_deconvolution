import argparse

from glob import glob
import re
import os.path
import numpy as np
import testFunction
import time


from deconvsolver import Simulation, MarkerPositions
import postprocess as pp


CLONE_DECAY = .8
MAX_COPY = 10

PROBLEM_INDEX = re.compile('simulate_data_(.*?).pkl')


def yaml_dump_list(lst, outfile):
    ncol = len(lst)
    print(' - [ ', file=outfile, end='')
    for j in range(ncol):
        if j != 0:
            print(', ', file=outfile, end='')
        print(lst[j], file=outfile, end='')
    print(' ]', file=outfile)

def yaml_dump_matrix(mat, outfile):
    mrows, ncol = mat.shape
    for i in range(mrows):
        print(' - [ ', file=outfile, end='')
        for j in range(ncol):
            if j != 0:
                print(', ', file=outfile, end='')
            print(mat[i, j], file=outfile, end='')
        print(' ]', file=outfile)


def find_simulated_datafiles(args):
    simulated_datafiles = list(
        sorted(glob('%s/*.pkl' % (args.simulation_directory,))))
    start = args.offset
    return simulated_datafiles[start:start + args.limit]


def print_problem_banner(index, of, named):
    print("=====")
    print("=")
    print("= Running %s of %s (%s)" % (index + 1, of, named))
    print("=")
    print("=====")


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Do a trial run')
    parser.add_argument('interval_file')
    parser.add_argument('simulation_directory')
    parser.add_argument('input_date', type=str)
    parser.add_argument('output_date', type=str)
    parser.add_argument('--alphaf', type=float, default=1.0)
    parser.add_argument('--alpha1', type=float, default=1.0)
    parser.add_argument('--alpha2', type=float, default=1.0)
    parser.add_argument('--noise', type=float, default=0.0)
    parser.add_argument('--offset', type=int, default=0)
    parser.add_argument('--limit', type=int, default=1000000)
    parser.add_argument('--solver', type=str, default='gurobi')
    parser.add_argument('--referenceNum', type=int, default=0)
    parser.add_argument('--mask', type=float, default=1.0)
    parser.add_argument('--clone-decay', type=float, default=CLONE_DECAY)
    parser.add_argument('--prefix', default='origProb')
    return parser.parse_args()

def getFISHInfer(C, index):
    k, _ = C.shape
    m = len(index)
    FISH = np.zeros((k, m))
    for i in range(m):
        FISH[:,i] = C[:, index[i]]
    return FISH

def run_one_sample(sample, prob_id, intervals, args):
    simulation = Simulation.read(sample)

    '''
    Load the necessary data from the simulated data
    '''
    # here, bulk data is simulate
    bulk_data = simulation.bulk_data
    clone_freq = simulation.reference_frequencies
    true_cells, true_freqs = simulation.true_cells, simulation.true_frequencies
    true_ploidy = simulation.true_cells_ploidy
    expandReferFISH = simulation.expand_uFISH
    referFISHploidy = simulation.reference_fish_ploidy
    
    pr = np.diag(referFISHploidy)
    corPositons = [item for sublist in simulation.correlate_pos for item in sublist] 
    '''only get proportion of genomic for reference'''
    mask = np.random.choice(range(bulk_data.shape[1]), int(bulk_data.shape[1]*args.mask), replace=False)
    '''select reference cell'''
    sample_cells = simulation.reference_cells
    selected = np.random.choice(range(sample_cells.shape[0]), size=args.referenceNum, replace=False)
    reference_cell = sample_cells[selected, :]
    #sample_cells[:] = 2    # assume there is no real SCS data
    '''rescale reference FISH'''
    expandReferFISH[expandReferFISH>40] = 40

    ''''initializa cell tree'''
    k = sample_cells.shape[0]
    ctotal = np.zeros((2*k, sample_cells.shape[1]), dtype=np.float)
    #1 use normalized cell
    ctotal[k:k*2, :] = sample_cells
    #2 use unnomralized cell
    #ctotal[k:k*2, :] = np.dot(np.diag(referPloidy).T, sample_cells) / 2

    ctotal = np.concatenate((ctotal, 2*np.ones((1, sample_cells.shape[1]))), axis=0)
    
    oldObj = [0, 0, 0, 0]
    thresholdI = 10 ** (-4)
    step = 1
    start_time = time.time()
    if args.solver == 'gurobi':
        import ILP_solver_v7 as gs
    elif args.solver == 'scip':
        import ILP_solver_SCIP as gs
    else:
        exit("Only Gurobi and SCIP available")

    while True:
        F, objVal1 = gs.updateProportion(bulk_data, ctotal, pr, clone_freq, args.alphaf, k, root=2*k)
        S1, objVal2 = gs.updateTree(ctotal, k, args.alpha1, root=2*k)
        step += 1
        C, objVal4 = gs.updateCopyNum(bulk_data, F, S1, None, sample_cells, pr, k, expandReferFISH, corPositons, mask=mask, alpha1=args.alpha1, alpha2=args.alpha2, Cap=10)
        
        FISHInfer = getFISHInfer(C, corPositons) #normalized C to get normalized FISH
        #FISHInfer = np.dot(pr/2, FISHInfer)
        H = expandReferFISH
        pr, objVal5 = gs.updatePloidy(bulk_data, C, F, H, FISHInfer, np.diag(referFISHploidy), args.alpha2)

        #1 use normalize cell
        ctotal[0:k, :] = C
        #2 use unnormalized cell
        #ctotal[0:k, :] = np.dot(pr/2, C)

        change = abs(oldObj[3] - objVal4)
        oldObj[3] = objVal4
        #oldObj[2] = objVal3
        oldObj[1] = objVal2
        oldObj[0] = objVal1
        if (change < thresholdI or step > 100):
            break

    print('Finished in', step, 'iterations.')
    print('Total time is %s minutes' % (round((time.time()-start_time)/60, 2)))
    print("The parameter for alphaf, alpha1, alpha2 are ",
        "%s, %s, %s" % (
            args.alphaf, args.alpha1, args.alpha2))
    
    post_fix = '_'.join(
        map(str,
            [prob_id, args.alphaf, args.alpha1, args.alpha2]))
    '''
    save files in txt format
    '''
    pr_result = []
    for i in range(pr.shape[0]):
        pr_result.append(pr[i, i])

    basename = '../results/' + args.output_date + '_%s' % args.mask + '/sample_output_' + post_fix
    with open(basename + '.txt', 'wt') as outfile:
        print('---', file=outfile)
        yaml_dump_matrix(C.T, outfile)
        print('---', file=outfile)
        yaml_dump_matrix(F.T, outfile)
        print('---', file=outfile)
        yaml_dump_list(pr_result, outfile)
        if args.alpha1 > 0:
            print('---', file=outfile)
            yaml_dump_matrix(S1, outfile)

    '''
    calculate the necessary correctiness
    '''
    #normalized_ploidy = [2, 2, 2, 2, 2, 2]
    pp.calculations(true_cells, C.T, true_freqs, F.T, true_ploidy, pr_result, args.output_date + '_%s' % args.mask, post_fix)


def main():
    args = parse_args()
    intervals = MarkerPositions.read(args.interval_file)
    simulated_datafiles = find_simulated_datafiles(args)
    testFunction.CheckDirectory('../results/' + args.output_date + '_%s' % args.mask)

    for i, sample in enumerate(simulated_datafiles):
        prob_id = os.path.basename(sample)
        m = PROBLEM_INDEX.match(prob_id)
        if m:
            prob_id = m.group(1)

        print_problem_banner(index=i, of=len(
            simulated_datafiles), named=sample)
        run_one_sample(sample, prob_id, intervals, args)

if __name__ == '__main__':
    main()
