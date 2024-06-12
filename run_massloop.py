import argparse
import os
import sys
import copy
import math
import ROOT

   
parser = argparse.ArgumentParser(description='run')

parser.add_argument('--none', action='store_true'  , help = 'none')
parser.add_argument('--dryrun', action='store_true'  , help = 'dry run')
parser.add_argument('--tag',   default='SmearRealistic' , help = 'algo')
parser.add_argument('--niter', dest = 'niter'  , type = int,  default=1, help='')

args = parser.parse_args()

if __name__ == '__main__':

    cmd_histo_iter0 = './massscales --lumi=-1 --firstIter=0 --lastIter=2 --skipReco '+\
        ' --tag='+args.tag+' '+\
        ' --run=Iter0 '+\
        ' --nRMSforGausFit=-1.0 '+\
        ' --biasResolution=0.1 ' +\
        ' --minNumEvents=100 --minNumEventsPerBin=30 '+\
        ' --minNumMassBins=4 '+\
        ' --rebin=2 '+\
        ' --fitNorm --fitWidth '    
    print(cmd_histo_iter0)
    if not args.dryrun:
        os.system(cmd_histo_iter0)
        #print()
    cmd_fit_iter0 = './massfit --nevents=1 --bias=-1 '+\
        '--tag='+args.tag+' '+\
        '--run=Iter0 '
    print(cmd_fit_iter0)
    if not args.dryrun:
        os.system(cmd_fit_iter0)
        #print()
    cmd_resol_iter0 = './resolfit --nevents=1 --bias=-1 '+\
        ' --tag='+args.tag+' '+\
        ' --run=Iter0 '+\
        ' --maxSigmaErr=0.1 '
    print(cmd_resol_iter0)
    if not args.dryrun:
        os.system(cmd_resol_iter0)
        #print()

    for iter in range(1, args.niter+1):
        cmd_histo_iteri = cmd_histo_iter0.replace('--run=Iter0', '--run=Iter'+str(iter))
        cmd_histo_iteri += ' --usePrevFit '+\
            ' --tagPrevFit='+args.tag+' '+\
            ' --runPrevFit=Iter'+str(iter-1)+' '
        cmd_histo_iteri += ' --useSmearFit '+\
            ' --tagSmearFit='+args.tag+' '+\
            ' --runSmearFit=Iter'+str(iter-1)+' '
        print(cmd_histo_iteri)
        if not args.dryrun:
            os.system(cmd_histo_iteri)
            #print()
        cmd_fit_iteri = cmd_fit_iter0.replace('--run=Iter0', '--run=Iter'+str(iter))
        print(cmd_fit_iteri)
        if not args.dryrun:
            os.system(cmd_fit_iteri)
            #print()
        cmd_resol_iteri = cmd_resol_iter0.replace('--run=Iter0', '--run=Iter'+str(iter))
        print(cmd_resol_iteri)
        if not args.dryrun:
            os.system(cmd_resol_iteri)
            #print()

