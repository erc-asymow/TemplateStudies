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

args = parser.parse_args()

if __name__ == '__main__':

    cmd_histo_iter0 = './massscales --lumi=-1 --firstIter=0 --lastIter=2 --skipReco '+\
        ' --tag='+args.tag+' '+\
        ' --run=Iter0 '+\
        ' --nRMSforGausFit=-1.0 '+\
        ' --biasResolution=-1.0 ' +\
        ' --minNumEvents=100 --minNumEventsPerBin=30 '+\
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
    
    cmd_histo_iter1 = cmd_histo_iter0.replace('--run=Iter0', '--run=Iter1')
    cmd_histo_iter1 += ' --usePrevFit '+\
        ' --tagPrevFit='+args.tag+' '+\
        ' --runPrevFit=Iter0 '
    print(cmd_histo_iter1)
    if not args.dryrun:
        os.system(cmd_histo_iter1)
        #print()
    cmd_fit_iter1 = cmd_fit_iter0.replace('--run=Iter0', '--run=Iter1')
    print(cmd_fit_iter1)
    if not args.dryrun:
        os.system(cmd_fit_iter1)
        #print()

    cmd_histo_iter2 = cmd_histo_iter0.replace('--run=Iter0', '--run=Iter2')
    cmd_histo_iter2 += ' --usePrevFit '+\
        ' --tagPrevFit='+args.tag+' '+\
        ' --runPrevFit=Iter1 '
    print(cmd_histo_iter2)
    if not args.dryrun:
        os.system(cmd_histo_iter2)
    cmd_fit_iter2 = cmd_fit_iter0.replace('--run=Iter0', '--run=Iter2')
    print(cmd_fit_iter2)
    if not args.dryrun:
        os.system(cmd_fit_iter2)

    
