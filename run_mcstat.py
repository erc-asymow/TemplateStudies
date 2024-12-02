import argparse
import os
import sys
import copy
import math

from multiprocessing import Process
   
parser = argparse.ArgumentParser(description='run')

parser.add_argument('--none', action='store_true'  , help = 'none')
parser.add_argument('--dryrun', action='store_true'  , help = 'dry run')
parser.add_argument('--tag', default=''  , help = '')
parser.add_argument('--val', dest = 'val'  , type = float,  default=0.0, help='')

args = parser.parse_args()

def run_one_opt(options, iopt):
    if args.dryrun:
        print(options[iopt])
    else:
        os.system(options[iopt])

def run_all( do_m5s=False ):

    # default
    options_base = []
    base_str = './mcstat --ntoys=10000 --ntoysFC=4000 --nbins=200 --nevents=2000000  --lumiscale=1  --asym=0.015  --doFC --decorrelate'
    #base_str = './mcstat --ntoys=10000 --ntoysFC=10 --nbins=200 --nevents=2000000  --lumiscale=1  --asym=0.015  --doFC --decorrelate'

    posttag = ''
    #posttag += 'adhoc'
    if do_m5s:
        posttag += 'm5s'
        base_str += ' --nsigmas=-5'
        
    options_base.append( base_str + '                                     --tag=1_200_0p015_decorr'+posttag )
    #options_base.append( base_str + ' --doFCcheat                         --tag=1_200_0p015_decorr'+posttag+'_FCCheat' )
    #options_base.append( base_str + ' --doBarlett                         --tag=1_200_0p015_decorr'+posttag+'_Barlett' )
    #options_base.append( base_str + ' --doPoisson                         --tag=1_200_0p015_decorr'+posttag+'_Poisson' )
    #options_base.append( base_str + ' --computeJtildeError                --tag=1_200_0p015_decorr'+posttag+'_Jtilde' )
    #options_base.append( base_str + ' --FCfixToTrue                       --tag=1_200_0p015_decorr'+posttag+'_FCfixToTrue' )
    #options_base.append( base_str + ' --doBarlett --FCfixToTrue           --tag=1_200_0p015_decorr'+posttag+'_BarlettFCfixToTrue' )
    #options_base.append( base_str + ' --computeJtildeError  --FCfixToTrue --tag=1_200_0p015_decorr'+posttag+'_JtildeFCfixToTrue' )

    options = []
    for iop in options_base:
        options.append( iop )
        opt = iop.replace( '--decorrelate', '             ' ).replace( '_decorr', '_corr' )
        options.append( opt )
        opt = iop.replace( '--asym=0.015', '--asym=0.030' ).replace( '_0p015_', '_0p03_' )
        options.append( opt )
        opt = iop.replace( '--asym=0.015', '--asym=0.060' ).replace( '_0p015_', '_0p06_' )
        options.append( opt )
        opt = iop.replace( '--nbins=200', '--nbins=100' ).replace( '1_200', '1_100' )
        options.append( opt )
        opt = iop.replace( '--nbins=200', '--nbins=20 ' ).replace( '1_200', '1_20' )
        options.append( opt )
        opt = iop.replace( '--lumiscale=1 ', '--lumiscale=10' ).replace( '1_200', '10_200' )
        options.append( opt )
        opt = iop.replace( '--lumiscale=1 ', '--lumiscale=40' ).replace( '1_200', '40_200' )
        options.append( opt )
    
    ps = []        
    counter = 0
    for i,ip in enumerate(options):
        p = Process(target=run_one_opt, args=(options, i))
        p.start()
        ps.append(p)
        counter += 1
    for p in ps: 
        p.join()
    return

if __name__ == '__main__':
    run_all( do_m5s=False )
