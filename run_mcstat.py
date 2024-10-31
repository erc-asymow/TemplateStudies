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

def run_all():
    options = []
    # default
    options.append( './mcstat --ntoys=10000 --ntoysFC=1000 --nbins=200 --nevents=2000000  --lumiscale=1 --asym=0.015  --doFC --decorrelate             --tag=1_200_0p015_decorr' )
    options.append( './mcstat --ntoys=10000 --ntoysFC=1000 --nbins=200 --nevents=2000000  --lumiscale=1 --asym=0.015  --doFC --decorrelate --doFCcheat --tag=1_200_0p015_decorr_FCCheat' )
    options.append( './mcstat --ntoys=10000 --ntoysFC=1000 --nbins=200 --nevents=2000000  --lumiscale=1 --asym=0.015  --doFC --decorrelate --doBarlett --tag=1_200_0p015_decorr_Barlett' )
    # asym = 0.03
    options.append( './mcstat --ntoys=10000 --ntoysFC=1000 --nbins=200 --nevents=2000000  --lumiscale=1 --asym=0.03   --doFC --decorrelate             --tag=1_200_0p03_decorr' )
    options.append( './mcstat --ntoys=10000 --ntoysFC=1000 --nbins=200 --nevents=2000000  --lumiscale=1 --asym=0.03   --doFC --decorrelate --doFCcheat --tag=1_200_0p03_decorr_FCCheat' )
    options.append( './mcstat --ntoys=10000 --ntoysFC=1000 --nbins=200 --nevents=2000000  --lumiscale=1 --asym=0.03   --doFC --decorrelate --doBarlett --tag=1_200_0p03_decorr_Barlett' )
    # asym = 0.05
    options.append( './mcstat --ntoys=10000 --ntoysFC=1000 --nbins=200 --nevents=2000000  --lumiscale=1 --asym=0.06   --doFC --decorrelate             --tag=1_200_0p06_decorr' )
    options.append( './mcstat --ntoys=10000 --ntoysFC=1000 --nbins=200 --nevents=2000000  --lumiscale=1 --asym=0.06   --doFC --decorrelate --doFCcheat --tag=1_200_0p06_decorr_FCCheat' )
    options.append( './mcstat --ntoys=10000 --ntoysFC=1000 --nbins=200 --nevents=2000000  --lumiscale=1 --asym=0.06   --doFC --decorrelate --doBarlett --tag=1_200_0p06_decorr_Barlett' )
    # nbins = 100
    options.append( './mcstat --ntoys=10000 --ntoysFC=1000 --nbins=100 --nevents=2000000  --lumiscale=1 --asym=0.015  --doFC --decorrelate             --tag=1_100_0p015_decorr' )
    options.append( './mcstat --ntoys=10000 --ntoysFC=1000 --nbins=100 --nevents=2000000  --lumiscale=1 --asym=0.015  --doFC --decorrelate --doFCcheat --tag=1_100_0p015_decorr_FCCheat' )
    options.append( './mcstat --ntoys=10000 --ntoysFC=1000 --nbins=100 --nevents=2000000  --lumiscale=1 --asym=0.015  --doFC --decorrelate --doBarlett --tag=1_100_0p015_decorr_Barlett' )
    # nbins = 20
    options.append( './mcstat --ntoys=10000 --ntoysFC=1000 --nbins=20  --nevents=2000000  --lumiscale=1 --asym=0.015  --doFC --decorrelate             --tag=1_20_0p015_decorr' )
    options.append( './mcstat --ntoys=10000 --ntoysFC=1000 --nbins=20  --nevents=2000000  --lumiscale=1 --asym=0.015  --doFC --decorrelate --doFCcheat --tag=1_20_0p015_decorr_FCCheat' )
    options.append( './mcstat --ntoys=10000 --ntoysFC=1000 --nbins=20  --nevents=2000000  --lumiscale=1 --asym=0.015  --doFC --decorrelate --doBarlett --tag=1_20_0p015_decorr_Barlett' )
    # lumiscale = 10
    options.append( './mcstat --ntoys=10000 --ntoysFC=1000 --nbins=200 --nevents=2000000  --lumiscale=10 --asym=0.015  --doFC --decorrelate             --tag=10_200_0p015_decorr' )
    options.append( './mcstat --ntoys=10000 --ntoysFC=1000 --nbins=200 --nevents=2000000  --lumiscale=10 --asym=0.015  --doFC --decorrelate --doFCcheat --tag=10_200_0p015_decorr_FCCheat' )
    options.append( './mcstat --ntoys=10000 --ntoysFC=1000 --nbins=200 --nevents=2000000  --lumiscale=10 --asym=0.015  --doFC --decorrelate --doBarlett --tag=10_200_0p015_decorr_Barlett' )
    # default correlated
    options.append( './mcstat --ntoys=10000 --ntoysFC=1000 --nbins=200 --nevents=2000000  --lumiscale=1 --asym=0.015  --doFC             --tag=1_200_0p015_corr' )
    options.append( './mcstat --ntoys=10000 --ntoysFC=1000 --nbins=200 --nevents=2000000  --lumiscale=1 --asym=0.015  --doFC --doFCcheat --tag=1_200_0p015_corr_FCCheat' )
    options.append( './mcstat --ntoys=10000 --ntoysFC=1000 --nbins=200 --nevents=2000000  --lumiscale=1 --asym=0.015  --doFC --doBarlett --tag=1_200_0p015_corr_Barlett' )
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
    run_all()
