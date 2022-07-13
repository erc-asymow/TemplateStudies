import argparse
import os
import sys
import copy

parser = argparse.ArgumentParser(description='run')

parser.add_argument('--none', action='store_true'  , help = 'none')
parser.add_argument('--tag',  default='test'       , help = 'tag')
parser.add_argument('--algo', default='all'        , help = 'algo')
parser.add_argument('--nevents', dest = 'nevents'  , type = int,  default=10000000, help='number of events')
parser.add_argument('--rebinX', dest = 'rebinX'  , type = int,  default=-1, help='')
parser.add_argument('--rebinY', dest = 'rebinY'  , type = int,  default=-1, help='')

args = parser.parse_args()

pol_default = {
    'corr'   : [10,4],
    'A0'     : [3,3],
    'A1'     : [3,3],
    'A2'     : [3,3],
    'A3'     : [3,3],
    'A4'     : [3,3]
}

if args.algo=='jac2':
    command  = './jac2 --nevents='+str(args.nevents) +' --tag='+args.tag+' --run=closure'
    for k in pol_default.keys():
        command += ' --degs_'+k+'_x='+str(pol_default[k][0])+' --degs_'+k+'_y='+str(pol_default[k][1])
    command += ' --toyTF2_corr'
    print command
    os.system(command)   

elif args.algo=='fit':
    command  = './fit --nevents='+str(args.nevents) +' --tag='+args.tag+' --run=closure'
    for k in pol_default.keys():
        command += ' --degs_'+k+'_x='+str(pol_default[k][0])+' --degs_'+k+'_y='+str(pol_default[k][1])
    if (args.rebinX>0 or args.rebinY>0):
        command += ' --rebinX='+str(args.rebinX)+' --rebinY='+str(args.rebinY)
    fit_opts = ['jUL', 'j0', 'j1', 'j2', 'j3', 'j4']
    for i in range(1, len(fit_opts)+1):
        new_command = copy.deepcopy(command)
        for j in range(i): 
            new_command += (' --'+fit_opts[j])
        print new_command
        os.system(new_command) 

elif args.algo=='fit_fast':
    command  = './fit --nevents='+str(args.nevents) +' --tag='+args.tag+' --run=closure'
    for k in pol_default.keys():
        command += ' --degs_'+k+'_x='+str(pol_default[k][0])+' --degs_'+k+'_y='+str(pol_default[k][1])
    if (args.rebinX>0 or args.rebinY>0):
        command += ' --rebinX='+str(args.rebinX)+' --rebinY='+str(args.rebinY)
    #fit_opts = ['jUL', 'j0', 'j1', 'j2', 'j3', 'j4']
    fit_opts = ['jUL', 'j0', 'j1']
    for i in range(len(fit_opts)):
        command += (' --'+fit_opts[i])
    print command
    os.system(command) 
    
