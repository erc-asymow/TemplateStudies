import argparse
import os
import sys
import copy

parser = argparse.ArgumentParser(description='run')

parser.add_argument('--none', action='store_true'  , help = 'none')
parser.add_argument('--dryrun', action='store_true'  , help = 'dry run')
parser.add_argument('--tag',  default='test'       , help = 'tag')
parser.add_argument('--post_tag',  default='test'  , help = 'post tag')
parser.add_argument('--algo',   default='all'        , help = 'algo')
parser.add_argument('--run',   default='full'        , help = 'run')
parser.add_argument('--nevents', dest = 'nevents'  , type = int,  default=10000000, help='number of events')
parser.add_argument('--rebinX', dest = 'rebinX'  , type = int,  default=-1, help='')
parser.add_argument('--rebinY', dest = 'rebinY'  , type = int,  default=-1, help='')
parser.add_argument('--scale0', action='store_true'  , help = 'scale A0(x,y)')
parser.add_argument('--scale1', action='store_true'  , help = 'scale A1(x,y)')
parser.add_argument('--scale2', action='store_true'  , help = 'scale A2(x,y)')
parser.add_argument('--scale3', action='store_true'  , help = 'scale A3(x,y)')
parser.add_argument('--scale4', action='store_true'  , help = 'scale A4(x,y)')
parser.add_argument('--jacmass', dest = 'jacmass'  , type = int,  default=-1, help='')

args = parser.parse_args()

pol_default = {
    'run'    : "full",
    'corr'   : [10,4],
    'A0'     : [3,2],
    'A1'     : [3,3],
    'A2'     : [3,2],
    'A3'     : [3,4],
    'A4'     : [3,3]
}

pol_systs = []

pol_systs.append( pol_default )

pol_syst = copy.deepcopy(pol_default)
pol_syst['A0'] = [2,2]
pol_syst['A2'] = [2,2]
pol_systs.append( pol_syst )

pol_syst= copy.deepcopy(pol_default)
pol_syst['A3'] = [2,4]
pol_systs.append( pol_syst )

pol_syst= copy.deepcopy(pol_default)
pol_syst['A3'] = [3,2]
pol_systs.append( pol_syst )

pol_syst= copy.deepcopy(pol_default)
pol_syst['A4'] = [1,3]
pol_systs.append( pol_syst )

pol_syst = copy.deepcopy(pol_default)
pol_syst['run'] = "corr"
pol_syst['corr'] = [3,2]
pol_syst['A0']   = [2,2]
pol_syst['A1']   = [2,2]
pol_syst['A2']   = [2,2]
pol_syst['A2']   = [2,2]
pol_syst['A3']   = [2,2]
pol_syst['A4']   = [2,2]
pol_systs.append( pol_syst )

pol_syst = copy.deepcopy(pol_default)
pol_syst['run'] = "corr"
pol_syst['corr'] = [5,2]
pol_syst['A0']   = [2,2]
pol_syst['A1']   = [2,2]
pol_syst['A2']   = [2,2]
pol_syst['A2']   = [2,2]
pol_syst['A3']   = [2,2]
pol_syst['A4']   = [2,2]
pol_systs.append( pol_syst )

pol_syst = copy.deepcopy(pol_default)
pol_syst['run'] = "corr"
pol_syst['corr'] = [3,2]
pol_syst['A0']   = [2,2]
pol_syst['A1']   = [2,2]
pol_syst['A2']   = [2,2]
pol_syst['A2']   = [2,2]
pol_syst['A3']   = [1,2]
pol_syst['A4']   = [2,2]
pol_systs.append( pol_syst )

pol_syst = copy.deepcopy(pol_default)
pol_syst['run'] = "corr"
pol_syst['corr'] = [3,2]
pol_syst['A0']   = [2,2]
pol_syst['A1']   = [2,2]
pol_syst['A2']   = [2,2]
pol_syst['A2']   = [2,2]
pol_syst['A3']   = [2,2]
pol_syst['A4']   = [1,2]
pol_systs.append( pol_syst )

pol_syst = copy.deepcopy(pol_default)
pol_syst['run'] = "corr"
pol_syst['corr'] = [3,2]
pol_syst['A0']   = [1,2]
pol_syst['A1']   = [1,2]
pol_syst['A2']   = [1,2]
pol_syst['A2']   = [1,2]
pol_syst['A3']   = [1,2]
pol_syst['A4']   = [1,2]
pol_systs.append( pol_syst )

pol_syst = copy.deepcopy(pol_default)
pol_syst['run'] = "grid"
pol_syst['corr'] = [8,6]
pol_syst['A0']   = [8,6]
pol_syst['A1']   = [8,6]
pol_syst['A2']   = [8,6]
pol_syst['A2']   = [8,6]
pol_syst['A3']   = [8,6]
pol_syst['A4']   = [8,6]
pol_systs.append( pol_syst )

if args.algo=='jac2':
    command  = './jac2 --nevents='+str(args.nevents) +' --tag='+args.tag+' --run='+pol_default['run']
    for k in pol_default.keys():
        if k=='run': continue
        command += ' --degs_'+k+'_x='+str(pol_default[k][0])+' --degs_'+k+'_y='+str(pol_default[k][1])
    print(command)
    os.system(command)   

elif args.algo=='jac2_systs':
    count = 0
    for syst in pol_systs:
        command  = './jac2 --nevents='+str(args.nevents) +' --tag='+args.tag
        for k in syst.keys():
            if k=='run':
                command += ' --run='+syst[k]
            else:
                command += ' --degs_'+k+'_x='+str(syst[k][0])+' --degs_'+k+'_y='+str(syst[k][1])
        print(command)
        if not args.dryrun:
            os.system(command)
        count += 1
    print('%.0f jobs will be submitted' % count)

elif args.algo=='jac2_systs_vsN':
    count = 0
    for syst in pol_systs:
        for n in [10000000,
                  100000000,
                  1000000000,
                  4000000000,
                  10000000000,
                  ]:            
            command  = './jac2 --nevents='+str(n) +' --tag='+args.tag
            for k in syst.keys():
                if k=='run':
                    command += ' --run='+syst[k]
                else:
                    command += ' --degs_'+k+'_x='+str(syst[k][0])+' --degs_'+k+'_y='+str(syst[k][1])
            print(command)
            if not args.dryrun:
                os.system(command)   
            count += 1
            
    print('%.0f jobs will be submitted' % count)
            
elif args.algo=='fit':
    command  = './fit --nevents='+str(args.nevents) +' --tag='+args.tag+' --run=closure --post_tag='+args.post_tag
    for k in pol_default.keys():
        command += ' --degs_'+k+'_x='+str(pol_default[k][0])+' --degs_'+k+'_y='+str(pol_default[k][1])
    if (args.rebinX>0 or args.rebinY>0):
        command += ' --rebinX='+str(args.rebinX)+' --rebinY='+str(args.rebinY)
    for hel in ['0','1','2','3','4']:
        if ( hasattr(args,"scale"+hel) and getattr(args,"scale"+hel) ):
            command += (' --scale'+hel)
    if args.jacmass>-1:
        command += (' --jacmass='+str(args.jacmass))
            
    fit_opts = [' --jUL', 
                ' --jUL --j0', 
                ' --jUL --j0 --j1', 
                ' --jUL --j0 --j1 --j2',
                ' --jUL --j0 --j1 --j2 --j3',
                ' --jUL --j0 --j1 --j2 --j4',
                ' --jUL --j0 --j1 --j2 --j3 --j4',
            ]
    for i in fit_opts:
        new_command = command+i
        print(new_command)
        os.system(new_command) 

elif args.algo=='fit_systs':
    for syst in pol_systs:
        command  = './fit --nevents='+str(args.nevents) +' --tag='+args.tag+' --run=closure --post_tag='+args.post_tag
        for k in syst.keys():
            command += ' --degs_'+k+'_x='+str(syst[k][0])+' --degs_'+k+'_y='+str(syst[k][1])
        if (args.rebinX>0 or args.rebinY>0):
            command += ' --rebinX='+str(args.rebinX)+' --rebinY='+str(args.rebinY)
        for hel in ['0','1','2','3','4']:
            if ( hasattr(args,"scale"+hel) and getattr(args,"scale"+hel) ):
                command += (' --scale'+hel)
        if args.jacmass>-1:
            command += (' --jacmass='+str(args.jacmass))

        fit_opts = [' --jUL', 
                    ' --jUL --j0', 
                    ' --jUL --j0 --j1', 
                    ' --jUL --j0 --j1 --j2',
                    ' --jUL --j0 --j1 --j2 --j3',
                    ' --jUL --j0 --j1 --j2 --j4',
                    ' --jUL --j0 --j1 --j2 --j3 --j4',
                ]
        for i in fit_opts:
            new_command = command+i
            print(new_command)
            os.system(new_command) 
        
elif args.algo=='fit_fast':
    command  = './fit --nevents='+str(args.nevents) +' --tag='+args.tag+' --run=closure --post_tag='+args.post_tag
    for k in pol_default.keys():
        command += ' --degs_'+k+'_x='+str(pol_default[k][0])+' --degs_'+k+'_y='+str(pol_default[k][1])
    if (args.rebinX>0 or args.rebinY>0):
        command += ' --rebinX='+str(args.rebinX)+' --rebinY='+str(args.rebinY)
    fit_opts = ['jUL', 'j0', 'j1', 'j2', 'j3', 'j4']
    #fit_opts = ['jUL', 'j0', 'j1']
    for i in range(len(fit_opts)):
        command += (' --'+fit_opts[i])
    print(command)
    os.system(command) 

elif args.algo=='fit_grid_fast':
    command  = './fit_grid --nevents='+str(args.nevents) +' --tag='+args.tag+' --run=grid --post_tag='+args.post_tag
    for k in pol_default.keys():
        if 'corr' not in k: continue
        command += ' --degs_'+k+'_x='+str(pol_default[k][0])+' --degs_'+k+'_y='+str(pol_default[k][1])
    if (args.rebinX>0 or args.rebinY>0):
        command += ' --rebinX='+str(args.rebinX)+' --rebinY='+str(args.rebinY)
    fit_opts = ['jUL', 'j0', 'j1', 'j2', 'j3', 'j4']
    #fit_opts = ['jUL', 'j0', 'j1']
    for i in range(len(fit_opts)):
        command += (' --'+fit_opts[i])
    print(command)
    os.system(command) 

elif args.algo=='fit_grid':
    command  = './fit_grid --nevents='+str(args.nevents) +' --tag='+args.tag+' --run=grid --post_tag='+args.post_tag
    for k in pol_default.keys():
        if 'corr' not in k: continue
        command += ' --degs_'+k+'_x='+str(pol_default[k][0])+' --degs_'+k+'_y='+str(pol_default[k][1])
    if (args.rebinX>0 or args.rebinY>0):
        command += ' --rebinX='+str(args.rebinX)+' --rebinY='+str(args.rebinY)
    fit_opts = [' --jUL', 
                ' --jUL --j0', 
                ' --jUL --j0 --j1', 
                ' --jUL --j0 --j1 --j2',
                ' --jUL --j0 --j1 --j2 --j3',
                ' --jUL --j0 --j1 --j2 --j4',
                ' --jUL --j0 --j1 --j2 --j3 --j4',
            ]
    for i in fit_opts:
        new_command = command+i
        print(new_command)
        os.system(new_command) 
