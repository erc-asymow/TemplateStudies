import argparse
import os
import sys
import copy

parser = argparse.ArgumentParser(description='run')

parser.add_argument('--none', action='store_true'  , help = 'none')
parser.add_argument('--tag',  default='test'       , help = 'tag')
parser.add_argument('--post_tag',  default='test'  , help = 'post tag')
parser.add_argument('--algo', default='all'        , help = 'algo')
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

'''
pol_default = {
    'corr'   : [10,4],
    'A0'     : [3,3],
    'A1'     : [3,3],
    'A2'     : [3,3],
    'A3'     : [3,3],
    'A4'     : [3,3]
}
'''
pol_default = {
    'corr'   : [10,4],
    'A0'     : [3,3],
    'A1'     : [3,3],
    'A2'     : [3,3],
    'A3'     : [3,3],
    'A4'     : [3,3]
}

pol_systs = []
pol_syst0 = copy.deepcopy(pol_default)
pol_syst0['corr'] = [8,4]
#pol_systs.append( pol_syst0 )
pol_syst1 = copy.deepcopy(pol_default)
pol_syst1['corr'] = [10,6]
pol_systs.append( pol_syst1 )
pol_syst2 = copy.deepcopy(pol_default)
pol_syst2['A3'] = [3,4]
#pol_systs.append( pol_syst2 )
pol_syst3 = copy.deepcopy(pol_default)
pol_syst3['A4'] = [3,4]
#pol_systs.append( pol_syst3 )

if args.algo=='jac2':
    command  = './jac2 --nevents='+str(args.nevents) +' --tag='+args.tag+' --run=closure'
    for k in pol_default.keys():
        command += ' --degs_'+k+'_x='+str(pol_default[k][0])+' --degs_'+k+'_y='+str(pol_default[k][1])
    command += ' --toyTF2_corr'
    print(command)
    os.system(command)   

elif args.algo=='jac2_systs':
    for syst in pol_systs:
        command  = './jac2 --nevents='+str(args.nevents) +' --tag='+args.tag+' --run=closure'
        for k in syst.keys():
            command += ' --degs_'+k+'_x='+str(syst[k][0])+' --degs_'+k+'_y='+str(syst[k][1])
        command += ' --toyTF2_corr'
        print(command)
        os.system(command)   

elif args.algo=='jac3':
    command  = './jac3 --nevents='+str(args.nevents) +' --tag='+args.tag+' --run=grid'
    for k in pol_default.keys():
        if 'corr' not in k: continue
        command += ' --degs_'+k+'_x='+str(pol_default[k][0])+' --degs_'+k+'_y='+str(pol_default[k][1])
    command += ' --toyTF2_corr'
    print(command)
    os.system(command)   

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
